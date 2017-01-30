// Platform check
#ifdef PLATFORM_CF1
#error ESTIMATOR = kalman is only compatible with the Crazyflie 2.0
#endif

#include "estimator_ME_dynamical.h"

#include "stm32f4xx.h"
#include "FreeRTOS.h"
#include "queue.h"
#include "task.h"
#include "sensors.h"
#include "log.h"
#include "param.h"
#include "math.h"
#include "arm_math.h"

/******************************************************************************
 * Flter settings, (1) can enable automatic resetting when infeasible NAN
 * values are present in the covariance matrices, (2) enforce symmetry of the
 * covariance matrices (3) run multiple estimators (EKF/UKF/GPF) enable various
 * means of multilateration through TOA WLS or by estimating the clock skew
 * by means of an independent kalman filter (as in the work of Mike Hamer)
 *****************************************************************************/
#define FILTER_NAN_CHECK
#define FILTER_SYMMETRY_CHECK
//#define USE_UKF_FILTER
//#define USE_ESKF_FILTER
#define USE_EKF_FILTER
#define USE_WLS_LATERATION           // TOA
//#define USE_LS_LATERATION          // TOA
//#define USE_KALMAN_SKEW_LATERATION // TDOA

// Sanity check to see that at least one of the filters have been defined
#if !defined(USE_EKF_FILTER) && !defined(USE_UKF_FILTER) && !defined(USE_ESKF_FILTER)
#error Either the EKF or UKF must be explicitly chosen
#endif

// Sanity check to see that at least one of the multilaterative algorithms have been defined
#if !defined(USE_WLS_LATERATION) && !defined(USE_LS_LATERATION) && !defined(USE_KALMAN_SKEW_LATERATION)
#error Either the multilateration algorithm must be explicitly chosen
#endif

// Tuning parameters
#define FILTER_RATE RATE_100_HZ // this is slower than the IMU update rate of 500Hz
// TODO include more...

// Constants used in the estimator
#define DEG_TO_RAD (PI/180.0f)
#define RAD_TO_DEG (180.0f/PI)
#define GRAVITY_MAGNITUDE (9.81f)
#define CRAZYFLIE_WEIGHT_grams (27.0f)
#define CONTROL_TO_ACC (GRAVITY_MAGNITUDE*60.0f/CRAZYFLIE_WEIGHT_grams/65536.0f)
#define SPEED_OF_LIGHT (299792458)
#define METERS_PER_TDOATICK (4.691763979e-3f)
#define SECONDS_PER_TDOATICK (15.650040064e-12f)

// Initial variances, uncertain of position, but know we're stationary and roughly flat
static const float stdDevInitialPosition_xy = 100;
static const float stdDevInitialPosition_z = 1;
static const float stdDevInitialVelocity = 0.01;
static const float stdDevInitialAttitude = 0.01;
static const float stdDevInitialAngularRate = 0.01;
static const float stdDevInitialSkew = 0.1;

// Process variances
static float procNoiseVel = 0;
static float procNoisePos = 0;
static float procNoiseAtt = 0;
static float procNoiseRate = 0;
static float procNoiseSkew = 10e-6f; // seconds per second^2 (is multiplied by dt to give skew noise)

/******************************************************************************
 * These funcitons have been appropriated and re-written the work of Hamer to
 * allow inclusion of positional measurements with the addition of a TOA queue
 * for static multilateration 
 *****************************************************************************/
static xQueueHandle distDataQueue;
#define DIST_QUEUE_LENGTH (10)

static void stateEstimatorMultilateration(TOAArray_t *TOAArray, Axis3f *measuredPsoition);

static inline bool stateEstimatorHasTOAArray(TOAArray_t *TOAArray) {
  return (pdTRUE == xQueueReceive(distDataQueue, TOAArray, 0));
}

/******************************************************************************
 * These funcitons have been appropriated from the work of Hamer to allow
 * inclusion of positional measurements with the addition of a TOA queue for
 * multilateration 
 *****************************************************************************/
// (A) - Predicting the current state forward
static void dynamicalFilterPrediction();
static void dynamicalFilterProcessNoise();

// (B) - Measurement update and correction based on sensors
static void dynamicalEstimatorCorrection();

// (C) - Rewrite filter state into the format expected by other modules */
static void dynamicalFilterConversion();


// State vector and covariance matrix
// The internally-estimated state is:
// - X, Y, Z: the quad's position in the global frame
// - PX, PY, PZ: the quad's velocity in the global frame
// - qw, qx, qy, qz: attitude quaternion body to global frame
// - wx, wy, wz: Angular rates in the body
// stored as a column vector
typedef enum
{
  IND_PX, IND_PY, IND_PZ, IND_VX, IND_VY, IND_VZ, IND_QW, IND_QX, IND_QY, IND_QZ, IND_WX, IND_WY, IND_WZ, STATE_DIMENSION
} stateVectorIdx_t;
static float S[STATE_DIMENSION];
static float P[STATE_DIMENSION][STATE_DIMENSION];
static arm_matrix_instance_f32 Pm = {STATE_DIMENSION, STATE_DIMENSION, (float *)P};

// Internal variables
static bool isInit = false;
static bool resetEstimation = true;
static float thrustAccumulator;
static Axis3f accAccumulator;
static Axis3f gyroAccumulator;
static Axis3f measuredPosition;
static uint32_t accAccumulatorCount;
static uint32_t thrustAccumulatorCount;
static uint32_t gyroAccumulatorCount;
static bool quadIsFlying = false;
static float stateSkew;
static float varSkew;
static uint32_t lastPrediction;
static uint32_t lastFlightCmd;
static uint32_t takeoffTime;
static uint32_t tdoaCount;
static int32_t lastTDOAUpdate;
static float LS_Pseudoinverse[3][3];

/******************************************************************************
 * Supporting utility functions for (i) matrix multiplication, (ii) solving
 * linear systems with Cholesky decomposition (A = LL*) or (ii) gaussian 
 * elimination with partial pivoting (A = LUP) followed by forward and backward
 * substitution. The method also contains a NaN check to see when the system
 * state or covariance is infeasible and should be reset.
 *****************************************************************************/
static inline void mat_trans(const arm_matrix_instance_f32 * pSrc, arm_matrix_instance_f32 * pDst)
{ configASSERT(ARM_MATH_SUCCESS == arm_mat_trans_f32(pSrc, pDst)); }
static inline void mat_inv(const arm_matrix_instance_f32 * pSrc, arm_matrix_instance_f32 * pDst)
{ configASSERT(ARM_MATH_SUCCESS == arm_mat_inverse_f32(pSrc, pDst)); }
static inline void mat_mult(const arm_matrix_instance_f32 * pSrcA, const arm_matrix_instance_f32 * pSrcB, arm_matrix_instance_f32 * pDst)
{ configASSERT(ARM_MATH_SUCCESS == arm_mat_mult_f32(pSrcA, pSrcB, pDst)); }
static inline float arm_sqrt(float32_t in)
{ float pOut = 0; arm_status result = arm_sqrt_f32(in, &pOut); configASSERT(ARM_MATH_SUCCESS == result); return pOut; }

#define CHOL_MAX_ELEMENTS 4 // Maximum dimasion allowed
#define LUP_MAX_ELEMENTS 4  // Maximum dimasion allowed

static void cholSolve(float (*ATA)[CHOL_MAX_ELEMENTS], float *b, float *x, int N)
{
  // LINEAR SYSTEM SOLVER BY CHOLESKY DECOMPOSITION
  // Written by Marcus Greiff 15/10/2016 based on Numerical Linear
  // Algebra, Lloyd N. Trefethen page 175.
  //
  // This is a fast and robust method of solving a linear system
  // ATA*x = b where ATA is a is a square, symmetric and positive
  // definite matrix (A^TA). The algorithm uses the Cholesky decomposition
  // in combination with forward and backward substitution, preserving ATA
  // but overwriting both b and x.
  //
  // No feasibility checks are made, so the user has to check/prove
  // positive definiteness of ATA beforehand and be sure to set 
  // CHOL_MAX_ELEMENTS larger than or equal to the highest intended N.
  float R[CHOL_MAX_ELEMENTS][CHOL_MAX_ELEMENTS] = {0};
  float y[CHOL_MAX_ELEMENTS] = {0};
  int i, j, k;
  for (i = 0; i < N; i++){
    for (j = i; j < N; j++){
      R[i][j] = ATA[i][j];
    }
  }
  
  // Decomposition
  for (k = 0; k < N; k++){
    for (j = k+1; j<N; j++){
      for (i = j; i < N; i++){
        R[j][i] -= R[k][i]*R[k][j]/R[k][k];
      }
    }
    for (i = N-1; i>=k; i--){
      R[k][i] /= sqrtf(R[k][k]);
    }
  }

  // Forward sub
  for (i = 0; i < N; i++){
    y[i] = b[i]/R[i][i];
    for (j = i + 1; j < N; j++){
      b[j] -= R[i][j]*y[i];
    }
  }

  // Backward sub
  for (i = N - 1; i >=0 ; i--){
    x[i] = y[i]/R[i][i];
    for (j = 0; j < i ; j++){
      y[j] -= R[j][i]*x[i];
    }
  }
}

#ifdef FILTER_NAN_CHECK
static void dynamicalEstimatorAssertNotNaN() {
  if ((isnan(S[IND_PX])) ||
      (isnan(S[IND_PY])) ||
      (isnan(S[IND_PZ])) ||
      (isnan(S[IND_VX])) ||
      (isnan(S[IND_VY])) ||
      (isnan(S[IND_VZ])) ||
      (isnan(S[IND_QW])) ||
      (isnan(S[IND_QX])) ||
      (isnan(S[IND_QY])) ||
      (isnan(S[IND_QZ])) ||
      (isnan(S[IND_WX])) ||
      (isnan(S[IND_WY])) ||
      (isnan(S[IND_WZ]))) {
      resetEstimation = true;
    }
  for(int i=0; i<STATE_DIMENSION; i++) {
    for(int j=0; j<STATE_DIMENSION; j++) {
      if (isnan(P[i][j])) {
        resetEstimation = true;
      }
    }
  }
}
#else
static void dynamicalEstimatorAssertNotNaN()
{
  return false;
}
#endif

void stateEstimatorUpdate(state_t *state, sensorData_t *sensors, control_t *control)
{
  /****************************************************************************
   * Main state estimator callback, with the general structure appropriated
   * from the work of Hamer
   ***************************************************************************/
   uint32_t tick=1;
  // If the client (via a parameter update) triggers an estimator reset
  if (resetEstimation) { stateEstimatorInit(); resetEstimation = false; }
  
  // Check if an update has been made and finalizes the current state
  bool doneUpdate = false;

  // Average the last IMU measurements.
  if (sensorsReadAcc(&sensors->acc)) {
    accAccumulator.x += sensors->acc.x;
    accAccumulator.y += sensors->acc.y;
    accAccumulator.z += sensors->acc.z;
    accAccumulatorCount++;
  }
  if (sensorsReadGyro(&sensors->gyro)) {
    gyroAccumulator.x += sensors->gyro.x;
    gyroAccumulator.y += sensors->gyro.y;
    gyroAccumulator.z += sensors->gyro.z;
    gyroAccumulatorCount++;
  }

  // Average the thrust command from the last timestep, generated externally by the controller
  thrustAccumulator += control->thrust * CONTROL_TO_ACC; // thrust is in grams, we need ms^-2
  thrustAccumulatorCount++;

  // Prediction and correction
  if (RATE_DO_EXECUTE(FILTER_RATE, tick)
      && gyroAccumulatorCount > 0
      && accAccumulatorCount > 0
      && thrustAccumulatorCount > 0)
  {
    // Accelerometer is in Gs but the estimator requires ms^-2
    gyroAccumulator.x *= (GRAVITY_MAGNITUDE / gyroAccumulatorCount);
    gyroAccumulator.y *= (GRAVITY_MAGNITUDE / gyroAccumulatorCount);
    gyroAccumulator.z *= (GRAVITY_MAGNITUDE / gyroAccumulatorCount);

    // Gyro is in deg/sec but the estimator requires rad/sec
    accAccumulator.x *= (DEG_TO_RAD / gyroAccumulatorCount);
    accAccumulator.y /= (DEG_TO_RAD / gyroAccumulatorCount);
    accAccumulator.z /= (DEG_TO_RAD / gyroAccumulatorCount);
    thrustAccumulator /= thrustAccumulatorCount;

    // If a new ranging measurement has been done, the position is determined
    // through a multilateration and the 
    TOAArray_t TOAArray;
    if (stateEstimatorHasTOAArray(&TOAArray))
    {
      stateEstimatorMultilateration(&TOAArray, &measuredPosition);
    }

    // Compute the filter time step
    float dt = 1.0/(float)FILTER_RATE;
    
    // State estimator prediction
    dynamicalFilterPrediction(&control);
    
    // State estimator update
    dynamicalEstimatorCorrection(); // TODO include more parameters

    // Set to true if update has been done successully
    doneUpdate = true;
  }
  
  // TODO Write method for adding process noise, should be used in the
  float dt = 1/RATE_MAIN_LOOP;
  dynamicalFilterProcessNoise();

  // TODO Handle incoming measurements

  // Finalize internal snante state if an update has been made
  if (doneUpdate)
  {
    dynamicalFilterConversion();
    dynamicalEstimatorAssertNotNaN();
  }
}

static void dynamicalFilterPrediction(state_t *state, control_t *control)
{
  /****************************************************************************
   * In this method, the quaternion euler dynamics are used to predict the
   * state and covariance matrix forward. This is done in three separate ways
   * depending on the filter employed and is solely based on the current state
   * and control signals.
   ***************************************************************************/
  dynamicalEstimatorAssertNotNaN();
}

static void dynamicalFilterProcessNoise()
{
  dynamicalEstimatorAssertNotNaN();
}


static void dynamicalEstimatorCorrection()
{
  dynamicalEstimatorAssertNotNaN();
}


static void dynamicalFilterConversion()
{
  dynamicalEstimatorAssertNotNaN();
}

void stateEstimatorInit(void) {

  if (!isInit)
  {
    distDataQueue = xQueueCreate(DIST_QUEUE_LENGTH, sizeof(TOAArray_t));
  }
  else
  {
    xQueueReset(distDataQueue);
  }

  lastPrediction = xTaskGetTickCount();
  lastTDOAUpdate = xTaskGetTickCount();

  accAccumulator = (Axis3f){.axis={0}};
  gyroAccumulator = (Axis3f){.axis={0}};
  thrustAccumulator = 0;

  accAccumulatorCount = 0;
  gyroAccumulatorCount = 0;
  thrustAccumulatorCount = 0;
  
  S[IND_PX] = 0.5;
  S[IND_PY] = 0.5;
  S[IND_PZ] = 0;
  S[IND_VX] = 0;
  S[IND_VY] = 0;
  S[IND_VZ] = 0;
  S[IND_QW] = 1.0;
  S[IND_QX] = 0;
  S[IND_QY] = 0;
  S[IND_QZ] = 0;
  S[IND_WX] = 0;
  S[IND_WY] = 0;
  S[IND_WZ] = 0;

  // Initialize covariance matrix
  for (int i=0; i< STATE_DIMENSION; i++) {
    for (int j=0; j < STATE_DIMENSION; j++) {
      P[i][j] = 0;
    }
  }

  P[IND_PX][IND_PX] = powf(stdDevInitialPosition_xy, 2);
  P[IND_PY][IND_PY] = powf(stdDevInitialPosition_xy, 2);
  P[IND_PZ][IND_PZ] = powf(stdDevInitialPosition_z, 2);
  
  P[IND_VX][IND_VX] = powf(stdDevInitialVelocity, 2);
  P[IND_VY][IND_VY] = powf(stdDevInitialVelocity, 2);
  P[IND_VZ][IND_VZ] = powf(stdDevInitialVelocity, 2);
  
  P[IND_QW][IND_QW] = powf(stdDevInitialPosition_xy, 2);
  P[IND_QX][IND_QX] = powf(stdDevInitialPosition_xy, 2);
  P[IND_QY][IND_QY] = powf(stdDevInitialPosition_xy, 2);
  P[IND_QZ][IND_QZ] = powf(stdDevInitialPosition_xy, 2);

  P[IND_WX][IND_WX] = powf(stdDevInitialVelocity, 2);
  P[IND_WY][IND_WY] = powf(stdDevInitialVelocity, 2);
  P[IND_WZ][IND_WZ] = powf(stdDevInitialVelocity, 2);

  varSkew = powf(stdDevInitialSkew, 2);
  
  tdoaCount = 0;
  isInit = true;
}

static void stateEstimatorMultilateration(TOAArray_t *TOAArray, Axis3f *measuredPosition)
{
#if defined(USE_WLS_LATERATION)
  // Weighted LS multilateration
  // Solve argmin[(A\tilde{p} - b)W(A\tilde{p} - b)], where
  //
  // Robot position: \mathbf{p} = [px,py,pz]\in\mathbb{R}^{3\times1}
  // Anchor (i) position: \mathbf{p}_i
  // Optimization vector: \tide{\mathbf{p}} = [px,py,pz,||\mathbf{p}||_2^2]
  // 
  // From which \mathbf{p} is extracted
  int WLS_N = 4;
  float WLS_tildeP[4];
  float WLS_A[6][4];
  float WLS_AWA[4][4] = {0};
  float WLS_B[6];
  float WLS_W[6];
  float WLS_AWB[4] = {0};
  
  // Temporary variables
  float WLS_di2;
  float WLS_pi[3];
  float WLS_pi2;
  int ii, jj, kk;
  
  for (ii = 0; ii < TOAArray->numberOfAnchors; ii++){
    WLS_pi[0] = TOAArray->anchorPosition[ii].x;
    WLS_pi[1] = TOAArray->anchorPosition[ii].y;
    WLS_pi[2] = TOAArray->anchorPosition[ii].z;
    WLS_pi2 = WLS_pi[0]*WLS_pi[0] + WLS_pi[1]*WLS_pi[1] + WLS_pi[2]*WLS_pi[2];
    WLS_di2 = powf(TOAArray->distances[ii], 2);
    if (WLS_di2 < 0.01){
      WLS_di2 = 0.01; // Condition for retaining feasibility in Cholesky, making AWA PSD
    }

    // Form a_i
    WLS_A[ii][0] = -2 * WLS_pi[0];
    WLS_A[ii][1] = -2 * WLS_pi[1];
    WLS_A[ii][2] = -2 * WLS_pi[2];
    WLS_A[ii][3] = 1.0;

    // form b_i
    WLS_B[ii] = WLS_di2 - WLS_pi2;

    // form w_i
    WLS_W[ii] = 1.0 / (4.0 * WLS_di2);
  }
  
  // form (A^T*W*A)
  for (ii = 0; ii < WLS_N; ii++){
    for (jj = 0; jj < WLS_N; jj++){
      for (kk = 0; kk < TOAArray->numberOfAnchors; kk++){
        WLS_AWA[ii][jj] += WLS_A[kk][ii] * WLS_W[kk] * WLS_A[kk][jj];
      }
    }
  }

  // form (A^T*W*B)
  for (ii = 0; ii < WLS_N; ii++){
    for (kk = 0; kk < TOAArray->numberOfAnchors; kk++){
      WLS_AWA[ii][jj] += WLS_A[kk][ii] * WLS_W[kk] * WLS_B[kk];
    }
  }
  
  // Solve linear system
  cholSolve(WLS_AWA, WLS_AWB, WLS_tildeP, WLS_N);
  
  // Extract position
  measuredPosition->x = WLS_tildeP[0];
  measuredPosition->y = WLS_tildeP[1];
  measuredPosition->z = WLS_tildeP[2];
#else
  // TODO: Regular LS multilateration
  measuredPosition->x = 1.0;
  measuredPosition->y = 2.0;
  measuredPosition->z = 3.0;
#endif
}

/******************************************************************************
 * Basic functionality for enqueueing external measurements, TOA, TDOA, and
 * MOCCAP for use in the dynamical filter
 *****************************************************************************/
static bool enqueueExternalMeasurement(xQueueHandle queue, void *measurement)
{
  portBASE_TYPE result;
  bool isInInterrupt = (SCB->ICSR & SCB_ICSR_VECTACTIVE_Msk) != 0;
  
  if (isInInterrupt) {
    portBASE_TYPE xHigherPriorityTaskWoken = pdFALSE;
    result = xQueueSendFromISR(queue, measurement, &xHigherPriorityTaskWoken);
    if(xHigherPriorityTaskWoken == pdTRUE)
    {
      portYIELD();
    }
  } else {
    result = xQueueSend(queue, measurement, 0);
  }
  return (result==pdTRUE);
}

bool stateEstimatorEnqueueTOAArray(TOAArray_t *TOAArray)
{
  return enqueueExternalMeasurement(distDataQueue, (void *)TOAArray);
}


/******************************************************************************
 * Functionality for testing the estimatior
 *****************************************************************************/
bool stateEstimatorTest(void)
{
  // TODO Include test for cholesky solver
  // TODO include test for LUP solver
  return isInit;
}
