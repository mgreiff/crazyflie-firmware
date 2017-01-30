// Note that FreeRTOS must be included to use the timers, and that the
// script fails if they are included in the wrong order.
#include "M_sequenceCommander.h"

#include "stm32f4xx.h"

#include "FreeRTOS.h"
#include "timers.h"
#include "queue.h"
#include "param.h"
#include "log.h"
#include "math.h"
#include "num.h"
#include "math.h"

#define REFERENCE_RATE RATE_25_HZ
#define REFERENCE_OFFSET 7
#define FLATNESS_RATE RATE_50_HZ
#define FLATNESS_OFFSET 14
#define LOAD_RATE RATE_10_HZ
#define READ_TRAJECTORY_QUEUE

#define LP_FILTER_LENGTH 5
#define POLYNOMIAL_DEGREE 5
#define MAX_TRAJECTORY_TYPE 5 // The total number of trajectory types (for assertion)

//#define USE_TEST_SEQUENCE_A // <- Uncomment to run a simple test sequence
//#define USE_TEST_SEQUENCE_B // <- Uncomment to run a simple test sequence
//#define USE_TEST_SEQUENCE_C // <- Uncomment to run a simple test sequence

static float flatOutput[4][5];          // Used for debugging and plotting
trajectory_t trajObj[NUM_FLAT_OUTPUTS]; // Stores the state trajectory

// Structure for preserving the LP filter history
static struct{
  float Y[FLAT_OUTPUT_DIMENSIONS][LP_FILTER_LENGTH];
  float U[FLAT_OUTPUT_DIMENSIONS][LP_FILTER_LENGTH];
  uint32_t ptr;
} filterHistory[NUM_FLAT_OUTPUTS];

// LP coefficients for derivative determination at 50 Hz
const float LP_denominator[LP_FILTER_LENGTH] = {
   1.0f                 , -3.441860465116279f  , 4.442401297998918f   , -2.548354232960619f  , 0.548192480346180f      // Denominator
};
const float LP_numerator[FLAT_OUTPUT_DIMENSIONS][LP_FILTER_LENGTH]   = {
  {0.000023692516762456f, 0.000094770067049822f, 0.000142155100574734f, 0.000094770067049822f, 0.000023692516762456f}, // Position
  {   0.002369251676246f,    0.004738503352491f,                  0.0f,  - 0.004738503352491f,   -0.002369251676246f}, // Velocity
  {   0.236925167624556f,                  0.0f,   -0.473850335249112f,                  0.0f,    0.236925167624556f}, // Acceleration
  {  23.692516762455604f,  -47.385033524911208f,                  0.0f,   47.385033524911208f,  -23.692516762455604f}, // Jerk
  {    2369.25167624556f,    -9477.00670498224f,    14215.51005747336f,    -9477.00670498224f,     2369.25167624556f}  // Snap (jerk derivative)
};

// Log variables for debugging (to be removed)
static uint32_t testCounter[6] = {0};
static float testVal[6];
static uint16_t nPackagesReceived;

// ~~~ Test sequence A data ~~~
const float testTrajectoryX[8][6] = {
  {     0.0f,       0.0f,       0.0f,   -0.1813f,    0.0594f,       0.0f},
  { -0.5000f,   -0.2748f,    0.3377f,   -0.0058f,   -0.0159f,       0.0f},
  {     0.0f,    0.4967f,   -0.0795f,   -0.0306f,    0.0043f,       0.0f},
  {  0.5000f,   -0.0497f,   -0.1589f,    0.0323f,   -0.0014f,       0.0f},
  {     0.0f,   -0.3444f,    0.0000f,    0.0207f,    0.0014f,       0.0f},
  { -0.5000f,   -0.0497f,    0.1589f,    0.0041f,   -0.0043f,       0.0f},
  {     0.0f,    0.4967f,    0.0795f,   -0.1333f,    0.0159f,       0.0f},
  {  0.5000f,   -0.2748f,   -0.3377f,    0.2939f,   -0.0594f,       0.0f}
};

const float testTrajectoryY[8][6] = {
  {     0.0f,       0.0f,       0.0f,    0.1304f,   -0.0339f,       0.0f},
  {  0.5000f,    0.4785f,   -0.0323f,   -0.0592f,    0.0091f,       0.0f},
  {  1.0000f,   -0.0695f,   -0.1689f,    0.0443f,   -0.0025f,       0.0f},
  {  0.5000f,   -0.2930f,    0.0373f,   -0.0095f,    0.0008f,       0.0f},
  {     0.0f,   -0.2318f,   -0.0000f,   -0.0029f,   -0.0008f,       0.0f},
  { -0.5000f,   -0.2930f,   -0.0373f,    0.0244f,    0.0025f,       0.0f},
  { -1.0000f,   -0.0695f,    0.1689f,    0.0137f,   -0.0091f,       0.0f},
  { -0.5000f,    0.4785f,    0.0323f,   -0.1411f,    0.0339f,       0.0f}
};
 
const float testTrajectoryZ[8][6] = {
  {     0.0f,       0.0f,       0.0f,   -0.0172f,    0.0086f,       0.0f},
  {     0.0f,    0.0687f,    0.1030f,   -0.0328f,   -0.0023f,       0.0f},
  {  0.2500f,    0.0133f,   -0.1491f,    0.0388f,    0.0006f,       0.0f},
  {     0.0f,   -0.0985f,    0.0980f,   -0.0241f,   -0.0001f,       0.0f},
  {     0.0f,   -0.0000f,   -0.0497f,    0.0251f,   -0.0001f,       0.0f},
  {     0.0f,    0.0985f,    0.0980f,   -0.0436f,    0.0006f,       0.0f},
  {  0.2500f,   -0.0133f,   -0.1491f,    0.0512f,   -0.0023f,       0.0f},
  {     0.0f,   -0.0687f,    0.1030f,   -0.0515f,    0.0086f,       0.0f}
};

const float testTrajectoryTime[8] = {2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f};

void sequeceCommanderTestA(){
  // Loads the test XYZ trajectory on the desired data fromat
  for (int r = 0; r < 8; r++){
    // Coefficients
    for (int c = 0; c < 6; c++){
      trajObj[0].data[r][c] = testTrajectoryX[r][c]; 
      trajObj[1].data[r][c] = testTrajectoryY[r][c]; 
      trajObj[2].data[r][c] = testTrajectoryZ[r][c];
    }
    // Times
    for (int dim = 0; dim < 3; dim++){
      trajObj[dim].time[r] = testTrajectoryTime[r];
      trajObj[dim].type[r] = 1; // polynomial
    }
  }
  // Ending with a spline will cause the trajecotry to diverge,
  // we must therefore en the trajectory with a hovering point
  for (int dim = 0; dim < 3; dim++){
    trajObj[dim].data[8][0] = 2.0f;
    trajObj[dim].time[8] = 5.0f;
    trajObj[dim].type[8] = 0; // LP
    trajObj[dim].numberOfEntries = 9;
    trajObj[dim].status = true;
    trajObj[dim].startTime = 5.0f;
  }
  // Loads the sinusodial yaw starting at time 0
  trajObj[3].data[0][0] = 0.2f; // Amplitude
  trajObj[3].data[0][1] = 0.0f;
  trajObj[3].data[0][2] = 0.5f; // Frequency
  trajObj[3].data[0][3] = 0.0f;
  trajObj[3].time[0] = 30.0f;
  trajObj[3].numberOfEntries = 1;
  trajObj[3].startTime = 0.0f;
  trajObj[3].type[0] = 2;    // function
  trajObj[3].status = true;
}

// ~~~ Test sequence B data ~~~
void sequeceCommanderTestB(){
  trajObj[0].data[0][0] = 1.0f;
  trajObj[0].data[1][0] = 0.0f;
  trajObj[0].data[2][0] = 2.0f;
  trajObj[0].data[3][0] = 0.5f;
  trajObj[0].time[0] = 5.0f;
  trajObj[0].time[1] = 5.0f;
  trajObj[0].time[2] = 5.0f;
  trajObj[0].time[3] = 5.0f;
  trajObj[0].type[0] = 0;
  trajObj[0].type[1] = 0;
  trajObj[0].type[2] = 0;
  trajObj[0].type[3] = 0;
  trajObj[0].numberOfEntries = 4;
  trajObj[0].status = true;
  trajObj[0].startTime = 5.0f;

  trajObj[1].data[0][0] = 0.1f; // Amplitude
  trajObj[1].data[0][1] = 0.0f;
  trajObj[1].data[0][2] = 1.0f; // Frequency
  trajObj[1].data[0][3] = 0.0f;
  trajObj[1].time[0] = 30.0f;
  trajObj[1].numberOfEntries = 1;
  trajObj[1].startTime = 0.0f;
  trajObj[1].type[0] = 2;       // function
  trajObj[1].status = true;

  trajObj[2].data[0][0] = 0.2f; // Amplitude
  trajObj[2].data[0][1] = 0.0f;
  trajObj[2].data[0][2] = 1.0f; // Frequency
  trajObj[2].data[0][3] = 0.0f;
  trajObj[2].time[0] = 30.0f;
  trajObj[2].numberOfEntries = 1;
  trajObj[2].startTime = 0.0f;
  trajObj[2].type[0] = 2;       // function
  trajObj[2].status = true;

  trajObj[3].data[0][0] = 0.3f; // Amplitude
  trajObj[3].data[0][1] = 0.0f;
  trajObj[3].data[0][2] = 1.0f; // Frequency
  trajObj[3].data[0][3] = 0.0f;
  trajObj[3].time[0] = 30.0f;
  trajObj[3].numberOfEntries = 1;
  trajObj[3].startTime = 0.0f;
  trajObj[3].type[0] = 2;       // function
  trajObj[3].status = true;
}

// ~~~ Test sequence C ~~~
void sequeceCommanderTestC(){
  // X
  trajObj[0].data[0][0] = 2.0f;
  trajObj[0].data[0][1] = 2.8f;
  trajObj[0].data[0][2] = 3.0f;
  trajObj[0].data[0][3] = 3.4f;
  trajObj[0].data[1][0] = 3.4f;
  trajObj[0].data[1][1] = 3.8f;
  trajObj[0].data[1][2] = 4.0f;
  trajObj[0].data[1][3] = 4.2f;
  trajObj[0].data[2][0] = 4.2f;
  trajObj[0].data[2][1] = 4.4f;
  trajObj[0].data[2][2] = 4.0f;
  trajObj[0].data[2][3] = 3.6f;

  // Y
  trajObj[1].data[0][0] = 1.5f;
  trajObj[1].data[0][1] = 1.5f;
  trajObj[1].data[0][2] = 2.1f;
  trajObj[1].data[0][3] = 1.5f;
  trajObj[1].data[1][0] = 1.5f;
  trajObj[1].data[1][1] = 0.9f;
  trajObj[1].data[1][2] = 0.9f;
  trajObj[1].data[1][3] = 1.5f;
  trajObj[1].data[2][0] = 1.5f;
  trajObj[1].data[2][1] = 2.1f;
  trajObj[1].data[2][2] = 1.7f;
  trajObj[1].data[2][3] = 1.7f;

  // Z
  trajObj[2].data[0][0] = 1.0f;
  trajObj[2].data[0][1] = 1.0f;
  trajObj[2].data[0][2] = 0.6f;
  trajObj[2].data[0][3] = 1.2f;
  trajObj[2].data[1][0] = 1.2f;
  trajObj[2].data[1][1] = 1.8f;
  trajObj[2].data[1][2] = 1.2f;
  trajObj[2].data[1][3] = 1.4f;
  trajObj[2].data[2][0] = 1.4f;
  trajObj[2].data[2][1] = 1.6f;
  trajObj[2].data[2][2] = 1.2f;
  trajObj[2].data[2][3] = 1.2f;

  // Times
  trajObj[0].time[0] = 5.0f;
  trajObj[1].time[0] = 5.0f;
  trajObj[2].time[0] = 5.0f;

  trajObj[0].time[1] = 6.0f;
  trajObj[1].time[1] = 6.0f;
  trajObj[2].time[1] = 6.0f;

  trajObj[0].time[2] = 4.0f;
  trajObj[1].time[2] = 4.0f;
  trajObj[2].time[2] = 4.0f;

  trajObj[1].startTime = 5.0f;
  trajObj[2].startTime = 5.0f;
  trajObj[3].startTime = 5.0f;

  // Number of entries
  trajObj[0].numberOfEntries = 3;
  trajObj[1].numberOfEntries = 3;
  trajObj[2].numberOfEntries = 3;

  // Type and status
  trajObj[0].type[0] = 3;
  trajObj[1].type[0] = 3;
  trajObj[2].type[0] = 3;
  trajObj[0].type[1] = 3;
  trajObj[1].type[1] = 3;
  trajObj[2].type[1] = 3;
  trajObj[0].type[2] = 3;
  trajObj[1].type[2] = 3;
  trajObj[2].type[2] = 3;
  trajObj[0].status = true;
  trajObj[1].status = true;
  trajObj[2].status = true;
  
  // Yaw reference
  trajObj[3].data[0][0] = 0.3f; // Amplitude
  trajObj[3].data[0][1] = 0.0f;
  trajObj[3].data[0][2] = 1.0f; // Frequency
  trajObj[3].data[0][3] = 0.0f;
  trajObj[3].time[0] = 30.0f;
  trajObj[3].numberOfEntries = 1;
  trajObj[3].startTime = 0.0f;
  trajObj[3].type[0] = 2;       // function
  trajObj[3].status = true;
}

// ~~~ The CRTP packet functions, structs and queques ~~~
#define TRAJECTORY_QUEUE_LENGTH (32)
#define CRTP_PORT_TRAJECTORY 0x08
#define CRTP_PORT_EMERGENCY 0x09
static xQueueHandle packetQueue;

static inline bool sequenceCommanderHasPacket(CRTPPacket *pk) {
  return (pdTRUE == xQueueReceive(packetQueue, pk, 0));
}

static void sequenceCommanderCRTPCallback(CRTPPacket* pk) {
  // Queues the incoming packet
  xQueueSend(packetQueue, pk, 0);
}

static void emergencyCRTPCallback(CRTPPacket* pk) {
  // Raises an assertion and shuts down the system, more encompassing and safe
  // than sendong zero in control reference and state given that setpoints and
  // states are updated at different times in the control system.
  configASSERT(false);
}

// TODO: Migrate to init function
static bool isInit = false;
void commanderInit_M() {
  if (isInit) {
    return;
  }
  packetQueue = xQueueCreate(TRAJECTORY_QUEUE_LENGTH, sizeof(CRTPPacket));
  xQueueReset(packetQueue);
  crtpRegisterPortCB(CRTP_PORT_TRAJECTORY, sequenceCommanderCRTPCallback);
  crtpRegisterPortCB(CRTP_PORT_EMERGENCY, emergencyCRTPCallback);
  isInit = true;
}

// ~~~ Commander functionality ~~~
void commanderGetSetpoint_M(M_setpoint_t *setpoint,
                            M_mode_t *mode,
                            M_state_t *state,
                            uint32_t currentTick)
{
  /****************************************************************************
   * This sequence commander loads trajectories thorugh the parameter framework,
   * allows the setting of controller modes, evaluate the trajectories and the
   * computation of the system states via the flatness equations.
   ***************************************************************************/

  // ~~~ Detect when switching modes ~~~
  if (mode->controller.current != mode->controller.previous){
    mode->controller.previous = mode->controller.current;
    mode->startTick = currentTick;
  }

#if defined(USE_TEST_SEQUENCE_A)
  // Used to test function/poly implementation and flatness
  sequeceCommanderTestA();
#elif defined(USE_TEST_SEQUENCE_B)
  // Used to test the LP filtering
  sequeceCommanderTestB();
#elif defined(USE_TEST_SEQUENCE_C)
  // Used to test the Bezier evaluation
  sequeceCommanderTestC();
#endif

  // ~~~ Decode and read CRTP packages based on the entered data type ~~~
  CRTPPacket packet;
  bool continueReading = true;
  while (sequenceCommanderHasPacket(&packet) == true && continueReading){
    // Load data into trajectory object and determine if the queue should be
    // or cleared based on the information in the packet
    continueReading = decodePacket(&packet, mode, setpoint);
    if (continueReading != true){
      // Clear queue and exit while loop
      xQueueReset(packetQueue);
    }
    // Accumulative package count for logging
    nPackagesReceived++;
  }
  
  // ~~~ Update references based on data ~~~
  if (trajObj[0].status == true &&
      trajObj[1].status == true &&
      trajObj[2].status == true &&
      trajObj[3].status == true)
    {
    if (RATE_DO_EXECUTE_OFFSET(REFERENCE_RATE, REFERENCE_OFFSET, currentTick)) {
      // Iterate over the flat outputs separately
      for (int dimension = 0; dimension < NUM_FLAT_OUTPUTS; dimension++){
        // Act if the trajectory is fully loaded and and has some entries
        if (trajObj[dimension].numberOfEntries > 0){

          // ~~~ Reads the cache and loads the trajectory if applicable ~~~
          // All trajectories are read fomr time 0 (at synchronization) and with a time lag
          uint32_t syncTick = trajObj[dimension].syncronizationTick;
          float startTime = trajObj[dimension].startTime;
          float currentTime = (1.0 / RATE_MAIN_LOOP) * (float)(currentTick - syncTick) + startTime;

          // Compute current index in data field
          int index = -1;
          float totalTime = trajObj[dimension].startTime;
          for (int ii = 0; ii < trajObj[dimension].numberOfEntries; ii++){
            if (totalTime < currentTime){
              totalTime += trajObj[dimension].time[ii];
              index++;
            }
          }

          // check if we are ahead of the specified starttime
          if (index > -1){
            // Time and index of the active trajectory was used
            float startTime = totalTime - trajObj[dimension].time[index];

            if (trajObj[dimension].type[index] == 0) {
              // LP filter point path
              eval_point_path(setpoint, index, dimension);
            } else if (trajObj[dimension].type[index] == 1) {
              // Polynomial spline evaluation
              eval_polynomial(setpoint, startTime, currentTime, index, dimension);
            } else if (trajObj[dimension].type[index] == 2) {
              // Function evaluation
              eval_function(setpoint, currentTime, index, dimension);
            } else if (trajObj[dimension].type[index] == 3) {
              // Bezier curve evaluation
              eval_bezier(setpoint, startTime, currentTime, index, dimension);
            } else {
              configASSERT(false); // A mode other than {0,1,2,3} should not be possible
            }
          }
        }
      }
    }
  }
  // Used for debugging
  for (int ii = 0; ii < 4; ii++ ){
    for (int jj = 0; jj < 5; jj++ ){
      flatOutput[ii][jj] = setpoint->gamma[ii][jj];
    }
  }
}

void eval_point_path(M_setpoint_t *setpoint, int index, int dim)
{
  /****************************************************************************
   * Simulates a fourth order system to find feasible derivative terms when 
   * using the linear path setpoints. Note that each derivative term is
   * simulated with a separate difference equation, but that the characteristic
   * polynomial of the system is invariant across the derivative terms on
   * on account of Tustin's approximation.
   ***************************************************************************/
  float Ycontrib[FLAT_OUTPUT_DIMENSIONS] = {0};
  float Ucontrib[FLAT_OUTPUT_DIMENSIONS] = {0};
  uint32_t p = filterHistory[dim].ptr;
  for (int n = 0; n < FLAT_OUTPUT_DIMENSIONS; n++){
    filterHistory[dim].U[n][p % FLAT_OUTPUT_DIMENSIONS] = trajObj[dim].data[index][0];
    for (int ii = 1; ii < FLAT_OUTPUT_DIMENSIONS; ii++){
      Ycontrib[n] -= LP_denominator[ii] * filterHistory[dim].Y[n][(ii + p) % FLAT_OUTPUT_DIMENSIONS];
    }
    for (int ii = 0; ii < FLAT_OUTPUT_DIMENSIONS; ii++){
      Ucontrib[n] += LP_numerator[n][ii] * filterHistory[dim].U[n][(ii + p) % FLAT_OUTPUT_DIMENSIONS];
    }
    filterHistory[dim].Y[n][p % FLAT_OUTPUT_DIMENSIONS] = Ycontrib[n] + Ucontrib[n];
    setpoint->gamma[dim][n] = Ycontrib[n] + Ucontrib[n];
  }
  filterHistory[dim].ptr--;
}

void eval_function(M_setpoint_t *setpoint, float t, int index, int dimension)
{
  /****************************************************************************
   * Evaluates a sinusoid function and it's derivative terms.
   ***************************************************************************/
  float A = trajObj[dimension].data[index][0];
  float B = trajObj[dimension].data[index][1];
  float w = trajObj[dimension].data[index][2];
  float p = trajObj[dimension].data[index][3];
  float vSin = sinf(w * t + p);
  float vCos = cosf(w * t + p);
  float w2 = w * w;
  float w3 = w2 * w;
  float w4 = w3 * w;
  setpoint->gamma[dimension][0] = A * vSin + B;
  setpoint->gamma[dimension][1] = A * w * vCos;
  setpoint->gamma[dimension][2] = -A * w2 * vSin;
  setpoint->gamma[dimension][3] = -A * w3 * vCos;
  setpoint->gamma[dimension][4] = A * w4 * vSin;
}

void eval_polynomial(M_setpoint_t *setpoint, float startTime, float currentTime, int index, int dimension)
{
  /****************************************************************************
   * Evaluates a polynomial function and it's derivative terms.
   ***************************************************************************/
  float t0 = 1.0f;
  float t1 = currentTime - startTime;
  float t2 = t1 * t1;
  float t3 = t2 * t1;
  float t4 = t3 * t1;
  float t5 = t4 * t1;
  float t[6] = {t0, t1, t2, t3, t4, t5};
  float c[6] = {trajObj[dimension].data[index][0], trajObj[dimension].data[index][1],
                trajObj[dimension].data[index][2], trajObj[dimension].data[index][3],
                trajObj[dimension].data[index][4], trajObj[dimension].data[index][5]};

  int maxorder = 5;
  // Iterates over the flat output derivatives
  for (int order = 0; order < maxorder; order++){
    setpoint->gamma[dimension][order] = 0.0f;
    // Write values
    for (int n = 0; n < (maxorder - order); n++){
      setpoint->gamma[dimension][order] += t[n] * c[n+order];
    }
    // Derive polynomial coefficients
    int count = 0;
    for (int ii = order; ii < maxorder; ii++){
        c[ii] *= count;
        count++;
    }
  }
}

void eval_bezier(M_setpoint_t *setpoint, float startTime, float currentTime, int index, int dimension){
  /****************************************************************************
   * Evaluates a bezier curve and it's derivative terms.
   ***************************************************************************/
  // TODO: Write out equations, currently does nothing.
  float t = (currentTime - startTime) / trajObj[dimension].time[index];
  float tp2 = t * t;
  float tp3 = t * tp2;
  float tm = 1.0f - t;
  float tmp2 = tm * tm;
  float tmp3 = tm * tmp2;
  float p0 = trajObj[dimension].data[index][0];
  float p1 = trajObj[dimension].data[index][1];
  float p2 = trajObj[dimension].data[index][2];
  float p3 = trajObj[dimension].data[index][3];
  setpoint->gamma[dimension][0] = tmp3 * p0 + 3.0f * tmp2 * t * p1 + 3.0f * tm * tp2 * p2 + tp3 * p3;
  setpoint->gamma[dimension][1] = 3.0f * tmp2 * (p1 - p0) + 6.0f * tm * t * ( p2 - p1 )  + 3.0f * tp2 * ( p3 - p2 );
  setpoint->gamma[dimension][2] = 6.0f * tm * (p2 - 2.0f * p1 + p0) + 6.0f * t * (p3 - 2.0f * p2 + p1);
  setpoint->gamma[dimension][3] = 6.0f * ( p3 - p0 ) + 18.0f * ( p1 - p2 );
  setpoint->gamma[dimension][4] = 0.0f;
}

bool decodePacket(CRTPPacket *pk, M_mode_t *mode, M_setpoint_t *setpoint){
  /****************************************************************************
   * Check how to cast data based on the first byte of data and fill the
   * trajectory structure based on the data contained in the CRTPPacket.
   *
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   *
   * packetType=0 - Synchronization (synchronization, clearing and settings of the trajectory)
   * 
   * |packetType|clear|circular0|..|circular3|number0|..|number3|time0|..|time3|
   * |int8      |int8 |int8     |..|int8     |int8   |..|int8   |int16|..|int16|
   *
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   *
   * packetType=1 - Trajectory data (part of a greater pre-loaded trajectory)
   *
   * |packetType|data0|..|data5|time |index|dimension|number|type|
   * |int8      |int16|..|int16|int16|int8 |int8     |int8  |int8|
   *
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   *
   * packetType=2 - Direct data (points in flat output space with derivatives)
   *
   * |packetType|enable|xpos |..|xjerk|ypos |..|yjerk|zpos |..|zjerk|yaw  |yawrate|
   * |int8      |int8  |int16|..|int16|int16|..|int16|int16|..|int16|int16|int16  |
   *
   ***************************************************************************/
  
  if (pk->data[0] == 0){
    // ~~~ Synchronization packet ~~~
    //
    // This packet is used to clear a trajectory object, set basic settings
    // and can also be used to synchronize the trajectory. By default, the
    // tick starts at 0 and time at 0 when initializing the quadcopter. By
    // setting a synchronization tick and specifying a starttime allows for
    // synchronization while executing a trajectory.
    //
    // When clearing the trajectory, the trajectory status is set to false
    // This can only be set to true (allowing setpoint evaluation) by
    // providing a starttime and a synchronization tick, at which point
    // a check is made to see if all trajectory entries have been set. If so,
    // the reference trajectory is followed.

    synchronizationPacket_t *packet = (synchronizationPacket_t*)(pk->data);

    // Fill the trajectory object
    for (int dim = 0; dim < NUM_FLAT_OUTPUTS; dim++){
      if (packet->synchronize == 1){
        // Fill trajectory settings and asserts that the loading is complete
        trajObj[dim].syncronizationTick = xTaskGetTickCount();
        trajObj[dim].status = assertTrajectoryLoaded(dim);
        testCounter[0]++;
      } else {
        // Clear the trajectory object
        trajObj[dim].status = false;
        trajObj[dim].circular = (packet->circular[dim] == 1);
        trajObj[dim].numberOfEntries = packet->number[dim];
        trajObj[dim].startTime = half2single(packet->time[dim]);
        clearTrajectory(dim);
        testCounter[1]++;
      }
    }
    // Clears queue and proceeds without evaluating any trajectory objects
    return true;  

  } else if (pk->data[0] == 1) {
    // ~~~ Trajectory data packet ~~~
    trajectoryPacket_t *packet = (trajectoryPacket_t*)(pk->data);
    
    testCounter[2]++;
    
    uint8_t index = packet->index;
    uint8_t number = packet->number;
    uint8_t dimension = packet->dimension;
    uint8_t type = packet->type;
    float time = half2single(packet->time);
    testVal[0] = time;

    // Sanity check for inputs, loading wonky data here would be very bad..!
    configASSERT(index <= number);              // The data is corrupted if the index > number
    configASSERT(index <= number);              // The data is corrupted if the index > number
    configASSERT(dimension < MAX_DIMENSION);    // A total of four possible dimensions
    configASSERT(type < MAX_TRAJECTORY_TYPE);   // A total of five possible types
    //configASSERT(time > 0.0f);                  // Time must always be greater than zero

    // Halt setpoint generation until all trajectories have been loaded
    trajObj[dimension].status = false;
    
    // Load data
    for (int ii = 0; ii < FLAT_OUTPUT_DIMENSIONS; ii++){
      trajObj[dimension].data[index][ii] = (float)half2single(packet->data[ii]);
    }
    trajObj[dimension].time[index] = time;
    trajObj[dimension].type[index] = type;
    trajObj[dimension].isset[index]  = true;
    
    // Takes the next element in the queue
    return true;

  } else if (pk->data[0] == 2) {
    // ~~~ Direct data packet ~~~
    pointPacket_t *packet = (pointPacket_t*)(pk->data);

    testCounter[3]++;
    
    // x
    setpoint->gamma[0][0] = half2single(packet->x.pos);
    setpoint->gamma[0][1] = half2single(packet->x.vel);
    setpoint->gamma[0][2] = half2single(packet->x.acc);
    setpoint->gamma[0][3] = half2single(packet->x.jerk);
    setpoint->gamma[0][4] = 0.0f;
    
    // y
    setpoint->gamma[1][0] = half2single(packet->y.pos);
    setpoint->gamma[1][1] = half2single(packet->y.vel);
    setpoint->gamma[1][2] = half2single(packet->y.acc);
    setpoint->gamma[1][3] = half2single(packet->y.jerk);
    setpoint->gamma[1][4] = 0.0f;
    
    // z
    setpoint->gamma[2][0] = half2single(packet->z.pos);
    setpoint->gamma[2][1] = half2single(packet->z.vel);
    setpoint->gamma[2][2] = half2single(packet->z.acc);
    setpoint->gamma[2][3] = half2single(packet->z.jerk);
    setpoint->gamma[2][4] = 0.0f;
    
    // yaw
    setpoint->gamma[3][0] = half2single(packet->yaw.pos);
    setpoint->gamma[3][1] = half2single(packet->yaw.vel);
    setpoint->gamma[3][2] = 0.0f;
    setpoint->gamma[3][3] = 0.0f;
    setpoint->gamma[3][4] = 0.0f;
  
    
    // mode
    mode->controller.current = packet->mode;
    mode->flatness.eulerAngles = (bool)(packet->mode);
    mode->flatness.quaternion = (bool)(packet->mode);
    mode->flatness.any = (bool)(packet->mode);
    mode->flatness.bodyacceleration = (bool)(packet->mode);
    mode->flatness.bodyrate = (bool)(packet->mode);
    mode->flatness.torque = (bool)(packet->mode);
    
    // Clears queue and proceeds without evaluating any trajectory objects
    // as the flat output have been specified, the data is used directly in
    // the flatness generator
    for (int dim = 0; dim < NUM_FLAT_OUTPUTS; dim++){
       trajObj[dim].status = false;
    }
    return false;
  } else {
    // Unsupported data packet, type > 2 is not defined.
    configASSERT(false); 
    return true;
  }
}

bool assertTrajectoryLoaded(int dimension){
  /****************************************************************************
   * Iterates over the elements in a trajectory along a flat dimension and
   * checks that all elements at indices up to *.numberOfEntries have been set
   ***************************************************************************/
  for (int index = 0; index < trajObj[dimension].numberOfEntries; index++){
    if (trajObj[dimension].isset[index] == false){
      return false;
    }
  }
  return true;
}

void clearTrajectory(int dimension){
  /****************************************************************************
   * Sets the isset flag to false for all elementss in all trajectories
   ***************************************************************************/
  for (int index = 0; index < MAX_TRAJECTORY_ENTRIES; index++){
    trajObj[dimension].isset[index] = false;
  }
}

LOG_GROUP_START(testA)
LOG_ADD(LOG_FLOAT, a, &flatOutput[0][0])
LOG_ADD(LOG_FLOAT, b, &flatOutput[1][0])
LOG_ADD(LOG_FLOAT, c, &flatOutput[2][0])
LOG_ADD(LOG_FLOAT, d, &flatOutput[0][1])
LOG_ADD(LOG_FLOAT, e, &flatOutput[1][1])
LOG_ADD(LOG_FLOAT, f, &flatOutput[2][1])
LOG_GROUP_STOP(testA)

LOG_GROUP_START(testB)
LOG_ADD(LOG_FLOAT, a, &testVal[0])
LOG_ADD(LOG_UINT32, b, &testCounter[0])
LOG_ADD(LOG_UINT32, c, &testCounter[1])
LOG_ADD(LOG_UINT32, d, &testCounter[2])
LOG_ADD(LOG_UINT32, e, &testCounter[3])
LOG_ADD(LOG_UINT32, f, &testCounter[4])
LOG_GROUP_STOP(testB)
