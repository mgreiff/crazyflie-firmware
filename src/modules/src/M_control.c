#include "FreeRTOS.h"
#include "param.h"
#include "log.h"
#include "math.h"

#include "M_control.h"
#include "arm_math.h"

// Loop rates
#define PD_ATTITUDE_RATE RATE_500_HZ
#define PD_ELEVATION_RATE RATE_500_HZ
#define PID_TRANSLATION_RATE RATE_500_HZ
#define LQR_HOVER_RATE RATE_500_HZ
#define SE3_ATTIDUDE_RATE RATE_500_HZ
#define SE3_TRANSLATION_RATE RATE_500_HZ

// LQR dimensions
#define LQR_N_CSIG  4
#define LQR_8_N_STATE 8
#define LQR_12_N_STATE 12

// Parameter for logging
control_error_t controlError;

// PID controller parameters for logging
static float attitude_KP_PR;    // Pitch and roll
static float attitude_KD_PR;    // Pitch and roll
static float attitude_KP_Y;     // Yaw
static float attitude_KD_Y;     // Yaw
static float elevation_KP;      // Z
static float elevation_KD;      // Z
static float translation_KP;    // XY
static float translation_KI;    // XY
static float translation_KD;    // XY
static float translation_N;     // XY
static float translation_alpha; // XY
static float translation_beta;  // XY

static float errorA;
static float errorB;
static float errorC;
static float errorD;
static float errorE;
static float errorF;

// SE(3) controller tuning parameters
//static float SE3_KP = 0.2f;       // X, Y ,Z
//static float SE3_KD = 0.08f;      // dX, dY, dZ
//static float SE3_KRX = 0.006f;     // R
//static float SE3_KWX = 0.000025f;  // Omega_B
//static float SE3_KRZ = 0.006f;     // R
//static float SE3_KWZ = 0.000025f;  // Omega_B

// SE(3) controller tuning parameters
static float SE3_KP = 0.2f;       // X, Y ,Z
static float SE3_KD = 0.08f;      // dX, dY, dZ
static float SE3_KRX = 0.0058f;     // R
static float SE3_KWX = 0.000020f;  // Omega_B
static float SE3_KRZ = 0.0058f;     // R
static float SE3_KWZ = 0.000020f;  // Omega_B

//static float SE3_KP = 0.2f;       // X, Y ,Z
//static float SE3_KD = 0.08f;      // dX, dY, dZ
//static float SE3_KRX = 0.0055f;     // R
//static float SE3_KWX = 0.000017f;  // Omega_B
//static float SE3_KRZ = 0.006f;     // R
//static float SE3_KWZ = 0.000025f;  // Omega_B

// Help function for ARM matrix multiplication
static inline void mat_mult(const arm_matrix_instance_f32 * pSrcA,
                            const arm_matrix_instance_f32 * pSrcB,
                            arm_matrix_instance_f32 * pDst)
{ configASSERT(ARM_MATH_SUCCESS == arm_mat_mult_f32(pSrcA, pSrcB, pDst)); }



// TILQR - LQR control signals
static float U[LQR_N_CSIG];
static arm_matrix_instance_f32 Um = {LQR_N_CSIG, 1, (float *)U};
// TILQR - 12 state gain matrix [z,y,z,dx,dy,dz, phi,theta,psi,dphi,dtheta,dpsi]
static float K_12[4][12];
static arm_matrix_instance_f32 Km_12 = {LQR_N_CSIG, LQR_12_N_STATE, (float *)K_12};
// TILQR - 8 state gain matrix [z, dz, phi,theta,psi,dphi,dtheta,dpsi]
static float K_8[4][8];
static arm_matrix_instance_f32 Km_8 = {LQR_N_CSIG, LQR_8_N_STATE, (float *)K_8};

void stateControllerInit_M(){
  attitude_KP_PR = 13.0f;   // Pitch and roll
  attitude_KD_PR = 7.0f;    // Pitch and roll
  attitude_KP_Y = 8.0f;     // Yaw
  attitude_KD_Y = 4.0f;     // Yaw
  
  elevation_KP = 10.0f;     // Z
  elevation_KD = 5.0f;      // Z
  
  translation_KP = 0.2f;    // XY
  translation_KI = 0.0f;    // XY
  translation_KD = 0.24f;   // XY
  translation_N = 100.0f;   // XY
  translation_alpha = 1.0f; // XY
  translation_beta = 1.0f;  // XY

  K_8[0][0] = 3.162277660168373f;
  K_8[0][1] = 1.082008019142515f;
  K_8[1][2] = 0.003214342856635f;
  K_8[1][5] = 1.082008019142515f;
  K_8[2][3] = 0.003214342856635f;
  K_8[2][6] = 1.082008019142515f;
  K_8[3][4] = 0.01f;
  K_8[3][7] = 0.003253613375925f;
}

void stateController_M(M_control_t *controlsignal,
                       M_setpoint_t *setpoint,
                       M_mode_t *mode,
                       M_state_t *state,
                       uint32_t tick)
{
  /****************************************************************************
  * High level method for updating the control signal object depending on the
  * specified controller mode. Simply calls the appropriate combination of
  * lower level controllers. Additional controllers can be implemented and
  * integrated by simply adding a new mode.
  ****************************************************************************/
  if (mode->controller.current == 0 || mode->controller.current == 0xFF){
    // IDLE state - do nothing and set references to zero
    controlsignal->thrust = 0;
    controlsignal->torque.x = 0;
    controlsignal->torque.y = 0;
    controlsignal->torque.z = 0;
  } else if (mode->controller.current == 1) {
    //
    // ~~~ PID hover with PID translation control ~~~ 
    //
    // Map error in xy to roll/pitch
    PID_translation_control(setpoint, state, tick);
    // Compute control signal thrust
    PD_elevation_control(controlsignal, setpoint, state, tick);
    // Compute control signal torque
    PD_attitude_control(controlsignal, setpoint, state, tick);  
  } else if (mode->controller.current == 2) {
    //
    // ~~~ LQR hover with PID translation control ~~~
    //
    // Map error in xy to roll/pitch 
    PID_translation_control(setpoint, state, tick); 
    LQR_hover_control(controlsignal, setpoint, state, tick); 
  } else if (mode->controller.current == 3) {
    //
    // ~~~ LQR full control ~~~ 
    //
    // Compute control signals (both thrust and torque)
    LQR_full_control(controlsignal, setpoint, state, tick);
  } else if (mode->controller.current == 4) {
    //
    // ~~~ SE(3) full control ~~~ 
    //
    // Compute control signals (both thrust and torque)
    SE3_full_control(controlsignal, setpoint, state, tick); 
  }
}

void PD_elevation_control(M_control_t *controlsignal,
                          M_setpoint_t *setpoint,
                          M_state_t *state,
                          uint32_t tick)
{ 
  if (RATE_DO_EXECUTE(PD_ATTITUDE_RATE, tick)) {
    controlError.position.z = setpoint->position.z - state->position.z;
    controlError.velocity.z = setpoint->velocity.z - state->velocity.z;

    controlsignal->thrust = (elevation_KP * controlError.position.z +
                             elevation_KD * controlError.velocity.z +
                             Q_MASS * Q_GRAVITY) / (setpoint->rotation[2][2]);
  }
}

void PD_attitude_control(M_control_t *controlsignal,
                         M_setpoint_t *setpoint,
                         M_state_t *state,
                         uint32_t tick)
{
  if (RATE_DO_EXECUTE(PD_ATTITUDE_RATE, tick)) {
    controlError.eulerAngles.roll = setpoint->eulerAngles.roll - state->eulerAngles.roll;
    controlError.eulerAngles.pitch = setpoint->eulerAngles.pitch - state->eulerAngles.pitch;
    controlError.eulerAngles.yaw = setpoint->eulerAngles.yaw - state->eulerAngles.yaw;

    controlError.eulerRates.roll = setpoint->eulerRates.roll - state->eulerRates.roll;
    controlError.eulerRates.pitch = setpoint->eulerRates.pitch - state->eulerRates.pitch;
    controlError.eulerRates.yaw = setpoint->eulerRates.yaw - state->eulerRates.yaw;
  
    controlsignal->torque.x = attitude_KP_PR * controlError.eulerAngles.roll +
                              attitude_KD_PR * controlError.eulerRates.roll;
    controlsignal->torque.y = attitude_KP_PR * controlError.eulerAngles.pitch +
                              attitude_KD_PR * controlError.eulerRates.pitch;
    controlsignal->torque.z = attitude_KP_Y * controlError.eulerAngles.yaw +
                              attitude_KP_Y * controlError.eulerRates.yaw;
  }
}

void PID_translation_control(M_setpoint_t *setpoint,
                             M_state_t *state,
                             uint32_t tick)
{
  if (RATE_DO_EXECUTE(PID_TRANSLATION_RATE, tick)) {
    // Control errors in the global frame
    float e_GX = setpoint->position.x - state->position.x;
    float e_GY = setpoint->position.y - state->position.y;
    float e_GZ = setpoint->position.z - state->position.z;

    float de_GX = setpoint->velocity.x - state->velocity.x;
    float de_GY = setpoint->velocity.y - state->velocity.y;
    float de_GZ = setpoint->velocity.z - state->velocity.z;

    //e_B = R^T*e_G to get control errors in the body frame
    float e_BX = state->rotation[0][0] * e_GX + state->rotation[1][0] * e_GY + state->rotation[2][0] * e_GZ;
    float e_BY = state->rotation[0][1] * e_GX + state->rotation[1][1] * e_GY + state->rotation[2][1] * e_GZ;

    float de_BX = state->rotation[0][0] * de_GX + state->rotation[1][0] * de_GY + state->rotation[2][0] * de_GZ;
    float de_BY = state->rotation[0][1] * de_GX + state->rotation[1][1] * de_GY + state->rotation[2][1] * de_GZ;
    
    // Map x-error to negative theta (roll)
    setpoint->eulerAngles.roll = -(translation_KP * e_BX + translation_KD * de_BX);
    
    // Map y-error to positive psi (pitch)
    setpoint->eulerAngles.pitch =  translation_KP * e_BY + translation_KD * de_BY;
  }
}

void LQR_hover_control(M_control_t *controlsignal,
                       M_setpoint_t *setpoint,
                       M_state_t *state,
                       uint32_t tick)
{
  if (RATE_DO_EXECUTE(LQR_HOVER_RATE, tick)) {
    // Compute control errors and store for logging
    controlError.position.z = setpoint->position.x - state->position.z;
    controlError.velocity.z = setpoint->velocity.x - state->velocity.z;
    controlError.eulerAngles.roll = setpoint->eulerAngles.roll - state->eulerAngles.roll;
    controlError.eulerAngles.pitch = setpoint->eulerAngles.pitch - state->eulerAngles.pitch;
    controlError.eulerAngles.yaw = setpoint->eulerAngles.yaw - state->eulerAngles.yaw;
    controlError.eulerRates.roll = setpoint->eulerRates.roll - state->eulerRates.roll;
    controlError.eulerRates.pitch = setpoint->eulerRates.pitch - state->eulerRates.pitch;
    controlError.eulerRates.yaw = setpoint->eulerRates.yaw - state->eulerRates.yaw;
    
    // Fill error vector
    float X [8] = {
      controlError.position.z, controlError.velocity.z,
      controlError.eulerAngles.roll, controlError.eulerAngles.pitch, controlError.eulerAngles.yaw,
      controlError.eulerRates.roll, controlError.eulerRates.pitch, controlError.eulerRates.yaw
    };
    arm_matrix_instance_f32 Xm = {LQR_8_N_STATE, 1, (float *)X};

    // Compute control signals from the LQR gain
    mat_mult(&Km_8, &Xm, &Um);
    
    //  note the linearization about (T, tau_x, tau_x, tau_z) = (mg, 0, 0, 0)
    controlsignal->thrust = Q_MASS * Q_GRAVITY + U[0];
    controlsignal->torque.x = U[1];
    controlsignal->torque.y = U[2];
    controlsignal->torque.z = U[3];
  }
}

void LQR_full_control(M_control_t *controlsignal,
                      M_setpoint_t *setpoint,
                      M_state_t *state,
                      uint32_t tick)
{
  if (RATE_DO_EXECUTE(LQR_HOVER_RATE, tick)) {
    // Compute control errors and store for logging
    controlError.position.x = setpoint->position.x - state->position.x;
    controlError.position.y = setpoint->position.x - state->position.y;
    controlError.position.z = setpoint->position.x - state->position.z;
    controlError.velocity.x = setpoint->velocity.x - state->velocity.x;
    controlError.velocity.y = setpoint->velocity.x - state->velocity.y;
    controlError.velocity.z = setpoint->velocity.x - state->velocity.z;
    controlError.eulerAngles.roll = setpoint->eulerAngles.roll - state->eulerAngles.roll;
    controlError.eulerAngles.pitch = setpoint->eulerAngles.pitch - state->eulerAngles.pitch;
    controlError.eulerAngles.yaw = setpoint->eulerAngles.yaw - state->eulerAngles.yaw;
    controlError.eulerRates.roll = setpoint->eulerRates.roll - state->eulerRates.roll;
    controlError.eulerRates.pitch = setpoint->eulerRates.pitch - state->eulerRates.pitch;
    controlError.eulerRates.yaw = setpoint->eulerRates.yaw - state->eulerRates.yaw;
    
    // Fill error vector
    float X [12] = {
      controlError.position.x, controlError.position.y, controlError.position.z,
      controlError.velocity.x, controlError.velocity.y, controlError.velocity.z,
      controlError.eulerAngles.roll, controlError.eulerAngles.pitch, controlError.eulerAngles.yaw,
      controlError.eulerRates.roll, controlError.eulerRates.pitch, controlError.eulerRates.yaw
    };
    arm_matrix_instance_f32 Xm = {LQR_12_N_STATE, 1, (float *)X};

    // Compute control signals from the LQR gain
    mat_mult(&Km_12, &Xm, &Um);
    
    //  note the linearization about (T, tau_x, tau_x, tau_z) = (mg, 0, 0, 0)
    controlsignal->thrust = Q_MASS * Q_GRAVITY + U[0];
    controlsignal->torque.x = U[1];
    controlsignal->torque.y = U[2];
    controlsignal->torque.z = U[3];
  }
}

void SE3_full_control(M_control_t *controlsignal,
                      M_setpoint_t *setpoint,
                      M_state_t *state,
                      uint32_t tick)
{
  if (RATE_DO_EXECUTE(SE3_TRANSLATION_RATE, tick)) {
    // ~~~ Compute control errors in translational position and velocity ~~~
    controlError.position.x = setpoint->position.x - state->position.x;
    controlError.position.y = setpoint->position.y - state->position.y;
    controlError.position.z = setpoint->position.z - state->position.z;

    controlError.velocity.x = setpoint->velocity.x - state->velocity.x;
    controlError.velocity.y = setpoint->velocity.y - state->velocity.y;
    controlError.velocity.z = setpoint->velocity.z - state->velocity.z;

    // ~~~ Desired forces (F_d) in the global frame ~~~
    float Fd[3];
    Fd[0] = SE3_KP * controlError.position.x +
            SE3_KD * controlError.velocity.x +
            Q_MASS * setpoint->acceleration.x;
    Fd[1] = SE3_KP * controlError.position.y +
            SE3_KD * controlError.velocity.y + 
            Q_MASS * setpoint->acceleration.y;
    Fd[2] = SE3_KP * controlError.position.z +
            SE3_KD * controlError.velocity.z +
            Q_MASS * setpoint->acceleration.z +
            Q_MASS * Q_GRAVITY;

    errorE = SE3_KP * controlError.position.z;
    errorF = SE3_KD * controlError.velocity.z;
    // ~~~ Normalize F_d to find the desired unit basis vectors ~~~
    // z-body axis
    float zBd [3];
    normalizeVector(Fd, zBd);
    
    // x-body axis
    float xBd [3];
    float yCd [3];
    float yCdcrosszBd [3];
    yCd[0] = cosf(setpoint->eulerAngles.yaw + (float)M_PI/2.0f);
    yCd[1] = sinf(setpoint->eulerAngles.yaw + (float)M_PI/2.0f);
    yCd[2] = 0.0;

    crossProd(yCd, zBd, yCdcrosszBd);
    normalizeVector(yCdcrosszBd, xBd);

    // y-body axis
    float yBd [3];
    crossProd(zBd, xBd, yBd);

    // ~~~ Compute control errors in rotation and body rate ~~~
    // Here, we know that R_d = [xB_d, yB_d, zB_d] and compute eR = (R_d'*R-R'*R_d)^v
    controlError.rotation.x = -0.5f * (  state->rotation[0][1]*zBd[0]
                                     + state->rotation[1][1]*zBd[1]
                                     + state->rotation[2][1]*zBd[2]
                                     - state->rotation[0][2]*yBd[0]
                                     - state->rotation[1][2]*yBd[1]
                                     - state->rotation[2][2]*yBd[2]);

    controlError.rotation.y = -0.5f * (  state->rotation[0][2]*xBd[0]
                                     + state->rotation[1][2]*xBd[1]
                                     + state->rotation[2][2]*xBd[2]
                                     - state->rotation[0][0]*zBd[0]
                                     - state->rotation[1][0]*zBd[1]
                                     - state->rotation[2][0]*zBd[2]);

    controlError.rotation.z = -0.5f * (  state->rotation[0][0]*yBd[0]
                                     + state->rotation[1][0]*yBd[1]
                                     + state->rotation[2][0]*yBd[2]
                                     - state->rotation[0][1]*xBd[0]
                                     - state->rotation[1][1]*xBd[1]
                                     - state->rotation[2][1]*xBd[2]);
    
    controlError.bodyrate.x = setpoint->bodyrate.x - state->bodyrate.x;
    controlError.bodyrate.y = setpoint->bodyrate.y - state->bodyrate.y;
    controlError.bodyrate.z = setpoint->bodyrate.z - state->bodyrate.z;

    // Compute the body frame z-axis, corresponding to z_B = R*z_G = R[:][2], where T = F \cdot z_B
    controlsignal->thrust = Fd[0] * state->rotation[0][2] +
                            Fd[1] * state->rotation[1][2] +
                            Fd[2] * state->rotation[2][2];       // [N]
    
    errorA = SE3_KRX * controlError.rotation.x;
    errorB = SE3_KWX * controlError.bodyrate.x;
    errorC = SE3_KRZ * controlError.rotation.z;
    errorD = SE3_KWZ * controlError.bodyrate.z;

    // Compute the body torques
    controlsignal->torque.x = SE3_KRX * controlError.rotation.x +
                              SE3_KWX * controlError.bodyrate.x;   // [Nm]
    controlsignal->torque.y = SE3_KRX * controlError.rotation.y +
                              SE3_KWX * controlError.bodyrate.y;   // [Nm]
    controlsignal->torque.z = SE3_KRZ * controlError.rotation.z +
                              SE3_KWZ * controlError.bodyrate.z;   // [Nm]
  }
}

void crossProd(float u[], float v[], float res[])
{
  /****************************************************************************
  * Computes the resulting cross product of u and v in res.
  *
  * ARGS:
  *   u - Vector in R^3 (float [3]).
  *   v - Vector in R^3 (float [3]).
  *   res - Vector in R^3 (float [3]). 
  ****************************************************************************/
  res[0] =   ( (u[1] * v[2]) - (u[2] * v[1]) );
  res[1] = - ( (u[0] * v[2]) - (u[2] * v[0]) );
  res[2] =   ( (u[0] * v[1]) - (u[1] * v[0]) );
}

void normalizeVector(float u[], float res[])
{
  /****************************************************************************
  * Computes the resulting cross product of u and v in res.
  *
  * ARGS:
  *   u - Vector in R^3 (float [3]).
  *   v - Vector in R^3 (float [3]).
  *   res - Vector in R^3 (float [3]). 
  ****************************************************************************/
  float u_two_norm = sqrtf(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
  res[0] =  u[0] / u_two_norm;
  res[1] =  u[1] / u_two_norm;
  res[2] =  u[2] / u_two_norm;
}

PARAM_GROUP_START(SE3)
PARAM_ADD(PARAM_FLOAT, KP, &SE3_KP)
PARAM_ADD(PARAM_FLOAT, KD, &SE3_KD)
PARAM_ADD(PARAM_FLOAT, KRX, &SE3_KRX)
PARAM_ADD(PARAM_FLOAT, KWX, &SE3_KWX)
PARAM_ADD(PARAM_FLOAT, KRZ, &SE3_KRZ)
PARAM_ADD(PARAM_FLOAT, KWZ, &SE3_KWZ)
PARAM_GROUP_STOP(SE3)

/*
PARAM_GROUP_START(tuning)
PARAM_ADD(PARAM_FLOAT, p_A, &tuningA)
PARAM_ADD(PARAM_FLOAT, p_B, &tuningB)
PARAM_ADD(PARAM_FLOAT, p_C, &tuningC)
PARAM_ADD(PARAM_FLOAT, p_D, &tuningD)
PARAM_GROUP_STOP(tuning)
*/

LOG_GROUP_START(errors)
LOG_ADD(LOG_FLOAT, x, &errorA)
LOG_ADD(LOG_FLOAT, y, &errorB)
LOG_ADD(LOG_FLOAT, z, &errorC)
LOG_ADD(LOG_FLOAT, dx, &errorD)
LOG_ADD(LOG_FLOAT, dy, &errorE)
LOG_ADD(LOG_FLOAT, dz, &errorF)
LOG_GROUP_STOP(errors)

