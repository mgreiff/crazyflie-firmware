#include "M_powerDistribution.h"
#include "param.h"
#include "math.h"
#include "motors.h"
#include "log.h"

#define POWER_X_CONFIGURATION

// SATURATION_EPSILON set low for speed, high for increased responsiveness
static float SATURATION_EPSILON = 0.2f;   //\in [0,0.5] (unitless)
static float LIM_MAX_OMEGA = 3500.0f;     // [rad/s]
static float LIM_MIN_OMEGA = 1300.0f;        // [rad/s]
static float Q_COEFF_K = 0.0000000128f;   // [N * m * s^2]
static float Q_COEFF_B = 0.0000000000765f;   // [N * s^2]
static float Q_COEFF_BETA_1 = 0.0877f;   // [.]
static float Q_COEFF_BETA_2 = 0.0660f;   // [.]

// Temporary variables for constants to improve efficiency
static float d4K;
static float d2KL;
static float d4B;
static float B1d2B2;
static float B1d2B2sq;
static float KdB2;

static float thrust_max_lim;
static float thrust_min_lim;
static float omegasq_max_lim;
static float omegasq_min_lim;

// Temporary variables for debugging
// TODO: remove
static float tempFloat[4];
static float tempInt[2];

void powerDistributionInit_M()
{
  /****************************************************************************
   * Initializes motor and pre-computes parameters used in the rotor speed maps
   ***************************************************************************/

  // Temporary parameters
  d4K = 1.0f / (4.0f * Q_COEFF_K);
  d2KL = 1.0f / (2.0f * Q_COEFF_K * Q_LENGTH);
  d4B = 1.0f / (4.0f * Q_COEFF_B);
  B1d2B2 = Q_COEFF_BETA_1 / (2.0f * Q_COEFF_BETA_2);
  B1d2B2sq = B1d2B2 * B1d2B2;
  KdB2 = Q_COEFF_K / Q_COEFF_BETA_2;
  
  // Saturation bounds
  omegasq_max_lim = LIM_MAX_OMEGA * LIM_MAX_OMEGA;
  omegasq_min_lim = LIM_MIN_OMEGA * LIM_MIN_OMEGA;
  thrust_max_lim = (1.0f - SATURATION_EPSILON) * Q_COEFF_K * 4.0f * omegasq_max_lim;
  thrust_min_lim = SATURATION_EPSILON * Q_COEFF_K * 4.0f * omegasq_min_lim;
  
  // Initializes motors
  motorsInit(motorMapDefaultBrushed);
}

void powerDistribution_M(M_mode_t *mode, M_control_t *controlsignal)
{
  /****************************************************************************
   * This fucntion saturates the control signals and writes the PWM duty cycle
   * To the motors based on thrust [N] and torque [Nm] input signals.
   *
   * The purpose of saturating the control signals in this way is to always
   * have the system controllable (in the case where we are saturated in
   * both thrust and rotor speeds simultaneously, we might not be able to
   * generate differing thrust between rotors, in which case we have no way
   * of affecting the attitude --> crash). With this sceme, we should be able
   * to follow loops in the positional state space.
   *
   * The power distribution can be run in two configurations.
   *
   *     “x”-configuration     forward        “+”-configuration
   *                              ^
   *    m4,w4 0    0 m1,w1                        m1,w1   
   *           \  /                                  0
   *            \/              x ^                  |
   *            /\                |       m4,w4  0—-—+-——0 m2,w2
   *           /  \          y    |                  |
   *    m3,w3 0    0 m2,w2   <-——-+ z (up)           0 m3,w3
   *
   * The body coordinate system (right handed)is always given with the x-axis
   * from the centre of  mass (origin) to the motor 1 (m1), and the y axis is
   * always  from origin to motor 4 (m4).
   *
   * In both configurations, thrust is given by the map
   *
   * T = k(w1^2 + w2^2 + w3^2 + w4^2)
   *
   * and torque about the z-axis (corresponding to yaw) is given by
   *
   * tau_z = b(w1^2 + w2^2 + w3^2 + w4^2)
   * 
   * where k and b are (very small) constants as given by blade theory.
   * The difference is in the torque induced about the x and y axes.
   *
   * ~~~ “+”-configuration ~~~
   *
   * tau_x = k*l*(-w2^2 + w4^2)
   * tau_y = k*l*(-w1^2 + w3^2)
   *
   * The map from thrusts and torques [T,tau_x, tau_z, tau_z]^T to rotor
   * speeds is then
   *
   *     | 1/(4*k),          0, -1/(2*k*l),  1/(4*b)|
   * M = | 1/(4*k), -1/(2*k*l),          0, -1/(4*b)|
   *     | 1/(4*k),          0,  1/(2*k*l),  1/(4*b)|
   *     | 1/(4*k),  1/(2*k*l),          0, -1/(4*b)|
   *
   * ~~~ “x”-configuration ~~~
   *
   * tau_x = (k*l/sqrt(2)) * (-w1^2 - w2^2 + w3^2 + w4^2)
   * tau_y = (k*l/sqrt(2)) * (-w1^2 + w2^2 + w3^2 - w4^2)
   *
   * The map from thrusts and torques [T,tau_x, tau_z, tau_z]^T to rotor
   * speeds is then
   *
   *     | 1/(4*k), -2^(1/2)/(4*k*l), -2^(1/2)/(4*k*l),  1/(4*b)|
   * M = | 1/(4*k), -2^(1/2)/(4*k*l),  2^(1/2)/(4*k*l), -1/(4*b)|
   *     | 1/(4*k),  2^(1/2)/(4*k*l),  2^(1/2)/(4*k*l),  1/(4*b)|
   *     | 1/(4*k),  2^(1/2)/(4*k*l), -2^(1/2)/(4*k*l), -1/(4*b)|
   *
   ***************************************************************************/

  // Saturate control signals and compute PWM
  saturate_control(mode, controlsignal);
  
  // Write to motors
  motorsSetRatio(MOTOR_M1, controlsignal->Power[0]);
  motorsSetRatio(MOTOR_M2, controlsignal->Power[1]);
  motorsSetRatio(MOTOR_M3, controlsignal->Power[2]);
  motorsSetRatio(MOTOR_M4, controlsignal->Power[3]);
}

void saturate_control(M_mode_t *mode, M_control_t *controlsignal)
{
  if (mode->controller.current == 0) { // idle
    controlsignal->Power[0] = 0;
    controlsignal->Power[1] = 0;
    controlsignal->Power[2] = 0;
    controlsignal->Power[3] = 0;
    return;
  }
  
  if (mode->controller.current == 0xFF) { // spin motors
    controlsignal->Power[0] = 25000; //arbitrary
    controlsignal->Power[1] = 25000;
    controlsignal->Power[2] = 25000;
    controlsignal->Power[3] = 25000;
    return;
  }
  
  // Saturates thrust
  if (controlsignal->thrust > thrust_max_lim){
    controlsignal->thrust = thrust_max_lim;                      // [N]
  } else if (controlsignal->thrust < thrust_min_lim){
    controlsignal->thrust = thrust_min_lim;                      // [N]
  }

  // Map thrusts and torques to rotor speeds squared
  float Tpart = d4K * controlsignal->thrust;
  float tauZpart = d4B * controlsignal->torque.z;
#if defined(POWER_PLUS_CONFIGURATION)
  float tauXpart = d2KL * controlsignal->torque.x;
  float tauYpart = d2KL * controlsignal->torque.y;
  controlsignal->omegasq[0] = fabs(Tpart - tauYpart + tauZpart); // [rad^2/s^2]
  controlsignal->omegasq[1] = fabs(Tpart - tauXpart - tauZpart); // [rad^2/s^2]
  controlsignal->omegasq[2] = fabs(Tpart + tauYpart + tauZpart); // [rad^2/s^2]
  controlsignal->omegasq[3] = fabs(Tpart + tauXpart - tauZpart); // [rad^2/s^2]
#else
  float tauXpart = (d2KL/sqrtf(2.0f)) * controlsignal->torque.x;
  float tauYpart = (d2KL/sqrtf(2.0f)) * controlsignal->torque.y;
  controlsignal->omegasq[0] = fabs(Tpart - tauXpart - tauYpart - tauZpart); // [rad^2/s^2]
  controlsignal->omegasq[1] = fabs(Tpart - tauXpart + tauYpart + tauZpart); // [rad^2/s^2]
  controlsignal->omegasq[2] = fabs(Tpart + tauXpart + tauYpart - tauZpart); // [rad^2/s^2]
  controlsignal->omegasq[3] = fabs(Tpart + tauXpart - tauYpart + tauZpart); // [rad^2/s^2]
#endif

  // Log variables
  tempFloat[0] = sqrtf(controlsignal->omegasq[0]);
  tempFloat[1] = sqrtf(controlsignal->omegasq[1]);

  for (int ii = 0; ii < 4; ii++){
    // Saturate rotor speed squared
    if (controlsignal->omegasq[ii] > omegasq_max_lim){
      controlsignal->omegasq[ii] = omegasq_max_lim;
    } else if (controlsignal->omegasq[ii] < omegasq_min_lim){
      controlsignal->omegasq[ii] = omegasq_min_lim;
    }
    
    // Find and saturate PWM duty cycle (this value should never be able to become negative if
    // the sarutation bounds are set correctly but is included as a safety measure)
    // TODO: something is wrong in this specific line
    controlsignal->PWM[ii] = -B1d2B2 + sqrt(B1d2B2sq + controlsignal->omegasq[ii] * KdB2); // \in[0,1] (float)
    if (controlsignal->PWM[ii] > 1.0f){
      controlsignal->PWM[ii] = 1.0f;
    } else if (controlsignal->PWM[ii] < 0.0f){
      controlsignal->PWM[ii] = 0.0f;
    }

    // Convert to integerpower setting
    controlsignal->Power[ii] = (uint16_t)(controlsignal->PWM[ii] * 65535.0f); // \in[0,65535] (uint16_t)
  }
  // Log variables
  tempFloat[2] = controlsignal->PWM[0];
  tempFloat[3] = controlsignal->PWM[0];
  tempInt[0] = controlsignal->Power[0];
  tempInt[1] = controlsignal->Power[1];
}

PARAM_GROUP_START(pDist)
PARAM_ADD(PARAM_FLOAT, SAT_EPSILON, &SATURATION_EPSILON)
PARAM_ADD(PARAM_FLOAT, MAX_OMEGA, &LIM_MAX_OMEGA)
PARAM_ADD(PARAM_FLOAT, MIN_OMEGA, &LIM_MIN_OMEGA)
PARAM_ADD(PARAM_FLOAT, K, &Q_COEFF_K)
PARAM_ADD(PARAM_FLOAT, B, &Q_COEFF_B)
PARAM_ADD(PARAM_FLOAT, BETA_1, &Q_COEFF_BETA_1)
PARAM_ADD(PARAM_FLOAT, BETA_2, &Q_COEFF_BETA_2)
PARAM_GROUP_STOP(pDist)


LOG_GROUP_START(powerTest)
LOG_ADD(LOG_FLOAT, aa, &tempFloat[0])
LOG_ADD(LOG_FLOAT, bb, &tempFloat[1])
LOG_ADD(LOG_FLOAT, cc, &tempFloat[2])
LOG_ADD(LOG_FLOAT, dd, &tempFloat[3])
LOG_ADD(LOG_UINT16, ee, &tempInt[0])
LOG_ADD(LOG_UINT16, ff, &tempInt[1])
LOG_GROUP_STOP(powerTest)
