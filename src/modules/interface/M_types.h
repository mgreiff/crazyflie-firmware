/******************************************************************************
 * M_types.h Contains definitions of structs used in the M_* control system.
 * The structure was written so as not to interfere with preexisting datatypes,
 * to allow easy integration with the stock system once complete. In general,
 * the structures haev been expanded to contain more information, and
 * the parameters of the quadcopter model have been provided for LQR/Flatness.
 *****************************************************************************/

// Avoid redefinitions
#ifndef __M_TYPES_H__
#define __M_TYPES_H__

// Included to get the main loop rate and control_t definitions
#include <stabilizer_types.h> 
#include <stdbool.h>

// Maximum sizes of the trjectories
#define NUM_FLAT_OUTPUTS 4
#define FLAT_OUTPUT_DIMENSIONS 5
#define MAX_TRAJECTORY_ENTRIES 20
#define MAX_DIMENSION 6

// Mathematical constants
#define PI 3.14159265359

// Quadcopter model parameters
#define Q_INERTIA_XX 0.0000166f  // [kg * m^2]
#define Q_INERTIA_YY 0.0000166f  // [kg * m^2]
#define Q_INERTIA_ZZ 0.00002930f // [kg * m^2]
#define Q_LENGTH 0.046f          // [m]
#define Q_MASS 0.03f            // [kg]
#define Q_GRAVITY 9.81f          // [m/s^2]

typedef struct M_quaternion_s
{
  float qw;     // Real part
  float qv [3]; // Imaginary part
} M_quaternion_t;

/******************************************************************************
 * The setpoint for the controller, contains all the flat outputs and
 * system states computed in the flatness equations. Many of which will
 * not be used in a final implementation, but are included for debugging
 * purposes.
 *****************************************************************************/
typedef struct M_setpoint_s
{
  float gamma[NUM_FLAT_OUTPUTS][FLAT_OUTPUT_DIMENSIONS]; // Flat outputs
  point_t position;            // Position in the global frame        [m]
  point_t velocity;            // Velocity in the global frame        [m/s]
  point_t acceleration;        // Acceleration in the global frame    [m/s]
  attitude_t eulerAngles;      // Euler angles                        [rad]
  attitude_t eulerRates;       // Euler angle rates                   [rad/s]
  M_quaternion_t quaternion;   // Quaternion attitude                 [rad]
  float rotation [3][3];       // SO(3) rotation                      [.]
  point_t bodyrate;            // Angular rates in the body frame     [rad/s]
  point_t bodyacceleration;    // Angular acc in the body frame       [rad/s^2]
  float thrust;                // Thrust along z_B (body frame)       [N]
  point_t torque;              // Torques                             [Nm]
} M_setpoint_t;

/******************************************************************************
 * The setpoint for the controller, contains all the flat outputs and
 * system states computed in the flatness equations. Many of which will
 * not be used in a final implementation, but are included for debugging
 * purposes.
 *****************************************************************************/
typedef struct M_state_s
{
  point_t position;            // Position                            [m]
  point_t velocity;            // Velocity                            [m/s]
  float rotation [3][3];       // SO(3) rotation                      [.]
  M_quaternion_t quaternion;   // Quaternion attitude                 [rad]
  attitude_t eulerAngles;      // Euler angles                        [rad]
  attitude_t eulerRates;       // Euler angle rates                   [rad/s]
  point_t bodyrate;            // Angular rates in the body frame     [rad/s]
} M_state_t;

/******************************************************************************
 * Trajecory for a single flat output.
 *****************************************************************************/
typedef struct trajectory_s
{
  bool circular;                                        // Set to true if the trajectory should restart upon terminating
  bool status;                                          // Determines if the trajectory is complete (being loaded or missing data)
  bool isset [MAX_TRAJECTORY_ENTRIES];                  // Used to check that a complete trajectory has been loaded
  float data [MAX_TRAJECTORY_ENTRIES][MAX_DIMENSION];   // trajectory data
  float time [MAX_TRAJECTORY_ENTRIES];                  // time of each trajectory component [s]
  uint8_t type [MAX_TRAJECTORY_ENTRIES];                // trajectory type
  float startTime;                                      // start time of the trajectory [s]
  uint32_t syncronizationTick;                          // time at which the object was last synchronized
  uint8_t numberOfEntries;                              // number of splines/points
} trajectory_t;

/******************************************************************************
 * Controller mode object containing an integer representation of the mode and
 * The tick at which the controller mode was most recently switched.
 *****************************************************************************/
typedef struct M_mode_s{
  uint32_t startTick;
  union {
    struct {
      uint8_t current;
      uint8_t previous;
    } controller;
    uint8_t controllerModes[2];
  };
  union {
    struct {
      bool any;                    // Set to 1 if the flatness is to be invoked
      bool eulerAngles;            // Compute an euler angle parametrization
      bool quaternion;             // Compute a quaternion parametrization
      bool bodyrate;               // Compute the body rates
      bool bodyacceleration;       // Compute the body accelerations
      bool torque;                 // Compute the body torque
    } flatness;
    bool flatnessModes[6];
  };
} M_mode_t;

/******************************************************************************
 * Control signal object, containing rotor speeds, PWM and thrust/torque vector
 * for debugging purposes.
 *****************************************************************************/
typedef struct M_control_s
{
  float PWM [4];               // PWM duty cycles (1,2,3,4)        \in[0,1]
  uint16_t Power [4];           // The motor power  (1,2,3,4)       \in[0,65535]
  float omegasq [4];           // Rotor speeds squared (1,2,3,4)   [rad^2/s^2]
  float thrust;                // Thrust along z_B (body frame)    [N]
  point_t torque;              // Torques                          [Nm]
} M_control_t;

/******************************************************************************
 * Control error object, used for debuggng controllers.
 *****************************************************************************/
typedef struct control_error_s
{
  point_t position;      // Error in translational position
  point_t velocity;      // Error in translational velcities
  attitude_t eulerAngles;// Error in euler angles  (if applicable)
  attitude_t eulerRates; // Error in euler angular rates (if applicable)
  point_t rotation;      // SO(3) rotation (if applicable)
  point_t bodyrate;      // Error in body frame rates (if applicable)
} control_error_t;


#endif //__ME_TYPES_H__
