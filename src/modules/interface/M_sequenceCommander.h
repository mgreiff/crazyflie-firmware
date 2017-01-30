/******************************************************************************
 * The sequence commander loads trajectories thorugh the parameter framework,
 * (should be revised to use packets instead) allows the loading and evaluation
 * of trajectories in the flat output space gamma=(x,y,z,yaw) \in\mathbb{R]^4.
 *
 * ~~~ Supported reference trajectory types ~~~
 * (0) Points in R^4 \in\mathbb{R]^4 corresponding to the flat outputs
 * (1) Polynomial trajectories in \in\mathbb{R]^4
 * (2) Basic sinusoid functions in \in\mathbb{R]^4
 * (3) Bezier cureves functions in \in\mathbb{R]^4 
 *
 * Each trajectory is represented as a trajectory object containing the a data
 * field which is specific for each trajectory type (see M_sequenceCommander.c
 * at the bottom of the script for clarifications as to the structure).
 *
 * In order to generate the flat output and their derivatives for flatness
 * generation, we require the fourth derivative of the trajectory. While simple
 * in the case of polynomial splines and sinusoidal functions, special care
 * must be taken when handling points. In the current implementation the
 * fourth order system
 *
 *          a^4 * s^i
 *   G(s) = ---------
 *          (a + s)^4
 *
 * is discretized with a = 7.5 at 50 Hz and simulated with (0), the point
 * trajectory, to provide well defined derivative terms.
 *
 * TODO: Rewrite loading of trajectory using packets
 * TODO: Support Bezier curves
 * TODO: Remove trajectory example
 *
 * ~~~ Changelog ~~~
 * Editor        | Date       | Description
 *-----------------------------------------------------------------------------
 * Marcus Greiff | 28-11-2016 | Initial commit with dummy example so that
 *               |            | may be generated without prioviding a
 *               |            | trajectory file.
 *-----------------------------------------------------------------------------
 * Marcus Greiff | 11-12-2016 | Updated with packets for synchronization,
 *               |            | trajectory data and direct data.
 *****************************************************************************/

#ifndef __M_SEQUENCECOMMANDER_H__
#define __M_SEQUENCECOMMANDER_H__

#include <stdbool.h>
#include <stdint.h>

#include "M_types.h"
#include "crtp.h"

// ~~~ Packet definitions ~~~

// Synchronization data
typedef struct {         // total size: 18
  uint8_t packetType;    // Type, used to distinguish between package types (= 0 here)
  uint8_t synchronize;   // Specifies if the packet should be used to clear or to synchronize the data
  uint8_t circular[4];   // If the trajectories should be circular i.e. starting over
  uint8_t number[4];     // The trajectory mode along a specific dimension
  uint16_t time [4];     // Time to start after receiving the synchronization package

} __attribute__((packed)) synchronizationPacket_t;

// Trajectory data
typedef struct {         // total size: 18
  uint8_t packetType;    // PacketType, used to distinguish between package types (= 1 here)
  uint16_t data [6];     // Trajectory data 
  uint16_t time;         // Trajectory time 
  uint8_t index;         // Contains index of current element in the complete trajectory 
  uint8_t dimension;     // Contains the dimension of the current trajectory 
  uint8_t number;        // Contains total number of elements in the complete trajectory 
  uint8_t type;          // Contains total type of the trajectory {0,1,2,3,4}
} __attribute__((packed)) trajectoryPacket_t;

typedef struct {
  uint16_t pos;          // use uint16_t to hold float16_t
  uint16_t vel;
  uint16_t acc;
  uint16_t jerk;
} __attribute__((packed)) crtpLinearControlReference_t;

typedef struct {
  uint16_t pos;          // use uint16_t to hold float16_t
  uint16_t vel;
} __attribute__((packed)) crtpAngularControlReference_t;

// Data of point in flat output space
typedef struct {                     // total size: 30 (should be OK?)
  uint8_t packetType;                // Type, used to distinguish between package types (= 2 here)
  uint8_t mode;                    // TODO: Don't quite know whet this does (possibly remove @Mike?)
  crtpLinearControlReference_t x;    // size 8
  crtpLinearControlReference_t y;    // size 8
  crtpLinearControlReference_t z;    // size 8
  crtpAngularControlReference_t yaw; // size 4
} __attribute__((packed)) pointPacket_t;

// ~~~ Functions ~~~

// Update setpoint based on the current controller mode
void commanderGetSetpoint_M(M_setpoint_t *setpoint,
                            M_mode_t *mode,
                            M_state_t *state,
                            uint32_t currentTick);

// Initialization
void commanderInit_M();

// Decoding of packets
bool decodePacket(CRTPPacket *pk, M_mode_t *mode, M_setpoint_t *setpoint);

// Evaluate point trajectory object at the current time
void eval_point_path(M_setpoint_t *setpoint,
                     int index,
                     int dimension);

// Evaluate polynomial trajectory object at the current time
void eval_polynomial(M_setpoint_t *setpoint,
                     float startTime,
                     float currentTime,
                     int index,
                     int dimension);

// Evaluate sinusodial trajectory object at the current time
void eval_function(M_setpoint_t *setpoint,
                   float currentTime,
                   int index,
                   int dimension);

// Evaluate Bezier trajectory object at the current time
void eval_bezier(M_setpoint_t *setpoint,
                 float startTime,
                 float currentTime,
                 int index,
                 int dimension);

// Evaluates polynomial
float polyval(float coeff[],
              float tp,
              int order);

// Assert that the trajectory has been fully loaded
bool assertTrajectoryLoaded(int dimension);

// Clear trajectry object
void clearTrajectory(int dimension);

#define RATE_DO_EXECUTE_OFFSET(RATE_HZ, OFFSET, TICK) (((TICK-OFFSET) % (RATE_MAIN_LOOP / RATE_HZ)) == 0)

#endif //__ME_SEQUENCECOMMANDER_H__
