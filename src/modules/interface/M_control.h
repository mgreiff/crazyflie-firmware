// Avoid redefinitions
#ifndef __M_CONTROLLER_H__
#define __M_CONTROLLER_H__

#include <stdbool.h>
#include <stdint.h>

#include "M_types.h"

// High level computation of control signals
void stateController_M(M_control_t *controlsignal,
                       M_setpoint_t *setpoint,
                       M_mode_t *mode,
                       M_state_t *state,
                       uint32_t tick);

// Attitude PD control (used to hover and in external/internal mode)
void PD_elevation_control(M_control_t *controlsignal,
                          M_setpoint_t *setpoint,
                          M_state_t *state,
                          uint32_t tick);

// Elevation PD control (used to hover and in external/internal mode)
void PD_attitude_control(M_control_t *controlsignal,
                         M_setpoint_t *setpoint,
                         M_state_t *state,
                         uint32_t tick);

// Position PID control (used to hover and follow trajectories in internal mode)
void PID_translation_control(M_setpoint_t *setpoint,
                             M_state_t *state,
                             uint32_t tick);

// TI-LQR control with euler angle parametrization (used to hover in internal
// mode)
void LQR_hover_control(M_control_t *controlsignal,
                       M_setpoint_t *setpoint,
                       M_state_t *state,
                       uint32_t tick);

// TI-LQR control with euler angle parametrization (used to follow trajectories
// in internal mode)
void LQR_full_control(M_control_t *controlsignal,
                      M_setpoint_t *setpoint,
                      M_state_t *state,
                      uint32_t tick);
                       
// Geometric translation tracking control on SE(3) (follow trajectories in
// internal mode)
void SE3_full_control(M_control_t *controlsignal,
                      M_setpoint_t *setpoint,
                      M_state_t *state,
                      uint32_t tick);

// Initialization
void stateControllerInit_M();

// Compute cross product
void crossProd(float u[], float v[], float res[]);

// Normalize vector
void normalizeVector(float u[], float res[]);

#endif //__M_CONTROLLER_H__
