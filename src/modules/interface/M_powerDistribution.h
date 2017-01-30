/******************************************************************************
 * This script allows for motor control using known SI-units and mappings, the
 * hope is to get intermediary readings with sensible units while also
 * saturating the controlsignals in a way that allows for a completely
 * controllable system, regardless of attitude and active controllers.
 *
 * ~~~ Changelog ~~~
 * Editor        | Date       | Description
 *-----------------------------------------------------------------------------
 * Marcus Greiff | 28-11-2016 | Initial commit (not yet tested for flying)
 *****************************************************************************/
 
#ifndef __M_POWER_DISTRIBUTION_H__
#define __M_POWER_DISTRIBUTION_H__

#include <stdbool.h>
#include <stdint.h>
#include <M_types.h>

// Compute and set PWM output
void powerDistribution_M(M_mode_t *mode, M_control_t *controlsignal);

// Saturate control signals
void saturate_control(M_mode_t *mode, M_control_t *controlsignal);

// Initialize motors
void powerDistributionInit_M();

#endif // __M_POWER_DISTRIBUTION_H__
