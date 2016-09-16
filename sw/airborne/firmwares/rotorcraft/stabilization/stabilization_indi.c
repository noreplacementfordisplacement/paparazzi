/*Copyright (C) Ewoud Smeur <ewoud_smeur@msn.com>
		 Daan Hoppener <daanhoppener@gmail.com> (WLS CONTROL)
 * MAVLab Delft University of Technology
 *
 * This file is part of paparazzi.
 *
 * paparazzi is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2, or (at your option)
 * any later version.
 *
 * paparazzi is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with paparazzi; see the file COPYING.  If not, write to
 * the Free Software Foundation, 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

/** @file stabilization_attitude_quat_indi.c
 * @brief MAVLab Delft University of Technology
 * This control algorithm is Incremental Nonlinear Dynamic Inversion (INDI)
 *
 * This is a simplified implementation of the publication in the
 * journal of Control Guidance and Dynamics: Adaptive Incremental Nonlinear
 * Dynamic Inversion for Attitude Control of Micro Aerial Vehicles
 * http://arc.aiaa.org/doi/pdf/10.2514/1.G001490
 */

#include "firmwares/rotorcraft/stabilization/stabilization_indi.h"
#include "firmwares/rotorcraft/stabilization/stabilization_attitude.h"
#include "firmwares/rotorcraft/stabilization/stabilization_attitude_rc_setpoint.h"
#include "firmwares/rotorcraft/stabilization/stabilization_attitude_quat_transformations.h"

#include "state.h"
#include "generated/airframe.h"
#include "paparazzi.h"
#include "subsystems/radio_control.h"

// Actuator feedback headers
#include "boards/bebop/actuators.h"
#include "subsystems/commands.h"

// Matrix Inverse headers
#include "math/pprz_algebra_float.h"
#include "math/pprz_algebra_int.h"


// WLS Control allocator libraries
 #include "wls/wls_alloc.h"
// #include "wls/qr_solve.h"
// #include "wls/r8lib.h"
#include <stdio.h>

// WLS defines
#define MAX_MOTOR_WLS 8600 //From experiment data we know the Bebop does not really go over 8100 RPM
#define MIN_MOTOR_WLS 3000 //Minimum RPM from the motors
#define MARINUS 4 //MxN matrix control effectiveness, MARINUS = M: the dimension of the control objective
#define NICO 4 //MxN matrix control effectiveness, NICO = N: the dimension of the actuators


#if !defined(STABILIZATION_INDI_ACT_DYN_P) && !defined(STABILIZATION_INDI_ACT_DYN_Q) && !defined(STABILIZATION_INDI_ACT_DYN_R)
#error You have to define the first order time constant of the actuator dynamics!
#endif

// these parameters are used in the filtering of the angular acceleration
// define them in the airframe file if different values are required
#ifndef STABILIZATION_INDI_FILT_OMEGA
#define STABILIZATION_INDI_FILT_OMEGA 50.0 //should be 50
#endif

#define STABILIZATION_INDI_FILT_OMEGA2 (STABILIZATION_INDI_FILT_OMEGA*STABILIZATION_INDI_FILT_OMEGA)

#ifndef STABILIZATION_INDI_FILT_ZETA
#define STABILIZATION_INDI_FILT_ZETA 0.55
#endif

//FIXME: Filter values for MANUAL WLS G1, G2 Identification
#define IDENTIFICATION_INDI_FILT_OMEGA 25
#define IDENTIFICATION_INDI_FILT_OMEGA2 625
#define IDENTIFICATION_FILT_ZETA 0.55


// the yaw sometimes requires more filtering
#ifndef STABILIZATION_INDI_FILT_OMEGA_R
#define STABILIZATION_INDI_FILT_OMEGA_R STABILIZATION_INDI_FILT_OMEGA
#endif

#ifndef STABILIZATION_INDI_FILT_ZETA_R
#define STABILIZATION_INDI_FILT_ZETA_R STABILIZATION_INDI_FILT_ZETA
#endif

#ifndef STABILIZATION_INDI_MAX_RATE
#define STABILIZATION_INDI_MAX_RATE 6.0
#endif

#if STABILIZATION_INDI_USE_ADAPTIVE
#warning "Use caution with adaptive indi. See the wiki for more info"
#endif

#ifndef STABILIZATION_INDI_MAX_R
#define STABILIZATION_INDI_MAX_R STABILIZATION_ATTITUDE_SP_MAX_R
#endif

#define VECT4_INTEGRATE(_a, _b, _c) { \
_a[0] = _a[0] + _b[0]/_c; \
_a[1] = _a[1] + _b[1]/_c; \
_a[2] = _a[2] + _b[2]/_c; \
_a[3] = _a[3] + _b[3]/_c; \
}

#define RATES_INTEGRATE(_a, _b, _c) { \
_a.p = _a.p + _b.p/_c; \
_a.q = _a.q + _b.q/_c; \
_a.r = _a.r + _b.r/_c; \
}

struct Int32Eulers stab_att_sp_euler;
struct Int32Quat   stab_att_sp_quat;

//Variable headers used in WLS control allocator
float u_actuators[NICO] = {0.0, 0.0, 0.0, 0.0};
float udot_actuators[NICO] = {0.0, 0.0, 0.0, 0.0};
float udotdot_actuators[NICO] = {0.0, 0.0, 0.0, 0.0};
float u_act_dyn_actuators[NICO] = {0.0, 0.0, 0.0, 0.0};
float u_cmd[NICO] = {0.0 , 0.0, 0.0, 0.0};
float umin[NICO] = {0.0 , 0.0, 0.0, 0.0};
float umax[NICO] = {0.0 , 0.0, 0.0, 0.0};

int32_t in_cmd_wls[NICO]; //FIXME: "Jerryrig" to communicate with motor_mixing

float wls_temp_thrust = 0; //Incremental thruststatic float Wv[MARINUS] = {10, 10, 1, 5}; //State prioritization {W Roll, W pitch, W yaw, TOTAL THRUST}
static float Wv[MARINUS] = {100, 100, 1, 5}; //State prioritization {W Roll, W pitch, W yaw, TOTAL THRUST}
bool regindi = true; // Boolean to indicate if regular INDI should be run (unconstrained)
bool wls_adaptive  = false; //Boolean to indicate if G1 and G2 are going to be adaptive
/*float B_tmp[MARINUS][NICO] = {{0.0210, -0.0210, -0.0210, 0.0210},
					{0.015, 0.015, -0.015, -0.015},
					{-0.081, 0.081, -0.081, 0.081},
					{0.25, 0.25, 0.25, 0.25}}; //G1 + G2 for WLS control allocator */

// Previously identified Pseudo Inverse of G1+G2
static float B_pinv[4][3] = {{ 13.0152,   19.0064,   -3.5344}, {-12.4792, 20.7766, 3.6274}, {-12.2737 ,-21.2221  ,-3.5606}, {11.8252, -20.0975 ,3.7451}};

float** Bwls;

// For estimating the actuator effectiveness
struct FloatRates rate_estimation = {0., 0., 0.};
struct FloatRates ratedot_estimation = {0., 0., 0.};
struct FloatRates ratedotdot_estimation = {0., 0., 0.};
float u_estimate[NICO] = {0.0, 0.0, 0.0, 0.0};
float udot_estimate[NICO] = {0.0, 0.0, 0.0, 0.0};
float udotdot_estimate[NICO] = {0.0, 0.0, 0.0, 0.0};

//FIXME: BSVARS to estimate the actuator effectiveness
float G1G2[3][4] = {{0.0 , 0.0, 0.0 , 0.0 },
{0.0 , 0.0, 0.0 , 0.0 },
{0.0, 0.0, 0.0, 0.0}};
float G1_new[3][4] = {{0.0 , 0.0, 0.0 , 0.0 },
{0.0 , 0.0, 0.0 , 0.0 },
{0.0, 0.0, 0.0, 0.0}};
float G2_new[4] = {0.0, 0.0, 0.0, 0.0};

// Logged data
float G1wls[3][4] = {{20.570812, -20.139616, -20.075964, 19.833261},
{11.753228, 12.406771, -12.40436, -12.717836},
{-1.774367, 1.461394, -0.34027, 0.596049}};

/*float G1wls[3][4] = {{20.570812, -20.139616, 0, 19.833261},
{11.753228, 12.406771, 0, -12.717836},
{-1.774367, 1.461394, 0, 0.596049}};*/


float G2wls[4] = {-64.577644, 63.09156, -66.577477, 73.646568};

/*float G2wls[4] = {-64.577644, 63.09156, 0, 73.646568};*/


float mu1 = 0.00001;
float mu2 = 0.00001*600.0;
float dx_estimation[3] = {0.0, 0.0, 0.0};
float du_estimation[4] = {0.0, 0.0, 0.0, 0.0};
float ddu_estimation[4] = {0.0, 0.0, 0.0, 0.0};
float dx_error_disp[3];

//Functions
static int32_t stabilization_att_indi_cmd[COMMANDS_NB];
static inline void stabilization_indi_calc_cmd(int32_t indi_commands[], struct Int32Quat *att_err, bool rate_control);
static void stabilization_indi_second_order_filter_init(struct IndiFilter *filter, float omega, float zeta, float omega_r);
static void stabilization_indi_second_order_filter(struct IndiFilter *filter, struct FloatRates *input);
static inline void lms_estimation(void);

//Control Effectiveness Estimator Functions
static void calc_g1_element(float du_norm, float dx_error, int8_t i, int8_t j, float mu_extra);
static void calc_g2_element(float dx_error, int8_t j, float mu_extra);
static void lmsg1_estimation(void);

  // input u (output of WLS controller)
  float u[NICO] = {0.0, 0.0, 0.0, 0.0};

#define INDI_EST_SCALE 0.001 //The G values are scaled to avoid numerical problems during the estimation
struct IndiVariables indi = {
  .max_rate = STABILIZATION_INDI_MAX_RATE,
  .attitude_max_yaw_rate = STABILIZATION_INDI_MAX_RATE,

  .g1 = {STABILIZATION_INDI_G1_P, STABILIZATION_INDI_G1_Q, STABILIZATION_INDI_G1_R},
  .g2 = STABILIZATION_INDI_G2_R,
  .reference_acceleration = {
    STABILIZATION_INDI_REF_ERR_P,
    STABILIZATION_INDI_REF_ERR_Q,
    STABILIZATION_INDI_REF_ERR_R,
    STABILIZATION_INDI_REF_RATE_P,
    STABILIZATION_INDI_REF_RATE_Q,
    STABILIZATION_INDI_REF_RATE_R},

  /* Estimation parameters for adaptive INDI */
  .est = {
    .g1 = {
      STABILIZATION_INDI_G1_P / INDI_EST_SCALE,
      STABILIZATION_INDI_G1_Q / INDI_EST_SCALE,
      STABILIZATION_INDI_G1_R / INDI_EST_SCALE},
    .g2 = STABILIZATION_INDI_G2_R / INDI_EST_SCALE,
    .mu = STABILIZATION_INDI_ADAPTIVE_MU,
  },

#if STABILIZATION_INDI_USE_ADAPTIVE
  .adaptive = TRUE,
#else
  .adaptive = FALSE,
#endif
};

#if PERIODIC_TELEMETRY
#include "subsystems/datalink/telemetry.h"

static void send_att_indi(struct transport_tx *trans, struct link_device *dev)
{
  //The estimated G values are scaled, so scale them back before sending
  struct FloatRates g1_disp;
  RATES_SMUL(g1_disp, indi.est.g1, INDI_EST_SCALE);
  float g2_disp = indi.est.g2 * INDI_EST_SCALE;

  pprz_msg_send_STAB_ATTITUDE_INDI(trans, dev, AC_ID,
                                   &indi.rate.dx.p,
                                   &indi.rate.dx.q,
                                   &indi.rate.dx.r,
                                   &indi.angular_accel_ref.p,
                                   &indi.angular_accel_ref.q,
                                   &indi.angular_accel_ref.r,
                                   &g1_disp.p,
                                   &g1_disp.q,
                                   &g1_disp.r,
                                   &g2_disp);
}
#endif

void stabilization_indi_init(void)
{
  // Initialize filters
  stabilization_indi_second_order_filter_init(&indi.rate, STABILIZATION_INDI_FILT_OMEGA, STABILIZATION_INDI_FILT_ZETA, STABILIZATION_INDI_FILT_OMEGA_R);
  stabilization_indi_second_order_filter_init(&indi.u, STABILIZATION_INDI_FILT_OMEGA, STABILIZATION_INDI_FILT_ZETA, STABILIZATION_INDI_FILT_OMEGA_R);
  stabilization_indi_second_order_filter_init(&indi.est.rate, 10.0, 0.8, 10.0); //FIXME: no magic number
  stabilization_indi_second_order_filter_init(&indi.est.u, 10.0, 0.8, 10.0); //FIXME: no magic number

#if PERIODIC_TELEMETRY
  register_periodic_telemetry(DefaultPeriodic, PPRZ_MSG_ID_STAB_ATTITUDE_INDI, send_att_indi);
#endif
}

void stabilization_indi_enter(void)
{
  /* reset psi setpoint to current psi angle */
  stab_att_sp_euler.psi = stabilization_attitude_get_heading_i();
 
  FLOAT_RATES_ZERO(indi.rate.x);
  FLOAT_RATES_ZERO(indi.rate.dx);
  FLOAT_RATES_ZERO(indi.rate.ddx);
  FLOAT_RATES_ZERO(indi.angular_accel_ref);
  FLOAT_RATES_ZERO(indi.du);
  FLOAT_RATES_ZERO(indi.u_act_dyn);
  FLOAT_RATES_ZERO(indi.u_in);
  FLOAT_RATES_ZERO(indi.u.x);
  FLOAT_RATES_ZERO(indi.u.dx);
  FLOAT_RATES_ZERO(indi.u.ddx);
}

void stabilization_indi_set_failsafe_setpoint(void)
{
  /* set failsafe to zero roll/pitch and current heading */
  int32_t heading2 = stabilization_attitude_get_heading_i() / 2;
  PPRZ_ITRIG_COS(stab_att_sp_quat.qi, heading2);
  stab_att_sp_quat.qx = 0;
  stab_att_sp_quat.qy = 0;
  PPRZ_ITRIG_SIN(stab_att_sp_quat.qz, heading2);
}

void stabilization_indi_set_rpy_setpoint_i(struct Int32Eulers *rpy)
{
  // stab_att_sp_euler.psi still used in ref..
  stab_att_sp_euler = *rpy;

  quat_from_rpy_cmd_i(&stab_att_sp_quat, &stab_att_sp_euler);
}

void stabilization_indi_set_earth_cmd_i(struct Int32Vect2 *cmd, int32_t heading)
{
  // stab_att_sp_euler.psi still used in ref..
  stab_att_sp_euler.psi = heading;

  // compute sp_euler phi/theta for debugging/telemetry
  /* Rotate horizontal commands to body frame by psi */
  int32_t psi = stateGetNedToBodyEulers_i()->psi;
  int32_t s_psi, c_psi;
  PPRZ_ITRIG_SIN(s_psi, psi);
  PPRZ_ITRIG_COS(c_psi, psi);
  stab_att_sp_euler.phi = (-s_psi * cmd->x + c_psi * cmd->y) >> INT32_TRIG_FRAC;
  stab_att_sp_euler.theta = -(c_psi * cmd->x + s_psi * cmd->y) >> INT32_TRIG_FRAC;

  quat_from_earth_cmd_i(&stab_att_sp_quat, cmd, heading);
}

static inline void stabilization_indi_calc_cmd(int32_t indi_commands[], struct Int32Quat *att_err, bool rate_control)
{
  /* Propagate the second order filter on the gyroscopes */
  struct FloatRates *body_rates = stateGetBodyRates_f();
  stabilization_indi_second_order_filter(&indi.rate, body_rates);

  //The rates used for feedback are by default the measured rates. If needed they can be filtered (see below)
  struct FloatRates rates_for_feedback;
  RATES_COPY(rates_for_feedback, (*body_rates));

  //If there is a lot of noise on the gyroscope, it might be good to use the filtered value for feedback.
  //Note that due to the delay, the PD controller can not be as aggressive.
#if STABILIZATION_INDI_FILTER_ROLL_RATE
  rates_for_feedback.p = indi.rate.x.p;
#endif
#if STABILIZATION_INDI_FILTER_PITCH_RATE
  rates_for_feedback.q = indi.rate.x.q;
#endif
#if STABILIZATION_INDI_FILTER_YAW_RATE
  rates_for_feedback.r = indi.rate.x.r;
#endif

  indi.angular_accel_ref.p = indi.reference_acceleration.err_p * QUAT1_FLOAT_OF_BFP(att_err->qx)
                             - indi.reference_acceleration.rate_p * rates_for_feedback.p;

  indi.angular_accel_ref.q = indi.reference_acceleration.err_q * QUAT1_FLOAT_OF_BFP(att_err->qy)
                             - indi.reference_acceleration.rate_q * rates_for_feedback.q;

  //This separates the P and D controller and lets you impose a maximum yaw rate.
  float rate_ref_r = indi.reference_acceleration.err_r * QUAT1_FLOAT_OF_BFP(att_err->qz)/indi.reference_acceleration.rate_r;
  BoundAbs(rate_ref_r, indi.attitude_max_yaw_rate);
  indi.angular_accel_ref.r = indi.reference_acceleration.rate_r * (rate_ref_r - rates_for_feedback.r);

  /* Check if we are running the rate controller and overwrite */
  if(rate_control) {
    indi.angular_accel_ref.p =  indi.reference_acceleration.rate_p * ((float)radio_control.values[RADIO_ROLL]  / MAX_PPRZ * indi.max_rate - body_rates->p);
    indi.angular_accel_ref.q =  indi.reference_acceleration.rate_q * ((float)radio_control.values[RADIO_PITCH] / MAX_PPRZ * indi.max_rate - body_rates->q);
    indi.angular_accel_ref.r =  indi.reference_acceleration.rate_r * ((float)radio_control.values[RADIO_YAW]   / MAX_PPRZ * indi.max_rate - body_rates->r);
  }

  //Increment in angular acceleration requires increment in control input
  //G1 is the control effectiveness. In the yaw axis, we need something additional: G2.
  //It takes care of the angular acceleration caused by the change in rotation rate of the propellers
  //(they have significant inertia, see the paper mentioned in the header for more explanation)
  indi.du.p = 1.0 / indi.g1.p * (indi.angular_accel_ref.p - indi.rate.dx.p);
  indi.du.q = 1.0 / indi.g1.q * (indi.angular_accel_ref.q - indi.rate.dx.q);
  indi.du.r = 1.0 / (indi.g1.r + indi.g2) * (indi.angular_accel_ref.r - indi.rate.dx.r + indi.g2 * indi.du.r);
  
  //---------------------------------------------------------------------
  //---------------------------------------------------------------------
  // 				0 0
  // WLS Control Allocator       x
  //                            0 0
  //---------------------------------------------------------------------
  //---------------------------------------------------------------------
   
   // Current actuator state (RPM FB) (IN MM coordinates)
   u_act_dyn_actuators[0] = (actuators_bebop.rpm_obs[0]- MIN_MOTOR_WLS)*(MAX_PPRZ/9000); 
   u_act_dyn_actuators[1] = (actuators_bebop.rpm_obs[1]- MIN_MOTOR_WLS)*(MAX_PPRZ/9000); 
   u_act_dyn_actuators[2] = (actuators_bebop.rpm_obs[2]- MIN_MOTOR_WLS)*(MAX_PPRZ/9000); 
   u_act_dyn_actuators[3] = (actuators_bebop.rpm_obs[3]- MIN_MOTOR_WLS)*(MAX_PPRZ/9000); 

    // Filter actuators
    VECT4_INTEGRATE(u_actuators, udot_actuators, 512.0);
    VECT4_INTEGRATE(udot_actuators, udotdot_actuators, 512.0);
    
   // Update filter states
   udotdot_actuators[0] = -udot_actuators[0] * 2*STABILIZATION_INDI_FILT_ZETA*STABILIZATION_INDI_FILT_OMEGA + (u_act_dyn_actuators[0] - u_actuators[0])*STABILIZATION_INDI_FILT_OMEGA2;
   udotdot_actuators[1] = -udot_actuators[1] * 2*STABILIZATION_INDI_FILT_ZETA*STABILIZATION_INDI_FILT_OMEGA + (u_act_dyn_actuators[1] - u_actuators[1])*STABILIZATION_INDI_FILT_OMEGA2;
   udotdot_actuators[2] = -udot_actuators[2] * 2*STABILIZATION_INDI_FILT_ZETA*STABILIZATION_INDI_FILT_OMEGA + (u_act_dyn_actuators[2] - u_actuators[2])*STABILIZATION_INDI_FILT_OMEGA2;
   udotdot_actuators[3] = -udot_actuators[3] * 2*STABILIZATION_INDI_FILT_ZETA*STABILIZATION_INDI_FILT_OMEGA + (u_act_dyn_actuators[3] - u_actuators[3])*STABILIZATION_INDI_FILT_OMEGA2;

   // Bound total input
   Bound(u_actuators[0], 0, MAX_MOTOR_WLS);
   Bound(u_actuators[1], 0, MAX_MOTOR_WLS);
   Bound(u_actuators[2], 0, MAX_MOTOR_WLS);
   Bound(u_actuators[3], 0, MAX_MOTOR_WLS);
 
  //MAX possible MINUMUM increment 
  umin[0] = -u_actuators[0]; 
  umin[1] = -u_actuators[1]; 
  umin[2] = -u_actuators[2]; 
  umin[3] = -u_actuators[3];
   
  //MAX possible MAXIMUM increment
  umax[0] = MAX_MOTOR_WLS - u_actuators[0];
  umax[1] = MAX_MOTOR_WLS - u_actuators[1];
  umax[2] = MAX_MOTOR_WLS - u_actuators[2];
  umax[3] = MAX_MOTOR_WLS - u_actuators[3];

//  static float** B; // Initialize **B
  
  //FIXME: Not really elegant
//  if (wls_adaptive == false){
  //Static control effectiveness matrix for the WLS control allocator, last row is not considered INDI, is only used for the Thrust command
	float B_tmp[MARINUS][NICO]  = {{G1wls[0][0]*INDI_EST_SCALE, G1wls[0][1]*INDI_EST_SCALE, G1wls[0][2]*INDI_EST_SCALE, G1wls[0][3]*INDI_EST_SCALE},{G1wls[1][0]*INDI_EST_SCALE, G1wls[1][1]*INDI_EST_SCALE, G1wls[1][2]*INDI_EST_SCALE, G1wls[1][3]*INDI_EST_SCALE}, {G1wls[2][0]*INDI_EST_SCALE + G2wls[0]*INDI_EST_SCALE, G1wls[2][1]*INDI_EST_SCALE + G2wls[1]*INDI_EST_SCALE, G1wls[2][2]*INDI_EST_SCALE + G2wls[2]*INDI_EST_SCALE, G1wls[2][3]*INDI_EST_SCALE + G2wls[3]*INDI_EST_SCALE}, {0.297628763, 0.2770495344, 0.1933141364, 0.2320075662}};
//	}
 // FIXME: Allocate and free **B each time... not very efficient
 	Bwls = (float**)calloc(MARINUS, sizeof(float*)); 
 	   for (int i = 0; i < MARINUS; i++) {
 	       Bwls[i] = (float*)calloc(NICO, sizeof(float*));
 	       for (int j = 0; j < NICO; j++) Bwls[i][j] = B_tmp[i][j];
 	   }

  //FIXME: Jerryrig "incremental Thrust"
  wls_temp_thrust = stabilization_cmd[COMMAND_THRUST] - (u_actuators[0] + u_actuators[1] + u_actuators[2] +u_actuators[3])/4; // incremental thrust

  float v[MARINUS];

  v[0] = (indi.angular_accel_ref.p - indi.rate.dx.p);
  v[1] = (indi.angular_accel_ref.q - indi.rate.dx.q);  
  v[2] = (indi.angular_accel_ref.r - indi.rate.dx.r + wlsg2_fb);
  v[3] = wls_temp_thrust;

if (regindi == false){
  // WLS Control Allocator
  // u: output incremental actuator commands, v: control objective, umin,umax: maximum possible positive and negative increments, B: control effectiveness matrix (needs to be of ** type), Wv: priority of WLS control allocator, gamma (1000): weight on control objective solution (should be => 1000), rmax (100): maximum number of iterations
  wls_alloc(u,v,umin,umax,Bwls,NICO,MARINUS,0,0,Wv,0,0,10000,100);
  }

//FIXME: Free Bwls
  free(Bwls);

if (regindi == true){
// Calculate increments with regular INDI
	u[0] = (B_pinv[0][0] * v[0]) + (B_pinv[0][1] * v[1]) + (B_pinv[0][2] * v[2]);
  	u[1] = (B_pinv[1][0] * v[0]) + (B_pinv[1][1] * v[1]) + (B_pinv[1][2] * v[2]);
  	u[2] = (B_pinv[2][0] * v[0]) + (B_pinv[2][1] * v[1]) + (B_pinv[2][2] * v[2]);
  	u[3] = (B_pinv[3][0] * v[0]) + (B_pinv[3][1] * v[1]) + (B_pinv[3][2] * v[2]);
}	

 // Compute G2 feedback for INDI on yaw axis
  wlsg2_fb = (G2wls[0]*INDI_EST_SCALE*u[0] + G2wls[1]*INDI_EST_SCALE*u[1] + G2wls[2]*INDI_EST_SCALE*u[2] + G2wls[3]*INDI_EST_SCALE*u[3]);

//========================================================================================================================================================
//FIXME: MANUAL ADAPTIVE G1 && G2 FILTERING (less elegent than the other in this document)
//Sensor dynamics and Filtering (same filter as on gyro measurements) MANUAL ACTUATOR EFFECTIVENESS
  RATES_INTEGRATE(rate_estimation,ratedot_estimation,512.0);
  RATES_INTEGRATE(ratedot_estimation,ratedotdot_estimation,512.0);

  ratedotdot_estimation.p = -ratedot_estimation.p * 2*IDENTIFICATION_FILT_ZETA*IDENTIFICATION_INDI_FILT_OMEGA + (stateGetBodyRates_f()->p - rate_estimation.p)*IDENTIFICATION_INDI_FILT_OMEGA2;
  ratedotdot_estimation.q = -ratedot_estimation.q * 2*IDENTIFICATION_FILT_ZETA*IDENTIFICATION_INDI_FILT_OMEGA + (stateGetBodyRates_f()->q - rate_estimation.q)*IDENTIFICATION_INDI_FILT_OMEGA2;
  ratedotdot_estimation.r = -ratedot_estimation.r * 2*IDENTIFICATION_FILT_ZETA*IDENTIFICATION_INDI_FILT_OMEGA + (stateGetBodyRates_f()->r - rate_estimation.r)*IDENTIFICATION_INDI_FILT_OMEGA2;

    // Filter actuators
    VECT4_INTEGRATE(u_estimate, udot_estimate, 512.0);
    VECT4_INTEGRATE(udot_estimate, udotdot_estimate, 512.0);
    
   // Update filter states
   udotdot_estimate[0] = -udot_estimate[0] * 2*IDENTIFICATION_FILT_ZETA*IDENTIFICATION_INDI_FILT_OMEGA + (u_act_dyn_actuators[0] - u_estimate[0])*IDENTIFICATION_INDI_FILT_OMEGA2;
   udotdot_estimate[1] = -udot_estimate[1] * 2*IDENTIFICATION_FILT_ZETA*IDENTIFICATION_INDI_FILT_OMEGA  + (u_act_dyn_actuators[1] - u_estimate[1])*IDENTIFICATION_INDI_FILT_OMEGA2;
   udotdot_estimate[2] = -udot_estimate[2] * 2*IDENTIFICATION_FILT_ZETA*IDENTIFICATION_INDI_FILT_OMEGA + (u_act_dyn_actuators[2] - u_estimate[2])*IDENTIFICATION_INDI_FILT_OMEGA2;
   udotdot_estimate[3] = -udot_estimate[3] * 2*IDENTIFICATION_FILT_ZETA*IDENTIFICATION_INDI_FILT_OMEGA + (u_act_dyn_actuators[3] - u_estimate[3])*IDENTIFICATION_INDI_FILT_OMEGA2;
//========================================================================================================================================================

  //FIXME: This is to kill the INDI/WLS control
  if (stabilization_cmd[COMMAND_THRUST] < 300) {
	u_cmd[0] = 0;
	u_cmd[1] = 0;
	u_cmd[2] = 0;
	u_cmd[3] = 0;
	}
  else{
  // If everything is correct this can be used as direct actuator input
  u_cmd[0] = u[0] + u_actuators[0];
  u_cmd[1] = u[1] + u_actuators[1]; 
  u_cmd[2] = u[2] + u_actuators[2]; 
  u_cmd[3] = u[3] + u_actuators[3];

 	if (regindi == true){ // Add thrust to REGULAR INDI
  	float avg_u_in = (u_cmd[0] + u_cmd[1] + u_cmd[2] + u_cmd[3])/4.0;

  	if(avg_u_in > 1.0) {
    		u_cmd[0] = u_cmd[0] /avg_u_in * stabilization_cmd[COMMAND_THRUST];
    		u_cmd[1] = u_cmd[1] /avg_u_in * stabilization_cmd[COMMAND_THRUST];
    		u_cmd[2] = u_cmd[2] /avg_u_in * stabilization_cmd[COMMAND_THRUST];
    		u_cmd[3] = u_cmd[3] /avg_u_in * stabilization_cmd[COMMAND_THRUST];
		 	}
		}

    if (wls_adaptive == true){
	lmsg1_estimation();
     }	

   // Bound total output
   Bound(u_cmd[0], 0, MAX_MOTOR_WLS); //should be MAX_MOTOR_WLS (REA
   Bound(u_cmd[1], 0, MAX_MOTOR_WLS);
   Bound(u_cmd[2], 0, MAX_MOTOR_WLS);
   Bound(u_cmd[3], 0, MAX_MOTOR_WLS);
  }
 
//@@@@@@@@@@@@@@@@@@@@ Get your log on! @@@@@@@@@@@@@@@@@@@@@@@@@@

 phim = stateGetNedToBodyEulers_f() -> phi;  // Get body euler angles
 thetam = stateGetNedToBodyEulers_f() -> theta;
 psim = stateGetNedToBodyEulers_f() -> psi;

 pratem =  body_rates->p; // Get body euler rates (unfiltered)
 qratem =  body_rates->q;
 rratem =  body_rates->r;

 prateref =  indi.angular_accel_ref.p; // Get angular acc. reference 
 qrateref =  indi.angular_accel_ref.q; 
 rrateref =  indi.angular_accel_ref.r;

 vlog[0] = v[0]; // Log inputs to WLS
 vlog[1] = v[1];
 vlog[2] = v[2];
 vlog[3] = v[3];
 vlog[4] = (indi.g1.r + indi.g2) * (indi.angular_accel_ref.r - indi.rate.dx.r + indi.g2 * indi.du.r); 

 u_cmd_log[0] = u_actuators[0]; //log percieved actuator feedback
 u_cmd_log[1] = u_actuators[1];
 u_cmd_log[2] = u_actuators[2];
 u_cmd_log[3] = u_actuators[3];

 accestlog[0] = ratedotdot_estimation.p;
 accestlog[1] = ratedotdot_estimation.q;
 accestlog[2] = ratedotdot_estimation.r;

 phiref = stab_att_sp_euler.phi;  // Get body euler reference angles
 thetaref = stab_att_sp_euler.theta;
 psiref = stab_att_sp_euler.psi; 

 //---------------------------------------------------------------------
 //---------------------------------------------------------------------
 // 		               0 0
 // END of wls control          x
 //                            0 0
 //---------------------------------------------------------------------
 //---------------------------------------------------------------------
  
  //add the increment to the total control input
  indi.u_in.p = indi.u.x.p + indi.du.p;
  indi.u_in.q = indi.u.x.q + indi.du.q;
  indi.u_in.r = indi.u.x.r + indi.du.r;

  //bound the total control input
  Bound(indi.u_in.p, -4500, 4500);
  Bound(indi.u_in.q, -4500, 4500);
  Bound(indi.u_in.r, -4500, 4500);

  //Propagate input filters
  //first order actuator dynamics
  indi.u_act_dyn.p = indi.u_act_dyn.p + STABILIZATION_INDI_ACT_DYN_P * (indi.u_in.p - indi.u_act_dyn.p);
  indi.u_act_dyn.q = indi.u_act_dyn.q + STABILIZATION_INDI_ACT_DYN_Q * (indi.u_in.q - indi.u_act_dyn.q);
  indi.u_act_dyn.r = indi.u_act_dyn.r + STABILIZATION_INDI_ACT_DYN_R * (indi.u_in.r - indi.u_act_dyn.r);

  //sensor filter
  stabilization_indi_second_order_filter(&indi.u, &indi.u_act_dyn);

  //Don't increment if thrust is off
  //TODO: this should be something more elegant, but without this the inputs will increment to the maximum before
  //even getting in the air.
  if (stabilization_cmd[COMMAND_THRUST] < 300) {
    FLOAT_RATES_ZERO(indi.du);
    FLOAT_RATES_ZERO(indi.u_act_dyn);
    FLOAT_RATES_ZERO(indi.u_in);
    FLOAT_RATES_ZERO(indi.u.x);
    FLOAT_RATES_ZERO(indi.u.dx);
    FLOAT_RATES_ZERO(indi.u.ddx);
  } else {
    // only run the estimation if the commands are not zero.
    lms_estimation();
  }
  /*  INDI feedback */
  indi_commands[COMMAND_ROLL] = indi.u_in.p;
  indi_commands[COMMAND_PITCH] = indi.u_in.q;
  indi_commands[COMMAND_YAW] = indi.u_in.r;
  indi_commands[COMMAND_WLS_1] = u_cmd[0];
  indi_commands[COMMAND_WLS_2] = u_cmd[1]; 
  indi_commands[COMMAND_WLS_3] = u_cmd[2];
  indi_commands[COMMAND_WLS_4] = u_cmd[3];
}
// =======================================================================================
// :::::::::::::::: WLS Actuator Effectiveness Estimator functions :::::::::::::::::::::::
static void calc_g1_element(float du_norm, float dx_error, int8_t i, int8_t j, float mu_extra) {
  G1_new[i][j] = G1wls[i][j] - du_estimation[j]*mu1*dx_error*mu_extra;
}

static void calc_g2_element(float dx_error, int8_t j, float mu_extra) {
  G2_new[j] = G2wls[j] - ddu_estimation[j]*mu2*dx_error*mu_extra;
}

static void lmsg1_estimation(void) {
  dx_estimation[0] = ratedotdot_estimation.p;
  dx_estimation[1] = ratedotdot_estimation.q;
  dx_estimation[2] = ratedotdot_estimation.r;
  du_estimation[0] = udot_estimate[0]/1000; 
  du_estimation[1] = udot_estimate[1]/1000; 
  du_estimation[2] = udot_estimate[2]/1000; 
  du_estimation[3] = udot_estimate[3]/1000; 
  ddu_estimation[0] = udotdot_estimate[0]/1000/512.0;
  ddu_estimation[1] = udotdot_estimate[1]/1000/512.0;
  ddu_estimation[2] = udotdot_estimate[2]/1000/512.0;
  ddu_estimation[3] = udotdot_estimate[3]/1000/512.0;

   //Estimation of G
  float du_norm = du_estimation[0]*du_estimation[0] + du_estimation[1]*du_estimation[1] +du_estimation[2]*du_estimation[2] + du_estimation[3]*du_estimation[3];
  float dx_norm = dx_estimation[0]*dx_estimation[0] + dx_estimation[1]*dx_estimation[1] + dx_estimation[2]*dx_estimation[2];
//   if((dx_norm > 40000.0)) {
    for(int8_t i=0; i<3; i++) {
      float mu_extra = 1.0; // adapt according to percieved error DEFAULT mu_extra = 1.0
      float dx_error = G1wls[i][0]*du_estimation[0] + G1wls[i][1]*du_estimation[1] + G1wls[i][2]*du_estimation[2] + G1wls[i][3]*du_estimation[3] - dx_estimation[i];
      dx_error_disp[i] = dx_error;
      //if the yaw axis, also use G2
      if(i==2) {
        mu_extra = 0.3; //default mu extra = 0.3
        dx_error = dx_error + G2wls[0]*ddu_estimation[0] + G2wls[1]*ddu_estimation[1] + G2wls[2]*ddu_estimation[2] + G2wls[3]*ddu_estimation[3];
        for(int8_t j=0; j<4; j++) {
          calc_g2_element(dx_error,j, mu_extra);
        }
      }
      for(int8_t j=0; j<4; j++) {
        calc_g1_element(du_norm, dx_error, i, j, mu_extra);
      }
    }

    for(int8_t j=0; j<4; j++) {
      G2wls[j] = G2_new[j];
      for(int8_t i=0; i<3; i++) {
        G1wls[i][j] = G1_new[i][j];
      }
    }
 // }
}

/* void calc_g1g2_pseudo_inv(void) {

  //sum of G1 and G2
  for(int8_t i=0; i<3; i++) {
    for(int8_t j=0; j<4; j++) {
      if(i<2)
        G1G2[i][j] = G1[i][j]/1000.0;
      else
        G1G2[i][j] = G1[i][j]/1000.0 + G2[j]/1000.0;
    }
  }

  //G1G2*transpose(G1G2)
  //calculate matrix multiplication of its transpose 3x4 x 4x3
  float element = 0;
  for(int8_t row=0; row<3; row++) {
    for(int8_t col=0; col<3; col++) {
      element = 0;
      for(int8_t i=0; i<4; i++) {
        element = element + G1G2[row][i]*G1G2[col][i];
      }
      MAT33_ELMT(G1G2_trans_mult,row,col) = element;
    }
  }

  //there are numerical errors if the scaling is not right.
  MAT33_MULT_SCALAR(G1G2_trans_mult,100.0);

  //inverse of 3x3 matrix
  MAT33_INV(G1G2inv,G1G2_trans_mult);

  //scale back
  MAT33_MULT_SCALAR(G1G2inv,100.0);

  //G1G2'*G1G2inv
  //calculate matrix multiplication 4x3 x 3x3
  for(int8_t row=0; row<4; row++) {
    for(int8_t col=0; col<3; col++) {
      element = 0;
      for(int8_t i=0; i<3; i++) {
        element = element + G1G2[i][row]*MAT33_ELMT(G1G2inv,col,i);
      }
      G1G2_pseudo_inv[row][col] = element;
    }
  }
} */


// ===========================================================================================


void stabilization_indi_run(bool enable_integrator __attribute__((unused)), bool rate_control)
{
  /* attitude error                          */
  struct Int32Quat att_err;
  struct Int32Quat *att_quat = stateGetNedToBodyQuat_i();
  int32_quat_inv_comp(&att_err, att_quat, &stab_att_sp_quat);
  /* wrap it in the shortest direction       */
  int32_quat_wrap_shortest(&att_err);
  int32_quat_normalize(&att_err);

  /* compute the INDI command */
  stabilization_indi_calc_cmd(stabilization_att_indi_cmd, &att_err, rate_control);

  /* copy the INDI command */
  stabilization_cmd[COMMAND_ROLL] =  stabilization_att_indi_cmd[COMMAND_ROLL];
  stabilization_cmd[COMMAND_PITCH] = stabilization_att_indi_cmd[COMMAND_PITCH];
  stabilization_cmd[COMMAND_YAW] = stabilization_att_indi_cmd[COMMAND_YAW];

 //  WLS control allocator functions
   stabilization_cmd[COMMAND_WLS_1] =  stabilization_att_indi_cmd[COMMAND_WLS_1]; 
   stabilization_cmd[COMMAND_WLS_2] =  stabilization_att_indi_cmd[COMMAND_WLS_2];
   stabilization_cmd[COMMAND_WLS_3] =  stabilization_att_indi_cmd[COMMAND_WLS_3];  
   stabilization_cmd[COMMAND_WLS_4] =  stabilization_att_indi_cmd[COMMAND_WLS_4];
 //  stabilization_cmd[COMMAND_WLS_4] =  stabilization_att_indi_cmd[COMMAND_WLS_3]; 

  //FIXME: Jerryrig to pass vars to motor_mixing.c
  in_cmd_wls[0]  =  stabilization_att_indi_cmd[COMMAND_WLS_1];
  in_cmd_wls[1]  =  stabilization_att_indi_cmd[COMMAND_WLS_2];
  in_cmd_wls[2]  =  stabilization_att_indi_cmd[COMMAND_WLS_3];
  in_cmd_wls[3]  =  stabilization_att_indi_cmd[COMMAND_WLS_4];

  /* bound the result */
  BoundAbs(stabilization_cmd[COMMAND_ROLL], MAX_PPRZ);
  BoundAbs(stabilization_cmd[COMMAND_PITCH], MAX_PPRZ);
  BoundAbs(stabilization_cmd[COMMAND_YAW], MAX_PPRZ); 

}

// This function reads rc commands
void stabilization_indi_read_rc(bool in_flight, bool in_carefree, bool coordinated_turn)
{
  struct FloatQuat q_sp;
#if USE_EARTH_BOUND_RC_SETPOINT
  stabilization_attitude_read_rc_setpoint_quat_earth_bound_f(&q_sp, in_flight, in_carefree, coordinated_turn);
#else
  stabilization_attitude_read_rc_setpoint_quat_f(&q_sp, in_flight, in_carefree, coordinated_turn);
#endif
  QUAT_BFP_OF_REAL(stab_att_sp_quat, q_sp);
}

// Initialize a second order low pass filter
static void stabilization_indi_second_order_filter_init(struct IndiFilter *filter, float omega, float zeta, float omega_r)
{
  filter->omega = omega;
  filter->omega2 = omega * omega;
  filter->zeta = zeta;
  filter->omega_r = omega_r;
  filter->omega2_r = omega_r * omega_r;
}

// This is a simple second order low pass filter
static void stabilization_indi_second_order_filter(struct IndiFilter *filter, struct FloatRates *input)
{
  float_rates_integrate_fi(&filter->x, &filter->dx, 1.0 / PERIODIC_FREQUENCY);
  float_rates_integrate_fi(&filter->dx, &filter->ddx, 1.0 / PERIODIC_FREQUENCY);

  filter->ddx.p = -filter->dx.p * 2 * filter->zeta * filter->omega   + (input->p - filter->x.p) * filter->omega2;
  filter->ddx.q = -filter->dx.q * 2 * filter->zeta * filter->omega   + (input->q - filter->x.q) * filter->omega2;
  filter->ddx.r = -filter->dx.r * 2 * filter->zeta * filter->omega_r + (input->r - filter->x.r) * filter->omega2_r;
}

// This is a Least Mean Squares adaptive filter
// It estiamtes the actuator effectiveness online by comparing the expected angular acceleration based on the inputs with the measured angular acceleration
static inline void lms_estimation(void)
{
  static struct IndiEstimation *est = &indi.est;
  // Only pass really low frequencies so you don't adapt to noise
  stabilization_indi_second_order_filter(&est->u, &indi.u_act_dyn);
  struct FloatRates *body_rates = stateGetBodyRates_f();
  stabilization_indi_second_order_filter(&est->rate, body_rates);

  // The inputs are scaled in order to avoid overflows
  float du = est->u.dx.p * INDI_EST_SCALE;
  est->g1.p = est->g1.p - (est->g1.p * du - est->rate.ddx.p) * du * est->mu;
  du = est->u.dx.q * INDI_EST_SCALE;
  est->g1.q = est->g1.q - (est->g1.q * du - est->rate.ddx.q) * du * est->mu;
  du = est->u.dx.r * INDI_EST_SCALE;
  float ddu = est->u.ddx.r * INDI_EST_SCALE / PERIODIC_FREQUENCY;
  float error = (est->g1.r * du + est->g2 * ddu - est->rate.ddx.r);
  est->g1.r = est->g1.r - error * du * est->mu / 3;
  est->g2 = est->g2 - error * 1000 * ddu * est->mu / 3;

  //the g values should be larger than zero, otherwise there is positive feedback, the command will go to max and there is nothing to learn anymore...
  if (est->g1.p < 0.01) { est->g1.p = 0.01; }
  if (est->g1.q < 0.01) { est->g1.q = 0.01; }
  if (est->g1.r < 0.01) { est->g1.r = 0.01; }
  if (est->g2   < 0.01) { est->g2 = 0.01; }

  if (indi.adaptive) {
    //Commit the estimated G values and apply the scaling
    indi.g1.p = est->g1.p * INDI_EST_SCALE;
    indi.g1.q = est->g1.q * INDI_EST_SCALE;
    indi.g1.r = est->g1.r * INDI_EST_SCALE;
    indi.g2   = est->g2 * INDI_EST_SCALE;
  }
}
