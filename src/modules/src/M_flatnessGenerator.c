#include "stm32f4xx.h"

#include "FreeRTOS.h"
#include "param.h"
#include "log.h"
#include "math.h" 

#include "M_flatnessGenerator.h"
#include "arm_math.h"

// Utility ARM functions
static inline void mat_inv(const arm_matrix_instance_f32 * pSrc, arm_matrix_instance_f32 * pDst)
{ configASSERT(ARM_MATH_SUCCESS == arm_mat_inverse_f32(pSrc, pDst)); }
static inline void mat_mult(const arm_matrix_instance_f32 * pSrcA, const arm_matrix_instance_f32 * pSrcB, arm_matrix_instance_f32 * pDst)
{ configASSERT(ARM_MATH_SUCCESS == arm_mat_mult_f32(pSrcA, pSrcB, pDst)); }


static float testOutput[3];

void flatnessExpansion_M(M_setpoint_t *setpoint, M_mode_t *mode)
{
  /****************************************************************************
   * Fills the setpoint structure from the flat outputs providing data only
   * where requested in order to yield better computational efficiency.
   *
   * Current status
   * Position          - Tested and functional
   * Velocity          - Tested and functional
   * Acceleration      - Tested and functional
   * Euler angles      - Tested and functional
   * Quaterinon        - Tested and functional
   * Rotation matrix   - Tested and functional
   * Euler rates       - Tested and functional
   * Body rates        - Tested and functional
   * Body acceleration - Tested and functional (*)
   * Thrust            - Tested and functional
   * Torques           - Tested and functional (*)
   *
   * (*) There are very subtle differences in angular accelerations between the
   *     Simulink and C implementations, these do not matter presently, but
   *     have to be investigated further at some point.
   ***************************************************************************/

  // Loads global position, velocity and acceleration
  setpoint->position.x = setpoint->gamma[0][0];
  setpoint->position.y = setpoint->gamma[1][0];
  setpoint->position.z = setpoint->gamma[2][0];
  
  setpoint->velocity.x = setpoint->gamma[0][1];
  setpoint->velocity.y = setpoint->gamma[1][1];
  setpoint->velocity.z = setpoint->gamma[2][1];

  setpoint->acceleration.x = setpoint->gamma[0][2];
  setpoint->acceleration.y = setpoint->gamma[1][2];
  setpoint->acceleration.z = setpoint->gamma[2][2];
  
  // Computes thrust (T) and the body unit axis (xB,yB,zB)
  float atilde[3] = {setpoint->acceleration.x,
                     setpoint->acceleration.y,
                     setpoint->acceleration.z + Q_GRAVITY};
  float normAtilde = twoNorm(atilde);                        // [m/s^2]
  float T = Q_MASS * normAtilde;                             // [N]
  float zB[3] = {atilde[0]/normAtilde, atilde[1]/normAtilde, atilde[2]/normAtilde};
  float yaw = setpoint->gamma[3][0];
  float yC[3] = {cosf(yaw + PI/2.0f), sinf(yaw + PI/2.0f), 0.0f};
  float xB[3];
  float xBhat[3];
  crossProduct(yC, zB, xBhat);
  normalize(xBhat, xB);
  float yB[3];
  crossProduct(zB, xB, yB);
  setpoint->thrust = T;

  // Computes the rotation matrix (R) with arm instance
  float R[3][3] = {{xB[0],yB[0],zB[0]},
                   {xB[1],yB[1],zB[1]},
                   {xB[2],yB[2],zB[2]}};
  for (int ii = 0; ii < 3; ii++){
    for (int jj = 0; jj < 3; jj++){
      setpoint->rotation[ii][jj] = R[ii][jj];
    }
  }
  arm_matrix_instance_f32 Rm = {3, 3, (float *)R};
  
  if (mode->flatness.eulerAngles == true){
    // R-P-Y Z-Y-X Tait-byrian parametrization of the rotation
    float roll = -asin(R[2][0]);
    float pitch = atan2(R[2][1]/cosf(roll), R[2][2]/cosf(roll));
    setpoint->eulerAngles.roll = roll;
    setpoint->eulerAngles.pitch = pitch;
    setpoint->eulerAngles.yaw = yaw;
  }

  if (mode->flatness.quaternion == true){
    // TODO Quaterinon parametrization of the rotation
    float qsq[4] = {
      0.25f * fabs( R[0][0] + R[1][1] + R[2][2] + 1.0f),
      0.25f * fabs( R[0][0] - R[1][1] - R[2][2] + 1.0f),
      0.25f * fabs(-R[0][0] + R[1][1] - R[2][2] + 1.0f), 
      0.25f * fabs(-R[0][0] - R[1][1] + R[2][2] + 1.0f)
    };
    setpoint->quaternion.qw = sqrtf(qsq[0]);
    setpoint->quaternion.qv[0] = signFunction(R[2][1] - R[1][2]) * sqrtf(qsq[1]);
    setpoint->quaternion.qv[1] = signFunction(R[0][2] - R[2][0]) * sqrtf(qsq[2]);
    setpoint->quaternion.qv[2] = signFunction(R[1][0] - R[0][1]) * sqrtf(qsq[3]);
  }

  mode->flatness.bodyrate = true;
  mode->flatness.bodyacceleration = true;
  mode->flatness.torque = true;
  
  if (mode->flatness.bodyrate == true){ 
    // ~~~ Compute body rates and auler angle rates ~~~
    float zG[3] = {0,0,1};
    float yawrate = setpoint->gamma[3][1];
    float jerk[3] = {setpoint->gamma[0][3], setpoint->gamma[1][3], setpoint->gamma[2][3]};
    float projectedW[3];
    float zBdjerk = dotProduct(zB, jerk);
    float mdT = Q_MASS / T;
    
    // Compute projection
    for (int dim = 0; dim < 3; dim++){
      projectedW[dim] = mdT * (jerk[dim] - zBdjerk * zB[dim]);
    }

    // Set up linear system
    float Wbasis[3][3] = {{xB[0],yC[0],zG[0]},
                          {xB[1],yC[1],zG[1]},
                          {xB[2],yC[2],zG[2]}};
                          
    // Matrix instanciation and multiplicaiton
    arm_matrix_instance_f32 Wbasism = {3, 3, (float *)Wbasis};
    float invWbasis[3][3];
    arm_matrix_instance_f32 invWbasism = {3, 3, (float *)invWbasis};
    float A[3][3];
    arm_matrix_instance_f32 Am = {3, 3, (float *)A};
    mat_inv(&Wbasism, &invWbasism);
    mat_mult(&invWbasism, &Rm, &Am);
    
    // Compute angular rates
    float omegax = -dotProduct(yB, projectedW);
    float omegay = dotProduct(xB, projectedW);
    float omegaz = (yawrate - A[2][0]*omegax - A[2][1]*omegay)/ A[2][2];
    float omega [3] = {omegax, omegay, omegaz};
    setpoint->bodyrate.x = omegax;
    setpoint->bodyrate.y = omegay;
    setpoint->bodyrate.z = omegaz;
    
    float rollrate = A[0][0] * omegax + A[0][1] * omegay + A[0][2] * omegaz;
    float pitchrate = A[1][0] * omegax + A[1][1] * omegay + A[1][2] * omegaz;
    setpoint->eulerRates.roll =rollrate;
    setpoint->eulerRates.pitch = pitchrate;
    setpoint->eulerRates.yaw = yawrate;

    if (mode->flatness.bodyacceleration == true) {
      // ~~~ Compute body accelerations ~~~
      float yawaccel = setpoint->gamma[3][2];
      float snap[3] = {setpoint->gamma[0][4], setpoint->gamma[1][4], setpoint->gamma[2][4]};
      
      float TMzB[3] = {T * zB[0], T * zB[1], T * zB[2]};
      float omegaCTMzB[3];
      crossProduct(omega, TMzB, omegaCTMzB);
      float omegaComegaCTMzB[3];
      crossProduct(omega, omegaCTMzB, omegaComegaCTMzB);
      float v[3] = {Q_MASS * snap[0] - omegaComegaCTMzB[0],
                    Q_MASS * snap[1] - omegaComegaCTMzB[1],
                    Q_MASS * snap[2] - omegaComegaCTMzB[2]};

      float projectedA[3];
      float zBDv = dotProduct(zB, v);
      float mMzBDjerkMzB[3] = {Q_MASS * zBdjerk * zB[0],
                               Q_MASS * zBdjerk * zB[1],
                               Q_MASS * zBdjerk * zB[2]};
      float omegaCmMzBDjerkMzB[3];
      crossProduct(omega, mMzBDjerkMzB, omegaCmMzBDjerkMzB);

      // Compute projection
      for (int ii = 0; ii < 3; ii++){
        projectedA[ii] = ( v[ii] - (zBDv * zB[ii]) - (2 * omegaCmMzBDjerkMzB[0]) ) / T;
      }

      // Set up linear system
      arm_matrix_instance_f32 om = {3, 1, (float *)omega};
      float oc [3][3] = {{0, -omegaz, omegay},{omegaz, 0, -omegax},{-omegay, omegax, 0}};
      arm_matrix_instance_f32 ocm = {3, 3, (float *)oc};
      float ocMR [3][3];
      arm_matrix_instance_f32 ocMRm = {3, 3, (float *)ocMR};
      float ocMRMo [3];
      arm_matrix_instance_f32 ocMRMom = {3, 1, (float *)ocMRMo};
      mat_mult(&ocm, &Rm, &ocMRm);              //  [omega]_cross R
      mat_mult(&ocMRm, &om, &ocMRMom);          // ([omega]_cross R) * omega
      
      float dphixb [3] = {rollrate * xB[0], rollrate * xB[1], rollrate * xB[2]};
      float omegaCdphixb [3];
      crossProduct(omega, dphixb, omegaCdphixb);

      float omegac [3] = {0, 0, omegaz};
      float dthetaMyc [3] = {pitchrate * yC[0], pitchrate * yC[1], pitchrate * yC[2]};
      float omegacCdthetaMyc [3];
      crossProduct(omegac, dthetaMyc, omegacCdthetaMyc);

      float X[3];
      for (int ii = 0; ii < 3; ii++){
        X[ii] = ocMRMo[ii] - omegaCdphixb[ii] - omegacCdthetaMyc[ii];
      }
      
      // Matrix instanciation and multiplicaiton
      arm_matrix_instance_f32 Xm = {3, 1, (float *)X};
      float B[3];
      arm_matrix_instance_f32 Bm = {3, 1, (float *)B};
      mat_mult(&invWbasism, &Xm, &Bm);

      // Compute angular accelerations
      float alphax = -dotProduct(yB, projectedA);
      float alphay = dotProduct(xB, projectedA);
      float alphaz = (yawaccel - B[2] - A[2][0] * alphax - A[2][1] * alphay) / A[2][2];
      setpoint->bodyacceleration.x = alphax;
      setpoint->bodyacceleration.y = alphay;
      setpoint->bodyacceleration.z = alphaz;
    
      if (mode->flatness.torque == true){
        // Compute feedforward torque terms
        float Ialpha [3] = {Q_INERTIA_XX * alphax, Q_INERTIA_YY * alphay, Q_INERTIA_ZZ * alphaz};
        float omegaCIalpha [3];
        crossProduct(omega, Ialpha, omegaCIalpha);
        setpoint->torque.x = Ialpha[0] + omegaCIalpha[0];
        setpoint->torque.y = Ialpha[1] + omegaCIalpha[1];
        setpoint->torque.z = Ialpha[2] + omegaCIalpha[2];

        testOutput[0] = setpoint->torque.x;
        testOutput[1] = setpoint->torque.y;
        testOutput[2] = setpoint->torque.z;
      
      }
    }
  }
}


void crossProduct(float u[], float v[], float res[])
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

float dotProduct(float u[], float v[])
{
  /****************************************************************************
  * Computes the resulting cross product of u and v.
  *
  * ARGS:
  *   u - Vector in R^3 (float [3]).
  *   v - Vector in R^3 (float [3]).
  * RETURNS:
  *   value - The dot product (float)
  ****************************************************************************/
  return v[0] * u[0] + v[1] * u[1] + v[2] * u[2];
}

void normalize(float u[], float res[])
{
  /****************************************************************************
  * Computes the resulting normalized vector of u and in res.
  *
  * ARGS:
  *   u - Vector in R^3 (float [3]).
  *   resulting vector - Vector in R^3 (float [3]).
  ****************************************************************************/
  float nrm = twoNorm(u);
  for (int i = 0; i < 3; i++){res[i] = u[i] / nrm;}
}

float twoNorm(float u[])
{
  /****************************************************************************
  * Computes the two norm of the array u.
  *
  * ARGS:
  *   u - Vector in R^3 (float [3]).
  *
  * RETURNS:
  *   value - Two norm of u
  ****************************************************************************/
  return sqrtf( powf(u[0], 2) + powf(u[1], 2) + powf(u[2], 2) );
}

float signFunction(float u){
  /****************************************************************************
  * Returns the sign of a float
  *
  * ARGS:
  *   u - floating point number
  *
  * RETURNS:
  *   sign - floating point number \in {-1,0,1}
  ****************************************************************************/
  if (u > 0.0f) return 1.0f;
  if (u < 0.0f) return -1.0f;
  return 0.0f;
}

LOG_GROUP_START(testC)
LOG_ADD(LOG_FLOAT, a, &testOutput[0])
LOG_ADD(LOG_FLOAT, b, &testOutput[1])
LOG_ADD(LOG_FLOAT, c, &testOutput[2])
LOG_GROUP_STOP(testC)
