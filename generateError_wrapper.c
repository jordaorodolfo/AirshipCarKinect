

/*
 * Include Files
 *
 */
#if defined(MATLAB_MEX_FILE)
#include "tmwtypes.h"
#include "simstruc_types.h"
#else
#include "rtwtypes.h"
#endif

/* %%%-SFUNWIZ_wrapper_includes_Changes_BEGIN --- EDIT HERE TO _END */
#include <math.h>
/* %%%-SFUNWIZ_wrapper_includes_Changes_END --- EDIT HERE TO _BEGIN */
#define u_width 1
#define y_width 1
/*
 * Create external references here.  
 *
 */
/* %%%-SFUNWIZ_wrapper_externs_Changes_BEGIN --- EDIT HERE TO _END */
/* extern double func(double a); */
/* %%%-SFUNWIZ_wrapper_externs_Changes_END --- EDIT HERE TO _BEGIN */

/*
 * Output functions
 *
 */
void generateError_Outputs_wrapper(real_T *rho,
			real_T *alpha,
			real_T *beta,
			real_T *flightMode,
			const real_T *xD)
{
/* %%%-SFUNWIZ_wrapper_Outputs_Changes_BEGIN --- EDIT HERE TO _END */
// The outputs are equal to the discrete states
rho[0]        = xD[0];   
alpha[0]      = xD[1];
beta[0]       = xD[2];
flightMode[0] = xD[3];
/* %%%-SFUNWIZ_wrapper_Outputs_Changes_END --- EDIT HERE TO _BEGIN */
}

/*
  * Updates function
  *
  */
void generateError_Update_wrapper(const real_T *y0,
			const real_T *x0,
			const real_T *theta0,
			const real_T *wind,
			const real_T *objective,
			const real_T *angOffset,
			real_T *rho,
			real_T *alpha,
			real_T *beta,
			real_T *flightMode,
			real_T *xD)
{
  /* %%%-SFUNWIZ_wrapper_Update_Changes_BEGIN --- EDIT HERE TO _END */
#define PI 3.14159265
float errorX, errorY, errorTheta;
float refX, refY, refTheta, flightFlag;

refX       = objective[1];     // Grabs destination points
refY       = objective[0];
refTheta   = objective[2];
flightFlag = objective[3];

xD[3] = 0;

if(flightFlag == 2){        // Defines hover mode
    refTheta = wind[0];
    xD[3]    = 1;
}
else if(flightFlag == 3){   // Defines follow mode
    refX = refX + 5*sin(angOffset[0] + refTheta);
    refY = refY + 5*cos(angOffset[0] + refTheta);
    xD[3] = 2;
}

errorX     = refX - x0[0];
errorY     = refY - y0[0];
errorTheta = atan2(errorX, errorY);

if((theta0[0] >= PI/4) && (errorTheta < -PI/4)){ // Fix for alpha error
    errorTheta = errorTheta + PI*2;
}
else if((theta0[0] < -PI/4) && (errorTheta >= PI/4)){
    errorTheta = errorTheta - PI*2;
}

if((theta0[0] >= PI/4) && (refTheta < -PI/4)){ // Fix for beta error
    refTheta = refTheta + PI*2;
}
else if((theta0[0] < -PI/4) && (refTheta >= PI/4)){
    refTheta = refTheta - PI*2;
}

xD[0] = sqrt(pow(errorX,2) + pow(errorY,2));
xD[1] = -theta0[0] + errorTheta;
xD[2] = -theta0[0] - xD[1] + refTheta;

if(xD[1] > (PI/2)){
    xD[0] = 6;
    xD[1] = 50;
    xD[2] = 0;
}
else if(xD[1] < (-PI/2)){
    xD[0] = 6;
    xD[1] = -50;
    xD[2] = 0;
//     if(xD[1] < (-PI+0.2618)){
//         xD[1] = 50;
//     }
}

if(flightFlag == 2){    // Hovering
    if(xD[0] <= 5){
        xD[0] = 0;
    }
    
//     if(xD[1] < (-PI/2) || xD[1] > (PI/2)){
//         xD[0] = -xD[0];    // Control law-converted error
//         xD[1] = -theta0[0] + atan2(-errorX, -errorY);
//         xD[2] = -theta0[0] - xD[1] + refTheta;
//     }
}

if(flightFlag == 3){    // Following target
    if(xD[0] <= 2){
        xD[0] = 0;
    }    
    
//     if(xD[1] < (-PI/2) || xD[1] > (PI/2)){
//         xD[0] = -xD[0];    // Control law-converted error
//         xD[1] = -theta0[0] + atan2(-errorX, -errorY);
//         xD[2] = -theta0[0] - xD[1] + refTheta;
//     }
}
/* %%%-SFUNWIZ_wrapper_Update_Changes_END --- EDIT HERE TO _BEGIN */
}
