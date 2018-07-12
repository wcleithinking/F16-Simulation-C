/*
    The nonlinear dynamic model for the F-16 fighting falcon
    Maintained by Wenchao Lei, e-mail: wcleithinking@gmail.com
*/

#include "math.h"
#include "lofi_F16_AeroData.c"  // the LOFI look-up table file
#include "hifi_F16_AeroData.c"  // the HIFI look-up table file

/*****************************************************************************/
/*                      Function Declaration                                 */
/*****************************************************************************/
void atmos(double,double,double*);  // atmosphere
void accels(double*,double*,double*);
void nlplant(double*,double*);
int fix(double);
int sign(double);

/*****************************************************************************/
/*                      mexFunction used in matlab environment               */
/*****************************************************************************/
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    #define XU prhs[0]  // the first input parameter
    #define XDOTY plhs[0]   // the first output parameter

    int i;
    double *xup, *xdotyp;    // creat two pointers

    if ( mxGetM(XU)==18 && mxGetN(XU)==1 ) {
        /* Call nlplant Function */
        xup = mxGetPr(XU);  // let the pointer xup point to the input XU 
        XDOTY = mxCreatDoubleMatrix(18,1,mxREAL); // initialize the first output parameter
        xdotyp = mxGetPr(XDOTY) // let the pointter xdotyp point to the output XDOTY
        nlplant(xup,xdotyp); // use the nlplant function to update the first output parameter
        /* Debug */
        
    }  
} // end mexFunction

