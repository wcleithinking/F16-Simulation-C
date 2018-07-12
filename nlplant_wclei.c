/*
    The nonlinear dynamic model for the F-16 fighting falcon
    Maintained by Wenchao Lei, e-mail: wcleithinking@gmail.com
    Note : This work focus on adding some necessary annotations for the 
    original file "nlplant.c" from the following website:
    https://www.aem.umn.edu/people/faculty/balas/darpa_sec/SEC.Software.html#F16Manual 
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
/*                  The mexFunction used in Matlab Environment               */
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
        for (i=0;i<=14;i++){
            printf("xdotp(%d) = %e\n",i+1,xdotp[i]);
        }   // end Debug

    }   // end if 
    else {
        mexErrMsgTxt("Input and/or output is wrong size.");
    }   // end else

}   // end mexFunction

/*****************************************************************************/
/*                   The Nonlinear Dynamics of F-16                          */
/*****************************************************************************/
void nlplant(double *xu,double *xdot){

    int fi_flag;    // the flag for choosing LOFI or HIFI

    /* The F-16 Constants */
    double g = 32.17    // gravity, ft/s^2 
    double m = 636.94   // mass, slugs


}