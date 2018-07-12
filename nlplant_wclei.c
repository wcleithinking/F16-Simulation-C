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
    #define XDOT plhs[0]   // the first output parameter

    int i;
    double *xup, *xdotp;    // creat two pointers

    if ( mxGetM(XU)==18 && mxGetN(XU)==1 ) {
        /* Call nlplant Function */
        xup = mxGetPr(XU);  // let the pointer xup point to the input XU 
        XDOT = mxCreatDoubleMatrix(18,1,mxREAL); // initialize the first output parameter
        xdotp = mxGetPr(XDOT) // let the pointter xdotp point to the output XDOT
        nlplant(xup,xdotp); // use the nlplant function to update the first output parameter

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
    double g = 32.17;    // gravity, ft/s^2 ( about 9.8 m/s^2 ) 
    double m = 636.94;   // mass, slugs ( about 636.94*14.593903 kg )
    double B = 30.0;     // span, ft ( about  30*0.3048 m)
    double S = 300.0;    // plantform area, ft^2 ( about 300*0.3048^2 m^2 )
    double cbar = 11.32; // mean aero chord, ft
    double xcgr = 0.35;  // reference center of gravity as a fraction of cbar
    double xcg = 0.30;   // center of gravity as a fraction of cbar
    double Heng = 0.0;   // turbine momentum along roll axis
    double pi = acos(-1);// use the acos function from "math.h" to define the value of pi
    double r2d = 180.0/pi; // radians to degrees
    double d2r = pi/180.0; // degrees to radians

    /* NASA Data ( translated via Equation 2.4--2.6 on Page 80 of Stevens and Lewis ) */
    double Jy = 55814.0;    // slug*ft^2 (about 55814*14.593903*0.3048^2 kg*m^2)
    double Jxz = 982.0;     // slug*ft^2
    double Jz = 63100.0     // slug*ft^2
    double Jx = 9496.0      // slug*ft^2

    double npos, epos, alt; // north position, east position and altitude 
    double phi, theta, psi; // the Euler angles
    double vt, alpha, beta; // total velocity, angle of attack and side-slip angle
    double P, Q, R; // the angular rate
    double sa, ca, sb, cb, tb, st, ct, tt, sphi, cphi, spsi, cpsi;  // the trigonometric functions
    double T, el, ail, rud, dail, drud, lef, dlef;  // control inputs
    double qbar, mach, ps;  // dynamic pressure, Mach number and the static pressure
    double U, V, W, Udot, Vdot, Wdot; // velocities
    double L_tot, M_tot, N_tot, denom;  // moments and the denominator

    double Cx_tot, Cx, delta_Cx_lef, dXdQ, Cxq, delta_Cxq_lef;  // for C_{X_T}
    double Cz_tot, Cz, delta_Cz_lef, dZdQ, Czq, delta_Czq_lef;  // for C_{Z_t}
    double Cm_tot, Cm, eta_el, delta_Cm_lef, dMdQ, Cmq, delta_Cmq_lef, delta_Cm, delta_Cm_ds;   // for C_{m_T}
    double Cy_tot, Cy, delta_Cy_lef, dYdail, delta_Cy_r30, dYdR, dYdP;  // for C_{Y_T}
    double delta_Cy_a20, delta_Cy_a20_lef, Cyr, delta_Cyr_lef, Cyp, delta_Cyp_lef;
    double Cn_tot, Cn, delta_Cn_lef, dNdail, delta_Cn_r30, dNdR, dNdP, delta_Cnbeta;    // for C_{n_T}
    double delta_Cn_a20, delta_Cn_a20_lef, Cnr, delta_Cnr_lef, Cnp, delta_Cnp_lef;
    double Cl_tot, Cl, delta_Cl_lef, dLdail, delta_Cl_r30, dLdR, dLdP, delta_Clbeta;    // for C_{l_T}
    double delta_Cl_a20, delta_Cl_a20_lef, Clr, delta_Clr_lef, Clp, delta_Clp_lef;

    /* A Temp */
    double *temp;                                // a temp pointer
    temp = (double *)malloc(9 * sizeof(double)); // memory allocation for the variable pointed by temp

    /* Get the States */
    npos  = xu[0];
    epos  = xu[1];
    alt   = xu[2];
    phi   = xu[3];  // rad    
    theta = xu[4];  // rad
    psi   = xu[5];  // rad
    vc    = xu[6];
    alpha = xu[7];  // rad
    beta  = xu[8];  // rad
    P     = xu[9];
    Q     = xu[10];
    R     = xu[11];
    sa    = sin(alpha);
    ca    = cos(alpha);
    sb    = sin(beta);
    cb    = cos(beta);
    tb    = tan(beta);
    st    = sin(theta);
    ct    = cos(theta);
    tt    = tan(theta);
    sphi  = sin(phi);
    cphi  = cos(phi);
    spsi  = sin(psi);
    cpsi  = cos(psi);

    if (vt<=0.01) {vt = 0.01;}  // the minimum velocity

    /* Get the Control Inputs */
    T     = xu[12]; // thrust
    el    = xu[13]; // elevator in degrees
    ail   = xu[14]; // ailerons in degrees
    rud   = xu[15]; // rudder in degrees
    lef   = xu[16]; // leading edge flap in degrees
    dail  = ail/21.5;   // max angle 21.5 degrees
    drud = rud/30.0;    // max angle 30.0 degrees
    dlef = (1-lef/25.0);    // max angle 25.0 degrees

    /* Get the FI Flag */
    fi_flag = xu[17];

    /* Get the Atmospheric Parameters*/
    atmos(alt,vt,temp);
    mach = temp[0];
    qbar = temp[1];
    ps   = temp[2]; 

    /*-----------------------------------------------------------------------*/
    /*                          The Dynamics                                 */
    /*-----------------------------------------------------------------------*/
    /* Navigation Equations */
    U = vt*ca*cb;
    V = vt*sb;
    W = vt*sa*cb;
    xdot[0] = U*(ct*cpsi) + V*(sphi*cpsi*st - cphi*spsi) + W*(cphi*st*cpsi + sphi*spsi);    // nposdot
    xdot[1] = U*(ct*spsi) + V*(sphi*spsi*st + cphi*cpsi) + W*(cphi*st*spsi - sphi*cpsi);    // eposdot
    xdot[2] = U*st - V*(sphi*ct) - W*(cphi*ct); // altdot
    /* Kinematic Equations */
    xdot[3] = P + tt*(Q*sphi + R*cphi); // phidot
    xdot[4] = Q*cphi - R*sphi;  // thetadot
    xdot[5] = (Q*sphi + R*cphi)/ct; //psidot


}