/*
    The nonlinear dynamic model for the F-16 fighting falcon
    Maintained by Wenchao Lei, e-mail: wcleithinking@gmail.com
    Note : This work focus on adding some necessary annotations for the 
    original file "nlplant.c" from the following website:
    https://www.aem.umn.edu/people/faculty/balas/darpa_sec/SEC.Software.html#F16Manual 
*/

#include <math.h>
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
    /*                          AeroDynamics                                 */
    /*-----------------------------------------------------------------------*/
    if (fi_flag == 1) // HIFI Aerodynamic Table
    {
        hifi_C(alpha, beta, el, temp);
        Cx = temp[0];
        Cz = temp[1];
        Cm = temp[2];
        Cy = temp[3];
        Cn = temp[4];
        Cl = temp[5];

        hifi_damping(alpha, temp);
        Cxq = temp[0];
        Cyr = temp[1];
        Cyp = temp[2];
        Czq = temp[3];
        Clr = temp[4];
        Clp = temp[5];
        Cmq = temp[6];
        Cnr = temp[7];
        Cnp = temp[8];

        hifi_C_lef(alpha, beta, temp);
        delta_Cx_lef = temp[0];
        delta_Cz_lef = temp[1];
        delta_Cm_lef = temp[2];
        delta_Cy_lef = temp[3];
        delta_Cn_lef = temp[4];
        delta_Cl_lef = temp[5];

        hifi_damping_lef(alpha, temp);
        delta_Cxq_lef = temp[0];
        delta_Cyr_lef = temp[1];
        delta_Cyp_lef = temp[2];
        delta_Czq_lef = temp[3];
        delta_Clr_lef = temp[4];
        delta_Clp_lef = temp[5];
        delta_Cmq_lef = temp[6];
        delta_Cnr_lef = temp[7];
        delta_Cnp_lef = temp[8];

        hifi_rudder(alpha, beta, temp);
        delta_Cy_r30 = temp[0];
        delta_Cn_r30 = temp[1];
        delta_Cl_r30 = temp[2];

        hifi_ailerons(alpha, beta, temp);
        delta_Cy_a20 = temp[0];
        delta_Cy_a20_lef = temp[1];
        delta_Cn_a20 = temp[2];
        delta_Cn_a20_lef = temp[3];
        delta_Cl_a20 = temp[4];
        delta_Cl_a20_lef = temp[5];

        hifi_other_coeffs(alpha, el, temp);
        delta_Cnbeta = temp[0];
        delta_Clbeta = temp[1];
        delta_Cm = temp[2];
        eta_el = temp[3];
        delta_Cm_ds = 0;   // ignore the deep-stall effect
    }                      // end if
    else if (fi_flag == 0) // LOFI Aerodynamic Table
    {
        dlef = 0.0;

        damping(alpha, temp);
        Cxq = temp[0];
        Cyr = temp[1];
        Cyp = temp[2];
        Czq = temp[3];
        Clr = temp[4];
        Clp = temp[5];
        Cmq = temp[6];
        Cnr = temp[7];
        Cnp = temp[8];

        dmomdcon(alpha, beta, temp);
        delta_Cl_a20 = temp[0];
        delta_Cl_r30 = temp[1];
        delta_Cn_a20 = temp[2];
        delta_Cn_r30 = temp[3];

        clcn(alpha, beta, temp);
        Cl = temp[0];
        Cn = temp[1];

        cxcm(alpha, el, temp);
        Cx = temp[0];
        Cm = temp[1];

        Cy = -0.02 * beta + 0.021 * dail + 0.086 * drud;

        cz(alpha, beta, el, temp);
        Cz = temp[0];

        /* set all higher order terms of hifi that are not applicable to lofi equal to zero. */
        delta_Cx_lef = 0.0;
        delta_Cz_lef = 0.0;
        delta_Cm_lef = 0.0;
        delta_Cy_lef = 0.0;
        delta_Cn_lef = 0.0;
        delta_Cl_lef = 0.0;
        delta_Cxq_lef = 0.0;
        delta_Cyr_lef = 0.0;
        delta_Cyp_lef = 0.0;
        delta_Czq_lef = 0.0;
        delta_Clr_lef = 0.0;
        delta_Clp_lef = 0.0;
        delta_Cmq_lef = 0.0;
        delta_Cnr_lef = 0.0;
        delta_Cnp_lef = 0.0;
        delta_Cy_r30 = 0.0;
        delta_Cy_a20 = 0.0;
        delta_Cy_a20_lef = 0.0;
        delta_Cn_a20_lef = 0.0;
        delta_Cl_a20_lef = 0.0;
        delta_Cnbeta = 0.0;
        delta_Clbeta = 0.0;
        delta_Cm = 0.0;
        eta_el = 1.0; // Needs to be one, see Equation for Cm_tot.
        delta_Cm_ds = 0.0;

    } // end else if

    /* Compute Cx_tot, Cz_tot, Cm_tot, Cy_tot, Cn_tot, and Cl_tot (based on Page 37--40 of NASA report 1538) */
    // Cx_tot
    dXdQ = (cbar / (2 * vt)) * (Cxq + delta_Cxq_lef * dlef);
    Cx_tot = Cx + delta_Cx_lef * dlef + dXdQ * Q;
    // Cz_tot
    dZdQ = (cbar / (2 * vt)) * (Czq + delta_Czq_lef * dlef); // original version: dZdQ = (cbar/(2*vt))*(Czq + delta_Cz_lef*dlef)
    Cz_tot = Cz + delta_Cz_lef * dlef + dZdQ * Q;
    // Cm_tot
    dMdQ = (cbar / (2 * vt)) * (Cmq + delta_Cmq_lef * dlef);
    Cm_tot = Cm * eta_el + Cz_tot * (xcgr - xcg) + delta_Cm_lef * dlef + dMdQ * Q + delta_Cm + delta_Cm_ds;
    // Cy_tot
    dYdail = delta_Cy_a20 + delta_Cy_a20_lef * dlef;
    dYdR = (B / (2 * vt)) * (Cyr + delta_Cyr_lef * dlef);
    dYdP = (B / (2 * vt)) * (Cyp + delta_Cyp_lef * dlef);
    Cy_tot = Cy + delta_Cy_lef * dlef + dYdail * dail + delta_Cy_r30 * drud + dYdR * R + dYdP * P;
    // Cn_tot
    dNdail = delta_Cn_a20 + delta_Cn_a20_lef * dlef;
    dNdR = (B / (2 * vt)) * (Cnr + delta_Cnr_lef * dlef);
    dNdP = (B / (2 * vt)) * (Cnp + delta_Cnp_lef * dlef);
    Cn_tot = Cn + delta_Cn_lef * dlef - Cy_tot * (xcgr - xcg) * (cbar / B) + dNdail * dail + delta_Cn_r30 * drud + dNdR * R + dNdP * P + delta_Cnbeta * beta;
    // Cl_tot
    dLdail = delta_Cl_a20 + delta_Cl_a20_lef * dlef;
    dLdR = (B / (2 * vt)) * (Clr + delta_Clr_lef * dlef);
    dLdP = (B / (2 * vt)) * (Clp + delta_Clp_lef * dlef);
    Cl_tot = Cl + delta_Cl_lef * dlef + dLdail * dail + delta_Cl_r30 * drud + dLdR * R + dLdP * P + delta_Clbeta * beta;

    /*-----------------------------------------------------------------------*/
    /*                              Dynamics                                 */
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

    /* Compute Udot, Vdot, and Wdot (based on Page 36 of NASA report 1538) */
    Udot = R*V - Q*W - g*st + qbar*S*Cx_tot/m + T/m;
    Vdot = P*W - R*U + g*ct*sphi + qbar*S*Cy_tot/m;
    Wdot = Q*U - P*V + g*ct*cphi + qbar*S*Cz_tot/m;
    /* vt_dot */
    xdot[6] = (U*Udot + V*Vdot + W*Wdot)/vt;
    /* alpha_dot */
    xdot[7] = (U*Wdot - W*Udot)/(U*U + W*W);
    /* beta_dot */
    xdot[8] = (Vdot*vt - V*xdot[6])/(vt*vt*cb);

    /* Compute Pdot, Qdot, and Rdot (based on Page 32 of Stevens and Lewis) */
    L_tot = Cl_tot*qbar*S*B;
    M_tot = Cm_tot*qbar*S*cbar;
    N_tot = Cn_tot*qbar*S*B;
    denom = Jx*Jz - Jxz*Jxz;
    xdot[9] = (Jz*L_tot + Jxz*N_tot - (Jz*(Jz-Jy)+Jxz*Jxz)*Q*R + Jxz*(Jx-Jy+Jz)*P*Q + Jxz*Q*Heng)/denom;    // Pdot
    xdot[10] = (M_tot + (Jz-Jx)*P*R - Jxz*(P*P-R*R) - R*Heng)/Jy;   // Qdot
    xdot[11] = (Jx*N_tot + Jxz*L_tot + (Jx*(Jx-Jy)+Jxz*Jxz)*P*Q - Jxz*(Jx-Jy+Jz)*Q*R + Jx*Q*Heng)/denom;    // Rdot

    /*  Create Accelerations */
    accels(xu,xdot,temp);
    xdot[12] = temp[0]; // anx_cg
    xdot[13] = temp[1]; // any_cg 
    xdot[14] = temp[2]; // anz_cg
    
    /* Record the Mach Number, Dynamic Presseure, and Static Pressure */
    xdot[15] = mach;
    xdot[16] = qbar;
    xdot[17] = ps;

    /* free the variable pointed by temp */
    free(temp);

};  // end function: nlplant()