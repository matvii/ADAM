#include"utils.h"
/*
 * This phase function is same as the one in
 * KOALA (by M. Kaasalainen,J. Durech and B. Carry)
 */
double phase_function(double *E,double *E0,double *params,double *dpdp)
{
    /*
     * E view direction
     * E0 sun direction
     * params=[e,c]
     * NOTE: Remember the derivatives wrt rotation angles
     */
    double a=acos(DOT(E,E0));
    double e=exp(-a/params[1]);
    double br=1+params[0]*e+params[2]*a;
    dpdp[0]=e;
    dpdp[1]=params[0]*e*a/pow(params[1],1);
    dpdp[2]=a;
    return br;
}
