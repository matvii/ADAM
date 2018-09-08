#include "utils.h"
/*
 * From Spherical Harmonic Lighting: The Gritty Details
 * by Robin Green
 */
double factorial(int n)
{
double k=1;
for(int i=1;i<=n;i++)
  k*=i;
return k;
}

double P(int l,int m,double x)
{

double pmm = 1.0;
if(m>0) {
double somx2=sqrt((1.0-x)*(1.0+x));
double fact=1.0;
for(int i=1; i<=m; i++) {
pmm*=(-fact)*somx2;
fact+=2.0;
}
}
if(l==m) return pmm;
double pmmp1=x*(2.0*m+1.0)*pmm;
if(l==m+1) return pmmp1;
double pll=0.0;
for(int ll=m+2;ll<=l;++ll) {
pll =((2.0*ll-1.0)*x*pmmp1-(ll+m-1.0)*pmm)/(ll-m);
pmm=pmmp1;
pmmp1=pll;
}
return pll;
}
double K(int l, int m)
{
double temp=((2.0*l+1.0)*factorial(l-m))/(4.0*PI*factorial(l+m));
return sqrt(temp);
}
double SH(int l, int m, double theta, double phi)
{
// return a point sample of a Spherical Harmonic basis function
// l is the band, range [0..N]
// m in the range [-l..l]
// theta in the range [0..Pi]
// phi in the range [0..2*Pi]

double sqrt2=0;
sqrt2=sqrt(2.0);
if(m==0) return K(l,0)*P(l,m,cos(theta));
else if(m>0) return sqrt2*K(l,m)*cos(m*phi)*P(l,m,cos(theta));
else return sqrt2*K(l,-m)*sin(-m*phi)*P(l,-m,cos(theta));
}

