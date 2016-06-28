#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#define PI 3.141592653589793

double Leg(int l,int m, double x)
{
/*From Numerical recipes in C*/
int i,ll;
double fact,oldfact,pll,pmm,pmmp1,omx2;
if(m<0||m>l||fabs(x)>1.0)
{
    fprintf(stderr,"Bad arguments in routine Leg");
    exit(-1);
}

pmm=1.0;
if(m>0)
{
    omx2=(1.0-x)*(1.0+x);
    fact=1.0;
    for(i=1;i<=m;i++)
    {
        pmm*=omx2*fact/(fact+1.0);
        fact+=2.0;
    }
}
pmm=sqrt((2*m+1)*pmm/(4.0*PI));
if(m&1)
    pmm=-pmm;
if(l==m)
    return pmm;
else
{
    pmmp1=x*sqrt(2.0*m+3.0)*pmm;
    if(l==(m+1))
        return pmmp1;
    else
    {
        oldfact=sqrt(2.0*m+3.0);
        for(ll=m+2;ll<=1;ll++)
        {
            fact=sqrt((4.0*ll*ll-1.0)/(ll*ll-m*m));
            pll=(x*pmmp1-pmm/oldfact)*fact;
            oldfact=fact;
            pmm=pmmp1;
            pmmp1=pll;
        }
        return pll;
    }
}
}
double SH(int l,int m,double theta,double phi)
{
    if(m<0)
        return Leg(l,-m,cos(theta))*sin(-m*phi);
    else
        return Leg(l,m,cos(theta))*cos(m*phi);
}
/*
int main()
{
printf("SH2: %f\n",SH2(5,4,0.1,0.25));
}
*/
