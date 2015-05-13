#include"mex.h"
#include<math.h>
#include<complex.h>
#define PI 3.141592653589793
#define PI2 9.869604401089358
inline double complex FTC(double u,double v,double a,double b,double c,double d);
inline double complex FTC2(double v,double a,double b,double c,double d);
inline double complex FTC3(double u,double v,double a,double b,double c,double d);
//Derivatives
inline double complex FTCda(double u,double v,double a,double b,double c,double d);
inline double complex FTCdb(double u,double v,double a,double b,double c,double d);
inline double complex FTCdc(double u,double v,double a,double b,double c,double d);
inline double complex FTCdd(double u,double v,double a,double b,double c,double d);
//Derivatives,pathological cases
inline double complex FTC0da(double v,double a,double b,double c,double d);
inline double complex FTC0db(double v,double a,double b,double c,double d);
inline double complex FTC0dc(double v,double a,double b,double c,double d);
inline double complex FTC0dd(double v,double a,double b,double c,double d);

inline double complex FTC3da(double u,double v,double a,double b,double c,double d);
inline double complex FTC3db(double u,double v,double a,double b,double c,double d);
inline double complex FTC3dc(double u,double v,double a,double b,double c,double d);
inline double complex FTC3dd(double u,double v,double a,double b,double c,double d);
void mexFunction(int nlhs, mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
  double *u1,*v1,u,v,a,b,c,d,g,h;
  int nf,deriv;
  double complex *F,*Fda,*Fdb,*Fdc,*Fdd,*Fdg,*Fdh;
  if(nrhs<8)
    mexErrMsgTxt("Input: u v a b c d h g and deriv");
  u1=mxGetPr(prhs[0]);
  v1=mxGetPr(prhs[1]);
  if(mxGetM(prhs[0])!=mxGetM(prhs[1]))
    mexErrMsgTxt("u and v should of same size");
  nf=mxGetM(prhs[0]);
  F=(double complex*)mxCalloc(nf,sizeof(double complex));
  
    
  a=(double)mxGetScalar(prhs[2]);
  b=(double)mxGetScalar(prhs[3]);
  c=(double)mxGetScalar(prhs[4]);
  d=(double)mxGetScalar(prhs[5]);
  g=(double)mxGetScalar(prhs[6]);
  h=(double)mxGetScalar(prhs[7]);
  deriv=(int)mxGetScalar(prhs[8]);
  plhs[0]=mxCreateDoubleMatrix(nf,1,mxCOMPLEX);
  double *Fr,*Fi;
  Fr=mxGetPr(plhs[0]);
  Fi=mxGetPi(plhs[0]);
  //Loop over all frequency pairs
  if(deriv==0)
  {
  for(int j=0;j<nf;j++)
  {
   u=u1[j];
   v=v1[j];
   if(u==0)
   {
     F[j]=FTC2(v,a,b,c,d)+FTC2(v,c,d,g,h)+FTC2(v,g,h,a,b);
     continue;
   }
   //We should replace ==0 with <eps
   if((c-a)*u+(d-b)*v==0)
     F[j]=FTC3(u,v,a,b,c,d);
   else
     F[j]=FTC(u,v,a,b,c,d);
   if((g-c)*u+(h-d)*v==0)
     F[j]=F[j]+FTC3(u,v,c,d,g,h);
   else
     F[j]=F[j]+FTC(u,v,c,d,g,h);
   if((a-g)*u+(b-h)*v==0)
     F[j]=F[j]+FTC3(u,v,g,h,a,b);
   else
     F[j]=F[j]+FTC(u,v,g,h,a,b);
  }

for(int j=0;j<nf;j++)
{
  Fr[j]=creal(F[j]);
  Fi[j]=cimag(F[j]);
}
return;
  }
  //Have to calculate the derivative
  Fda=(double complex*)mxCalloc(nf,sizeof(double complex));
  Fdb=(double complex*)mxCalloc(nf,sizeof(double complex));
  Fdc=(double complex*)mxCalloc(nf,sizeof(double complex));
  Fdd=(double complex*)mxCalloc(nf,sizeof(double complex));
  Fdg=(double complex*)mxCalloc(nf,sizeof(double complex));
  Fdh=(double complex*)mxCalloc(nf,sizeof(double complex));  
  plhs[1]=mxCreateDoubleMatrix(nf,1,mxCOMPLEX);
  plhs[2]=mxCreateDoubleMatrix(nf,1,mxCOMPLEX);
  plhs[3]=mxCreateDoubleMatrix(nf,1,mxCOMPLEX);
  plhs[4]=mxCreateDoubleMatrix(nf,1,mxCOMPLEX);
  plhs[5]=mxCreateDoubleMatrix(nf,1,mxCOMPLEX);
  plhs[6]=mxCreateDoubleMatrix(nf,1,mxCOMPLEX);
  
  double *Fdar,*Fdai,*Fdbr,*Fdbi,*Fdci,*Fdcr,*Fddi,*Fddr,*Fdgi,*Fdgr,*Fdhr,*Fdhi;
  
  Fr=mxGetPr(plhs[0]);
  Fi=mxGetPi(plhs[0]);
  Fdar=mxGetPr(plhs[1]);
  Fdai=mxGetPi(plhs[1]);
  Fdbr=mxGetPr(plhs[2]);
  Fdbi=mxGetPi(plhs[2]);
  Fdcr=mxGetPr(plhs[3]);
  Fdci=mxGetPi(plhs[3]);
  Fddr=mxGetPr(plhs[4]);
  Fddi=mxGetPi(plhs[4]);
  Fdgr=mxGetPr(plhs[5]);
  Fdgi=mxGetPi(plhs[5]);
  Fdhr=mxGetPr(plhs[6]);
  Fdhi=mxGetPi(plhs[6]);
  
  
  for(int j=0;j<nf;j++)
  {
    u=u1[j];
    v=v1[j];
    if(u==0)
    {
      F[j]=FTC2(v,a,b,c,d)+FTC2(v,c,d,g,h)+FTC2(v,g,h,a,b);
      Fda[j]=FTC0da(v,a,b,c,d)+FTC0dc(v,g,h,a,b);
      Fdb[j]=FTC0db(v,a,b,c,d)+FTC0dd(v,g,h,a,b);
      Fdc[j]=FTC0da(v,c,d,g,h)+FTC0dc(v,a,b,c,d);
      Fdd[j]=FTC0db(v,c,d,g,h)+FTC0dd(v,a,b,c,d);
      Fdg[j]=FTC0da(v,g,h,a,b)+FTC0dc(v,c,d,g,h);
      Fdh[j]=FTC0db(v,g,h,a,b)+FTC0dd(v,c,d,g,h);
      continue;
    }
    //We should replace ==0 with <eps
    if((c-a)*u+(d-b)*v==0)
    {
      F[j]=FTC3(u,v,a,b,c,d);
      Fda[j]=FTC3da(u,v,a,b,c,d);
      Fdb[j]=FTC3db(u,v,a,b,c,d);
      Fdc[j]=FTC3dc(u,v,a,b,c,d);
      Fdd[j]=FTC3dd(u,v,a,b,c,d);
      Fdg[j]=0;
      Fdh[j]=0;
      
    }
    
    else
    {
      F[j]=FTC(u,v,a,b,c,d);
      Fda[j]=FTCda(u,v,a,b,c,d);
      Fdb[j]=FTCdb(u,v,a,b,c,d);
      Fdc[j]=FTCdc(u,v,a,b,c,d);
      Fdd[j]=FTCdd(u,v,a,b,c,d);
      Fdg[j]=0;
      Fdh[j]=0;
    }
    if((g-c)*u+(h-d)*v==0)
    {
      F[j]=F[j]+FTC3(u,v,c,d,g,h);
      Fdc[j]=Fdc[j]+FTC3da(u,v,c,d,g,h);
      Fdd[j]=Fdd[j]+FTC3db(u,v,c,d,g,h);
      Fdg[j]=FTC3dc(u,v,c,d,g,h);
      Fdh[j]=FTC3dd(u,v,c,d,g,h);
    }
      
    else
    {
      F[j]=F[j]+FTC(u,v,c,d,g,h);
      Fdc[j]=Fdc[j]+FTCda(u,v,c,d,g,h);
      Fdd[j]=Fdd[j]+FTCdb(u,v,c,d,g,h);
      Fdg[j]=FTCdc(u,v,c,d,g,h);
      Fdh[j]=FTCdd(u,v,c,d,g,h);
    }
    if((a-g)*u+(b-h)*v==0)
    {
      F[j]=F[j]+FTC3(u,v,g,h,a,b);
      Fda[j]=Fda[j]+FTC3dc(u,v,g,h,a,b);
      Fdb[j]=Fdb[j]+FTC3dd(u,v,g,h,a,b);
      Fdg[j]=Fdg[j]+FTC3da(u,v,g,h,a,b);
      Fdh[j]=Fdh[j]+FTC3db(u,v,g,h,a,b);
    }
    
      else
      {
      F[j]=F[j]+FTC(u,v,g,h,a,b);
      Fda[j]=Fda[j]+FTCdc(u,v,g,h,a,b);
      Fdb[j]=Fdb[j]+FTCdd(u,v,g,h,a,b);
      Fdg[j]=Fdg[j]+FTCda(u,v,g,h,a,b);
      Fdh[j]=Fdh[j]+FTCdb(u,v,g,h,a,b);
      }
  }
  for(int j=0;j<nf;j++)
  {
    Fr[j]=creal(F[j]);
    Fi[j]=cimag(F[j]);
    Fdar[j]=creal(Fda[j]);
    Fdai[j]=cimag(Fda[j]);
    Fdbr[j]=creal(Fdb[j]);
    Fdbi[j]=cimag(Fdb[j]);
    Fdcr[j]=creal(Fdc[j]);
    Fdci[j]=cimag(Fdc[j]);
    Fddr[j]=creal(Fdd[j]);
    Fddi[j]=cimag(Fdd[j]);
    Fdgr[j]=creal(Fdg[j]);
    Fdgi[j]=cimag(Fdg[j]);
    Fdhr[j]=creal(Fdh[j]);
    Fdhi[j]=cimag(Fdh[j]);
  } 
  
}


 inline double complex FTC(double u,double v,double a,double b,double c,double d)
 {
  double complex Eab=cexp(-2.0*PI*I*(a*u+b*v));
  double complex Ecd=cexp(-2.0*PI*I*(c*u+d*v));
  return  -(d-b)/(4.0*PI2*u*((c-a)*u+(d-b)*v))*(Ecd-Eab);
 }
 inline double complex FTC2(double v,double a,double b,double c,double d)
 {
   return (c-a)/(4.0*PI2*pow(v,2)*(d-b))*(cexp(-2.0*PI*I*v*d)-cexp(-2.0*PI*I*b*v));
 }
 inline double complex FTC3(double u,double v,double a,double b,double c,double d)
 {
   return -(d-b)/(2.0*PI*I*u)*cexp(-2.0*PI*I*(a*u+b*v));
 }
 //Derivatives
 inline double complex FTCda(double u,double v,double a,double b,double c,double d)
 {
   double complex Eab=cexp(-2.0*PI*I*(a*u+b*v));
   double complex Ecd=cexp(-2.0*PI*I*(c*u+d*v));
   double complex cadb=(c-a)*u+(d-b)*v;
   return -(d-b)/(4.0*PI2*pow(cadb,2))*(Ecd-Eab)-(d-b)/(2.0*PI*(cadb))*(I*Eab);
 }
 inline double complex FTCdb(double u,double v,double a,double b,double c,double d)
 {
   double complex Eab=cexp(-2.0*PI*I*(a*u+b*v));
   double complex Ecd=cexp(-2.0*PI*I*(c*u+d*v)); 
   double complex cadb=(c-a)*u+(d-b)*v;
   return (-Eab+Ecd)/(4.0*PI2*u*cadb)*((b-d)*v/(cadb)+1)+I*(b-d)*Eab*v/(2.0*PI*u*cadb);
 }
 inline double complex FTCdc(double u,double v,double a,double b,double c,double d)
 {
   double complex Eab=cexp(-2.0*PI*I*(a*u+b*v));
   double complex Ecd=cexp(-2.0*PI*I*(c*u+d*v)); 
   double complex cadb=(c-a)*u+(d-b)*v;
   return (d-b)/(4.0*PI2*pow(cadb,2))*(Ecd-Eab)+(d-b)/(2.0*PI*cadb)*(I*Ecd);
 }
 inline double complex FTCdd(double u,double v,double a,double b,double c,double d)
 {
   double complex Eab=cexp(-2.0*PI*I*(a*u+b*v));
   double complex Ecd=cexp(-2.0*PI*I*(c*u+d*v));
   double complex cadb=(c-a)*u+(d-b)*v;
   return -1.0/(4.0*PI2)*(c-a)/pow(cadb,2)*(Ecd-Eab)+(d-b)/(2.0*PI*u*(cadb))*Ecd*I*v;
 }
 inline double complex FTC0da(double v,double a,double b,double c,double d)
 {
   double complex Eab=cexp(-2.0*PI*I*(b*v));
   double complex Ecd=cexp(-2.0*PI*I*(d*v)); 
   return (Eab-Ecd)/(4.0*(d-b)*PI2*pow(v,2));
 }
 inline double complex FTC0db(double v,double a,double b,double c,double d)
 {
   double complex Eab=cexp(-2.0*PI*I*(b*v));
   double complex Ecd=cexp(-2.0*PI*I*(d*v)); 
   return (c-a)*((Ecd-Eab)/pow(((d-b)*2.0*PI*v),2)+I*Eab/(2.0*(d-b)*PI*v));
 }
 inline double complex FTC0dc(double v,double a,double b,double c,double d)
 {
   return -FTC0da(v,a,b,c,d);
 }
 inline double complex FTC0dd(double v,double a,double b,double c,double d)
 {
   double complex Eab=cexp(-2.0*PI*I*(b*v));
   double complex Ecd=cexp(-2.0*PI*I*(d*v)); 
   return (c-a)*(-(Ecd-Eab)/pow((d-b)*2.0*PI*v,2)-I*Ecd/(2.0*(d-b)*PI*v));
 }
 inline double complex FTC3da(double u,double v,double a,double b,double c,double d)
 {
   return -(d-b)*cexp(-2.0*PI*I*(b*v));
 }
 inline double complex FTC3db(double u,double v,double a,double b,double c,double d)
 {
   double complex Eab=cexp(-2.0*PI*I*(b*v));
   return 1.0/(2.0*PI*I*u)*Eab+(d-b)/u*Eab*v;
 }
 inline double complex FTC3dc(double u,double v,double a,double b,double c,double d)
 {
   return 0;
 }
 inline double complex FTC3dd(double u,double v,double a,double b,double c,double d)
 {
   return -1.0/(2*PI*1i*u)*cexp(-2.0*PI*I*(b*v));
 }
 
 
   
   
    
  
  