#include"prepare.h"
#define PI2 9.869604401089358
 double complex FTC(double u,double v,double a,double b,double c,double d); 
 double complex FTC2(double v,double a,double b,double c,double d);
 double complex FTC3(double u,double v,double a,double b,double c,double d);
//Derivatives
 double complex FTCda(double u,double v,double a,double b,double c,double d);
 double complex FTCdb(double u,double v,double a,double b,double c,double d);
 double complex FTCdc(double u,double v,double a,double b,double c,double d);
 double complex FTCdd(double u,double v,double a,double b,double c,double d);
//Derivatives,pathological cases
 

 double complex FTC3da(double u,double v,double a,double b,double c,double d);
 double complex FTC3db(double u,double v,double a,double b,double c,double d);
 double complex FTC3dc(double u,double v,double a,double b,double c,double d);
 double complex FTC3dd(double u,double v,double a,double b,double c,double d);
void Calc_FTC(double *u1,double *v1,int nf,double a,double b,double c,double d,double g,double h,double complex *F)
{
  /*INPUT
   * nf frequency pairs u1,v1
   * Triangle vertices (a,b) (c,d) (g,h)
   * deriv==1 if derivatives are to be calculated
   */
  /*OUTPUT
   * F double complex array
   * Derivatives of F wrt a,b,c,d,g,h
   */
 // double u,v,a,b,c,d,g,h;
  double u,v;
  double eps=1e-14;
 
  //F=(double complex*)mxCalloc(nf,sizeof(double complex));
  
    
  /*
  double *Fr,*Fi;
 
  Fr=(double complex*)mxCalloc(nf,sizeof(complex));
  Fi=(double complex*)mxCalloc(nf,sizeof(complex));
  */
  //Loop over all frequency pairs
 
  for(int j=0;j<nf;j++)
  {
   u=u1[j];
   v=v1[j];
  
   //We should replace ==0 with <eps
   if(fabs((c-a)*u+(d-b)*v)>eps)
     F[j]=FTC(u,v,a,b,c,d);
  
   else
     F[j]=FTC3(u,v,a,b,c,d);
   if(fabs((g-c)*u+(h-d)*v)>eps)
   
     F[j]=F[j]+FTC(u,v,c,d,g,h);
   
   else
     F[j]=F[j]+FTC3(u,v,c,d,g,h);
   if(fabs((a-g)*u+(b-h)*v)>eps)
   
     F[j]=F[j]+FTC(u,v,g,h,a,b);
   
   else
     F[j]=F[j]+FTC3(u,v,g,h,a,b);
  }

/*for(int j=0;j<nf;j++)
{
  Fr[j]=creal(F[j]);
  Fi[j]=cimag(F[j]);
}
*/
  
}
void Calc_FTC_deriv(double *u1,double *v1,int nf,double a,double b,double c,double d,double g,double h,double complex *F,double complex *Fda,double complex *Fdb,double complex *Fdc,double complex *Fdd,double complex *Fdg,double complex *Fdh)
{
  /*INPUT
   * nf frequency pairs u1,v1
   * Triangle vertices (a,b) (c,d) (g,h)
   * deriv==1 if derivatives are to be calculated
   */
  /*OUTPUT
   * F double complex array
   * Derivatives of F wrt a,b,c,d,g,h
   */
  // double u,v,a,b,c,d,g,h;
  double eps=1e-18;
  double u,v;
  // int nf,deriv;
  //double complex *F,*Fda,*Fdb,*Fdc,*Fdd,*Fdg,*Fdh;
  
  
  //F=(double complex*)mxCalloc(nf,sizeof(double complex));
  
  
  /*
   *  double *Fr,*Fi;
   * 
   *  Fr=(double complex*)mxCalloc(nf,sizeof(complex));
   *  Fi=(double complex*)mxCalloc(nf,sizeof(complex));
   */
  //Loop over all frequency pairs
 
   
    /*for(int j=0;j<nf;j++)
     * {
     *  Fr[j]=creal(F[j]);
     *  Fi[j]=cimag(F[j]);
  }
  */
   
  /*
     //Have to calculate the derivative
     Fda=(double complex*)mxCalloc(nf,sizeof(double complex));
     Fdb=(double complex*)mxCalloc(nf,sizeof(double complex));
     Fdc=(double complex*)mxCalloc(nf,sizeof(double complex));
   *  Fdd=(double complex*)mxCalloc(nf,sizeof(double complex));
   *  Fdg=(double complex*)mxCalloc(nf,sizeof(double complex));
   *  Fdh=(double complex*)mxCalloc(nf,sizeof(double complex));  
   *  plhs[1]=mxCreateDoubleMatrix(nf,1,mxCOMPLEX);
   *  plhs[2]=mxCreateDoubleMatrix(nf,1,mxCOMPLEX);
   *  plhs[3]=mxCreateDoubleMatrix(nf,1,mxCOMPLEX);
   *  plhs[4]=mxCreateDoubleMatrix(nf,1,mxCOMPLEX);
   *  plhs[5]=mxCreateDoubleMatrix(nf,1,mxCOMPLEX);
   *  plhs[6]=mxCreateDoubleMatrix(nf,1,mxCOMPLEX);
   *  
   *  double *Fdar,*Fdai,*Fdbr,*Fdbi,*Fdci,*Fdcr,*Fddi,*Fddr,*Fdgi,*Fdgr,*Fdhr,*Fdhi;
   *  
   *  Fr=mxGetPr(plhs[0]);
   *  Fi=mxGetPi(plhs[0]);
   *  Fdar=mxGetPr(plhs[1]);
   *  Fdai=mxGetPi(plhs[1]);
   *  Fdbr=mxGetPr(plhs[2]);
   *  Fdbi=mxGetPi(plhs[2]);
   *  Fdcr=mxGetPr(plhs[3]);
   *  Fdci=mxGetPi(plhs[3]);
   *  Fddr=mxGetPr(plhs[4]);
   *  Fddi=mxGetPi(plhs[4]);
   *  Fdgr=mxGetPr(plhs[5]);
   *  Fdgi=mxGetPi(plhs[5]);
   *  Fdhr=mxGetPr(plhs[6]);
   *  Fdhi=mxGetPi(plhs[6]);
   */  
    // printf("eps:\n");
     //printf("%d %d %d \n",fabs((c-a)*u1[0]+(d-b)*v1[0])>eps,fabs((g-c)*u1[0]+(h-d)*v1[0])>eps,fabs((a-g)*u1[0]+(b-h)*v1[0])>eps);
     for(int j=0;j<nf;j++)
     {
       u=u1[j];
       v=v1[j];
      
//We should replace ==0 with <eps
if(fabs((c-a)*u+(d-b)*v)>eps)
{
F[j]=FTC(u,v,a,b,c,d);
Fda[j]=FTCda(u,v,a,b,c,d);
Fdb[j]=FTCdb(u,v,a,b,c,d);
Fdc[j]=FTCdc(u,v,a,b,c,d);
Fdd[j]=FTCdd(u,v,a,b,c,d);
Fdg[j]=0;
Fdh[j]=0;
//printf("Fdd1: %f %f\n",creal(FTCdd(u,v,a,b,c,d)),cimag(FTCdd(u,v,a,b,c,d)));
}

else
{
F[j]=FTC3(u,v,a,b,c,d);
Fda[j]=FTC3da(u,v,a,b,c,d);
Fdb[j]=FTC3db(u,v,a,b,c,d);
Fdc[j]=FTC3dc(u,v,a,b,c,d);
Fdd[j]=FTC3dd(u,v,a,b,c,d);
Fdg[j]=0;
Fdh[j]=0;
}
if(fabs((g-c)*u+(h-d)*v)>eps)
{
F[j]=F[j]+FTC(u,v,c,d,g,h);
Fdc[j]=Fdc[j]+FTCda(u,v,c,d,g,h);
Fdd[j]=Fdd[j]+FTCdb(u,v,c,d,g,h);
Fdg[j]=FTCdc(u,v,c,d,g,h);
Fdh[j]=FTCdd(u,v,c,d,g,h);
//printf("Fdd1: %f %f\n",creal(FTCdb(u,v,c,d,g,h)),cimag(FTCdb(u,v,c,d,g,h)));
}

else
{
F[j]=F[j]+FTC3(u,v,c,d,g,h);
Fdc[j]=Fdc[j]+FTC3da(u,v,c,d,g,h);
Fdd[j]=Fdd[j]+FTC3db(u,v,c,d,g,h);
Fdg[j]=FTC3dc(u,v,c,d,g,h);
Fdh[j]=FTC3dd(u,v,c,d,g,h);
}
if(fabs((a-g)*u+(b-h)*v)>eps)
{
F[j]=F[j]+FTC(u,v,g,h,a,b);
Fda[j]=Fda[j]+FTCdc(u,v,g,h,a,b);
Fdb[j]=Fdb[j]+FTCdd(u,v,g,h,a,b);
Fdg[j]=Fdg[j]+FTCda(u,v,g,h,a,b);
Fdh[j]=Fdh[j]+FTCdb(u,v,g,h,a,b);
}

else
{
F[j]=F[j]+FTC3(u,v,g,h,a,b);
Fda[j]=Fda[j]+FTC3dc(u,v,g,h,a,b);
Fdb[j]=Fdb[j]+FTC3dd(u,v,g,h,a,b);
Fdg[j]=Fdg[j]+FTC3da(u,v,g,h,a,b);
Fdh[j]=Fdh[j]+FTC3db(u,v,g,h,a,b);
}
//printf("Fdd: %f %f\n",creal(Fdd[0]),cimag(Fdd[0]));
}
/*
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
*/
  
}


 inline double complex FTC(double u,double v,double a,double b,double c,double d)
 {
  double complex Eab=cexp(-2.0*PI*I*(a*u+b*v));
  double complex Ecd=cexp(-2.0*PI*I*(c*u+d*v));
  double u2v2=1.0/(pow(u,2)+pow(v,2));
  
  return  1/(PI2*4.0)*((b-d)*u-(a-c)*v)/((a-c)*u+(b-d)*v)*u2v2*(Eab-Ecd);
 }
 
 inline double complex FTC3(double u,double v,double a,double b,double c,double d)
 {
   double u2v2=1.0/(pow(u,2)+pow(v,2));
   return -I/(2*PI)*((b-d)*u-(a-c)*v)*u2v2*cexp(-2.0*PI*I*(a*u+b*v));
 }
 //Derivatives
 inline double complex FTCda(double u,double v,double a,double b,double c,double d)
 {
   double complex Eab=cexp(-2.0*PI*I*(a*u+b*v));
   double complex Ecd=cexp(-2.0*PI*I*(c*u+d*v));
   
   double  cadb=(c-a)*u+(d-b)*v;
   double u2v2=1.0/(pow(u,2)+pow(v,2));
   double  H=-((b-d)*u-(a-c)*v)/cadb*u2v2;
   
   return 1/(4*PI2)*((d-b)/pow(cadb,2)*(Eab-Ecd)+H*(-2.0*PI*I*u*Eab));
 }
 inline double complex FTCdb(double u,double v,double a,double b,double c,double d)
 {
   double complex Eab=cexp(-2.0*PI*I*(a*u+b*v));
   double complex Ecd=cexp(-2.0*PI*I*(c*u+d*v)); 
   double cadb=(c-a)*u+(d-b)*v;
   double u2v2=1.0/(pow(u,2)+pow(v,2));
   double complex H=-((b-d)*u-(a-c)*v)/cadb*u2v2;
   return 1/(4*PI2)*((a-c)/pow(cadb,2)*(Eab-Ecd)+H*(-2.0*PI*I*v*Eab));
 }
 inline double complex FTCdc(double u,double v,double a,double b,double c,double d)
 {
   double complex Eab=cexp(-2.0*PI*I*(a*u+b*v));
   double complex Ecd=cexp(-2.0*PI*I*(c*u+d*v)); 
   double cadb=(c-a)*u+(d-b)*v;
   double u2v2=1.0/(pow(u,2)+pow(v,2));
   double complex H=-((b-d)*u-(a-c)*v)/cadb*u2v2;
   return 1/(4*PI2)*((b-d)/pow(cadb,2)*(Eab-Ecd)+H*(2.0*PI*I*u*Ecd));
 }
 inline double complex FTCdd(double u,double v,double a,double b,double c,double d)
 {
   double complex Eab=cexp(-2.0*PI*I*(a*u+b*v));
   double complex Ecd=cexp(-2.0*PI*I*(c*u+d*v));
   double cadb=(c-a)*u+(d-b)*v;
   double u2v2=1.0/(pow(u,2)+pow(v,2));
   double complex H=-((b-d)*u-(a-c)*v)/cadb*u2v2;
    return 1/(4*PI2)*((c-a)/pow(cadb,2)*(Eab-Ecd)+H*(2.0*PI*I*v*Ecd));
 }
 
 
 inline double complex FTC3da(double u,double v,double a,double b,double c,double d)
 {
   return I/(2*PI)*v*1.0/(pow(u,2)+pow(v,2))*cexp(-2*PI*I*(c*u+d*v));
 }
 inline double complex FTC3db(double u,double v,double a,double b,double c,double d)
 {

   return -I/(2*PI)*u*1.0/(pow(u,2)+pow(v,2))*cexp(-2*PI*I*(c*u+d*v));
 }
 inline double complex FTC3dc(double u,double v,double a,double b,double c,double d)
 {
   double u2v2=1.0/(pow(u,2)+pow(v,2));
   return ((d-b)*pow(u,2)*u2v2+(-I+2*(a-c)*PI*u)*v*1/(2*PI)*u2v2)*cexp(-2*PI*I*(c*u+d*v));
 }
 inline double complex FTC3dd(double u,double v,double a,double b,double c,double d)
 { double u2v2=1.0/(pow(u,2)+pow(v,2));
   return ((a-c)*pow(v,2)*u2v2+(I+2*(d-b)*PI*v)*u*1/(2*PI)*u2v2)*cexp(-2*PI*I*(c*u+d*v));
 }
 
 
   
   
    
  
  