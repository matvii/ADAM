#include"prepare.h"
//void Calculate_HF(int *tlist,double *vlist,int nfac,int nvert,double *angles,double *Eo,double *E0o,double *up,double TIME,double dist,double Gamma,double A,double R,int N,double WL,double *freqx,double *freqy,int nfreq,double *offset,double complex* F);
//void Calculate_Radiance(int *tlist, double *vlist,int nfac,int nvert,double* angles,double *CE,double* CE0,double t0,double Gamma, double A,double R,double L,int N,double *Flux,double *Fldx,double *Fldy,double *Fldz,double *FldA,int deriv);
//void Calculate_Temp(int *tlist, double *vlist,int nfac,int nvert,double* angles,double* E0,double t0,double Gamma, double A,double R,int N,double *Tres);
//void Calculate_Temp_deriv(int *tlist, double *vlist,int nfac,int nvert,double* angles,double* E0,double t0,double Gamma, double A,double R,int N,double *Tres,double *Tdx,double *Tdy,double *Tdz,double *TdA);
void Calculate_HF_deriv(int *tlist,double *vlist,int nfac,int nvert,double *angles,double *Eo,double *E0o,double *up,double TIME,double dist,double Gamma,double A,double Hdist,int N,double WL,double *freqx,double *freqy,int nfreq,double *offset,double complex *F,double complex *dFdx,double complex *dFdy,double complex *dFdz,double complex *dFdA,double complex *dFdoff);
void main()
{
 int tlist[6]={1,2,3,1,4,2};
 double vlist[12]={0,0,100,100,0,0,0,100,0,100,-200,100};
 double freqx[1]={0.1};
 double freqy[1]={0.05};
 double E0[3]={1/sqrt(3),1/sqrt(3),1/sqrt(3)};
 double E[3]={1/sqrt(2),1/sqrt(2),0};
 double up[3]={0,0,1};
 double complex *F;
 double dist=1;
 double angles[3]={0.1,0.2,21.54};
 double offset[2]={0.0,0.0};
 double t0=0.01;
 double Tres1[2],Tres2[2],Tres3;
 double Tdx[8],Tdy[8],Tdz[8],TdA[8];
 double Flux;
 double complex *Fdx,*Fdy,*Fdz,*FdA,*Fdoff;
 int nfac=2;
 int nvert=4;
 int nfreq=1;
 double Gamma=100.0;
 double A=0.1;
 double R=1;
 int N=16;
 
 F=(double complex*)calloc(nfreq,sizeof(double complex));
 Fdx=(double complex*)calloc(nfreq*nvert,sizeof(double complex));
 Fdy=(double complex*)calloc(nfreq*nvert,sizeof(double complex));
 Fdz=(double complex*)calloc(nfreq*nvert,sizeof(double complex));
 FdA=(double complex*)calloc(nfreq*3,sizeof(double complex));
 Fdoff=(double complex*)calloc(nfreq*2,sizeof(double complex));
 
 //Calculate_Temp(tlist,vlist,nfac,nvert,angles,E0,t0,Gamma, A,R,N,Tres2);
 //printf("Temp: %f %f\n",Tres2[0],Tres2[1]);

 //Calculate_Temp_deriv(tlist,vlist,nfac,nvert,angles,E0,t0,Gamma, A,R,N,Tres1,Tdx,Tdy,Tdz,TdA);
// printf("Temp: %f %f\n",Tres1[0],Tres1[1]);
 /*
//printf("\n %f \n",Tres1);
Calculate_Temp_deriv(tlist,vlist,nfac,nvert,angles,E0,t0,Gamma, A,R,N,&Tres2,Tdx,Tdy,Tdz,TdA);
Calculate_Temp(tlist,vlist,nfac,nvert,angles,E0,t0,Gamma, A,R,N,&Tres2);
//printf("\n %f \n",Tres2);
Calculate_Temp_deriv(tlist,vlist,nfac,nvert,angles,E0,t0,Gamma, A,R,N,&Tres3,Tdx,Tdy,Tdz,TdA);
Calculate_Temp(tlist,vlist,nfac,nvert,angles,E0,t0,Gamma, A,R,N,&Tres3);
//printf("\n %f \n",Tres3);
Calculate_Temp_deriv(tlist,vlist,nfac,nvert,angles,E0,t0,Gamma, A,R,N,&Tres2,Tdx,Tdy,Tdz,TdA);
Calculate_Temp(tlist,vlist,nfac,nvert,angles,E0,t0,Gamma, A,R,N,&Tres2);
//printf("\n %f \n",Tres2);
//Tres[0]=0;
*/
 //Calculate_Radiance(tlist,vlist,nfac,nvert,angles,E0,E0,t0,Gamma,A,R,1E-5,N,&Tres1,Tdx,Tdy,Tdz,TdA,1);
 
//printf("\n %.10e \n",Flux);
//Calculate_HF(tlist,vlist,nfac, nvert,angles,E0,E0,up,t0,dist,Gamma,A,R,N,1E-4,freqx,freqy,nfreq,offset,F);
//printf("F: %f %f\n",creal(F[0]),cimag(F[0]));
//F[0]=0.0;
Calculate_HF_deriv(tlist,vlist,nfac, nvert,angles,E,E0,up,t0,dist,Gamma,A,R,N,1E-5,freqx,freqy,nfreq,offset,F,Fdx,Fdy,Fdz,FdA,Fdoff);
//F[0]=0.0;
//Calculate_HF(tlist,vlist,nfac, nvert,angles,E0,E0,up,t0,dist,Gamma,A,R,N,1E-4,freqx,freqy,nfreq,offset,F);

printf("F: %f %f\n",creal(F[0]),cimag(F[0]));

printf("Derivs:\n");
 printf("%.10e  %.10e %.10e %.10e %.10e %.10e %.10e %.10e\n",creal(Fdx[0]),cimag(Fdx[0]),creal(Fdx[1]),cimag(Fdx[1]),creal(Fdx[2]),cimag(Fdx[2]),creal(Fdx[3]),cimag(Fdx[3]));
 printf("%.10e  %.10e %.10e %.10e %.10e %.10e %.10e %.10e\n",creal(Fdy[0]),cimag(Fdy[0]),creal(Fdy[1]),cimag(Fdy[1]),creal(Fdy[2]),cimag(Fdy[2]),creal(Fdy[3]),cimag(Fdy[3]));
 printf("%.10e  %.10e %.10e %.10e %.10e %.10e %.10e %.10e\n",creal(Fdz[0]),cimag(Fdz[0]),creal(Fdz[1]),cimag(Fdz[1]),creal(Fdz[2]),cimag(Fdz[2]),creal(Fdz[3]),cimag(Fdz[3]));
 printf("%.10e  %.10e %.10e %.10e %.10e %.10e\n",creal(FdA[0]),cimag(FdA[0]),creal(FdA[1]),cimag(FdA[1]),creal(FdA[2]),cimag(FdA[2]));
 printf("%.10e  %.10e %.10e  %.10e\n",creal(Fdoff[0]),cimag(Fdoff[0]),creal(Fdoff[1]),cimag(Fdoff[1]));

 
/*
Calculate_Radiance(tlist,vlist,nfac,nvert,angles,E0,E0,t0,Gamma,A,R,1E-5,N,Tres1,Tdx,Tdy,Tdz,TdA,1);
printf("Flux: %.10e %.10e\n",Tres1[0],Tres1[1]);
printf("Derivs:\n");
printf("%.10e %.10e %.10e %.10e\n",Tdx[0],Tdx[1],Tdx[2],Tdx[3]);
printf("%.10e %.10e %.10e %.10e\n",Tdy[0],Tdy[1],Tdy[2],Tdy[3]);
printf("%.10e %.10e %.10e %.10e\n",Tdz[0],Tdz[1],Tdz[2],Tdz[3]);
printf("%.10e %.10e %.10e\n",TdA[0],TdA[1],TdA[2]);

printf("%.10e %.10e %.10e %.10e\n",Tdx[4],Tdx[5],Tdx[6],Tdx[7]);
printf("%.10e %.10e %.10e %.10e\n",Tdy[4],Tdy[5],Tdy[6],Tdy[7]);
printf("%.10e %.10e %.10e %.10e\n",Tdz[4],Tdz[5],Tdz[6],Tdz[7]);
printf("%.10e %.10e %.10e\n",TdA[3],TdA[4],TdA[5]);
*/
  
}


