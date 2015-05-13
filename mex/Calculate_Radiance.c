#include"prepare.h"
void Calculate_Radiance(int *tlist, double *vlist,int nfac,int nvert,double* angles,double *CE,double* CE0,double t0,double Gamma, double A,double R,double L,int N,double *Flux,double *Fldx,double *Fldy,double *Fldz,double *FldA,int deriv)
{
  double beta,lambda,omega;
  double M[3][3],dMb[3][3],dMo[3][3],dMl[3][3];
  double E[3];
  double *normal,*mu;
  int *visible;
  int t1,t2,t3;
  double *v1,*v2,*v3;
  double planck;
  double hc=6.62606957e-34; //Planck's
double cc=299792458; //speed of light
double kc=1.3806488e-23; //Boltzmann
double *Tres;
double nu=cc/L;

double C1=2.0*hc*pow(nu,3.0)/pow(cc,2.0);
double C2=hc*nu/kc;

Tres=(double*)calloc(nfac,sizeof(double));
visible=(int*)calloc(nfac,sizeof(int));
beta=angles[0];
lambda=angles[1];
omega=angles[2];
if(deriv==0)
{
Calculate_Temp(tlist, vlist, nfac,nvert, angles,CE0,t0,Gamma, A, R, N,Tres);
rotate(beta,lambda,omega,0.0,t0,M,dMb,dMl,dMo);
mult_vector(M,CE,E);
 FindActualBlockers(tlist,vlist,nfac,nvert,E,E,1,visible);
 
 

for(int j=0;j<nfac;j++)
{
  
  if(visible[j]==0)
    continue;
 
 
     planck=C1/(exp(C2/Tres[j])-1);
     Flux[j]=planck;
     if(Tres[j]==0)
       Flux[j]=0;
}
free(Tres);
free(visible);
}
else
{
  double *Tdx,*Tdy,*Tdz,*TdA;
  Tdx=(double*)calloc(nfac*nvert,sizeof(double));
  Tdy=(double*)calloc(nfac*nvert,sizeof(double));
  Tdz=(double*)calloc(nfac*nvert,sizeof(double));
  TdA=(double*)calloc(nfac*3,sizeof(double));
  Calculate_Temp_deriv(tlist, vlist, nfac,nvert, angles,CE0,t0,Gamma, A, R, N,Tres,Tdx,Tdy,Tdz,TdA);

rotate(beta,lambda,omega,0.0,t0,M,dMb,dMl,dMo);
mult_vector(M,CE,E);
 FindActualBlockers(tlist,vlist,nfac,nvert,E,E,1,visible);
 
for(int j=0;j<nfac;j++)
{
  if(visible[j]==0)
    continue;
   t1=tlist[3*(j)]-1; //Indexing starts from 0, t1 starts from one
     t2=tlist[3*(j)+1]-1;
     t3=tlist[3*(j)+2]-1;
     

     planck=C1/(exp(C2/Tres[j])-1);
     Flux[j]=planck;
      
     if(Tres[j]==0)
       Flux[j]=0;
     Fldx[nvert*j+t1]=C1*1/pow(exp(C2/Tres[j])-1,2)*exp(C2/(Tres[j]))*C2/pow(Tres[j],2)*Tdx[nvert*j+t1];
     Fldx[nvert*j+t2]=C1*1/pow(exp(C2/Tres[j])-1,2)*exp(C2/(Tres[j]))*C2/pow(Tres[j],2)*Tdx[nvert*j+t2];
     Fldx[nvert*j+t3]=C1*1/pow(exp(C2/Tres[j])-1,2)*exp(C2/(Tres[j]))*C2/pow(Tres[j],2)*Tdx[nvert*j+t3];
     
     Fldy[nvert*j+t1]=C1*1/pow(exp(C2/Tres[j])-1,2)*exp(C2/(Tres[j]))*C2/pow(Tres[j],2)*Tdy[nvert*j+t1];
     Fldy[nvert*j+t2]=C1*1/pow(exp(C2/Tres[j])-1,2)*exp(C2/(Tres[j]))*C2/pow(Tres[j],2)*Tdy[nvert*j+t2];
     Fldy[nvert*j+t3]=C1*1/pow(exp(C2/Tres[j])-1,2)*exp(C2/(Tres[j]))*C2/pow(Tres[j],2)*Tdy[nvert*j+t3];
     
     Fldz[nvert*j+t1]=C1*1/pow(exp(C2/Tres[j])-1,2)*exp(C2/(Tres[j]))*C2/pow(Tres[j],2)*Tdz[nvert*j+t1];
     Fldz[nvert*j+t2]=C1*1/pow(exp(C2/Tres[j])-1,2)*exp(C2/(Tres[j]))*C2/pow(Tres[j],2)*Tdz[nvert*j+t2];
     Fldz[nvert*j+t3]=C1*1/pow(exp(C2/Tres[j])-1,2)*exp(C2/(Tres[j]))*C2/pow(Tres[j],2)*Tdz[nvert*j+t3];
     
     FldA[3*j]=C1*1/pow(exp(C2/Tres[j])-1,2)*exp(C2/(Tres[j]))*C2/pow(Tres[j],2)*TdA[3*j];
     FldA[3*j+1]=C1*1/pow(exp(C2/Tres[j])-1,2)*exp(C2/(Tres[j]))*C2/pow(Tres[j],2)*TdA[3*j+1];
     FldA[3*j+2]=C1*1/pow(exp(C2/Tres[j])-1,2)*exp(C2/(Tres[j]))*C2/pow(Tres[j],2)*TdA[3*j+2];
     if(Flux[j]==0)
     {
	Fldx[nvert*j+t1]=0;
	Fldx[nvert*j+t2]=0;
	Fldx[nvert*j+t3]=0;
	Fldy[nvert*j+t1]=0;
	Fldy[nvert*j+t2]=0;
	Fldy[nvert*j+t3]=0;
	Fldz[nvert*j+t1]=0;
	Fldz[nvert*j+t2]=0;
	Fldz[nvert*j+t3]=0;
	FldA[3*j]=0;
	FldA[3*j+1]=0;
	FldA[3*j+2]=0;
     }
     
}
free(Tdx);
free(Tdy);
free(Tdz);
free(TdA);
free(Tres);
free(visible);
}
}