#include"prepare.h"
#include<fftw3.h>


void Calculate_Temp(int *tlist, double *vlist,int nfac,int nvert,double* angles,double* E0,double t0,double Gamma, double A,double R,int N,double *Tres)
{
  double PHI,ep,sigma,omega,beta,lambda,dphi,offset;
  double M[3][3],dMb[3][3],dMo[3][3],dMl[3][3];
  double *mu0,*normal,*CE0,*u0;
  double res;
  double complex *dn;
  double complex tempres=0;
  double d0;
  double T0;
  double thetan;
  double delphi;
  double psin;
  
  int *visible;
 
  dn=(double complex*)calloc(N,sizeof(double complex));
  mu0=(double*)malloc(N*nfac*sizeof(double));
  u0=(double*)calloc(N*nfac,sizeof(double));
  normal=(double*)malloc(3*nfac*sizeof(double));
  visible=(int*)calloc(N*nfac,sizeof(int));
  CE0=(double*)malloc(3*N*sizeof(double));
  
  PHI=1373.0/(R*R);
  ep=0.9;
  sigma=5.670373E-8;
  dphi=2.0*PI/N;
  omega=angles[2];
  beta=angles[0];
  lambda=angles[1];
  offset=omega*t0;
 //FindActualBlockers(tlist,vlist,nfac,nvert,E0,E0,N,visible);
  for(int j=0;j<N;j++)
  {
    rotate(beta,lambda,omega,0.0,j*dphi/omega,M,dMb,dMl,dMo); //Calculate Rotation matrix
    mult_vector(M,E0,&(CE0[3*j])); //Rotated sun direction
  }
  FindActualBlockers(tlist,vlist,nfac,nvert,CE0,CE0,N,visible); //We need only visibility
   
  Calculate_Normals(tlist,vlist,nfac,nvert,normal); //normal nfac x 3 matrix
 
  real_matrix_multiplyT(normal,CE0,nfac,3,N,u0); //normal is nfac x 3, CE0 is Nx3, we calculate normal*CE0'
  for(int k=0;k<N;k++)
       printf("u0: %f\n",u0[k]);
  //result is u0=normal*CE0' , nfacxN matrix
  real_matrix_multiplyT_ele(u0,visible,nfac,N,mu0);
 
  free(u0);

  fftw_plan p;
  
  for(int facet=0;facet<nfac;facet++)
  {
    //Remember to divide the fft with N!!!!
    //Calculate fft for current facet
     p = fftw_plan_dft_r2c_1d(N,&(mu0[facet*N]),dn, FFTW_ESTIMATE);
     fftw_execute(p);
     
       //printf("dn: %f %f\n",creal(dn[k]),cimag(dn[k]));
     //dn contains N Fourier coefficients, we only need first N/2
     d0=creal(dn[0])/N; //This should be real in any case
    
     T0=pow((1-A)*PHI*d0/(ep*sigma),0.25);
      tempres=0;
     for(int i=0;i<N/2;i++)
     {
       thetan=Gamma/(4*ep*sigma*pow(T0,3))*sqrt(0.5*i*omega/(24.0*3600.0)); //is sqrt double? check this!
       delphi=atan(thetan/(thetan+1.0)); //is it ok to use atan?
       psin=pow((1+2*thetan+2*thetan*thetan),-0.5);
    
       
     if(i==0)   
       res=psin*d0;
    else
    {
    
       tempres+=psin*(dn[i])/N*cexp(I*(i*offset-delphi));
    }
  }
     res+=2.0*creal(tempres);
     //check for some pathological cases
     if(d0<1e-16)
       res=0;
     res=fabs(res);
     Tres[facet]=pow(((1.0-A)*PHI/(ep*sigma)*res),0.25);
  
  }
 
  free(dn);
  free(mu0);
  free(normal);
  free(visible);
  free(CE0);
   fftw_destroy_plan(p);
}
     
 