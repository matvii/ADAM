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
     {
       res=psin*d0;
     }
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
void Calculate_Temp_deriv(int *tlist, double *vlist,int nfac,int nvert,double* angles,double* E0,double t0,double Gamma, double A,double R,int N,double *Tres,double *Tresdx,double *Tresdy,double *Tresdz, double *TresdA)
{
  double PHI,ep,sigma,omega,beta,lambda,dphi,offset;
  double M[3][3],dMb[3][3],dMo[3][3],dMl[3][3];
  double *mu0,*normal,*CE0,*u0;
  double *CE0db,*CE0dl,*CE0do;
  double res;
  double complex *dn;
  double complex tempres=0;
  double d0;
  double T0;
  double thetan;
  double delphi;
  double psin;
  double dnodx1[3],dnodx2[3],dnodx3[3];
  double dnody1[3],dnody2[3],dnody3[3];
  double dnodz1[3],dnodz2[3],dnodz3[3];
  double *dmudx1,*dmudx2,*dmudx3;
  double *dmudy1,*dmudy2,*dmudy3;
  double *dmudz1,*dmudz2,*dmudz3;
  double *dmudb,*dmudl,*dmudo;
  double complex *dndx1,*dndy1,*dndz1;
  double complex *dndx2,*dndy2,*dndz2;
  double complex *dndx3,*dndy3,*dndz3;
  double complex *dndb,*dndl,*dndo;
  double complex dtresdx1,dtresdx2,dtresdx3;
  double complex dtresdy1,dtresdy2,dtresdy3;
  double complex dtresdz1,dtresdz2,dtresdz3;
  double complex dtresdb,dtresdl,dtresdo;
  double DT0d0;
  double DthetandT0,Dpsindth,Ddelphidth,Dthetando;
  double Dpsidx1,Dpsidx2,Dpsidx3;
  double Dpsidy1,Dpsidy2,Dpsidy3;
  double Dpsidz1,Dpsidz2,Dpsidz3;
  
  double Ddelphidx1,Ddelphidx2,Ddelphidx3;
  double Ddelphidy1,Ddelphidy2,Ddelphidy3;
  double Ddelphidz1,Ddelphidz2,Ddelphidz3;
  
  double Dpsidb,Dpsidl,Dpsido;
  double Ddelphidb,Ddelphidl,Ddelphido;
   double dresdx1=0.0,dresdx2=0.0,dresdx3=0.0;
      double dresdy1=0.0,dresdy2=0.0,dresdy3=0.0;
      double dresdz1=0.0,dresdz2=0.0,dresdz3=0.0;
      double dresdb=0.0,dresdl=0.0,dresdo=0.0;
  int t1,t2,t3;
  double *v1,*v2,*v3;
  
  int *visible;
  
  dn=(double complex*)calloc(N,sizeof(double complex));
  mu0=(double*)malloc(N*nfac*sizeof(double));
  u0=(double*)calloc(N*nfac,sizeof(double));
  normal=(double*)malloc(3*nfac*sizeof(double));
  visible=(int*)calloc(N*nfac,sizeof(int));
  CE0=(double*)malloc(3*N*sizeof(double));
  CE0db=(double*)malloc(3*N*sizeof(double));
  CE0dl=(double*)malloc(3*N*sizeof(double));
  CE0do=(double*)malloc(3*N*sizeof(double));
  dmudx1=(double*)malloc(nfac*N*sizeof(double));
  dmudx2=(double*)malloc(nfac*N*sizeof(double));
  dmudx3=(double*)malloc(nfac*N*sizeof(double));
  dmudy1=(double*)malloc(nfac*N*sizeof(double));
  dmudy2=(double*)malloc(nfac*N*sizeof(double));
  dmudy3=(double*)malloc(nfac*N*sizeof(double));
  dmudz1=(double*)malloc(nfac*N*sizeof(double));
  dmudz2=(double*)malloc(nfac*N*sizeof(double));
  dmudz3=(double*)malloc(nfac*N*sizeof(double));
  dmudb=(double*)malloc(nfac*N*sizeof(double));
  dmudl=(double*)malloc(nfac*N*sizeof(double));
  dmudo=(double*)malloc(nfac*N*sizeof(double));
  dndx1=(double complex*)malloc(N*sizeof(double complex));
  dndx2=(double complex*)malloc(N*sizeof(double complex));
  dndx3=(double complex*)malloc(N*sizeof(double complex));
  dndy1=(double complex*)malloc(N*sizeof(double complex));
  dndy2=(double complex*)malloc(N*sizeof(double complex));
  dndy3=(double complex*)malloc(N*sizeof(double complex));
  dndz1=(double complex*)malloc(N*sizeof(double complex));
  dndz2=(double complex*)malloc(N*sizeof(double complex));
  dndz3=(double complex*)malloc(N*sizeof(double complex));
  dndb=(double complex*)malloc(N*sizeof(double complex));
  dndl=(double complex*)malloc(N*sizeof(double complex));
  
  
  //Tres=(double*)malloc(nfac*sizeof(double));
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
    mult_vector(dMb,E0,&(CE0db[3*j])); //Derivatives wrt angles
    mult_vector(dMl,E0,&(CE0dl[3*j]));
    mult_vector(dMo,E0,&(CE0do[3*j]));
  }
  FindActualBlockers(tlist,vlist,nfac,nvert,CE0,CE0,N,visible); //We need only visibility
   
  Calculate_Normals(tlist,vlist,nfac,nvert,normal); //normal nfac x 3 matrix
 
  real_matrix_multiplyT(normal,CE0,nfac,3,N,u0); //normal is nfac x 3, CE0 is Nx3, we calculate normal*CE0'
  
  //result is u0=normal*CE0' , nfacxN matrix
  real_matrix_multiplyT_ele(u0,visible,nfac,N,mu0);
 
  free(u0);

  fftw_plan p;
  
  
   //Calculate derivatives
   
    //derivatives wrt angles
    real_matrix_multiplyT(normal,CE0db,nfac,3,N,dmudb); //normal is nfac x 3, CE0db is Nx3, we calculate normal*CE0db'
    real_matrix_multiplyT(normal,CE0dl,nfac,3,N,dmudl);
    real_matrix_multiplyT(normal,CE0do,nfac,3,N,dmudo);
    

   for(int j=0;j<nfac;j++)
   {
     t1=tlist[3*(j)]; //Indexing starts from 0, t1 starts from one
     t2=tlist[3*(j)+1];
     t3=tlist[3*(j)+2];
     v1=&(vlist[3*(t1-1)]);
     v2=&vlist[3*(t2-1)];
     v3=&vlist[3*(t3-1)];
     //Calculate Normal derivatives
     Calculate_Normal_Derivative(v1,v2,v3,dnodx1,dnodx2,dnodx3,dnody1,dnody2,dnody3,dnodz1,dnodz2,dnodz3);
     
    for(int k=0;k<N;k++)
    {
      dmudx1[j*N+k]=dnodx1[0]*CE0[3*k]+dnodx1[1]*CE0[3*k+1]+dnodx1[2]*CE0[3*k+2];
      dmudx2[j*N+k]=dnodx2[0]*CE0[3*k]+dnodx2[1]*CE0[3*k+1]+dnodx2[2]*CE0[3*k+2];
      dmudx3[j*N+k]=dnodx3[0]*CE0[3*k]+dnodx3[1]*CE0[3*k+1]+dnodx3[2]*CE0[3*k+2];
      
      dmudy1[j*N+k]=dnody1[0]*CE0[3*k]+dnody1[1]*CE0[3*k+1]+dnody1[2]*CE0[3*k+2];
      dmudy2[j*N+k]=dnody2[0]*CE0[3*k]+dnody2[1]*CE0[3*k+1]+dnody2[2]*CE0[3*k+2];
      dmudy3[j*N+k]=dnody3[0]*CE0[3*k]+dnody3[1]*CE0[3*k+1]+dnody3[2]*CE0[3*k+2];
      
      dmudz1[j*N+k]=dnodz1[0]*CE0[3*k]+dnodz1[1]*CE0[3*k+1]+dnodz1[2]*CE0[3*k+2];
      dmudz2[j*N+k]=dnodz2[0]*CE0[3*k]+dnodz2[1]*CE0[3*k+1]+dnodz2[2]*CE0[3*k+2];
      dmudz3[j*N+k]=dnodz3[0]*CE0[3*k]+dnodz3[1]*CE0[3*k+1]+dnodz3[2]*CE0[3*k+2];
    }
     
   }
   
   real_matrix_multiplyT_ele(dmudx1,visible,nfac,N,dmudx1);
   real_matrix_multiplyT_ele(dmudx2,visible,nfac,N,dmudx2);
   real_matrix_multiplyT_ele(dmudx3,visible,nfac,N,dmudx3);
   
    real_matrix_multiplyT_ele(dmudy1,visible,nfac,N,dmudy1);
   real_matrix_multiplyT_ele(dmudy2,visible,nfac,N,dmudy2);
   real_matrix_multiplyT_ele(dmudy3,visible,nfac,N,dmudy3);
   
    real_matrix_multiplyT_ele(dmudz1,visible,nfac,N,dmudz1);
   real_matrix_multiplyT_ele(dmudz2,visible,nfac,N,dmudz2);
   real_matrix_multiplyT_ele(dmudz3,visible,nfac,N,dmudz3);
   
    real_matrix_multiplyT_ele(dmudb,visible,nfac,N,dmudb);
   real_matrix_multiplyT_ele(dmudl,visible,nfac,N,dmudl);
   real_matrix_multiplyT_ele(dmudo,visible,nfac,N,dmudo);
 
 for(int facet=0;facet<nfac;facet++)
 {
     t1=tlist[3*(facet)]-1; //Indexing starts from 0, t1 starts from one
     t2=tlist[3*(facet)+1]-1;
     t3=tlist[3*(facet)+2]-1;
     
      p = fftw_plan_dft_r2c_1d(N,&(dmudx1[facet*N]),dndx1, FFTW_ESTIMATE);
     fftw_execute(p);
      p = fftw_plan_dft_r2c_1d(N,&(dmudx2[facet*N]),dndx2, FFTW_ESTIMATE);
     fftw_execute(p);
      p = fftw_plan_dft_r2c_1d(N,&(dmudx3[facet*N]),dndx3, FFTW_ESTIMATE);
     fftw_execute(p);
     
     p = fftw_plan_dft_r2c_1d(N,&(dmudy1[facet*N]),dndy1, FFTW_ESTIMATE);
     fftw_execute(p);
      p = fftw_plan_dft_r2c_1d(N,&(dmudy2[facet*N]),dndy2, FFTW_ESTIMATE);
     fftw_execute(p);
      p = fftw_plan_dft_r2c_1d(N,&(dmudy3[facet*N]),dndy3, FFTW_ESTIMATE);
     fftw_execute(p);
     
     p = fftw_plan_dft_r2c_1d(N,&(dmudz1[facet*N]),dndz1, FFTW_ESTIMATE);
     fftw_execute(p);
      p = fftw_plan_dft_r2c_1d(N,&(dmudz2[facet*N]),dndz2, FFTW_ESTIMATE);
     fftw_execute(p);
      p = fftw_plan_dft_r2c_1d(N,&(dmudz3[facet*N]),dndz3, FFTW_ESTIMATE);
     fftw_execute(p);
     p = fftw_plan_dft_r2c_1d(N,&(dmudb[facet*N]),dndb, FFTW_ESTIMATE);
     fftw_execute(p);
      p = fftw_plan_dft_r2c_1d(N,&(dmudl[facet*N]),dndl, FFTW_ESTIMATE);
     fftw_execute(p);
     
   
	  //Remember to divide the fft with N!!!!
    //Calculate fft for current facet
     p = fftw_plan_dft_r2c_1d(N,&(mu0[facet*N]),dn, FFTW_ESTIMATE);
     fftw_execute(p);
     
     //dn contains N Fourier coefficients, we only need first N/2
     d0=creal(dn[0])/N; //This should be real in any case
    
     T0=pow((1-A)*PHI*d0/(ep*sigma),0.25);
       DT0d0=pow((1-A)*PHI/(ep*sigma),0.25)*pow(d0,-0.75)*0.25;
 
 tempres=0;
      dtresdx1=0;
      dtresdx2=0;
      dtresdx3=0;
      dtresdy1=0;
      dtresdy2=0;
      dtresdy3=0;
      dtresdz1=0;
      dtresdz2=0;
      dtresdz3=0;
      dtresdb=0;
      dtresdl=0;
      dtresdo=0;
      dresdx1=0.0,dresdx2=0.0,dresdx3=0.0;
      dresdy1=0.0,dresdy2=0.0,dresdy3=0.0;
      dresdz1=0.0,dresdz2=0.0,dresdz3=0.0;
      dresdb=0.0,dresdl=0.0,dresdo=0.0;
     for(int i=0;i<N/2;i++)
     {
       thetan=Gamma/(4*ep*sigma*pow(T0,3))*sqrt(0.5*i*omega/(24.0*3600.0)); //is sqrt double? check this!
       delphi=atan(thetan/(thetan+1.0)); //is it ok to use atan?
       psin=pow((1+2*thetan+2*thetan*thetan),-0.5);
       
       DthetandT0=-3*pow(T0,-4)*Gamma/(4.0*ep*sigma)*sqrt(0.5*i*omega/(24.0*3600.0));
       
       Dpsindth=-0.5*pow((1+2*thetan+2*thetan*thetan),-1.5)*(2+4*thetan);
       Ddelphidth=1.0/(pow(thetan+1,2)+pow(thetan,2));
       
       Dthetando=Gamma/(4*ep*sigma*pow(T0,3))*pow(0.5*i*omega/(24*3600),-0.5)*0.5*0.5*i/(24.0*3600.0);
       
       
      Dpsidx1=Dpsindth*DthetandT0*DT0d0*creal(dndx1[0])/N;
      Dpsidx2=Dpsindth*DthetandT0*DT0d0*creal(dndx2[0])/N;
      Dpsidx3=Dpsindth*DthetandT0*DT0d0*creal(dndx3[0])/N;
      
      Dpsidy1=Dpsindth*DthetandT0*DT0d0*creal(dndy1[0])/N;
      Dpsidy2=Dpsindth*DthetandT0*DT0d0*creal(dndy2[0])/N;
      Dpsidy3=Dpsindth*DthetandT0*DT0d0*creal(dndy3[0])/N;
      
      Dpsidz1=Dpsindth*DthetandT0*DT0d0*creal(dndz1[0])/N;
      Dpsidz2=Dpsindth*DthetandT0*DT0d0*creal(dndz2[0])/N;
      Dpsidz3=Dpsindth*DthetandT0*DT0d0*creal(dndz3[0])/N;
      
      Ddelphidx1=Ddelphidth*DthetandT0*DT0d0*creal(dndx1[0])/N;
      Ddelphidx2=Ddelphidth*DthetandT0*DT0d0*creal(dndx2[0])/N;
      Ddelphidx3=Ddelphidth*DthetandT0*DT0d0*creal(dndx3[0])/N;
      
      Ddelphidy1=Ddelphidth*DthetandT0*DT0d0*creal(dndy1[0])/N;
      Ddelphidy2=Ddelphidth*DthetandT0*DT0d0*creal(dndy2[0])/N;
      Ddelphidy3=Ddelphidth*DthetandT0*DT0d0*creal(dndy3[0])/N;
      
      Ddelphidz1=Ddelphidth*DthetandT0*DT0d0*creal(dndz1[0])/N;
      Ddelphidz2=Ddelphidth*DthetandT0*DT0d0*creal(dndz2[0])/N;
      Ddelphidz3=Ddelphidth*DthetandT0*DT0d0*creal(dndz3[0])/N;
      
      Dpsidb=Dpsindth*DthetandT0*DT0d0*creal(dndb[0])/N;
      Dpsidl=Dpsindth*DthetandT0*DT0d0*creal(dndl[0])/N;
      Dpsido=Dpsindth*Dthetando;
      
      Ddelphidb=Ddelphidth*DthetandT0*DT0d0*creal(dndb[0])/N;
      Ddelphidl=Ddelphidth*DthetandT0*DT0d0*creal(dndl[0])/N;
      Ddelphido=Ddelphidth*Dthetando;
      if(i==0)
      {
	Dpsido=0;
	Ddelphido=0;
	dresdx1=Dpsidx1*d0+psin*creal(dndx1[0])/N;
	dresdx2=Dpsidx2*d0+psin*creal(dndx2[0])/N;
	dresdx3=Dpsidx3*d0+psin*creal(dndx3[0])/N;
	dresdy1=Dpsidy1*d0+psin*creal(dndy1[0])/N;
	dresdy2=Dpsidy2*d0+psin*creal(dndy2[0])/N;
	dresdy3=Dpsidy3*d0+psin*creal(dndy3[0])/N;
	dresdz1=Dpsidz1*d0+psin*creal(dndz1[0])/N;
	dresdz2=Dpsidz2*d0+psin*creal(dndz2[0])/N;
	dresdz3=Dpsidz3*d0+psin*creal(dndz3[0])/N;
	dresdb=psin*creal(dndb[0]/N)+Dpsidb*d0;
	dresdl=psin*creal(dndl[0]/N)+Dpsidl*d0;
	dresdo=0;
	res=psin*d0;
      }
	else
	{
      tempres+=psin*(dn[i])/N*cexp(I*(i*offset-delphi));
	dtresdx1+=(Dpsidx1*dn[i]/N+psin*dndx1[i]/N+psin*(dn[i])/N*(-I*Ddelphidx1))*cexp(I*(i*offset-delphi));
	
	dtresdx2+=(Dpsidx2*dn[i]/N+psin*dndx2[i]/N+psin*(dn[i])/N*(-I*Ddelphidx2))*cexp(I*(i*offset-delphi));
	dtresdx3+=(Dpsidx3*dn[i]/N+psin*dndx3[i]/N+psin*(dn[i])/N*(-I*Ddelphidx3))*cexp(I*(i*offset-delphi));
	
	dtresdy1+=(Dpsidy1*dn[i]/N+psin*dndy1[i]/N+psin*(dn[i])/N*(-I*Ddelphidy1))*cexp(I*(i*offset-delphi));
	dtresdy2+=(Dpsidy2*dn[i]/N+psin*dndy2[i]/N+psin*(dn[i])/N*(-I*Ddelphidy2))*cexp(I*(i*offset-delphi));
	dtresdy3+=(Dpsidy3*dn[i]/N+psin*dndy3[i]/N+psin*(dn[i])/N*(-I*Ddelphidy3))*cexp(I*(i*offset-delphi));
	
	dtresdz1+=(Dpsidz1*dn[i]/N+psin*dndz1[i]/N+psin*(dn[i])/N*(-I*Ddelphidz1))*cexp(I*(i*offset-delphi));
	dtresdz2+=(Dpsidz2*dn[i]/N+psin*dndz2[i]/N+psin*(dn[i])/N*(-I*Ddelphidz2))*cexp(I*(i*offset-delphi));
	dtresdz3+=(Dpsidz3*dn[i]/N+psin*dndz3[i]/N+psin*(dn[i])/N*(-I*Ddelphidz3))*cexp(I*(i*offset-delphi));
	
	dtresdb+=(Dpsidb*dn[i]/N+psin*dndb[i]/N+psin*dn[i]/N*(-I*Ddelphidb))*cexp(I*(i*offset-delphi));
	dtresdl+=(Dpsidl*dn[i]/N+psin*dndl[i]/N+psin*dn[i]/N*(-I*Ddelphidl))*cexp(I*(i*offset-delphi));
	dtresdo+=(Dpsido*dn[i]/N+psin*dn[i]/N*(-I*Ddelphido+I*i*t0))*cexp(I*(i*offset-delphi));
	
	}
    }
    // Dthetando[0]=1;
       


     
     
   
     
 
    
     dresdx1+=2*creal(dtresdx1);
     dresdx2+=2*creal(dtresdx2);
     dresdx3+=2*creal(dtresdx3);
     dresdy1+=2*creal(dtresdy1);
     dresdy2+=2*creal(dtresdy2);
     dresdy3+=2*creal(dtresdy3);
     dresdz1+=2*creal(dtresdz1);
     dresdz2+=2*creal(dtresdz2);
     dresdz3+=2*creal(dtresdz3);
     dresdb+=2*creal(dtresdb);
     dresdl+=2*creal(dtresdl);
     dresdo+=2*creal(dtresdo);
    
     res+=2.0*creal(tempres);
     //check for some pathological cases
     if(d0<1e-16)
       res=0;
     res=fabs(res);
     Tres[facet]=pow(((1.0-A)*PHI/(ep*sigma)*res),0.25);
     Tresdx[nvert*facet+t1]=pow((1.0-A)*PHI/(ep*sigma),0.25)*0.25*pow(res,-0.75)*dresdx1;
     Tresdx[nvert*facet+t2]=pow((1.0-A)*PHI/(ep*sigma),0.25)*0.25*pow(res,-0.75)*dresdx2;
     Tresdx[nvert*facet+t3]=pow((1.0-A)*PHI/(ep*sigma),0.25)*0.25*pow(res,-0.75)*dresdx3;
     
     Tresdy[nvert*facet+t1]=pow((1.0-A)*PHI/(ep*sigma),0.25)*0.25*pow(res,-0.75)*dresdy1;
     Tresdy[nvert*facet+t2]=pow((1.0-A)*PHI/(ep*sigma),0.25)*0.25*pow(res,-0.75)*dresdy2;
     Tresdy[nvert*facet+t3]=pow((1.0-A)*PHI/(ep*sigma),0.25)*0.25*pow(res,-0.75)*dresdy3;
     
     Tresdz[nvert*facet+t1]=pow((1.0-A)*PHI/(ep*sigma),0.25)*0.25*pow(res,-0.75)*dresdz1;
     Tresdz[nvert*facet+t2]=pow((1.0-A)*PHI/(ep*sigma),0.25)*0.25*pow(res,-0.75)*dresdz2;
     Tresdz[nvert*facet+t3]=pow((1.0-A)*PHI/(ep*sigma),0.25)*0.25*pow(res,-0.75)*dresdz3;
     
     TresdA[3*facet]=pow((1.0-A)*PHI/(ep*sigma),0.25)*0.25*pow(res,-0.75)*dresdb;
     TresdA[3*facet+1]=pow((1.0-A)*PHI/(ep*sigma),0.25)*0.25*pow(res,-0.75)*dresdl;
     TresdA[3*facet+2]=pow((1.0-A)*PHI/(ep*sigma),0.25)*0.25*pow(res,-0.75)*dresdo;
 }
     fftw_destroy_plan(p);    
     free(CE0);
     free(CE0db);
     free(CE0dl);
     free(CE0do);
     
     free(normal);
     
     free(dn);
     free(mu0);
     free(u0);
     free(visible);
     
     free(dmudx1);
     free(dmudx2);
     free(dmudx3);
     free(dmudy1);
     free(dmudy2);
     free(dmudy3);
     free(dmudz1);
     free(dmudz2);
     free(dmudz3);
     free(dmudb);
     free(dmudl);
     free(dmudo);
     free(dndx1);
     free(dndx2);
     free(dndx3);
     free(dndy1);
     free(dndy2);
     free(dndy3);
     free(dndz1);
     free(dndz2);
     free(dndz3);
     free(dndb);
     free(dndl);
     
     
     
     
     
     
     
     }
 