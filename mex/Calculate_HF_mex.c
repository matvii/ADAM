#include "prepare.h"
#include "num_of_threads.h"
#include<omp.h>
void Calculate_HF(int *tlist,double *vlist,int nfac,int nvert,double *angles,double *Eo,double *E0o,double *up,double TIME,double dist,double Gamma,double A,double Hdist,int N,double WL,double *freqx,double *freqy,int nfreq,double *offset,double complex *F);
void Calculate_HF_deriv(int *tlist,double *vlist,int nfac,int nvert,double *angles,double *Eo,double *E0o,double *up,double TIME,double dist,double Gamma,double A,double Hdist,int N,double WL,double *freqx,double *freqy,int nfreq,double *offset,double complex *F,double complex *dFdx,double complex *dFdy,double complex *dFdz,double complex *dFdA,double complex *dFdoff);
//void Calculate_AO(int *tlist,double *vlist,int nfac,int nvert,double *angles,double *Eo,double *E0o,double *up,double TIME,double dist,double *freqx,double *freqy,int nfreq,double *offset,double complex* F);
void Matrix_Multiply(complex double *A, double *B,complex double *C,int n,int m,int l);
//void Calculate_AO_deriv(int *tlist,double *vlist,int nfac,int nvert,double *angles,double *Eo,double *E0o,double *up,double TIME,double dist,double *freqx,double *freqy,int nfreq,double *offset,double complex* F,double complex* dFdx,double complex* dFdy,double complex *dFdz,double complex* dFdA,double complex* dFdoff);
void mexFunction(int nlhs, mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
 /*Same as the original, only exception is the inclusion of matrix D (For effective memory usage)
  */
 //OUTPUT: FT FTdx FTdy FTdz FTdA FTdoff
  if(nrhs!=16)
    mexErrMsgTxt("16 input arguments required:a list of  triangles, vertices, angles,FT.E,FT.E0,FT.TIME,FT.freq,FT.up,FT.dist,FT.Hdist,Gamma,A,N,FT.WL,offset,deriv\n");
  //INPUT
  int nvert,nfac;
  nvert=mxGetM(prhs[1]);
  nfac=mxGetM(prhs[0]);
  if(mxGetN(prhs[1])!=3)
    mexErrMsgTxt("vlist size:nvertx3");
  if(mxGetN(prhs[0])!=3)
    mexErrMsgTxt("tlist size:nvertx3");
  double *vlist,*vl,*tl;
  int *tlist;
  vlist=mxCalloc(3*nvert,sizeof(double));
  tlist=mxCalloc(3*nfac,sizeof(int));
  //NB: MATLAB uses column first ordering
  //We have to change this
  
  tl=mxGetPr(prhs[0]);
  vl=mxGetPr(prhs[1]);
  for(int j=0;j<nvert;j++)
  {
    vlist[3*j]=vl[j];
    vlist[3*j+1]=vl[j+nvert];
    vlist[3*j+2]=vl[j+2*nvert];
  }
  for(int j=0;j<nfac;j++)
  {
    tlist[3*j]=(int)tl[j];
    tlist[3*j+1]=(int)tl[j+nfac];
    tlist[3*j+2]=(int)tl[j+2*nfac];
  }
  
  double *angles;
  angles=mxGetPr(prhs[2]);
  double *E,*E0;
 int Nfft;
  double TIME;
 
  
  int nfreq;
 
  double *offset;
  offset=mxGetPr(prhs[14]);
  
 
  int deriv;
  deriv=(int)mxGetScalar(prhs[15]);
  
  double Gamma;
  double A;
  Gamma=(double)mxGetScalar(prhs[10]);
  A=(double)mxGetScalar(prhs[11]);
  Nfft=(int)mxGetScalar(prhs[12]);
  //Cell handling
  mwSize nobs;
  mwIndex index;
  mxArray *Freq_ptr;
  mxArray *FTEcell,*FTE0cell,*FTTIMEcell,*FTfreqcell,*FTupcell,*FTdistcell,*FTHdistcell,*FTWLcell;
  //double *FTE,*FTTIME,*FTfreq,*FTrfreq;
  int M,N;
  int *nopoints,*cumpoints,ntpoints;
  ntpoints=0;
  //double *freqx,*freqy,*FTfreq;
  FTEcell=prhs[3];
  FTE0cell=prhs[4];
  nobs = mxGetNumberOfElements(FTEcell); 
  
  FTTIMEcell=prhs[5];
  FTfreqcell=prhs[6];
  FTupcell=prhs[7];
  FTdistcell=prhs[8];
  FTHdistcell=prhs[9];
  FTWLcell=prhs[13];
  nopoints=mxCalloc(nobs,sizeof(int));
  cumpoints=mxCalloc(nobs+1,sizeof(int));
  for(int j=0;j<nobs;j++)
  {
   // Freq_ptr = mxGetCell(FTfreqcell, j);
   // FTfreq=mxGetPr(mxGetCell(FTfreqcell,j));
    M=mxGetM(mxGetCell(FTfreqcell, j));
    N=mxGetN(mxGetCell(FTfreqcell, j));
    nopoints[j]=M;
    ntpoints+=M;
    cumpoints[j+1]=cumpoints[j]+nopoints[j];
  // mexPrintf("M: %d N: %d\n",M,N);
  }
 
 //Get data corresponding to one observation
 int obsind=0;
 //FTE=mxGetPr(mxGetCell(FTEcell,obsind));
 //FTTIME=mxGetPr(mxGetCell(FTTIMEcell,obsind));
 //FTfreq=mxGetPr(mxGetCell(FTfreqcell,obsind));
 //FTrfreq=mxGetPr(mxGetCell(FTrfreqcell,obsind));
 //mexPrintf("radar freq: %f\n",*FTrfreq);
  //Allocate output variables
//  FT=(double complex*)mxCalloc(nfreq,sizeof(double complex));
   plhs[0]=mxCreateDoubleMatrix(ntpoints,1,mxCOMPLEX);
   double *FTr,*FTi;
   FTr=mxGetPr(plhs[0]);
    FTi=mxGetPi(plhs[0]);

  if(deriv==0)
  {  
omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
for(obsind=0;obsind<nobs;obsind++)
  {
    double *FTE,*FTE0,*FTTIME,*FTfreq,*FTup,*FTdist,*FTHdist,*FTWL;
    double complex *FT;
   //  obsind=omp_get_thread_num();
    FT=(double complex*)calloc(nopoints[obsind],sizeof(double complex));
   
   
    FTE=mxGetPr(mxGetCell(FTEcell,obsind));
    FTE0=mxGetPr(mxGetCell(FTE0cell,obsind));
    FTup=mxGetPr(mxGetCell(FTupcell,obsind));
    FTTIME=mxGetPr(mxGetCell(FTTIMEcell,obsind));
    FTfreq=mxGetPr(mxGetCell(FTfreqcell,obsind));
    FTdist=mxGetPr(mxGetCell(FTdistcell,obsind));
    FTHdist=mxGetPr(mxGetCell(FTHdistcell,obsind));
    FTWL=mxGetPr(mxGetCell(FTWLcell,obsind));
  //  printf("obsind: %d nopoints: %d cumpoints: %d total: %d radarfreq: %f\n",obsind,nopoints[obsind],cumpoints[obsind],ntpoints,*FTrfreq);
  //printf("TIME: %f freqs: %f %f nfreq: %d offset: %f %f\n",*FTTIME,*FTfreq,*(FTfreq+nopoints[obsind]),nopoints[obsind],*offset,*(offset+1));
   // Calculate_AO(tlist,vlist,nfac,nvert,angles,FTE,FTE0,FTup,*FTTIME,*FTdist,FTfreq,FTfreq+nopoints[obsind],nopoints[obsind],offset+2*obsind,FT);
    Calculate_HF(tlist,vlist,nfac,nvert,angles,FTE,FTE0,FTup,*FTTIME,*FTdist,Gamma,A,*FTHdist,Nfft,*FTWL,FTfreq,FTfreq+nopoints[obsind],nopoints[obsind],offset+2*obsind,FT);
  
   
     for(int j=0;j<nopoints[obsind];j++)
  {
    FTr[j+cumpoints[obsind]]=creal(FT[j]);
    FTi[j+cumpoints[obsind]]=cimag(FT[j]);
  }
  
  //return;
    free(FT);
  }
 
 return; 
}


int nvertf;

  nvertf=nvert;


  
  
  plhs[0]=mxCreateDoubleMatrix(ntpoints,1,mxCOMPLEX);
  plhs[1]=mxCreateDoubleMatrix(ntpoints,nvertf,mxCOMPLEX);
  plhs[2]=mxCreateDoubleMatrix(ntpoints,nvertf,mxCOMPLEX);
  plhs[3]=mxCreateDoubleMatrix(ntpoints,nvertf,mxCOMPLEX);
  plhs[4]=mxCreateDoubleMatrix(ntpoints,3,mxCOMPLEX);
  plhs[5]=mxCreateDoubleMatrix(ntpoints,2*nobs,mxCOMPLEX);
  plhs[6]=mxCreateDoubleMatrix(ntpoints,1,mxCOMPLEX);
  
  //Aux variables to handle complex C->Matlab
  double *FTdxr,*FTdxi,*FTdyr,*FTdyi,*FTdzr,*FTdzi,*FTdAr,*FTdAi,*FTdoffr,*FTdoffi,*FTdexpr,*FTdexpi;
  FTr=mxGetPr(plhs[0]);
  FTi=mxGetPi(plhs[0]);
  FTdxr=mxGetPr(plhs[1]);
  FTdxi=mxGetPi(plhs[1]);
  FTdyr=mxGetPr(plhs[2]);
  FTdyi=mxGetPi(plhs[2]);
  FTdzr=mxGetPr(plhs[3]);
  FTdzi=mxGetPi(plhs[3]);
  FTdAr=mxGetPr(plhs[4]);
  FTdAi=mxGetPi(plhs[4]);
  FTdoffr=mxGetPr(plhs[5]);
  FTdoffi=mxGetPi(plhs[5]);
  FTdexpr=mxGetPr(plhs[6]);
  FTdexpi=mxGetPi(plhs[6]);
  omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
for(obsind=0;obsind<nobs;obsind++)
  {
    int cind=0;
    int oind=0;
    double complex *FTdx,*FTdy,*FTdz,*FTdA,*FTdoff,*FTdxf,*FTdyf,*FTdzf;
     double *FTE,*FTE0,*FTTIME,*FTfreq,*FTdist,*FTup,*FTHdist,*FTWL;
    double complex *FT;
   //  obsind=omp_get_thread_num();
    FT=(double complex*)calloc(nopoints[obsind],sizeof(double complex));
   FTdx=(double complex*)calloc(nopoints[obsind]*nvertf,sizeof(double complex));
  FTdy=(double complex*)calloc(nopoints[obsind]*nvertf,sizeof(double complex));
  FTdz=(double complex*)calloc(nopoints[obsind]*nvertf,sizeof(double complex));
  FTdA=(double complex*)calloc(nopoints[obsind]*3,sizeof(double complex));
  FTdoff=(double complex*)calloc(nopoints[obsind]*2,sizeof(double complex));

  FTdist=mxGetPr(mxGetCell(FTdistcell,obsind));
  FTE0=mxGetPr(mxGetCell(FTE0cell,obsind));
    FTup=mxGetPr(mxGetCell(FTupcell,obsind));
    FTE=mxGetPr(mxGetCell(FTEcell,obsind));
    FTTIME=mxGetPr(mxGetCell(FTTIMEcell,obsind));
    FTfreq=mxGetPr(mxGetCell(FTfreqcell,obsind));
     FTHdist=mxGetPr(mxGetCell(FTHdistcell,obsind));
    FTWL=mxGetPr(mxGetCell(FTWLcell,obsind));
   
      //Calculate_AO_deriv(tlist,vlist,nfac,nvert,angles,FTE,FTE0,FTup, *FTTIME,*FTdist,FTfreq,FTfreq+nopoints[obsind],nopoints[obsind],offset+2*obsind,FT,FTdx,FTdy,FTdz,FTdA,FTdoff);
 Calculate_HF_deriv(tlist,vlist,nfac,nvert,angles,FTE,FTE0,FTup,*FTTIME,*FTdist,Gamma,A,*FTHdist,Nfft,*FTWL,FTfreq,FTfreq+nopoints[obsind],nopoints[obsind],offset+2*obsind,FT,FTdx,FTdy,FTdz,FTdA,FTdoff);
  
  //Copy variables to matlab
  cind=cumpoints[obsind];
  oind=nopoints[obsind];
  for(int j=0;j<nopoints[obsind];j++)
  {
    FTr[j+cind]=creal(FT[j]);
    FTi[j+cind]=cimag(FT[j]);
   
    FTdAr[j+cind]=creal(FTdA[j*3]);
    FTdAi[j+cind]=cimag(FTdA[j*3]);
    FTdAr[j+ntpoints+cind]=creal(FTdA[j*3+1]);
    FTdAi[j+ntpoints+cind]=cimag(FTdA[j*3+1]);
    FTdAr[j+2*ntpoints+cind]=creal(FTdA[j*3+2]);
    FTdAi[j+2*ntpoints+cind]=cimag(FTdA[j*3+2]);
    FTdoffr[j+cind+2*obsind*ntpoints]=creal(FTdoff[j*2]);
    FTdoffi[j+cind+2*obsind*ntpoints]=cimag(FTdoff[j*2]);
    FTdoffr[j+cind+(2*obsind+1)*ntpoints]=creal(FTdoff[j*2+1]);
    FTdoffi[j+cind+(2*obsind+1)*ntpoints]=cimag(FTdoff[j*2+1]);
    //Now the derivatives wrt vertices
    for(int k=0;k<nvertf;k++)
    {
      FTdxr[j+k*ntpoints+cind]=creal(FTdx[j*nvertf+k]);
      FTdxi[j+k*ntpoints+cind]=cimag(FTdx[j*nvertf+k]);
      FTdyr[j+k*ntpoints+cind]=creal(FTdy[j*nvertf+k]);
      FTdyi[j+k*ntpoints+cind]=cimag(FTdy[j*nvertf+k]);
      FTdzr[j+k*ntpoints+cind]=creal(FTdz[j*nvertf+k]);
      FTdzi[j+k*ntpoints+cind]=cimag(FTdz[j*nvertf+k]);
    }
  }
  
  free(FT);
  free(FTdx);
  free(FTdy);
  free(FTdz);
  free(FTdA);
  free(FTdoff);
 

}

}
complex double Vector_Multiply(complex double *A, double *B,int n)
{
  complex double res=0;
  for(int i=0;i<n;i++)
    res+=A[i]*B[i];
  return res;
}
void Matrix_Multiply(complex double *A,double *B,complex double *C,int m,int n,int l)
{
  //NB: It is assumed that the matrix B is in MATLAB format, ie columns first
  //However, matrices A,C use C type indexing
  //This algorithm does not work for usual matrices
 //NB: We could openMP for fast implementation:
  // omp_set_num_threads(8);
//#pragma omp parallel for shared(A,B,C)
  for(int i=0;i<m;i++)
    for(int j=0;j<l;j++)
    {
      C[i*l+j]=Vector_Multiply(&A[n*i],&B[n*j],n);
    }
}