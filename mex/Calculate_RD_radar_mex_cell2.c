#include "prepare.h"
#include<omp.h>
#include"num_of_threads.h"
void Calculate_Range_Doppler(int *tlist,double *vlist,int nfac,int nvert,double *angles,double *Eo,double TIME,double *freqx,double *freqy,int nfreq,double rfreq,double *offset,double scal,double rexpe,double complex* F);
void Calculate_Range_Doppler_deriv(int *tlist,double *vlist,int nfac,int nvert,double *angles,double *Eo,double TIME,double *freqx,double *freqy,int nfreq,double rfreq,double *offset,double scal,double rexpe,double complex* F,double complex *FTdx,double complex *FTdy,double complex *FTdz,double complex *FTdA,double complex *FTdoff,double complex *FTdexp);
void Matrix_Multiply(complex double *A, double *B,complex double *C,int n,int m,int l);
void mexFunction(int nlhs, mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
 /*Same as the original, only exception is the inclusion of matrix D (For effective memory usage)
  */
  if(nrhs!=12)
    mexErrMsgTxt("12 input arguments required:a list of  triangles, vertices, angles,FT.E,FT.TIME,FT.freq,FT.radarfreq,offset,scale,rexp,matrix D,deriv\n");
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
  double *E;
 // E=mxGetPr(prhs[3]);
  double TIME;
  //TIME=(double)mxGetScalar(prhs[4]);
  
  int nfreq;
  //nfreq=mxGetM(prhs[5]);
  //if(mxGetM(prhs[5])!=mxGetM(prhs[6]))
   // mexErrMsgTxt("freqx and freqy must be of equal length");
  
//  freqx=mxGetPr(prhs[5]);
 // freqy=mxGetPr(prhs[6]);
  double *offset;
  offset=mxGetPr(prhs[7]);
  double *scale;
  double f0,rexp;
  scale=mxGetPr(prhs[8]);
  //f0=(double)mxGetScalar(prhs[9]);
  rexp=(double)mxGetScalar(prhs[9]);
  int deriv;
  deriv=(int)mxGetScalar(prhs[11]);
  //OUTPUT: FT FTdx FTdy FTdz FTdA FTdoff FTdexp
  //double complex *FT, *FTdx,*FTdy,*FTdz,*FTdA,*FTdoff,*FTdexp;
  //Cell handling
  mwSize nobs;
  mwIndex index;
  mxArray *Freq_ptr;
  mxArray *FTEcell,*FTTIMEcell,*FTfreqcell,*FTrfreqcell;
  //double *FTE,*FTTIME,*FTfreq,*FTrfreq;
  int M,N;
  int *nopoints,*cumpoints,ntpoints;
  ntpoints=0;
  //double *freqx,*freqy,*FTfreq;
  FTEcell=prhs[3];
  nobs = mxGetNumberOfElements(FTEcell); 
  FTTIMEcell=prhs[4];
  FTfreqcell=prhs[5];
  FTrfreqcell=prhs[6];
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
    double *FTE,*FTTIME,*FTfreq,*FTrfreq;
    double complex *FT;
   //  obsind=omp_get_thread_num();
    FT=(double complex*)calloc(nopoints[obsind],sizeof(double complex));
   
   //plhs[0]=mxCreateDoubleMatrix(nfreq,1,mxCOMPLEX);
   //double *FTr,*FTi;
   //FTr=mxGetPr(plhs[0]);
    //FTi=mxGetPi(plhs[0]);
    FTE=mxGetPr(mxGetCell(FTEcell,obsind));
    FTTIME=mxGetPr(mxGetCell(FTTIMEcell,obsind));
    FTfreq=mxGetPr(mxGetCell(FTfreqcell,obsind));
    FTrfreq=mxGetPr(mxGetCell(FTrfreqcell,obsind));
  //  printf("obsind: %d nopoints: %d cumpoints: %d total: %d radarfreq: %f\n",obsind,nopoints[obsind],cumpoints[obsind],ntpoints,*FTrfreq);
  //printf("TIME: %f freqs: %f %f nfreq: %d offset: %f %f\n",*FTTIME,*FTfreq,*(FTfreq+nopoints[obsind]),nopoints[obsind],*offset,*(offset+1));
    Calculate_Range_Doppler(tlist,vlist,nfac,nvert,angles,FTE,*FTTIME,FTfreq,FTfreq+nopoints[obsind],nopoints[obsind],*FTrfreq,offset+2*obsind,*(scale+obsind),rexp,FT);
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

double *D;
int Ds;
int nvertf;
Ds=mxGetM(prhs[10]);
if(Ds>0)
{
  D=mxGetPr(prhs[10]);
  nvertf=mxGetN(prhs[10]);
}
else
  nvertf=nvert;


  
  
 // plhs[0]=mxCreateDoubleMatrix(ntpoints,1,mxCOMPLEX);
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
    double complex *FTdx,*FTdy,*FTdz,*FTdA,*FTdoff,*FTdexp,*FTdxf,*FTdyf,*FTdzf;
     double *FTE,*FTTIME,*FTfreq,*FTrfreq;
    double complex *FT;
   //  obsind=omp_get_thread_num();
    FT=(double complex*)calloc(nopoints[obsind],sizeof(double complex));
   FTdx=(double complex*)calloc(nopoints[obsind]*nvertf,sizeof(double complex));
  FTdy=(double complex*)calloc(nopoints[obsind]*nvertf,sizeof(double complex));
  FTdz=(double complex*)calloc(nopoints[obsind]*nvertf,sizeof(double complex));
  FTdA=(double complex*)calloc(nopoints[obsind]*3,sizeof(double complex));
  FTdoff=(double complex*)calloc(nopoints[obsind]*2,sizeof(double complex));
  FTdexp=(double complex*)calloc(nopoints[obsind],sizeof(double complex));
  
  
    FTE=mxGetPr(mxGetCell(FTEcell,obsind));
    FTTIME=mxGetPr(mxGetCell(FTTIMEcell,obsind));
    FTfreq=mxGetPr(mxGetCell(FTfreqcell,obsind));
    FTrfreq=mxGetPr(mxGetCell(FTrfreqcell,obsind));
    if(Ds>0)
    {
      FTdxf=(double complex*)calloc(nopoints[obsind]*nvert,sizeof(double complex));
      FTdyf=(double complex*)calloc(nopoints[obsind]*nvert,sizeof(double complex));
      FTdzf=(double complex*)calloc(nopoints[obsind]*nvert,sizeof(double complex));
      Calculate_Range_Doppler_deriv(tlist,vlist,nfac,nvert,angles,FTE,*FTTIME,FTfreq,FTfreq+nopoints[obsind],nopoints[obsind],*FTrfreq,offset+2*obsind,*(scale+obsind),rexp,FT,FTdxf,FTdyf,FTdzf,FTdA,FTdoff,FTdexp);
      //Convert from vlistn->vlist. Only because we want to minimize memory usage
      Matrix_Multiply(FTdxf,D,FTdx,nopoints[obsind],nvert,nvertf);
      free(FTdxf);
      Matrix_Multiply(FTdyf,D,FTdy,nopoints[obsind],nvert,nvertf);
      free(FTdyf);
      Matrix_Multiply(FTdzf,D,FTdz,nopoints[obsind],nvert,nvertf);
      free(FTdzf);
    }
    else
      Calculate_Range_Doppler_deriv(tlist,vlist,nfac,nvert,angles,FTE,*FTTIME,FTfreq,FTfreq+nopoints[obsind],nopoints[obsind],*FTrfreq,offset+2*obsind,*(scale+obsind),rexp,FT,FTdx,FTdy,FTdz,FTdA,FTdoff,FTdexp);
      
 
  
  //Copy variables to matlab
  cind=cumpoints[obsind];
  oind=nopoints[obsind];
  for(int j=0;j<nopoints[obsind];j++)
  {
    FTr[j+cind]=creal(FT[j]);
    FTi[j+cind]=cimag(FT[j]);
    FTdexpr[j+cind]=creal(FTdexp[j]);
    FTdexpi[j+cind]=cimag(FTdexp[j]);
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
  free(FTdexp);

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
  for(int i=0;i<m;i++)
    for(int j=0;j<l;j++)
    {
      C[i*l+j]=Vector_Multiply(&A[n*i],&B[n*j],n);
    }
}