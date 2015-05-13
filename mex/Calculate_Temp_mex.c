#include "prepare.h"
#include "num_of_threads.h"
#include<omp.h>
void mexFunction(int nlhs, mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
 /*Same as the original, only exception is the inclusion of matrix D (For effective memory usage)
  */
 //OUTPUT: FT FTdx FTdy FTdz FTdA FTdoff
  if(nrhs!=10)
    mexErrMsgTxt("10 input arguments required:a list of  triangles, vertices, angles,E0,TIME,Hdist,Wavelength,A,Gamma,N\n");
  //INPUT
  double A,R,Gamma;
  int N;
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
 E0=mxGetPr(prhs[3]);
  double TIME;
  TIME=(double)mxGetScalar(prhs[4]);
  
  double Hdist;
  Hdist=(double)mxGetScalar(5);
  double WL;
  WL=(double)mxGetScalar(prhs[6]);
  A=(double)mxGetScalar(prhs[7]);
  Gamma=(double)mxGetScalar(prhs[8]);
  N=(int)mxGetScalar(prhs[9]);
  double *Temp;
  plhs[0]=mxCreateDoubleMatrix(nfac,1,mxREAL);
  Temp=mxGetPr(plhs[0]);
  //Temp=mxCalloc(nfac,sizeof(double));
  Calculate_Temp(tlist,vlist,nfac,nvert,angles,E0,TIME,Gamma,A,Hdist,N,Temp,0);
  
}