void mexFunction(int nlhs, mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
 /*Same as the original, only exception is the inclusion of matrix D (For effective memory usage)
  */
 //OUTPUT: FT FTdx FTdy FTdz FTdA FTdoff
  if(nrhs!=17)
    mexErrMsgTxt("17 input arguments required:a list of  triangles, vertices, angles,FT.E,FT.E0,FT.TIME,FT.freq,FT.up,FT.dist,FT.Hdist,Wavelength,A,R,Gamma,N,offset,deriv\n");
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
 // E=mxGetPr(prhs[3]);
  double TIME;
  //TIME=(double)mxGetScalar(prhs[4]);
  
  int nfreq;
  //nfreq=mxGetM(prhs[5]);
  //if(mxGetM(prhs[5])!=mxGetM(prhs[6]))
   // mexErrMsgTxt("freqx and freqy must be of equal length");
  
//  freqx=mxGetPr(prhs[5]);
 // freqy=mxGetPr(prhs[6]);
  A=(double)mxGetScalar(prhs[11]);
  R=(double)mxGetScalar(prhs[12]);
  Gamma=(double)mxGetScalar(prhs[13]);
  N=(int)mxGetScalar(prhs[14]);
  double *offset;
  offset=mxGetPr(prhs[15]);
  
  //f0=(double)mxGetScalar(prhs[9]);
  int deriv;
  deriv=(int)mxGetScalar(prhs[16]);
  
  //double complex *FT, *FTdx,*FTdy,*FTdz,*FTdA,*FTdoff;
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
  FTWLcell=prhs[10];
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