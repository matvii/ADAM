#include "prepare.h"
void mexFunction(int nlhs, mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
  /*
   * INPUT:
   * tlist,vlist
   */
  /*
   * OUTPUT:
   * normal,centroid,numberofblockers,IndexofBlockers
   */
  int *tlist;
  double *vlist,*vl,*tl;
  int nfac,nvert;
  double *normal,*centroid;
  int *NumofBlocks,*IndexofBlocks;
  double *normal_out,*centroid_out;
  int *IoB_out;
  
  if(nrhs!=2)
    mexErrMsgTxt("2 input arguments required: tlist and vlist");
  if(mxGetN(prhs[0])!=3 || mxGetN(prhs[1])!=3)
    mexErrMsgTxt("Matrix dimensions must be nfac x 3 and nvert x 3");
  nfac=mxGetM(prhs[0]);
  nvert=mxGetM(prhs[1]);
   
  tl=mxGetPr(prhs[0]);
  vl=mxGetPr(prhs[1]);
  vlist=mxCalloc(3*nvert,sizeof(double));
  tlist=mxCalloc(3*nfac,sizeof(int));
  
  //Convert from row first to column first
  for(int j=0;j<nvert;j++)
  {
    for (int k=0;k<3;k++)
    {
      vlist[3*j+k]=vl[k*nvert+j];
    }
  }
  for(int j=0;j<nfac;j++)
  {
    for (int k=0;k<3;k++)
    {
      tlist[3*j+k]=(int)tl[k*nfac+j];
    }
  }
  //Allocate for memory
  normal=mxCalloc(3*nfac,sizeof(double));
  centroid=mxCalloc(3*nfac,sizeof(double));
  IndexofBlocks=mxCalloc(nfac*nfac,sizeof(int));
  
  //Prepare output variables
  plhs[0]=mxCreateDoubleMatrix(nfac,3,mxREAL);
  plhs[1]=mxCreateDoubleMatrix(nfac,3,mxREAL);
  plhs[2]=mxCreateNumericMatrix(1,nfac,mxINT32_CLASS,mxREAL);
  plhs[3]=mxCreateNumericMatrix(nfac,nfac,mxINT32_CLASS,mxREAL);
  normal_out=mxGetPr(plhs[0]);
  centroid_out=mxGetPr(plhs[1]);
  NumofBlocks=mxGetData(plhs[2]);
  IoB_out=mxGetData(plhs[3]);
  
  FacetsOverHorizon(tlist,vlist,nfac,nvert,normal,centroid,NumofBlocks,IndexofBlocks);
  
  
  //Copy to output variables
  for(int j=0;j<nfac;j++)
  {
    for (int k=0;k<3;k++)
    {
      normal_out[k*nfac+j]=normal[3*j+k];
      centroid_out[k*nfac+j]=centroid[3*j+k];
    }
  }
  for(int j=0;j<nfac;j++)
  {
    for (int k=0;k<nfac;k++)
    {
      IoB_out[k*nfac+j]=IndexofBlocks[j*nfac+k];
    }
  }
  
  
}
  
    
  