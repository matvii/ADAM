#include "prepare.h"
void mexFunction(int nlhs, mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
  /*
   * INPUT:
   * tlist,vlist,E,E0
   */
  /*
   * OUTPUT:
   * visibility
   */
  int *tlist;
  double *vlist,*vl,*tl,*Ei,*E0i,*E,*E0;
  int nfac,nvert,nE;
  int *vis_out,*visible;
 
  if(nrhs!=4)
    mexErrMsgTxt("4 input arguments required: tlist,vlist,E,E0");
  if(mxGetN(prhs[0])!=3 || mxGetN(prhs[1])!=3 || mxGetN(prhs[2])!=3 ||mxGetN(prhs[3])!=3)
    mexErrMsgTxt("Matrix dimensions must be nfac x 3 and nvert x 3,nE x 3,nE x 3");
  if(mxGetM(prhs[2])!=mxGetM(prhs[3]))
    mexErrMsgTxt("Number of view and light directions must be equal");
  nfac=mxGetM(prhs[0]);
  nvert=mxGetM(prhs[1]);
  nE=mxGetM(prhs[2]);
  
  
  tl=mxGetPr(prhs[0]);
  vl=mxGetPr(prhs[1]);
  Ei=mxGetPr(prhs[2]);
  E0i=mxGetPr(prhs[3]);
  vlist=mxCalloc(3*nvert,sizeof(double));
  tlist=mxCalloc(3*nfac,sizeof(int));
  visible=mxCalloc(nE*nfac,sizeof(int));
  E=mxCalloc(3*nE,sizeof(double));
  E0=mxCalloc(3*nE,sizeof(double));
  //mexPrintf("nfac: %d nvert: %d nE: %d\n",nfac,nvert,nE);
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
  
  for(int j=0;j<nE;j++)
  {
    for(int k=0;k<3;k++)
    {
  
      E[3*j+k]=Ei[k*nE+j];
      E0[3*j+k]=E0i[k*nE+j];
    }
  }
  
 //Allocate for output
 plhs[0]=mxCreateNumericMatrix(nE,nfac,mxINT32_CLASS,mxREAL);
 vis_out=mxGetData(plhs[0]);

FindActualBlockers(tlist,vlist,nfac,nvert,E,E0,nE,visible);
 //Copy to output
 
 for(int j=0;j<nE;j++)
 {
   for (int k=0;k<nfac;k++)
   {
     vis_out[k*nE+j]=visible[nfac*j+k];
   }
 }

}