#include "prepare.h"

void mexFunction(int nlhs, mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
  double *e,*e0;
 int numfac,numvert,nE;
 double *n,*c, *vl, *tl;
 numvert=mxGetM(prhs[0]);
 numfac=mxGetM(prhs[1]);
 nE=mxGetM(prhs[2]);
 //double vlist[numvert][3],normal[numfac][3],centroid[numfac][3],E[nE][3],E0[nE][3];
 double **vlist,**normal,**centroid,**E,**E0;
 vlist=mxCalloc(numvert,sizeof(double*));
 normal=mxCalloc(numfac,sizeof(double*));
 centroid=mxCalloc(numfac,sizeof(double*));
 E=mxCalloc(nE,sizeof(double*));
 E0=mxCalloc(nE,sizeof(double*));
 
 for(int i=0;i<numvert;i++)
   vlist[i]=mxCalloc(3,sizeof(double));
 for(int i=0;i<numfac;i++)
 {
   normal[i]=mxCalloc(3,sizeof(double));
   centroid[i]=mxCalloc(3,sizeof(double));
 }
 for(int i=0;i<nE;i++)
 {
   E[i]=mxCalloc(3,sizeof(double));
   E0[i]=mxCalloc(3,sizeof(double));
 }
 
//  int tlist[numfac][3];
 int **tlist;
 tlist=mxCalloc(numfac,sizeof(int*));
 for(int i=0;i<numfac;i++)
   tlist[i]=mxCalloc(3,sizeof(int));
 /*OUTPUT*/
  plhs[0]=mxCreateDoubleMatrix(numfac,3,mxREAL);
  plhs[1]=mxCreateDoubleMatrix(numfac,3,mxREAL);
plhs[2]=mxCreateNumericMatrix(1,numfac,mxINT32_CLASS,mxREAL);
plhs[3]=mxCreateNumericMatrix(numfac,numfac,mxINT32_CLASS,mxREAL);
plhs[4]=mxCreateNumericMatrix(nE,numfac,mxINT32_CLASS,mxREAL);

if(nrhs!=4)
    mexErrMsgTxt("4 input arguments required:a list of vertices, a list of triangles, view directions and light directions.\n");
if((!mxIsNumeric(prhs[0]))||(!mxIsNumeric(prhs[1]))||(!mxIsNumeric(prhs[2]))||(!mxIsNumeric(prhs[3])))
    mexErrMsgTxt("Input must be numeric.");
if((mxGetN(prhs[0])!=3)||(mxGetN(prhs[1])!=3)||(mxGetN(prhs[2])!=3)||(mxGetN(prhs[2])!=3))
    mexErrMsgTxt("Size of input matrices must be: nvert x 3,nfac x 3, nE x 3 and nE x 3");
if(nE!=mxGetM(prhs[3]))
    mexErrMsgTxt("Number of view and light directions must be equal");


n=mxGetPr(plhs[0]);
c=mxGetPr(plhs[1]);

/* Some variables for prepare function*/

int *nbl;
/*nbl=mxCalloc(numfac,sizeof(int));*/
nbl=mxGetData(plhs[2]);
int *ibll;
ibll=mxGetData(plhs[3]);
/*ibll=mxCalloc(numfac*numfac,sizeof(int));
/*ibl=mxCalloc(numfac,sizeof(int*));
for(int j=0;j<numfac;j++)
  ibl[j]=&ibll[j*numfac]; /*for 2d indexing*/


  
vl=mxGetPr(prhs[0]);
tl=mxGetPr(prhs[1]);
e=mxGetPr(prhs[2]);
e0=mxGetPr(prhs[3]);
for(int j=0;j<3;j++)
{
  for(int i=0;i<numvert;i++)
  {
    vlist[i][j]=vl[i+j*numvert];
  }
}

for(int j=0;j<3;j++)
{
  for(int i=0;i<numfac;i++)
    tlist[i][j]=*(tl+i+numfac*j);
}

for(int j=0;j<3;j++)
  for(int i=0;i<nE;i++)
  {
    E[i][j]=e[i+nE*j];
   E0[i][j]=e0[i+nE*j];
  }
 
/*mexPrintf("nvert: %d nfac: %d\n",nvert,nfac);
mexPrintf("%d %d %d\n",tlist[nfac-1][0],tlist[nfac-1][1],tlist[nfac-1][2]);
mexPrintf("%f %f %f",vlist[nvert-1][0],vlist[nvert-1][1],vlist[nvert-1][2]);
*/

prepare(numfac,numvert,vlist,tlist,normal,centroid,nbl,ibll);


for(int j=0;j<3;j++)
{
  for(int i=0;i<numfac;i++)
  {
    n[i+j*numfac]=normal[i][j];
    c[i+j*numfac]=centroid[i][j];
  }
}
int *visiblel,*visiblel2,**visible2;


visiblel2=mxGetData(plhs[4]);
visiblel=mxCalloc(nE*numfac,sizeof(int));
visible2=mxCalloc(nE,sizeof(int*));
findblockers(numfac,numvert,nE,E,E0,normal,centroid,nbl,ibll,vlist,tlist,visiblel);
for(int j=0;j<nE;j++)
visible2[j]=&visiblel[j*numfac];
for(int j=0;j<numfac;j++)
{
  for(int i=0;i<nE;i++)
  {
    visiblel2[i+j*nE]=visible2[i][j];
    
  }
}
}