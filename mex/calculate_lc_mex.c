#include "prepare.h"
void calculate_lcurve(int nE,int numfac,int numvert,int **tlist,double **vlist,double *angles,double **Eo,double **E0o,double *TIME,double *bright,double *dbrightx,double *dbrighty,double *dbrightz,double *dbrightb,double *dbrightl,double *drighto);
void mexFunction(int nlhs, mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
  double *e,*e0,*angles,*TIME;
 int numfac,numvert,nE;
 double  *vl, *tl,*bright,*dbrightx2,*dbrighty2,*dbrightz2,*dbrightb,*dbrightl,*dbrighto;
 double *dbrightx,*dbrighty,*dbrightz;
 numvert=mxGetM(prhs[0]);
 numfac=mxGetM(prhs[1]);
 double nvert=numvert;
 nE=mxGetM(prhs[2]);
 
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
 

 int **tlist;
 tlist=mxCalloc(numfac,sizeof(int*));
 for(int i=0;i<numfac;i++)
   tlist[i]=mxCalloc(3,sizeof(int));
 

 plhs[0]=mxCreateDoubleMatrix(1,nE,mxREAL);
  plhs[1]=mxCreateDoubleMatrix(nE,numvert,mxREAL);
  plhs[2]=mxCreateDoubleMatrix(nE,numvert,mxREAL);
  plhs[3]=mxCreateDoubleMatrix(nE,numvert,mxREAL);
  
  plhs[4]=mxCreateDoubleMatrix(1,nE,mxREAL);
  plhs[5]=mxCreateDoubleMatrix(1,nE,mxREAL);
  plhs[6]=mxCreateDoubleMatrix(1,nE,mxREAL);

bright=mxGetPr(plhs[0]);
dbrightx2=mxGetPr(plhs[1]);
dbrighty2=mxGetPr(plhs[2]);
dbrightz2=mxGetPr(plhs[3]);
dbrightb=mxGetPr(plhs[4]);
dbrightl=mxGetPr(plhs[5]);
dbrighto=mxGetPr(plhs[6]);
dbrightx=mxCalloc(nE*nvert,sizeof(double));
dbrighty=mxCalloc(nE*nvert,sizeof(double));
dbrightz=mxCalloc(nE*nvert,sizeof(double));
if(nrhs!=6)
  
    mexErrMsgTxt("6 input arguments required:a list of vertices, a list of triangles, view directions and light directions,TIME,angles\n");
if((!mxIsNumeric(prhs[0]))||(!mxIsNumeric(prhs[1]))||(!mxIsNumeric(prhs[2]))||(!mxIsNumeric(prhs[3])))
    mexErrMsgTxt("Input must be numeric.");
if((mxGetN(prhs[0])!=3)||(mxGetN(prhs[1])!=3)||(mxGetN(prhs[2])!=3)||(mxGetN(prhs[2])!=3))
    mexErrMsgTxt("Size of input matrices must be: nvert x 3,nfac x 3, nE x 3 and nE x 3");
if(nE!=mxGetM(prhs[3]))
    mexErrMsgTxt("Number of view and light directions must be equal");



  
vl=mxGetPr(prhs[0]);
tl=mxGetPr(prhs[1]);
e=mxGetPr(prhs[2]);
e0=mxGetPr(prhs[3]);
TIME=mxGetPr(prhs[4]);
angles=mxGetPr(prhs[5]);

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
 


 calculate_lcurve(nE,numfac,numvert,tlist,vlist,angles,E,E0,TIME,bright,dbrightx,dbrighty,dbrightz,dbrightb,dbrightl,dbrighto);

for(int i=0;i<numvert;i++)
  for(int j=0;j<nE;j++)
  {
    dbrightx2[j+i*nE]=dbrightx[j*numvert+i];
    dbrighty2[j+i*nE]=dbrighty[j*numvert+i];
    dbrightz2[j+i*nE]=dbrightz[j*numvert+i];
  }

}