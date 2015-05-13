#include "prepare.h"


void prepare(int numfac,int numvert,double **vlist,int **tlist,double **normal,double **centroid,int *nbl,int *ibll)
{
  double side1[3],side2[3],normN;


int *nfac;
nfac=mxCalloc(numvert,sizeof(int));
int *idone;
idone=mxCalloc(numfac,sizeof(int));

int **ibl;

ibl=mxCalloc(numfac,sizeof(int*));
for(int j=0;j<numfac;j++)
  ibl[j]=&ibll[j*numfac]; /*for 2d indexing*/
double vn[3],vec[3],vec2[3],vn2[3];
double altd;
int *ifacl;
int **ifac;
ifacl=mxCalloc(numvert*numfac,sizeof(int));
ifac=mxCalloc(numvert,sizeof(int*));
for(int j=0;j<numvert;j++)
  ifac[j]=&ifacl[j*numfac]; /*for 2d indexing*/


/*printf("sizeof: %d sizeof0: %d",sizeof(ifacl),sizeof(ifacl[0]));*/
for(int j=0;j<numfac;j++)
{
  idone[j]=0;
  nbl[j]=0;
 
}
for(int j=0;j<numvert;j++)
  nfac[j]=0;

  for(int j=0;j<numfac;j++)
  {
    
    for (int k=0;k<3;k++)
    {
      side1[k]=vlist[tlist[j][1]-1][k]-vlist[tlist[j][0]-1][k];
       side2[k]=vlist[tlist[j][2]-1][k]-vlist[tlist[j][0]-1][k];
    }
    cross(side1,side2,&normal[j][0]); /*cross product*/
    normN=NORM(normal[j]);
    normal[j][0]=normal[j][0]/normN;
    normal[j][1]=normal[j][1]/normN;
    normal[j][2]=normal[j][2]/normN;
    
    
    centroid[j][0]=(vlist[tlist[j][0]-1][0]+vlist[tlist[j][1]-1][0]+vlist[tlist[j][2]-1][0])/3;
    centroid[j][1]=(vlist[tlist[j][0]-1][1]+vlist[tlist[j][1]-1][1]+vlist[tlist[j][2]-1][1])/3;
    centroid[j][2]=(vlist[tlist[j][0]-1][2]+vlist[tlist[j][1]-1][2]+vlist[tlist[j][2]-1][2])/3;
    
  }
  /*mexPrintf("\n numvert:%d numfac: %d\n",numvert,numfac);*/
  for(int j=0;j<numfac;j++)
  {
    nfac[tlist[j][0]-1]=nfac[tlist[j][0]-1]+1;
   
 ifac[tlist[j][0]-1][nfac[tlist[j][0]-1]-1]=j+1;
    nfac[tlist[j][1]-1]=nfac[tlist[j][1]-1]+1;
    ifac[tlist[j][1]-1][nfac[tlist[j][1]-1]-1]=j+1;
    nfac[tlist[j][2]-1]=nfac[tlist[j][2]-1]+1;
    ifac[tlist[j][2]-1][nfac[tlist[j][2]-1]-1]=j+1;
  }
  
  for(int i=0;i<numfac;i++)
  {
    for(int k=0;k<3;k++)
	vn[k]=normal[i][k];
   
    for(int j=0;j<numvert;j++)
    {
      if(j==(tlist[i][0]-1)||j==(tlist[i][1]-1)||j==(tlist[i][2]-1))
	continue;
      for(int k=0;k<3;k++)
	vec[k]=vlist[j][k]-centroid[i][k];
      if(DOT(vec,vn)>EP)
      {
	for(int l=0;l<nfac[j];l++)
	  if(idone[ifac[j][l]-1]!=i+1)
	  {
	    idone[ifac[j][l]-1]=i+1;
	    for(int k=0;k<3;k++)
	      vn2[k]=normal[ifac[j][l]-1][k];
	    altd=DOT(vec,vn2);
	    
	    if(altd<0) 
	    {
	      nbl[i]=nbl[i]+1;
	      ibl[i][nbl[i]-1]=ifac[j][l];
	    }
	  }
      }
    }
  }
	    


      
    
    
}

