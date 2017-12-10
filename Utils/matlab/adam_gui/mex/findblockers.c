#include "prepare.h"


void findblockers(int nfac,int nvert,int nE,double **E,double **E0,double **normal,double **centroid,int *nbl,int *ibll, double **vlist,int **tlist,int *visiblel)
{
  double temp_a,temp_b;
  /*Input
   * E,E0 view and sun directions, nE x 3
   * normal, centroid double nfac x 3
   * nbl, nfac int array, number of blockers per triangle
   * ibl nfac x nfac array, list of blocker triangles
   * vlist, tlist nvertx3 and nfacx3 arrays
   * Output:
   
   
   * visible nE x nfac array, visible(i,j)=1 is triangle j is visible and not blocked when viewed from E(i), E0(i)*/
  
  /* Do some memory magic to allow familiar matrix indexing. Check if there is cleaner way to do this*/
  int **seen, **lit, **sandl, **visible,**ibl,*sandll,*seenl,*litl,*nseen,*nlit,*nsandl;
/*ifacl=mxCalloc(numvert*numfac,sizeof(int))*/
nseen=mxCalloc(nE,sizeof(int));
nsandl=mxCalloc(nE,sizeof(int));
nlit=mxCalloc(nE,sizeof(int));
sandll=mxCalloc(nE*nfac,sizeof(int));
seenl=mxCalloc(nE*nfac,sizeof(int));
litl=mxCalloc(nE*nfac,sizeof(int));
seen=mxCalloc(nE,sizeof(int*));
lit=mxCalloc(nE,sizeof(int*));
sandl=mxCalloc(nE,sizeof(int*));
visible=mxCalloc(nE,sizeof(int*));
ibl=mxCalloc(nfac,sizeof(int*));
for(int j=0;j<nfac;j++)
  ibl[j]=&ibll[j*nfac]; /*for 2d indexing*/

for(int j=0;j<nE;j++)
{
  seen[j]=&seenl[j*nfac]; /*for 2d indexing*/
  lit[j]=&litl[j*nfac];
  sandl[j]=&sandll[j*nfac];
  visible[j]=&visiblel[j*nfac];
}

  /*First determine which triangles are seen and lit:*/
  for(int j=0;j<nE;j++)
  {
    for(int k=0;k<nfac;k++)
    {
      temp_a=DOT(E[j],normal[k]);
      temp_b=DOT(E0[j],normal[k]);
      if((temp_a>EP)&&(temp_b>EP))
      {

	seen[j][k]=1;
	nseen[j]=nseen[j]+1;

	lit[j][k]=1;
	nlit[j]=nlit[j]+1;
	sandl[j][nsandl[j]]=k+1;
	nsandl[j]=nsandl[j]+1;
	visible[j][k]=1;
      }
	else if(temp_a>EP)
	{
	  seen[j][k]=1;
	  nseen[j]=nseen[j]+1;
	}
	else if(temp_b>EP)
	{
	  lit[j][k]=1;
	  nlit[j]=nlit[j]+1;
	}
	  
	  
	
   
    }
  }
int isequal=0;
    for(int j=0;j<nE;j++) /*for each view direction*/
    {
      isequal=0;
    if(pow(E[j][0]-E0[j][0],2)+pow(E[j][1]-E0[j][1],2)+pow(E[j][2]-E0[j][2],2)<1e-20)
      isequal=1;
      for(int jk=0;jk<nsandl[j];jk++) /*for each triangle which is seen and lit */
      {
	
	
	for(int jjk=0;jjk<nbl[sandl[j][jk]-1];jjk++)/*number of possible blockers for triangle sandl[j][jk]. Note that we started indexing from 1 */
	{
	 /*NB We could replace this seen/lit mess with one-sided triangle test*/
	  if((seen[j][ibl[sandl[j][jk]-1][jjk]-1]==0)&&(is_in_triangle(&centroid[sandl[j][jk]-1][0],&E[j][0],&vlist[tlist[ibl[sandl[j][jk]-1][jjk]-1][0]-1][0],&vlist[tlist[ibl[sandl[j][jk]-1][jjk]-1][1]-1][0],&vlist[tlist[ibl[sandl[j][jk]-1][jjk]-1][2]-1][0])==1))
	  {
	    visible[j][sandl[j][jk]-1]=0;
	
	    break;
	  }
	}
	

	if(visible[j][sandl[j][jk]-1]==1 && isequal==0)
	{
	 
	  for(int jjk=0;jjk<nbl[sandl[j][jk]-1];jjk++)
	  { 
	    
	   
	    if((lit[j][ibl[sandl[j][jk]-1][jjk]-1]==0)&&(is_in_triangle(&centroid[sandl[j][jk]-1][0],&E0[j][0],&vlist[tlist[ibl[sandl[j][jk]-1][jjk]-1][0]-1][0],&vlist[tlist[ibl[sandl[j][jk]-1][jjk]-1][1]-1][0],&vlist[tlist[ibl[sandl[j][jk]-1][jjk]-1][2]-1][0])==1))
	  {
	    
	    visible[j][sandl[j][jk]-1]=0;
	    
	    break;
	  }
	  }
	}
	
      }
    }
}


