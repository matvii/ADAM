#include"prepare.h"
void FindActualBlockers(int *tlist,double *vlist,int nfac,int nvert,double *E,double *E0,int nE,int *visible)
{
  /*Find actual blocked facets, taking into account the view and illumination directions
   * INPUT
   * E,E0 nEx3 vectors, view and illumination directions
   * OUTPUT
   * visible nE x nfac matrix, ==1 if facet is visible 
   */
  double *normal,*centroid;
  double *cnormal,*ccent;
  double *bnormal;
  double mu,mu0,bmu,bmu0;
  double *cE,*cE0;
  int i1,i2,i3;
  double *v1,*v2,*v3;
  int blockindex;
  int *IndexofBlocks;
  int *NumofBlocks;
  int blocked=0;
  normal=calloc(nfac*3,sizeof(double));
  centroid=calloc(nfac*3,sizeof(double));
  NumofBlocks=calloc(nfac,sizeof(int));
  IndexofBlocks=calloc(nfac*nfac,sizeof(int));
  //Find potential blockers, independent of view and ill directions
  FacetsOverHorizon(tlist,vlist,nfac,nvert,normal,centroid,NumofBlocks,IndexofBlocks);
  //Loop over observations
  for(int j=0;j<nE;j++)
  {
    cE=E+j*3;
    cE0=E0+j*3;
    //Loop over facets
    for(int k=0;k<nfac;k++)
    {
      blocked=0;
      //Take facet normal
      cnormal=normal+k*3;
      ccent=centroid+k*3;
      //Is the current facet illuminated and visible?
      mu=DOT(cE,cnormal);
      mu0=DOT(cE0,cnormal);
      
	
      if(mu<0 || mu0<0)
	continue;
      //facet is seen, test for blockers
      if(NumofBlocks[k]==0)
      {
	visible[nfac*j+k]=1;
	continue;
      }
      //Ok, so there are facets above local horizon,
      //test for blockers
      for(int l=0;l<NumofBlocks[k];l++)
      {
	blockindex=IndexofBlocks[k*nfac+l]-1;
	
	  //CONSIDER ONLY BACKFACING FACETS, IE NOT ILLUMINATED OR NOT VISIBLE
	  //IF THE MESH IS OK, THERE WILL BE NO PROBLEMS
	  bnormal=normal+3*blockindex;
	  bmu=DOT(cE,bnormal);
	  bmu0=DOT(cE0,bnormal);
	  //Is the blocking facet visible? If it is skip it, since there will be an invisible blocking facet too
	 if(bmu>0 && bmu0>0)
	    continue;
	  
	  i1=tlist[3*blockindex]-1;
	  i2=tlist[3*blockindex+1]-1;
	  i3=tlist[3*blockindex+2]-1;
	  v1=vlist+3*i1;
	  v2=vlist+3*i2;
	  v3=vlist+3*i3;
	  //Test for eye obstruction
	  //if(k==26 && l==67)
	 //mexPrintf("facet:%d blocker: %d, is_in_triangle: %d\n block: %d\n",k+1,l+1,is_in_triangle(ccent,cE,v1,v2,v3),k*nfac+l); 
	  
	  if(bmu<0 && is_in_triangle(ccent,cE,v1,v2,v3))
	  {
	   // if(k==26)
	    //mexPrintf("in if: k:%d blockid: %d\n",k+1,blockindex+1);
	    blocked=1;
	    break;
	  }
	  //Test for light obstruction
	  if(bmu0<0 && is_in_triangle(ccent,cE0,v1,v2,v3))
	  {
	    blocked=1;
	    break;
	  }
	
      }
      
      //No blockers, so visible
      if(blocked==0)
      {
	//mexPrintf("Not blocked: %d\n",k+1);
	visible[nfac*j+k]=1;
      }
    }
  }
  free(normal);
  free(centroid);
  free(NumofBlocks);
  free(IndexofBlocks);
}
      
      
    
  
  