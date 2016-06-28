#include"utils.h"
void FacetsOverHorizon(int *tlist,double *vlist,int nfac,int nvert,double *normal,double *centroid,int *NumofBlocks,int *IndexofBlocks)
{
  /*INPUT:
   * tlist 
   * vlist
   * NOTE: MATLAB uses columns first ordering, C uses rows first. We assume rows first because it is convenient
   * nfac number of facets
   * nvert number of vertices
   */
  /*OUTPUT:
   * Normal List of facet normals, rows first ordering since these will not leave C program
   * centroid, List of facet centroids, rows first
   * NumofBlocks number of facets above current facet, 1xnfac vector
   * IndexofBlocks 1xnfac^2  vector, indexes the facets blocking the current facet NOT BINARY MATRIX
   * NOTE INDEXING FROM 1!!
   */
 /*NOTE: FACET A IS ABOVE FACET B IF 1) ONE OF ITS VERTICES IS ABOVE THE LOCAL HORIZON OF B AND
  *2) FACET A IS FACING TOWARDS FACET B
  *CONDITION 2 IS SOMEWHAT CONTROVERSIAL. CHECK THIS!
  */ 
 
 
  int j1,j2,j3;
  int vindex;
  double *v1,*v2,*v3;
  double*w;
  double side1[3],side2[3];
  double cvvec[3];
  double norm;
  double cnormal[3];
  double *jcent,*jnormal,*knormal;
  double cent[3];
 /*
  for(int j=0;j<nfac*nfac;j++)
    IndexofBlocks[j]=0;
  for(int j=0;j<nfac;j++)
    NumofBlocks[j]=0;
  */
  //For each facet
  for(int j=0;j<nfac;j++)
  {
    //Vertex indices of the current facet
    //Note that C indices from 0, matlab from 1
    j1=tlist[j*3]-1;
    j2=tlist[j*3+1]-1;
    j3=tlist[j*3+2]-1;
    //Current vertices
    v1=vlist+j1*3;
    v2=vlist+j2*3;
    v3=vlist+j3*3;
    //Calculate normals and centroids
    for(int i=0;i<3;i++)
    {
      side1[i]=v2[i]-v1[i];
      side2[i]=v3[i]-v1[i];
      cent[i]=(v1[i]+v2[i]+v3[i])/3.0;
    }
    cross(side1,side2,cnormal);
    norm=NORM(cnormal);
    
    //Store centroids and normals
    for(int i=0;i<3;i++)
    {
      cnormal[i]=cnormal[i]/norm;
      normal[3*j+i]=cnormal[i];
      centroid[3*j+i]=cent[i];
    }
  }
  for(int j=0;j<nfac;j++)
  {
    //Vertex indices of the current facet
    //Note that C indices from 0, matlab from 1
    j1=tlist[j*3]-1;
    j2=tlist[j*3+1]-1;
    j3=tlist[j*3+2]-1;
    //Current vertices
    v1=vlist+j1*3;
    v2=vlist+j2*3;
    v3=vlist+j3*3;
    jcent=centroid+3*j;
    jnormal=normal+3*j;
    //Find potential blockers
    for(int k=0;k<nfac;k++)
    {
      if(k==j)
	continue;
      //Vertex indices of the possible blocker facet
     
      //Facet is above local horizon of the current facet if dot product of facet normal with the vector from facet centroid to
      //a vertex is positive.
      //Check each vertex
      //Facet normal
      knormal=normal+3*k;
      
      for(int i=0;i<3;i++)
      {
	vindex=tlist[k*3+i]-1;
	
	if((vindex==j1) || (vindex==j2) ||( vindex==j3))
	  continue;
	w=vlist+vindex*3; //Points to the current vertex
	cvvec[0]=w[0]-jcent[0];
	cvvec[1]=w[1]-jcent[1];
	cvvec[2]=w[2]-jcent[2];
	
	if((DOT(cvvec,jnormal))>0 && (DOT(cvvec,knormal)<0)) //Vertex of the facet must be above local horizon, and the blocking facet must be facing the current facet
	{
	  IndexofBlocks[nfac*j+(NumofBlocks[j]++)]=k+1;
	  //NumofBlocks[j]++;
	  break;
	}
      }
    }
  }
}
	  
	  
	  
	  
	
	  
	
	
      
      
      
      
