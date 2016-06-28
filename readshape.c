#include<stdio.h>
#include<stdlib.h>
void read_shape(char* filename,int **facets2,double** vertices2,int *nfac2,int *nvert2,int type3)
{
    /* Read shape and return vertices and facets. If type3 is nonzero, then file is assumed to contain '3' between
     * every facet */
  int nvert,nfac;
  FILE *fid;
  fid=fopen(filename,"r");
  fscanf(fid,"%d %d",&nvert,&nfac);
  //printf("fid%d\n",(int)fid);
  //printf("facets: %d vertices %d\n",nfac,nvert);
  int *facets;
  double *vertices;
  facets=malloc(3*nfac*sizeof(int));
  vertices=malloc(3*nvert*sizeof(double));
  for(int i=0;i<nvert;i++)
    fscanf(fid,"%lf %lf %lf",&vertices[i*3],&vertices[i*3+1],&vertices[i*3+2]);
  if(type3)
  for(int i=0;i<nfac;i++)
  {
      fscanf(fid,"%d");
      fscanf(fid,"%d %d %d",&facets[i*3],&facets[i*3+1],&facets[i*3+2]);
  }
  else
      for(int i=0;i<nfac;i++)
          fscanf(fid,"%d %d %d",&facets[i*3],&facets[i*3+1],&facets[i*3+2]);
  fclose(fid);
  *nvert2=nvert;
  *nfac2=nfac;
  *facets2=facets;
  *vertices2=vertices;
}

/*
int main()
{
    char filen[]="shape.txt";
    double *vertices;
    int *facets;
    int nfac,nvert;
    read_shape(filen,&facets,&vertices,&nvert,&nfac,0);
    printf("facets: %d vertices %d\n",nfac,nvert);
    printf("First facet: %d %d %d\n",facets[0],facets[1],facets[2]);
    printf("Last facet: %d %d %d\n",facets[3*nfac-3],facets[3*nfac-2],facets[3*nfac-1]);
    printf("First vertex: %f %f %f\n",vertices[0],vertices[1],vertices[2]);
    printf("Last vertex: %f %f %f\n",vertices[3*nvert-3],vertices[3*nvert-2],vertices[3*nvert-1]);
    
}
*/
