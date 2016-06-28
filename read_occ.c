#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include"structs.h"
#include"utils.h"

OCstruct  *read_occ(char* filename,double min_tim,double *up)
{
  OCstruct *OC=malloc(sizeof(OCstruct));
  FILE *fid;
  fid=fopen(filename,"r");
  double temp;
  double E[3];
  int ntotal=0;
  int noc;
  double sumx=0,sumy=0;
  int nobs;
  double *data,*TIME,*etime,dist;
  int *type;
  if(fid==NULL)
  {
      perror("Cannot open OCC file!");
      exit(-1);
  }
  if(fscanf(fid,"%d",&noc)!=1)
  {
      perror("Error reading OCC file!");
      exit(-1);
  }
  
  OC->noc=noc;
  OC->nobs=calloc(noc,sizeof(int));
  OC->data=calloc(noc,sizeof(double*));
  OC->TIME=calloc(noc,sizeof(double*));
  OC->E=calloc(3*noc,sizeof(double));
  OC->V=calloc(3*noc,sizeof(double));
  OC->up=calloc(3*noc,sizeof(double));
  OC->distance=calloc(noc,sizeof(double));
  OC->type=calloc(noc,sizeof(int*));
  OC->etime=calloc(noc,sizeof(double*));
  
  for(int j=0;j<noc;j++)
  {
      //read direction vector and velocity
      fscanf(fid,"%lf %lf %lf",E,E+1,E+2);
      fscanf(fid,"%lf %lf %lf",OC->V+3*j,OC->V+3*j+1,OC->V+3*j+2);
      fscanf(fid,"%lf",&temp);
      fscanf(fid,"%d",&nobs);
      OC->nobs[j]=nobs;
      OC->data[j]=calloc(4*nobs,sizeof(double));
      OC->TIME[j]=calloc(2*nobs,sizeof(double));
      OC->type[j]=calloc(nobs,sizeof(int));
      OC->etime[j]=calloc(2*nobs,sizeof(double));
        data=OC->data[j];
        TIME=OC->TIME[j];
        etime=OC->etime[j];
        type=OC->type[j];
        dist=sqrt(pow(E[0],2)+pow(E[1],2)+pow(E[2],2));
        OC->E[3*j]=E[0]/dist;
        OC->E[3*j+1]=E[1]/dist;
        OC->E[3*j+2]=E[2]/dist;
        OC->distance[j]=dist;
        OC->up[3*j]=up[0];
        OC->up[3*j+1]=up[1];   
        OC->up[3*j+2]=up[2]; //Orientation of  camera in ecl
        sumx=0.0;
        sumy=0.0;
        ntotal+=nobs;
      for(int k=0;k<nobs;k++)
      {
          if(fscanf(fid,"%lf %lf %lf %lf %lf %lf %lf %lf %d",TIME+2*k,etime+2*k,data+4*k,data+4*k+1,TIME+2*k+1,etime+2*k+1,data+4*k+2,data+4*k+3,type+k)!=9)
          {
              perror("Error reading a line in OCC file!");
              exit(-1);
          }
          
          //substract min_tim and LT correction from OCC TIME
          TIME[2*k]=TIME[2*k]-min_tim-(dist/299792.458)/60/60/24;
          TIME[2*k+1]=TIME[2*k+1]-min_tim-(dist/299792.458)/60/60/24;
          //Calculate sum of coordinates
          sumx+=data[4*k]+data[4*k+2];
          sumy+=data[4*k+1]+data[4*k+3];
          
          
      }
     
      //Center the coordinates
      sumx=sumx/(2*nobs);
      sumy=sumy/(2*nobs);
      for(int k=0;k<nobs;k++)
      {
          data[4*k]-=sumx;
          data[4*k+2]-=sumx;
          data[4*k+1]-=sumy;
          data[4*k+3]-=sumy;
      }
          
  }
              
  OC->ntotal=ntotal;
  fclose(fid);
  return OC;
      
  
}
/*
void main()
{
    char filename[]="/home/matvii/matlab/Hertha.occ";
    OCstruct *OC;
    OC=read_occ(filename,2443846.0);
    print_matrix(OC->E,1,3);
    print_matrix(OC->V,1,3);
    printf("%d %d %d\n",OC->noc,OC->nobs[0],OC->ntotal);
    print_matrix(OC->data[0],18,4);
    print_matrix(OC->TIME[0],18,2);
}
*/
