#include<stdio.h>
#include<stdlib.h>
#include"structs.h"
#include"utils.h"
#include<math.h>

struct LC  *read_lcurve2(char* filename,double min_tim)
{
  struct LC *tlc;
  char *buffer=calloc(2048,sizeof(char));
  tlc=malloc(sizeof(struct LC));
  FILE *fid;
  fid=fopen(filename,"r");
  if(fid==NULL)
  {
      perror("Cannot open lcurve file!");
      exit(-1);
  }
  int cal,nlc,nobs;
  double *br,cumsum;
  double E[3],E0[3];
  double lE,lE0;
  double time;
  int total=0;
  fgets(buffer,2048,fid);
  
  if(sscanf(buffer,"%d",&nlc)==0) /*Total number of lightcurves*/
  {
    fprintf(stderr,"Error reading lcurve file (total number of lightcurves)!");
      exit(-1);
  }
  
  (*tlc).nlc=nlc;
  (*tlc).nobs=malloc(nlc*sizeof(int));
  (*tlc).lcs=malloc(nlc*sizeof(double*));
  (*tlc).E=malloc(nlc*sizeof(double*));
  (*tlc).E0=malloc(nlc*sizeof(double*));
  (*tlc).TIME=malloc(nlc*sizeof(double*));
  (*tlc).rel=calloc(nlc,sizeof(int));
  
  tlc->calib=0;
  for(int i=0;i<nlc;i++)
  {
      fgets(buffer,2048,fid);
     
      if(sscanf(buffer,"%d %d",&nobs,&cal)!=2)
          {
    fprintf(stderr,"Error reading lcurve file!(nobs: %d,cal: %d), lc: %d",nobs,cal,i+1);
      exit(-1);
  }
          /*Only relative for now, remember to fix this*/
      /*TODO: Also remember to remove comments*/
      (*tlc).nobs[i]=nobs;
      (*tlc).lcs[i]=malloc(nobs*sizeof(double));
      (*tlc).E[i]=malloc(3*nobs*sizeof(double));
      (*tlc).E0[i]=malloc(3*nobs*sizeof(double));
      (*tlc).TIME[i]=malloc(nobs*sizeof(double));
      br=malloc(nobs*sizeof(double));
      cumsum=0;
      if(cal==0)
          (*tlc).rel[i]=1;
      for(int j=0;j<nobs;j++)
      {
          fgets(buffer,2048,fid);
         
         if(sscanf(buffer,"%lf %lf %lf %lf %lf %lf %lf %lf",&time,&(br[j]),E0,E0+1,E0+2,E,E+1,E+2)!=8)
             {
    perror("Error reading lcurve file!");
      exit(-1);
             }
          cumsum+=br[j];
          lE=NORM(E);
          lE0=NORM(E0);
          (*tlc).E0[i][3*j]=E0[0]/lE0;
          (*tlc).E0[i][3*j+1]=E0[1]/lE0;
          (*tlc).E0[i][3*j+2]=E0[2]/lE0;
          
          (*tlc).E[i][3*j]=E[0]/lE;
          (*tlc).E[i][3*j+1]=E[1]/lE;
          (*tlc).E[i][3*j+2]=E[2]/lE;
          
          (*tlc).TIME[i][j]=time-min_tim;
      }
      /*Copy relative brightness values to array*/
      /////////////////////////////////////////////////////////7
      //DO NOT USE CALIBRATED
      cal=0;
      //////////////////////////////////////////////////////////
      if(cal==0)
      {
          for(int k=0;k<nobs;k++)
              (*tlc).lcs[i][k]=br[k]*nobs/cumsum;
      }
      else
      {
          for(int k=0;k<nobs;k++)
              (*tlc).lcs[i][k]=br[k];
          tlc->calib=1;
          printf("calibrated!\n");
      }
              
          
      free(br);
      total+=nobs;
      
  }
  fclose(fid);
  tlc->ntotal=total;
  return tlc;
}

int main()
{
    struct LC *tlc;
    tlc=read_lcurve2("juno/Juno_LC2.txt",0);
    printf("Number of lightcurves: %d obs: %d E: %f %f %f\n",(*tlc).nlc,(*tlc).nobs[37],(*tlc).E0[37][0],(*tlc).E0[37][1],(*tlc).E0[37][2]);
    printf("Brightness:\n");
    for(int i=0;i<10;i++)
        printf("%f\n",(*tlc).lcs[44][i]);
    printf("Total number of points: %d\n",tlc->ntotal);
}

