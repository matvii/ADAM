
#include"utils.h"
struct CNTR *read_contour(char *filename,double min_tim,int type)
{
    /*
     * Read contours from a file
     * min_tim is substracted from observation times
     * type=0 if contours are given as x-y coordinates
     * type=1 if contours are radius-angle
     * ARE TIMES LT-CORRECTED??
     */
    struct CNTR *C;
    C=malloc(sizeof(struct CNTR));
    FILE *fid;
    int ncont;
    fid=fopen(filename,"r");
    double deg2rad=PI/180.0;
    double angle,radius;
     char delims[]=" \t\r\n\f\v,";
     char *token;
     char *buffer,*filebuff;
     double time,tilt,pixscale;
     double scale;
     double lE,lE0;
     double x,y;
     int nobs;
     buffer=malloc(10000);
     double E[3],E0[3];
     int count=0;
     if(fid==NULL)
    {
        perror("Error opening the contour file\n");
        exit(-1);
    }
    fgets(buffer,2048,fid);
     if(sscanf(buffer,"%d",&ncont)==0) /*Total number of contours*/
  {
    fprintf(stderr,"Error reading contour file (total number of contours)!");
      exit(-1);
  }
  (*C).ncont=ncont;
  (*C).nobs=calloc(ncont,sizeof(int));
  (*C).datax=calloc(ncont,sizeof(double*));
  (*C).datay=calloc(ncont,sizeof(double*));
  (*C).TIME=calloc(ncont,sizeof(double));
  (*C).E=calloc(3*ncont,sizeof(double));
  (*C).E0=calloc(3*ncont,sizeof(double));
  (*C).up=calloc(3*ncont,sizeof(double));
  (*C).distance=calloc(ncont,sizeof(double));
  //Loop over contours
  for(int j=0;j<ncont;j++)
  {
    fgets(buffer,2048,fid);
     
      if(sscanf(buffer,"%lf %lf %lf",E0,E0+1,E0+2)!=3)
          {
            fprintf(stderr,"Error reading E0 in the contour file!\n");
            exit(-1);
          }
      fgets(buffer,2048,fid);    
    if(sscanf(buffer,"%lf %lf %lf",E,E+1,E+2)!=3)
          {
            fprintf(stderr,"Error reading E in the contour file!\n");
            exit(-1);
        }
       fgets(buffer,2048,fid);
    
    if(sscanf(buffer,"%lf %lf",&time,&tilt)!=2)
          {
            fprintf(stderr,"Error reading (time,tilt) in the contour file!\n");
            exit(-1);
        }
         fgets(buffer,2048,fid);
        if(sscanf(buffer,"%lf",&pixscale)!=1)
          {
            fprintf(stderr,"Error reading pixscale in the contour file!\n");
            exit(-1);
        }
         fgets(buffer,2048,fid);
        if(sscanf(buffer,"%d",&nobs)!=1)
          {
            fprintf(stderr,"Error reading the number of points in the contour file!\n");
            exit(-1);
        }
        lE=NORM(E);
        lE0=NORM(E0);
        (*C).datax[j]=calloc(nobs,sizeof(double));
        (*C).datay[j]=calloc(nobs,sizeof(double));
        (*C).nobs[j]=nobs;
        (*C).E[3*j]=E[0]/lE;
        (*C).E[3*j+1]=E[1]/lE;
        (*C).E[3*j+2]=E[2]/lE;
        (*C).E0[3*j]=E0[0]/lE0;
        (*C).E0[3*j+1]=E0[1]/lE0;
        (*C).E0[3*j+2]=E0[2]/lE0;
        (*C).up[3*j]=0;
        (*C).up[3*j+1]=0;
        (*C).up[3*j+2]=0;
        (*C).TIME[j]=time-min_tim;
        (*C).distance[j]=lE;
        count=count+nobs;
        //Now we will read the actual contour points
        
        //Scale to km
        scale=1.0/(lE*pixscale*149597871)*180/PI*3600;
        
        if(type==0) //We are reading x-y coordinate pairs
        {
        for(int k=0;k<nobs;k++)
        {
            fgets(buffer,2048,fid);
            
           if(sscanf(buffer,"%lf %lf",&x,&y)!=2)
          {
            fprintf(stderr,"Error parsing contour point  in the contour file!(contour %d, index %d)\n",j,k);
            exit(-1);
          }
          (*C).datax[j][k]=(cos(tilt*deg2rad)*x-sin(tilt*deg2rad)*y)/scale;
          (*C).datay[j][k]=(cos(tilt*deg2rad)*y+sin(tilt*deg2rad)*x)/scale;
        }
        }
        else
        {
            for(int k=0;k<nobs;k++)
        {
            fgets(buffer,2048,fid);
           if(sscanf(buffer,"%lf %lf",&angle,&radius)!=2)
          {
            fprintf(stderr,"Error parsing (angle,radius) point  in the contour file!(contour %d, index %d)\n",j,k);
            exit(-1);
          }
          x=radius*cos(angle); //IN RADIANS?
          y=radius*sin(angle);
          (*C).datax[j][k]=(cos(tilt*deg2rad)*x-sin(tilt*deg2rad)*y)/scale;
          (*C).datay[j][k]=(cos(tilt*deg2rad)*y+sin(tilt*deg2rad)*x)/scale;
        }
        }
            
  }
  free(buffer);
  fclose(fid);
  (*C).ntotal=count;
  return C;
}
     
int main()
{
    char filename[]="/tmp/Hermione_AO_contours";
    struct CNTR *C;
    
    C=read_contour(filename,0,1);
    int last=C->ncont-1;
    printf("\n");
   for(int j=0;j<10;j++)
       printf(" %4.2f \n",(*C).datax[last][j]);
   printf("\n");
   for(int j=0;j<10;j++)
       printf(" %4.2f \n",(*C).datay[last][j]);

    printf("\n");
   for(int j=0;j<3;j++)
       printf(" %4.2f \n",(*C).E[3*last+j]);
    printf("\n");
   for(int j=0;j<3;j++)
       printf(" %4.2f \n",(*C).E0[3*last+j]);
   printf("\n");
   printf("TIME: %7.4f\n",C->TIME[last]);
}

