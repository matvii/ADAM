#include"utils.h"
#include"structs.h"
void write_lcs(char *lcfile,double *LCout,char *outputlc)
{
    FILE *fid,*fidout;
    fid=fopen(lcfile,"r");
    fidout=fopen(outputlc,"w");
    int rel,nlc,nobs;
  double br;
  double E[3],E0[3];
  double lE,lE0;
  double time;
  int totalcount=0;
  int total=0;
  if(fscanf(fid,"%d",&nlc)==0) /*Total number of lightcurves*/
  {
    perror("Error reading lcurve file!");
      exit(-1);
  }
  fprintf(fidout,"%d\n",nlc);
  for(int i=0;i<nlc;i++)
  {
      if(fscanf(fid,"%d %d",&nobs,&rel)!=2)
      {
          perror("Error reading lcurve file!");
          exit(-1);
      }
        fprintf(fidout,"%d %d\n",nobs,rel);
          /*Only relative for now, remember to fix this*/
      /*TODO: Also remember to remove comments*/
      
      for(int j=0;j<nobs;j++)
      {
         if(fscanf(fid,"%lf %lf %lf %lf %lf %lf %lf %lf",&time,&br,E0,E0+1,E0+2,E,E+1,E+2)!=8)
             {
    perror("Error reading lcurve file!");
      exit(-1);
             }
             
         fprintf(fidout,"%lf %lf %lf %lf %lf %lf %lf %lf\n",time,-LCout[totalcount++],E0[0],E0[1],E0[2],E[0],E[1],E[2]);
          
      }
      /*Copy relative brightness values to array*/
      
      
  }
  fclose(fid);
  fclose(fidout);
}
void main(int argc, char *argv[])
{
    /*
     * Usage:
     * -l lc.txt -a beta lambda omega -t T0 -s shapefile -o outputlc
     */
    int argcount=1;
    char *arg;
    char *lcfile=calloc(100,sizeof(char));
    char *shapefile=calloc(100,sizeof(char));
    char *outputlc=calloc(100,sizeof(char));
    double angles[4];
    double T0=0;
    LCstruct *LC;
    
    
    while(argcount<=argc)
    {
        if(strcmp(argv[argcount],"-l")==0)
        {
            strcpy(lcfile,argv[argcount+1]);
            argcount=argcount+2;
            continue;
        }
        if(strcmp(argv[argcount],"-a")==0)
        {
            angles[0]=(90-atof(argv[argcount+1]))*PI/180;
            angles[1]=atof(argv[argcount+2])*PI/180;
            angles[2]=24*2*PI/atof(argv[argcount+3]);
            angles[3]=0;
            argcount=argcount+4;
            continue;
        }
        if(strcmp(argv[argcount],"-t")==0)
        {
            T0=atof(argv[argcount+1]);
            argcount=argcount+2;
            continue;
        }
        if(strcmp(argv[argcount],"-s")==0)
        {
            strcpy(shapefile,argv[argcount+1]);
            argcount=argcount+2;
            continue;
        }
        if(strcmp(argv[argcount],"-o")==0)
        {
            strcpy(outputlc,argv[argcount+1]);
            argcount=argcount+2;
            continue;
        }
        fprintf(stderr,"Unknown option %s, exiting\n",argv[argcount]);
        fprintf(stderr,"Usage: -l lc.txt -a beta lambda P -t T0 -s shapefile -o outputlc\n");
        exit(-1);
    }
    LC=read_lcurve(lcfile,T0);
    int *tlist,nfac,nvert;
    double *vlist;
    int nLCtotal=LC->ntotal;
    double *LCout=calloc(nLCtotal,sizeof(double));
    read_shape(shapefile,&tlist,&vlist,&nfac,&nvert,0);
    for(int j=0;j<LC->ntotal;j++)
        zero_array(LC->lcs[j],LC-nobs[j]);
    calculate_lcs(tlist,vlist,nfac,nvert,angles,LC,NULL,nvert,nvert,LCout,NULL,NULL,NULL,NULL,0);
    write_lcs(lcfile,LCout,outputlc);
}
            
            
            
    
        
