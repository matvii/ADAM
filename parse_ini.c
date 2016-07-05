#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include"iniparser.h"
#include<math.h>
#include"structs.h"
#include"utils.h"
#include"globals.h"
//Build command gcc test_adam_ini.c -o test_adam_ini -Iiniparser/src -Liniparser -liniparser
//Here are global variables
int     INI_HAVE_LC=0;
int     INI_HAVE_AO=0;
int     INI_HAVE_OC=0;
int     INI_HAVE_HF=0;
int     INI_HAVE_RD=0;
double INI_ANGLE_B=NAN;
double INI_ANGLE_L=NAN;
double INI_ANGLE_P=NAN;
double INI_ANGLE_PHI0=0;
int INI_LMAX=0;
int INI_SD_LEVEL=0;
char *INI_SHAPE_FILE=NULL;
char *INI_INPUT_AO_OFFSET=NULL;
char *INI_OUTPUT_AO_OFFSET=NULL;
char *INI_INPUT_RD_OFFSET=NULL;
char *INI_OUTPUT_RD_OFFSET=NULL;
char *OUT_SHAPE_FILE=NULL;
char *OUT_PARAM_FILE=NULL;
char *OUT_SHAPE_PARAM_FILE=NULL;
double *INI_OC_OFFSET=NULL;
char *OUT_LC_FILE=NULL;
int USE_ELLIPSOID=0;
double ELLIPSOID_SEMI_A=0;
double ELLIPSOID_SEMI_B=0;
double ELLIPSOID_SEMI_C=0;
double INI_MIN_TIM=0;
int INI_NROWS=0;
int NUM_OF_ROUNDS=0;
double INI_LCW=0;
double INI_AOW=0;
double INI_OCW=0;
double INI_CW=0;
double INI_AW=0;
double INI_DW=0;
double INI_OW=0;
double INI_RW=0;
double *INI_AO_WEIGHT=NULL;
double *INI_OC_WEIGHT=NULL;
double *INI_RD_WEIGHT=NULL;
double INI_LAMBDA=1;
double INI_RDEXP=0.59;
int INI_AO_SCALING=0;
double *INI_PARAMS;
AOstruct *INI_AO;
LCstruct *INI_LC;
OCstruct *INI_OC;
RDstruct *INI_RD;
int parse_ini(char *filename)
{
    char *ephmfile;
    int nAO=0;
    int nRD=0;
    char *OCfile;
    char **RDfiles;
    double *AOWeight;
    int UseAOWeight=0;
    dictionary *ini;
    ini=iniparser_load(filename);
    
    if(ini==NULL)
    {
        perror("Cannot load ini file\n");
        exit(-1);
    }
     char *s;
    char *results;
    s=iniparser_getstring(ini,"Shape:InitEllipsoid",NULL);
    if(s!=NULL)
    {
        USE_ELLIPSOID=1;
        ELLIPSOID_SEMI_A=atof(strtok(s,","));
        ELLIPSOID_SEMI_B=atof(strtok(NULL,","));
        ELLIPSOID_SEMI_C=atof(strtok(NULL,","));
    }
    s=iniparser_getstring(ini,"Shape:Nrows",NULL);
    if(s!=NULL)
        INI_NROWS=atoi(s);
    s=iniparser_getstring(ini,"Shape:InitShapeFile",NULL);
    if(s!=NULL)
    {
        int slength=strlen(s);
        INI_SHAPE_FILE=calloc(slength+1,sizeof(char));
        strcpy(INI_SHAPE_FILE,s);
    }
    s=iniparser_getstring(ini,"Shape:SDLevel",NULL);
    if(s!=NULL)
        INI_SD_LEVEL=atoi(s);
    s=iniparser_getstring(ini,"Shape:LMAX",NULL);
    if(s!=NULL)
        INI_LMAX=atoi(s);
    s=iniparser_getstring(ini,"Shape:MinTim",NULL);
    if(s!=NULL)
        INI_MIN_TIM=atof(s);
    else
        printf("Zero time not set. Assuming 0\n");
    s=iniparser_getstring(ini,"Shape:Angles",NULL);
    if(s!=NULL)
    {
        INI_ANGLE_B=atof(strtok(s,","));
        INI_ANGLE_L=atof(strtok(NULL,","));
        INI_ANGLE_P=atof(strtok(NULL,","));
        if(fabs(INI_ANGLE_P)<1e-10)
        {
            fprintf(stderr,"Rotation period is initialized to zero (or very small). Check that Angles are properly set in the ini file\n");
            exit(1);
        }
        char *phi0=strtok(NULL,",");
        if(phi0!=NULL)
            INI_ANGLE_PHI0=atof(phi0);
    }
    else
    {
        fprintf(stderr,"Rotation angles must be set in the ini file\n");
        exit(1);
    }
    //Parse optimization
    s=iniparser_getstring(ini,"Optimization:NumberofRounds","50");
    NUM_OF_ROUNDS=atoi(s);
    s=iniparser_getstring(ini,"Optimization:UseAOScaling","0");
    INI_AO_SCALING=atoi(s);
    s=iniparser_getstring(ini,"Optimization:LCWeight","2");
    INI_LCW=atof(s);
    s=iniparser_getstring(ini,"Optimization:AOWeight","2");
    INI_AOW=atof(s);
    s=iniparser_getstring(ini,"Optimization:OCWeight","1");
    INI_OCW=atof(s);
    s=iniparser_getstring(ini,"Optimization:ConvexWeight","20");
    INI_CW=atof(s);
    s=iniparser_getstring(ini,"Optimization:AreaWeight","20");
    INI_AW=atof(s);
    s=iniparser_getstring(ini,"Optimization:DiAWeight","2");
    INI_DW=atof(s);
    s=iniparser_getstring(ini,"Optimization:OctWeight","20");
    INI_OW=atof(s);
    s=iniparser_getstring(ini,"Optimization:RDWeight","1");
    INI_RW=atof(s);
    s=iniparser_getstring(ini,"Optimization:Lambda","1");
    INI_LAMBDA=atof(s);
    s=iniparser_getstring(ini,"Optimization:RDexp","2");
    INI_RDEXP=log(atoi(s));
    //Parse Data
    s=iniparser_getstring(ini,"Data:UseLC","1");
    INI_HAVE_LC=atoi(s);
    s=iniparser_getstring(ini,"Data:UseAO","0");
    if(atoi(s)>0)
    {
        INI_HAVE_AO=1;
        nAO=atoi(s);
    }
    s=iniparser_getstring(ini,"Data:UseOC","0");
    INI_HAVE_OC=atoi(s);
    s=iniparser_getstring(ini,"Data:UseRD","0");
    if(atoi(s)>0)
    {
        INI_HAVE_RD=1;
        nRD=atoi(s);
    }
    s=iniparser_getstring(ini,"LC:LCFile",NULL);
    //Prepare LC data
    INI_LC=read_lcurve(s,INI_MIN_TIM);
    s=iniparser_getstring(ini,"LC:params",NULL);
    if(s!=NULL)
    {
        INI_PARAMS=calloc(4,sizeof(double));
        INI_PARAMS[0]=atof(strtok(s,","));
        INI_PARAMS[1]=atof(strtok(NULL,","));
        INI_PARAMS[2]=atof(strtok(NULL,","));
        INI_PARAMS[3]=atof(strtok(NULL,","));
    }
    //If AO data exists, read the images
    char **AOfiles;
    char **PSFfiles;
    int *x0;
    int *y0;
    int *nx;
    int *ny;
    int *px0;
    int *py0;
    int *npx;
    int *npy;
    double *date;
    double *pixscalex,*pixscaley;
    double *up;
    double OCup[3];
    //Ephm section
    s=iniparser_getstring(ini,"Ephm:EphFile",NULL);
    if(s==NULL)
    {
         perror("Ephemeris information required. Set file in Ephm");
        exit(-1);
    }
    ephmfile=calloc(strlen(s)+1,sizeof(char));
    strcpy(ephmfile,s);
   
    s=iniparser_getstring(ini,"Output:LCOutputFile",NULL);
    if(s!=NULL)
    {
        OUT_LC_FILE=calloc(strlen(s)+1,sizeof(char));
        strcpy(OUT_LC_FILE,s);
    }
        
     s=iniparser_getstring(ini,"Output:ShapeFile",NULL);
    if(s==NULL)
    {
         perror("Warning: Output shapefile is not set");
    }
    OUT_SHAPE_FILE=calloc(strlen(s)+1,sizeof(char));
    strcpy(OUT_SHAPE_FILE,s);
     s=iniparser_getstring(ini,"Output:AnglesFile",NULL);
    if(s==NULL)
    {
         perror("Warning: Output angles file is not set");
    }
    OUT_PARAM_FILE=calloc(strlen(s)+1,sizeof(char));
    strcpy(OUT_PARAM_FILE,s);
    s=iniparser_getstring(ini,"Output:ShapeParamFile",NULL);
    if(s!=NULL)
    {
    OUT_SHAPE_PARAM_FILE=calloc(strlen(s)+1,sizeof(char));
    strcpy(OUT_SHAPE_PARAM_FILE,s);
    }
    
    
    if(nAO>0)
    {
    int *LowFreq=calloc(nAO,sizeof(int));
    AOfiles=calloc(nAO,sizeof(char*));
    AOWeight=calloc(nAO,sizeof(double));
    PSFfiles=calloc(nAO,sizeof(char*));
    x0=calloc(nAO,sizeof(int));
    y0=calloc(nAO,sizeof(int));
    nx=calloc(nAO,sizeof(int));
    ny=calloc(nAO,sizeof(int));
    px0=calloc(nAO,sizeof(int));
    py0=calloc(nAO,sizeof(int));
    npx=calloc(nAO,sizeof(int));
    npy=calloc(nAO,sizeof(int));
    date=calloc(nAO,sizeof(double));
    pixscalex=calloc(nAO,sizeof(double));
    pixscaley=calloc(nAO,sizeof(double));
    up=calloc(3*nAO,sizeof(double));
    for(int k=0;k<nAO;k++)
        AOWeight[k]=1.0;
    s=iniparser_getstring(ini,"Optimization:AOOffsetFile",NULL);
    if(s!=NULL)
    {
        INI_INPUT_AO_OFFSET=calloc(strlen(s)+1,sizeof(char));
        strcpy(INI_INPUT_AO_OFFSET,s);
    }
    s=iniparser_getstring(ini,"Output:AOOffsetFile",NULL);
    if(s!=NULL)
    {
        INI_OUTPUT_AO_OFFSET=calloc(strlen(s)+1,sizeof(char));
        strcpy(INI_OUTPUT_AO_OFFSET,s);
    }
        char sect[20];
        for(int j=0;j<nAO;j++)
        {
            snprintf(sect,20,"AO%d:AOFile",j+1); //AO file name
            
            s=iniparser_getstring(ini,sect,NULL);
            
            AOfiles[j]=calloc(strlen(s)+1,sizeof(char));
            strcpy(AOfiles[j],s);
            
            snprintf(sect,20,"AO%d:PSFFile",j+1); 
            s=iniparser_getstring(ini,sect,NULL);
            if(s!=NULL)
            {
            PSFfiles[j]=calloc(strlen(s)+1,sizeof(char));
            strcpy(PSFfiles[j],s);
            }
             snprintf(sect,20,"AO%d:AORect",j+1); 
            s=iniparser_getstring(ini,sect,NULL);
            if(s!=NULL)
            {
            x0[j]=atoi(strtok(s,","));
            y0[j]=atoi(strtok(NULL,","));
            }
            snprintf(sect,20,"AO%d:AOSize",j+1); 
            s=iniparser_getstring(ini,sect,NULL);
            if(s!=NULL)
            {
            nx[j]=atoi(strtok(s,","));
            ny[j]=atoi(strtok(NULL,","));
            }
            snprintf(sect,20,"AO%d:Date",j+1); 
            s=iniparser_getstring(ini,sect,NULL);
           
            if(s!=NULL)
            {
           date[j]=atof(s);
            }
            else
                date[j]=NAN;
            
             snprintf(sect,20,"AO%d:PSFRect",j+1); 
            s=iniparser_getstring(ini,sect,NULL);
            if(s!=NULL)
            {
            px0[j]=atoi(strtok(s,","));
            py0[j]=atoi(strtok(NULL,","));
            }
            snprintf(sect,20,"AO%d:PSFSize",j+1); 
            s=iniparser_getstring(ini,sect,NULL);
            if(s!=NULL)
            {
            npx[j]=atoi(strtok(s,","));
            npy[j]=atoi(strtok(NULL,","));
            }
             snprintf(sect,20,"AO%d:PixScale",j+1); 
            s=iniparser_getstring(ini,sect,NULL);
            
            if(s==NULL)
            {
                perror("Pixel scale must be set in ini!");
                exit(-1);
            }
            
            pixscalex[j]=atof(strtok(s,","));
            pixscaley[j]=atof(strtok(NULL,","));
            snprintf(sect,20,"AO%d:SetCamera",j+1);
             
            s=iniparser_getstring(ini,sect,"Unset");
            
            if(strcmp(s,"Equ")==0)
            {
                up[3*j+0]=0;
                up[3*j+1]=0.397748474527011;
                up[3*j+2]=0.917494496447491;
            }
            else if(strcmp(s,"Equ")==0)
            {
                up[3*j+0]=0;
                up[3*j+1]=0;
                up[3*j+2]=1;
            }
            else
            {
               up[3*j+0]=0;
                up[3*j+1]=0;
                up[3*j+2]=0;
                
            }
            snprintf(sect,20,"AO%d:SetCamUp",j+1);
            s=iniparser_getstring(ini,sect,NULL);
            if(s!=NULL)
            {
                up[3*j+0]=atof(strtok(s,","));
                up[3*j+1]=atof(strtok(NULL,","));
                up[3*j+2]=atof(strtok(NULL,","));
            }
           snprintf(sect,20,"AO%d:LowFreq",j+1);
            s=iniparser_getstring(ini,sect,NULL);
            if(s!=NULL)
            LowFreq[j]=atoi(s);
            snprintf(sect,20,"AO%d:Weight",j+1);
            s=iniparser_getstring(ini,sect,NULL);
            if(s!=NULL)
            {
                UseAOWeight=1;
                AOWeight[j]=atof(s);
            }
        }
        if(UseAOWeight==1)
            INI_AO_WEIGHT=AOWeight;
        else
            free(AOWeight);
        
    
    int Large=100;
    
    
      double *E=calloc(Large*3,sizeof(double));
    double *E0=calloc(Large*3,sizeof(double));
    double *TIME=calloc(Large,sizeof(double));  
    int nephm=read_ephm_data(ephmfile,TIME,E,E0);
    INI_AO=process_ao_images(AOfiles,PSFfiles,nAO,x0,y0,nx,ny,px0,py0,npx,npy,pixscalex,pixscaley,date,INI_MIN_TIM,E,E0,up,TIME,nephm,LowFreq);
    free(E);
    free(E0);
    free(TIME);
    for(int j=0;j<nAO;j++)
    {
    free(AOfiles[j]);
    if(PSFfiles[j]!=NULL)
    free(PSFfiles[j]);
    }
    free(x0);
    free(y0);
    free(nx);
    free(ny);
    free(pixscalex);
    free(pixscaley);
    }
    if(INI_HAVE_OC>0)
    {
        OCfile=iniparser_getstring(ini,"OC:OCFile",NULL);
        if(OCfile==NULL)
        {
            perror("UseOC set but filename is not!");
            exit(-1);
        }
        s=iniparser_getstring(ini,"OC:SetCamera",NULL);   
        if(strcmp(s,"Equ")==0)
        {
            OCup[0]=0;
            OCup[1]=0.397748474527011;
            OCup[2]=0.917494496447491;
        }
        else
        {
            OCup[0]=0;
            OCup[1]=0;
            OCup[2]=1;
        }
        s=iniparser_getstring(ini,"OC:SetCamUp",NULL);
        if(s!=NULL)
        {
            OCup[0]=atof(strtok(s,","));
            OCup[1]=atof(strtok(NULL,","));
            OCup[2]=atof(strtok(NULL,","));
        }  
        INI_OC=read_occ(OCfile,INI_MIN_TIM,OCup);
        s=iniparser_getstring(ini,"OC:OCCOffset",NULL);
        if(s!=NULL)
        {
            INI_OC_OFFSET=calloc(2*(INI_OC->noc),sizeof(double));
            INI_OC_OFFSET[0]=atof(strtok(s,","));
            for(int k=1;k<2*(INI_OC->noc);k++)
            INI_OC_OFFSET[k]=atof(strtok(NULL,","));
         
            
           
        }  
        s=iniparser_getstring(ini,"OC:ChordWeightFile",NULL);
        if(s!=NULL)
        {
            INI_OC_WEIGHT=calloc(INI_OC->ntotal,sizeof(double));
            for(int l=0;l<INI_OC->ntotal;l++)
                INI_OC_WEIGHT[l]=1.0;
            read_vector_file(s,INI_OC_WEIGHT,INI_OC->ntotal);
        }
        
            
            
    }
        if(nRD>0)
    {
    int *RDLowFreq=calloc(nRD,sizeof(int));
    RDfiles=calloc(nRD,sizeof(char*));
   
    //PSFfiles=calloc(nRD,sizeof(char*));
    x0=calloc(nRD,sizeof(int));
    y0=calloc(nRD,sizeof(int));
    nx=calloc(nRD,sizeof(int));
    ny=calloc(nRD,sizeof(int));
    //px0=calloc(nAO,sizeof(int));
    //py0=calloc(nAO,sizeof(int));
    //npx=calloc(nAO,sizeof(int));
    //npy=calloc(nAO,sizeof(int));
    date=calloc(nRD,sizeof(double));
    pixscalex=calloc(nRD,sizeof(double));
    pixscaley=calloc(nRD,sizeof(double));
    double *RadarFreq=calloc(nRD,sizeof(double));
    //up=calloc(3*nAO,sizeof(double));
    s=iniparser_getstring(ini,"Optimization:RDOffsetFile",NULL);
    if(s!=NULL)
    {
        INI_INPUT_RD_OFFSET=calloc(strlen(s)+1,sizeof(char));
        strcpy(INI_INPUT_RD_OFFSET,s);
    }
    s=iniparser_getstring(ini,"Output:RDOffsetFile",NULL);
    if(s!=NULL)
    {
        INI_OUTPUT_RD_OFFSET=calloc(strlen(s)+1,sizeof(char));
        strcpy(INI_OUTPUT_RD_OFFSET,s);
    }
        char sect[20];
        for(int j=0;j<nRD;j++)
        {
            snprintf(sect,20,"RD%d:RDFile",j+1); //AO file name
            
            s=iniparser_getstring(ini,sect,NULL);
            
            RDfiles[j]=calloc(strlen(s)+1,sizeof(char));
            strcpy(RDfiles[j],s);
            
           // snprintf(sect,20,"AO%d:PSFFile",j+1); 
           // s=iniparser_getstring(ini,sect,NULL);
           // if(s!=NULL)
          //  {
          //  PSFfiles[j]=calloc(strlen(s)+1,sizeof(char));
         //   strcpy(PSFfiles[j],s);
          //  }
             snprintf(sect,20,"RD%d:RDRect",j+1); 
            s=iniparser_getstring(ini,sect,NULL);
            if(s!=NULL)
            {
            x0[j]=atoi(strtok(s,","));
            y0[j]=atoi(strtok(NULL,","));
            }
            snprintf(sect,20,"RD%d:RDSize",j+1); 
            s=iniparser_getstring(ini,sect,NULL);
            if(s!=NULL)
            {
            nx[j]=atoi(strtok(s,","));
            ny[j]=atoi(strtok(NULL,","));
            }
            snprintf(sect,20,"RD%d:Date",j+1); 
            s=iniparser_getstring(ini,sect,NULL);
           
            if(s!=NULL)
            {
           date[j]=atof(s);
            }
            else
                date[j]=NAN;
            
        
             snprintf(sect,20,"RD%d:PixScale",j+1); 
            s=iniparser_getstring(ini,sect,NULL);
            
            if(s==NULL)
            {
                fprintf(stderr,"Pixel scale of radar image is not set in ini, trying to read from fits\n");
                pixscalex[j]=0; //THIS IS FREQUENCY SPACING
                pixscaley[j]=0; //PIXEL SIZE IN us, delay
            }
            else
            {
            pixscalex[j]=atof(strtok(s,","));
            pixscaley[j]=atof(strtok(NULL,","));
            }
            
      snprintf(sect,20,"RD%d:RadarFreq",j+1); 
            s=iniparser_getstring(ini,sect,NULL);
            
            if(s==NULL)
            {
                fprintf(stderr,"Radar frequency of radar image must be set in ini\n");
                exit(-1);
                
            }
            else
            {
            RadarFreq[j]=atof(strtok(s,","));
            
            }
            snprintf(sect,20,"RD%d:LowFreq",j+1);
            s=iniparser_getstring(ini,sect,NULL);
            if(s!=NULL)
            RDLowFreq[j]=atoi(s);
             
       }
        
    
    int Large=100;
    
    
      double *E=calloc(Large*3,sizeof(double));
    double *E0=calloc(Large*3,sizeof(double));
    double *TIME=calloc(Large,sizeof(double));  
    int nephm=read_ephm_data(ephmfile,TIME,E,E0);
    INI_RD=process_rd_images(RDfiles,nRD,x0,y0,nx,ny,pixscalex,pixscaley,date,RadarFreq,INI_MIN_TIM,E,TIME,nephm,RDLowFreq);
    free(E);
    free(E0);
    free(TIME);
    for(int j=0;j<nRD;j++)
    {
    free(RDfiles[j]);
   
    }
    free(x0);
    free(y0);
    free(nx);
    free(ny);
    free(pixscalex);
    free(pixscaley);
    }       
    iniparser_freedict(ini);
    

}
