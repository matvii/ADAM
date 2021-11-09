#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include"iniparser/src/iniparser.h"
#include<math.h>
#include"structs.h"
#include"utils.h"
#include"globals.h"
//Build command gcc test_adam_ini.c -o test_adam_ini -Iiniparser/src -Liniparser -liniparser
//Here are global variables
double INI_MAX_RD_ANGLE=90;
int INI_RADAR_FLIP=0;
int     *INI_AO_FLIP;
int     INI_SUBDIV_TYPE_BUTTERFLY=0;
int     INI_IGNORE_AO_ALBEDO=0;
int     INI_HAVE_LC=0;
int     INI_HAVE_AO=0;
int     INI_HAVE_OC=0;
int     INI_HAVE_HF=0;
int     INI_HAVE_RD=0;
int     INI_HAVE_CNTR=0;
int     INI_CNTR_IS_SPARSE=0;
int     INI_CNTR_RAD=0;
int     INI_ALBEDO_FIT_ONLY=0;
int INI_FIT_AO_ALBEDO=0;
int INI_RESTORE_ALBEDO=1;
double INI_SET_RD_ZERO=-10;
double INI_ANGLE_B=NAN;
double INI_ANGLE_L=NAN;
double INI_ANGLE_P=NAN;
double INI_ANGLE_PHI0=0;
int INI_LMAX=0;
int INI_SD_LEVEL=0;
char *INI_ALBEDO_FILE=NULL;
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
int *INI_SPARSE_LCS=NULL;
int INI_SPARSE_LC_NMR=0;
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
double INI_DW_DEC=1.0;
double INI_DW=0;
double INI_OW=0;
double INI_RW=0;
double INI_ZMAX_WEIGHT=0;
double INI_ZMAX=10;
double INI_CHRDW=0;
double INI_CNTR_WEIGHT=0;
double INI_INER_WEIGHT=0;
double INI_LAMBDAINC=10;
double INI_LAMBDADEC=10;
double INI_LAMBDAMAX=1e6;
double INI_MINDEC=0.1;
double *INI_AO_WEIGHT=NULL;
double *INI_OC_WEIGHT=NULL;
double *INI_RD_WEIGHT=NULL;
double *INI_LC_WEIGHTS=NULL;
int *INI_PHASE_MASK=NULL;
double INI_LAMBDA=1;
double INI_RDEXP=0.59;
int INI_AO_SCALING=1;
double *INI_AO_ROTANGLE=NULL;
double *INI_PHASE_PARAMS=NULL;
double *INI_HAPKE=NULL;
AOstruct *INI_AO;
LCstruct *INI_LC;
OCstruct *INI_OC;
RDstruct *INI_RD;
CNTRstruct *INI_CNTR;
double *INI_CHORD_OFFSET=NULL;
int INI_MASK_SET=0;
int *INI_FREE_CHORD_LIST=NULL;
int INI_FREE_CHORD_NMR=0;
int INI_FIX_SHAPE=0;
int INI_FIX_ANGLES=0;
int INI_FIX_PERIOD=0;
int INI_FIX_A1=0;
int INI_STAR_SHAPED=0;
int INI_LC_ARE_RELATIVE=0;
int INI_FIT_ALBEDO=0;
double INI_ALBREGW=1;
double INI_ALBEDO_MAX=1;
double INI_ALBEDO_MIN=0;
double INI_NDCHORD_WEIGHT=1000;
double INI_COM_WEIGHT=0.0;
int INI_LOGEXP=4;
char *INI_WRITE_STATE_FILE=NULL;
char *INI_ALBEDO_OUT_FILE=NULL;
char *OUT_OBJSHAPE_FILE=NULL;
char *INI_RESTORE_STATE=NULL;
char *OUT_LOG_FILE=NULL;
double *INI_AO_REDUCE_ZERO=NULL;
double *INI_DIA_PARAMS=NULL;
int MAX_VERTICES=10000;
int INI_FIX_VERTEX_NBR=0;
int *INI_FIX_VERTEX_LIST;
int INI_VERTEX_NORMAL=0;
int *INI_RD_CENTER;
int *INI_RD_DIM;
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
    int GlobalAOLowFreq=0;
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
    s=iniparser_getstring(ini,"Shape:RestoreState",NULL);
    if(s!=NULL)
    {
        INI_RESTORE_STATE=calloc(strlen(s)+1,sizeof(char));
        strcpy(INI_RESTORE_STATE,s);
    }
    s=iniparser_getstring(ini,"Shape:SDLevel",NULL);
    if(s!=NULL)
        INI_SD_LEVEL=atoi(s);
    s=iniparser_getstring(ini,"Shape:ButterflySubdiv",NULL);
    if(s!=NULL)
        INI_SUBDIV_TYPE_BUTTERFLY=atoi(s);
    
    s=iniparser_getstring(ini,"Shape:LMAX",NULL);
    if(s!=NULL)
        INI_LMAX=atoi(s);
    s=iniparser_getstring(ini,"Shape:AlbedoFitOnly",0);
    if(s!=NULL)
        INI_ALBEDO_FIT_ONLY=atoi(s);
    s=iniparser_getstring(ini,"Shape:InitAlbedoFile",NULL);
    if(s!=NULL)
    {
        int slength=strlen(s);
        INI_ALBEDO_FILE=calloc(slength+1,sizeof(char));
        strcpy(INI_ALBEDO_FILE,s);
    }
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
    s=iniparser_getstring(ini,"Shape:FixShape",NULL);
    if(s!=NULL)
        INI_FIX_SHAPE=atoi(s);
    s=iniparser_getstring(ini,"Shape:StarShaped",NULL);
    if(s!=NULL)
        INI_STAR_SHAPED=atoi(s);
    s=iniparser_getstring(ini,"Shape:FixAngles",NULL);
    if(s!=NULL)
        INI_FIX_ANGLES=atoi(s);
    s=iniparser_getstring(ini,"Shape:FixPeriod",NULL);
    if(s!=NULL)
        INI_FIX_PERIOD=atoi(s);
    s=iniparser_getstring(ini,"Shape:FixA1",NULL);
    if(s!=NULL)
        INI_FIX_A1=atoi(s);
    s=iniparser_getstring(ini,"Shape:FixShapeVertexFile",NULL);
    if(s!=NULL)
    {
        INI_FIX_VERTEX_NBR=read_vector_fileI_alloc(s,&INI_FIX_VERTEX_LIST,MAX_VERTICES);
        if(INI_FIX_VERTEX_NBR<1)
        {
            perror("Failed to process file in Shape:FixShapeVertexFile");
            INI_FIX_VERTEX_NBR=0;
        }
        INI_MASK_SET=1;
       
    }
    s=iniparser_getstring(ini,"Shape:RestoreAlbedo","1");
    INI_RESTORE_ALBEDO=atoi(s);
    if(INI_FIX_SHAPE==1 || INI_FIX_ANGLES==1 ||INI_FIX_A1==1 || INI_STAR_SHAPED==1||INI_FIX_PERIOD==1)
        INI_MASK_SET=1;
    s=iniparser_getstring(ini,"Shape:VertexNormalFit","0");
    INI_VERTEX_NORMAL=atoi(s);
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
    s=iniparser_getstring(ini,"Optimization:DiAParams",NULL);
    if(s!=NULL)
    {
        INI_DIA_PARAMS=calloc(3,sizeof(double));
        INI_DIA_PARAMS[0]=atof(strtok(s,","));
        INI_DIA_PARAMS[1]=atof(strtok(NULL,","));
        INI_DIA_PARAMS[2]=atof(strtok(NULL,","));
        
       
    }
    s=iniparser_getstring(ini,"Optimization:DiAWeightDec","1");
    INI_DW_DEC=atof(s);
    s=iniparser_getstring(ini,"Optimization:OctWeight","20");
    INI_OW=atof(s);
    s=iniparser_getstring(ini,"Optimization:RDWeight","1");
    INI_RW=atof(s);
    s=iniparser_getstring(ini,"Optimization:RDFlip","0");
    INI_RADAR_FLIP=atoi(s);
    s=iniparser_getstring(ini,"Optimization:Lambda","1");
    INI_LAMBDA=atof(s);
    s=iniparser_getstring(ini,"Optimization:LambdaInc","10");
    INI_LAMBDAINC=atof(s);
    s=iniparser_getstring(ini,"Optimization:LambdaDec","10");
    INI_LAMBDADEC=atof(s);
    s=iniparser_getstring(ini,"Optimization:LambdaMax","1000000");
    INI_LAMBDAMAX=atof(s);
    s=iniparser_getstring(ini,"Optimization:MinDec","0.1");
    INI_MINDEC=atof(s);
    s=iniparser_getstring(ini,"Optimization:RDexp","2");
    INI_RDEXP=log(atof(s));
    s=iniparser_getstring(ini,"Optimization:ChordWeight","1");
    INI_CHRDW=atof(s);
    s=iniparser_getstring(ini,"Optimization:AlbRegWeight","1");
    INI_ALBREGW=atof(s);
    s=iniparser_getstring(ini,"Optimization:CNTRWeight","1");
    INI_CNTR_WEIGHT=atof(s);
    //Parse Data
    s=iniparser_getstring(ini,"Optimization:RestrictZCoord",NULL);
        if(s!=NULL)
            INI_ZMAX=atof(s);
    s=iniparser_getstring(ini,"Optimization:RestrictZcoordWeight",NULL);
    if(s!=NULL)
        INI_ZMAX_WEIGHT=atof(s);
    s=iniparser_getstring(ini,"Optimization:COMWeight","0");
    if(s!=NULL)
        INI_COM_WEIGHT=atof(s);
    s=iniparser_getstring(ini,"Optimization:INERWeight","0");
    if(s!=NULL)
        INI_INER_WEIGHT=atof(s);
    s=iniparser_getstring(ini,"Optimization:IgnoreAOAlbedo","0");
    if(s!=NULL)
        INI_IGNORE_AO_ALBEDO=atoi(s);
    s=iniparser_getstring(ini,"Optimization:FitAOAlbedo","0");
    if(s!=NULL)
        INI_FIT_AO_ALBEDO=atoi(s);
    if(INI_FIT_AO_ALBEDO)
        INI_IGNORE_AO_ALBEDO=0;
    s=iniparser_getstring(ini,"Optimization:AOLowFreq",NULL);
    if(s!=NULL)
        GlobalAOLowFreq=atoi(s);
    s=iniparser_getstring(ini,"Optimization:SetRDZero","-10");
    INI_SET_RD_ZERO=atof(s);
    s=iniparser_getstring(ini,"Optimization:MaxRDAngle",NULL);
    if(s!=NULL)
        INI_MAX_RD_ANGLE=atof(s);
    s=iniparser_getstring(ini,"Data:UseLC","1");
    INI_HAVE_LC=atoi(s);
    s=iniparser_getstring(ini,"Data:UseAO","0");
    if(atoi(s)>0)
    {
        INI_HAVE_AO=1;
        nAO=atoi(s);
    }
    s=iniparser_getstring(ini,"Data:UseCNTR","0");
    if(atoi(s)>0)
        INI_HAVE_CNTR=1;
    s=iniparser_getstring(ini,"Data:UseOC","0");
    INI_HAVE_OC=atoi(s);
    s=iniparser_getstring(ini,"Data:UseRD","0");
    if(atoi(s)>0)
    {
        INI_HAVE_RD=1;
        nRD=atoi(s);
    }
    s=iniparser_getstring(ini,"LC:AllLCRelative","0");
    INI_LC_ARE_RELATIVE=atoi(s);
    s=iniparser_getstring(ini,"LC:LCFile",NULL);
    //Prepare LC data
    if(INI_HAVE_LC)
    {
    INI_LC=read_lcurve(s,INI_MIN_TIM);
    //Read LC Weight file
    s=iniparser_getstring(ini,"LC:LCWeightFile",NULL);
    INI_LC_WEIGHTS=calloc(INI_LC->nlc,sizeof(double));
    for(int jk=0;jk<INI_LC->nlc;jk++)
        INI_LC_WEIGHTS[jk]=1.0;
    if(s!=NULL)
        read_weight_file(s,INI_LC_WEIGHTS,INI_LC->nlc);
    }
    s=iniparser_getstring(ini,"LC:PhaseParams",NULL);
    if(s!=NULL)
    {
        
        INI_PHASE_PARAMS=calloc(3,sizeof(double));
        INI_PHASE_PARAMS[0]=atof(strtok(s,","));
        INI_PHASE_PARAMS[1]=atof(strtok(NULL,","));
        INI_PHASE_PARAMS[2]=atof(strtok(NULL,","));
        
       
    }
   
        
    s=iniparser_getstring(ini,"LC:HapkeParams",NULL);
    if(s!=NULL)
    {
        INI_HAPKE=calloc(5,sizeof(double));
        INI_HAPKE[0]=atof(strtok(s,","));
        INI_HAPKE[1]=atof(strtok(NULL,","));
        INI_HAPKE[2]=atof(strtok(NULL,","));
        INI_HAPKE[3]=atof(strtok(NULL,","));
        INI_HAPKE[4]=atof(strtok(NULL,","));
        if(INI_PHASE_PARAMS!=NULL)
        {
          fprintf(stderr,"Both PhaseParams and HapkeParams is set; ignoring PhaseParams.\n");
          INI_PHASE_PARAMS=NULL;
        }
          
    }
    if(INI_PHASE_PARAMS!=NULL)
    {
    s=iniparser_getstring(ini,"LC:FixedParams",NULL);
    if(s!=NULL)
    {
        INI_PHASE_MASK=calloc(4,sizeof(int));
        INI_PHASE_MASK[0]=atoi(strtok(s,","));
        INI_PHASE_MASK[1]=atoi(strtok(NULL,","));
        INI_PHASE_MASK[2]=atoi(strtok(NULL,","));
        INI_PHASE_MASK[3]=atoi(strtok(NULL,","));
         INI_MASK_SET=1;
    }
    }
    if(INI_HAVE_LC==1 && INI_LC->ncalib>0 && INI_PHASE_PARAMS==NULL && INI_HAPKE==NULL)
    {
        fprintf(stderr,"WARNING: There are calibrated lightcurves, but phase params or Hapke are not set. Either set AllLCRelative=1 or PhaseParams=... or HapkeParams\n");
        
    }
    s=iniparser_getstring(ini,"LC:FitAlbedo","0");
    INI_FIT_ALBEDO=atoi(s);
   
    s=iniparser_getstring(ini,"LC:AlbedoMax","1");
    INI_ALBEDO_MAX=atof(s);
    s=iniparser_getstring(ini,"LC:AlbedoMin","0");
    INI_ALBEDO_MIN=atof(s);
    s=iniparser_getstring(ini,"Optimization:FitAOAlbedo","0");
    if(s!=NULL)
        INI_FIT_AO_ALBEDO=atoi(s);
    if(INI_FIT_AO_ALBEDO==1 && INI_FIT_ALBEDO==0)
    {
        fprintf(stderr,"FitAOALbedo is set but LC:FitAlbedo is not. Setting LC:FitAlbedo to 1. Remember also set LC:ALbedoMax and AlbedoMin. to reasonable values\n");
        INI_FIT_ALBEDO=1;
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
    if(s!=NULL)
    {
    ephmfile=calloc(strlen(s)+1,sizeof(char));
    strcpy(ephmfile,s);
    }
    s=iniparser_getstring(ini,"Output:LCOutputFile",NULL);
    if(s!=NULL)
    {
        OUT_LC_FILE=calloc(strlen(s)+1,sizeof(char));
        strcpy(OUT_LC_FILE,s);
    }
     s=iniparser_getstring(ini,"Output:LogFile",NULL);
     if(s!=NULL)
     {
         OUT_LOG_FILE=calloc(strlen(s)+1,sizeof(char));
         strcpy(OUT_LOG_FILE,s);
     }
     s=iniparser_getstring(ini,"Output:ShapeFile",NULL);
    if(s==NULL)
    {
         perror("Warning: Output shapefile is not set");
    }
    OUT_SHAPE_FILE=calloc(strlen(s)+1,sizeof(char));
    strcpy(OUT_SHAPE_FILE,s);
    s=iniparser_getstring(ini,"Output:ShapeObjFile",NULL);
    if(s!=NULL)
    {
        OUT_OBJSHAPE_FILE=calloc(strlen(s)+1,sizeof(char));
        strcpy(OUT_OBJSHAPE_FILE,s);
    }
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
    s=iniparser_getstring(ini,"Output:StateFile",NULL);
    if(s!=NULL)
    {
        INI_WRITE_STATE_FILE=calloc(strlen(s)+1,sizeof(char));
        strcpy(INI_WRITE_STATE_FILE,s);
    }
    if(INI_CHECKFIT==1)
    {
        if(s==NULL)
            perror("Error: checkfit passed, but StateFile not set in the inifile");
        INI_RESTORE_STATE=calloc(strlen(INI_WRITE_STATE_FILE)+1,sizeof(char));
        strcpy(INI_RESTORE_STATE,INI_WRITE_STATE_FILE);
    }
    s=iniparser_getstring(ini,"Output:AlbedoFile",NULL);
    if(s!=NULL)
    {
        INI_ALBEDO_OUT_FILE=calloc(strlen(s)+1,sizeof(char));
        strcpy(INI_ALBEDO_OUT_FILE,s);
    }
    if(INI_HAVE_CNTR)
    {
        s=iniparser_getstring(ini,"CNTR:Type","Cart");
        int type=0;
        int rotate=0;
        if(strcmp(s,"Ang")==0)
            type=1;
        s=iniparser_getstring(ini,"CNTR:Rotate","1");
        rotate=atoi(s);
        s=iniparser_getstring(ini,"CNTR:Sparse","1");
        INI_CNTR_IS_SPARSE=atoi(s);
        s=iniparser_getstring(ini,"CNTR:CRFile",NULL);
        INI_CNTR=read_contour(s,INI_MIN_TIM,type,rotate);
        s=iniparser_getstring(ini,"CNTR:FitAlg","Cart");
        if(strcmp(s,"Rad")==0)
        {
            INI_CNTR_RAD=1;
            INI_CNTR_IS_SPARSE=1;
        }
            
    }
    
    if(nAO>0)
    {
    INI_AO_REDUCE_ZERO=calloc(nAO,sizeof(double));
    INI_AO_ROTANGLE=calloc(nAO,sizeof(double));
    INI_AO_FLIP=calloc(nAO,sizeof(int));
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
            snprintf(sect,20,"AO%d:SetZero",j+1);
            s=iniparser_getstring(ini,sect,NULL);
            if(s!=NULL)
                INI_AO_REDUCE_ZERO[j]=atof(s);
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
            else if(strcmp(s,"Ecc")==0)
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
            snprintf(sect,20,"AO%d:RotAngle",j+1);
            s=iniparser_getstring(ini,sect,NULL);
            if(s!=NULL)
                INI_AO_ROTANGLE[j]=atof(s);
            snprintf(sect,20,"AO%d:FlipHor",j+1);
            s=iniparser_getstring(ini,sect,NULL);
            if(s!=NULL)
                INI_AO_FLIP[j]=atoi(s);
               
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
        
    
    int Large=150;
    if(GlobalAOLowFreq==1)
        for(int jk=0;jk<nAO;jk++)
            LowFreq[jk]=1;
        
    
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
            parse_vector(s,INI_OC_OFFSET,2*(INI_OC->noc));
         
            
           
        } 
        s=iniparser_getstring(ini,"OC:ChordWeight",NULL);
        if(s!=NULL)
        {
            INI_OC_WEIGHT=calloc(INI_OC->ntotal,sizeof(double));
            for(int l=0;l<INI_OC->ntotal;l++)
                INI_OC_WEIGHT[l]=1.0;
             parse_vector(s,INI_OC_WEIGHT,INI_OC->ntotal);
        }
        s=iniparser_getstring(ini,"OC:ChordWeightFile",NULL);
        if(s!=NULL)
        {
            INI_OC_WEIGHT=calloc(INI_OC->ntotal,sizeof(double));
            for(int l=0;l<INI_OC->ntotal;l++)
                INI_OC_WEIGHT[l]=1.0;
            read_vector_file(s,INI_OC_WEIGHT,INI_OC->ntotal);
        }
        s=iniparser_getstring(ini,"OC:FreeChords",NULL);
        if(s!=NULL)
        {
            INI_MASK_SET=1;
            INI_FREE_CHORD_LIST=calloc(INI_OC->ntotal,sizeof(double));
            INI_FREE_CHORD_NMR=parse_vectorI(s,INI_FREE_CHORD_LIST,INI_OC->ntotal);
        }
        s=iniparser_getstring(ini,"OC:ChordOffset",NULL);
        if(s!=NULL)
        {
            INI_CHORD_OFFSET=calloc(INI_OC->ntotal,sizeof(double));
            parse_vector(s,INI_CHORD_OFFSET,INI_OC->ntotal);
        }
         s=iniparser_getstring(ini,"OC:NDChordWeight",NULL);
         if(s!=NULL)
             INI_NDCHORD_WEIGHT=atof(s);
         s=iniparser_getstring(ini,"OC:LogExp",NULL);
         if(s!=NULL)
             INI_LOGEXP=atoi(s);
        
            
            
    }
        if(nRD>0)
    {
    double* RD_iWeight=calloc(nRD,sizeof(double));
    for(int k=0;k<nRD;k++)
        RD_iWeight[k]=1.0;
    int UseRDWeight=0;
    int *RDLowFreq=calloc(nRD,sizeof(int));
    RDfiles=calloc(nRD,sizeof(char*));
    int *cx,*cy,*cxdim,*cydim;
   cx=calloc(nRD,sizeof(int));
   cy=calloc(nRD,sizeof(int));
   cxdim=calloc(nRD,sizeof(int));
   cydim=calloc(nRD,sizeof(int));
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
            snprintf(sect,20,"RD%d:RDcenter",j+1);
            s=iniparser_getstring(ini,sect,NULL);
            if(s!=NULL)
            {
                cx[j]=atoi(strtok(s,","));
                cy[j]=atoi(strtok(NULL,","));
            }
            snprintf(sect,20,"RD%d:RDFromCenter",j+1);
            s=iniparser_getstring(ini,sect,NULL);
            if(s!=NULL)
            {
                cxdim[j]=atoi(strtok(s,","));
                cydim[j]=atoi(strtok(NULL,","));
                //Rectx=cx[j]-cxdim[j]
                //size=2*cxdim[j]
                //[cx[j]-cxdim[j]:2*cxdim[j]+cx[j]-1
            }
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
            snprintf(sect,20,"RD%d:Weight",j+1);
            s=iniparser_getstring(ini,sect,NULL);
            if(s!=NULL)
            {
                UseRDWeight=1;
                RD_iWeight[j]=atof(s);
            }
             
       }
       if(UseRDWeight==1)
        INI_RD_WEIGHT=RD_iWeight;
       else
           free(RD_iWeight);
    
    int Large=500;
    
    
      double *E=calloc(Large*3,sizeof(double));
    double *E0=calloc(Large*3,sizeof(double));
    double *TIME=calloc(Large,sizeof(double));  
    int nephm=read_ephm_data(ephmfile,TIME,E,E0);
    INI_RD=process_rd_images(RDfiles,nRD,x0,y0,nx,ny,cx,cy,cxdim,cydim,pixscalex,pixscaley,date,RadarFreq,INI_MIN_TIM,E,TIME,nephm,RDLowFreq);
    free(E);
    free(E0);
    free(TIME);
    for(int j=0;j<nRD;j++)
    {
    free(RDfiles[j]);
   
    }
    free(cx);
    free(cy);
    free(cxdim);
    free(cydim);
    free(x0);
    free(y0);
    free(nx);
    free(ny);
    free(pixscalex);
    free(pixscaley);
    }       
    iniparser_freedict(ini);
    

}
