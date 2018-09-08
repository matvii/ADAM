#include"utils.h"
#include"globals.h"
int INI_CHECKFIT=0;
char *INI_INPUTOBJSHAPE=NULL;
char *INI_INPUTSHAPE=NULL;
int INI_VERBOSE=0;
double *INI_SET_AO_ALBEDO=NULL;
void fit_model_albedo(LCstruct *LC,AOstruct *AO,double** feAlbedo,int* alblength);
void check_model_fit(LCstruct *LC,AOstruct *AO,OCstruct *OC,RDstruct *RD,CNTRstruct *CR);
int main(int argc, char *argv[])
{
    int argcount=1;
    double *feAlbedo;
    int alblength;
    while(argcount<argc)
    {
    if(strncmp("--checkfit",argv[argcount],10)==0)
        INI_CHECKFIT=1;
    if(strncmp("--objshapefile",argv[argcount],14)==0)
    {
        INI_INPUTOBJSHAPE=calloc(100,sizeof(char));
        strcpy(INI_INPUTOBJSHAPE,argv[argcount+1]);
    }
    if(strncmp("--shapefile",argv[argcount],11)==0)
    {
        INI_INPUTSHAPE=calloc(100,sizeof(char));
        strcpy(INI_INPUTSHAPE,argv[argcount+1]);
    }
    if(strncmp("--AOalbfile",argv[argcount],11)==0)
    {
        int count=0;
      count=read_values_from_file(argv[argcount+1],&INI_SET_AO_ALBEDO);
      printf("Read %d facet albedo values\n",count);
    }
    argcount++;
    }
   
        
    if(argc==1)
        parse_ini("adam.ini");
    else
        parse_ini(argv[1]);
   
   if(INI_CHECKFIT==1)
       {
           check_model_fit(INI_LC,INI_AO,INI_OC,INI_RD,INI_CNTR);
           exit(0);
       }
    if(INI_VERTEX_NORMAL==1)
    {
        printf("Vertex normal fitting enabled\n");
        fit_vertex_model_to_LC_AO(INI_LC,INI_AO,INI_OC,INI_RD,INI_CNTR);
        exit(0);
    }
    double Alblimits[2];
    Alblimits[0]=INI_ALBEDO_MIN;
    Alblimits[1]=INI_ALBEDO_MAX;
    
    if(INI_ALBEDO_FIT_ONLY==1)
    {
        printf("Albedo only fitting\n");
        INI_FIT_ALBEDO=1;
        INI_FIT_AO_ALBEDO=1;
        fit_model_albedo(INI_LC,INI_AO,&feAlbedo,&alblength);
        if(INI_ALBEDO_OUT_FILE!=NULL)
            {
                double *falbedo=calloc(alblength,sizeof(double));
                for(int j=0;j<alblength;j++)
                    falbedo[j]=(Alblimits[0]+Alblimits[1])/2.0+(Alblimits[1]-Alblimits[0])/2.0*tanh(feAlbedo[j]);
                write_matrix_file(INI_ALBEDO_OUT_FILE,falbedo,1,alblength);
                free(falbedo);
            }
            else
                fprintf(stderr,"Albedo file unset\n");
            exit(0);
    }
    if(INI_LMAX>0)
       fit_oct_model_to_LC_AO(INI_LC,INI_AO,INI_OC,INI_RD,INI_CNTR);
    else
       fit_subdiv_model_to_LC_AO(INI_LC,INI_AO,INI_OC,INI_RD,INI_CNTR);
}
    
        
