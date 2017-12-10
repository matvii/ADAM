#include"utils.h"
#include"globals.h"
int INI_CHECKFIT=0;
char *INI_INPUTOBJSHAPE=NULL;
char *INI_INPUTSHAPE=NULL;
int INI_VERBOSE=0;
double *INI_SET_AO_ALBEDO=NULL;
void check_model_fit(LCstruct *LC,AOstruct *AO,OCstruct *OC,RDstruct *RD,CNTRstruct *CR);
int main(int argc, char *argv[])
{
    int argcount=1;
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
    if(INI_LMAX>0)
       fit_oct_model_to_LC_AO(INI_LC,INI_AO,INI_OC,INI_RD,INI_CNTR);
    else
       fit_subdiv_model_to_LC_AO(INI_LC,INI_AO,INI_OC,INI_RD,INI_CNTR);
}
    
        
