#include"utils.h"
#include"globals.h"

int main(int argc, char *argv[])
{
    char *inifile;
    if(argc==1)
        parse_ini("adam.ini");
    else
        parse_ini(argv[1]);
   
    //If INI_LMAX is set, we assume octantoids
   // printf("INI_LMAX %d INI_CW %f INI_OW %f\n",INI_LMAX,INI_CW,INI_OW);
    //OC debugging:
    
    if(INI_LMAX>0)
       fit_oct_model_to_LC_AO(INI_LC,INI_AO,INI_OC,INI_RD);
    else
       fit_subdiv_model_to_LC_AO(INI_LC,INI_AO,INI_OC,INI_RD);
}
    
        
