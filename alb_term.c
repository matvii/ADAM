#include"utils.h"
double Albedo_Term(int *tlist,double *vlist,int nfac,int nvert, double *Alimit,double *Alb,int index,double* dAlbv)
/*Calculate albedo values for the each facet
 * Here albedo value is mean of vertex values
 */
{
    int v1,v2,v3;
    int talb;
    double la,ha;
    double Albedo;
    la=Alimit[0];
    ha=Alimit[1];
    v1=tlist[3*index]-1;
    v2=tlist[3*index+1]-1;
    v3=tlist[3*index+2]-1;
    Albedo=(la+ha)/2+(ha-la)/2*(tanh(Alb[v1])+tanh(Alb[v2])+tanh(Alb[v3]))/3;
    if(dAlbv!=NULL)
    {
        dAlbv[0]=(ha-la)/2*(1-pow(tanh(Alb[v1]),2))/3;
        dAlbv[1]=(ha-la)/2*(1-pow(tanh(Alb[v2]),2))/3;
        dAlbv[2]=(ha-la)/2*(1-pow(tanh(Alb[v3]),2))/3;
    }
    return Albedo;
}

            
