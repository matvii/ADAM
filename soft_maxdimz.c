#include"utils.h"
#include"matrix_ops.h"
void find_maxz(double *vlist,int nvert, double *maxz,int *maxv)
{
    double maxzval=-1e6;
    double maxzind=-1;
    for(int j=0;j<nvert;j++)
        if(maxzval<vlist[3*j+2])
        {
            maxzval=vlist[3*j+2];
            maxzind=j;
        }
    *maxz=maxzval;
    *maxv=maxzind;
}
 void find_minz(double *vlist,int nvert, double *minz,int *minv)
{
    double minzval=1e6;
    double minzind=-1;
    for(int j=0;j<nvert;j++)
        if(minzval>vlist[3*j+2])
        {
            minzval=vlist[3*j+2];
            minzind=j;
        }
    *minz=minzval;
    *minv=minzind;
}   
        
void soft_maxdimz(int *tlist,double *vlist,int nfacn,int nvertn,double *D,int nvert,double zmax,double c,double *zm,double *dz)
{
    /*
     * Find whether z dimension of the polyhedron is larger than the zmax
     * Calculate distance using the logistic function
     * i
     * OUTPUT:
     * zm: distance
     * dz: 1xnvert vector, only 2 entries are nonzero [if D is nvertnxnvert vector, otherwise dz is 1xnvertn vector]
     * 
     */
    int minv,maxv;
    double maxz,minz;
    double *dzn=calloc(nvertn,sizeof(double));
    find_maxz(vlist,nvertn,&maxz,&maxv);
    find_minz(vlist,nvertn,&minz,&minv);
    double zmv=(maxz-minz)-zmax;
    *zm=1.0/(1.0+exp(-c*zmv));
   
    dzn[maxv]=*zm*(1-*zm)*c;
    dzn[minv]=-*zm*(1-*zm)*c;
    if(D!=NULL)
        matrix_prod(dzn,1,nvertn,D,nvert,dz);
    else
        memcpy(dz,dzn,sizeof(double)*nvertn);
    free(dzn);
}

