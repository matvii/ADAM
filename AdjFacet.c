#include"utils.h"

void AdjFacet(int *tlist,double *vlist,int nfac,int nvert,int *A)
{
    /*A is nvertxnvert matrix, facets with common edge (i1,i2) are A(i1,i2), A(i2,i1)*/
    int i1,i2,i3;
    for(int j=0;j<nfac;j++)
    {
        i1=tlist[3*j]-1;
        i2=tlist[3*j+1]-1;
        i3=tlist[3*j+2]-1;
        set_elI(A,nvert,nvert,j+1,i1,i2);
        set_elI(A,nvert,nvert,j+1,i2,i3);
        set_elI(A,nvert,nvert,j+1,i3,i1);
    }
}
        
