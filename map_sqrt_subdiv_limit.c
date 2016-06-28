#include"utils.h"

void map_sqrt_subdiv_limit(int *tlist,double *vlist,int nfac,int nvert,double **D2)
{
    /*
     * Map vertices to limit vertices.
     * Here it is assumed that tlist,vlist results from sqrt(3) subdivision
     * D2 will contain the derivative matrix wrt vertices
     * NOTE: vlist is overwritten
     */
    double *vlistn=calloc(3*nvert,sizeof(double));
    int *tlistn=calloc(3*nfac,sizeof(int));
    double *D=calloc(nvert*nvert,sizeof(double));
    int *vMo=calloc(nvert*nvert,sizeof(int));
   
   int i1,i2,i3;
    for(int j=0;j<nfac;j++)
    {
        i1=tlist[3*j+0]-1;
        i2=tlist[3*j+1]-1;
        i3=tlist[3*j+2]-1;
        set_elI(vMo,nvert,nvert,1,i1,i2);
        set_elI(vMo,nvert,nvert,1,i1,i3);
        set_elI(vMo,nvert,nvert,1,i2,i1);
        set_elI(vMo,nvert,nvert,1,i2,i3);
        set_elI(vMo,nvert,nvert,1,i3,i1);
        set_elI(vMo,nvert,nvert,1,i3,i2);
    }
    int *Neigh=calloc(nvert,sizeof(int));
    
    int n=0;
    double vxsum;
    double vysum;
    double vzsum;
    double an;
    double bn;
    for(int j=0;j<nvert;j++)
    {
       
       n=0;
       vxsum=0;
       vysum=0;
       vzsum=0;
       for(int k=0;k<nvert;k++)
       {
           if(get_elI(vMo,nvert,nvert,j,k))
           {
               Neigh[n++]=k;
               vxsum+=vlist[3*k];
               vysum+=vlist[3*k+1];
               vzsum+=vlist[3*k+2];
           }
       }
       an=1.0/9.0*(4-2*cos(2*PI/n));
       bn=3.0*an/(1.0+3*an);
       vlistn[3*j+0]=(1-bn)*vlist[3*j+0]+bn*1.0/n*vxsum;
       vlistn[3*j+1]=(1-bn)*vlist[3*j+1]+bn*1.0/n*vysum;
       vlistn[3*j+2]=(1-bn)*vlist[3*j+2]+bn*1.0/n*vzsum;
       
       set_el(D,nvert,nvert,1-bn,j,j);
       for(int k=0;k<n;k++)
        set_el(D,nvert,nvert,bn/n,j,Neigh[k]);
    }
    *D2=D;
    memcpy(vlist,vlistn,sizeof(double)*3*nvert);
    
}
        
 int main()
 {
     int *tlist;
    double *vlist;
    double *D;
    int nfac,nvert;
    char file[]="mshape.txt";
    read_shape(file,&tlist,&vlist,&nfac,&nvert,0);
    map_sqrt_subdiv_limit(tlist,vlist,nfac,nvert,&D);
    write_matrix_file("/tmp/D.txt",D,nvert,nvert);
    write_matrix_file("/tmp/vlist.txt",vlist,nvert,3);
    write_matrix_fileI("/tmp/tlist.txt",tlist,nfac,3);
 }
 
