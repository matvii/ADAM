#include"utils.h"
/*#include<stdlib.h>
#include<math.h>
#include<stdio.h>
#include<cblas.h>*/
#include"matrix_ops.h"
#define PI 3.141592653589793
void sqrt3_subdiv(int *tlist,double *vlist,int nfac,int nvert,int **tlistn2,double **vlistn2,int *nfacn2,int *nvertn2,double **D2)
{
    /*Subdivide polyhedron given by facet list tlist and vertex list vlist using the Sqrt3 subdivision. 
     * MEMORY IS ALLOCATED HERE */
    /*Number of vertices in subdivided triangle:
     * We know that nvert=2+1/2*nfac*/
    int nfacn=3*nfac;
    int nvertn=2+nfacn/2;
    int * tlistn=calloc(3*nfacn,sizeof(int));
    double * vlistn=calloc(3*nvertn,sizeof(double));
    double *D=calloc(nvertn*nvert,sizeof(double));
   
    //print_matrixI(tlistn,nfacn,3);
    /*Generate adjacency matrix: Facets with common edge (i1,i2) are A(i1,i2) and A(i2,i1);*/
    /*Note that facet numbering in A starts from 1*/
    int *A=calloc(nvert*nvert,sizeof(int));
    int *vMo=calloc(nvert*nvert,sizeof(int));
    int i0,i1,i2;
    double v0x,v0y,v0z;
    double v1x,v1y,v1z;
    double v2x,v2y,v2z;
    double cx,cy,cz;
    int *centindex=calloc(nfac,sizeof(int));
    for(int j=0;j<nfac;j++)
    {
        i0=get_elI(tlist,nfac,3,j,0)-1;
        i1=get_elI(tlist,nfac,3,j,1)-1;
        i2=get_elI(tlist,nfac,3,j,2)-1;
        set_elI(vMo,nvert,nvert,1,i0,i1);
        set_elI(vMo,nvert,nvert,1,i0,i2);
        set_elI(vMo,nvert,nvert,1,i1,i0);
        set_elI(vMo,nvert,nvert,1,i1,i2);
        set_elI(vMo,nvert,nvert,1,i2,i0);
        set_elI(vMo,nvert,nvert,1,i2,i1);
        
        set_elI(A,nvert,nvert,j+1,i0,i1);
        set_elI(A,nvert,nvert,j+1,i1,i2);
        set_elI(A,nvert,nvert,j+1,i2,i0);
    }
    /* Update old vertices*/
    int nbv=0;
    double bn1=0,bn2x,bn2y,bn2z;
    int *V;
    int vindex=0;
    int count=0;
    for(int j=0;j<nvert;j++)
    {
        nbv=ind2vec(vMo,nvert,&V,j);
        bn1=1.0/(9.0*nbv)*(4.0-2.0*cos(2.0*PI/nbv));
        bn2x=(1.0-nbv*bn1)*get_el(vlist,nvert,3,j,0)+bn1*sum_matelC(vlist,nvert,3,V,nbv,0);
        bn2y=(1.0-nbv*bn1)*get_el(vlist,nvert,3,j,1)+bn1*sum_matelC(vlist,nvert,3,V,nbv,1);
        bn2z=(1.0-nbv*bn1)*get_el(vlist,nvert,3,j,2)+bn1*sum_matelC(vlist,nvert,3,V,nbv,2);
       
        set_el(vlistn,nvertn,3,bn2x,j,0);
        set_el(vlistn,nvertn,3,bn2y,j,1);
        set_el(vlistn,nvertn,3,bn2z,j,2);
        set_el(D,nvertn,nvert,1.0-nbv*bn1,j,j);
        for(int k=0;k<nbv;k++)
            set_el(D,nvertn,nvert,bn1,j,V[k]);
        free(V);
    }
    /* Add new vertices and triangles, minus one triangle*/
    
    for(int j=0;j<nfac;j++)
    {
        i0=get_elI(tlist,nfac,3,j,0)-1;
        i1=get_elI(tlist,nfac,3,j,1)-1;
        i2=get_elI(tlist,nfac,3,j,2)-1;
        v0x=get_el(vlist,nvert,3,i0,0);
        v0y=get_el(vlist,nvert,3,i0,1);
        v0z=get_el(vlist,nvert,3,i0,2);
        v1x=get_el(vlist,nvert,3,i1,0);
        v1y=get_el(vlist,nvert,3,i1,1);
        v1z=get_el(vlist,nvert,3,i1,2);
        v2x=get_el(vlist,nvert,3,i2,0);
        v2y=get_el(vlist,nvert,3,i2,1);
        v2z=get_el(vlist,nvert,3,i2,2);
        
        cx=(v0x+v1x+v2x)/3.0;
        cy=(v0y+v1y+v2y)/3.0;
        cz=(v0z+v1z+v2z)/3.0;
        vindex=count+nvert;
        count++;
        set_el(vlistn,nvertn,3,cx,vindex,0);
        set_el(vlistn,nvertn,3,cy,vindex,1);
        set_el(vlistn,nvertn,3,cz,vindex,2);
        set_elI(centindex,1,nfac,vindex,0,j);
        set_el(D,nvertn,nvert,1.0/3.0, vindex,i0);
        set_el(D,nvertn,nvert,1.0/3.0, vindex,i1);
        set_el(D,nvertn,nvert,1.0/3.0, vindex,i2);
    }
    /*Now we flip the edge between facets*/
    int tcount=0;
    int nf1,nf2,nf3;
    for(int j=0;j<nfac;j++)
    {
       i0=get_elI(tlist,nfac,3,j,0)-1;
       i1=get_elI(tlist,nfac,3,j,1)-1;
       i2=get_elI(tlist,nfac,3,j,2)-1; 
       /*Take the neighbor facets*/
       nf1=get_elI(A,nvert,nvert,i1,i0);
       nf2=get_elI(A,nvert,nvert,i2,i1);
       nf3=get_elI(A,nvert,nvert,i0,i2);
       if(nf1>0)
       {
           /*tlistn(tcount,:)=[centindex(j) i1 centindex(nf1)];*/
        //   printf("nf1: %d %d %d %d\n",nf1,centindex[j]+1,i0+1,centindex[nf1-1]+1);
           set_elI(tlistn,nfacn,3,centindex[j]+1,tcount,0); /*centindex+1 or not? CHECK THIS!*/
           set_elI(tlistn,nfacn,3,i0+1,tcount,1);
           set_elI(tlistn,nfacn,3,centindex[nf1-1]+1,tcount,2);
           tcount++;
           set_elI(tlistn,nfacn,3,centindex[nf1-1]+1,tcount,0);
           set_elI(tlistn,nfacn,3,i1+1,tcount,1);
           set_elI(tlistn,nfacn,3,centindex[j]+1,tcount,2);
           tcount++;
           set_elI(A,nvert,nvert,0,i1,i0);
           set_elI(A,nvert,nvert,0,i0,i1);
       }
       if(nf2>0) 
       {
           set_elI(tlistn,nfacn,3,centindex[j]+1,tcount,0); /*centindex+1 or not? CHECK THIS!*/
           set_elI(tlistn,nfacn,3,i1+1,tcount,1);
           set_elI(tlistn,nfacn,3,centindex[nf2-1]+1,tcount,2);
           tcount++;
           set_elI(tlistn,nfacn,3,centindex[nf2-1]+1,tcount,0); /*centindex+1 or not? CHECK THIS!*/
           set_elI(tlistn,nfacn,3,i2+1,tcount,1);
           set_elI(tlistn,nfacn,3,centindex[j]+1,tcount,2);
           tcount++;
           set_elI(A,nvert,nvert,0,i2,i1);
           set_elI(A,nvert,nvert,0,i1,i2);
       }
       if(nf3>0)
       {
            set_elI(tlistn,nfacn,3,centindex[j]+1,tcount,0); /*centindex+1 or not? CHECK THIS!*/
           set_elI(tlistn,nfacn,3,i2+1,tcount,1);
           set_elI(tlistn,nfacn,3,centindex[nf3-1]+1,tcount,2);
           tcount++;
            set_elI(tlistn,nfacn,3,centindex[nf3-1]+1,tcount,0); /*centindex+1 or not? CHECK THIS!*/
           set_elI(tlistn,nfacn,3,i0+1,tcount,1);
           set_elI(tlistn,nfacn,3,centindex[j]+1,tcount,2);
           tcount++;
           set_elI(A,nvert,nvert,0,i2,i0);
           set_elI(A,nvert,nvert,0,i0,i2);
       }
    }
    *vlistn2=vlistn;
    *tlistn2=tlistn;
    *nfacn2=nfacn;
    *nvertn2=nvertn;
    *D2=D;
    free(A);
    free(vMo);
    free(centindex);
   
}
void map_sqrt_subdiv_limit(int *tlist,double *vlist,int nfac,int nvert,double **D2)
{
    /*
     * Map vertices to limit vertices.
     * Here it is assumed that tlist,vlist results from sqrt(3) subdivision
     * D2 will contain the derivative matrix wrt vertices
     * NOTE: vlist is overwritten
     */
    double *vlistn=calloc(3*nvert,sizeof(double));
    //int *tlistn=calloc(3*nfac,sizeof(int));
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
    free(vlistn);
    free(vMo);
    free(Neigh);
    
}
void Sqrt3_Subdiv(int *tlist,double* vlist,int nfac,int nvert,int **tlistn,double **vlistn,int *nfacn,int *nvertn,double **D,int sdstep)
/*Do subdivision*/
{
    int *tlistn2,*tlistn3;
    double *vlistn2,*vlistn3,*D2,*D3,*DL,*DT;
    int nvertn2,nfacn2;
    if(sdstep==2)
    {
   sqrt3_subdiv(tlist,vlist,nfac,nvert,&tlistn2,&vlistn2,&nfacn2,&nvertn2,&D2); 
   sqrt3_subdiv(tlistn2,vlistn2,nfacn2,nvertn2,tlistn,vlistn,nfacn,nvertn,&D3);
   map_sqrt_subdiv_limit(*tlistn,*vlistn,*nfacn,*nvertn,&DL);
   free(vlistn2);
   free(tlistn2);
   DT=calloc(*nvertn*nvert,sizeof(double));
   *D=calloc(*nvertn*nvert,sizeof(double));
   matrix_prod(D3,*nvertn,nvertn2,D2,nvert,DT);
   matrix_prod(DL,*nvertn,*nvertn,DT,nvert,*D);
       free(D2);
       free(D3);
       free(DL);
       free(DT);
    }
    else if(sdstep==1)
    {
        
        sqrt3_subdiv(tlist,vlist,nfac,nvert,tlistn,vlistn,nfacn,nvertn,&DT);
        map_sqrt_subdiv_limit(*tlistn,*vlistn,*nfacn,*nvertn,&DL);
        *D=calloc(*nvertn*nvert,sizeof(double));
        matrix_prod(DL,*nvertn,*nvertn,DT,nvert,*D);
        free(DT);
        free(DL);
        
    }
    else if(sdstep==0)
    {
        *tlistn=calloc(nfac*3,sizeof(int));
        *vlistn=calloc(nvert*3,sizeof(double));
        
        memcpy(*tlistn,tlist,sizeof(int)*3*nfac);
        memcpy(*vlistn,vlist,sizeof(double)*3*nvert);
         *nfacn=nfac;
        *nvertn=nvert;
        *D=calloc(*nvertn*nvert,sizeof(double));
        map_sqrt_subdiv_limit(*tlistn,*vlistn,*nfacn,*nvertn,D);
       
    }
    else if(sdstep==-1)
    {
        *tlistn=calloc(nfac*3,sizeof(int));
        *vlistn=calloc(nvert*3,sizeof(double));
        *D=NULL;
        memcpy(*tlistn,tlist,sizeof(int)*3*nfac);
        memcpy(*vlistn,vlist,sizeof(double)*3*nvert);
        *nfacn=nfac;
        *nvertn=nvert;
    }
    else
    {
        fprintf(stderr,"sdstep value %d not currently supported\n",sdstep);
        exit(-1);
    }
        
    
        
}
/*
  int main()
  {
     char  filename[]="mshape.txt";
      int *tlist,*tlistn,*tlistn2;
      double *vlist,*vlistn,*vlistn2,*D2,*D;
      int nvert,nfac,nvertn,nfacn,nvertn2,nfacn2;
      read_shape(filename,&tlist,&vlist,&nfac,&nvert,0);
      printf("nfac: %d nvert: %d\n",nfac,nvert);
     // sqrt3_subdiv(tlist,vlist,nfac,nvert,&tlistn,&vlistn,&nfacn,&nvertn,&D); 
     // printf("nfac: %d nvert: %d\n",nfacn,nvertn);
     // sqrt3_subdiv(tlistn,vlistn,nfacn,nvertn,&tlistn2,&vlistn2,&nfacn2,&nvertn2,&D2);
      //printf("nfac: %d nvert: %d\n",nfacn2,nvertn2);
      //print_matrixI(tlist,nfac,3);
      //print_matrix(vlist,nvert,3);
      //sqrt3_subdiv(tlist,nfac,vlist,nvert,&tlistn,&nfacn,&vlistn,&nvertn,&D);
      Sqrt3_Subdiv(tlist,vlist,nfac,nvert,&tlistn,&vlistn,&nfacn,&nvertn,&D,1);
     printf("nfacn: %d nvertn: %d\n",nfacn,nvertn);
      //print_matrix(D,nvertn,nvert);
    //print_matrixI(tlistn,10,3);
    //print_matrixI(tlist,8,3);
     write_shape_file("/tmp/sshape.txt",tlistn,vlistn,nfacn,nvertn);
     write_matrix_file("/tmp/D.txt",D,nvertn,nvert);
     
     free(tlistn);
     free(vlistn);
     
     free(tlist);
     free(vlist);
     */
   /*
   int tlist[]={1,2,     3,
     1,     3,     4,
     1,     4,     5,
     1,     5,     2,
     6,     3,     2,
     6,     4,     3,
     6,     5,     4,
     6,     2,     5};
   double vlist[]={  0,        0,    1.0000,
    1.0000,         0,    0.0000,
    0.0000,    1.0000,    0.0000,
   -1.0000,    0.0000,    0.0000,
   -0.0000,   -1.0000,    0.0000,
    0.0000,         0,   -1.0000};
   double *vlistn,*D;
   int *tlistn,nvertn,nvert,nfac,nfacn;
   nvert=6;
   nfac=8;
   sqrt3_subdiv(tlist,vlist,nfac,nvert,&tlistn,&vlistn,&nfacn,&nvertn,&D);
   printf("nfacn: %d nvertn: %d\n",nfacn,nvertn);
    write_shape_file("/tmp/sshape.txt",tlistn,vlistn,nfacn,nvertn);
     write_matrix_file("/tmp/D.txt",D,nvertn,nvert);
  }
       
  */     

