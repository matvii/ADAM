#include"utils.h"
#include"matrix_ops.h"
#define PI 3.141592653589793
void butterfly_subdiv_step(int *tlist,double *vlist,int nfac,int nvert,int **tlistn2,double **vlistn2,int *nfacn2,int *nvertn2,double **D2);
void butterfly_triangle(int* tlist,double* vlist,int nfac,int nvert,unsigned char* V,short* VT,int facind,int vi1,int vi2,double * w23,int * verts,double* weights,int *nverts);
void butterfly_subdiv(int *tlist,double *vlist,int nfac,int nvert,int **tlistn,double **vlistn,int *nfacn,int *nvertn,double **D,int sdstep)
{
    int *tlistn2,*tlistn3;
    double *vlistn2,*vlistn3,*D2,*D3,*DL,*DT;
    int nvertn2,nfacn2;
    if(sdstep==2)
    {
   butterfly_subdiv_step(tlist,vlist,nfac,nvert,&tlistn2,&vlistn2,&nfacn2,&nvertn2,&D2); 
   butterfly_subdiv_step(tlistn2,vlistn2,nfacn2,nvertn2,tlistn,vlistn,nfacn,nvertn,&D3);
  
   free(vlistn2);
   free(tlistn2);
   *D=calloc(*nvertn*nvert,sizeof(double));
   matrix_prod(D3,*nvertn,nvertn2,D2,nvert,*D);
  
       free(D2);
       free(D3);
      
    }
   
    else if(sdstep==1)
    {
        
        butterfly_subdiv_step(tlist,vlist,nfac,nvert,tlistn,vlistn,nfacn,nvertn,D);
       
        
    }
    else
    {
        fprintf(stderr,"sdstep value %d not currently supported\n",sdstep);
        exit(-1);
    }
        
}
void butterfly_subdiv_step(int *tlist,double *vlist,int nfac,int nvert,int **tlistn2,double **vlistn2,int *nfacn2,int *nvertn2,double **D2)
{
/*Subdivide polyhedron given by facet list tlist and vertex list vlist using the butterfly subdivision. 
     * MEMORY IS ALLOCATED HERE */
    /* VT is a matrix where VT(i,j)=t, if edge (i,j) is in triangle t
     * V is a valance vector, V(i)=valence of vertex i
     */
    int nfacn=4*nfac;
    int nvertn=nfac+nvert-2+nvert;
    int * tlistn=calloc(3*nfacn,sizeof(int));
    double * vlistn=calloc(3*nvertn,sizeof(double));
    double *D=calloc(nvertn*nvert,sizeof(double));
    for(int j=0;j<nvert;j++)
        set_el(D,nvertn,nvert,1,j,j);
    memcpy(vlistn,vlist,3*nvert*sizeof(double));
    short *VT=calloc(nvert*nvert,sizeof(short));
    unsigned char *V=calloc(nvert,sizeof(char));
    short  *E=calloc(nvert*nvert,sizeof(short));
    int i1,i2,i3,v1,v2,v3,vind;
    int w12i,w23i,w31i;
    double w12[3],w23[3],w31[3];
    int verts[10];
    double weights[10];
    int nverts=0;
    for(int j=0;j<nfac;j++)
    {
        i1=tlist[3*j]-1;
        i2=tlist[3*j+1]-1;
        i3=tlist[3*j+2]-1;
        set_elS(VT,nvert,nvert,j,i1,i2);
        set_elS(VT,nvert,nvert,j,i2,i3);
        set_elS(VT,nvert,nvert,j,i3,i1);
        
        //NOTE: VT[i,j] goes  from 0 to nfac-1 
        V[i1]++;
        V[i2]++;
        V[i3]++;
    }
    
    vind=nvert;
    for(int j=0;j<nfac;j++)
    {
        v1=get_elI(tlist,nfac,3,j,0)-1;
        v2=get_elI(tlist,nfac,3,j,1)-1;
        v3=get_elI(tlist,nfac,3,j,2)-1;
        
        if(get_elS(E,nvert,nvert,v1,v2)==0)
        {
            butterfly_triangle(tlist,vlist,nfac,nvert,V,VT,j,0,1,w12,verts,weights,&nverts);
            for(int k=0;k<nverts;k++)
                set_el(D,nvertn,nvert,weights[k],vind,verts[k]);
            for(int k=0;k<3;k++)
                set_el(vlistn,nvertn,3,w12[k],vind,k);
            w12i=vind;
            vind++;
            set_elS(E,nvert,nvert,w12i,v1,v2);
            set_elS(E,nvert,nvert,w12i,v2,v1);
        }
        else
            w12i=get_elS(E,nvert,nvert,v1,v2);
        
        if(get_elS(E,nvert,nvert,v2,v3)==0)
        {
            butterfly_triangle(tlist,vlist,nfac,nvert,V,VT,j,1,2,w23,verts,weights,&nverts);
            for(int k=0;k<nverts;k++)
                set_el(D,nvertn,nvert,weights[k],vind,verts[k]);
            for(int k=0;k<3;k++)
                set_el(vlistn,nvertn,3,w23[k],vind,k);
            w23i=vind;
            vind++;
            set_elS(E,nvert,nvert,w23i,v2,v3);
            set_elS(E,nvert,nvert,w23i,v3,v2);
        }
        else
            w23i=get_elS(E,nvert,nvert,v2,v3);
        
        if(get_elS(E,nvert,nvert,v3,v1)==0)
        {
            butterfly_triangle(tlist,vlist,nfac,nvert,V,VT,j,2,0,w31,verts,weights,&nverts);
            for(int k=0;k<nverts;k++)
                set_el(D,nvertn,nvert,weights[k],vind,verts[k]);
            for(int k=0;k<3;k++)
                set_el(vlistn,nvertn,3,w31[k],vind,k);
            w31i=vind;
            vind++;
            set_elS(E,nvert,nvert,w31i,v3,v1);
            set_elS(E,nvert,nvert,w31i,v1,v3);
        }
        else
            w31i=get_elS(E,nvert,nvert,v3,v1);
        
        //Add 4 new triangles
        set_elI(tlistn,nfacn,3,v1+1,4*j,0);
        set_elI(tlistn,nfacn,3,w12i+1,4*j,1);
        set_elI(tlistn,nfacn,3,w31i+1,4*j,2);
        
        set_elI(tlistn,nfacn,3,v2+1,4*j+1,0);
        set_elI(tlistn,nfacn,3,w23i+1,4*j+1,1);
        set_elI(tlistn,nfacn,3,w12i+1,4*j+1,2);
        
        set_elI(tlistn,nfacn,3,v3+1,4*j+2,0);
        set_elI(tlistn,nfacn,3,w31i+1,4*j+2,1);
        set_elI(tlistn,nfacn,3,w23i+1,4*j+2,2);
        
        set_elI(tlistn,nfacn,3,w12i+1,4*j+3,0);
        set_elI(tlistn,nfacn,3,w23i+1,4*j+3,1);
        set_elI(tlistn,nfacn,3,w31i+1,4*j+3,2);
    }
    free(VT);
    free(V);
    free(E);
    *vlistn2=vlistn;
    *tlistn2=tlistn;
    *nfacn2=nfacn;
    *nvertn2=nvertn;
    *D2=D;
}
double V5(int n)
{
    return (0.25*cos(2*PI*n/5.0)+0.5*cos(4.0*PI*n/5.0))/5.0;
}
int setdiff(int *V1,int v2,int v3)
{
    /*Assume V1 1x3 int vector,
     * Return element of V1 !=v2 and v3
     */
    if(V1[0]!=v2 && V1[0]!=v3)
        return V1[0];
    if(V1[1]!=v2 && V1[1]!=v3)
        return V1[1];
    if(V1[2]!=v2 && V1[2]!=v3)
        return V1[2];
    puts("Error: Failure in setdiff in butterfly subdivision");
        exit(1);
}
        
   void butterfly_triangle(int* tlist,double* vlist,int nfac,int nvert,unsigned char* V,short* VT,int facind,int vi1,int vi2,double * w,int * verts,double* weights,int *nverts)
   {
       int v1,v2,v3;
       int t1,t2,t3,t4,t5,t6,t7,t8,t9,t10;
       int v,v4,v5,v6,v7,v8,v9,v10;
       int w0,w1,w2,w3,w4;
       double wc,a,b,c,d;
       v1=get_elI(tlist,nfac,3,facind,vi1)-1;
       v2=get_elI(tlist,nfac,3,facind,vi2)-1;
       v3=setdiff(tlist+3*facind,v1+1,v2+1)-1;
       wc=-1.0/32.0;
       a=0.5-wc;
       b=1.0/8.0+2*wc;
       c=-1.0/16.0-wc;
       d=wc;
       if(V[v1]==6 && V[v2]==6)
       {
           t1=facind;
           t2=get_elS(VT,nvert,nvert,v2,v1);
           v4=setdiff(tlist+3*t2,v1+1,v2+1)-1;
           
           t4=get_elS(VT,nvert,nvert,v1,v3);
           v5=setdiff(tlist+3*t4,v3+1,v1+1)-1;
           t3=get_elS(VT,nvert,nvert,v4,v1);
           v6=setdiff(tlist+3*t3,v1+1,v4+1)-1;
           t5=get_elS(VT,nvert,nvert,v1,v5);
           v7=setdiff(tlist+3*t5,v1+1,v5+1)-1;
           
           //Lower half of the butterfly mask, from vertex v2
           
           t8=get_elS(VT,nvert,nvert,v3,v2);
           t7=get_elS(VT,nvert,nvert,v2,v4);
           v9=setdiff(tlist+3*t7,v2+1,v4+1)-1;
           v8=setdiff(tlist+3*t8,v3+1,v2+1)-1;
           t9=get_elS(VT,nvert,nvert,v2,v9);
           v10=setdiff(tlist+3*t9,v2+1,v9+1)-1;
           
           for(int k=0;k<3;k++)
               w[k]=a*vlist[3*v1+k]+a*vlist[3*v2+k]+b*vlist[3*v4+k]+b*vlist[3*v3+k]+c*vlist[3*v5+k]+c*vlist[3*v6+k]+c*vlist[3*v9+k]+c*vlist[3*v8+k]+d*vlist[3*v7+k]+d*vlist[3*v10+k];
           *nverts=10;
           verts[0]=v1;
           verts[1]=v2;
           verts[2]=v4;
           verts[3]=v3;
           verts[4]=v5;
           verts[5]=v6;
           verts[6]=v9;
           verts[7]=v8;
           verts[8]=v7;
           verts[9]=v10;
           
           weights[0]=a;
           weights[1]=a;
           weights[2]=b;
           weights[3]=b;
           weights[4]=c;
           weights[5]=c;
           weights[6]=c;
           weights[7]=c;
           weights[8]=d;
           weights[9]=d;
       }
       else if(V[v1]==4)
       {
           t2=get_elS(VT,nvert,nvert,v1,v3);
           v4=setdiff(tlist+3*t2,v1+1,v3+1)-1;
           for(int k=0;k<3;k++)
               w[k]=3.0/4*vlist[3*v1+k]+3.0/8*vlist[3*v2+k]-1.0/8*vlist[3*v4+k];
           verts[0]=v1;
           verts[1]=v2;
           verts[2]=v4;
           *nverts=3;
           weights[0]=3.0/4;
           weights[1]=3.0/8;
           weights[2]=-1.0/8;
       }
       else if(V[v2]==4) 
       {
           t2=get_elS(VT,nvert,nvert,v2,v1);
           v4=setdiff(tlist+3*t2,v2+1,v1+1)-1;
           t3=get_elS(VT,nvert,nvert,v2,v4);
           v5=setdiff(tlist+3*t3,v2+1,v4+1)-1;
           for(int k=0;k<3;k++)
               w[k]=3.0/4*vlist[3*v2+k]+3.0/8*vlist[3*v1+k]-1.0/8*vlist[3*v5+k];
           verts[0]=v2;
           verts[1]=v1;
           verts[2]=v5;
           *nverts=3;
           weights[0]=3.0/4;
           weights[1]=3.0/8;
           weights[2]=-1.0/8;
       }
       else if(V[v1]==5)
       {
           v=v1;
           w0=v2;
           w1=v3;
           t2=get_elS(VT,nvert,nvert,v,w1);
           w2=setdiff(tlist+3*t2,w1+1,v+1)-1;
           t3=get_elS(VT,nvert,nvert,v,w2);
           w3=setdiff(tlist+3*t3,v+1,w2+1)-1;
           t4=get_elS(VT,nvert,nvert,v,w3);
           w4=setdiff(tlist+3*t4,v+1,w3+1)-1;
           for(int k=0;k<3;k++)
               w[k]=3.0/4*vlist[3*v+k]+V5(0)*vlist[3*w0+k]+V5(1)*vlist[3*w1+k]+V5(2)*vlist[3*w2+k]+V5(3)*vlist[3*w3+k]+V5(4)*vlist[3*w4+k];
           verts[0]=v;
           verts[1]=w0;
           verts[2]=w1;
           verts[3]=w2;
           verts[4]=w3;
           verts[5]=w4;
           *nverts=6;
           weights[0]=3.0/4;
           weights[1]=V5(0);
           weights[2]=V5(1);
           weights[3]=V5(2);
           weights[4]=V5(3);
           weights[5]=V5(4);
       }
       else if(V[v2]==5)
       {
           v=v2;
           w0=v1;
           w1=v3;
           t2=get_elS(VT,nvert,nvert,v,w1);
           w2=setdiff(tlist+3*t2,w1+1,v+1)-1;
           t3=get_elS(VT,nvert,nvert,v,w2);
           w3=setdiff(tlist+3*t3,v+1,w2+1)-1;
           t4=get_elS(VT,nvert,nvert,v,w3);
           w4=setdiff(tlist+3*t4,v+1,w3+1);
           for(int k=0;k<3;k++)
               w[k]=3.0/4*vlist[3*v+k]+V5(0)*vlist[3*w0+k]+V5(1)*vlist[3*w1+k]+V5(2)*vlist[3*w2+k]+V5(3)*vlist[3*w3+k]+V5(4)*vlist[3*w4+k];
           verts[0]=v;
           verts[1]=w0;
           verts[2]=w1;
           verts[3]=w2;
           verts[4]=w3;
           verts[5]=w4;
           *nverts=6;
           weights[0]=3.0/4;
           weights[1]=V5(0);
           weights[2]=V5(1);
           weights[3]=V5(2);
           weights[4]=V5(3);
           weights[5]=V5(4);
       }
   }
               
       
  
