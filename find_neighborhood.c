#include"utils.h"
void find_neighborhood(int *tlist,double *vlist,int nfac,int nvert,int *E,int *N,int *E2,int *A)
{
    /*E(i,j)=1 if there is and edge from vertex i to j, nvertxnvert matrix
     * N(i,j)=1, if vertex i belongs to facet j nvertxnfac matrix
     * E2(i,j)=k if edge(i,j) is in the triangle k, nvert x nvert matrix
     * A(i,j)=1 if facet i and j are adjacent, nfacxnfac matrix
     * ALL INDEXING IS FROM ZERO
     * NB: THIS WORKS ONLY FOR SHAPES WITHOUT BOUNDARY. IF BOUNDARY EXISTS, WE SHOULD TEST IF E2[i2*nvert+i1]>0 BEFORE SETTING A*/
    /*NB2: E2 is a bit problematic, maybe we should index from 1?*/
    int i1,i2,i3;
    for(int j=0;j<nfac;j++)
    {
        i1=tlist[3*j]-1;
        i2=tlist[3*j+1]-1;
        i3=tlist[3*j+2]-1;
        E[i1*nvert+i2]=1;
        E[i2*nvert+i1]=1;
        
        E[i2*nvert+i3]=1;
        E[i3*nvert+i2]=1;
        
        E[i1*nvert+i3]=1;
        E[i3*nvert+i1]=1;
        
        N[i1*nfac+j]=1;
        N[i2*nfac+j]=1;
        N[i3*nfac+j]=1;
        
        E2[i1*nvert+i2]=j;
        E2[i2*nvert+i3]=j;
        E2[i3*nvert+i1]=j;
    }
    for(int j=0;j<nfac;j++)
    {
        i1=tlist[3*j]-1;
        i2=tlist[3*j+1]-1;
        i3=tlist[3*j+2]-1;
        A[j*nfac+E2[i2*nvert+i1]]=1;
        A[j*nfac+E2[i3*nvert+i2]]=1;
        A[j*nfac+E2[i1*nvert+i3]]=1;
    }
}
/*
int main()
 {
     //int tlist[]={1,2,3,1,3,4};//2x3
     //double vlist[]={-1,0,1,0,-1,0,0,1,0,0,10,10}; //4x3
     int tlist[]={1,2,3,1,3,4,1,4,5,1,5,2,6,3,2,6,4,3,6,5,4,6,2,5}; //8x3
     double vlist[]={0,0,1,1,0,0,0,1,0,-1,0,0,0,-1,0,0,0,-1}; //6x3
     int nfac=8;
     int nvert=6;
     int *E,*N,*E2,*A;
     E=calloc(nvert*nvert,sizeof(int));
     N=calloc(nvert*nfac,sizeof(int));
     E2=calloc(nvert*nvert,sizeof(int));
     A=calloc(nfac*nfac,sizeof(int));
    find_neighborhood(tlist,vlist,nfac,nvert,E,N,E2,A);
    print_matrixI(E,nvert,nvert);
    print_matrixI(N,nvert,nfac);
    print_matrixI(E2,nvert,nvert);
    print_matrixI(A,nfac,nfac);
 }
 */
