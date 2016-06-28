#include"utils.h"
#include"matrix_ops.h"
void convex_reg(int* tlist,double* vlist,int nfac,int nvert,double *D,int dm,int dn,double *res,double *drdv)
{
    /*Calculate convex reg term
     * OUTPUT:
     * res a double
     * drdv 3*nvert+3 array (dn array if D!=NULL)
     */
    
    int nedge;
    double result=0;
    int *E,*N,*E2,*A;
    double *v1,*v2,*v3,n1[3],n2[3];
    double *w1,*w2,*w3;
    double cangle;
    double TA=0;
    double *dTAdx,*dTAdy,*dTAdz;
    dTAdx=calloc(nvert,sizeof(double));
    dTAdy=calloc(nvert,sizeof(double));
    dTAdz=calloc(nvert,sizeof(double));
    
    double dn1dx1[3],dn1dx2[3],dn1dx3[3];
    double dn1dy1[3],dn1dy2[3],dn1dy3[3];
    double dn1dz1[3],dn1dz2[3],dn1dz3[3];
    
    double dn2dx1[3],dn2dx2[3],dn2dx3[3];
    double dn2dy1[3],dn2dy2[3],dn2dy3[3];
    double dn2dz1[3],dn2dz2[3],dn2dz3[3];
    double area1,area2;
    double dA1dx[3],dA1dy[3],dA1dz[3];
    double dA2dx[3],dA2dy[3],dA2dz[3];
    double *dresdx,*dresdy,*dresdz;
    dresdx=calloc(nvert,sizeof(double));
    dresdy=calloc(nvert,sizeof(double));
    dresdz=calloc(nvert,sizeof(double));
    E=calloc(nvert*nvert,sizeof(int));
    N=calloc(nvert*nfac,sizeof(int));
    E2=calloc(nvert*nvert,sizeof(int));
    A=calloc(nfac*nfac,sizeof(int));
    
    if(D!=NULL && dm!=nvert)
    {
        puts("Error: Number of vertex coordinates is not equal to the number of rows in D.");
        exit(1);
    }
    if(D==NULL)
        dn=nvert;
    double *drdx,*drdy,*drdz;
    drdx=calloc(dn,sizeof(double));
    drdy=calloc(dn,sizeof(double));
    drdz=calloc(dn,sizeof(double));
    
    find_neighborhood(tlist,vlist,nfac,nvert,E,N,E2,A);
    
    free(E);
    free(E2);
    free(N);
    nedge=nfac+nvert-2;
    int ind;
    int CB=0;
    int i1,i2,i3;
    int j1,j2,j3;
    int count=0;
    int *NumofBlocks;
    int *IndexofBlocks;
    double *normal,*centroid;
    NumofBlocks=calloc(nfac,sizeof(int));
    IndexofBlocks=calloc(nfac*nfac,sizeof(int));
    normal=malloc(3*nfac*sizeof(double)); //We don't really need these
    centroid=malloc(3*nfac*sizeof(double));
    FacetsOverHorizon(tlist,vlist,nfac,nvert,normal,centroid,NumofBlocks,IndexofBlocks);
    free(normal);
    free(centroid);
    
    
    for(int k=0;k<nfac;k++)
    {
        
        i1=tlist[3*k]-1;
        i2=tlist[3*k+1]-1;
        i3=tlist[3*k+2]-1;
        v1=vlist+3*i1;
        v2=vlist+3*i2;
        v3=vlist+3*i3;
        
        Calculate_Area_and_Normal_Derivative(v1,v2,v3,n1,dn1dx1,dn1dx2,dn1dx3,dn1dy1,dn1dy2,dn1dy3,dn1dz1,dn1dz2,dn1dz3,&area1,dA1dx,dA1dy,dA1dz);
        
        TA+=area1;
        dTAdx[i1]+=dA1dx[0];
        dTAdx[i2]+=dA1dx[1];
        dTAdx[i3]+=dA1dx[2];
        
        dTAdy[i1]+=dA1dy[0];
        dTAdy[i2]+=dA1dy[1];
        dTAdy[i3]+=dA1dy[2];
        
        dTAdz[i1]+=dA1dz[0];
        dTAdz[i2]+=dA1dz[1];
        dTAdz[i3]+=dA1dz[2];
        if(NumofBlocks[k]==0)
            continue; //There are no facets above the local horizon of current facets
            
            for(int j=0;j<NumofBlocks[k];j++)
            {
                CB=IndexofBlocks[nfac*k+j]-1; //Current potential blocker
                
                if(A[nfac*k+CB]==0)
                    continue; //Facets are not adjacent
                    
                    j1=tlist[3*CB]-1;
                j2=tlist[3*CB+1]-1;
                j3=tlist[3*CB+2]-1;
                
                w1=vlist+3*j1;
                w2=vlist+3*j2;
                w3=vlist+3*j3;
                Calculate_Area_and_Normal_Derivative(w1,w2,w3,n2,dn2dx1,dn2dx2,dn2dx3,dn2dy1,dn2dy2,dn2dy3,dn2dz1,dn2dz2,dn2dz3,&area2,dA2dx,dA2dy,dA2dz);
                
                cangle=DOT(n1,n2);
                result+=(area2)*(1-cangle);
                dresdx[i1]+=(area2)*(-DOT(dn1dx1,n2));
                dresdx[i2]+=(area2)*(-DOT(dn1dx2,n2));
                dresdx[i3]+=(area2)*(-DOT(dn1dx3,n2));
                
                dresdy[i1]+=(area2)*(-DOT(dn1dy1,n2));
                dresdy[i2]+=(area2)*(-DOT(dn1dy2,n2));
                dresdy[i3]+=(area2)*(-DOT(dn1dy3,n2));
                
                dresdz[i1]+=(area2)*(-DOT(dn1dz1,n2));
                dresdz[i2]+=(area2)*(-DOT(dn1dz2,n2));
                dresdz[i3]+=(area2)*(-DOT(dn1dz3,n2));
                
                dresdx[j1]+=dA2dx[0]*(1-cangle)+(area2)*(-DOT(dn2dx1,n1));
                dresdx[j2]+=dA2dx[1]*(1-cangle)+(area2)*(-DOT(dn2dx2,n1));
                dresdx[j3]+=dA2dx[2]*(1-cangle)+(area2)*(-DOT(dn2dx3,n1));
                
                dresdy[j1]+=dA2dy[0]*(1-cangle)+(area2)*(-DOT(dn2dy1,n1));
                dresdy[j2]+=dA2dy[1]*(1-cangle)+(area2)*(-DOT(dn2dy2,n1));
                dresdy[j3]+=dA2dy[2]*(1-cangle)+(area2)*(-DOT(dn2dy3,n1));
                
                dresdz[j1]+=dA2dz[0]*(1-cangle)+(area2)*(-DOT(dn2dz1,n1));
                dresdz[j2]+=dA2dz[1]*(1-cangle)+(area2)*(-DOT(dn2dz2,n1));
                dresdz[j3]+=dA2dz[2]*(1-cangle)+(area2)*(-DOT(dn2dz3,n1));
                
                
            }
            
            
    }
    
    free(A);
    (*res)=result/TA;
    if(D==NULL)
    {
        for(int j=0;j<nvert;j++)
        {
           
            drdx[j]=(dresdx[j]*TA-result*dTAdx[j])/pow(TA,2);
            drdy[j]=(dresdy[j]*TA-result*dTAdy[j])/pow(TA,2);
            drdz[j]=(dresdz[j]*TA-result*dTAdz[j])/pow(TA,2);
        }
        set_submatrix(drdv,1,3*dn+3,drdx,1,dn,0,0);
        set_submatrix(drdv,1,3*dn+3,drdy,1,dn,0,dn);
        set_submatrix(drdv,1,3*dn+3,drdz,1,dn,0,2*dn);
    }
    else
    {
        for(int j=0;j<nvert;j++)
        {
           
            dresdx[j]=(dresdx[j]*TA-result*dTAdx[j])/pow(TA,2);
            dresdy[j]=(dresdy[j]*TA-result*dTAdy[j])/pow(TA,2);
            dresdz[j]=(dresdz[j]*TA-result*dTAdz[j])/pow(TA,2);
        }
        matrix_prod(dresdx,1,nvert,D,dn,drdx);
        matrix_prod(dresdy,1,nvert,D,dn,drdy);
        matrix_prod(dresdz,1,nvert,D,dn,drdz);
        set_submatrix(drdv,1,3*dn+3,drdx,1,dn,0,0);
        set_submatrix(drdv,1,3*dn+3,drdy,1,dn,0,dn);
        set_submatrix(drdv,1,3*dn+3,drdz,1,dn,0,2*dn);
    } 
    free(drdx);
    free(drdy);
    free(drdz);
    free(dTAdx);
    free(dTAdy);
    free(dTAdz);
    free(dresdx);
    free(dresdy);
    free(dresdz);
    free(NumofBlocks);
    free(IndexofBlocks);
    
    
}
/*
int main()
 {
     //int tlist[]={1,2,3,1,3,4};//2x3
     //double vlist[]={-1,0,1,0,-1,0,0,1,0,0,10,10}; //4x3
     int tlist[]={1,2,3,1,3,4,1,4,5,1,5,2,6,3,2,6,4,3,6,5,4,6,2,5}; //8x3
     double vlist[]={10,0,1,1,0,0,0,1,0,-1,0,0,0,-1,0,0,0,-1}; //6x3
     int nfac=8;
     int nvert=6;
     int nedge=nfac+nvert-2;
     double res,*dresdv;
     int *EV;
     double *D;
     D=calloc(6*6,sizeof(double));
     for(int j=0;j<6;j++)
         D[j+6*j]=2;
     
     print_matrix(D,6,6);
     dresdv=calloc(3*nvert+3,sizeof(double));
   
    
     convex_reg(tlist,vlist,nfac,nvert,D,6,6,&res,dresdv);
    printf("res: %f\n",res);
     print_matrix(dresdv,1,3*nvert+3);
     free(D);
    free(dresdv);
   

 }
 
*/
