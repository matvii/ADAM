#include"utils.h"
#include<stdlib.h>
#include<math.h>
#include"matrix_ops.h"
void area_reg(int *tlist,double *vlist,int nfac,int nvert,double *D,int dm,int dn,double *area,double *dArdv)
{
    /*Calculate area regularization, penalize divergence from mean area. area(j)=A_j/Sum(A)-1/n*/
    /*area: nfac array
     * dArdv: nfac x 3*nvert+3 array */
    /*NB ALLOCATE MEMORY BEFOREHAND*/
    double *v1,*v2,*v3,n[3];
    double dndx1[3],dndx2[3],dndx3[3];
    double dndy1[3],dndy2[3],dndy3[3];
    double dndz1[3],dndz2[3],dndz3[3];
    double *areadx,*aready,*areadz,*dArdx,*dArdy,*dArdz;
    areadx=calloc(nvert*nfac,sizeof(double));
    aready=calloc(nvert*nfac,sizeof(double));
    areadz=calloc(nvert*nfac,sizeof(double));
    double A;
    double dAdx[3],dAdy[3],dAdz[3];
    double TA=0;
    int j1,j2,j3;
    if(D!=NULL && dm!=nvert)
    {
        puts("Error: Number of vertex coordinates is not equal to the number of rows in D.");
        exit(1);
    }
    if(D==NULL)
        dn=nvert;
    for(int j=0;j<nfac;j++)
    {
        j1=tlist[3*j]-1; /*NB: WE INDEX FROM ZERO*/
        j2=tlist[3*j+1]-1;
        j3=tlist[3*j+2]-1;
        v1=vlist+3*j1;
        v2=vlist+3*j2;
        v3=vlist+3*j3;
        Calculate_Area_and_Normal_Derivative(v1,v2,v3,n,dndx1,dndx2,dndx3,dndy1,dndy2,dndy3,dndz1,dndz2,dndz3,&A,dAdx,dAdy,dAdz);
        TA=TA+A;
        area[j]=A;
        areadx[j*nvert+j1]=dAdx[0];
        areadx[j*nvert+j2]=dAdx[1];
        areadx[j*nvert+j3]=dAdx[2];
        
        aready[j*nvert+j1]=dAdy[0];
        aready[j*nvert+j2]=dAdy[1];
        aready[j*nvert+j3]=dAdy[2];
        
        areadz[j*nvert+j1]=dAdz[0];
        areadz[j*nvert+j2]=dAdz[1];
        areadz[j*nvert+j3]=dAdz[2];
    }
    //Divide with total area and calculate final derivatives
    //Derivative of total area wrt vertices
    double *dTAdx,*dTAdy,*dTAdz;
    dTAdx=malloc(nvert*sizeof(double));
    dTAdy=malloc(nvert*sizeof(double));
    dTAdz=malloc(nvert*sizeof(double));
    for(int k=0;k<nvert;k++)
    {
        dTAdx[k]=0;
        dTAdy[k]=0;
        dTAdz[k]=0;
        for(int j=0;j<nfac;j++)
        {
          dTAdx[k]+=areadx[j*nvert+k];
          dTAdy[k]+=aready[j*nvert+k];  
          dTAdz[k]+=areadz[j*nvert+k];
        }
    }
    if(D==NULL)
    {
        
    for(int j=0;j<nfac;j++)
    {
        
        for(int k=0;k<nvert;k++)
        {
            dArdv[j*(3*nvert+3)+k]=(areadx[j*nvert+k]*TA-area[j]*dTAdx[k])/pow(TA,2);
            dArdv[j*(3*nvert+3)+k+nvert]=(aready[j*nvert+k]*TA-area[j]*dTAdy[k])/pow(TA,2);
            dArdv[j*(3*nvert+3)+k+2*nvert]=(areadz[j*nvert+k]*TA-area[j]*dTAdz[k])/pow(TA,2);
        }
        area[j]=area[j]/TA-1.0/nfac;
    }
    
    }
    else
    {
        dArdx=calloc(nvert*nfac,sizeof(double));
        dArdy=calloc(nvert*nfac,sizeof(double));
        dArdz=calloc(nvert*nfac,sizeof(double));
        for(int j=0;j<nfac;j++)
    {
        
        for(int k=0;k<nvert;k++)
        {
            areadx[j*nvert+k]=(areadx[j*nvert+k]*TA-area[j]*dTAdx[k])/pow(TA,2);
            aready[j*nvert+k]=(aready[j*nvert+k]*TA-area[j]*dTAdy[k])/pow(TA,2);
            areadz[j*nvert+k]=(areadz[j*nvert+k]*TA-area[j]*dTAdz[k])/pow(TA,2);
        }
        area[j]=area[j]/TA-1.0/nfac;
    }
    matrix_prod(areadx,nfac,nvert,D,dn,dArdx);
    matrix_prod(aready,nfac,nvert,D,dn,dArdy);
    matrix_prod(areadz,nfac,nvert,D,dn,dArdz);
    
    set_submatrix(dArdv,nfac,3*dn+3,dArdx,nfac,dn,0,0);
    set_submatrix(dArdv,nfac,3*dn+3,dArdy,nfac,dn,0,dn);
    set_submatrix(dArdv,nfac,3*dn+3,dArdz,nfac,dn,0,2*dn);
     free(dArdx);
    free(dArdy);
    free(dArdz);
    }
   
    free(areadx);
    free(aready);
    free(areadz);
    free(dTAdx);
    free(dTAdy);
    free(dTAdz);
}
/* 
 int main()
 {
      int nfac=8,nfacn;
      int nvert=6,nvertn;
     int *tlist;
     double *vlist;
     double *vlistn;
     int *tlistn;
     double *D;
     char shapefile[]="shape2.txt";
     read_shape(shapefile,&tlist,&vlist,&nfac,&nvert,0);
      double ini_dims[]={90,90,70};
      mul_cols(vlist,nvert,3,ini_dims);
     Sqrt3_Subdiv(tlist,vlist,nfac,nvert,&tlistn,&vlistn,&nfacn,&nvertn,&D);
     int nedgen=nfacn+nvertn-2;
     int nedge=nfac+nvert-2;
     double res,*dresdv;
   
   
     double *area,*areadv;
     area=calloc(nfacn,sizeof(double));
     areadv=calloc(nfacn*(3*nvert+3),sizeof(double));
     
     area_reg(tlistn,vlistn,nfacn,nvertn,D,nvertn,nvert,area,areadv);
     //print_matrix(area,1,nfac);
    // print_matrix(areadv,nfac,3*nvert+3);
     
    // free(D);
     write_matrix_file("/tmp/area2.txt",area,1,nfacn);
    write_matrix_file("/tmp/areadv2.txt",areadv,nfacn,3*nvert+3);
    free(area);
     free(areadv);
 }
 
*/
