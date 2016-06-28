#include"utils.h"
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include"matrix_ops.h"

void octantoid_to_trimesh(double *a,int LMAX,int nrows,int *tlist,double *vlist,double *dvda,int padding)
{
    /*
     * Convert octantoid representation a=[ax by cz], where ax by ca are 1 x (LMAX+1)^2 vectors,
     * to the standard tlist,vlist representation.
     * also outputs the derivative matrix dx/da
     * vlist is 4*nrows^2+2 x 3 matrix
     * tlist 8nrows^2 x 3 matrix
     * Pad the derivative matrix with additional "padding" zero columns and rows, with identity matrix
     * [da 0
     *  0 I]
     */
    double *theta,*phi;
    int nvert=4*pow(nrows,2)+2;
    theta=calloc(nvert,sizeof(double));
    phi=calloc(nvert,sizeof(double));
    triangulate_sphere(nrows,theta,phi,tlist);
    int al=pow(LMAX+1,2);
    double *B=calloc(nvert*al,sizeof(double));
    double *Bx=calloc(nvert,sizeof(double));
    double *Bxy=calloc(nvert,sizeof(double));
    double *Bxz=calloc(nvert,sizeof(double));
    double *shxy=calloc(al,sizeof(double));
    double *shxz=calloc(al,sizeof(double));
    for(int j=0;j<al;j++)
    {
        shxy[j]=a[j]+a[j+al];
        shxz[j]=a[j]+a[j+2*al];
    }
    for(int j=0;j<=LMAX;j++)
        for(int k=-j;k<=j;k++)
            for(int i=0;i<nvert;i++)
                B[i*al+j*(j+1)+k+1-1]=SH(j,k,theta[i],phi[i]);
            matrix_vectorprod(B,nvert,al,a,Bx,0);
        matrix_vectorprod(B,nvert,al,shxy,Bxy,0);
    matrix_vectorprod(B,nvert,al,shxz,Bxz,0);
   
    for(int j=0;j<nvert;j++)
    {
        vlist[3*j]=exp(Bx[j])*sin(theta[j])*cos(phi[j]);
        vlist[3*j+1]=exp(Bxy[j])*sin(theta[j])*sin(phi[j]);
        vlist[3*j+2]=exp(Bxz[j])*cos(theta[j]);
    }
    zero_array(dvda,3*nvert*3*al);
    double *dxda=calloc(nvert*al,sizeof(double));
    double *dyda=calloc(nvert*al,sizeof(double));
    double *dzda=calloc(nvert*al,sizeof(double));
    for(int i=0;i<nvert;i++)
        for(int j=0;j<al;j++)
        {
            dxda[i*al+j]=vlist[3*i]*B[i*al+j];
            dyda[i*al+j]=vlist[3*i+1]*B[i*al+j];
            dzda[i*al+j]=vlist[3*i+2]*B[i*al+j];
        }
        if(padding==0)
        {
    set_submatrix(dvda,3*nvert,3*al,dxda,nvert,al,0,0);
    set_submatrix(dvda,3*nvert,3*al,dyda,nvert,al,nvert,0);
    set_submatrix(dvda,3*nvert,3*al,dyda,nvert,al,nvert,al);
    set_submatrix(dvda,3*nvert,3*al,dzda,nvert,al,2*nvert,0);
    set_submatrix(dvda,3*nvert,3*al,dzda,nvert,al,2*nvert,2*al);
        }
        else
        {
    set_submatrix(dvda,3*nvert+padding,3*al+padding,dxda,nvert,al,0,0);
    set_submatrix(dvda,3*nvert+padding,3*al+padding,dyda,nvert,al,nvert,0);
    set_submatrix(dvda,3*nvert+padding,3*al+padding,dyda,nvert,al,nvert,al);
    set_submatrix(dvda,3*nvert+padding,3*al+padding,dzda,nvert,al,2*nvert,0);
    set_submatrix(dvda,3*nvert+padding,3*al+padding,dzda,nvert,al,2*nvert,2*al);
    for(int j=0;j<padding;j++)
        set_el(dvda,3*nvert+padding,3*al+padding,1,3*nvert+j,3*al+j);
        }
    free(B);
    free(Bx);
    free(theta);
    free(phi);
    free(Bxy);
    free(Bxz);
    free(shxy);
    free(shxz);
    free(dxda);
    free(dyda);
    free(dzda);
}
/*
int main()
{
    int LMAX=1;
    int al=pow(LMAX+1,2);
    int *tlist;
    int nrows=1;
    double *vlist;
    int nfac=8*pow(nrows,2);
    int nvert=4*pow(nrows,2)+2;
    tlist=calloc(3*nfac,sizeof(int));
    vlist=calloc(4*nvert,sizeof(double));
    double *a=calloc(3*al,sizeof(double));
    
    //double a[]={15.635,0.013691,0.22937,0.1341,0.10937,0.0024299,0.10371,-0.024498,0.03098,0.063325,0.010701,0.022324,-0.097043,0.080485,0.09755,0.12292,0.068109,-0.08211,0.014059,-0.0087267,-0.038375,-0.037296,0.0047849,-0.0040609,-0.012102,0.026923,0.0071193,0.014678,0.001805,0.0022437,-0.011298,0.048475,0.03679,-0.049613,0.033315,-0.00017863,-0.29452,0.013917,-0.073025,0.012105,0.10192,0.031868,0.023307,-0.011363,0.036648,0.020415,-0.014803,-0.015465,0.0091759,0.011069,0.030681,0.01112,0.0077766,0.019238,0.013775,0.0086287,-0.011219,-0.046121,0.008848,0.012846,0.026349,-0.01996,0.023934,-0.00027727,-0.010195,-0.030403,0.0040127,-0.022541,0.0039028,-0.0049176,-0.0065048,-0.004749,-2.2149,0.0034772,-0.00035969,0.044307,-0.026401,-0.026718,-0.0010601,-0.058143,0.013648,-0.0084772,0.025036,0.011575,0.0012623,0.054822,0.011374,-0.0087401,0.0055519,0.0044708,-0.066513,-0.035699,0.020633,0.043237,-0.023028,0.060712,-0.0033598,0.0057378,-0.014311,0.027262,0.053754,0.0052758,0.0072281,0.042371,0.0064283,-0.011575,-0.0013839,0.010978,1.206,3.1454,29.689};
    a[0]=log(90)/0.2821;
    int padding=3;
    a[al]=log(80)/0.2821-a[0];
    a[2*al]=log(70)/0.2821-a[0];
    double *dvda=calloc((3*al+padding)*(3*nvert+padding),sizeof(double));
    octantoid_to_trimesh(a,LMAX,nrows,tlist,vlist,dvda,padding);
    print_matrix(dvda,3*nvert+padding,3*al+padding);
    //write_shape_file("/tmp/ashape.txt",tlist,vlist,nfac,nvert);
   // write_matrix_file("/tmp/dvda.txt",dvda,3*nvert,3*al);
}
*/
