#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"utils.h"

void generate_sphere(int nrows,int *tlist,double *vlist)
{
    /*
     * divides a unit sphere into triangles
     * vlist is (4*nrows^2+2)x3 matrix
     * tlist is 8*nrows^2 x 3 matrix
     */
    int *ifp;
    double x,y,z;
    
    double *t,*f;
    t=calloc(4*pow(nrows,2)+2,sizeof(double));
    f=calloc(4*pow(nrows,2)+2,sizeof(double));
    triangulate_sphere(nrows,t,f,tlist);
    for(int j=0;j<4*pow(nrows,2)+2;j++)
    {
    x=sin(t[j])*cos(f[j]);
    y=sin(t[j])*sin(f[j]);
    z=cos(t[j]);
    vlist[3*j]=x;
    vlist[3*j+1]=y;
    vlist[3*j+2]=z;
    }
    free(t);
    free(f);
}
void generate_ellipsoid(int nrows,double a,double b,double c,int *tlist,double *vlist)
{
    /*
     * divides a unit sphere into triangles
     * vlist is (4*nrows^2+2)x3 matrix
     * tlist is 8*nrows^2 x 3 matrix
     */
    int *ifp;
    double x,y,z;
    
    double *t,*f;
    t=calloc(4*pow(nrows,2)+2,sizeof(double));
    f=calloc(4*pow(nrows,2)+2,sizeof(double));
    triangulate_sphere(nrows,t,f,tlist);
    for(int j=0;j<4*pow(nrows,2)+2;j++)
    {
    x=sin(t[j])*cos(f[j]);
    y=sin(t[j])*sin(f[j]);
    z=cos(t[j]);
    vlist[3*j]=a*x;
    vlist[3*j+1]=b*y;
    vlist[3*j+2]=c*z;
    }
    free(t);
    free(f);
}
/*
void main()
{
    int *tlist;
    int nrows=4;
    double *vlist;
    int nfac=8*pow(nrows,2);
    int nvert=4*pow(nrows,2)+2;
    tlist=calloc(3*nfac,sizeof(int));
    vlist=calloc(4*nvert,sizeof(double));
   generate_ellipsoid(nrows,3,2,1,tlist,vlist);
   write_shape_file("/tmp/shape.txt",tlist,vlist,nfac,nvert);
}
  */  
