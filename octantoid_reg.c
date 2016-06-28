#include"utils.h"
#include<stdlib.h>
#include<stdio.h>

void octantoid_reg(double *a,int LMAX,double *oreg,double *doreg)
{
    /*
     * Octantoid regularization
     * INPUT:
     * a octantoid representation
     * oreg pointer to double
     * doreg 3*(LMAX+1)^2 array
     */
    double reg=0;
    int al=pow(LMAX+1,2);
    int count=0;
    double *vec=calloc(al,sizeof(double));
    for(int i=0;i<=LMAX;i++)
        for(int j=-i;j<=i;j++)
        {
            vec[count]=(double)i/LMAX;
            count++;
        }
    double *xc=a;
    double *yc=a+al;
    double *zc=a+2*al;
    for(int i=0;i<al;i++)
    {
        reg+=vec[i]*yc[i]*yc[i]+vec[i]*zc[i]*zc[i];
        doreg[al+i]=2*vec[i]*yc[i];
        doreg[2*al+i]=2*vec[i]*zc[i];
    }
    *oreg=reg;
    free(vec);
}
/*
void main()
{
   int LMAX=2;
   int al=pow(LMAX+1,2);
   double *a;
   a=calloc(3*al,sizeof(double));
   for(int j=0;j<27;j++)
       a[j]=j;
   
   double oreg;
   
   double *doreg=calloc(3*al,sizeof(double));
   octantoid_reg(a,LMAX,&oreg,doreg);
   printf("oreg: %f\n",oreg);
   write_matrix_file("/tmp/doreg.txt",doreg,1,3*al);
}
*/
