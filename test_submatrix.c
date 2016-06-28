#include"utils.h"
#include<stdio.h>
#include<stdlib.h>
void print_submatrix(double *M,int m,int n,int k1,int l1,int k2,int l2)
{
    /*Print submatrix of mxn matrix*/
    for(int j=0;j<k2;j++)
    {
        printf("\n");
        for(int i=0;i<l2;i++)
            printf(" %.4f ",M[(j+k1)*n+i+l1]);
    }
    printf("\n");
}
int main()
{
    double *D;
    D=calloc(5*5,sizeof(double));
    for(int j=0;j<25;j++)
        D[j]=j;
    print_matrix(D,5,5);
    print_submatrix(D,5,5,0,0,3,3);
    print_submatrix(D,5,5,2,3,3,2);
}
