#include<stdlib.h>
#include<stdio.h>

void write_matrix_file(char * str,double *M,int m,int n)
{
    FILE *fp;
    fp=fopen(str,"w");
    if(fp==NULL)
        exit(-1);
    for(int j=0;j<m;j++)
    {
        for(int k=0;k<n;k++)
        {
            fprintf(fp,"%.5f ",M[j*n+k]);
        }
        fprintf(fp,"\n");
    }
}
int main()
{
    char fn[]="/tmp/test.txt";
    double M[]={1,2,3,4,5,6,7,8,9};
    write_matrix_file(fn,M,3,3);
}
