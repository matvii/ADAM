#include"utils.h"
#include<stdio.h>
#include<stdlib.h>
int read_ephm_data(char *filename,double *TIME,double *E,double *E0)
{
    
    double cao=1.731446326742403e+02;
    FILE *fid;
    char *buffer=calloc(2048,sizeof(char));
    double E2[3];
    double E02[3];
    double time;
    double lttime;
    int count=0;
    int nreads=0;
    fid=fopen(filename,"r");
     if(fid==NULL)
    {
        fprintf(stderr,"Error opening ephm file %s\n",filename);
        exit(-1);
    }
    while(fgets(buffer,2048,fid))
    {
        nreads=sscanf(buffer,"%lf %lf %lf %lf %lf %lf %lf",&time,E02,E02+1,E02+2,E2,E2+1,E2+2);
    
        if(nreads==7)
        {
            E0[3*count]=E02[0];
            E0[3*count+1]=E02[1];
            E0[3*count+2]=E02[2];
            E[3*count]=E2[0];
            E[3*count+1]=E2[1];
            E[3*count+2]=E2[2];
            TIME[count]=time;
        }
        else if(nreads==4)
        {
            E[3*count]=E02[0];
            E[3*count+1]=E02[1];
            E[3*count+2]=E02[2];
            TIME[count]=time;  
        }
        else
        {
           break;
        }
        count++;
    }
    
    fclose(fid);
    free(buffer);
    return count;
}
            
