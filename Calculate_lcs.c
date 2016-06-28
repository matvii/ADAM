#include"utils.h"
#include"matrix_ops.h"
#include"structs.h"
#include "num_of_threads.h"
#include<omp.h>
void calculate_lcurve2(int *tlist,double *vlist,int numfac,int numvert,double *angles,double *Eo,double *E0o,int nE,double *TIME,double *bright,double *dbrightx,double *dbrighty,double *dbrightz,double *dbrightb,double *dbrightl,double *dbrighto,double *A,double *Alimit,double *dA);

void calculate_lcs(int *tlist,double *vlist,int nfac,int nvert,double *angles,LCstruct *LC,double *D,int dm,int dn,double *LCout,double *dLCdv,double *Albedo,double *Alimit,double *dAlb,int deriv)
{
    /*Calculates the lightcurves corresponding to geometries described in LCstruct
     * Optionally Alb contains facet albedos, Alimit albedo limits.
     * OUTPUT:
     * LCout ntpoints array, where nlcp is the total number of lightcurve points
     * NOTE: LCout=lc_data-lc_model
     * dLCdv ntpoints x nvertf+3 matrix, (=dLCdv=dLCdx*D dLCdy*D dLCdz*D dLCdA)
     * dAlb derivatives wrt albedo, ntpoints x nfac array */
    
    
    
    if(D!=NULL && nvert!=dm)
    {
        puts("Error: Number of vertex coordinates is not equal to the number of rows in D.");
        exit(1);
    }
    int nvertf=dn;
    if(D==NULL)
        nvertf=nvert;
    int cind;
    int nlc,ntpoints,pinlc;
    int *cumpoints;
    int *nobs;
    nlc=LC->nlc;
    nobs=LC->nobs;
    cumpoints=malloc((nlc+1)*sizeof(int));
    cumpoints[0]=0;
    for(int i=1;i<=nlc;i++)
        cumpoints[i]=cumpoints[i-1]+nobs[i-1]; //Cumulative sum of points
        ntpoints=cumpoints[nlc]; //Total number of observed points
        omp_set_num_threads(1);
    printf("IN LCS, nlc is %d\n",nlc);
    #pragma omp parallel for
    for(int j=0;j<nlc;j++)
    {
         printf("HERE: j=%d\n",j);
        cind=cumpoints[j]; //total number of points in previous lightcurves
        pinlc=nobs[j]; //points in current lightcurve
        double *E,*E0,*TIME;
        double *bright;
        double *dbrightx;
        double *dbrighty;
        double *dbrightz;
        double *dbrightxf;
        double *dbrightyf;
        double *dbrightzf;
        double *dbrightb;
        double *dbrightl;
        double *dbrighto;
        double *dA;
        double *lcs;
        
        bright=calloc(pinlc,sizeof(double));
        dbrightx=calloc(pinlc*nvertf,sizeof(double));
        dbrighty=calloc(pinlc*nvertf,sizeof(double));
        dbrightz=calloc(pinlc*nvertf,sizeof(double));
        dbrightb=calloc(pinlc,sizeof(double));
        dbrightl=calloc(pinlc,sizeof(double));
        dbrighto=calloc(pinlc,sizeof(double));
        if(Albedo!=NULL)
            dA=calloc(pinlc*nfac,sizeof(double)); //If no albedo, then no albedo derivatives
            E=LC->E[j];
        E0=LC->E0[j];
        TIME=LC->TIME[j];
        lcs=LC->lcs[j];
        if(D!=NULL)
        {
            dbrightxf=calloc(pinlc*nvert,sizeof(double));
            dbrightyf=calloc(pinlc*nvert,sizeof(double));  
            dbrightyf=calloc(pinlc*nvert,sizeof(double));  
            calculate_lcurve2(tlist,vlist,nfac,nvert,angles,E,E0,pinlc,TIME,bright,dbrightxf,dbrightyf,dbrightzf,dbrightb,dbrightl,dbrighto,Albedo,Alimit,dA);
             printf("HERE: j=%d\n",j);
            matrix_prod(dbrightxf,pinlc,nvert,D,nvertf,dbrightx);
            matrix_prod(dbrightyf,pinlc,nvert,D,nvertf,dbrighty);
            matrix_prod(dbrightzf,pinlc,nvert,D,nvertf,dbrightz);
            free(dbrightxf);
            free(dbrightyf);
            free(dbrightzf);
        }
        else
            calculate_lcurve2(tlist,vlist,nfac,nvert,angles,E,E0,pinlc,TIME,bright,dbrightx,dbrighty,dbrightz,dbrightb,dbrightl,dbrighto,Albedo,Alimit,dA);
        /*Copy stuff to correct places*/
       
        for(int k=0;k<pinlc;k++)
            LCout[k+cumpoints[j]]=lcs[k]-bright[k];
        if(deriv==1)
        {
            /*Copy derivatives*/
            
            set_submatrix(dLCdv,ntpoints,3*nvertf+3,dbrightx,pinlc,nvertf,cind,0);
           
            set_submatrix(dLCdv,ntpoints,3*nvertf+3,dbrighty,pinlc,nvertf,cind,nvertf);
            set_submatrix(dLCdv,ntpoints,3*nvertf+3,dbrightz,pinlc,nvertf,cind,2*nvertf);
            
            set_submatrix(dLCdv,ntpoints,3*nvertf+3,dbrightb,pinlc,1,cind,3*nvertf);
            set_submatrix(dLCdv,ntpoints,3*nvertf+3,dbrightl,pinlc,1,cind,3*nvertf+1);
            set_submatrix(dLCdv,ntpoints,3*nvertf+3,dbrighto,pinlc,1,cind,3*nvertf+2);
            free(dbrightx);
            free(dbrighty);
            free(dbrightz);
            free(dbrightb);
            free(dbrightl);
            free(dbrighto);
            
            if(Albedo!=NULL)
            {
                set_submatrix(dAlb,ntpoints,nfac,dA,pinlc,nfac,cind,0);
                free(dA);
            }
            
        }
       
    }
    free(cumpoints);
}
/*
int main()
{
    
    int tlist[]={1,2,3};
    double vlist[]={0.0,-2.0,0.0,0.5,0.0,-1.0,0.0,1.0,1.0};
    double angles[]={0,0,30.15};
    int nvert=3;
    int nfac=1;
    int nobs[]={2,2};
    int nkcs=2;
    int ntpoints=4;
    //double E[]={1,0,0};
    double E2[]={1,0.1,0.1};
    double E[6];
    E[0]=1;
    E[1]=0;
    E[2]=0;
    double norm=NORM(E2);
    //printf("norm: %f\n",norm);
    for(int j=0;j<3;j++)
        E[j+3]=E2[j]/norm;
    double E0[6]={1,0,0,1,0,0};
    double TIME2[]={-0.03,0};
    double TIME1[]={0,0.03};
    double zeros[]={0.0,0.0};
    double *LCout;
    double *dLCdv;
    LCout=malloc(4*sizeof(double));
    dLCdv=calloc(4*(3*3+3),sizeof(double));
    LCstruct LC;
    LC.nlc=2;
    LC.nobs=nobs;
    LC.E=malloc(2*sizeof(double*));
    LC.E0=malloc(2*sizeof(double*));
    LC.lcs=malloc(2*sizeof(double*));
    LC.TIME=malloc(2*sizeof(double*));
    LC.E[0]=E;
    LC.E[1]=E;
    LC.E0[0]=E0;
    LC.E0[1]=E0;
    LC.TIME[0]=TIME2;
    LC.TIME[1]=TIME1;
    LC.lcs[0]=zeros;
    LC.lcs[1]=zeros;
    calculate_lcs(tlist,vlist,nfac,nvert,angles,&LC,NULL,0,0,LCout,dLCdv,NULL,NULL,NULL,1);
    print_matrix(LCout,1,4);
    print_matrix(dLCdv,4,3*nvert+3);
}

*/
