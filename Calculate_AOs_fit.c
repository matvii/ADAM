#include "num_of_threads.h"
#include<omp.h>
#include"utils.h"
#include"structs.h"
#include"matrix_ops.h"
#include"globals.h"

void Calculate_AOs_fit(int *tlist,double *vlist,int nfac,int nvert,double *angles,AOstruct *AOs,double *offset,double *D,int dm,int dn,double *Weight,double *scale,double *FT,double *FTdv,double* FTdS,double *Albedo,double *Alimit,double *dA,int deriv)
{
 /*tlist,vlist,angles -the asteroid shape
  * AO struct contains the AO data
  * offset naox2 vector, offsets,
  * D is the derivative matrix (dm x dn), derivatives of vertex coordinates wrt parameters 
  * Weight is additional weighting terms for individual AO images, 1xnao vector (not implemented yet)
  * Scale additional scaling terms for each ao image.
  * Albedo nfac vector, optional
  * Alimit, albedo limit 2-vector
  * deriv==1, then the derivatives will be calculated
  * OUTPUT:
  * FTr,FTi real and imaginary results
  * Derivative matrix FTdvr (real) FTdvi (imag)
  */
 /* Denote the total number of data points by ntpoints. Then
  * FT is 2*ntpoints vector
  * FTdv is 2*ntpoints x (3*dn+3+2*nao) matrix =[real(FTdx)*D real(FTdy)*D real(FTdz)*D real(FTdA) real(FTdoff);
  * imag(FTdx)*D imag(FTdy)*D imag(FTdz)*D imag(FTdA) imag(FTdoff);...]
  * FTdS is an optional matrix for Scaling terms
  * dAr, dAi 2*ntpointsxnfac matrices, only if albedo is to be fitted
  * NOTE THAT FTdv is assumed to be initialized to zero
  
  */
 /*TBD: Combine real and complex matrices here*/
 int DisNULL=0;
 int D1V=0;
 int D3V=0;
 int UseScale=0;
 int UseWeight=0;
 if(scale!=NULL)
     UseScale=1;
 int nao;
  nao=AOs->nao; //Number of AO images 
 /*First some sanity checking*/
 if(D==NULL)
     DisNULL=1;
 
 if(!DisNULL && nvert!=dm)
 {
     puts("Error: nvert is not equal dm.");
     exit(1);
 }
 double TB=0;
 if(Weight!=NULL)
     UseWeight=1;

  
  int M,N;
  int *nopoints,*cumpoints,ntpoints;
  nopoints=AOs->nobs; //Array, number of samples in each AO image
  cumpoints=malloc((nao+1)*sizeof(int));
  cumpoints[0]=0;
  for(int i=1;i<=nao;i++)
      cumpoints[i]=cumpoints[i-1]+nopoints[i-1]; //cumpoints is  the cumulative sum of all observation points, used for indexing
  
    ntpoints=cumpoints[nao];//Total number of points
  
     
omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
for(int obsind=0;obsind<nao;obsind++)
  {
      double Scale=1;
     if(UseScale==1)
         Scale=exp(scale[obsind]);
    double *FTE,*FTE0,*FTTIME,*FTfreqx,*FTfreqy,*FTup,*FTdist,*datar,*datai;
    double *FTr;
    double *FTi;
    double *psfi;
    double *psfr;
   double W;
   if(UseWeight==1)
       W=Weight[obsind];
   else
       W=1;
    FTr=calloc(nopoints[obsind],sizeof(double));
   FTi=calloc(nopoints[obsind],sizeof(double));
   FTE=AOs->E+3*obsind;
   FTE0=AOs->E0+3*obsind;
   FTup=AOs->up+3*obsind;
   FTTIME=AOs->TIME+obsind;
   FTfreqx=AOs->freqx[obsind];
   FTfreqy=AOs->freqy[obsind];
   FTdist=AOs->distance+obsind;
   datar=AOs->datar[obsind];
   datai=AOs->datai[obsind];
   psfr=AOs->psfr[obsind];
   psfi=AOs->psfi[obsind];
   
  // double time=omp_get_wtime();
    TB=Calculate_AO(tlist,vlist,nfac,nvert,angles,FTE,FTE0,FTup,*FTTIME,*FTdist,FTfreqx,FTfreqy,nopoints[obsind],offset+2*obsind,FTr,FTi,Albedo,Alimit);
   // printf("Time taken: %f\n",omp_get_wtime()-time);
  if(psfr==NULL || psfi==NULL)
  {
     for(int j=0;j<nopoints[obsind];j++)
  {
    
         FT[j+cumpoints[obsind]]=W*(datar[j]-Scale*FTr[j]*TB/INI_AO_TOTAL_BRIGHT[obsind]);
    FT[j+cumpoints[obsind]+ntpoints]=W*(datai[j]-Scale*FTi[j]*TB/INI_AO_TOTAL_BRIGHT[obsind]);
   
  }
  }
  else
  {
     for(int j=0;j<nopoints[obsind];j++)
  {
    
        FT[j+cumpoints[obsind]]=W*(datar[j]-Scale*(psfr[j]*FTr[j]-psfi[j]*FTi[j])*TB/INI_AO_TOTAL_BRIGHT[obsind]); 
    FT[j+cumpoints[obsind]+ntpoints]=W*(datai[j]-Scale*(psfi[j]*FTr[j]+psfr[j]*FTi[j])*TB/INI_AO_TOTAL_BRIGHT[obsind]);
  }
  }
 
    free(FTr);
    free(FTi);
  }
  free(cumpoints);
 
}

