#include "num_of_threads.h"
#include<omp.h>
#include"utils.h"
#include"structs.h"
#include"matrix_ops.h"


void Calculate_AOs(int *tlist,double *vlist,int nfac,int nvert,double *angles,AOstruct *AOs,double *offset,double *D,int dm,int dn,double *Weight,double *scale,double *FT,double *FTdv,double* FTdS,int deriv)
{
 /*tlist,vlist,angles -the asteroid shape
  * AO struct contains the AO data
  * offset naox2 vector, offsets,
  * D is the derivative matrix (dm x dn), derivatives of vertex coordinates wrt parameters 
  * Weight is additional weighting terms for individual AO images, 1xnao vector (not implemented yet)
  * Scale additional scaling terms for each ao image.
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
  * NOTE THAT FTdv is assumed to be initialized to zero
  
  */
 /*TBD: Combine real and complex matrices here*/
 /*IN THIS VERSION D can be either nvert x dn or 3*nvert x dn matrix*/
 int DisNULL=0;
 int D1V=0;
 int D3V=0;
 int UseScale=0;
 if(scale!=NULL)
     UseScale=1;
 int nao;
  nao=AOs->nao; //Number of AO images 
 /*First some sanity checking*/
 if(D==NULL)
     DisNULL=1;
 int nparams;
 if(DisNULL==0)
 {
    if(dm==3*nvert)
    {
     D3V=1;
     nparams=dn;
    }
    else if(dm==nvert)
    {
        D1V=1;
        nparams=3*dn;
    }
    else
    {
        perror("if D is not null, then dm must be equal to either nvert or 3*nvert\n");
        exit(1);
    }
 }
     
 
 
 
  
  int M,N;
  int *nopoints,*cumpoints,ntpoints;
  nopoints=AOs->nobs; //Array, number of samples in each AO image
  cumpoints=malloc((nao+1)*sizeof(int));
  cumpoints[0]=0;
  for(int i=1;i<=nao;i++)
      cumpoints[i]=cumpoints[i-1]+nopoints[i-1]; //cumpoints is  the cumulative sum of all observation points, used for indexing
  
    ntpoints=cumpoints[nao];//Total number of points
  if(deriv==0)
  {
     
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
    Calculate_AO(tlist,vlist,nfac,nvert,angles,FTE,FTE0,FTup,*FTTIME,*FTdist,FTfreqx,FTfreqy,nopoints[obsind],offset+2*obsind,FTr,FTi);
   // printf("Time taken: %f\n",omp_get_wtime()-time);
  if(psfr==NULL || psfi==NULL)
  {
     for(int j=0;j<nopoints[obsind];j++)
  {
    
         FT[j+cumpoints[obsind]]=(datar[j]-Scale*FTr[j]);
    FT[j+cumpoints[obsind]+ntpoints]=(datai[j]-Scale*FTi[j]);
   
  }
  }
  else
  {
     for(int j=0;j<nopoints[obsind];j++)
  {
    
        FT[j+cumpoints[obsind]]=(datar[j]-Scale*(psfr[j]*FTr[j]-psfi[j]*FTi[j])); 
    FT[j+cumpoints[obsind]+ntpoints]=(datai[j]-Scale*(psfi[j]*FTr[j]+psfr[j]*FTi[j]));
  }
  }
 
    free(FTr);
    free(FTi);
  }
  free(cumpoints);
 return; 
}

int nvertf;
if(D!=NULL)
    nvertf=dn;
else
{
  nvertf=nvert;
  dn=nvert;
}

  
zero_array(FTdv,2*ntpoints*(3*dn+3+2*nao));
if(UseScale==1)
zero_array(FTdS,2*ntpoints*nao);
  
  

omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
for(int obsind=0;obsind<nao;obsind++)
  {
      
      double Scale=1;
     if(UseScale==1)
         Scale=exp(scale[obsind]);
      
    int cind=0;
    int oind=0;
    double *FTdxyz,*FTdar,*FTdai;
    double *FTdxfr,*FTdxfi,*FTdyfr,*FTdyfi,*FTdzfr,*FTdzfi;
     double *FTE,*FTE0,*FTTIME,*FTfreqx,*FTfreqy,*FTdist,*FTup,*datar,*datai;
     double *FTdAr,*FTdAi,*FTdoffr,*FTdoffi,*FTdxr,*FTdxi,*FTdyr,*FTdyi,*FTdzr,*FTdzi;
    double  *FTr,*FTi;
    double *psfr,*psfi;
    psfr=AOs->psfr[obsind];
    psfi=AOs->psfi[obsind];
   //  obsind=omp_get_thread_num();
    FTr=calloc(nopoints[obsind],sizeof(double));
    FTi=calloc(nopoints[obsind],sizeof(double));
  
  //TBD: This is a temporary solution, fix this!
    FTdAr=calloc(nopoints[obsind]*3,sizeof(double));
    FTdAi=calloc(nopoints[obsind]*3,sizeof(double));
    FTdoffr=calloc(nopoints[obsind]*2,sizeof(double));
    FTdoffi=calloc(nopoints[obsind]*2,sizeof(double));
    FTdxr=calloc(nopoints[obsind]*nvertf,sizeof(double));
    FTdxi=calloc(nopoints[obsind]*nvertf,sizeof(double));
    FTdyr=calloc(nopoints[obsind]*nvertf,sizeof(double));
    FTdyi=calloc(nopoints[obsind]*nvertf,sizeof(double));
    FTdzr=calloc(nopoints[obsind]*nvertf,sizeof(double));
    FTdzi=calloc(nopoints[obsind]*nvertf,sizeof(double));
   
 datar=AOs->datar[obsind];
   datai=AOs->datai[obsind];
   
   FTE=AOs->E+3*obsind;
   FTE0=AOs->E0+3*obsind;
   FTup=AOs->up+3*obsind;
   FTTIME=AOs->TIME+obsind;
   FTfreqx=AOs->freqx[obsind];
   FTfreqy=AOs->freqy[obsind];
   FTdist=AOs->distance+obsind;
    if(D1V==1)
    {
      FTdxfr=calloc(nopoints[obsind]*nvert,sizeof(double));
      FTdyfr=calloc(nopoints[obsind]*nvert,sizeof(double));
      FTdzfr=calloc(nopoints[obsind]*nvert,sizeof(double));
     FTdxfi=calloc(nopoints[obsind]*nvert,sizeof(double));
      FTdyfi=calloc(nopoints[obsind]*nvert,sizeof(double));
      FTdzfi=calloc(nopoints[obsind]*nvert,sizeof(double));
      //double time=omp_get_wtime();
      Calculate_AO_deriv(tlist,vlist,nfac,nvert,angles,FTE,FTE0,FTup, *FTTIME,*FTdist,FTfreqx,FTfreqy,nopoints[obsind],offset+2*obsind,FTr,FTi,FTdxfr,FTdxfi,FTdyfr,FTdyfi,FTdzfr,FTdzfi,FTdAr,FTdAi,FTdoffr,FTdoffi);
     // printf("Time taken: %f\n",omp_get_wtime()-time);
     
     
      //Convert from vlistn->vlist by multiplying with D
      matrix_prod(FTdxfr,nopoints[obsind],nvert,D,nvertf,FTdxr);
      matrix_prod(FTdxfi,nopoints[obsind],nvert,D,nvertf,FTdxi);
      free(FTdxfr);
      free(FTdxfi);
      matrix_prod(FTdyfr,nopoints[obsind],nvert,D,nvertf,FTdyr);
      matrix_prod(FTdyfi,nopoints[obsind],nvert,D,nvertf,FTdyi);
      free(FTdyfr);
      free(FTdyfi);
      matrix_prod(FTdzfr,nopoints[obsind],nvert,D,nvertf,FTdzr);
      matrix_prod(FTdzfi,nopoints[obsind],nvert,D,nvertf,FTdzi);
      free(FTdzfr);
      free(FTdzfi);
    }
    else if(D3V==1)
    {
         FTdxfr=calloc(nopoints[obsind]*nvert,sizeof(double));
      FTdyfr=calloc(nopoints[obsind]*nvert,sizeof(double));
      FTdzfr=calloc(nopoints[obsind]*nvert,sizeof(double));
     FTdxfi=calloc(nopoints[obsind]*nvert,sizeof(double));
      FTdyfi=calloc(nopoints[obsind]*nvert,sizeof(double));
      FTdzfi=calloc(nopoints[obsind]*nvert,sizeof(double));
      FTdxyz=calloc(nopoints[obsind]*3*nvert,sizeof(double));
      FTda=calloc(nopoints[obsind]*nparam,sizeof(double));
      //double time=omp_get_wtime();
      Calculate_AO_deriv(tlist,vlist,nfac,nvert,angles,FTE,FTE0,FTup, *FTTIME,*FTdist,FTfreqx,FTfreqy,nopoints[obsind],offset+2*obsind,FTr,FTi,FTdxfr,FTdxfi,FTdyfr,FTdyfi,FTdzfr,FTdzfi,FTdAr,FTdAi,FTdoffr,FTdoffi);
     // printf("Time taken: %f\n",omp_get_wtime()-time);
     
     
      //Convert from vertex derivatives to derivatives wrt params by multiplying with D
      set_submatrix(FTdxyz,nopoints[obsind],3*nvert,FTdxfr,nopoints[obsind],nvert,0,0);
      set_submatrix(FTdxyz,nopoints[obsind],3*nvert,FTdyfr,nopoints[obsind],nvert,0,nvert);
      set_submatrix(FTdxyz,nopoints[obsind],3*nvert,FTdzfr,nopoints[obsind],nvert,0,2*nvert);
      matrix_prod(FTdxyz,nopoints[obsind],3*nvert,D,nparams,FTdar);
      set_submatrix(FTdxyz,nopoints[obsind],3*nvert,FTdxfi,nopoints[obsind],nvert,0,0);
      set_submatrix(FTdxyzv,nopoints[obsind],3*nvert,FTdyfi,nopoints[obsind],nvert,0,nvert);
      set_submatrix(FTdxyz,nopoints[obsind],3*nvert,FTdzfi,nopoints[obsind],nvert,0,2*nvert);
      matrix_prod(FTdxyz,nopoints[obsind],3*nvert,D,nparams,FTdai);
      free(FTdxyz);
      free(FTdxfr);
      free(FTdxfi);
      free(FTdyfr);
      free(FTdyfi);
      free(FTdzfr);
      free(FTdzfi);
    }
    else
      Calculate_AO_deriv(tlist,vlist,nfac,nvert,angles,FTE,FTE0,FTup, *FTTIME,*FTdist,FTfreqx,FTfreqy,nopoints[obsind],offset+2*obsind,FTr,FTi,FTdxr,FTdxi,FTdyr,FTdyi,FTdzr,FTdzi,FTdAr,FTdAi,FTdoffr,FTdoffi);
 
  

  cind=cumpoints[obsind];
  oind=nopoints[obsind];
  if(psfr==NULL || psfi==NULL)
  {
  for(int j=0;j<oind;j++)
  {
    FTr[j]=FTr[j]*Scale;
    FTi[j]=FTi[j]*Scale;
    FT[j+cind]=(datar[j]-FTr[j]); 
    FT[j+cind+ntpoints]=(datai[j]-FTi[j]);
   
    
    
  }
  }
  else
  {
      
      double temp;
      for(int j=0;j<oind;j++)
      {
          temp=FTr[j];
          FTr[j]=Scale*(psfr[j]*FTr[j]-psfi[j]*FTi[j]);
          FTi[j]=Scale*(psfi[j]*temp+psfr[j]*FTi[j]);
          
          FT[j+cind]=(datar[j]-FTr[j]); 
          FT[j+cind+ntpoints]=(datai[j]-FTi[j]);
          
      }
      
      //Multiply derivatives with the psf
      if(D3V==0)
      {
          for(int k1=0;k1<oind;k1++)
          {
              for(int k2=0;k2<nvertf;k2++)
              {
                  temp=FTdxr[k1*nvertf+k2];
                  FTdxr[k1*nvertf+k2]=FTdxr[k1*nvertf+k2]*psfr[k1]-psfi[k1]*FTdxi[k1*nvertf+k2];
                  FTdxi[k1*nvertf+k2]=FTdxi[k1*nvertf+k2]*psfr[k1]+psfi[k1]*temp;
                  temp=FTdyr[k1*nvertf+k2];
                  FTdyr[k1*nvertf+k2]=FTdyr[k1*nvertf+k2]*psfr[k1]-psfi[k1]*FTdyi[k1*nvertf+k2];
                  FTdyi[k1*nvertf+k2]=FTdyi[k1*nvertf+k2]*psfr[k1]+psfi[k1]*temp;
                  temp=FTdzr[k1*nvertf+k2];
                  FTdzr[k1*nvertf+k2]=FTdzr[k1*nvertf+k2]*psfr[k1]-psfi[k1]*FTdzi[k1*nvertf+k2];
                  FTdzi[k1*nvertf+k2]=FTdzi[k1*nvertf+k2]*psfr[k1]+psfi[k1]*temp;
              }
              for(int k2=0;k2<3;k2++)
              {
                  temp=FTdAr[k1*3+k2];
                  FTdAr[k1*3+k2]=FTdAr[k1*3+k2]*psfr[k1]-psfi[k1]*FTdAi[k1*3+k2];
                  FTdAi[k1*3+k2]=FTdAi[k1*3+k2]*psfr[k1]+psfi[k1]*temp;
              }
              temp=FTdoffr[k1*2];
              FTdoffr[k1*2]=FTdoffr[k1*2]*psfr[k1]-psfi[k1]*FTdoffi[k1*2];
              FTdoffi[k1*2]=FTdoffi[k1*2]*psfr[k1]+psfi[k1]*temp;
              temp=FTdoffr[k1*2+1];
              FTdoffr[k1*2+1]=FTdoffr[k1*2+1]*psfr[k1]-psfi[k1]*FTdoffi[k1*2+1];
              FTdoffi[k1*2+1]=FTdoffi[k1*2+1]*psfr[k1]+psfi[k1]*temp;
          }
      }
      else
      {
         for(int k1=0;k1<oind;k1++)
      {
          for(int k2=0;k2<nparams;k2++)
          {
              temp=FTdar[k1*nparams+k2];
              FTdar[k1*nparams+k2]=FTdxr[k1*nparams+k2]*psfr[k1]-psfi[k1]*FTdxi[k1*nparams+k2];
              FTdai[k1*nparams+k2]=FTdxi[k1*nparams+k2]*psfr[k1]+psfi[k1]*temp;
             
          }
          for(int k2=0;k2<3;k2++)
          {
              temp=FTdAr[k1*3+k2];
              FTdAr[k1*3+k2]=FTdAr[k1*3+k2]*psfr[k1]-psfi[k1]*FTdAi[k1*3+k2];
              FTdAi[k1*3+k2]=FTdAi[k1*3+k2]*psfr[k1]+psfi[k1]*temp;
          }
          temp=FTdoffr[k1*2];
          FTdoffr[k1*2]=FTdoffr[k1*2]*psfr[k1]-psfi[k1]*FTdoffi[k1*2];
          FTdoffi[k1*2]=FTdoffi[k1*2]*psfr[k1]+psfi[k1]*temp;
          temp=FTdoffr[k1*2+1];
          FTdoffr[k1*2+1]=FTdoffr[k1*2+1]*psfr[k1]-psfi[k1]*FTdoffi[k1*2+1];
          FTdoffi[k1*2+1]=FTdoffi[k1*2+1]*psfr[k1]+psfi[k1]*temp;
      } 
      }
          
  }
  
  
  
 
/*Copy submatrices to the final matrix. This is a temporary solution. Streamline this to avoid unnecessary copying
 * FTdv is is 2*ntpoints x 3*dn+3+2*nao matrix
 */
if(UseScale==1)
{
    if(D3V==1)
    {
       mult_with_cons(FTdar,oind,nparams,Scale);
    mult_with_cons(FTdai,oind,nparams,Scale);
    }
    else
    {
      mult_with_cons(FTdxr,oind,nvertf,Scale);
    mult_with_cons(FTdxi,oind,nvertf,Scale);
    mult_with_cons(FTdyr,oind,nvertf,Scale);
    mult_with_cons(FTdyi,oind,nvertf,Scale);
    mult_with_cons(FTdzr,oind,nvertf,Scale);
    mult_with_cons(FTdzi,oind,nvertf,Scale);
    }
    mult_with_cons(FTdAr,oind,3,Scale);
    mult_with_cons(FTdAi,oind,3,Scale); 
    
    
    mult_with_cons(FTdAr,oind,3,Scale);
    mult_with_cons(FTdAi,oind,3,Scale);
    
    mult_with_cons(FTdoffr,oind,2,Scale);
    mult_with_cons(FTdoffi,oind,2,Scale);
    //derivatives wrt Scale
    
    set_submatrix(FTdS,2*ntpoints,nao,FTr,oind,1,cind,obsind);
    set_submatrix(FTdS,2*ntpoints,nao,FTi,oind,1,cind+ntpoints,obsind);
}  
free(FTr);
  free(FTi);
  if(D3V==1)
  {
      set_submatrix(FTdv,2*ntpoints,nparams+3+2*nao,FTdar,oind,nparams,cind,0);
set_submatrix(FTdv,2*ntpoints,nparams+3+2*nao,FTdai,oind,nparams,cind+ntpoints,0);


set_submatrix(FTdv,2*ntpoints,nparams+3+2*nao,FTdAr,oind,3,cind,nparams);
set_submatrix(FTdv,2*ntpoints,nparams+3+2*nao,FTdAi,oind,3,cind+ntpoints,nparams);

set_submatrix(FTdv,2*ntpoints,nparams+3+2*nao,FTdoffr,oind,2,cind,nparams+3+2*obsind);
set_submatrix(FTdv,2*ntpoints,nparams+3+2*nao,FTdoffi,oind,2,cind+ntpoints,nparams+3+2*obsind);

free(FTdar);
free(FTdai);

free(FTdAr);
free(FTdAi);
free(FTdoffr);
free(FTdoffi);
  }
  else
  {
set_submatrix(FTdv,2*ntpoints,3*dn+3+2*nao,FTdxr,oind,nvertf,cind,0);
set_submatrix(FTdv,2*ntpoints,3*dn+3+2*nao,FTdxi,oind,nvertf,cind+ntpoints,0);

set_submatrix(FTdv,2*ntpoints,3*dn+3+2*nao,FTdyr,oind,nvertf,cind,nvertf);
set_submatrix(FTdv,2*ntpoints,3*dn+3+2*nao,FTdyi,oind,nvertf,cind+ntpoints,nvertf);

set_submatrix(FTdv,2*ntpoints,3*dn+3+2*nao,FTdzr,oind,nvertf,cind,2*nvertf);
set_submatrix(FTdv,2*ntpoints,3*dn+3+2*nao,FTdzi,oind,nvertf,cind+ntpoints,2*nvertf);

set_submatrix(FTdv,2*ntpoints,3*dn+3+2*nao,FTdAr,oind,3,cind,3*nvertf);
set_submatrix(FTdv,2*ntpoints,3*dn+3+2*nao,FTdAi,oind,3,cind+ntpoints,3*nvertf);

set_submatrix(FTdv,2*ntpoints,3*dn+3+2*nao,FTdoffr,oind,2,cind,3*nvertf+3+2*obsind);
set_submatrix(FTdv,2*ntpoints,3*dn+3+2*nao,FTdoffi,oind,2,cind+ntpoints,3*nvertf+3+2*obsind);

free(FTdxr);
free(FTdxi);
free(FTdyr);
free(FTdyi);
free(FTdzr);
free(FTdzi);
free(FTdAr);
free(FTdAi);
free(FTdoffr);
free(FTdoffi);
  }
}
free(cumpoints);
}

/*
void main()
{
//     int tlist[]={1,2,3,
//         1,3,4,
//         2,4,3,
//         1,2,4}; //4 facets
//     double vlist[]={0.0,-2.0,0.0
//         ,0.5,0.0,-1.0,
//         0.0,1.0,1.0,
//         -3,1,4};
//     
//     int nvert=4;
//     int nfac=4;
    int *tlist;
    double *vlist;
    int nfac,nvert;
    char file[]="mshape.txt";
    read_shape(file,&tlist,&vlist,&nfac,&nvert,0);
    int nobs[]={29,29,29};
    int nao=3;
    int ntpoints=3*29;
    //double E[]={1,0,0};
    double E2[]={1,0.1,0.1};
    double E[9];
    E[0]=1;
    E[1]=0;
    E[2]=0;
    E[6]=1;
    E[7]=0;
    E[8]=0;
    double norm=NORM(E2);
    //printf("norm: %f\n",norm);
    for(int j=0;j<3;j++)
        E[j+3]=E2[j]/norm;
    double E0[]={1,0,0,1,0,0,1,0,0};
    double TIME[]={0.1,0.2,-0.1};
    double distance[]={0.00137879506,0.00137879506,0.00137879506};
    double scale[]={1,1,1};
    double up[]={0,0,1,0,0,1,0,0,1};
    double *datar=calloc(29,sizeof(double));
    double *datai=calloc(29,sizeof(double));
    double freqx[]={-1.0000,   -0.9300,   -0.8600,   -0.7900,   -0.7200,   -0.6500,   -0.5800,   -0.5100,   -0.4400,   -0.3700,   -0.3000,
        -0.2300,   -0.1600,   -0.0900,   -0.0200,    0.0500,    0.1200,    0.1900,    0.2600,    0.3300,    0.4000,
        0.4700,    0.5400,    0.6100,    0.6800,    0.7500,    0.8200,    0.8900,    0.9600};
    double freqy[]={1.2900,    1.2200,    1.1500,    1.0800,    1.0100,    0.9400,    0.8700,    0.8000,    0.7300,    0.6600,
        0.5900,    0.5200,    0.4500,    0.3800,    0.3100,    0.2400,    0.1700,    0.1000,    0.0300,   -0.0400,   -0.1100,
        -0.1800,   -0.2500,   -0.3200,   -0.3900,   -0.4600,   -0.5300,   -0.6000,   -0.6700,
    };
    double freqy2[]={-0.3,0.05};
    double freqy3[]={-0.5,-0.1};
    double freqx2[]={0.1,0.15};
    double angles[]={0.1,0.3,30,0};
    double offset[]={0.1,0.2,0.5,-0.1,0,-0.3};
    double D[]={1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1};
    double psf[]={0.5377,1.8339,-2.2588,    0.8622    ,0.3188  , -1.3077,   -0.4336,    0.3426 ,   3.5784,    2.7694,   -1.3499,    3.0349,    0.7254,   -0.0631  ,  0.7147,   -0.2050,   -0.1241,    1.4897 ,   1.4090,    1.4172,    0.6715 ,  -1.2075,    0.7172 ,   1.6302,    0.4889 ,   1.0347,    0.7269, -0.3034,    0.2939};
    double *Weight;
    double *FT,*FTdv;
    FT=calloc(2*ntpoints,sizeof(double));
    FTdv=calloc(2*ntpoints*(3*nvert+2*nao+3),sizeof(double));
    AOstruct AO;
    AO.nao=3;
    AO.nobs=nobs;
    AO.datar=calloc(nao,sizeof(double*));
    AO.datai=calloc(nao,sizeof(double*));
    AO.freqx=calloc(nao,sizeof(double*));
    AO.freqy=calloc(nao,sizeof(double*));
    AO.psfi=calloc(nao,sizeof(double*));
    AO.psfr=calloc(nao,sizeof(double*));
    AO.datar[0]=datar;
    AO.datai[0]=datai;
    AO.freqx[0]=freqx;
    AO.freqy[0]=freqy;
    AO.datar[1]=datar;
    AO.datai[1]=datai;
    AO.freqx[1]=freqx;
    AO.freqy[1]=freqy;
    AO.datar[2]=datar;
    AO.datai[2]=datai;
    AO.freqx[2]=freqx;
    AO.freqy[2]=freqy;
    AO.psfr[0]=psf;
    AO.psfi[0]=psf;
    AO.psfr[1]=psf;
    AO.psfi[1]=psf;
    AO.psfi[2]=psf;
    AO.psfr[2]=psf;
    AO.E=E;
    AO.E0=E0;
    AO.TIME=TIME;
    AO.distance=distance;
    AO.scalex=scale;
    AO.scaley=scale;
    AO.up=up;
   
    
    Calculate_AOs(tlist,vlist,nfac,nvert,angles,&AO,offset,NULL,nvert,nvert,Weight,NULL,FT,FTdv,NULL,1);
    
    
    //print_matrix(FT,1,2*ntpoints);
    //print_matrix(FTdv,2*ntpoints,3*nvert+2*nao+3);
    write_matrix_file("/tmp/FT.txt",FT,2*ntpoints,1);
  write_matrix_file("/tmp/FTdv.txt",FTdv,2*ntpoints,3*nvert+2*nao+3);
    free(FT);
    free(FTdv);
    free(AO.datar);
    free(AO.datai);
    free(AO.freqx);
    free(AO.freqy);
}

*/