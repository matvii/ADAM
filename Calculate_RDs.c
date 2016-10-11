#include "utils.h"
#include"num_of_threads.h"
#include"matrix_ops.h"
void Calculate_Range_Doppler(int *tlist,double *vlist,int nfac,int nvert,double *angles,double *Eo,double TIME,double *freqx,double *freqy,int nfreq,double rfreq,double *offset,double scal,double rexpe,double * Fr,double *Fi);
void Calculate_Range_Doppler_deriv(int *tlist,double *vlist,int nfac,int nvert,double *angles,double *Eo,double TIME,double *freqx,double *freqy,int nfreq,double rfreq,double *offset,double scal,double rexpe,double * Fr,double  *Fi,double *FTdxr,double *FTdxi,double *FTdyr,double *FTdyi,double *FTdzr,double *FTdzi,double *FTdAr,double *FTdAi,double *FTdoffr,double *FTdoffi,double *FTdexpr,double *FTdexpi);


void Calculate_RDs(int *tlist,double *vlist,int nfac,int nvert,double *angles,RDstruct  *RDs,double *offset,double *D,int dm,int dn,double *Weight,double *scale,double rexp,double *FT,double *FTdv,double *FTdoff,double *FTdsc,double *FTdxp,int deriv)
{
 /*Same as the original, only exception is the inclusion of matrix D (For effective memory usage)
  */
 int DisNULL=0;
 int D1V=0;
 int D3V=0;
 int UseScale=0;
 int UseWeight=0;
 if(scale!=NULL)
     UseScale=1;
 int nRD;
  nRD=RDs->nRD; //Number of RD images 
 /*First some sanity checking*/
 if(D==NULL)
     DisNULL=1;
 
 if(!DisNULL && nvert!=dm)
 {
     puts("Error: nvert is not equal dm.");
     exit(1);
 }
 if(Weight!=NULL)
     UseWeight=1;
 
  int *nopoints,*cumpoints,ntpoints;
  nopoints=RDs->nobs; //Array, number of samples in each RD image
  cumpoints=malloc((nRD+1)*sizeof(int));
  cumpoints[0]=0;
  for(int i=1;i<=nRD;i++)
      cumpoints[i]=cumpoints[i-1]+nopoints[i-1]; //cumpoints is  the cumulative sum of all observation points, used for indexing
  
    ntpoints=cumpoints[nRD];//Total number of points
 
  
  
 
  if(deriv==0)
  {  
omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
for(int obsind=0;obsind<nRD;obsind++)
  {
    double *FTE,*FTTIME,*FTfreqx,*FTfreqy,*FTrfreq,*datar,*datai;
    double  *FTr,*FTi;
    double W;
    if(UseWeight==1)
        W=Weight[obsind];
    else
        W=1;
    FTr=calloc(nopoints[obsind],sizeof(double));
   FTi=calloc(nopoints[obsind],sizeof(double));
   
    FTE=RDs->E+3*obsind;
    FTTIME=RDs->TIME+obsind;
    FTfreqx=RDs->freqx[obsind];
    FTfreqy=RDs->freqy[obsind];
    FTrfreq=RDs->rfreq+obsind;
    datar=RDs->datar[obsind];
    datai=RDs->datai[obsind];
 
    Calculate_Range_Doppler(tlist,vlist,nfac,nvert,angles,FTE,*FTTIME,FTfreqx,FTfreqy,nopoints[obsind],*FTrfreq,offset+2*obsind,scale[obsind],rexp,FTr,FTi);
     for(int j=0;j<nopoints[obsind];j++)
  {
    FT[j+cumpoints[obsind]]=W*(datar[j]-FTr[j]);
    FT[j+cumpoints[obsind]+ntpoints]=W*(datai[j]-FTi[j]);
  }
  //return;
    free(FTr);
    free(FTi);
  }
  
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

zero_array(FTdv,2*ntpoints*(3*nvertf+3));
zero_array(FTdoff,2*ntpoints*2*nRD);
zero_array(FTdsc,2*ntpoints*nRD);
zero_array(FTdxp,2*ntpoints);


  
  

  omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
for(int obsind=0;obsind<nRD;obsind++)
  {
    int cind=0;
    int oind=0;
    double *FTdxr,*FTdxi,*FTdyr,*FTdyi,*FTdzr,*FTdzi,*FTdAr,*FTdAi,*FTdoffr,*FTdoffi,*FTdexpr,*FTdexpi,*FTdxfr,*FTdxfi,*FTdyfr,*FTdyfi,*FTdzfr,*FTdzfi;
     double *FTE,*FTTIME,*FTfreqx,*FTfreqy,*FTrfreq;
    double *FTr,*FTi,*datar,*datai;
    double W;
    if(UseWeight==1)
        W=Weight[obsind];
    else
        W=1;
   //  obsind=omp_get_thread_num();
    FTr=calloc(nopoints[obsind],sizeof(double));
    FTi=calloc(nopoints[obsind],sizeof(double));
   FTdxr=calloc(nopoints[obsind]*nvertf,sizeof(double));
   FTdxi=calloc(nopoints[obsind]*nvertf,sizeof(double));
  FTdyr=calloc(nopoints[obsind]*nvertf,sizeof(double));
  FTdyi=calloc(nopoints[obsind]*nvertf,sizeof(double));
  FTdzr=calloc(nopoints[obsind]*nvertf,sizeof(double));
  FTdzi=calloc(nopoints[obsind]*nvertf,sizeof(double));
  FTdAr=calloc(nopoints[obsind]*3,sizeof(double));
  FTdAi=calloc(nopoints[obsind]*3,sizeof(double));
  FTdoffr=calloc(nopoints[obsind]*2,sizeof(double));
  FTdoffi=calloc(nopoints[obsind]*2,sizeof(double));
  FTdexpr=calloc(nopoints[obsind],sizeof(double));
  FTdexpi=calloc(nopoints[obsind],sizeof(double));
  
   FTE=RDs->E+3*obsind;
    FTTIME=RDs->TIME+obsind;
    FTfreqx=RDs->freqx[obsind];
    FTfreqy=RDs->freqy[obsind];
    FTrfreq=RDs->rfreq+obsind;
    datar=RDs->datar[obsind];
    datai=RDs->datai[obsind];
    if(D!=NULL)
    {
      FTdxfr=calloc(nopoints[obsind]*nvert,sizeof(double));
      FTdyfr=calloc(nopoints[obsind]*nvert,sizeof(double));
      FTdzfr=calloc(nopoints[obsind]*nvert,sizeof(double));
      FTdxfi=calloc(nopoints[obsind]*nvert,sizeof(double));
      FTdyfi=calloc(nopoints[obsind]*nvert,sizeof(double));
      FTdzfi=calloc(nopoints[obsind]*nvert,sizeof(double));
      
      Calculate_Range_Doppler_deriv(tlist,vlist,nfac,nvert,angles,FTE,*FTTIME,FTfreqx,FTfreqy,nopoints[obsind],*FTrfreq,offset+2*obsind,scale[obsind],rexp,FTr,FTi,FTdxfr,FTdxfi,FTdyfr,FTdyfi,FTdzfr,FTdzfi,FTdAr,FTdAi,FTdoffr,FTdoffi,FTdexpr,FTdexpi);
      //Convert from vlistn->vlist. Only because we want to minimize memory usage
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
    else
      Calculate_Range_Doppler_deriv(tlist,vlist,nfac,nvert,angles,FTE,*FTTIME,FTfreqx,FTfreqy,nopoints[obsind],*FTrfreq,offset+2*obsind,scale[obsind],rexp,FTr,FTi,FTdxr,FTdxi,FTdyr,FTdyi,FTdzr,FTdzi,FTdAr,FTdAi,FTdoffr,FTdoffi,FTdexpr,FTdexpi);
     for(int j=0;j<nopoints[obsind];j++)
  {
    FT[j+cumpoints[obsind]]=W*(datar[j]-FTr[j]);
    FT[j+cumpoints[obsind]+ntpoints]=W*(datai[j]-FTi[j]);
  }
//    print_matrix(vlist,10,3);
// print_matrix(angles,1,3);
// print_matrix(offset,1,2);
// printf("Scale: %f rexp:%f TIME: %f, E: %f %f %f\n",scale[0],rexp,*FTTIME,FTE[0],FTE[1],FTE[2]);
// write_matrix_file("/tmp/FTfreqx.txt",FTfreqx,1,nopoints[obsind]);
// write_matrix_file("/tmp/FTfreqy.txt",FTfreqy,1,nopoints[obsind]);
 // write_matrix_file("/tmp/FTr.txt",FTr,1,nopoints[obsind]);
  // write_matrix_file("/tmp/FTi.txt",FTi,1,nopoints[obsind]);
  //Copy variables to matlab
  cind=cumpoints[obsind];
  oind=nopoints[obsind];
  if(UseWeight==1)
  {
      mult_with_cons(FTdxr,oind,dn,W);
      mult_with_cons(FTdxi,oind,dn,W);
      mult_with_cons(FTdyr,oind,dn,W);
      mult_with_cons(FTdyi,oind,dn,W);
      mult_with_cons(FTdzr,oind,dn,W);
      mult_with_cons(FTdzi,oind,dn,W);
      mult_with_cons(FTdAr,oind,3,W);
      mult_with_cons(FTdAi,oind,3,W);
      mult_with_cons(FTdoffr,oind,2,W);
      mult_with_cons(FTdoffi,oind,2,W);
      mult_with_cons(FTr,oind,1,W);
    mult_with_cons(FTi,oind,1,W);
     mult_with_cons(FTdexpr,oind,1,W);
     mult_with_cons(FTdexpi,oind,1,W);
  }
  set_submatrix(FTdv,2*ntpoints,3*dn+3,FTdxr,oind,dn,cind,0);
  set_submatrix(FTdv,2*ntpoints,3*dn+3,FTdxi,oind,dn,cind+ntpoints,0);
  
  set_submatrix(FTdv,2*ntpoints,3*dn+3,FTdyr,oind,dn,cind,dn);
  set_submatrix(FTdv,2*ntpoints,3*dn+3,FTdyi,oind,dn,cind+ntpoints,dn);
  
  set_submatrix(FTdv,2*ntpoints,3*dn+3,FTdzr,oind,dn,cind,2*dn);
  set_submatrix(FTdv,2*ntpoints,3*dn+3,FTdzi,oind,dn,cind+ntpoints,2*dn);
  
  set_submatrix(FTdv,2*ntpoints,3*dn+3,FTdAr,oind,3,cind,3*dn);
  set_submatrix(FTdv,2*ntpoints,3*dn+3,FTdAi,oind,3,cind+ntpoints,3*dn);
  
  set_submatrix(FTdoff,2*ntpoints,2*nRD,FTdoffr,oind,2,cind,2*obsind);
  set_submatrix(FTdoff,2*ntpoints,2*nRD,FTdoffi,oind,2,cind+ntpoints,2*obsind);
  
  set_submatrix(FTdsc,2*ntpoints,nRD,FTr,oind,1,cind,obsind);
  set_submatrix(FTdsc,2*ntpoints,nRD,FTi,oind,1,cind+ntpoints,obsind);
  
  set_submatrix(FTdxp,2*ntpoints,1,FTdexpr,oind,1,cind,0);
  set_submatrix(FTdxp,2*ntpoints,1,FTdexpi,oind,1,cind+ntpoints,0);
  free(FTr);
  free(FTi);
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
free(FTdexpr);
free(FTdexpi);
  

}
}
/*
int main()
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
    int nRD=3;
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
    double rfreq[]={90e9,90e9,90e9};
    double norm=NORM(E2);
    //printf("norm: %f\n",norm);
    for(int j=0;j<3;j++)
        E[j+3]=E2[j]/norm;
    
    double TIME[]={0.1,0.2,-0.1};
    double distance[]={0.00137879506,0.00137879506,0.00137879506};
    //double scale[]={0.1,0.1,0.1};
    
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
    double angles[]={0,0,30,0};
    double offset[]={0,0,0,0,0,0};
    double D[]={1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1};
    double psf[]={0.5377,1.8339,-2.2588,    0.8622    ,0.3188  , -1.3077,   -0.4336,    0.3426 ,   3.5784,    2.7694,   -1.3499,    3.0349,    0.7254,   -0.0631  ,  0.7147,   -0.2050,   -0.1241,    1.4897 ,   1.4090,    1.4172,    0.6715 ,  -1.2075,    0.7172 ,   1.6302,    0.4889 ,   1.0347,    0.7269, -0.3034,    0.2939};
    double *Weight;
    double *FT,*FTdv;
    double *FTdoff=calloc(2*ntpoints*2*nRD,sizeof(double));
    double *FTdsc=calloc(2*ntpoints*nRD,sizeof(double));
    double *FTdxp=calloc(2*ntpoints,sizeof(double));
    FT=calloc(2*ntpoints,sizeof(double));
    FTdv=calloc(2*ntpoints*(3*nvert+3),sizeof(double));
    RDstruct RD;
    RD.nRD=3;
    RD.nobs=nobs;
    RD.datar=calloc(nRD,sizeof(double*));
    RD.datai=calloc(nRD,sizeof(double*));
    RD.freqx=calloc(nRD,sizeof(double*));
    RD.freqy=calloc(nRD,sizeof(double*));
    
    RD.datar[0]=datar;
    RD.datai[0]=datai;
    RD.freqx[0]=freqx;
    RD.freqy[0]=freqy;
    RD.datar[1]=datar;
    RD.datai[1]=datai;
    RD.freqx[1]=freqx;
    RD.freqy[1]=freqy;
    RD.datar[2]=datar;
    RD.datai[2]=datai;
    RD.freqx[2]=freqx;
    RD.freqy[2]=freqy;
    RD.rfreq=rfreq;
    RD.E=E;
    
    RD.TIME=TIME;
    RD.distance=distance;
   printf("nfac: %d nvert: %d\n",nfac,nvert);
    double scale[]={0.0,0.0,0.0};
   Calculate_RDs(tlist,vlist,nfac,nvert,angles,&RD,offset,NULL,nvert,nvert,scale,0.69,FT,FTdv,FTdoff,FTdsc,FTdxp,1);
    
    //Calculate_AOs(tlist,vlist,nfac,nvert,angles,&AO,offset,NULL,nvert,nvert,Weight,NULL,FT,FTdv,NULL,1);
    
    
    //print_matrix(FT,1,2*ntpoints);
    //print_matrix(FTdv,2*ntpoints,3*nvert+2*nao+3);
    write_matrix_file("/tmp/FT.txt",FT,2*ntpoints,1);
  write_matrix_file("/tmp/FTdv.txt",FTdv,2*ntpoints,3*nvert+3);
  write_matrix_file("/tmp/FTdoff.txt",FTdoff,2*ntpoints,2*nRD);
  write_matrix_file("/tmp/FTdsc.txt",FTdsc,2*ntpoints,nRD);
    free(FT);
    free(FTdv);
   
}
*/
