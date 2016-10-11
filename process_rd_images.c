#include"utils.h"
#include<stdio.h>
#include"structs.h"

RDstruct * process_rd_images(char **filenames,int nRD,int *x0,int *y0,int *nx,int *ny, double *dx,double *dy,double *dates,double *RadarFreq,double min_tim,double *E,double *TIME,int nephm,int *LowFreq)
{
    /*INPUT:
     * filenames: filenames of AO images 
     * psfnames: filenames of psf images (optional)
     * nao: number of images
     * x0,y0,nx,ny define a rectangle in current image 
     * dx,dy pixel sizes
     * dates: observation dates corresponding to images. If nan, use date from image
     * E observation directions
     * TIME times corresponding to directions
     */
   
    RDstruct *RD=calloc(1,sizeof(RDstruct));
    RD->nRD=nRD;
    RD->nobs=calloc(nRD,sizeof(int));
    RD->datar=calloc(nRD,sizeof(double*));
    RD->datai=calloc(nRD,sizeof(double*));
//     AO->psfr=calloc(nao,sizeof(double*));
//     AO->psfi=calloc(nao,sizeof(double*));
    RD->freqx=calloc(nRD,sizeof(double*));
    RD->freqy=calloc(nRD,sizeof(double*));
    RD->E=calloc(3*nRD,sizeof(double));
    
    RD->TIME=calloc(nRD,sizeof(double));
    RD->distance=calloc(nRD,sizeof(double));
    RD->scalex=calloc(nRD,sizeof(double));
    RD->scaley=calloc(nRD,sizeof(double));
    RD->rfreq=calloc(nRD,sizeof(double));
    memcpy(RD->rfreq,RadarFreq,nRD*sizeof(double));
    double *datar,*datai,*freqx,*freqy;
    double cdelt1=0,cdelt2=0;
    double cao=1.731446326742403e+02; // c in AU/days
    double *buffer;
    double MaxFreqx,MaxFreqy;
//     double *psfbuffer;
    
    double date;
    int xsize,ysize;
    int FreqCount=0;
    int ntotal=0;
    double *zMr,*zMi,*Fx,*Fy;
    double norm1,norm2;
//     double upr[3];
    int index=0;
    //For each image, fill RD structure
    for(int j=0;j<nRD;j++)
    {
      
        
//         if(psfnames[j]!=NULL)
//         {
//             readfits(psfnames[j],&buffer,xp0[j],yp0[j],npx[j],npy[j],&date,&xsize,&ysize);
//             AO->psfr[j]=calloc(xsize*ysize/2,sizeof(double));
//             AO->psfi[j]=calloc(xsize*ysize/2,sizeof(double));
//             freqx=calloc(xsize*ysize/2,sizeof(double));
//             freqy=calloc(xsize*ysize/2,sizeof(double));
//             calc_image_fft(buffer,xsize,ysize,dx[j],dy[j],AO->psfr[j],AO->psfi[j],freqx,freqy);
//             free(freqx);
//             free(freqy);
//             free(buffer);
//         }
//         else
//         {
//             AO->psfr[j]=NULL;
//             AO->psfi[j]=NULL;
//         }
        readfits_rd(filenames[j],&buffer,x0[j],y0[j],nx[j],ny[j],&date,&xsize,&ysize,&cdelt1,&cdelt2); 
       // printf("%7.4f\n",date);
        if(dx[j]==0)
        {
            if(cdelt1>0)
                dx[j]=cdelt1;
            else
            {
                fprintf(stderr,"pixel size of radar image %d is undefined\n",j);
                exit(1);
            }
        }
         if(dy[j]==0)
        {
            if(cdelt2>0)
                dy[j]=cdelt2;
            else
            {
                fprintf(stderr,"pixel size of radar image %d is undefined\n",j);
                exit(1);
            }
        } 
          RD->scalex[j]=dx[j];
        RD->scaley[j]=dy[j];
        datai=calloc(xsize*ysize/2,sizeof(double));
        datar=calloc(xsize*ysize/2,sizeof(double));
        freqx=calloc(xsize*ysize/2,sizeof(double));
            freqy=calloc(xsize*ysize/2,sizeof(double));
        RD->datar[j]=calloc(xsize*ysize/2,sizeof(double));
        RD->datai[j]=calloc(xsize*ysize/2,sizeof(double));
        RD->freqx[j]=calloc(xsize*ysize/2,sizeof(double));
        RD->freqy[j]=calloc(xsize*ysize/2,sizeof(double));
       //Note that we reverse order of dy and dx here
        calc_image_fft_unnormed(buffer,xsize,ysize,dy[j],dx[j],datar,datai,freqx,freqy);
//         printf("xsize:%d ysize %d dx: %f dy: %f\n",xsize,ysize,dx[j],dy[j]);
       // printf("freqx: %f %f freqy: %f %f\n",RD->freqx[0][0],RD->freqx[0][1],RD->freqy[0][0],RD->freqy[0][1]);
        free(buffer);
        if(LowFreq!=NULL && LowFreq[j]==1)
        {
            FreqCount=0;
            MaxFreqx=1.0/(4*dy[j]); //Freqx corresponds to frequency here
            MaxFreqy=1.0/(4*dx[j]);
            for(int k=0;k<xsize*ysize/2-1;k++)
            {
                if(fabs(freqx[k])<=MaxFreqx && fabs(freqy[k])<=MaxFreqy)
                {
                    RD->freqx[j][FreqCount]=freqx[k];
                    RD->freqy[j][FreqCount]=freqy[k];
                    RD->datar[j][FreqCount]=datar[k];
                    RD->datai[j][FreqCount]=datai[k];
                     
                    FreqCount++;
                }
            }
             RD->nobs[j]=FreqCount;
        ntotal+=FreqCount;
        }
        else
        {
        RD->nobs[j]=xsize*ysize/2-1;
        ntotal+=xsize*ysize/2-1;
        memcpy(RD->datar[j],datar,(xsize*ysize/2-1)*sizeof(double));
        memcpy(RD->datai[j],datai,(xsize*ysize/2-1)*sizeof(double));
        memcpy(RD->freqx[j],freqx,(xsize*ysize/2-1)*sizeof(double));
        memcpy(RD->freqy[j],freqy,(xsize*ysize/2-1)*sizeof(double));
        }
        if(!isnan(dates[j]))
        {
            RD->TIME[j]=dates[j];
        }
        else if(!isnan(date))
            RD->TIME[j]=date;
        else
        {
            fprintf(stderr,"NO  JDMEAN in fits file and valid DATE was not given in ini\n");
            exit(-1);
        }
        //Find TIME and observation directions corresponding to current AO image
       
        index=find_index(TIME,nephm,RD->TIME[j]);
        norm1=NORM((E+3*index));
        
        //Correct for LT
       // RD->TIME[j]-=(norm1/cao+min_tim);
        //NO LT CORRECTION:
         RD->TIME[j]-=(min_tim);
        RD->distance[j]=norm1;
        for(int k=0;k<3;k++)
        {
            RD->E[3*j+k]=E[3*index+k]/norm1;
           
            
        }
       
       free(datar);
        free(datai);
        free(freqx);
        free(freqy); 
    }
    RD->ntotal=ntotal;
    return RD;
}
/*
int main()
{
    char file1[]="Metis/2.fits";
    char file2[]="Metis/1.fits";
    int nx[]={150,150};
    int ny[]={150,150};
    char **files=calloc(2,sizeof(char*));
    files[0]=file1;
    files[1]=file2;
    int zerovec[]={0,0};
    int fiftyvec[]={51,51};
    double dates[2];
    dates[0]=NAN;
    dates[1]=NAN;
    double pixscale[]={0.009942,0.009942};
    char **psfnames=calloc(2,sizeof(char*));
    char ephmfile[]="Metis/ephm.dat";
    int nao=2;
    AOstruct *result;
    int nobs=2,reads;
    double *E=calloc(nobs*3,sizeof(double));
    double *E0=calloc(nobs*3,sizeof(double));
    double *TIME=calloc(nobs,sizeof(double));
    double up[]={0,0,1,0,0,1};
    reads=read_ephm_data(ephmfile,TIME,E,E0);
      printf("dates: %f %f\n",TIME[0],TIME[1]);
    result=process_ao_images(files,psfnames,nao,fiftyvec,fiftyvec,nx,ny,zerovec,zerovec,nx,ny,pixscale,pixscale,dates,0,E,E0,up,TIME);
  //  printf("dates: %f %f\n",result->TIME[0],result->TIME[1]);
   // printf("distance: %f %f\n",result->distance[0],result->distance[1]);
   // printf("Number of points: %d\n",result->ntotal);
    write_matrix_file("/tmp/datar.txt",result->datar[0],1,result->nobs[0]);
    write_matrix_file("/tmp/datai.txt",result->datai[0],1,result->nobs[0]);
    write_matrix_file("/tmp/Fx.txt",result->freqx[0],1,result->nobs[0]);
    write_matrix_file("/tmp/Fy.txt",result->freqy[0],1,result->nobs[0]);
   // write_matrix_file("/tmp/psfr.txt",result->psfr[0],1,result->nobs[0]);
   //write_matrix_file("/tmp/psfi.txt",result->psfr[0],1,result->nobs[0]);
    
   char shapefile[]="shape2.txt";
    
    double ini_dims[]={90,90,70};
    double angles[]={66*PI/180,185*PI/180,24*2*PI*1.0/5.079177,0};
    double min_tim=2449830.78281;
    //double *offset=calloc(2,sizeof(double));
    int *tlist,*tlistn;
    double *vlist,*vlist2,*vlistn,*vlistn2,*D;
    double *angles2;
    int nvert,nfac,nfacn,nvertn;
    int nfreq=result->nobs[0];
    read_shape(shapefile,&tlist,&vlist,&nfac,&nvert,0); 
    Sqrt3_Subdiv(tlist,vlist,nfac,nvert,&tlistn,&vlistn,&nfacn,&nvertn,&D);
double complex *F0,*dFdx,*dFdy,*dFdz,*dFdA,*dFdoff;
F0=calloc(nfreq,sizeof(double complex));
 dFdx=(double complex*)calloc(nfreq*nvertn,sizeof(double complex));
 dFdy=(double complex*)calloc(nfreq*nvertn,sizeof(double complex));
 dFdz=(double complex*)calloc(nfreq*nvertn,sizeof(double complex));
 dFdA=(double complex*)calloc(nfreq*3,sizeof(double complex));
 dFdoff=(double complex*)calloc(nfreq*3,sizeof(double complex));
    Calculate_AO_deriv(tlistn,vlistn,nfacn,nvertn,angles,result->E,result->E0,result->up,result->TIME[0],result->distance[0],result->freqx[0],result->freqy[0],result->nobs[0],offset,F0,dFdx,dFdy,dFdz,dFdA,dFdoff);
AOstruct *AO=result;
 int nAOtotal=2*(AO->ntotal);
    int nAOcols=3*nvert+3+2*AO->nao;
    int nAOoffsets=2*(AO->nao);
    double *offset=calloc(2*(AO->nao),sizeof(double));
    
    double *AOout=calloc(nAOtotal,sizeof(double));
    double *AOdv=calloc(nAOtotal*nAOcols,sizeof(double));
    Calculate_AOs(tlistn,vlistn,nfacn,nvertn,angles,AO,offset,D,nvertn,nvert,NULL,AOout,AOdv,1);
    
}
*/
