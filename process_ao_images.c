#include"utils.h"
#include<stdio.h>
#include"structs.h"

AOstruct * process_ao_images(char **filenames,char **psfnames,int nao,int *x0,int *y0,int *nx,int *ny,int *xp0,int *yp0,int *npx,int *npy, double *dx,double *dy,double *dates,double min_tim,double *E,double *E0,double *up,double *TIME,int nephm,int *LowFreq)
{
    /*INPUT:
     * filenames: filenames of AO images 
     * psfnames: filenames of psf images (optional)
     * nao: number of images
     * x0,y0,nx,ny define a rectangle in current image 
     * dx,dy pixel sizes
     * dates: observation dates corresponding to images. If nan, use date from image
     * E,E0 illumination and observation directions
     * up camera directions
     * TIME times corresponding to directions
     */
   
    AOstruct *AO=calloc(1,sizeof(AOstruct));
    AO->nao=nao;
    AO->nobs=calloc(nao,sizeof(int));
    AO->datar=calloc(nao,sizeof(double*));
    AO->datai=calloc(nao,sizeof(double*));
    AO->psfr=calloc(nao,sizeof(double*));
    AO->psfi=calloc(nao,sizeof(double*));
    AO->freqx=calloc(nao,sizeof(double*));
    AO->freqy=calloc(nao,sizeof(double*));
    AO->E=calloc(3*nao,sizeof(double));
    AO->E0=calloc(3*nao,sizeof(double));
    AO->TIME=calloc(nao,sizeof(double));
    AO->distance=calloc(nao,sizeof(double));
    AO->scalex=calloc(nao,sizeof(double));
    AO->scaley=calloc(nao,sizeof(double));
    AO->up=calloc(3*nao,sizeof(double));
    double *datar,*datai,*psfr,*psfi;
    int count=0;
    double cao=1.731446326742403e+02; // c in AO/days
    double *buffer;
    double *psfbuffer;
    double *freqx,*freqy;
    double date;
    int xsize,ysize;
    int ntotal=0;
    double *zMr,*zMi,*Fx,*Fy;
    double norm1,norm2;
    double upr[3];
    double MaxFreqy=0;
    double MaxFreqx=0;
    int FreqCount=0;
    int index=0;
    //For each image, fill AO structure
    for(int j=0;j<nao;j++)
    {
        AO->scalex[j]=dx[j];
        AO->scaley[j]=dy[j];
        count=0;
      
   
        if(psfnames[j]!=NULL)
        {
            readfits(psfnames[j],&buffer,xp0[j],yp0[j],npx[j],npy[j],&date,&xsize,&ysize);
            AO->psfr[j]=calloc(xsize*ysize/2,sizeof(double));
            
            AO->psfi[j]=calloc(xsize*ysize/2,sizeof(double));
             psfr=calloc(xsize*ysize/2,sizeof(double));
            psfi=calloc(xsize*ysize/2,sizeof(double));
            freqx=calloc(xsize*ysize/2,sizeof(double));
            freqy=calloc(xsize*ysize/2,sizeof(double));
            calc_image_fft(buffer,xsize,ysize,dx[j],dy[j],psfr,psfi,freqx,freqy);
            free(freqx);
            free(freqy);
            free(buffer);
        }
        else
        {
            AO->psfr[j]=NULL;
            AO->psfi[j]=NULL;
        }
        readfits(filenames[j],&buffer,x0[j],y0[j],nx[j],ny[j],&date,&xsize,&ysize); 
       datai=calloc(xsize*ysize/2,sizeof(double));
        datar=calloc(xsize*ysize/2,sizeof(double));
        freqx=calloc(xsize*ysize/2,sizeof(double));
            freqy=calloc(xsize*ysize/2,sizeof(double));
        AO->datar[j]=calloc(xsize*ysize/2,sizeof(double));
        AO->datai[j]=calloc(xsize*ysize/2,sizeof(double));
        AO->freqx[j]=calloc(xsize*ysize/2,sizeof(double));
        AO->freqy[j]=calloc(xsize*ysize/2,sizeof(double));
       
        calc_image_fft(buffer,xsize,ysize,dx[j],dy[j],datar,datai,freqx,freqy);
        
        free(buffer);
        if(LowFreq!=NULL && LowFreq[j]==1)
        {
            FreqCount=0;
            MaxFreqx=1.0/(4*dx[j]);
            MaxFreqy=1.0/(4*dy[j]);
            for(int k=0;k<xsize*ysize/2-1;k++)
            {
                if(fabs(freqx[k])<=MaxFreqx && fabs(freqy[k])<=MaxFreqy)
                {
                    AO->freqx[j][FreqCount]=freqx[k];
                    AO->freqy[j][FreqCount]=freqy[k];
                    AO->datar[j][FreqCount]=datar[k];
                    AO->datai[j][FreqCount]=datai[k];
                     if(psfnames[j]!=NULL);
                     {
                     AO->psfr[j][FreqCount]=psfr[k];
                     AO->psfi[j][FreqCount]=psfi[k];
                     }
                    FreqCount++;
                }
            }
            AO->nobs[j]=FreqCount;
        ntotal+=FreqCount;
        }
        else
        {
        AO->nobs[j]=xsize*ysize/2-1;
        ntotal+=xsize*ysize/2-1;
        memcpy(AO->datar[j],datar,(xsize*ysize/2-1)*sizeof(double));
        memcpy(AO->datai[j],datai,(xsize*ysize/2-1)*sizeof(double));
        if(psfnames[j]!=NULL)
        {
        memcpy(AO->psfr[j],psfr,(xsize*ysize/2-1)*sizeof(double));
        memcpy(AO->psfi[j],psfi,(xsize*ysize/2-1)*sizeof(double));
        
        }
        memcpy(AO->freqx[j],freqx,(xsize*ysize/2-1)*sizeof(double));
        memcpy(AO->freqy[j],freqy,(xsize*ysize/2-1)*sizeof(double));
        }
        if(!isnan(dates[j]))
        {
            AO->TIME[j]=dates[j];
        }
        else if(!isnan(date))
            AO->TIME[j]=date+2400000.5;
        else
        {
            fprintf(stderr,"NO VALID MJD-OBS in fits file and valid DATE was not given in ini");
            exit(-1);
        }
        //Find TIME and observation directions corresponding to current AO image
       
        index=find_index(TIME,nephm,AO->TIME[j]);
        norm1=NORM((E+3*index));
        norm2=NORM((E0+3*index));
        //Correct for LT
        AO->TIME[j]-=(norm1/cao+min_tim);
        AO->distance[j]=norm1;
        for(int k=0;k<3;k++)
        {
            AO->E[3*j+k]=E[3*index+k]/norm1;
            AO->E0[3*j+k]=E0[3*index+k]/norm2;
            
        }
        if(up[3*j]==0 && up[3*j+1]==0 && up[3*j+2]==0)
        {
            calc_rot_frame(filenames[j],AO->E+3*j,upr);
            
            AO->up[3*j]=upr[0];
            AO->up[3*j+1]=upr[1];
            AO->up[3*j+2]=upr[2];
            
           // printf("Calculating up direction from file\n");
        }
        else
        {
            AO->up[3*j]=up[3*j];
            AO->up[3*j+1]=up[3*j+1];
            AO->up[3*j+2]=up[3*j+2];
        }
        if(psfnames[j]!=NULL)
        {
        free(psfi);
        free(psfr);
        }
        free(datar);
        free(datai);
        free(freqx);
        free(freqy);
        
    }
    AO->ntotal=ntotal;
    return AO;
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
    int *LowFreq=calloc(2,sizeof(double));
    LowFreq[0]=1;
    double up[]={0,0,1,0,0,1};
    reads=read_ephm_data(ephmfile,TIME,E,E0);
      printf("dates: %f %f\n",TIME[0],TIME[1]);
    result=process_ao_images(files,psfnames,nao,fiftyvec,fiftyvec,nx,ny,zerovec,zerovec,nx,ny,pixscale,pixscale,dates,0,E,E0,up,TIME,nao,LowFreq);
  //  printf("dates: %f %f\n",result->TIME[0],result->TIME[1]);
   // printf("distance: %f %f\n",result->distance[0],result->distance[1]);
   // printf("Number of points: %d\n",result->ntotal);
    write_matrix_file("/tmp/datar2.txt",result->datar[0],1,result->nobs[0]);
    write_matrix_file("/tmp/datai2.txt",result->datai[0],1,result->nobs[0]);
    write_matrix_file("/tmp/Fx2.txt",result->freqx[0],1,result->nobs[0]);
    write_matrix_file("/tmp/Fy2.txt",result->freqy[0],1,result->nobs[0]);
   // write_matrix_file("/tmp/psfr.txt",result->psfr[0],1,result->nobs[0]);
   //write_matrix_file("/tmp/psfi.txt",result->psfr[0],1,result->nobs[0]);
}    
//    char shapefile[]="shape2.txt";
//     
//     double ini_dims[]={90,90,70};
//     double angles[]={66*PI/180,185*PI/180,24*2*PI*1.0/5.079177,0};
//     double min_tim=2449830.78281;
//     //double *offset=calloc(2,sizeof(double));
//     int *tlist,*tlistn;
//     double *vlist,*vlist2,*vlistn,*vlistn2,*D;
//     double *angles2;
//     int nvert,nfac,nfacn,nvertn;
//     int nfreq=result->nobs[0];
//     read_shape(shapefile,&tlist,&vlist,&nfac,&nvert,0); 
//     Sqrt3_Subdiv(tlist,vlist,nfac,nvert,&tlistn,&vlistn,&nfacn,&nvertn,&D);
// double complex *F0,*dFdx,*dFdy,*dFdz,*dFdA,*dFdoff;
// F0=calloc(nfreq,sizeof(double complex));
//  dFdx=(double complex*)calloc(nfreq*nvertn,sizeof(double complex));
//  dFdy=(double complex*)calloc(nfreq*nvertn,sizeof(double complex));
//  dFdz=(double complex*)calloc(nfreq*nvertn,sizeof(double complex));
//  dFdA=(double complex*)calloc(nfreq*3,sizeof(double complex));
//  dFdoff=(double complex*)calloc(nfreq*3,sizeof(double complex));
//     Calculate_AO_deriv(tlistn,vlistn,nfacn,nvertn,angles,result->E,result->E0,result->up,result->TIME[0],result->distance[0],result->freqx[0],result->freqy[0],result->nobs[0],offset,F0,dFdx,dFdy,dFdz,dFdA,dFdoff);
// AOstruct *AO=result;
//  int nAOtotal=2*(AO->ntotal);
//     int nAOcols=3*nvert+3+2*AO->nao;
//     int nAOoffsets=2*(AO->nao);
//     double *offset=calloc(2*(AO->nao),sizeof(double));
//     
//     double *AOout=calloc(nAOtotal,sizeof(double));
//     double *AOdv=calloc(nAOtotal*nAOcols,sizeof(double));
//     Calculate_AOs(tlistn,vlistn,nfacn,nvertn,angles,AO,offset,D,nvertn,nvert,NULL,AOout,AOdv,1);
//     
// }
*/
