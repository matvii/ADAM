#include"wcstools-3.9.2/libwcs/fitsfile.h"
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"utils.h"
#include"matrix_ops.h"
#include"globals.h"
void readfits_rd(char* filename,double **buffer,int cx,int cy,int cxdim,int cydim,double *date,int *xsize,int *ysize,double *cdelt1,double *cdelt2)
{
    /*Routine to read range-Doppler radar fits and information specific to it*/
    /*NOTE: xsize is HORIZONTAL axis,ysize is VERTICAL, so resulting image is ysize x xsize matrix*/
    /*NOTE 2: In fact xsize is naxis1, so it frequency axis.*/
    char *head;
    char *image;
    double *buf;
    double crpix1,crpix2;
    int imsize=0;
    int nx=0,ny=0,x0=0,y0=0;
    int nlog;
    int naxis1,naxis2;
    float *fimage;
    double *dimage;
    int bitpix=0;
    char keyword[]="JDMEAN";
    int lhead,nbhead;
    head=fitsrhead(filename,&lhead,&nbhead);
    if(!hgetr8(head,keyword,date))
        *date=NAN;
    if(!hgetr8(head,"CDELT1",cdelt1))
        *cdelt1=0;
    if(!hgetr8(head,"CDELT2",cdelt2))
        *cdelt2=0;
    if(!hgetr8(head,"CRPIX1",&crpix1))
        crpix1=0;
    if(!hgetr8(head,"CRPIX2",&crpix2))
        crpix2=0;
    if(crpix1>0 && crpix2>0)
        printf("Image %s, Doppler res %f, Delay res %f, Doppler center %f, Delay center %f\n",filename,*cdelt1,*cdelt2,crpix1,crpix2);
    
    if(cx>0 && cy>0 && crpix1>0)
        printf("RDcenter set, overriding CRPIX from the %s\n",filename); 
    
    if(cx>0 && cy>0 && cxdim>0 &&cydim>0)
    {
       
        nx=2*cxdim;
        ny=2*cydim;
        x0=cx-cxdim;
        y0=cy-cydim;
    }
    else if(cxdim>0 && cydim>0)
    {
       
        cx=round(crpix1);
        cy=round(crpix2);
        nx=2*cxdim;
        ny=2*cydim;
        x0=cx-cxdim;
        y0=cy-cydim;
    }
    
    
   
    hgeti4(head,"NAXIS1",&naxis1);
    hgeti4(head,"NAXIS2",&naxis2);
    hgeti4(head,"BITPIX",&bitpix);
// printf("file: %s cx %d cy %d nx %d ny %d x0: %d y0: %d cxdim: %d  cydim: %d naxis1 %d naxis2 %d\n",filename,cx,cy,nx,ny,x0,y0,cxdim,cydim,naxis1,naxis2);
    if(x0+nx>naxis1)
    {
        printf("in image %s, set image size in ini is larger than physical dimension, using full image in x direction\n",filename);
        x0=1;
        nx=naxis1-naxis1%2;
    }
    if(y0+ny>naxis2)
    {
        printf("in image %s, set image size in ini is larger than physical dimension, using full image in y direction\n",filename);
        y0=1;
        ny=naxis2-naxis2%2;
    }
    
   if ((x0==0 || y0==0 || nx==0 || ny==0 || x0+nx>naxis1 || y0+ny>naxis2))
   {
       x0=1;
       y0=1;
       nx=naxis1-naxis1%2;
       ny=naxis2-naxis2%2;
  //  printf("x0 %d y0 %d nx: %d ny: %d %d %d\n",x0,y0,nx,ny,x0+nx,y0+ny);
   //    printf("Using the whole image %s, but size is odd. Fixing\n",filename);
   }
   
    
   
    *xsize=nx;
    *ysize=ny;
    imsize=nx*ny;
    image=fitsrsect(filename,head,nbhead,x0,y0,nx,ny,nlog);
//    printf("x0:%d y0: %d\n",x0,y0);
  //  printf("Input readfits: %d %d %d %d\n",x0,y0,nx,ny);
    if(bitpix!=-64)
        fimage=(float *)image;
    else
        dimage=(double *)image;
        
    
    buf=calloc(imsize,sizeof(double));
      if(bitpix!=-64)
      for(int j=0;j<imsize;j++)
         //buf[j]=fmax(fimage[j],0);
          buf[j]=fmax(fimage[j],INI_SET_RD_ZERO);
    else
        for(int j=0;j<imsize;j++)
            buf[j]=fmax(dimage[j],INI_SET_RD_ZERO);
      //NOTE: HERE WE flip and transpose the image
       if(INI_RADAR_FLIP==1) 
     flip_dim(buf,*ysize,*xsize,2);
    //  matrix_transpose(buf,*ysize,*xsize);
     
    free(head);
    free(image);
    *buffer=buf;
//    printf("xsize: %d ysize: %d\n",*xsize,*ysize);
  // write_matrix_file("/tmp/test.txt",*buffer,*ysize,*xsize);
  //exit(1);
}
/*
int main()
{
    char filename[]="/tmp/test.fits";
    double *buffer;
    double date;
    int xsize,ysize;
    double cdelt1=0,cdelt2=0;
    readfits_rd(filename,&buffer,0,0,0,0,&date,&xsize,&ysize,&cdelt1,&cdelt2);
    printf("date %f xsize: %d ysize: %d cdelt1: %f cdelt2: %f\n",date,xsize,ysize,cdelt1,cdelt2);
    write_matrix_file("/tmp/B.txt",buffer,xsize,ysize);
    
    //Then transpose
}
*/
