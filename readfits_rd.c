#include"wcstools-3.9.2/libwcs/fitsfile.h"
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"utils.h"
#include"matrix_ops.h"
#include"globals.h"
void readfits_rd(char* filename,double **buffer,int x0,int y0,int nx,int ny,double *date,int *xsize,int *ysize,double *cdelt1,double *cdelt2)
{
    /*Routine to read range-Doppler radar fits and information specific to it*/
    /*NOTE: xsize is HORIZONTAL axis,ysize is VERTICAL, so resulting image is ysize x xsize matrix*/
    /*NOTE 2: In fact xsize is naxis1, so it frequency axis.*/
    char *head;
    char *image;
    double *buf;
    int imsize=0;
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
    hgeti4(head,"NAXIS1",&naxis1);
    hgeti4(head,"NAXIS2",&naxis2);
    hgeti4(head,"BITPIX",&bitpix);
    if(nx==0 || ny==0 || nx>naxis1 || ny>naxis2)
    {
    *xsize=naxis1;
    *ysize=naxis2;
    imsize=naxis1*naxis2;
    image=fitsrimage(filename, nbhead, head);
    if(bitpix!=-64)
        fimage=(float *)image;
    else
        dimage=(double *)image;
    
    }
    else
    {
    *xsize=nx;
    *ysize=ny;
    imsize=nx*ny;
    image=fitsrsect(filename,head,nbhead,x0,y0,nx,ny,nlog);
    //printf("x0:%d y0: %d\n",x0,y0);
    //printf("Input readfits: %d %d %d %d\n",x0,y0,nx,ny);
    if(bitpix!=-64)
        fimage=(float *)image;
    else
        dimage=(double *)image;
        
    }
    buf=calloc(imsize,sizeof(double));
      if(bitpix!=-64)
      for(int j=0;j<imsize;j++)
         //buf[j]=fmax(fimage[j],0);
          buf[j]=fmax(fimage[j],INI_SET_RD_ZERO);
    else
        for(int j=0;j<imsize;j++)
            buf[j]=fmax(dimage[j],INI_SET_RD_ZERO);
      //NOTE: HERE WE flip and transpose the image
     // flip_dim(buf,*ysize,*xsize,2);
    //  matrix_transpose(buf,*ysize,*xsize);
     
    free(head);
    free(image);
    *buffer=buf;
  
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
