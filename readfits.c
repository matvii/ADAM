#include"wcstools-3.9.2/libwcs/fitsfile.h"
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
void readfits(char* filename,double **buffer,int x0,int y0,int nx,int ny,double *date,int *xsize,int *ysize)
{
    char *head;
    char *image;
    double *buf;
    int imsize=0;
    int nlog;
    int naxis1,naxis2;
    int bitpix;
    float *fimage;
    double *dimage;
    char keyword[]="MJD-OBS";
    int lhead,nbhead;
    head=fitsrhead(filename,&lhead,&nbhead);
    if(!hgetr8(head,keyword,date))
        *date=NAN;
    hgeti4(head,"NAXIS1",&naxis1);
    hgeti4(head,"NAXIS2",&naxis2);
    hgeti4(head,"BITPIX",&bitpix);
    if ((x0<=0 || y0<=0 || nx<=0 || ny<=0 || x0+nx>naxis1 || y0+ny>naxis2) && (naxis1%2!=0 || naxis2%2!=0))
   {
       x0=1;
       y0=1;
       nx=naxis1-naxis1%2;
       ny=naxis2-naxis2%2;
       printf("Using the whole image %s, but size is odd. Fixing\n",filename);
   }
    if(x0<=0 || y0<=0 ||nx<=0 || ny<=0 || x0+nx>naxis1 || y0+ny>naxis2)
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
         buf[j]=fmax(fimage[j],0);
    else
        for(int j=0;j<imsize;j++)
            buf[j]=fmax(dimage[j],0);
    free(head);
    free(image);
    *buffer=buf;
    //printf("xsize: %d ysize: %d\n",*xsize,*ysize);
    //write_matrix_file("/tmp/test.txt",*buffer,*ysize,*xsize);
  //exit(1);
}
 /*       
void main()
{
    char *head,test;
    char *image;
    float *fimage;
    int xsize,ysize;
    int lhead,nbhead,ival;
    char filename[]="test.fits";
    char keyword[]="MJD-OBS";
    int imsize;
    double date;
    double *buffer;
    double *buffer2;
    //head=fitsrhead(filename,&lhead,&nbhead);
    //printf("ldhead: %d nbhead: %d\n",lhead,nbhead);
    //printf("%s\n",head);
    //if(hgetr8(head,keyword,&date))
    //    printf("DATE: %f\n",date);
 // image=fitsrimage(filename, nbhead, head);
//  fimage=(float *)image;
  
//  printf("\n %f \n",fimage[256*255+255]);
    readfits(filename,&buffer,0,0,0,0,&date,&xsize,&ysize);
    printf("xsize: %d ysize: %d\n",xsize,ysize);
    printf("date: %f\n",date);
    for(int j=0;j<10;j++)
     printf(" %f ",buffer[j+2*256]);
    
    date=0;
    xsize=0;
    ysize=0;
    readfits(filename,&buffer2,0,0,10,10,&date,&xsize,&ysize);
    printf("xsize: %d ysize: %d\n",xsize,ysize);
    printf("date: %f\n",date);
    for(int j=0;j<10;j++)
     printf(" %f ",buffer2[j+2*10]);
    
}
*/
