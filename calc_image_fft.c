
#include"Kissfft/tools/kiss_fftndr.h"
//#include"Kissfft/kiss_fft.h"
//#include"Kissfft/tools/kiss_fftr.h"
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<complex.h>
#include"utils.h"
void calc_image_fft(double *M,int m,int n,double dx,double dy,double *zMr,double *zMi,double *Fx,double *Fy)
{
    /*M image matrix mxn, m and n EVEN. Horizontal is x direction, vertical y
     * Zero coordinate at (m/2+1,n/2+1)
     * Lx=n*dx, Ly=m*dy
     * Returns 2d fft transform of M,m x (n/2) matrix, with reals in in zMr and imaginary part in zMi 
     * corresponding frequencies in x-direction
     * 0,1/Lx,2/Lx,...(n/2-1)*1/Lx
     * and in y direction
     * 0,1/Ly,..,(m/2-1)*1/Ly,-m/2*Ly,...,-1/Ly
     * So first n/2 entries in zMr+i*zMi correspond to spatial frequencies fx=0,1/Lx,...(n/2-1)*1/Lx and fy=0
     * next n/2 entries correspond to fx=0,1/Lx,...(n/2-1)*1/Lx and fy=1/Ly
     * NOTE: We normalize with DC term and then drop it*/
    /*Output: zMr,zMi m*(n/2)-1 arrays, Fx,Fy m*(n/2)-1 arrays. The kth term of zMr+zMi corresponds to frequency Fx(k),Fy(k). */
    /*TBD: PSF!!*/  
    const int dim[2]={m,n};
    const int dimcount=2;
    double dc;
    double Lx,Ly;
    Lx=n*dx;
    Ly=m*dy;
    kiss_fft_cpx *fft2_out;
    fft2_out=malloc(dim[0]*dim[1]*sizeof(kiss_fft_cpx));
    kiss_fftndr_cfg stf = kiss_fftndr_alloc(dim, dimcount, 0, 0, 0);
    kiss_fftndr(stf, M, fft2_out);
    /*fft2_out is mx(n/2+1) complex matrix matrix*/
    dc=fft2_out[0].r; /*dc term*/
    /*Origo is assumed to at the center of image, so we have to compensate for phase shift
     * caused by fft*/
    for(int i=0;i<m;i++)
    {
        for(int j=0;j<n/2;j++)
        {
        
            
            zMr[i*n/2+j]=pow(-1.0,i+j)*fft2_out[i*(n/2+1)+j].r/dc;
            zMi[i*n/2+j]=pow(-1.0,i+j)*fft2_out[i*(n/2+1)+j].i/dc;
            Fx[i*n/2+j]=j*1/Lx;
            Fy[i*n/2+j]=i*1/Ly;
            if(i>m/2-1)
                Fy[i*n/2+j]=(i-m)*1/Ly;
        }
    }
    free(fft2_out);
    free(stf);
}
 /*           
   int main()
{
    double a[]={0.16218,0.79428,0.31122,0.52853,0.16565,0.60198,0.26297,0.65408,0.68921,0.74815,0.45054,0.083821,0.22898,0.91334,0.15238,0.82582,0.53834,0.99613,0.078176,0.44268,0.10665,0.9619,0.0046342,0.77491,0.8173,0.86869,0.084436,0.39978,0.25987,0.80007,0.43141,0.91065,0.18185,0.2638,0.14554,0.13607,0.86929,0.5797,0.54986,0.14495,0.85303,0.62206,0.35095,0.51325,0.40181,0.075967,0.23992,0.12332,0.18391,0.23995,0.41727,0.049654,0.90272,0.94479,0.49086,0.48925,0.33772,0.90005,0.36925,0.1112,0.78025,0.38974,0.24169,0.40391,0.096455,0.13197,0.94205,0.95613,0.57521,0.05978,0.23478,0.35316,0.82119,0.015403,0.043024,0.16899,0.64912,0.73172,0.64775,0.45092};
    double *Mr,*Mi,*Fx,*Fy;
    Mr=calloc(8*5,sizeof(double));
    Mi=calloc(8*5,sizeof(double));
    Fx=calloc(8*5,sizeof(double));
    Fy=calloc(8*5,sizeof(double));
    calc_image_fft(a,8,10,0.5,0.5,Mr,Mi,Fx,Fy);
    //Zeros frequency here is rubbish since it is discarded
    print_matrix(Fx,8,5);
    print_matrix(Fy,8,5);
    print_matrix(Mr,8,5);
    print_matrix(Mi,8,5);
    write_matrix_file("/tmp/Fx2.txt",Fx,1,40-1);
    write_matrix_file("/tmp/Fy2.txt",Fy,1,40-1);
    write_matrix_file("/tmp/Mr2.txt",Mr,1,40-1);
    write_matrix_file("/tmp/Mi2.txt",Mi,1,40-1);
//     free(Mr);
//     free(Mi);
//     free(Fx);
//     free(Fy);
//     char  filename[]="Hebe/jittered_result0.fits";
//     double *buffer;
//     int xsize;
//     int ysize;
//     double date=NAN;
//     readfits(filename,&buffer,79,79,100,100,&date,&xsize,&ysize);
//     printf("xsize: %d ysize: %d First:%f Last:%f\n",xsize,ysize,buffer[0],buffer[xsize*ysize-1]);
//     //print_matrix(buffer,ysize,xsize);
//         double *Mr,*Mi,*Fx,*Fy;
//         int imsize=xsize*ysize;
//      Mr=calloc(imsize/2,sizeof(double));
//      Mi=calloc(imsize/2,sizeof(double));
//      Fx=calloc(imsize/2,sizeof(double));
//      Fy=calloc(imsize/2,sizeof(double));
//     calc_image_fft(buffer,xsize,ysize,0.009942,0.009942,Mr,Mi,Fx,Fy);
//     write_matrix_file("/tmp/buffer.txt",buffer,ysize,xsize);
//     write_matrix_file("/tmp/Mr.txt",Mr,1,imsize/2);
//     write_matrix_file("/tmp/Mi.txt",Mi,1,imsize/2);
//     write_matrix_file("/tmp/Fx.txt",Fx,1,imsize/2);
//     write_matrix_file("/tmp/Fy.txt",Fx,1,imsize/2);
}

*/
