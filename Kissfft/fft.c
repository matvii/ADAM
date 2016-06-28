#include"/tmp/kiss_fft130/tools/kiss_fftr.h"
#include<stdio.h>
int main()
{
  //complex out[8];
  int size=8;
  kiss_fft_cpx fft_out[size];
  kiss_fftr_cfg cfg=kiss_fftr_alloc(size,0,0,0);
  double in[8]={0.1,-1.1,2.0,-4.1,3.1,2.0,-4.4,1.0};
  //kiss_fftr_cfg cfg = kiss_fftr_alloc( 8,0 ,0,0 );
kiss_fftr( cfg,in,fft_out);
	 printf("\n");
	 for(int j=0;j<8;j++)
	   printf("%f I%f\n",fft_out[j].r,fft_out[j].i);
	 printf("\n");
}



//gcc -std=c99 fft.c  kiss_fftr.c ../kiss_fft.c -lm -I/tmp/kiss_fft130/tools -I/tmp/kiss_fft130 -o fft
