#include<math.h>
#include<stdio.h>
#define PI 3.141592654
double sinc(double x)
{
double limit=sqrt(sqrt(2.22e-16));
if(fabs(x*PI)>limit)
    return sin(x*PI)/(PI*x);
else
    return 1-pow(x*PI,2)/6.0+pow(x*PI,4)/120.0;
}

int main()
{
    for(int i=-5;i<5;i++)
    {
        printf("%f\n",sinc(i/100.0));
    }
}
