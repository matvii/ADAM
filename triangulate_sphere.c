#include<stdio.h>
#include<math.h>
#include"utils.h"
#include<stdlib.h>
void triangulate_sphere(int nrows,double *t,double *f,int *ifp)
{
    /*
     * t and f are 4*pow(nrows,2)+2 arrays
     * ifp is 8*pow(nrows,2) x 3 matrix
     */
    double dth=PI/(2*nrows);
    double dph;
    int k=0;
    t[0]=0;
    f[0]=0;
    for(int i=1;i<=nrows;i++)
    {
        dph=PI/(2*i);
        for(int j=0;j<=(4*i-1);j++)
        {
            k++;
            t[k]=i*dth;
            f[k]=j*dph;
        }
    }
    for(int i=nrows-1;i>=1;i--)
    {
        dph=PI/(2*i);
        for(int j=0;j<=(4*i-1);j++)
        {
            k++;
            t[k]=PI-i*dth;
            f[k]=j*dph;
        }
    }
    t[k+1]=PI;
    f[k+1]=0;
    
    int *nod=calloc((2*nrows+1)*(4*nrows+1),sizeof(int));
    int nnod=1;
    nod[0]=1;
    
    for(int i=1;i<=nrows;i++)
    {
        for(int j=0;j<=(4*i-1);j++)
        {
            nnod++;
            set_elI(nod,2*nrows+1,4*nrows+1,nnod,i,j);
            if(j==0)
                set_elI(nod,2*nrows+1,4*nrows+1,nnod,i,4*i);
        }
    }
    
    
    for(int i=nrows-1;i>=1;i--)
    {
        for(int j=0;j<=4*i-1;j++)
        {
            nnod++;
            set_elI(nod,2*nrows+1,4*nrows+1,nnod,2*nrows-i,j);
            if(j==0)
                set_elI(nod,2*nrows+1,4*nrows+1,nnod,2*nrows-i,4*i);
        }
    }
     
    set_elI(nod,2*nrows+1,4*nrows+1,nnod+1,2*nrows,0);
    
    int ntri=0,j0,el;
   
    for(int j1=1;j1<=nrows;j1++)
    {
        for(int j3=1;j3<=4;j3++)
        {
            j0=(j3-1)*j1;
            ntri++;
            el=get_elI(nod,2*nrows+1,4*nrows+1,j1-1,j0-(j3-1));
            set_elI(ifp,8*pow(nrows,2),3,el,ntri-1,0);
            
            el=get_elI(nod,2*nrows+1,4*nrows+1,j1,j0);
            set_elI(ifp,8*pow(nrows,2),3,el,ntri-1,1);
            
            el=get_elI(nod,2*nrows+1,4*nrows+1,j1,j0+1);
            set_elI(ifp,8*pow(nrows,2),3,el,ntri-1,2);
            for(int j2=j0+1;j2<=j0+j1-1;j2++)
            {
                ntri++;
                el=get_elI(nod,2*nrows+1,4*nrows+1,j1,j2);
                set_elI(ifp,8*pow(nrows,2),3,el,ntri-1,0);
            
                el=get_elI(nod,2*nrows+1,4*nrows+1,j1-1,j2-(j3-1));
                set_elI(ifp,8*pow(nrows,2),3,el,ntri-1,1);
            
                el=get_elI(nod,2*nrows+1,4*nrows+1,j1-1,j2-(j3-1)-1);
                set_elI(ifp,8*pow(nrows,2),3,el,ntri-1,2);
                
                ntri++;
                
                el=get_elI(nod,2*nrows+1,4*nrows+1,j1-1,j2-(j3-1));
                set_elI(ifp,8*pow(nrows,2),3,el,ntri-1,0);
            
                el=get_elI(nod,2*nrows+1,4*nrows+1,j1,j2);
                set_elI(ifp,8*pow(nrows,2),3,el,ntri-1,1);
            
                el=get_elI(nod,2*nrows+1,4*nrows+1,j1,j2+1);
                set_elI(ifp,8*pow(nrows,2),3,el,ntri-1,2);
            }
        }
    }
    
    
 /*Lower hemisphere*/
    for(int j1=nrows+1;j1<=2*nrows;j1++)
    {
        for(int j3=1;j3<=4;j3++)
        {
            j0=(j3-1)*(2*nrows-j1);
            ntri++;
            el=get_elI(nod,2*nrows+1,4*nrows+1,j1,j0);
            set_elI(ifp,8*pow(nrows,2),3,el,ntri-1,0);
            el=get_elI(nod,2*nrows+1,4*nrows+1,j1-1,j0+(j3-1)+1);
            set_elI(ifp,8*pow(nrows,2),3,el,ntri-1,1);
            el=get_elI(nod,2*nrows+1,4*nrows+1,j1-1,j0+(j3-1));
            set_elI(ifp,8*pow(nrows,2),3,el,ntri-1,2);
            
            for(int j2=j0+1;j2<=j0+(2*nrows-j1);j2++)
            {
                
                ntri++;
                el=get_elI(nod,2*nrows+1,4*nrows+1,j1,j2);
                set_elI(ifp,8*pow(nrows,2),3,el,ntri-1,0);
                el=get_elI(nod,2*nrows+1,4*nrows+1,j1-1,j2+(j3-1));
                set_elI(ifp,8*pow(nrows,2),3,el,ntri-1,1);
                el=get_elI(nod,2*nrows+1,4*nrows+1,j1,j2-1);
                set_elI(ifp,8*pow(nrows,2),3,el,ntri-1,2);
                
                ntri++;
                el=get_elI(nod,2*nrows+1,4*nrows+1,j1,j2);
               
                set_elI(ifp,8*pow(nrows,2),3,el,ntri-1,0);
            
                el=get_elI(nod,2*nrows+1,4*nrows+1,j1-1,j2+1+(j3-1));
                set_elI(ifp,8*pow(nrows,2),3,el,ntri-1,1);
            
                el=get_elI(nod,2*nrows+1,4*nrows+1,j1-1,j2+(j3-1));
                set_elI(ifp,8*pow(nrows,2),3,el,ntri-1,2);
            }
        }
    }
  
}
/*   
void main()
{
    int *ifp;
    int nrows=2;
    double *t,*f;
    ifp=calloc(3*8*pow(nrows,2),sizeof(int));
    t=calloc(4*pow(nrows,2)+2,sizeof(double));
    f=calloc(4*pow(nrows,2)+2,sizeof(double));
     triangulate_sphere(nrows,t,f,ifp);
     print_matrix(t,1,4*pow(nrows,2)+2);
     print_matrix(f,1,4*pow(nrows,2)+2);
     print_matrixI(ifp,8*pow(nrows,2),3);
}
 

*/

        
