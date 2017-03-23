#include<string.h>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"utils.h"


int linesegment_intersect_noderiv(double *a,double *b,double *x,double *y,double *p)
{
    /*
     * Check whether line determined by points a,b
     *intersects the line segment determined by x,y
     *If intersection, inters=true and the intersection is x+t(y-x)
     *p=(px,py)
     *dpx=[dpx/dx0,dpx/dx1,dpx/dy0,dpx/dy1] IF x=(x0,y0), y=(x1,y1). NOTE THAT
     *HERE x=(x0,x1) y=(y0,y1). SO dpx is derivative wrt 1.point x-coord,
     *1.point y-coord, 2.point x-coord 2.point y-coord
     *dtdx=[dt/dx0 dt/dy0], dtdy=[dt/dx1,dt/dy1];Checks whether line determined by point a,b 
     */
    double eps=1E-10;
    int inters=0;
    double t=0;
    double b0=b[0];
    double b1=b[1];
    double a0=a[0];
    double a1=a[1];
    double x0=x[0];
    double x1=x[1];
    double y0=y[0];
    double y1=y[1];
    double dtdx[2],dtdy[2];
    double denom=(b0-a0)*(x1-y1)-(b1-a1)*(x0-y0);
    
    double s;
    if(fabs(denom)<eps)
        return 0;
    t=((a0-x0)*(b1-a1)-(a1-x1)*(b0-a0))/denom;
    if(t<0 || t>1)
        return 0;
    
    s=-((x0-a0)*(y1-x1)-(x1-a1)*(y0-x0))/denom;
    if(s<0 || s>1)
        return 0;
   p[0]=x[0]+t*(y[0]-x[0]);
    p[1]=x[1]+t*(y[1]-x[1]);
   
    return 1;
}



void find_edge_noderiv(int *tlist,double *vlist2,int nfac,int nvert,int *visible,double *a,double *b,int *cledge,double *clpoint,int *inters)
{
    /*
     * *Find outer edge that line a->b intersects
    *We assume the shape is already rotated and projected to plane
    *determined by the z axis.
    *Adj is the nvertxnvert adjavency matrix, where Adj(i,j)!=0 if vertices i
    *and j are connected by an edge
    *OUTPUT
    *cledge closest edge to a, corresponds to vertices cledge[0] and cledge[1]
    * clt is t value corresponding to the intersection point,
    *clpoint=vlist(cledge[0],:)+clt*(vlist(cledge[1],:)-vlist(cledge[0],:))
    */
    double eps=1E-10;
    int *A=calloc(nvert*nvert,sizeof(int));
    int i1,i2,i3;
    double u1,u2,w1,w2,n3;
    double *v1,*v2,*v3;
    double interp[2],ip[2];
    
    double dist,cldist=1E9;
    double *vlist=calloc(3*nvert,sizeof(double));
    memcpy(vlist,vlist2,sizeof(double)*nvert*3);
    
    double bbx[2],bby[2];
    
    bbx[0]=minv(vlist,nvert,0);
    bbx[1]=maxv(vlist,nvert,0);
    bby[0]=minv(vlist,nvert,1);
    bby[1]=maxv(vlist,nvert,1);
    int t=1;
    double p[2];
    double b11,b12,b21,b22;
    double bv11,bv12,bv21,bv22;
    p[0]=a[0]+t*(b[0]-a[0]);
    p[1]=a[1]+t*(b[1]-a[1]);
    while(p[0]<bbx[1] && p[0]>bbx[0] && p[1]<bby[1] && p[1]>bby[0])
    {
        t++;
        p[0]=a[0]+t*(b[0]-a[0]);
        p[1]=a[1]+t*(b[1]-a[1]);
    }
   
    b11=fmin(0,p[0]);
    b12=fmin(0,p[1]);
    b21=fmax(0,p[0]);
    b22=fmax(0,p[1]);
    int sinters=0;
      for(int j=0;j<nfac;j++)
    {
        i1=tlist[3*j]-1; 
        i2=tlist[3*j+1]-1;
        i3=tlist[3*j+2]-1;
        
        if(visible[j]==0)
            continue;
        v1=vlist+3*i1;
        v2=vlist+3*i2;
        v3=vlist+3*i3;
        
        bv11=fmin(v1[0],fmin(v2[0],v3[0]));
        bv12=fmin(v1[1],fmin(v2[1],v3[1]));
        bv21=fmax(v1[0],fmax(v2[0],v3[0]));
        bv22=fmax(v1[1],fmax(v2[1],v3[1]));
        if(b12>bv22 || b22<bv12)
            continue;
        if(b21<bv11 || bv21<b11)
            continue;
            if(linesegment_intersect_noderiv(a,p,v1,v2,ip))
            {
                sinters=1;
                 //ip  is the intersection point of line a->b and the current edge
                dist=sqrt(pow(ip[0]-p[0],2)+pow(ip[1]-p[1],2));
                if(dist<cldist)
                {
                    cldist=dist;
                    cledge[0]=i1;
                    cledge[1]=i2;
                    clpoint[0]=ip[0];
                    clpoint[1]=ip[1];
                   
                }
               
            }
          
      
       
            if(linesegment_intersect_noderiv(a,p,v2,v3,ip))
            {
                
                 sinters=1;
                dist=sqrt(pow(ip[0]-p[0],2)+pow(ip[1]-p[1],2));
                 if(dist<cldist)
                {
                    cldist=dist;
                    cledge[0]=i2;
                    cledge[1]=i3;
                    clpoint[0]=ip[0];
                    clpoint[1]=ip[1];
                   
                }
                
            }
          
        
            if(linesegment_intersect_noderiv(a,p,v3,v1,ip))
            {
                
                 sinters=1;
                dist=sqrt(pow(ip[0]-p[0],2)+pow(ip[1]-p[1],2));
                 if(dist<cldist)
                {
                    cldist=dist;
                    cledge[0]=i3;
                    cledge[1]=i1;
                    clpoint[0]=ip[0];
                    clpoint[1]=ip[1];
                   
                }
               
            }
            
            
    }
   
    free(A);
    free(vlist);
    *inters=sinters;
}
            
            
