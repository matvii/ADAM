#include<string.h>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"utils.h"

int line_intersect(double *a,double *b,double *x,double *y,double *p,double *dpx,double *dpy)
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
    if(fabs(denom)<eps)
        return 0;
    t=((a0-x0)*(b1-a1)-(a1-x1)*(b0-a0))/denom;
    if(t>=0 && t<=1)
    {
        inters=1;
        p[0]=x[0]+t*(y[0]-x[0]);
        p[1]=x[1]+t*(y[1]-x[1]);
    }
    dtdx[0]=(a1-b1)*(a1*(b0-y0)+b1*y0-b0*y1+a0*(-b1+y1))/pow(denom,2);
    dtdx[1]=(a1 - b1)*(-b1*x0 + a1*(-b0 + x0) + a0*(b1 - x1) + b0*x1)/pow(denom,2);
    dtdy[0]=(a0 - b0)*(-b1*y0 + a1*(-b0 + y0) + a0*(b1 - y1) + b0*y1)/pow(denom,2);
    dtdy[1]=-(a0 - b0)*(-b1*x0 + a1*(-b0 + x0) + a0*(b1 - x1) + b0*x1)/pow(denom,2);
    dpx[0]=1+dtdx[0]*(y0-x0)-t;
    dpx[1]=dtdy[0]*(y0-x0);
    dpx[2]=dtdx[1]*(y0-x0)+t;
    dpx[3]=dtdy[1]*(y0-x0);
    
    dpy[0]=dtdx[0]*(y1-x1);
    dpy[1]=1+dtdy[0]*(y1-x1)-t;
    dpy[2]=dtdx[1]*(y1-x1);
    dpy[3]=dtdy[1]*(y1-x1)+t;
    return inters;
}



void find_chord(int *tlist,double *vlist2,int nfac,int nvert,double *offset,double *a,double *b,int *Adj,int *cledge,int *faedge,double *clpoint,double *fapoint,int *inters,double *dclpx,double *dclpy,double *dfapx,double *dfapy)
{
    /*
     * *Find two outer edges that line a->b intersects
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
    double ip[2];
    double dpx[4],dpy[4];
    double dist,cldist=1E9,fadist=0;
    double *vlist=calloc(3*nvert,sizeof(double));
    memcpy(vlist,vlist2,sizeof(double)*nvert*3);
    memcpy(A,Adj,sizeof(int)*nvert*nvert);
    double bbx[2],bby[2];
    double offx=offset[0];
    double offy=offset[1];
    //Add offsets to vlist
    for(int j=0;j<nvert;j++)
    {
        vlist[3*j]+=offx;
        vlist[3*j+1]+=offy;
    }
    bbx[0]=minv(vlist,nvert,0);
    bbx[1]=maxv(vlist,nvert,0);
    bby[0]=minv(vlist,nvert,1);
    bby[1]=maxv(vlist,nvert,1);
    int t=0;
    double p[2];
    p[0]=a[0]+t*(b[0]-a[0]);
    p[1]=a[1]+t*(b[1]-a[1]);
    while(p[0]<bbx[1] && p[0]>bbx[0] && p[1]<bby[1] && p[1]>bby[0])
    {
        t--;
        p[0]=a[0]+t*(b[0]-a[0]);
        p[1]=a[1]+t*(b[1]-a[1]);
    }
    // p is a point on line a-b outside the shape. We use p to find the closest
    //and farthest intersection points of model and and chord a->b. 
    clpoint[0]=0;
    clpoint[1]=0;
    fapoint[0]=0;
    fapoint[1]=0;
    int sinters=0;
      for(int j=0;j<nfac;j++)
    {
        i1=tlist[3*j]-1; 
        i2=tlist[3*j+1]-1;
        i3=tlist[3*j+2]-1;
        u1=vlist[3*i2]-vlist[3*i1];
        u2=vlist[3*i2+1]-vlist[3*i1+1];
        w1=vlist[3*i3]-vlist[3*i1];
        w2=vlist[3*i3+1]-vlist[3*i1+1];
        n3=u1*w2-u2*w1;
        if(n3<=0)
            continue;
        v1=vlist+3*i1;
        v2=vlist+3*i2;
        v3=vlist+3*i3;
        
        if(get_elI(A,nvert,nvert,i1,i2))
        {
           
            if(line_intersect(a,b,v1,v2,ip,dpx,dpy))
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
                    memcpy(dclpx,dpx,sizeof(double)*4);
                    memcpy(dclpy,dpy,sizeof(double)*4);
                }
                if(dist>fadist)
                {
                    fadist=dist;
                    faedge[0]=i1;
                    faedge[1]=i2;
                    fapoint[0]=ip[0];
                    fapoint[1]=ip[1];
                    
                    memcpy(dfapx,dpx,sizeof(double)*4);
                    memcpy(dfapy,dpy,sizeof(double)*4);
                }
            }
            set_elI(A,nvert,nvert,0,i2,i1);
            set_elI(A,nvert,nvert,0,i1,i2);
        }
      
        if(get_elI(A,nvert,nvert,i2,i3))
        {
            if(line_intersect(a,b,v2,v3,ip,dpx,dpy))
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
                    memcpy(dclpx,dpx,sizeof(double)*4);
                    memcpy(dclpy,dpy,sizeof(double)*4);
                }
                if(dist>fadist)
                {
                    fadist=dist;
                    faedge[0]=i2;
                    faedge[1]=i3;
                    fapoint[0]=ip[0];
                    fapoint[1]=ip[1];
                    
                    memcpy(dfapx,dpx,sizeof(double)*4);
                  
                    memcpy(dfapy,dpy,sizeof(double)*4);
                }
            }
            set_elI(A,nvert,nvert,0,i2,i3);
            set_elI(A,nvert,nvert,0,i3,i2);
        }
       
        if(get_elI(A,nvert,nvert,i3,i1))
        {
            if(line_intersect(a,b,v3,v1,ip,dpx,dpy))
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
                    memcpy(dclpx,dpx,sizeof(double)*4);
                    memcpy(dclpy,dpy,sizeof(double)*4);
                }
                if(dist>fadist)
                {
                    fadist=dist;
                    faedge[0]=i3;
                    faedge[1]=i1;
                    fapoint[0]=ip[0];
                    fapoint[1]=ip[1];
                    
                    memcpy(dfapx,dpx,sizeof(double)*4);
                  
                    memcpy(dfapy,dpy,sizeof(double)*4);
                }
            }
            set_elI(A,nvert,nvert,0,i3,i1);
            set_elI(A,nvert,nvert,0,i1,i3);
        }
    }
    if(pow(clpoint[0]-fapoint[0],2)+pow(clpoint[1]-fapoint[1],2)>eps && sinters==1)
        *inters=1;
    else
        *inters=0;
    free(A);
    free(vlist);
}
            
 
