#include<math.h>
int line_intersect(double *a,double *b,double *x,double *y,double *p,double *dpx,double *dpy)
{
    /*
     * Check whether line determined by points a,b
     *intersects the line segment determined by x,y
     *If intersection, inters=true and the intersection is x+t(y-x)
     *p=(px,py)
     *dpx=[dpx/dx0,dpx/dy0,dpx/dx1,dpx/dy1] IF x=(x0,y0), y=(x1,y1). NOTE THAT
     *HERE x=(x0,x1) y=(y0,y1). SO dpx is derivative wrt 1.point x-coord,
     *1.point y-coord, 2.point x-coord 2.point y-coord
     *dtdx=[dt/dx0 dtdx1], dtdy=[dt/dy0,dtdy1];Checks whether line determined by point a,b 
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
