#include"utils.h"
#include"matrix_ops.h"
#include"globals.h"

 void find_edge(int *tlist,double *vlist2,int nfac,int nvert,int *visible,double *offset,double *a,double *b,int *cledge,double *clpoint,int *inters,double *dclpx,double *dclpy,double *doffx,double *doffy);     
void Fit_Contour2(int *tlist,double *vlist,int nfac,int nvert,double *angles,double *up,double *E,double *E0,double TIME,double *offset,double *datax,double *datay,int npoints,double *dist,double *dx,double *dy,double *dz,double *dangles,double *dtox,double *dtoy)
{
    /*INPUT:
     * TIMES nchordsx2 disappearance and appearance times
     * offset=[x,y] offset in plane
     * chords nchordsx4 coordinates, disappearance and appearance
     * Chordoffset nchords vector,containing chord offset in seconds
     * OUTPUT:
     * dist nchordsx4 matrix, x and y distance ofclosest model points to disappearance and appearance points
     * dx,dy,dz 4nchordsxnvert matrix
     * dangles 4nchordsx3 matrix
     * dtox,dtoy 2*npointsx1, derivatives wrt offsets
     * dCOdoff 2*npointsxnchords matrix, derivatives for chord offsets
     * */
    
    double M[3][3],MRT[3][3],MRdb[3][3],MRdl[3][3],MRdo[3][3],MRTT[9];
    double Er[3],E0r[3];
    double MRdbT[9],MRdlT[9],MRdoT[9];
    double R[3][3],Rb[3][3],Rl[3][3],Ro[3][3],Rt[3][3],RbT[3][3],RlT[3][3],RoT[3][3];
    int cledge[2],inters;
    int edget[2];
    int ic1,ic2,fa1,fa2;
    double LARGE_VALUE=1e4;
    double clpoint[2],dclpx[4],dclpy[4];
    double temp4[4];
    double V[3];
    double off[]={0,0};
     double nd_dist,nd_dx,nd_dy;
     int nd_vert;
   
   // int *Adj=calloc(nvert*nvert,sizeof(int));
   // AdjFacet(tlist,vlist,nfac,nvert,Adj);
    Calculate_Frame_Matrix(E,up,M);
    rotate(angles[0],angles[1],angles[2],angles[3],TIME,R,Rb,Rl,Ro);
   //Here we determine which facets are visible
    //We rotate view/sun directions
    int *visible;
   visible=calloc(nfac,sizeof(int));
   mult_vector(R,E,Er);
    mult_vector(R,E0,E0r);
   FindActualBlockers(tlist,vlist,nfac,nvert,Er,E0r,1,visible);
    //visible[j]==1 if jth facet is visible
   
    transpose(R,Rt); //Transpose, since we rotate the model, not view directions
 transpose(Rb,RbT);
 transpose(Rl,RlT);
 transpose(Ro,RoT);
 mult_mat(M,Rt,MRT); 
double mr11,mr12,mr13,mr21,mr22,mr23;
mr11=MRT[0][0];
mr12=MRT[0][1];
mr13=MRT[0][2];
mr21=MRT[1][0];
mr22=MRT[1][1];
mr23=MRT[1][2];
 mult_mat(M,RbT,MRdb);
 mult_mat(M,RlT,MRdl);
 mult_mat(M,RoT,MRdo);
 zero_array(dist,npoints);
zero_array(dx,(npoints)*nvert);
zero_array(dy,(npoints)*nvert);
zero_array(dz,(npoints)*nvert);
zero_array(dangles,(npoints)*3);
zero_array(dtox,npoints);
zero_array(dtox,npoints);

transpose2(MRT,MRTT);
transpose2(MRdb,MRdbT);
transpose2(MRdl,MRdlT);
transpose2(MRdo,MRdoT);
double a[2],b[2];
 double *vlist2=calloc(nvert*3,sizeof(double));
 double *vlist2b=calloc(nvert*3,sizeof(double));
 double *vlist2l=calloc(nvert*3,sizeof(double));
 double *vlist2o=calloc(nvert*3,sizeof(double));
 matrix_prod(vlist,nvert,3,MRTT,3,vlist2); //vlist2 contains rotated and projected model
  matrix_prod(vlist,nvert,3,MRdbT,3,vlist2b);
  matrix_prod(vlist,nvert,3,MRdlT,3,vlist2l);
  matrix_prod(vlist,nvert,3,MRdoT,3,vlist2o);
  double w;
 double el=0;
 double d2=0;
 double distx;
 double disty;
 double distxdx1,distxdy1,distxdz1,distxdx2,distxdy2,distxdz2;
 double distydx1,distydy1,distydz1,distydx2,distydy2,distydz2;
 double distxdox,distxdoy,distydox,distydoy;
 double distxdb,distxdl,distxdo;
 double distydb,distydl,distydo;
 double doffx[2],doffy[2];
 a[0]=0;
 a[1]=0;
 double *dist2=calloc(nvert,sizeof(double));
 int *vertices=calloc(nvert,sizeof(int));
 for(int j=0;j<npoints;j++)
 {
     b[0]=datax[j];
     b[1]=datay[j];
     
     
    
     find_edge(tlist,vlist2,nfac,nvert,visible,offset,a,b,cledge,clpoint,&inters,dclpx,dclpy,doffx,doffy);
     
     if(inters==1)
     {
         
             
         distx=clpoint[0]-datax[j];
         disty=clpoint[1]-datay[j];
         distxdx1=dclpx[0]*mr11+dclpx[1]*mr21; //distx wrt ic1
         distxdy1=dclpx[0]*mr12+dclpx[1]*mr22;
         distxdz1=dclpx[0]*mr13+dclpx[1]*mr23;
         
         distxdx2=dclpx[2]*mr11+dclpx[3]*mr21;//distx wrt ic2
         distxdy2=dclpx[2]*mr12+dclpx[3]*mr22;
         distxdz2=dclpx[2]*mr13+dclpx[3]*mr23;
        
         
         distydx1=dclpy[0]*mr11+dclpy[1]*mr21; //disty wrt ic1
         distydy1=dclpy[0]*mr12+dclpy[1]*mr22;
         distydz1=dclpy[0]*mr13+dclpy[1]*mr23;
         
         distydx2=dclpy[2]*mr11+dclpy[3]*mr21;
         distydy2=dclpy[2]*mr12+dclpy[3]*mr22;
         distydz2=dclpy[2]*mr13+dclpy[3]*mr23;
        d2=pow(distx,2)+pow(disty,2);
         dist[j]=sqrt(d2); //distance in x
        
        
         ic1=cledge[0];
         ic2=cledge[1];
          vertices[ic1]=1;
         vertices[ic2]=1;
        set_el(dx,nvert+npoints,nvert,pow(d2,-0.5)*(distx*distxdx1+disty*distydx1),j,ic1);
        set_el(dx,nvert+npoints,nvert,pow(d2,-0.5)*(distx*distxdx2+disty*distydx2),j,ic2);
        
        set_el(dy,nvert+npoints,nvert,pow(d2,-0.5)*(distx*distxdy1+disty*distydy1),j,ic1);
        set_el(dy,nvert+npoints,nvert,pow(d2,-0.5)*(distx*distxdy2+disty*distydy2),j,ic2);
        
        set_el(dz,nvert+npoints,nvert,pow(d2,-0.5)*(distx*distxdz1+disty*distydz1),j,ic1);
        set_el(dz,nvert+npoints,nvert,pow(d2,-0.5)*(distx*distxdz2+disty*distydz2),j,ic2);
         /*
         el=dclpx[0]*mr11+dclpx[1]*mr21;
         
         set_el(dx,2*npoints,nvert,el,2*j,ic1);
          el=dclpx[0]*mr12+dclpx[1]*mr22;
         
         set_el(dy,2*npoints,nvert,el,2*j,ic1);
          el=dclpx[0]*mr13+dclpx[1]*mr23;
         
         set_el(dz,2*npoints,nvert,el,2*j,ic1);
         
          el=dclpx[2]*mr11+dclpx[3]*mr21;
        
         set_el(dx,2*npoints,nvert,el,2*j,ic2);
        el=dclpx[2]*mr12+dclpx[3]*mr22;
        
         set_el(dy,2*npoints,nvert,el,2*j,ic2);
          el=dclpx[2]*mr13+dclpx[3]*mr23;
          
         set_el(dz,2*npoints,nvert,el,2*j,ic2);
         
         el=dclpy[0]*mr11+dclpy[1]*mr21;
         
         set_el(dx,2*npoints,nvert,el,2*j+1,ic1);
          el=dclpy[0]*mr12+dclpy[1]*mr22;
         
         set_el(dy,2*npoints,nvert,el,2*j+1,ic1);
          el=dclpy[0]*mr13+dclpy[1]*mr23;
         
         set_el(dz,2*npoints,nvert,el,2*j+1,ic1);
         
         el=dclpy[2]*mr11+dclpy[3]*mr21;
      
         set_el(dx,2*npoints,nvert,el,2*j+1,ic2);
        el=dclpy[2]*mr12+dclpy[3]*mr22;
       
         set_el(dy,2*npoints,nvert,el,2*j+1,ic2);
          el=dclpy[2]*mr13+dclpy[3]*mr23;
          
         set_el(dz,2*npoints,nvert,el,2*j+1,ic2);
         */
        
         //Chords offset derivatives
        
         //Offset derivatives
         distxdox=(dclpx[0]+dclpx[2]);
         distydox=(dclpy[0]+dclpy[2]);
         distxdoy=(dclpx[1]+dclpx[3]);
         distydoy=(dclpy[1]+dclpy[3]);
         dtox[j]=pow(d2,-0.5)*(distx*doffx[0]+disty*doffx[1]);
         dtoy[j]=pow(d2,-0.5)*(distx*doffy[0]+disty*doffy[1]);
         
         /*
         dtox[2*j]=(dclpx[0]+dclpx[2]);
         dtox[2*j+1]=(dclpy[0]+dclpy[2]);
        
         
         dtoy[2*j]=(dclpx[1]+dclpx[3]);
         dtoy[2*j+1]=(dclpy[1]+dclpy[3]);
        */
         distxdb=(dclpx[0]*get_el(vlist2b,nvert,3,ic1,0)+dclpx[1]*get_el(vlist2b,nvert,3,ic1,1)+dclpx[2]*get_el(vlist2b,nvert,3,ic2,0)+dclpx[3]*get_el(vlist2b,nvert,3,ic2,1));
         distydb=(dclpy[0]*get_el(vlist2b,nvert,3,ic1,0)+dclpy[1]*get_el(vlist2b,nvert,3,ic1,1)+dclpy[2]*get_el(vlist2b,nvert,3,ic2,0)+dclpy[3]*get_el(vlist2b,nvert,3,ic2,1));
         distxdl=(dclpx[0]*get_el(vlist2l,nvert,3,ic1,0)+dclpx[1]*get_el(vlist2l,nvert,3,ic1,1)+dclpx[2]*get_el(vlist2l,nvert,3,ic2,0)+dclpx[3]*get_el(vlist2l,nvert,3,ic2,1));
         distydl=(dclpy[0]*get_el(vlist2l,nvert,3,ic1,0)+dclpy[1]*get_el(vlist2l,nvert,3,ic1,1)+dclpy[2]*get_el(vlist2l,nvert,3,ic2,0)+dclpy[3]*get_el(vlist2l,nvert,3,ic2,1));
         distxdo=(dclpx[0]*get_el(vlist2o,nvert,3,ic1,0)+dclpx[1]*get_el(vlist2o,nvert,3,ic1,1)+dclpx[2]*get_el(vlist2o,nvert,3,ic2,0)+dclpx[3]*get_el(vlist2o,nvert,3,ic2,1));
         distydo=(dclpy[0]*get_el(vlist2o,nvert,3,ic1,0)+dclpy[1]*get_el(vlist2o,nvert,3,ic1,1)+dclpy[2]*get_el(vlist2o,nvert,3,ic2,0)+dclpy[3]*get_el(vlist2o,nvert,3,ic2,1));
         /*
         dangles[2*j*3]=(dclpx[0]*get_el(vlist2b,nvert,3,ic1,0)+dclpx[1]*get_el(vlist2b,nvert,3,ic1,1)+dclpx[2]*get_el(vlist2b,nvert,3,ic2,0)+dclpx[3]*get_el(vlist2b,nvert,3,ic2,1));
         dangles[2*j*3+1]=(dclpx[0]*get_el(vlist2l,nvert,3,ic1,0)+dclpx[1]*get_el(vlist2l,nvert,3,ic1,1)+dclpx[2]*get_el(vlist2l,nvert,3,ic2,0)+dclpx[3]*get_el(vlist2l,nvert,3,ic2,1));
         dangles[2*j*3+2]=(dclpx[0]*get_el(vlist2o,nvert,3,ic1,0)+dclpx[1]*get_el(vlist2o,nvert,3,ic1,1)+dclpx[2]*get_el(vlist2o,nvert,3,ic2,0)+dclpx[3]*get_el(vlist2o,nvert,3,ic2,1));
         
         dangles[2*j*3+3+0]=(dclpy[0]*get_el(vlist2b,nvert,3,ic1,0)+dclpy[1]*get_el(vlist2b,nvert,3,ic1,1)+dclpy[2]*get_el(vlist2b,nvert,3,ic2,0)+dclpy[3]*get_el(vlist2b,nvert,3,ic2,1));
         dangles[2*j*3+3+1]=(dclpy[0]*get_el(vlist2l,nvert,3,ic1,0)+dclpy[1]*get_el(vlist2l,nvert,3,ic1,1)+dclpy[2]*get_el(vlist2l,nvert,3,ic2,0)+dclpy[3]*get_el(vlist2l,nvert,3,ic2,1));
         dangles[2*j*3+3+2]=(dclpy[0]*get_el(vlist2o,nvert,3,ic1,0)+dclpy[1]*get_el(vlist2o,nvert,3,ic1,1)+dclpy[2]*get_el(vlist2o,nvert,3,ic2,0)+dclpy[3]*get_el(vlist2o,nvert,3,ic2,1));
         */
         dangles[3*j]=pow(d2,-0.5)*(distx*distxdb+disty*distydb);
         dangles[3*j+1]=pow(d2,-0.5)*(distx*distxdl+disty*distydl);
         dangles[3*j+2]=pow(d2,-0.5)*(distx*distxdo+disty*distydo);
     }
     else
     {
         printf("No intersection with Boundary contour. j :%d \n",j);
        dist[j]=LARGE_VALUE;
       
     }
 }

 free(vlist2);
 free(vlist2b);
 free(vlist2l);
 free(vlist2o);
 
}
