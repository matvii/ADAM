#include"utils.h"
#include"matrix_ops.h"
#include"globals.h"

 void find_edge(int *tlist,double *vlist2,int nfac,int nvert,int *visible,double *offset,double *a,double *b,int *cledge,double *clpoint,int *inters,double *dclpx,double *dclpy,double *doffx,double *doffy);
 void Find_Boundary(int *tlist,double *vlist,int nfac,int nvert,double *angles,double *up,double *E,double *E0,double TIME,int *edges,int* nedges,int *Bvert);
void Fit_Contour(int *tlist,double *vlist,int nfac,int nvert,double *angles,double *up,double *E,double *E0,double TIME,double *offset,double *datax,double *datay,int npoints,double *dist,double *dx,double *dy,double *dz,double *dangles,double *dtox,double *dtoy)
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
    int vind;
    int ic1,ic2,fa1,fa2;
    double LARGE_VALUE=1e4;
    double clpoint[2],dclpx[4],dclpy[4];
    double temp4[4];
    double V[3];
    double off[]={0,0};
     double nd_dist,nd_dx,nd_dy;
     int nd_vert;
     int *Edges=calloc(3*nfac/2*2,sizeof(int));
     int nedges;
     int *Bvert=calloc(nvert,sizeof(int));
   Find_Boundary(tlist,vlist,nfac,nvert,angles,up,E,E0,TIME,Edges,&nedges,Bvert);
 // printf("Vertices:\n");
 // print_matrixI(Bvert,1,nvert);
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
   //FindActualBlockers(tlist,vlist,nfac,nvert,Er,E0r,1,visible);
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
 zero_array(dist,npoints+nvert);
zero_array(dx,(npoints+nvert)*nvert);
zero_array(dy,(npoints+nvert)*nvert);
zero_array(dz,(npoints+nvert)*nvert);
zero_array(dangles,(npoints+nvert)*3);
zero_array(dtox,npoints+nvert);
zero_array(dtox,npoints+nvert);

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
 double distdx1,distdy1,distdz1,distdx2,distdy2,distdz2;

 double distxdox,distxdoy,distydox,distydoy;
 double distxdb,distxdl,distxdo;
 double distydb,distydl,distydo;
 double doffx[2],doffy[2];
 a[0]=0;
 a[1]=0;
 double *dist2=calloc(nvert,sizeof(double));
 int *vertices=calloc(nvert,sizeof(int));
 
 //Model boundary edges are contained in Edges list
 //Model boundary vertices in Bvert
 //For each data point, we find the closest boundary edge. 
 //We add offset to the data
 int cedge=-1;
 double edge_dist=1e9;
 double cdist;
 int *closest_edge=calloc(npoints,sizeof(int));
 double *v0,*v1;
 double v[2];
 double N0,N1;
 double nN;
 double t;
 int i1,i2;
double P0,P1;
double v00,v01,v10,v11,b0,b1;
double dP0dx0,dP0dx1,dP0dy0,dP0dy1;
double dP1dx0,dP1dx1,dP1dy0,dP1dy1;

double dtdx0,dtdx1,dtdy0,dtdy1,dtdz0,dtdz1;
double distdx0,distdy0;
double dtdb0,dtdb1;
int cl_vertex;
double vertex_dist,vdist;
//For each data point, find closest vertex, and closest edge.
//Choose whichever is closest
 for(int j=0;j<npoints;j++)
 {
     b0=datax[j]+offset[0];
     b1=datay[j]+offset[1];
     
     //First find closest vertex
     vertex_dist=1e9;
     for(int h=0;h<nvert;h++)
     {
         if(Bvert[h]==0) //Not a boundary vertex
             continue;
         v0=vlist2+3*h;
         vdist=sqrt(pow(b0-v0[0],2)+pow(b1-v0[1],2));
         if(vdist<vertex_dist)
         {
             vertex_dist=vdist;
             cl_vertex=h;
         }
     }
     //Find closest edege
     edge_dist=1e9;
     for(int k=0;k<nedges;k++)
     {
         v0=vlist2+Edges[2*k]*3;
         v1=vlist2+Edges[2*k+1]*3;
         v00=v0[0];
        v01=v0[1];
        v10=v1[0];
        v11=v1[1];
         N0=v1[0]-v0[0];
         N1=v1[1]-v0[1];
         nN=pow(N0,2)+pow(N1,2);
        t=((b0-v00)*N0+(b1-v01)*N1)/nN;
       
         if(t<0 || t>1) //Point projection is  outside of the edge
             continue;
          P0=v00+t*N0;
     P1=v01+t*N1;
     cdist=sqrt(pow((b0-P0),2)+pow(b1-P1,2));
    
         if(cdist<edge_dist)
         {
            
             edge_dist=cdist;
             closest_edge[j]=k; //Closest edge is determined by the vertex pair(Edges[2*j]->Edges[2*j+1])
         }
         
         //TBD: REMEMBER TO CHECK THE PROJECTION ROUTINE AND SIGN OF t
     }
     //For each point, calculate final distances and derivatives   
     if(edge_dist<vertex_dist) //So edge is closest
     {
       //  printf("Point %d, Edge %d is closest\n",j,closest_edge[j]);
     cedge=closest_edge[j];
     ic1=Edges[2*cedge];
     ic2=Edges[2*cedge+1];
     v0=vlist2+ic1*3;
     v1=vlist2+ic2*3;
     v00=v0[0];
     v01=v0[1];
     v10=v1[0];
     v11=v1[1];
     
     N0=v10-v00;
     N1=v11-v01;
     nN=pow(N0,2)+pow(N1,2);
     t=((b0-v00)*N0+(b1-v01)*N1)/nN;
     P0=v00+t*N0;
     P1=v01+t*N1;
     cdist=pow((b0-P0),2)+pow(b1-P1,2);
     dist[j]=sqrt(cdist);        
  
     //Derivatives
     dtdx0=((-b0 + 2*v00 - v10)*(pow(v00 - v10,2) + pow(v01 - v11,2)) + 2*(-v00 + v10)*((b0 - v00)*(-v00 + v10) + (b1 - v01)*(-v01 + v11)))/pow(nN,2);
     dtdx1=((b0 - v00)*(pow(v00 - v10,2) + pow(v01 - v11,2)) - 2*(-v00 + v10)*((b0 - v00)*(-v00 + v10) + (b1 - v01)*(-v01 + v11)))/pow(nN,2);
     dtdy0=((pow(v00 - v10,2) + pow(v01 - v11,2))*(-b1 + 2*v01 - v11) + 2*(-v01 + v11)*((b0 - v00)*(-v00 + v10) + (b1 - v01)*(-v01 + v11)))/pow(nN,2);
     dtdy1=((b1 - v01)*(pow(v00 - v10,2) +pow(v01 - v11,2)) - 2*(-v01 + v11)*((b0 - v00)*(-v00 + v10) + (b1 - v01)*(-v01 + v11)))/pow(nN,2);
     dtdb0=N0/nN;
     dtdb1=N1/nN;
     dP0dx0=1+dtdx0*N0-t;
     dP0dx1=0+dtdx1*N0+t;
     dP1dx0=0+dtdx0*N1+0;
     dP1dx1=0+dtdx1*N1+0;
     dP0dy0=0+dtdy0*N0+0;
     dP0dy1=0+dtdy1*N0+0;
     dP1dy0=1+dtdy0*N1-t;
     dP1dy1=0+dtdy1*N1+t;
     distdx0=1/dist[j]*(-(b0-P0)*dP0dx0-(b1-P1)*dP1dx0);
     distdy0=1/dist[j]*(-(b0-P0)*dP0dy0-(b1-P1)*dP1dy0);
     distdx1=1/dist[j]*(-(b0-P0)*dP0dx1-(b1-P1)*dP1dx1);
     distdy1=1/dist[j]*(-(b0-P0)*dP0dy1-(b1-P1)*dP1dy1);
     
     //dx is npoints x nvert
     dx[j*nvert+ic1]=distdx0*mr11+distdy0*mr21;
     dy[j*nvert+ic1]=distdx0*mr12+distdy0*mr22;
     dz[j*nvert+ic1]=distdx0*mr13+distdy0*mr23;
     
     dx[j*nvert+ic2]=distdx1*mr11+distdy1*mr21;
     dy[j*nvert+ic2]=distdx1*mr12+distdy1*mr22;
     dz[j*nvert+ic2]=distdx1*mr13+distdy1*mr23;
     
     //Derivatives wrt angles. //npointsx3 matrix
     
     dangles[j*3+0]=distdx0*vlist2b[3*ic1+0]+distdy0*vlist2b[3*ic1+1]+distdx1*vlist2b[3*ic2+0]+distdy1*vlist2b[3*ic2+1];
     dangles[j*3+1]=distdx0*vlist2l[3*ic1+0]+distdy0*vlist2l[3*ic1+1]+distdx1*vlist2l[3*ic2+0]+distdy1*vlist2l[3*ic2+1];
     dangles[j*3+2]=distdx0*vlist2o[3*ic1+0]+distdy0*vlist2o[3*ic1+1]+distdx1*vlist2o[3*ic2+0]+distdy1*vlist2o[3*ic2+1];
     //Derivatives wrt offsets. //npoints x 1 matrix
     dtox[j]=1/dist[j]*((b0-P0)*(1-N0*dtdb0)+(b1-P1)*(-N1*dtdb0));
     dtoy[j]=1/dist[j]*((b0-P0)*(-N0*dtdb1)+(b1-P1)*(1-N1*dtdb1));
     }
     else //Vertex is closest
     {
       //  printf("Point %d, vertex %d is closest\n",j,cl_vertex);
         ic1=cl_vertex;
         v0=vlist2+3*ic1;
         dist[j]=vertex_dist;
         //sqrt(pow(b0-v0[0],2)+pow(b1-v0[1],2));
         //Derivatives
         dx[j*nvert+ic1]=-1/vertex_dist*((b0-v0[0])*mr11+(b1-v0[1])*mr21);
         dy[j*nvert+ic1]=-1/vertex_dist*((b0-v0[0])*mr12+(b1-v0[1])*mr22);
         dz[j*nvert+ic1]=-1/vertex_dist*((b0-v0[0])*mr13+(b1-v0[1])*mr23);
         
         dangles[j*3+0]=-1/vertex_dist*((b0-v0[0])*vlist2b[3*ic1]+(b1-v0[1])*vlist2b[3*ic1+1]);
         dangles[j*3+1]=-1/vertex_dist*((b0-v0[0])*vlist2l[3*ic1]+(b1-v0[1])*vlist2l[3*ic1+1]);
         dangles[j*3+2]=-1/vertex_dist*((b0-v0[0])*vlist2o[3*ic1]+(b1-v0[1])*vlist2o[3*ic1+1]);
         
         dtox[j]=1/vertex_dist*(b0-v0[0]);
         dtoy[j]=1/vertex_dist*(b1-v0[1]);
     }
}
//Now we find closest data point for each point on boundary

for(int j=0;j<nvert;j++)
 {
    // printf("Vertex %d (%7.2f,%7.2f), Boundary: %d\n",j,vlist2[3*j],vlist2[3*j+1],Bvert[j]);
     if(Bvert[j]==0)
         continue;
     v[0]=vlist2[3*j]-offset[0];
     v[1]=vlist2[3*j+1]-offset[1];
     vind=find_closest(v,datax,datay,npoints);
     distx=v[0]-datax[vind];
     disty=v[1]-datay[vind];
     
     //printf("For vertex %d, (%7.2f,%7.2f) %d (%7.2f,%7.2f) is closest\n",j,vlist2[3*j],vlist2[3*j+1],vind,datax[vind],datay[vind]);   
    d2=pow(distx,2)+pow(disty,2);
    dist[j+npoints]=sqrt(d2);
    dx[(npoints+j)*nvert+j]=pow(d2,-0.5)*(distx*mr11+disty*mr21);
    dy[(npoints+j)*nvert+j]=pow(d2,-0.5)*(distx*mr12+disty*mr22);
    dz[(npoints+j)*nvert+j]=pow(d2,-0.5)*(distx*mr13+disty*mr23);
    
    dtox[npoints+j]=-pow(d2,-0.5)*distx;
    dtoy[npoints+j]=-pow(d2,-0.5)*disty;
    
    dangles[3*(j+npoints)]=pow(d2,-0.5)*(distx*vlist2b[3*j]+disty*vlist2b[3*j+1]);
    dangles[3*(j+npoints)+1]=pow(d2,-0.5)*(distx*vlist2l[3*j]+disty*vlist2l[3*j+1]);
    dangles[3*(j+npoints)+2]=pow(d2,-0.5)*(distx*vlist2o[3*j]+disty*vlist2o[3*j+1]);
    //Derivatives
 }
 free(vlist2);
 free(vlist2b);
 free(vlist2l);
 free(vlist2o);
 free(closest_edge);
 free(dist2);
 free(vertices);
 free(Edges);
 free(Bvert);
 free(visible);
}
/*
int main()
{
    int *tlist,*tlistn,nfac,nvert,nfacn,nvertn;
    double *vlist,*vlist2,*vlistn,*D;
    CNTRstruct *CR;
    double angles[]={0.5,0.6,5,0};
    double angles2[]={0.5,0.6,5,0};
    double TIME=1.2;
    //angles[0]=0;
    //angles[1]=0;
    double offset[]={0.2,-0.1};
     double offset2[]={0.2,-0.1};
    CR=read_contour("test_cont",0,0);
    int nrows=2;
    nfac=8*pow(nrows,2);
        nvert=4*pow(nrows,2)+2;
        vlist=calloc(3*nvert,sizeof(double));
        
        tlist=calloc(3*nfac,sizeof(int));
        generate_ellipsoid(nrows,3,2,1,tlist,vlist);
        vlist2=calloc(3*nvert,sizeof(double));
        memcpy(vlist2,vlist,3*nvert*sizeof(double));
       int j=0;
       int *nobs=CR->nobs;
        double *dx=calloc((nvert+nobs[j])*nvert,sizeof(double));
        double *dy=calloc((nvert+nobs[j])*nvert,sizeof(double));
        double *dz=calloc((nvert+nobs[j])*nvert,sizeof(double));
        double *dangles=calloc((nvert+nobs[j])*3,sizeof(double));
        double *dtox=calloc((nvert+nobs[j]),sizeof(double));
        double *dtoy=calloc((nvert+nobs[j]),sizeof(double));
        double *dist=calloc(CR->ntotal+nvert,sizeof(double));
        
        double *dist2=calloc(CR->ntotal+nvert,sizeof(double));
        double eps=1e-7;
       
     Fit_Contour(tlist,vlist,nfac,nvert,angles,CR->up,CR->E,CR->E0,TIME,offset,CR->datax[0],CR->datay[0],CR->nobs[0],dist,dx,dy,dz,dangles,dtox,dtoy);
      print_matrix(dangles,nobs[j]+nvert,3);
//      
       offset2[0]=offset[0];
       offset2[1]=offset[1]+eps;
       vlist2[6*3+2]+=eps; 
    printf("\n");  
      angles2[2]+=eps;
       Fit_Contour(tlist,vlist,nfac,nvert,angles2,CR->up,CR->E,CR->E0,TIME,offset,CR->datax[0],CR->datay[0],CR->nobs[0],dist2,dx,dy,dz,dangles,dtox,dtoy);
// //       
        for(int j=0;j<nobs[0]+nvert;j++)
            printf(" %7.4f ",(dist2[j]-dist[j])/eps);
            printf("\n");
}

*/
