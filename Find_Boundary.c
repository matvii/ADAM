#include"utils.h"
#include"matrix_ops.h"
#include"globals.h"

 void find_edge_noderiv(int *tlist,double *vlist2,int nfac,int nvert,int *visible,double *a,double *b,int *cledge,double *clpoint,int *inters);     
void Find_Boundary(int *tlist,double *vlist,int nfac,int nvert,double *angles,double *up,double *E,double *E0,double TIME,int *edges,int* nedges,int *Bvert)
{
    /*INPUT:
     * Allocated edges to 2*(3*nfac/2) mult_vector
     * Bvert 1xnvert vector
     */
    /*OUTPUT:
     * Boundary edges: edges[2*i]->edges[2*i+1]
     * Bvert[j]==1 if vertex is on the boundary
     */
    
    double M[3][3],MRT[3][3],MRdb[3][3],MRdl[3][3],MRdo[3][3],MRTT[9];
    double Er[3],E0r[3];
    double MRdbT[9],MRdlT[9],MRdoT[9];
    double R[3][3],Rb[3][3],Rl[3][3],Ro[3][3],Rt[3][3],RbT[3][3],RlT[3][3],RoT[3][3];
    int cledge[2],inters;
    int edget[2];
    int ic1,ic2,fa1,fa2;
    double LARGE_VALUE=1e4;
    double clpoint[2];
    double temp4[4];
    double V[3];
    double off[]={0,0};
     double nd_dist,nd_dx,nd_dy;
     int nd_vert;
   int NDIV=100; //Number of lines used for determining boundary
   int npoints=NDIV;
   double *datax=calloc(NDIV,sizeof(double));
   double *datay=calloc(NDIV,sizeof(double));
   for(int j=0;j<NDIV;j++)
   {
       datax[j]=cos(2*PI*j/NDIV);
       datay[j]=sin(2*PI*j/NDIV);
   }
   int *Bmat=calloc(nvert*nvert,sizeof(int)); //something for indexing edges 
    
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
 
 mult_mat(M,Rt,MRT); 
double mr11,mr12,mr13,mr21,mr22,mr23;

 


transpose2(MRT,MRTT);

double a[2],b[2];
 double *vlist2=calloc(nvert*3,sizeof(double));
 
 matrix_prod(vlist,nvert,3,MRTT,3,vlist2); //vlist2 contains rotated and projected model
 
  double w;
 double el=0;
 double d2=0;
 double distx;
 double disty;
 
 a[0]=0;
 a[1]=0;
 
 int count=0;
 for(int j=0;j<npoints;j++)
 {
     b[0]=datax[j];
     b[1]=datay[j];
     
     
    
     find_edge_noderiv(tlist,vlist2,nfac,nvert,visible,a,b,cledge,clpoint,&inters);
     
     if(inters==1)
     {
         
             
        
         
         ic1=cledge[0];
         ic2=cledge[1];
         if(Bmat[ic1*nvert+ic2]==0)
         {
             Bmat[ic1*nvert+ic2]=1;
             Bmat[ic2*nvert+ic1]=1;
             edges[2*count]=ic1;
             edges[2*count+1]=ic2;
             count++;
         }
         Bvert[ic1]=1; //Save the boundary
         Bvert[ic2]=1;
        
     }
 }
 *nedges=count;
 
 
 free(vlist2);
 free(datax);
 free(datay);
 free(Bmat);
 free(visible);
}

// int main()
// {
//     int *tlist,*tlistn,nfac,nvert,nfacn,nvertn;
//     double *vlist,*vlist2,*vlistn,*D;
//     CNTRstruct *CR;
//     double angles[]={0.1,0.3,5,0};
//     double angles2[]={0.1,0.3,5,0};
//    
//     double offset[]={0.1,0.2};
//      double offset2[]={0,0};
//     CR=read_contour("test_cont",0,0);
//     CR->E0[0]=1/sqrt(3);
//     CR->E0[1]=1/sqrt(3);
//     CR->E0[2]=1/sqrt(3);
//     int nrows=2;
//     nfac=8*pow(nrows,2);
//         nvert=4*pow(nrows,2)+2;
//         vlist=calloc(3*nvert,sizeof(double));
//         
//         tlist=calloc(3*nfac,sizeof(int));
//         generate_ellipsoid(nrows,3,2,1,tlist,vlist);
//         vlist2=calloc(3*nvert,sizeof(double));
//         memcpy(vlist2,vlist,3*nvert*sizeof(double));
//        int j=0;
//        int *nobs=CR->nobs;
//        /* double *dx=calloc(nobs[j]*nvert,sizeof(double));
//         double *dy=calloc(nobs[j]*nvert,sizeof(double));
//         double *dz=calloc(nobs[j]*nvert,sizeof(double));
//         double *dangles=calloc(nobs[j]*3,sizeof(double));
//         double *dtox=calloc(nobs[j],sizeof(double));
//         double *dtoy=calloc(nobs[j],sizeof(double));
//         double *dist=calloc(CR->ntotal,sizeof(double));
//         
//         double *dist2=calloc(CR->ntotal,sizeof(double));
//         double eps=1e-7;
//        */
//        int *Edges=calloc(3*nfac/2*2,sizeof(int));
//        int nedges;
//        int *Bvert=calloc(nvert,sizeof(int));
//      Find_Boundary(tlist,vlist,nfac,nvert,angles,CR->up,CR->E,CR->E0,CR->TIME[0],Edges,&nedges,Bvert);
//      print_matrixI(Edges,1,2*nedges);
//      printf("Boundary vertices:\n");
//      for(int j=0;j<nvert;j++)
//          if(Bvert[j]==1)
//              printf(" %d ",j);
//     printf("\n");
//      
//      
// }
