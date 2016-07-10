#include"utils.h"
#include"matrix_ops.h"
void AdjFacet(int *tlist,double *vlist,int nfac,int nvert,int *A)
{
    /*A is nvertxnvert matrix, facets with common edge (i1,i2) are A(i1,i2), A(i2,i1)*/
    int i1,i2,i3;
    for(int j=0;j<nfac;j++)
    {
        i1=tlist[3*j]-1;
        i2=tlist[3*j+1]-1;
        i3=tlist[3*j+2]-1;
        set_elI(A,nvert,nvert,j+1,i1,i2);
        set_elI(A,nvert,nvert,j+1,i2,i3);
        set_elI(A,nvert,nvert,j+1,i3,i1);
    }
}
      
void Fit_Occ(int *tlist,double *vlist,int nfac,int nvert,double *angles,double *up,double *E,double *V0,double *TIME,double *offset,double *chords,int *type,int nchords,double *W,double *Chordoffset,double *dist,double *dx,double *dy,double *dz,double *dangles,double *dtox,double *dtoy,double* dCOdoff)
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
     * dtox,dtoy 4*nchordsx1, derivatives wrt offsets
     * dCOdoff 4*nchordsxnchords matrix, derivatives for chord offsets
     * */
    double time=0;
    double M[3][3],MRT[3][3],MRdb[3][3],MRdl[3][3],MRdo[3][3],MRTT[9];
    double MRdbT[9],MRdlT[9],MRdoT[9];
    double R[3][3],Rb[3][3],Rl[3][3],Ro[3][3],Rt[3][3],RbT[3][3],RlT[3][3],RoT[3][3];
    int faedge[2],cledge[2],inters;
    int edget[2];
    int ic1,ic2,fa1,fa2;
    double clpoint[2],fapoint[2],dclpx[4],dclpy[4],dfapx[4],dfapy[4];
    double temp4[4];
    double V[3];
    double dLx=0,dLy=0; //Chord offset
    for(int j=0;j<2*nchords;j++)
        time+=TIME[j];
    time=time/(2*nchords);
    int *Adj=calloc(nvert*nvert,sizeof(int));
    AdjFacet(tlist,vlist,nfac,nvert,Adj);
    Calculate_Frame_Matrix(E,up,M);
    rotate(angles[0],angles[1],angles[2],angles[3],time,R,Rb,Rl,Ro);
    //Calculate velocity in camera frame 
    V[0]=M[0][0]*V0[0]+M[0][1]*V0[1]+M[0][2]*V0[2];
    V[1]=M[1][0]*V0[0]+M[1][1]*V0[1]+M[1][2]*V0[2];
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
 
zero_array(dx,4*nchords*nvert);
zero_array(dy,4*nchords*nvert);
zero_array(dz,4*nchords*nvert);
zero_array(dangles,4*nchords*3);
zero_array(dtox,4*nchords);
zero_array(dtox,4*nchords);
zero_array(dCOdoff,4*nchords*nchords);
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
 for(int j=0;j<nchords;j++)
 {
     if(Chordoffset!=NULL)
     {
     dLx=V[0]*Chordoffset[j];
    dLy=V[1]*Chordoffset[j];
     }
     else
     {
         dLx=0;
         dLy=0;
     }
     if(W!=NULL)
         w=W[j];
     else
         w=1;
     a[0]=chords[4*j]+dLx;
     a[1]=chords[4*j+1]+dLy;
     b[0]=chords[4*j+2]+dLx;
     b[1]=chords[4*j+3]+dLy;
    
     find_chord(tlist,vlist2,nfac,nvert,offset,a,b,Adj,cledge,faedge,clpoint,fapoint,&inters,dclpx,dclpy,dfapx,dfapy);
     if(type[j]==-1)
     {
         if(inters==1)
         {
           dist[4*j]=w*1e3;
            dist[4*j+1]=w*1e3;
            dist[4*j+2]=w*1e3;
            dist[4*j+3]=w*1e3;
            printf("Intersection with no-detection chord %d (of %d)\n",j,nchords);
           
         }
         else
         {
         dist[4*j]=0;
            dist[4*j+1]=0;
            dist[4*j+2]=0;
            dist[4*j+3]=0;
         }
         continue;
       
     }
         
     if(inters==1)
     {
         //We have check that closest point is actually closest and not farthest. In that case, switch points
         if((fapoint[0]-clpoint[0])*(b[0]-a[0])+(fapoint[1]-clpoint[1])*(b[1]-a[1])<0)
         {
             memcpy(edget,cledge,sizeof(int)*2);
             memcpy(cledge,faedge,sizeof(int)*2);
             memcpy(faedge,edget,sizeof(int)*2);
             
             memcpy(temp4,clpoint,sizeof(double)*2);
             memcpy(clpoint,fapoint,sizeof(double)*2);
             memcpy(fapoint,temp4,sizeof(double)*2);
             
             memcpy(temp4,dclpx,sizeof(double)*4);
             memcpy(dclpx,dfapx,sizeof(double)*4);
             memcpy(dfapx,temp4,sizeof(double)*4);
             
             memcpy(temp4,dclpy,sizeof(double)*4);
             memcpy(dclpy,dfapy,sizeof(double)*4);
             memcpy(dfapy,temp4,sizeof(double)*4);
         }
             
             
         
         dist[4*j]=w*(clpoint[0]-chords[4*j]-dLx); //distance in x
         dist[4*j+1]=w*(clpoint[1]-chords[4*j+1]-dLy); //dist in y
         dist[4*j+2]=w*(fapoint[0]-chords[4*j+2]-dLx);
         dist[4*j+3]=w*(fapoint[1]-chords[4*j+3]-dLy);
         ic1=cledge[0];
         ic2=cledge[1];
         fa1=faedge[0];
         fa2=faedge[1];
         
         el=dclpx[0]*mr11+dclpx[1]*mr21;
         el=w*el;
         set_el(dx,4*nchords,nvert,el,4*j,ic1);
          el=dclpx[0]*mr12+dclpx[1]*mr22;
          el=w*el;
         set_el(dy,4*nchords,nvert,el,4*j,ic1);
          el=dclpx[0]*mr13+dclpx[1]*mr23;
          el=w*el;
         set_el(dz,4*nchords,nvert,el,4*j,ic1);
         
          el=dclpx[2]*mr11+dclpx[3]*mr21;
          el=w*el;
         set_el(dx,4*nchords,nvert,el,4*j,ic2);
        el=dclpx[2]*mr12+dclpx[3]*mr22;
        el=w*el;
         set_el(dy,4*nchords,nvert,el,4*j,ic2);
          el=dclpx[2]*mr13+dclpx[3]*mr23;
          el=w*el;
         set_el(dz,4*nchords,nvert,el,4*j,ic2);
         
         el=dclpy[0]*mr11+dclpy[1]*mr21;
         el=w*el;
         set_el(dx,4*nchords,nvert,el,4*j+1,ic1);
          el=dclpy[0]*mr12+dclpy[1]*mr22;
          el=w*el;
         set_el(dy,4*nchords,nvert,el,4*j+1,ic1);
          el=dclpy[0]*mr13+dclpy[1]*mr23;
          el=w*el;
         set_el(dz,4*nchords,nvert,el,4*j+1,ic1);
         
         el=dclpy[2]*mr11+dclpy[3]*mr21;
         el=w*el;
         set_el(dx,4*nchords,nvert,el,4*j+1,ic2);
        el=dclpy[2]*mr12+dclpy[3]*mr22;
        el=w*el;
         set_el(dy,4*nchords,nvert,el,4*j+1,ic2);
          el=dclpy[2]*mr13+dclpy[3]*mr23;
          el=w*el;
         set_el(dz,4*nchords,nvert,el,4*j+1,ic2);
         
          el=dfapx[0]*mr11+dfapx[1]*mr21;
          el=w*el;
         set_el(dx,4*nchords,nvert,el,4*j+2,fa1);
          el=dfapx[0]*mr12+dfapx[1]*mr22;
          el=w*el;
         set_el(dy,4*nchords,nvert,el,4*j+2,fa1);
          el=dfapx[0]*mr13+dfapx[1]*mr23;
          el=w*el;
         set_el(dz,4*nchords,nvert,el,4*j+2,fa1);
         
         el=dfapx[2]*mr11+dfapx[3]*mr21;
         el=w*el;
         set_el(dx,4*nchords,nvert,el,4*j+2,fa2);
          el=dfapx[2]*mr12+dfapx[3]*mr22;
          el=w*el;
         set_el(dy,4*nchords,nvert,el,4*j+2,fa2);
          el=dfapx[2]*mr13+dfapx[3]*mr23;
          el=w*el;
         set_el(dz,4*nchords,nvert,el,4*j+2,fa2);
         
           el=dfapy[0]*mr11+dfapy[1]*mr21;
           el=w*el;
         set_el(dx,4*nchords,nvert,el,4*j+3,fa1);
          el=dfapy[0]*mr12+dfapy[1]*mr22;
          el=w*el;
         set_el(dy,4*nchords,nvert,el,4*j+3,fa1);
          el=dfapy[0]*mr13+dfapy[1]*mr23;
          el=w*el;
         set_el(dz,4*nchords,nvert,el,4*j+3,fa1);
         
          el=dfapy[2]*mr11+dfapy[3]*mr21;
          el=w*el;
         set_el(dx,4*nchords,nvert,el,4*j+3,fa2);
          el=dfapy[2]*mr12+dfapy[3]*mr22;
          el=w*el;
         set_el(dy,4*nchords,nvert,el,4*j+3,fa2);
          el=dfapy[2]*mr13+dfapy[3]*mr23;
          el=w*el;
         set_el(dz,4*nchords,nvert,el,4*j+3,fa2);
         //Chords offset derivatives
         set_el(dCOdoff,4*nchords,nchords,-V[0],4*j,j);
         set_el(dCOdoff,4*nchords,nchords,-V[1],4*j+1,j);
         set_el(dCOdoff,4*nchords,nchords,-V[0],4*j+2,j);
         set_el(dCOdoff,4*nchords,nchords,-V[1],4*j+3,j);
         //Offset derivatives
         dtox[4*j]=w*(dclpx[0]+dclpx[2]);
         dtox[4*j+1]=w*(dclpy[0]+dclpy[2]);
         dtox[4*j+2]=w*(dfapx[0]+dfapx[2]);
         dtox[4*j+3]=w*(dfapy[0]+dfapy[2]);
         
         dtoy[4*j]=w*(dclpx[1]+dclpx[3]);
         dtoy[4*j+1]=w*(dclpy[1]+dclpy[3]);
         dtoy[4*j+2]=w*(dfapx[1]+dfapx[3]);
         dtoy[4*j+3]=w*(dfapy[1]+dfapy[3]);
         
         dangles[4*j*3]=w*(dclpx[0]*get_el(vlist2b,nvert,3,ic1,0)+dclpx[1]*get_el(vlist2b,nvert,3,ic1,1)+dclpx[2]*get_el(vlist2b,nvert,3,ic2,0)+dclpx[3]*get_el(vlist2b,nvert,3,ic2,1));
         dangles[4*j*3+1]=w*(dclpx[0]*get_el(vlist2l,nvert,3,ic1,0)+dclpx[1]*get_el(vlist2l,nvert,3,ic1,1)+dclpx[2]*get_el(vlist2l,nvert,3,ic2,0)+dclpx[3]*get_el(vlist2l,nvert,3,ic2,1));
         dangles[4*j*3+2]=w*(dclpx[0]*get_el(vlist2o,nvert,3,ic1,0)+dclpx[1]*get_el(vlist2o,nvert,3,ic1,1)+dclpx[2]*get_el(vlist2o,nvert,3,ic2,0)+dclpx[3]*get_el(vlist2o,nvert,3,ic2,1));
         
         dangles[4*j*3+3+0]=w*(dclpy[0]*get_el(vlist2b,nvert,3,ic1,0)+dclpy[1]*get_el(vlist2b,nvert,3,ic1,1)+dclpy[2]*get_el(vlist2b,nvert,3,ic2,0)+dclpy[3]*get_el(vlist2b,nvert,3,ic2,1));
         dangles[4*j*3+3+1]=w*(dclpy[0]*get_el(vlist2l,nvert,3,ic1,0)+dclpy[1]*get_el(vlist2l,nvert,3,ic1,1)+dclpy[2]*get_el(vlist2l,nvert,3,ic2,0)+dclpy[3]*get_el(vlist2l,nvert,3,ic2,1));
         dangles[4*j*3+3+2]=w*(dclpy[0]*get_el(vlist2o,nvert,3,ic1,0)+dclpy[1]*get_el(vlist2o,nvert,3,ic1,1)+dclpy[2]*get_el(vlist2o,nvert,3,ic2,0)+dclpy[3]*get_el(vlist2o,nvert,3,ic2,1));
         
         dangles[4*j*3+2*3+0]=w*(dfapx[0]*get_el(vlist2b,nvert,3,fa1,0)+dfapx[1]*get_el(vlist2b,nvert,3,fa1,1)+dfapx[2]*get_el(vlist2b,nvert,3,fa2,0)+dfapx[3]*get_el(vlist2b,nvert,3,fa2,1));
         dangles[4*j*3+2*3+1]=w*(dfapx[0]*get_el(vlist2l,nvert,3,fa1,0)+dfapx[1]*get_el(vlist2l,nvert,3,fa1,1)+dfapx[2]*get_el(vlist2l,nvert,3,fa2,0)+dfapx[3]*get_el(vlist2l,nvert,3,fa2,1));
         dangles[4*j*3+2*3+2]=w*(dfapx[0]*get_el(vlist2o,nvert,3,fa1,0)+dfapx[1]*get_el(vlist2o,nvert,3,fa1,1)+dfapx[2]*get_el(vlist2o,nvert,3,fa2,0)+dfapx[3]*get_el(vlist2o,nvert,3,fa2,1));
         
         dangles[4*j*3+3*3+0]=w*(dfapy[0]*get_el(vlist2b,nvert,3,fa1,0)+dfapy[1]*get_el(vlist2b,nvert,3,fa1,1)+dfapy[2]*get_el(vlist2b,nvert,3,fa2,0)+dfapy[3]*get_el(vlist2b,nvert,3,fa2,1));
         dangles[4*j*3+3*3+1]=w*(dfapy[0]*get_el(vlist2l,nvert,3,fa1,0)+dfapy[1]*get_el(vlist2l,nvert,3,fa1,1)+dfapy[2]*get_el(vlist2l,nvert,3,fa2,0)+dfapy[3]*get_el(vlist2l,nvert,3,fa2,1));
         dangles[4*j*3+3*3+2]=w*(dfapy[0]*get_el(vlist2o,nvert,3,fa1,0)+dfapy[1]*get_el(vlist2o,nvert,3,fa1,1)+dfapy[2]*get_el(vlist2o,nvert,3,fa2,0)+dfapy[3]*get_el(vlist2o,nvert,3,fa2,1));
     }
     else
     {
        
         dist[4*j]=w*chords[j*4];
         dist[4*j+1]=w*chords[j*4+1];
         dist[4*j+2]=w*chords[j*4+2];
         dist[4*j+3]=w*chords[j*4+3];
     }
 }
 free(vlist2);
 free(vlist2b);
 free(vlist2l);
 free(vlist2o);
 free(Adj);
}
