#include"prepare.h"
void Calculate_Frame_Matrix(double *E,double *up,double R[3][3]);
void Calculate_Normal_Derivative(double *v1,double *v2,double *v3,double *n,double *n1dx,double *n2dx,double *n3dx,double *n1dy,double *n2dy,double *n3dy,double *n1dz,double *n2dz,double *n3dz,double *area,double *dAdx,double *dAdy,double *dAdz);
void Calculate_AO_deriv(int *tlist,double *vlist,int nfac,int nvert,double *angles,double *Eo,double *E0o,double *up,double TIME,double dist,double *freqx,double *freqy,int nfreq,double *offset,double complex *F,double complex *dFdx,double complex *dFdy,double complex *dFdz,double complex *dFdA,double complex *dFdoff)
{
 double complex *F0,*FTda,*FTdb,*FTdc,*FTdd,*FTdh,*FTdg,*FTdx,*FTdy,*FTdz,*FTdA;
 FTdx=(double complex*)calloc(nfreq*nvert,sizeof(double complex));
 FTdy=(double complex*)calloc(nfreq*nvert,sizeof(double complex));
 FTdz=(double complex*)calloc(nfreq*nvert,sizeof(double complex));
 FTdA=(double complex*)calloc(nfreq*3,sizeof(double complex));
 F0=(double complex*)calloc(nfreq,sizeof(double complex));
 FTda=(double complex*)calloc(nfreq,sizeof(double complex));
 FTdb=(double complex*)calloc(nfreq,sizeof(double complex));
 FTdc=(double complex*)calloc(nfreq,sizeof(double complex));
 FTdd=(double complex*)calloc(nfreq,sizeof(double complex));
 FTdh=(double complex*)calloc(nfreq,sizeof(double complex));
 FTdg=(double complex*)calloc(nfreq,sizeof(double complex));
 
 double M[3][3],dMb[3][3],dMo[3][3],dMl[3][3],Mt[3][3],dMbT[3][3],dMoT[3][3],dMlT[3][3]; //Rotation matrices and their derivatives and transposes
 double R[3][3],Rdb[3][3],Rdl[3][3],Rdo[3][3],RT[3][3]; //Projection matrix, and derivatives of combined projection+rotation matrix
 double dEdb[3],dEdl[3],dEdo[3],dE0db[3],dE0dl[3],dE0do[3];
 double dndx1[3],dndx2[3],dndx3[3],dndy1[3],dndy2[3],dndy3[3],dndz1[3],dndz2[3],dndz3[3]; //Derivatives of the facet normal vector
 double dBdx1,dBdy1,dBdz1,dBdx2,dBdy2,dBdz2,dBdx3,dBdy3,dBdz3,dBdb,dBdl,dBdo; //Derivatives of facet brightness
 double *dTBdx,*dTBdy,*dTBdz,dTBdA[3]; //Derivatives of total brightness, allocating memory
 dTBdx=(double*)calloc(nvert,sizeof(double));
 dTBdy=(double*)calloc(nvert,sizeof(double));
 dTBdz=(double*)calloc(nvert,sizeof(double));
 double dmudx1,dmudy1,dmudz1,dmudx2,dmudy2,dmudz2,dmudx3,dmudy3,dmudz3;
 double dmu0dx1,dmu0dy1,dmu0dz1,dmu0dx2,dmu0dy2,dmu0dz2,dmu0dx3,dmu0dy3,dmu0dz3;
 double dmudl,dmudb,dmudo,dmu0dl,dmu0db,dmu0do; //Derivatives of mu and mu0
 double dAdx[3],dAdy[3],dAdz[3]; //Facet area derivatives
 double dadx,dady,dadz,dbdx,dbdy,dbdz; //Derivatives of projected vertices
 double E[3],E0[3]; //Earth and Sun direction, rotated
 double side1[3],side2[3];
 double n[3];
 double v1db[3],v2db[3],v3db[3],v1dl[3],v2dl[3],v3dl[3],v1do[3],v2do[3],v3do[3]; //Derivatives of 2d vertices wrt angles
 double vr1[3],vr2[3],vr3[3];
 double *v1,*v2,*v3;
 double complex scale;

 double dp;
 double B,TB=0;
 double norm;
 double mut,mu0t;
 int t1,t2,t3;
int j1,j2,j3;
 double mu,mu0,area;
 double *normal;
 int *visible;
 int tb1,tb2,tb3; //Indices to the vertices of possible blocker facet
 int blocked=0;
 //Distance km->arcsec
 dp=1/(dist*149597871.0)*180.0/PI*3600.0;
 visible=(int*)calloc(nfac,sizeof(int));
 //Calculate_Frame_Matrix_Derivatives(Eo,angles,TIME,rfreq,R,Rdb,Rdl,Rdo);
  dadx=R[0][0];
  dady=R[0][1];
  dadz=R[0][2];
  dbdx=R[1][0];
  dbdy=R[1][1];
  dbdz=R[1][2];
  
  //Calculate frame change matrix
  Calculate_Frame_Matrix(Eo,up,R);
 
 //FacetsOverHorizon(tlist,vlist,nfac,nvert,normal,centroid,NumofBlocks,IndexofBlocks);
 
 rotate(angles[0],angles[1],angles[2],0.0,TIME,M,dMb,dMl,dMo);
 //Construct asteroid->Camera frame matrix, which is
 //asteroid->world frame->camera frame
 transpose(M,Mt); //Transpose, since we rotate the model, not view directions
 transpose(dMb,dMbT);
 transpose(dMl,dMlT);
 transpose(dMo,dMoT);
 mult_mat(R,Mt,RT); 
 mult_vector(M,Eo,E);
 mult_vector(M,E0o,E0);
 
 mult_mat(R,dMbT,Rdb);
 mult_mat(R,dMlT,Rdl);
 mult_mat(R,dMoT,Rdo);
//Derivatives of E,E0 wrt beta,lambda,omega
 mult_vector(dMb,Eo,dEdb);
 mult_vector(dMl,Eo,dEdl);
 mult_vector(dMo,Eo,dEdo);
 mult_vector(dMb,E0o,dE0db);
 mult_vector(dMl,E0o,dE0dl);
 mult_vector(dMo,E0o,dE0do);
 /*For each facet,
  * 1)Check if facet is visible
  * 2) Calculate echo
  * 3) Convert triangle to range-Doppler frame
  * 4) Calculate FT
  */
//Find actual blockers
FindActualBlockers(tlist,vlist,nfac,nvert,E,E0,1,visible);
 //visible is nfac vector, visible[j]=1 if facet (j+1)th facet is visible
//NOTE INDEXING
mexPrintf("%f %f %f\n",vlist[0],vlist[1],vlist[2]);
//Main loop
/*
 for(int j=0;j<nfac;j++)
 {
   if(visible[j]==0)
     continue;
  //Calculate normal from facet vertices
   //Vertex indices of the current facet
   //Note that C indices from 0, matlab from 1
   j1=tlist[j*3]-1;
   j2=tlist[j*3+1]-1;
   j3=tlist[j*3+2]-1;
   //Current vertices
   for(int i=0;i<3;i++)
   {
   v1[i]=*(vlist+j1*3+i)*dp; //convert km->arcsec
   v2[i]=*(vlist+j2*3+i)*dp;
   v3[i]=*(vlist+j3*3+i)*dp;
   }
    //Calculate Normal derivatives (in the original frame)
     Calculate_Normal_Derivative(v1,v2,v3,n,dndx1,dndx2,dndx3,dndy1,dndy2,dndy3,dndz1,dndz2,dndz3,&area,dAdx,dAdy,dAdz);
   
   //Calculate normals and centroids
   
   mu=DOT(E,n);
   mu0=DOT(E0,n);
  //Convert to camera frame
   mult_vector(RT,v1,vr1);
   mult_vector(RT,v2,vr2);
   mult_vector(RT,v3,vr3);
   for(int i=0;i<3;i++)
   {
     vr1[i]=dp*vr1[i];
     vr2[i]=dp*vr2[i];
     vr3[i]=dp*vr3[i];
   }
    
     //Now we should convert to frequency domain, ie calculate the contribution of each facet
  //   Calc_FTC(freqx,freqy,nfreq,vr1[0],vr1[1],vr2[0],vr2[1],vr3[0],vr3[1],F0);
     Calc_FTC_deriv(freqx,freqy,nfreq,vr1[0],vr1[1],vr2[0],vr2[1],vr3[0],vr3[1],F0,FTda,FTdb,FTdc,FTdd,FTdg,FTdh);
   
   area=0.5*norm;
    //Derivatives wrt angles
     mult_vector(Rdb,v1,v1db);
     mult_vector(Rdb,v2,v2db);
     mult_vector(Rdb,v3,v3db);
     mult_vector(Rdl,v1,v1dl);
     mult_vector(Rdl,v2,v2dl);
     mult_vector(Rdl,v3,v3dl);
     mult_vector(Rdo,v1,v1do);
     mult_vector(Rdo,v2,v2do);
     mult_vector(Rdo,v3,v3do);
     //Derivatives of mu,mu0
     dmudx1=DOT(E,dndx1);
     dmudx2=DOT(E,dndx2);
     dmudx3=DOT(E,dndx3);
     dmudy1=DOT(E,dndy1);
     dmudy2=DOT(E,dndy2);
     dmudy3=DOT(E,dndy3);
     dmudz1=DOT(E,dndz1);
     dmudz2=DOT(E,dndz2);
     dmudz3=DOT(E,dndz3);
     dmudb=DOT(dEdb,n);
     dmudl=DOT(dEdl,n);
     dmudo=DOT(dEdo,n);
     dmu0dx1=DOT(E0,dndx1);
     dmu0dx2=DOT(E0,dndx2);
     dmu0dx3=DOT(E0,dndx3);
     dmu0dy1=DOT(E0,dndy1);
     dmu0dy2=DOT(E0,dndy2);
     dmu0dy3=DOT(E0,dndy3);
     dmu0dz1=DOT(E0,dndz1);
     dmu0dz2=DOT(E0,dndz2);
     dmu0dz3=DOT(E0,dndz3);
     dmu0db=DOT(dE0db,n);
     dmu0dl=DOT(dE0dl,n);
     dmu0do=DOT(dE0do,n);
   B=mu0*(1.0/(mu+mu0)+0.1); //mu removed here
   //Derivatives of B
   mu0t=(mu/pow(mu+mu0,2)+0.1);
   mut=-mu0/pow(mu+mu0,2);
   dBdx1=mu0t*dmu0dx1+mut*dmudx1;
   dBdx2=mu0t*dmu0dx2+mut*dmudx2;
   dBdx3=mu0t*dmu0dx3+mut*dmudx3;
   
   dBdy1=mu0t*dmu0dy1+mut*dmudy1;
   dBdy2=mu0t*dmu0dy2+mut*dmudy2;
   dBdy3=mu0t*dmu0dy3+mut*dmudy3;
   
   dBdz1=mu0t*dmu0dz1+mut*dmudz1;
   dBdz2=mu0t*dmu0dz2+mut*dmudz2;
   dBdz3=mu0t*dmu0dz3+mut*dmudz3;
   dBdb=mu0t*dmu0db+mut*dmudb;
   dBdl=mu0t*dmu0dl+mut*dmudl;
   dBdo=mu0t*dmu0do+mut*dmudo;
   //Derivative of total brightness
   dTBdx[j1]+=dBdx1*area*mu+B*dAdx[0]*mu+B*area*dmudx1;
   dTBdx[j2]+=dBdx2*area*mu+B*dAdx[1]*mu+B*area*dmudx2;
   dTBdx[j3]+=dBdx3*area*mu+B*dAdx[2]*mu+B*area*dmudx3;
   dTBdy[j1]+=dBdy1*area*mu+B*dAdy[0]*mu+B*area*dmudy1;
   dTBdy[j2]+=dBdy2*area*mu+B*dAdy[1]*mu+B*area*dmudy2;
   dTBdy[j3]+=dBdy3*area*mu+B*dAdy[2]*mu+B*area*dmudy3;
   dTBdz[j1]+=dBdz1*area*mu+B*dAdz[0]*mu+B*area*dmudz1;
   dTBdz[j2]+=dBdz2*area*mu+B*dAdz[1]*mu+B*area*dmudz2;
   dTBdx[j3]+=dBdz3*area*mu+B*dAdz[2]*mu+B*area*dmudz3;
   dTBdA[0]+=dBdb*area*mu+B*area*dmudb;
   dTBdA[1]+=dBdl*area*mu+B*area*dmudl;
   dTBdA[2]+=dBdo*area*mu+B*area*dmudo;
     for(int jf=0;jf<nfreq;jf++)
     {
     F[jf]+=B*F0[jf];
     
       
       FTdx[jf*nvert+j1]+=dBdx1*F0[jf]+B*(FTda[jf]*dadx+FTdb[jf]*dbdx);
       FTdx[jf*nvert+j2]+=dBdx2*F0[jf]+B*(FTdc[jf]*dadx+FTdd[jf]*dbdx);
       FTdx[jf*nvert+j3]+=dBdx3*F0[jf]+B*(FTdg[jf]*dadx+FTdh[jf]*dbdx);
       
       FTdy[jf*nvert+j1]+=dBdy1*F0[jf]+B*(FTda[jf]*dady+FTdb[jf]*dbdy);
       FTdy[jf*nvert+j2]+=dBdy2*F0[jf]+B*(FTdc[jf]*dady+FTdd[jf]*dbdy);
       FTdy[jf*nvert+j3]+=dBdy3*F0[jf]+B*(FTdg[jf]*dady+FTdh[jf]*dbdy);
       
       FTdz[jf*nvert+j1]+=dBdz1*F0[jf]+B*(FTda[jf]*dadz+FTdb[jf]*dbdz);
       FTdz[jf*nvert+j2]+=dBdz2*F0[jf]+B*(FTdc[jf]*dadz+FTdd[jf]*dbdz);
       FTdz[jf*nvert+j3]+=dBdz3*F0[jf]+B*(FTdg[jf]*dadz+FTdh[jf]*dbdz);
        //angle derivatives
       
       FTdA[jf*3+0]+=dBdb*F0[jf]+B*(FTda[jf]*v1db[0]+FTdb[jf]*v1db[1]+FTdc[jf]*v2db[0]+FTdd[jf]*v2db[1]+FTdg[jf]*v3db[0]+FTdh[jf]*v3db[1]);
       FTdA[jf*3+1]+=dBdl*F0[jf]+B*(FTda[jf]*v1dl[0]+FTdb[jf]*v1dl[1]+FTdc[jf]*v2dl[0]+FTdd[jf]*v2dl[1]+FTdg[jf]*v3dl[0]+FTdh[jf]*v3dl[1]);
       FTdA[jf*3+2]+=dBdo*F0[jf]+B*(FTda[jf]*v1do[0]+FTdb[jf]*v1do[1]+FTdc[jf]*v2do[0]+FTdd[jf]*v2do[1]+FTdg[jf]*v3do[0]+FTdh[jf]*v3do[1]);
     }
      TB=TB+B*area*mu;
}

//Normalize with total brightness

for(int j=0;j<nfreq;j++)
{
  scale=cexp(2.0*PI*I*(offset[0]*freqx[j]+offset[1]*freqy[j]));
  for(int k=0;k<nvert;k++)
  {
    dFdx[j*nvert+k]=scale*(FTdx[j*nvert+k]*TB-F[j]*dTBdx[k])/pow(TB,2);
    dFdy[j*nvert+k]=scale*(FTdy[j*nvert+k]*TB-F[j]*dTBdy[k])/pow(TB,2);
    dFdz[j*nvert+k]=scale*(FTdz[j*nvert+k]*TB-F[j]*dTBdz[k])/pow(TB,2);
  }
  dFdA[j*3+0]=scale*(FTdA[j*3+0]*TB-F[j]*dTBdA[0])/pow(TB,2);
  dFdA[j*3+1]=scale*(FTdA[j*3+1]*TB-F[j]*dTBdA[1])/pow(TB,2);
  dFdA[j*3+2]=scale*(FTdA[j*3+2]*TB-F[j]*dTBdA[2])/pow(TB,2);
  
  F[j]=cexp(2.0*PI*I*(offset[0]*freqx[j]+offset[1]*freqy[j]))*F[j]/TB;
  dFdoff[j*2+0]=2.0*PI*I*freqx[j]*F[j];
  dFdoff[j*2+1]=2.0*PI*I*freqy[j]*F[j];
}
*/
free(FTdx);
free(FTdy);
free(FTdz);
free(FTdA);
free(dTBdx);
free(dTBdy);
free(dTBdz);

free(FTda);
free(FTdb);
free(FTdc);
free(FTdd);
free(FTdg);
free(FTdh);
free(F0);
free(visible);
}
void Calculate_Normal_Derivative(double *w1,double *w2,double *w3,double *n,double *dndx1,double *dndx2,double *dndx3,double *dndy1,double *dndy2,double *dndy3,double *dndz1,double *dndz2,double *dndz3,double *area,double *dAdx,double *dAdy,double *dAdz)
{
  //Calculate derivatives of normal of triangle (w1 w2 w3) wrt triangle vertex coordinates
  double x1,x2,x3,y1,y2,y3,z1,z2,z3;
  double v1[3],v2[3],cv[3];
  double normc;
  double dcx1[3],dcx2[3],dcx3[3],dcy1[3],dcy2[3],dcy3[3],dcz1[3],dcz2[3],dcz3[3];
  double *dareax,*dareay,*dareaz;
  dareax=dAdx;
  dareay=dAdy;
  dareaz=dAdz;
  for(int i=0;i<3;i++)
  {
    v1[i]=w2[i]-w1[i];
    v2[i]=w3[i]-w1[i];
  }
  cross(v1,v2,cv);
  normc=NORM(cv);
  area[0]=normc;
  n[0]=cv[0]/normc;
  n[1]=cv[1]/normc;
  n[2]=cv[2]/normc;
  x1=w1[0];
  y1=w1[1];
  z1=w1[2];
  x2=w2[0];
  y2=w2[1];
  z2=w2[2];
  x3=w3[0];
  y3=w3[1];
  z3=w3[2];
  dcx1[0]=0;
  dcx1[1]=z3-z2;
  dcx1[2]=y2-y3;
  
  dcx2[0]=0;
  dcx2[1]=z1-z3;
  dcx2[2]=y3-y1;
  
  dcx3[0]=0;
  dcx3[1]=z2-z1;
  dcx3[2]=y1-y2;
  
  dcy1[0]=z2-z3;
  dcy1[1]=0;
  dcy1[2]=x3-x2;
  
  dcy2[0]=z3-z1;
  dcy2[1]=0;
  dcy2[2]=x1-x3;
  
  dcy3[0]=z1-z2;
  dcy3[1]=0;
  dcy3[2]=x2-x1;
  
  dcz1[0]=y3-y2;
  dcz1[1]=x2-x3;
  dcz1[2]=0;
  
  dcz2[0]=y1-y3;
  dcz2[1]=x3-x1;
  dcz2[2]=0;
  
  dcz3[0]=y2-y1;
  dcz3[1]=x1-x2;
  dcz3[2]=0;
  
  dareax[0]=1/(2*normc)*(cv[0]*dcx1[0]+cv[1]*dcx1[1]+cv[2]*dcx1[2]);
  dareax[1]=1/(2*normc)*(cv[0]*dcx2[0]+cv[1]*dcx2[1]+cv[2]*dcx2[2]);
  dareax[2]=1/(2*normc)*(cv[0]*dcx3[0]+cv[1]*dcx3[1]+cv[2]*dcx3[2]);
  
  dareay[0]=1/(2*normc)*(cv[0]*dcy1[0]+cv[1]*dcy1[1]+cv[2]*dcy1[2]);
  dareay[1]=1/(2*normc)*(cv[0]*dcy2[0]+cv[1]*dcy2[1]+cv[2]*dcy2[2]);
  dareay[2]=1/(2*normc)*(cv[0]*dcy3[0]+cv[1]*dcy3[1]+cv[2]*dcy3[2]);
  
  dareaz[0]=1/(2*normc)*(cv[0]*dcz1[0]+cv[1]*dcz1[1]+cv[2]*dcz1[2]);
  dareaz[1]=1/(2*normc)*(cv[0]*dcz2[0]+cv[1]*dcz2[1]+cv[2]*dcz2[2]);
  dareaz[2]=1/(2*normc)*(cv[0]*dcz3[0]+cv[1]*dcz3[1]+cv[2]*dcz3[2]);
  
  for(int jk=0;jk<3;jk++)
  {
    dndx1[jk]=(normc*dcx1[jk]-cv[jk]*2*dareax[0])/(normc*normc);
    dndx2[jk]=(normc*dcx2[jk]-cv[jk]*2*dareax[1])/(normc*normc);
    dndx3[jk]=(normc*dcx3[jk]-cv[jk]*2*dareax[2])/(normc*normc);
    
    dndy1[jk]=(normc*dcy1[jk]-cv[jk]*2*dareay[0])/(normc*normc);
    dndy2[jk]=(normc*dcy2[jk]-cv[jk]*2*dareay[1])/(normc*normc);
    dndy3[jk]=(normc*dcy3[jk]-cv[jk]*2*dareay[2])/(normc*normc);
    
    dndz1[jk]=(normc*dcz1[jk]-cv[jk]*2*dareaz[0])/(normc*normc);
    dndz2[jk]=(normc*dcz2[jk]-cv[jk]*2*dareaz[1])/(normc*normc);
    dndz3[jk]=(normc*dcz3[jk]-cv[jk]*2*dareaz[2])/(normc*normc);
  }
}
void Calculate_Frame_Matrix(double *E,double *up,double R[3][3])
{
  //E is the vector pointing to the observer (unrotated)
  // up is the up vector, ie camera orientation
  
  //R is the output, 3x3 matrix, world frame -> Camera frame
 double x[3],y[3],z[3];
 double nx,ny,nz;
 nz=NORM(E);
 z[0]=E[0]/nz;
 z[1]=E[1]/nz;
 z[2]=E[2]/nz;
 cross(up,z,x);
 nx=NORM(x);
 x[0]=x[0]/nx;
 x[1]=x[1]/nx;
 x[2]=x[2]/nx;
 cross(z,x,y);
 ny=NORM(y);
 y[0]=y[0]/ny;
 y[1]=y[1]/ny;
 y[2]=y[2]/ny;
 Convert_to_Matrix(x,y,z,R);
}
  
