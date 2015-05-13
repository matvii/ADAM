#include"prepare.h"
void Calculate_HF(int *tlist,double *vlist,int nfac,int nvert,double *angles,double *Eo,double *E0o,double *up,double TIME,double dist,double Gamma,double A,double Hdist,int N,double WL,double *freqx,double *freqy,int nfreq,double *offset,double complex *F)
{
 double complex *F0;
 F0=(double complex*)calloc(nfreq,sizeof(double complex));
 double *Flux,*Fldx,*Fldy,*Fldz,*FldA;
 Flux=(double*)calloc(nfac,sizeof(double));
 double M[3][3],dMb[3][3],dMo[3][3],dMl[3][3],Mt[3][3];
 double R[3][3],Rdb[3][3],Rdl[3][3],Rdo[3][3],RT[3][3];
 double E[3],E0[3];
 double normalr[3],side1[3],side2[3];
 double dechdx[3],dechdy[3],dechdz[3],dechdA[3];
 double n[3],*nb,*cent;
 double *vb1,*vb2,*vb3;
 double vr1[3],vr2[3],vr3[3];
 double *v1,*v2,*v3;
 double scale;
 double complex tscale,FTC;
 double dp;
 double B,TB=0;
 double norm;
 int t1,t2,t3,blocker,sign;
int j1,j2,j3;
 double mu,mu0,area,mub,ech,rexp;
 double *normal,*centroid;
 int *visible;
 int tb1,tb2,tb3; //Indices to the vertices of possible blocker facet
 int blocked=0;
 //Distance km->arcsec
 dp=1/(dist*149597871.0)*180.0/PI*3600.0;
 visible=(int*)calloc(nfac,sizeof(int));
  //Allocate for memory
// normal=(double*)malloc(3*nfac*sizeof(double));
 // centroid=(double*)mxCalloc(3*nfac,sizeof(double));
 // IndexofBlocks=(int*)mxCalloc(nfac,sizeof(int));
// NumofBlocks=(int*)mxCalloc(nfac,sizeof(int));
  //Calculate frame change matrix
  Calculate_Frame_Matrix(Eo,up,R);
 
 //FacetsOverHorizon(tlist,vlist,nfac,nvert,normal,centroid,NumofBlocks,IndexofBlocks);
 
 rotate(angles[0],angles[1],angles[2],0.0,TIME,M,dMb,dMl,dMo);
 //Construct asteroid->Camera frame matrix, which is
 //asteroid->world frame->camera frame
 transpose(M,Mt); //Transpose, since we rotate the model, not view directions
 mult_mat(R,Mt,RT); 
 mult_vector(M,Eo,E);
 mult_vector(M,E0o,E0);

 /*For each facet,
  * 1)Check if facet is visible
  * 2) Calculate echo
  * 3) Convert triangle to range-Doppler frame
  * 4) Calculate FT
  */
//Find actual blockers
FindActualBlockers(tlist,vlist,nfac,nvert,E,E,1,visible);
Calculate_Radiance(tlist,vlist,nfac,nvert,angles,Eo,E0o,TIME,Gamma, A,Hdist,WL,N,Flux,Fldx,Fldy,Fldz,FldA,0);
//for(int i=27;i<nfac;i++)
 // mexPrintf("fl%d: %.10e\n",i+1, Flux[i]);

 //visible is nfac vector, visible[j]=1 if facet (j+1)th facet is visible
//NOTE INDEXING
//mexPrintf("%f %f %f\n",vlist[0],vlist[1],vlist[2]);
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
   
   v1=vlist+j1*3;
   v2=vlist+j2*3;
   v3=vlist+j3*3;
   
   
   //Calculate normals and centroids
   for(int i=0;i<3;i++)
   {
     
     side1[i]=dp*(v2[i]-v1[i]); //Convert km->arcsec
     side2[i]=dp*(v3[i]-v1[i]);
    
   }
   cross(side1,side2,n);
   norm=NORM(n);
   n[0]=n[0]/norm;
   n[1]=n[1]/norm;
   n[2]=n[2]/norm;
   
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
     Calc_FTC(freqx,freqy,nfreq,vr1[0],vr1[1],vr2[0],vr2[1],vr3[0],vr3[1],F0);
    // printf("Fdd: %f %f\n",creal(FTdd[0]),cimag(FTdd[0]));
     //Note that we sum to F at each round, does not work, we need to multiply with the echo
     //Derivatives wrt angles
   area=0.5*norm;
  //if(j==0)
   //{
    // mexPrintf("area: %f mu: %f F0: %f\n",area,mu,F0[0]);
   //}
   //mexPrintf("area: %f normal: %f %f %f mu: %f mu0: %f\n",area,n[0],n[1],n[2],mu,mu0); 
     
   B=Flux[j];
     for(int jf=0;jf<nfreq;jf++)
     {
       //This should be taken outside of the loop
       
    //   mexPrintf("scale:%f offset: %f %f\n",scale,creal(cexp(2*PI*I*(offset[0]*freqx[jf]+offset[1]*freqy[jf]))),cimag(cexp(2*PI*I*(offset[0]*freqx[jf]+offset[1]*freqy[jf]))));
       //FTC=tscale*F0[jf];
       F[jf]+=B*F0[jf];
      }
      TB=TB+B*area*mu;
//  printf("Flux: %f area: %f mu: %f\n",B,area,mu);

}
//printf("Total brightness: %f\n",TB);  
//Normalize with total brightness
for(int j=0;j<nfreq;j++)
  F[j]=cexp(2.0*PI*I*(offset[0]*freqx[j]+offset[1]*freqy[j]))*F[j]/TB;
free(visible);
free(Flux);
free(F0);

}

void Calculate_HF_deriv(int *tlist,double *vlist,int nfac,int nvert,double *angles,double *Eo,double *E0o,double *up,double TIME,double dist,double Gamma,double A,double Hdist,int N,double WL,double *freqx,double *freqy,int nfreq,double *offset,double complex *F,double complex *dFdx,double complex *dFdy,double complex *dFdz,double complex *dFdA,double complex *dFdoff)
{
 double complex *F0,*FTda,*FTdb,*FTdc,*FTdd,*FTdh,*FTdg,*FTdx,*FTdy,*FTdz,*FTdA;
 FTdx=(double complex*)calloc(nfreq*nvert,sizeof(double complex));
 FTdy=(double complex*)calloc(nfreq*nvert,sizeof(double complex));
 FTdz=(double complex*)calloc(nfreq*nvert,sizeof(double complex));
 FTdA=(double complex*)calloc(nfreq*3,sizeof(double complex));
 F0=(double complex*)calloc(nfreq,sizeof(double complex));
 FTda=(double complex*)malloc(nfreq*sizeof(double complex));
 FTdb=(double complex*)malloc(nfreq*sizeof(double complex));
 FTdc=(double complex*)malloc(nfreq*sizeof(double complex));
 FTdd=(double complex*)malloc(nfreq*sizeof(double complex));
 FTdh=(double complex*)malloc(nfreq*sizeof(double complex));
 FTdg=(double complex*)malloc(nfreq*sizeof(double complex));
 
 double M[3][3],dMb[3][3],dMo[3][3],dMl[3][3],Mt[3][3],dMbT[3][3],dMoT[3][3],dMlT[3][3]; //Rotation matrices and their derivatives and transposes
 double R[3][3],Rdb[3][3],Rdl[3][3],Rdo[3][3],RT[3][3]; //Projection matrix, and derivatives of combined projection+rotation matrix
 double dEdb[3],dEdl[3],dEdo[3],dE0db[3],dE0dl[3],dE0do[3];
 double dndx1[3],dndx2[3],dndx3[3],dndy1[3],dndy2[3],dndy3[3],dndz1[3],dndz2[3],dndz3[3]; //Derivatives of the facet normal vector
 double dBdx1,dBdy1,dBdz1,dBdx2,dBdy2,dBdz2,dBdx3,dBdy3,dBdz3,dBdb,dBdl,dBdo; //Derivatives of facet brightness
 double *dTBdx,*dTBdy,*dTBdz; //Derivatives of total brightness, allocating memory
 double *Flux,*Fldx,*Fldy,*Fldz,*FldA;
 Flux=(double*)calloc(nfac,sizeof(double));
 Fldx=(double*)calloc(nfac*nvert,sizeof(double));
 Fldy=(double*)calloc(nfac*nvert,sizeof(double));
 Fldz=(double*)calloc(nfac*nvert,sizeof(double));
 FldA=(double*)calloc(nfac*3,sizeof(double));
 double dTBdA[3]={0};
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
 double v1[3],v2[3],v3[3];
 double complex scale;

 double dp;
 double B,TB=0.0;
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
  Calculate_Frame_Matrix(Eo,up,R);
  
  
  //Calculate frame change matrix
 
 
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
 dadx=RT[0][0];
  dady=RT[0][1];
  dadz=RT[0][2];
  dbdx=RT[1][0];
  dbdy=RT[1][1];
  dbdz=RT[1][2];
 /*For each facet,
  * 1)Check if facet is visible
  * 2) Calculate echo
  * 3) Convert triangle to range-Doppler frame
  * 4) Calculate FT
  */
//Find actual blockers
FindActualBlockers(tlist,vlist,nfac,nvert,E,E,1,visible);
 //visible is nfac vector, visible[j]=1 if facet (j+1)th facet is visible
//NOTE INDEXING
Calculate_Radiance(tlist,vlist,nfac,nvert,angles,Eo,E0o,TIME,Gamma, A,Hdist,WL,N,Flux,Fldx,Fldy,Fldz,FldA,1);

 //for(int j=0;j<nfac;j++)
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
    Calculate_Area_and_Normal_Derivative(v1,v2,v3,n,dndx1,dndx2,dndx3,dndy1,dndy2,dndy3,dndz1,dndz2,dndz3,&area,dAdx,dAdy,dAdz);
    
   //Calculate normals and centroids
   
   mu=DOT(E,n);
   
  //Convert to camera frame
   mult_vector(RT,v1,vr1);
   mult_vector(RT,v2,vr2);
   mult_vector(RT,v3,vr3);
  
    
     //Now we should convert to frequency domain, ie calculate the contribution of each facet
  //   Calc_FTC(freqx,freqy,nfreq,vr1[0],vr1[1],vr2[0],vr2[1],vr3[0],vr3[1],F0);
     Calc_FTC_deriv(freqx,freqy,nfreq,vr1[0],vr1[1],vr2[0],vr2[1],vr3[0],vr3[1],F0,FTda,FTdb,FTdc,FTdd,FTdg,FTdh);
  
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
     B=Flux[j];
   //Derivatives of B
  
   dBdx1=Fldx[j*nvert+j1]/dp;
   dBdx2=Fldx[j*nvert+j2]/dp;
   dBdx3=Fldx[j*nvert+j3]/dp;
   
   dBdy1=Fldy[j*nvert+j1]/dp;
   dBdy2=Fldy[j*nvert+j2]/dp;
   dBdy3=Fldy[j*nvert+j3]/dp;
   
   dBdz1=Fldz[j*nvert+j1]/dp;
   dBdz2=Fldz[j*nvert+j2]/dp;
   dBdz3=Fldz[j*nvert+j3]/dp;
   dBdb=FldA[3*j];
   dBdl=FldA[3*j+1];
   dBdo=FldA[3*j+2];
   //Derivative of total brightness
   dTBdx[j1]+=dBdx1*area*mu+B*dAdx[0]*mu+B*area*dmudx1;
   dTBdx[j2]+=dBdx2*area*mu+B*dAdx[1]*mu+B*area*dmudx2;
   dTBdx[j3]+=dBdx3*area*mu+B*dAdx[2]*mu+B*area*dmudx3;
   dTBdy[j1]+=dBdy1*area*mu+B*dAdy[0]*mu+B*area*dmudy1;
   dTBdy[j2]+=dBdy2*area*mu+B*dAdy[1]*mu+B*area*dmudy2;
   dTBdy[j3]+=dBdy3*area*mu+B*dAdy[2]*mu+B*area*dmudy3;
   dTBdz[j1]+=dBdz1*area*mu+B*dAdz[0]*mu+B*area*dmudz1;
   dTBdz[j2]+=dBdz2*area*mu+B*dAdz[1]*mu+B*area*dmudz2;
   dTBdz[j3]+=dBdz3*area*mu+B*dAdz[2]*mu+B*area*dmudz3;
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
    dFdx[j*nvert+k]=dp*scale*(FTdx[j*nvert+k]*TB-F[j]*dTBdx[k])/pow(TB,2);
    dFdy[j*nvert+k]=dp*scale*(FTdy[j*nvert+k]*TB-F[j]*dTBdy[k])/pow(TB,2);
    dFdz[j*nvert+k]=dp*scale*(FTdz[j*nvert+k]*TB-F[j]*dTBdz[k])/pow(TB,2);
  }
  dFdA[j*3+0]=scale*(FTdA[j*3+0]*TB-F[j]*dTBdA[0])/pow(TB,2);
  dFdA[j*3+1]=scale*(FTdA[j*3+1]*TB-F[j]*dTBdA[1])/pow(TB,2);
  dFdA[j*3+2]=scale*(FTdA[j*3+2]*TB-F[j]*dTBdA[2])/pow(TB,2);
  
  F[j]=cexp(2.0*PI*I*(offset[0]*freqx[j]+offset[1]*freqy[j]))*F[j]/TB;
  dFdoff[j*2+0]=2.0*PI*I*freqx[j]*F[j];
  dFdoff[j*2+1]=2.0*PI*I*freqy[j]*F[j];
}


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
free(Flux);
free(Fldx);
free(Fldy);
free(Fldz);
free(FldA);
}