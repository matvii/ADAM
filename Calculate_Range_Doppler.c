#include"utils.h"
#include"matrix_ops.h"


#define INI_MAX_RD_ANGLE 60
void multmat(double A[3][3],double B[3][3],double C[3][3])
{
  for(int i=0;i<3;i++)
  {
    for(int j=0;j<3;j++)
    {
  C[i][j]=A[i][0]*B[0][j]+A[i][1]*B[1][j]+A[i][2]*B[2][j];
    }
  }
} 


void multmat3t(double A[3][3],double B[3][3],double C[3][3],double D[3][3])
{
  //Third matrix is transposed
  double T[3][3],Ct[3][3];
  multmat(A,B,T);
  transpose(C,Ct);
  multmat(T,Ct,D);
}
void Calculate_Frame_Matrix_Derivatives(double *E,double *angles,double TIME,double freq,double R[3][3],double Rdb[3][3],double Rdl[3][3],double Rdo[3][3]);
void Calculate_Frame_Matrix2(double *E,double *angles,double TIME,double freq,double R[3][3]);
void Calculate_Normal_Derivative(double *v1,double *v2,double *v3,double *n1dx,double *n2dx,double *n3dx,double *n1dy,double *n2dy,double *n3dy,double *n1dz,double *n2dz,double *n3dz);

void Calculate_Range_Doppler_deriv(int *tlist,double *vlist,int nfac,int nvert,double *angles,double *Eo,double TIME,double *freqx,double *freqy,int nfreq,double rfreq,double *offset,double scal,double rexpe,double* Fr,double *Fi,double *FTdxrf,double *FTdxif,double *FTdyrf,double *FTdyif,double *FTdzrf,double *FTdzif,double *FTdArf,double *FTdAif,double  *FTdoffr,double *FTdoffi,double * FTdexpr,double *FTdexpi)
{
//Map triangle to the Range-Doppler frame
 double complex *F0,*FTda,*FTdb,*FTdc,*FTdd,*FTdh,*FTdg;
 double complex temp;
// double complex *dFda,*dFdb,*dFdc,*dFdd,*dFdh,*dFdg;
 F0=calloc(nfreq,sizeof(double complex));
 FTda=calloc(nfreq,sizeof(double complex));
 FTdb=calloc(nfreq,sizeof(double complex));
 FTdc=calloc(nfreq,sizeof(double complex));
 FTdd=calloc(nfreq,sizeof(double complex));
 FTdh=calloc(nfreq,sizeof(double complex));
 FTdg=calloc(nfreq,sizeof(double complex));
 double dmudl,dmudb,dmudo;
 double complex fscale;
 double *FTdxr=calloc(nfreq*nvert,sizeof(double));
 double *FTdxi=calloc(nfreq*nvert,sizeof(double));
 double *FTdyr=calloc(nfreq*nvert,sizeof(double));
 double *FTdyi=calloc(nfreq*nvert,sizeof(double));
 double *FTdzr=calloc(nfreq*nvert,sizeof(double));
 double *FTdzi=calloc(nfreq*nvert,sizeof(double));
 double *FTdAr=calloc(nfreq*3,sizeof(double));
 double *FTdAi=calloc(nfreq*3,sizeof(double));
 int j1,j2,j3;
 double max_cos_angle=cos(INI_MAX_RD_ANGLE*PI/180);
 double dAdx[3],dAdy[3],dAdz[3];
 double M[3][3],dMb[3][3],dMo[3][3],dMl[3][3];
 double R[3][3],Rdb[3][3],Rdl[3][3],Rdo[3][3];
 double E[3],dEdb[3],dEdl[3],dEdo[3];
 double normalr[3],side1[3],side2[3];
 double dechdx[3],dechdy[3],dechdz[3],dechdA[3];
 double *n,*nb,*cent;
 double *vb1,*vb2,*vb3;
 double vr1[3],vr2[3],vr3[3];
 double *v1,*v2,*v3;
 double scale;
 double complex tscale,FTC;
 double TB=0,area=0;
 double dadx,dady,dadz,dbdx,dbdy,dbdz,dcdx,dcdy,dcdz;
 scale=exp(scal);
 int t1,t2,t3,blocker,sign;
 double mu,mub,ech,rexp;
 double *normal,*centroid;
 int *IndexofBlocks,*NumofBlocks;
 int tb1,tb2,tb3; //Indices to the vertices of possible blocker facet
 int blocked=0;
 //Allocate for derivatives
 double dech=0.0; //placeholder of derivative of echo
 double dndx1[3],dndx2[3],dndx3[3],dndy1[3],dndy2[3],dndy3[3],dndz1[3],dndz2[3],dndz3[3];
double *dTBdx,*dTBdy,*dTBdz; //Derivatives of total brightness, allocating memory
 double dTBdA[3]={0.0,0.0,0.0};
 dTBdx=calloc(nvert,sizeof(double));
 dTBdy=calloc(nvert,sizeof(double));
 dTBdz=calloc(nvert,sizeof(double));
 double v1db[3],v1dl[3],v1do[3],v2db[3],v2dl[3],v2do[3],v3db[3],v3dl[3],v3do[3];
 double dmudx1,dmudy1,dmudz1,dmudx2,dmudy2,dmudz2,dmudx3,dmudy3,dmudz3;
  //Allocate for memory
  normal=calloc(3*nfac,sizeof(double));
  centroid=calloc(3*nfac,sizeof(double));
  IndexofBlocks=calloc(nfac*nfac,sizeof(int));
  NumofBlocks=calloc(nfac,sizeof(int));
  //Calculate frame change matrix
  Calculate_Frame_Matrix_Derivatives(Eo,angles,TIME,rfreq,R,Rdb,Rdl,Rdo);
  dadx=R[0][0];
  dady=R[0][1];
  dadz=R[0][2];
  dbdx=R[1][0];
  dbdy=R[1][1];
  dbdz=R[1][2];
  dcdx=R[2][0];
  dcdy=R[2][1];
  dcdz=R[2][2];
 
  //Ok, R is the conversion matrix, rotation angles included, so use Eo
  //Find possible blockers
 FacetsOverHorizon(tlist,vlist,nfac,nvert,normal,centroid,NumofBlocks,IndexofBlocks);
 
 rotate(angles[0],angles[1],angles[2],angles[3],TIME,M,dMb,dMl,dMo);
 
 mult_vector(M,Eo,E);
 //Derivatives of E wrt beta,lambda,omega
 mult_vector(dMb,Eo,dEdb);
 mult_vector(dMl,Eo,dEdl);
 mult_vector(dMo,Eo,dEdo);

 /*For each facet,
  * 1)Check if facet is visible
  * 2) Calculate echo
  * 3) Convert triangle to range-Doppler frame
  * 4) Calculate FT
  */
 
 //printf("nfac: %d\n",nfac);
// printf("offset: %f %f\n",offset[0],offset[1]);
 for(int j=0;j<nfac;j++)
 {
   
   n=&normal[3*j];
   mu=DOT(E,n); //Check if facet is even pointed towards the radar
  j1=tlist[j*3]-1;
   j2=tlist[j*3+1]-1;
   j3=tlist[j*3+2]-1;
   if(mu<EP)
     continue;
   if(mu<max_cos_angle)
       continue;
   //Test for blocking facets
   //Facet centroid
   
   cent=&centroid[3*j]; //j starts from 0
   blocked=0;
   
     for(int k=0;k<NumofBlocks[j];k++)
     {
       blocked=0;
       blocker=IndexofBlocks[nfac*j+k]; //This indexing starts from 1
       //Index to tlist, blocker facet
       nb=&normal[3*(blocker-1)];
       mub=DOT(E,nb);
       tb1=tlist[3*(blocker-1)]; //Indexing starts from 1
       tb2=tlist[3*(blocker-1)+1];
       tb3=tlist[3*(blocker-1)+2];
       
       //Blocker vertices
       vb1=&vlist[3*(tb1-1)];
       vb2=&vlist[3*(tb2-1)];
       vb3=&vlist[3*(tb3-1)];
       
       //Possible blocker should not be visible to the radar
       //Check this condition.
    
       if(mub<EP)
	 if(is_in_triangle(cent,E,vb1,vb2,vb3)==1)
	 {
	   blocked=1;
	   break;
	 }
	 
	
     }
    
     if(blocked==1) //Is this point in the current facet blocked?
       continue; //Blocked, skip it
     //Calculate normal derivatives
     
     //Ok, so this facet is visible. Calculate echo
     rexp=exp(rexpe)-1;
     //mu=DOT(E,n);
     ech=pow(mu,rexp);
    // printf("ech: %f\n",ech);
    
     t1=tlist[3*(j)]; //Indexing starts from 0, t1 starts from one
     t2=tlist[3*(j)+1];
     t3=tlist[3*(j)+2];
     //Now comes the hard part, we transfer to the radar frame
     v1=&(vlist[3*(t1-1)]);
     v2=&vlist[3*(t2-1)];
     v3=&vlist[3*(t3-1)];
     //Calculate Normal derivatives (in the original frame)
     //Calculate_Normal_Derivative(v1,v2,v3,dndx1,dndx2,dndx3,dndy1,dndy2,dndy3,dndz1,dndz2,dndz3);
     Calculate_Area_and_Normal_Derivative(v1,v2,v3,n,dndx1,dndx2,dndx3,dndy1,dndy2,dndy3,dndz1,dndz2,dndz3,&area,dAdx,dAdy,dAdz);
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
     //Calculate echo derivative
     dech=ech*log(mu)*rexp;
     dechdx[0]=rexp*pow(mu,rexp-1)*(DOT(E,dndx1));
     dechdx[1]=rexp*pow(mu,rexp-1)*(DOT(E,dndx2));
     dechdx[2]=rexp*pow(mu,rexp-1)*(DOT(E,dndx3));
     dechdy[0]=rexp*pow(mu,rexp-1)*(DOT(E,dndy1));
     dechdy[1]=rexp*pow(mu,rexp-1)*(DOT(E,dndy2));
     dechdy[2]=rexp*pow(mu,rexp-1)*(DOT(E,dndy3));
     dechdz[0]=rexp*pow(mu,rexp-1)*(DOT(E,dndz1));
     dechdz[1]=rexp*pow(mu,rexp-1)*(DOT(E,dndz2));
     dechdz[2]=rexp*pow(mu,rexp-1)*(DOT(E,dndz3));
     dechdA[0]=rexp*pow(mu,rexp-1)*(DOT(dEdb,n));
     dechdA[1]=rexp*pow(mu,rexp-1)*(DOT(dEdl,n));
     dechdA[2]=rexp*pow(mu,rexp-1)*(DOT(dEdo,n));
     
     mult_vector(R,v1,vr1);
     
     mult_vector(R,v2,vr2);
     mult_vector(R,v3,vr3);
     //Check if sign change is needed:
     for(int i=0;i<3;i++)
     {
      side1[i]=vr2[i]-vr1[i];
      side2[i]=vr3[i]-vr1[i];
     }
     cross(side1,side2,normalr);
     if(normalr[2]<0)
       sign=-1;
     else
       sign=1;
    Calculate_Area_and_Normal_Derivative(vr1,vr2,vr3,n,dndx1,dndx2,dndx3,dndy1,dndy2,dndy3,dndz1,dndz2,dndz3,&area,dAdx,dAdy,dAdz);
     //x coordinate is the delay, y frequency
     //Now we should convert to frequency domain, ie calculate the contribution of each facet
     Calc_FTC_deriv(freqx,freqy,nfreq,vr1[0],vr1[1],vr2[0],vr2[1],vr3[0],vr3[1],F0,FTda,FTdb,FTdc,FTdd,FTdg,FTdh);
    // printf("Fdd: %f %f\n",creal(FTdd[0]),cimag(FTdd[0]));
     //Note that we sum to F at each round, does not work, we need to multiply with the echo
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
    // printf("dechdA: %f %f %f\n",dechdA[0],dechdA[1],dechdA[2]);
     //printf("dEdb: %f %f %f\n",dEdb[0],dEdb[1],dEdb[2]);
     //REMEMBER SCALING FOR FTda...
    // printf("v1d0:\n");
    // printf("%f %f %f %f %f %f\n",v1do[0],v1do[1],v2do[0],v2do[1],v3do[0],v3do[1]);
    // printf("Matrix:\n");
     //for(int db=0;db<3;db++)
   //printf("%f %f %f\n",R[db][0],R[db][1],R[db][2]);
     
     for(int jf=0;jf<nfreq;jf++)
     {
       tscale=scale*cexp(2*PI*I*(offset[0]*freqx[jf]+offset[1]*freqy[jf]));
       
       FTC=tscale*F0[jf];
       //F[jf]+=sign*ech*FTC;
       Fr[jf]+=creal(sign*ech*FTC);
       Fi[jf]+=cimag(sign*ech*FTC);
       temp=tscale*sign*(F0[jf]*dechdx[0]+ech*(FTda[jf]*dadx+FTdb[jf]*dbdx));
      // FTdx[jf*nvert+t1-1]+=tscale*sign*(F0[jf]*dechdx[0]+ech*(FTda[jf]*dadx+FTdb[jf]*dbdx));
       FTdxr[jf*nvert+t1-1]+=creal(temp);
       FTdxi[jf*nvert+t1-1]+=cimag(temp);
       
       temp=tscale*sign*(F0[jf]*dechdx[1]+ech*(FTdc[jf]*dadx+FTdd[jf]*dbdx));
       //FTdx[jf*nvert+t2-1]+=tscale*sign*(F0[jf]*dechdx[1]+ech*(FTdc[jf]*dadx+FTdd[jf]*dbdx));
       FTdxr[jf*nvert+t2-1]+=creal(temp);
       FTdxi[jf*nvert+t2-1]+=cimag(temp);
       
       temp=tscale*sign*(F0[jf]*dechdx[2]+ech*(FTdg[jf]*dadx+FTdh[jf]*dbdx));
       //FTdx[jf*nvert+t3-1]+=tscale*sign*(F0[jf]*dechdx[2]+ech*(FTdg[jf]*dadx+FTdh[jf]*dbdx));
       FTdxr[jf*nvert+t3-1]+=creal(temp);
       FTdxi[jf*nvert+t3-1]+=cimag(temp);
       
       temp=tscale*sign*(F0[jf]*dechdy[0]+ech*(FTda[jf]*dady+FTdb[jf]*dbdy));
       //FTdy[jf*nvert+t1-1]+=tscale*sign*(F0[jf]*dechdy[0]+ech*(FTda[jf]*dady+FTdb[jf]*dbdy));
       FTdyr[jf*nvert+t1-1]+=creal(temp);
       FTdyi[jf*nvert+t1-1]+=cimag(temp);
       
       temp=tscale*sign*(F0[jf]*dechdy[1]+ech*(FTdc[jf]*dady+FTdd[jf]*dbdy));
       //FTdy[jf*nvert+t2-1]+=tscale*sign*(F0[jf]*dechdy[1]+ech*(FTdc[jf]*dady+FTdd[jf]*dbdy));
       FTdyr[jf*nvert+t2-1]+=creal(temp);
       FTdyi[jf*nvert+t2-1]+=cimag(temp);
       
       temp=tscale*sign*(F0[jf]*dechdy[2]+ech*(FTdg[jf]*dady+FTdh[jf]*dbdy));
       //FTdy[jf*nvert+t3-1]+=tscale*sign*(F0[jf]*dechdy[2]+ech*(FTdg[jf]*dady+FTdh[jf]*dbdy));
       FTdyr[jf*nvert+t3-1]+=creal(temp);
       FTdyi[jf*nvert+t3-1]+=cimag(temp);
       
       temp=tscale*sign*(F0[jf]*dechdz[0]+ech*(FTda[jf]*dadz+FTdb[jf]*dbdz));
       //FTdz[jf*nvert+t1-1]+=tscale*sign*(F0[jf]*dechdz[0]+ech*(FTda[jf]*dadz+FTdb[jf]*dbdz));
       FTdzr[jf*nvert+t1-1]+=creal(temp);
       FTdzi[jf*nvert+t1-1]+=cimag(temp);
       
       temp=tscale*sign*(F0[jf]*dechdz[1]+ech*(FTdc[jf]*dadz+FTdd[jf]*dbdz));
       //FTdz[jf*nvert+t2-1]+=tscale*sign*(F0[jf]*dechdz[1]+ech*(FTdc[jf]*dadz+FTdd[jf]*dbdz));
       FTdzr[jf*nvert+t2-1]+=creal(temp);
       FTdzi[jf*nvert+t2-1]+=cimag(temp);
       temp=tscale*sign*(F0[jf]*dechdz[2]+ech*(FTdg[jf]*dadz+FTdh[jf]*dbdz));
       //FTdz[jf*nvert+t3-1]+=tscale*sign*(F0[jf]*dechdz[2]+ech*(FTdg[jf]*dadz+FTdh[jf]*dbdz));
       FTdzr[jf*nvert+t3-1]+=creal(temp);
       FTdzi[jf*nvert+t3-1]+=cimag(temp);
       
       //angle derivatives
       temp=tscale*sign*(F0[jf]*dechdA[0]+ech*(FTda[jf]*v1db[0]+FTdb[jf]*v1db[1]+FTdc[jf]*v2db[0]+FTdd[jf]*v2db[1]+FTdg[jf]*v3db[0]+FTdh[jf]*v3db[1]));
       //FTdA[jf*3+0]+=tscale*sign*(F0[jf]*dechdA[0]+ech*(FTda[jf]*v1db[0]+FTdb[jf]*v1db[1]+FTdc[jf]*v2db[0]+FTdd[jf]*v2db[1]+FTdg[jf]*v3db[0]+FTdh[jf]*v3db[1]));
       FTdAr[jf*3+0]+=creal(temp);
       FTdAi[jf*3+0]+=cimag(temp);
       temp=tscale*sign*(F0[jf]*dechdA[1]+ech*(FTda[jf]*v1dl[0]+FTdb[jf]*v1dl[1]+FTdc[jf]*v2dl[0]+FTdd[jf]*v2dl[1]+FTdg[jf]*v3dl[0]+FTdh[jf]*v3dl[1]));
       //FTdA[jf*3+1]+=tscale*sign*(F0[jf]*dechdA[1]+ech*(FTda[jf]*v1dl[0]+FTdb[jf]*v1dl[1]+FTdc[jf]*v2dl[0]+FTdd[jf]*v2dl[1]+FTdg[jf]*v3dl[0]+FTdh[jf]*v3dl[1]));
       FTdAr[jf*3+1]+=creal(temp);
       FTdAi[jf*3+1]+=cimag(temp);
       temp=tscale*sign*(F0[jf]*dechdA[2]+ech*(FTda[jf]*v1do[0]+FTdb[jf]*v1do[1]+FTdc[jf]*v2do[0]+FTdd[jf]*v2do[1]+FTdg[jf]*v3do[0]+FTdh[jf]*v3do[1]));
       //FTdA[jf*3+2]+=tscale*sign*(F0[jf]*dechdA[2]+ech*(FTda[jf]*v1do[0]+FTdb[jf]*v1do[1]+FTdc[jf]*v2do[0]+FTdd[jf]*v2do[1]+FTdg[jf]*v3do[0]+FTdh[jf]*v3do[1]));
       FTdAr[jf*3+2]+=creal(temp);
       FTdAi[jf*3+2]+=cimag(temp);
      // printf("FTdA:%f %f\n",creal(FTdd[j]*v2do[1]),cimag(FTdd[j]*v2do[1]));
       //FTdexp[jf]+=sign*dech*FTC;
       FTdexpr[jf]+=creal(sign*dech*FTC);
       FTdexpi[jf]+=cimag(sign*dech*FTC);
       //FTdoff[jf*2+0]+=sign*(2*PI*ech*I*freqx[jf]*FTC);
       //FTdoffr[jf*2+0]+=creal(sign*(2*PI*ech*I*freqx[jf]*FTC));
       //FTdoffi[jf*2+0]+=cimag(sign*(2*PI*ech*I*freqx[jf]*FTC));
       //FTdoff[jf*2+1]+=sign*(2*PI*ech*I*freqy[jf]*FTC);
       //FTdoffr[jf*2+1]+=creal(sign*(2*PI*ech*I*freqy[jf]*FTC));
      // FTdoffi[jf*2+1]+=cimag(sign*(2*PI*ech*I*freqy[jf]*FTC));
     
  //printf("%f %f %f %f %f %f\n",vr1[0],vr1[1],vr2[0],vr2[1],vr3[0],vr3[1]);
       
 }
// printf("j:%d tscale: %f %f scale: %f sign: %d ech: %f\n",j,creal(tscale),cimag(tscale),scale,sign,ech);
 TB+=mu*ech*area;
 dTBdx[j1]+=dmudx1*ech*area+mu*dechdx[0]*area+mu*ech*(dAdx[0]*dadx+dAdy[0]*dbdx+dAdz[0]*dcdx);
 dTBdx[j2]+=dmudx2*ech*area+mu*dechdx[1]*area+mu*ech*(dAdx[1]*dadx+dAdy[1]*dbdx+dAdz[1]*dcdx);
 dTBdx[j3]+=dmudx3*ech*area+mu*dechdx[2]*area+mu*ech*(dAdx[2]*dadx+dAdy[2]*dbdx+dAdz[2]*dcdx);
 dTBdy[j1]+=dmudy1*ech*area+mu*dechdy[0]*area+mu*ech*(dAdx[0]*dady+dAdy[0]*dbdy+dAdz[0]*dcdy);
 dTBdy[j2]+=dmudy2*ech*area+mu*dechdy[1]*area+mu*ech*(dAdx[1]*dady+dAdy[1]*dbdy+dAdz[1]*dcdy);
 dTBdy[j3]+=dmudy3*ech*area+mu*dechdy[2]*area+mu*ech*(dAdx[2]*dady+dAdy[2]*dbdy+dAdz[2]*dcdy);
 dTBdz[j1]+=dmudz1*ech*area+mu*dechdz[0]*area+mu*ech*(dAdx[0]*dadz+dAdy[0]*dbdz+dAdz[0]*dcdz);
 dTBdz[j2]+=dmudz2*ech*area+mu*dechdz[1]*area+mu*ech*(dAdx[1]*dadz+dAdy[1]*dbdz+dAdz[1]*dcdz);
 dTBdz[j3]+=dmudz3*ech*area+mu*dechdz[2]*area+mu*ech*(dAdx[2]*dadz+dAdy[2]*dbdz+dAdz[2]*dcdz);
 
 dTBdA[0]+=dmudb*ech*area+mu*dechdA[0]*area;
 dTBdA[1]+=dmudl*ech*area+mu*dechdA[1]*area;
 dTBdA[2]+=dmudo*ech*area+mu*dechdA[2]*area;
}
// printf("TB: %7.8f\n",TB);
// write_matrix_file("/tmp/dTBdA.txt",dTBdA,1,3);
// write_matrix_file("/tmp/dTBdx.txt",dTBdx,1,nvert);
// write_matrix_file("/tmp/dTBdy.txt",dTBdy,1,nvert);
// write_matrix_file("/tmp/dTBdz.txt",dTBdz,1,nvert);
double TB2=pow(TB,2);

double tempr,tempi;
double fscalei,fscaler;
for(int j=0;j<nfreq;j++)
{
   // fscale=scale*cexp(2.0*PI*I*(offset[0]*freqx[j]+offset[1]*freqy[j]));
  //  fscaler=creal(fscale);
   // fscalei=cimag(fscale);
    for(int k=0;k<nvert;k++)
    {
        FTdxrf[j*nvert+k]=(FTdxr[j*nvert+k]*TB-Fr[j]*dTBdx[k])/TB2;
        FTdxif[j*nvert+k]=(FTdxi[j*nvert+k]*TB-Fi[j]*dTBdx[k])/TB2;
        FTdyrf[j*nvert+k]=(FTdyr[j*nvert+k]*TB-Fr[j]*dTBdy[k])/TB2;
        FTdyif[j*nvert+k]=(FTdyi[j*nvert+k]*TB-Fi[j]*dTBdy[k])/TB2;
        FTdzrf[j*nvert+k]=(FTdzr[j*nvert+k]*TB-Fr[j]*dTBdz[k])/TB2;
        FTdzif[j*nvert+k]=(FTdzi[j*nvert+k]*TB-Fi[j]*dTBdz[k])/TB2;
    }
    FTdArf[j*3+0]=(FTdAr[j*3+0]*TB-Fr[j]*dTBdA[0])/TB2;
    FTdAif[j*3+0]=(FTdAi[j*3+0]*TB-Fi[j]*dTBdA[0])/TB2; 
    FTdArf[j*3+1]=(FTdAr[j*3+1]*TB-Fr[j]*dTBdA[1])/TB2;
    FTdAif[j*3+1]=(FTdAi[j*3+1]*TB-Fi[j]*dTBdA[1])/TB2;     
    FTdArf[j*3+2]=(FTdAr[j*3+2]*TB-Fr[j]*dTBdA[2])/TB2;
    FTdAif[j*3+2]=(FTdAi[j*3+2]*TB-Fi[j]*dTBdA[2])/TB2;
    temp=(Fr[j]+I*Fi[j])/TB;
    Fr[j]=Fr[j]/TB;
    Fi[j]=Fi[j]/TB;
    FTdoffr[j*2+0]=creal(temp*2*PI*I*freqx[j]);
    FTdoffi[j*2+0]=cimag(temp*2*PI*I*freqx[j]);
    FTdoffr[j*2+1]=creal(temp*2*PI*I*freqy[j]);
    FTdoffi[j*2+1]=cimag(temp*2*PI*I*freqy[j]);
}
//   write_matrix_file("/tmp/FTr1.txt",Fr,1,nfreq);
//   write_matrix_file("/tmp/FTi1.txt",Fi,1,nfreq);  
  free(normal);
  free(centroid);
  free(IndexofBlocks);
  free(NumofBlocks);
  free(F0);
  free(FTda);
  free(FTdb);
  free(FTdc);
  free(FTdd);
  free(FTdg);
  free(FTdh);
  free(FTdxr);
  free(FTdxi);
  free(FTdyr);
  free(FTdyi);
  free(FTdzr);
  free(FTdzi);
  free(FTdAr);
  free(FTdAi);
  free(dTBdx);
  free(dTBdy);
  free(dTBdz);
 
}
void Calculate_Range_Doppler(int *tlist,double *vlist,int nfac,int nvert,double *angles,double *Eo,double TIME,double *freqx,double *freqy,int nfreq,double rfreq,double *offset,double scal,double rexpe,double * Fr,double *Fi)
{
//Map triangle to the Range-Doppler frame
 double complex *F0;
// double complex *dFda,*dFdb,*dFdc,*dFdd,*dFdh,*dFdg;
 F0=calloc(nfreq,sizeof(double complex));
 double max_cos_angle=cos(INI_MAX_RD_ANGLE*PI/180);
 double M[3][3],dMb[3][3],dMo[3][3],dMl[3][3];
 double R[3][3],Rdb[3][3],Rdl[3][3],Rdo[3][3];
 double E[3],dEdb[3],dEdl[3],dEdo[3];
 double normalr[3],side1[3],side2[3];
 double dechdx[3],dechdy[3],dechdz[3],dechdA[3];
 double *n,*nb,*cent;
 double *vb1,*vb2,*vb3;
 double vr1[3],vr2[3],vr3[3];
 double *v1,*v2,*v3;
 double complex scale,tscale;
 double complex fscale,FTC;
 
 scale=exp(scal);
 int t1,t2,t3,blocker,sign;
 double mu,mub,ech,rexp;
 double *normal,*centroid;
 int *IndexofBlocks,*NumofBlocks;
 int tb1,tb2,tb3; //Indices to the vertices of possible blocker facet
 int blocked=0;
 double TB=0,area;
 //Allocate for derivatives
 double dech=0.0; //placeholder of derivative of echo
 double dndx1[3],dndx2[3],dndx3[3],dndy1[3],dndy2[3],dndy3[3],dndz1[3],dndz2[3],dndz3[3];

 double v1db[3],v1dl[3],v1do[3],v2db[3],v2dl[3],v2do[3],v3db[3],v3dl[3],v3do[3];
  //Allocate for memory
  normal=calloc(3*nfac,sizeof(double));
  centroid=calloc(3*nfac,sizeof(double));
  IndexofBlocks=calloc(nfac*nfac,sizeof(int));
  NumofBlocks=calloc(nfac,sizeof(int));
  //Calculate frame change matrix
  Calculate_Frame_Matrix_Derivatives(Eo,angles,TIME,rfreq,R,Rdb,Rdl,Rdo);
 
 FacetsOverHorizon(tlist,vlist,nfac,nvert,normal,centroid,NumofBlocks,IndexofBlocks);
 
 rotate(angles[0],angles[1],angles[2],angles[3],TIME,M,dMb,dMl,dMo);
 
 mult_vector(M,Eo,E);
 

 /*For each facet,
  * 1)Check if facet is visible
  * 2) Calculate echo
  * 3) Convert triangle to range-Doppler frame
  * 4) Calculate FT
  */
 
 //printf("nfac: %d\n",nfac);
// printf("offset: %f %f\n",offset[0],offset[1]);
 for(int j=0;j<nfac;j++)
 {
   
   n=&normal[3*j];
   mu=DOT(E,n); //Check if facet is even pointed towards the radar
  
   if(mu<EP)
     continue;
   if(mu<max_cos_angle)
       continue;
   //Test for blocking facets
   //Facet centroid
   cent=&centroid[3*j]; //j starts from 0
   blocked=0;
   
     for(int k=0;k<NumofBlocks[j];k++)
     {
       blocked=0;
       blocker=IndexofBlocks[nfac*j+k]; //This indexing starts from 1
       //Index to tlist, blocker facet
       nb=&normal[3*(blocker-1)];
       mub=DOT(E,nb);
       tb1=tlist[3*(blocker-1)]; //Indexing starts from 1
       tb2=tlist[3*(blocker-1)+1];
       tb3=tlist[3*(blocker-1)+2];
       
       //Blocker vertices
       vb1=&vlist[3*(tb1-1)];
       vb2=&vlist[3*(tb2-1)];
       vb3=&vlist[3*(tb3-1)];
       
       //Possible blocker should not be visible to the radar
       //Check this condition.
    
       if(mub<EP)
	 if(is_in_triangle(cent,E,vb1,vb2,vb3)==1)
	 {
	   blocked=1;
	   break;
	 }
	 
	
     }
    
     if(blocked==1) //Is this point in the current facet blocked?
       continue; //Blocked, skip it
     //Calculate normal derivatives
     
     //Ok, so this facet is visible. Calculate echo
     rexp=exp(rexpe)-1;
     //mu=DOT(E,n);
     ech=pow(mu,rexp);
    // printf("ech: %f\n",ech);
    
     t1=tlist[3*(j)]; //Indexing starts from 0, t1 starts from one
     t2=tlist[3*(j)+1];
     t3=tlist[3*(j)+2];
     //Now comes the hard part, we transfer to the radar frame
     v1=&(vlist[3*(t1-1)]);
     v2=&vlist[3*(t2-1)];
     v3=&vlist[3*(t3-1)];
     //Calculate Normal derivatives (in the original frame)
     Calculate_Normal_Derivative(v1,v2,v3,dndx1,dndx2,dndx3,dndy1,dndy2,dndy3,dndz1,dndz2,dndz3);
     //Calculate echo derivative
     dech=ech*log(mu)*rexp;
     
     
     mult_vector(R,v1,vr1);
     
     mult_vector(R,v2,vr2);
     mult_vector(R,v3,vr3);
     
    // mexPrintf("%f %f %f\n",R[0][0],R[0][1],R[0][2]);
   //  mexPrintf("%f %f %f\n",R[1][0],R[1][1],R[1][2]);
   //  mexPrintf("%f %f %f\n",R[2][0],R[2][1],R[2][2]);
     //Check if sign change is needed:
     for(int i=0;i<3;i++)
     {
      side1[i]=vr2[i]-vr1[i];
      side2[i]=vr3[i]-vr1[i];
     }
     cross(side1,side2,normalr);
     if(normalr[2]<0)
       sign=-1;
     else
       sign=1;
     area=NORM(normalr)*0.5;
     //x coordinate is the delay, y frequency
     //Now we should convert to frequency domain, ie calculate the contribution of each facet
     Calc_FTC(freqx,freqy,nfreq,vr1[0],vr1[1],vr2[0],vr2[1],vr3[0],vr3[1],F0);
    // printf("Fdd: %f %f\n",creal(FTdd[0]),cimag(FTdd[0]));
     //Note that we sum to F at each round, does not work, we need to multiply with the echo
     //Derivatives wrt angles
   
     for(int jf=0;jf<nfreq;jf++)
     {
       tscale=scale*cexp(2*PI*I*(offset[0]*freqx[jf]+offset[1]*freqy[jf]));
    //   mexPrintf("scale:%f offset: %f %f\n",scale,creal(cexp(2*PI*I*(offset[0]*freqx[jf]+offset[1]*freqy[jf]))),cimag(cexp(2*PI*I*(offset[0]*freqx[jf]+offset[1]*freqy[jf]))));
       FTC=tscale*F0[jf];
       Fr[jf]+=creal(sign*ech*FTC);
       Fi[jf]+=cimag(sign*ech*FTC);
      
 }
// printf("j:%d tscale: %f %f scale: %f sign: %d ech: %f\n",j,creal(tscale),cimag(tscale),scale,sign,ech);
 TB+=mu*ech*area;
}
//printf("TB: %7.8f\n",TB);
for(int j=0;j<nfreq;j++)
{
    Fr[j]=Fr[j]/TB;
    Fi[j]=Fi[j]/TB;
}
// write_matrix_file("/tmp/FTr2.txt",Fr,1,nfreq);
//   write_matrix_file("/tmp/FTi2.txt",Fi,1,nfreq);  
free(F0);
free(normal);
free(centroid);
free(IndexofBlocks);
free(NumofBlocks);
}

void Calculate_Frame_Matrix2(double *E,double *angles,double TIME,double freq,double R[3][3])
{
  //E is the vector pointing to the radar (unrotated)
  //angles =beta,lambda,omega
  //R is the output, 3x3 matrix
  double xr[3],yr[3],nyr,zr[3],nzr;
  double w[3],w0[3];
  double M[3][3],dMb[3][3],dMl[3][3],dMo[3][3],A[3][3],Mt[3][3];
  double omegas=angles[2]/(60.0*60.0*24); //rads/sec
  w0[0]=0;
  w0[1]=0;
  w0[2]=omegas;
  
  rotate(angles[0],angles[1],angles[2],angles[3],TIME,M,dMb,dMl,dMo);
  //Transpose Should transpose derivat
  transpose(M,Mt);
  mult_vector(Mt,w0,w); //Here is rotation vector in the global frame
  cross(E,w,yr); //y vector, corresponds to frequency
  nyr=NORM(yr);
  yr[0]=yr[0]/nyr;
  yr[1]=yr[1]/nyr;
  yr[2]=yr[2]/nyr;
  cross(E,yr,zr); //z vector, projection direction
  nzr=NORM(zr);
  zr[0]=zr[0]/nzr;
  zr[1]=zr[1]/nzr;
  zr[2]=zr[2]/nzr;
  Convert_to_Matrix(E,yr,zr,A); //A converts x,y,z coordinates to delay, freq and projdir
  A[0][0]=A[0][0]*(-2.0)/LTS*1.0e6; //In usec
  A[1][1]=A[1][1]*2.0*freq/LTS*nyr;
  //Finally the rotation corresponding to angles
  multmat(A,Mt,R);
}
void Calculate_Frame_Matrix_Derivatives(double *E,double *angles,double TIME,double freq,double R[3][3],double Rdb[3][3],double Rdl[3][3],double Rdo[3][3])
{
  //E is the vector pointing to the radar (unrotated)
  //angles =beta,lambda,omega
  //R is the output, 3x3 matrix
  double xr[3],yr[3],nxr,zr[3],nzr;
  double w[3],w0[3];
  double M[3][3],Mdb[3][3],Mdl[3][3],Mdo[3][3],A[3][3],Mt[3][3];
  double omegas=angles[2]/(60.0*60.0*24); //rads/sec
  //double Erdb[3],Erdl[3],Erdo[3];
  double dwdb[3],dwdl[3],dwdo[3];
  double xru[3],xrudb[3],xrudl[3],xrudo[3];
  double nxrdb,nxrdo,nxrdl;
  double xrdb[3],xrdl[3],xrdo[3];
  double Adb[3][3],Adl[3][3],Ado[3][3];
  double dzdb[3],dzdl[3],dzdo[3];
  double Cdb[3][3]={0,0,0,0,0,0,0,0,0};
  double Cdl[3][3]={0,0,0,0,0,0,0,0,0};
  double Cdo[3][3]={0,0,0,0,0,0,0,0,0};
  double C[3][3]={0,0,0,0,0,0,0,0,0};
  double zerovec[3];
  zerovec[0]=0.0;
  zerovec[1]=0.0;
  zerovec[2]=0.0;
  w0[0]=0;
  w0[1]=0;
  w0[2]=omegas;
  
  yr[0]=-E[0]; //yr points to asteroid
  yr[1]=-E[1];
  yr[2]=-E[2];
  rotate(angles[0],angles[1],angles[2],angles[3],TIME,M,Mdb,Mdl,Mdo);
  //Transpose Should transpose derivat
  transpose(M,Mt);
  mult_vector(Mt,w0,w); //Here is rotation vector in the global frame
  
  cross(E,w,xr); //x vector, corresponds to frequency
  nxr=NORM(xr);
  xr[0]=xr[0]/nxr;
  xr[1]=xr[1]/nxr;
  xr[2]=xr[2]/nxr;
  cross(xr,yr,zr); //z vector, projection direction
  nzr=NORM(zr);
  zr[0]=zr[0]/nzr;
  zr[1]=zr[1]/nzr;
  zr[2]=zr[2]/nzr;
  Convert_to_Matrix(xr,yr,zr,A); //A converts x,y,z coordinates to freq, delay and projdir
  //A[0][0]=A[0][0]*(-2.0)/LTS*1.0e6; //In usec
  //A[1][1]=A[1][1]*2.0*freq/LTS*nyr;
  //Finally the rotation corresponding to angles
  //multmat(A,Mt,R);
  //printf("xr: %f %f %f\n",E[0],E[1],E[2]);
  //printf("yr: %f %f %f\n",yr[0],yr[1],yr[2]);
  //printf("zr: %f %f %f\n",zr[0],zr[1],zr[2]);
  //Derivatives
  //mult_vector(Mdb,E,Erdb);
  //mult_vector(Mdl,E,Erdl);
  //mult_vecotr(Mdo,E,Erdo);
  //Derivative of unnormalized yr
  
  for(int i=0;i<3;i++)
  {
    dwdb[i]=Mdb[2][i]*omegas;
    dwdl[i]=Mdl[2][i]*omegas;
    dwdo[i]=Mdo[2][i]*omegas+M[2][i]*1.0/86400.0;
  }
  cross(E,w,xru);
  cross(E,dwdb,xrudb);
  cross(E,dwdl,xrudl);
  cross(E,dwdo,xrudo);
  //Derivative of the norm
  nxrdb=(DOT(xru,xrudb))*1.0/nxr;
  nxrdl=(DOT(xru,xrudl))*1.0/nxr;
  nxrdo=(DOT(xru,xrudo))*1.0/nxr;
 // printf("nyr: \n %f %f %f\n",nyrdb,nyrdl,nyrdo);
  //Derivative of normalized vector
  for(int i=0;i<3;i++)
  {
    xrdb[i]=(xrudb[i]*nxr-xru[i]*nxrdb)/pow(nxr,2);
    xrdl[i]=(xrudl[i]*nxr-xru[i]*nxrdl)/pow(nxr,2);
    xrdo[i]=(xrudo[i]*nxr-xru[i]*nxrdo)/pow(nxr,2);
  }
  cross(xrdb,yr,dzdb);
  cross(xrdl,yr,dzdl);
  cross(xrdo,yr,dzdo);
  Convert_to_Matrix(xrdb,zerovec,dzdb,Adb);
  Convert_to_Matrix(xrdl,zerovec,dzdl,Adl);
  Convert_to_Matrix(xrdo,zerovec,dzdo,Ado);
  C[0][0]=2.0*freq/LTS*nxr;
  C[1][1]=2.0/LTS*1e6;
  C[2][2]=1;
  Cdb[0][0]=2*freq/LTS*nxrdb;
  Cdl[0][0]=2*freq/LTS*nxrdl;
  Cdo[0][0]=2*freq/LTS*nxrdo;
  //Conversion matrices
  double Rdb1[3][3],Rdb2[3][3],Rdb3[3][3];
  double Rdl1[3][3],Rdl2[3][3],Rdl3[3][3];
  double Rdo1[3][3],Rdo2[3][3],Rdo3[3][3];
  multmat3t(Cdb,A,M,Rdb1);
  multmat3t(C,Adb,M,Rdb2);
  multmat3t(C,A,Mdb,Rdb3);
  
  multmat3t(Cdl,A,M,Rdl1);
  multmat3t(C,Adl,M,Rdl2);
  multmat3t(C,A,Mdl,Rdl3);
  multmat3t(Cdo,A,M,Rdo1);
  multmat3t(C,Ado,M,Rdo2);
  multmat3t(C,A,Mdo,Rdo3);
  multmat3t(C,A,M,R);
  //printf("nyrdo:%f\n",nyrdo);
  //printf("w:%f %f %f\n",w[0],w[1],w[2]);
   //printf("C:\n");
    
  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
  {
    Rdb[i][j]=Rdb1[i][j]+Rdb2[i][j]+Rdb3[i][j];
    Rdl[i][j]=Rdl1[i][j]+Rdl2[i][j]+Rdl3[i][j];
    Rdo[i][j]=Rdo1[i][j]+Rdo2[i][j]+Rdo3[i][j];
  }
 // for(int db=0;db<3;db++)
  // printf("%f %f %f\n",Cdo[db][0],Cdo[db][1],Cdo[db][2]);
   
}


  
  
