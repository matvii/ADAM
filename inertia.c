#include"utils.h"
#include"matrix_ops.h"
double Pderiv(int h,int k,double vol,double *D,double *E,double *F,double *dP);
void inertia(int *tlist,double* vlist,int nfac,int nvert,double *D1,int m,int n,double *result,double *angle,double *dres)
{
    /*
     * Calculate angle between the inertia axis corresponding to the maximum moment of inertia, and the z-axis.
     * OUTPUT: result=(1-cos(angle)^2)
     * dres=1x3*nvert derivative wrt vertices [or 1x3*n is D1\=NULL]
     */
    double *dres2=calloc(3*nvert,sizeof(double));
   double CM[3]={0.0,0.0,0.0};
   double Pxx=0,Pyy=0.0,Pzz=0.0,Pyz=0.0,Pxz=0.0,Pxy=0.0;
   int v1,v2,v3;
   int h,k;
   double *D,*E,*F;
   double G[3],H[3],C[3],S[3];
   double vol=0.0;
   double tvol=0.0;
   double *dCM;
   dCM=calloc(3*3*nvert,sizeof(double));
   double *dPxx,*dPyy,*dPzz,*dPxy,*dPxz,*dPyz,*dP;
   dPxx=calloc(3*nvert,sizeof(double));
   dPyy=calloc(3*nvert,sizeof(double));
   dPzz=calloc(3*nvert,sizeof(double));
   dPxy=calloc(3*nvert,sizeof(double));
   dPxz=calloc(3*nvert,sizeof(double));
   dPyz=calloc(3*nvert,sizeof(double));
   dP=calloc(9,sizeof(double));
   double *dV=calloc(3*nvert,sizeof(double));
   double *dcm1,*dcm2,*dcm3;
   dcm1=calloc(3*nvert,sizeof(double));
   dcm2=calloc(3*nvert,sizeof(double));
   dcm3=calloc(3*nvert,sizeof(double));
   double dIc0,dIc1,dIc2,dIc4,dIc5,dIc8;
   double x0,y0,z0,x1,y1,z1,x2,y2,z2;
   double cm1,cm2,cm3,cm12,cm22,cm32;
   double *Im,*Ic,*M,*L;
   double P,U,T,eig;
   double vz,res;
   double ctheta=0.0,theta;
   double vec1[3],vec2[3],vec3[3];
   double max;
   double n1,n2,n3;
   double dT,dU;
   double dtheta,dctheta;
   double tvol2;
   double deig;
   double dL0,dL1,dL2,dL3,dL4,dL5,dL6,dL7,dL8;
   double dvzdL0,dvzdL1,dvzdL2,dvzdL4,dvzdL5,dvzdL8;
   int mindex;
   L=calloc(9,sizeof(double));
   Im=calloc(9,sizeof(double));
   Ic=calloc(9,sizeof(double));
   M=calloc(9,sizeof(double));
   double L0,L1,L2,L3,L4,L5,L6,L7,L8;
   for(int j=0;j<nfac;j++)
   {
       v1=tlist[3*j]-1;
       v2=tlist[3*j+1]-1;
       v3=tlist[3*j+2]-1;
       D=vlist+3*v1;
       E=vlist+3*v2;
       F=vlist+3*v3;
       x0=D[0];
       y0=D[1];
       z0=D[2];
       x1=E[0];
       y1=E[1];
       z1=E[2];
       x2=F[0];
       y2=F[1];
       z2=F[2];
       G[0]=E[0]-D[0];
       G[1]=E[1]-D[1];
       G[2]=E[2]-D[2];
       
       H[0]=F[0]-D[0];
       H[1]=F[1]-D[1];
       H[2]=F[2]-D[2];
       cross(G,H,C);
       vol=1.0/6.0*(D[0]*C[0]+D[1]*C[1]+D[2]*C[2]);
       
       S[0]=(D[0]+E[0]+F[0])/4;
       S[1]=(D[1]+E[1]+F[1])/4;
       S[2]=(D[2]+E[2]+F[2])/4;
       
       CM[0]=CM[0]+vol*S[0];
       CM[1]=CM[1]+vol*S[1];
       CM[2]=CM[2]+vol*S[2];
       //Derivatives of CM wrt vertex coordinates
       //x-coordinates, CM[0]
       dCM[v1]+=(x1*y2*(z0 - z1) + x2*(y0 - y2)*z1 - 2*x0*y2*z1 + 2*x0*y1*z2 + x1*(-y0 + y1)*z2 + x2*y1*(-z0 + z2))/24.0;
       dCM[v2]+=(2*x1*(y2*z0 - y0*z2) + x2*(-(y1*z0) + y2*z0 + y0*z1 - y0*z2) + x0*(y2*z0 - y2*z1 - y0*z2 + y1*z2))/24.0;
       dCM[v3]+=(-2*x2*y1*z0 + 2*x2*y0*z1 + x0*(y0 - y2)*z1 + x1*(-(y1*z0) + y2*z0 + y0*(z1 - z2)) + x0*y1*(-z0 + z2))/24.0;
       //y-coordinates, CM[0]
       dCM[nvert+v1]+=-((x0 + x1 + x2)*(-(x2*z1) + x1*z2))/24.0;
       dCM[nvert+v2]+=((x0 + x1 + x2)*(-(x2*z0) + x0*z2))/24.0;
       dCM[nvert+v3]+=-((x0 + x1 + x2)*(-(x1*z0) + x0*z1))/24.0;
       //z-coordinates, CM[0]
       dCM[2*nvert+v1]+=((x0 + x1 + x2)*(-(x2*y1) + x1*y2))/24.0;
       dCM[2*nvert+v2]+=-((x0 + x1 + x2)*(-(x2*y0) + x0*y2))/24.0;
       dCM[2*nvert+v3]+=((x0 + x1 + x2)*(-(x1*y0) + x0*y1))/24.0;
       //CM[1]
       dCM[v1+3*nvert]+=((y0 + y1 + y2)*(-(y2*z1) + y1*z2))/24.0;
       dCM[v2+3*nvert]+=-((y0 + y1 + y2)*(-(y2*z0) + y0*z2))/24.0;
       dCM[v3+3*nvert]+=((y0 + y1 + y2)*(-(y1*z0) + y0*z1))/24.0;
       
       dCM[nvert+v1+3*nvert]+=(-(x0*y2*z1) + x2*((2*y0 + y2)*z1 + y1*(-z0 + z1)) + x0*y1*z2 + x1*(y2*(z0 - z2) - (2*y0 + y1)*z2))/24.0;
       dCM[nvert+v2+3*nvert]+=(-(x2*((2*y1 + y2)*z0 + y0*(z0 - z1))) + x1*(y2*z0 - y0*z2) + x0*(-(y2*z1) + y0*z2 + 2*y1*z2 + y2*z2))/24.0;
       dCM[nvert+v3+3*nvert]+=(x2*(-(y1*z0) + y0*z1) + x1*((y1 + 2*y2)*z0 + y0*(z0 - z2)) - x0*(y0*z1 + y1*z1 + 2*y2*z1 - y1*z2))/24.0;
       
       dCM[2*nvert+v1+3*nvert]+=-((y0 + y1 + y2)*(x2*y1 - x1*y2))/24.0;
       dCM[2*nvert+v2+3*nvert]+=((y0 + y1 + y2)*(x2*y0 - x0*y2))/24.0;
       dCM[2*nvert+v3+3*nvert]+=((-(x1*y0) + x0*y1)*(y0 + y1 + y2))/24.0;
       
       //CM[2]
       dCM[v1+2*3*nvert]+=-((z0 + z1 + z2)*(y2*z1 - y1*z2))/24.0;
       dCM[v2+2*3*nvert]+=((z0 + z1 + z2)*(y2*z0 - y0*z2))/24.0;
       dCM[v3+2*3*nvert]+=((-(y1*z0) + y0*z1)*(z0 + z1 + z2))/24.0;
       
       dCM[v1+nvert+2*3*nvert]+=((z0 + z1 + z2)*(x2*z1 - x1*z2))/24.0;
       dCM[v2+nvert+2*3*nvert]+=-((z0 + z1 + z2)*(x2*z0 - x0*z2))/24.0;
       dCM[v3+nvert+2*3*nvert]+=((x1*z0 - x0*z1)*(z0 + z1 + z2))/24.0;
       
       dCM[v1+2*nvert+2*3*nvert]+=(-(x0*y2*z1) + x0*y1*z2 - x2*(-(y0*z1) + y1*(2*z0 + z1 + z2)) + x1*(-(y0*z2) + y2*(2*z0 + z1 + z2)))/24.0;
       dCM[v2+2*nvert+2*3*nvert]+=(x1*(y2*z0 - y0*z2) + x2*(-(y1*z0) + y0*(z0 + 2*z1 + z2)) - x0*(-(y1*z2) + y2*(z0 + 2*z1 + z2)))/24.0;
       dCM[v3+2*nvert+2*3*nvert]+=(-(x2*y1*z0) + x2*y0*z1 - x1*(-(y2*z0) + y0*(z0 + z1 + 2*z2)) + x0*(-(y2*z1) + y1*(z0 + z1 + 2*z2)))/24.0;
       
       Pxx=Pxx+Pderiv(0,0,vol,D,E,F,dP);
       dPxx[v1]+=dP[0];
       dPxx[v2]+=dP[1];
       dPxx[v3]+=dP[2];
       dPxx[v1+nvert]+=dP[3];
       dPxx[v2+nvert]+=dP[4];
       dPxx[v3+nvert]+=dP[5];
       dPxx[v1+2*nvert]+=dP[6];
       dPxx[v2+2*nvert]+=dP[7];
       dPxx[v3+2*nvert]+=dP[8];
       Pyy=Pyy+Pderiv(1,1,vol,D,E,F,dP);
       dPyy[v1]+=dP[0];
       dPyy[v2]+=dP[1];
       dPyy[v3]+=dP[2];
       dPyy[v1+nvert]+=dP[3];
       dPyy[v2+nvert]+=dP[4];
       dPyy[v3+nvert]+=dP[5];
       dPyy[v1+2*nvert]+=dP[6];
       dPyy[v2+2*nvert]+=dP[7];
       dPyy[v3+2*nvert]+=dP[8];
       Pzz=Pzz+Pderiv(2,2,vol,D,E,F,dP);
       dPzz[v1]+=dP[0];
       dPzz[v2]+=dP[1];
       dPzz[v3]+=dP[2];
       dPzz[v1+nvert]+=dP[3];
       dPzz[v2+nvert]+=dP[4];
       dPzz[v3+nvert]+=dP[5];
       dPzz[v1+2*nvert]+=dP[6];
       dPzz[v2+2*nvert]+=dP[7];
       dPzz[v3+2*nvert]+=dP[8];
       Pxz=Pxz+Pderiv(0,2,vol,D,E,F,dP);
       dPxz[v1]+=dP[0];
       dPxz[v2]+=dP[1];
       dPxz[v3]+=dP[2];
       dPxz[v1+nvert]+=dP[3];
       dPxz[v2+nvert]+=dP[4];
       dPxz[v3+nvert]+=dP[5];
       dPxz[v1+2*nvert]+=dP[6];
       dPxz[v2+2*nvert]+=dP[7];
       dPxz[v3+2*nvert]+=dP[8];
       Pxy=Pxy+Pderiv(0,1,vol,D,E,F,dP);
       dPxy[v1]+=dP[0];
       dPxy[v2]+=dP[1];
       dPxy[v3]+=dP[2];
       dPxy[v1+nvert]+=dP[3];
       dPxy[v2+nvert]+=dP[4];
       dPxy[v3+nvert]+=dP[5];
       dPxy[v1+2*nvert]+=dP[6];
       dPxy[v2+2*nvert]+=dP[7];
       dPxy[v3+2*nvert]+=dP[8];
       Pyz=Pyz+Pderiv(1,2,vol,D,E,F,dP);
       dPyz[v1]+=dP[0];
       dPyz[v2]+=dP[1];
       dPyz[v3]+=dP[2];
       dPyz[v1+nvert]+=dP[3];
       dPyz[v2+nvert]+=dP[4];
       dPyz[v3+nvert]+=dP[5];
       dPyz[v1+2*nvert]+=dP[6];
       dPyz[v2+2*nvert]+=dP[7];
       dPyz[v3+2*nvert]+=dP[8];
    
       tvol=tvol+vol;
       dV[v1]+=(-(y2*z1) + y1*z2)/6.0;
       dV[v2]+=(y2*z0 - y0*z2)/6.0;
       dV[v3]+=(-(y1*z0) + y0*z1)/6.0;
       
       dV[v1+nvert]+=(x2*z1 - x1*z2)/6.0;
       dV[v2+nvert]+=(-(x2*z0) + x0*z2)/6.0;
       dV[v3+nvert]+=(x1*z0 - x0*z1)/6.0;
       
       dV[v1+2*nvert]+=(-(x2*y1) + x1*y2)/6.0;
       dV[v2+2*nvert]+=(x2*y0 - x0*y2)/6.0;
       dV[v3+2*nvert]+=(-(x1*y0) + x0*y1)/6.0;
       
   }
      
       CM[0]=CM[0]/tvol;
       CM[1]=CM[1]/tvol;
       CM[2]=CM[2]/tvol;
       
       cm1=CM[0];
       cm2=CM[1];
       cm3=CM[2];
       cm12=pow(cm1,2);
       cm22=pow(cm2,2);
       cm32=pow(cm3,2);
       tvol2=pow(tvol,2);
       for(int j=0;j<nvert;j++)
       {
           dcm1[j]=(dCM[j]*tvol-CM[0]*dV[j])/tvol2;
           dcm1[j+nvert]=(dCM[j+nvert]*tvol-CM[0]*dV[j+nvert])/tvol2;
           dcm1[j+2*nvert]=(dCM[j+2*nvert]*tvol-CM[0]*dV[j+2*nvert])/tvol2;
           
           dcm2[j]=(dCM[j+3*nvert]*tvol-CM[1]*dV[j])/tvol2;
           dcm2[j+nvert]=(dCM[j+nvert+3*nvert]*tvol-CM[1]*dV[j+nvert])/tvol2;
           dcm2[j+2*nvert]=(dCM[j+2*nvert+3*nvert]*tvol-CM[1]*dV[j+2*nvert])/tvol2;
           
           dcm3[j]=(dCM[j+2*3*nvert]*tvol-CM[2]*dV[j])/tvol2;
           dcm3[j+nvert]=(dCM[j+nvert+2*3*nvert]*tvol-CM[2]*dV[j+nvert])/tvol2;
           dcm3[j+2*nvert]=(dCM[j+2*nvert+2*3*nvert]*tvol-CM[2]*dV[j+2*nvert])/tvol2;
       }
       free(dCM);
       free(dP);
       double dP0;
       Im[0]=Pyy+Pzz;
       Im[1]=-Pxy;
       Im[2]=-Pxz;
       Im[3]=-Pxy;
       Im[4]=Pxx+Pzz;
       Im[5]=-Pyz;
       Im[6]=-Pxz;
       Im[7]=-Pyz;
       Im[8]=Pxx+Pyy;
       
       Ic[0]=Im[0]-tvol*(cm22+cm32);
       Ic[1]=Im[1]-tvol*(-cm1*cm2);
       Ic[2]=Im[2]-tvol*(-cm1*cm3);
       Ic[3]=Im[3]-tvol*(-cm1*cm2);
       Ic[4]=Im[4]-tvol*(cm12+cm32);
       Ic[5]=Im[5]-tvol*(-cm2*cm3);
       Ic[6]=Im[6]-tvol*(-cm1*cm3);
       Ic[7]=Im[7]-tvol*(-cm2*cm3);
       Ic[8]=Im[8]-tvol*(cm12+cm22);
       double Ic0,Ic1,Ic2,Ic4,Ic5,Ic8;
       Ic0=Ic[0];
       Ic1=Ic[1];
       Ic2=Ic[2];
       Ic4=Ic[4];
       Ic5=Ic[5];
       Ic8=Ic[8];
       
       
      T=Ic[0]+Ic[4]+Ic[8];
       P=Ic[0]*Ic[4]+Ic[0]*Ic[8]+Ic[4]*Ic[8]-Ic[1]*Ic[1]-Ic[2]*Ic[2]-Ic[5]*Ic[5];
       U=sqrt(pow(T,2)-3*P)/3.0;
       
       ctheta=(-2.0*pow(T,3)+9*T*P-27*DET(Ic))/(54.0*pow(U,3));
       if(ctheta<-1)
           ctheta=-1;
       if(ctheta>1)
           ctheta=1;
       theta=acos(ctheta);
       
       double eig1,eig2,eig3;
       int eindex=1,vindex=-1;
       eig1=T/3.0-2*U*cos(theta/3.0);
       eig2=T/3.0-2*U*cos(theta/3.0-2.0*PI/3.0);
       eig3=T/3.0-2*U*cos(theta/3.0+2.0*PI/3.0);
       eig=eig1;
       if(eig<eig2)
       {
           eig=eig2;
           eindex=2;
       }
       if(eig<eig3)
       {
           eig=eig3;
           eindex=3;
       }
       
      
       L[0]=Ic[0]-eig;
       L[1]=Ic[1];
       L[2]=Ic[2];
       L[3]=Ic[3];
       L[4]=Ic[4]-eig;
       L[5]=Ic[5];
       L[6]=Ic[6];
       L[7]=Ic[7];
       L[8]=Ic[8]-eig;
       
       L0=L[0];
       L1=L[1];
       L2=L[2];
       L3=L[3];
       L4=L[4];
       L5=L[5];
       L6=L[6];
       L7=L[7];
       L8=L[8];
       cross(L,L+3,vec1);
       cross(L,L+6,vec2);
       cross(L+3,L+6,vec3);
       
       n1=NORM(vec1);
       n2=NORM(vec2);
       n3=NORM(vec3);
       
        if(n1>EP)
       {
           vz=vec1[2]/n1;
           *result=1-vz*vz;
           *angle=acos(vz)*180/PI;
           vindex=1;
       }
       else if(n2>EP)
       {
           vz=vec2[2]/n2;
           *result=1-vz*vz;
           *angle=acos(vz)*180/PI;
           vindex=2;
       }
       else if(n3>EP)
       {
           vz=vec3[2]/n3;
           *result=1-vz*vz;
           *angle=acos(vz)*180/PI;
           vindex=3;
       }
       else
       {
           fprintf(stderr,"Something went wrong in inertia calculation, all n1,n2 and n3 <EP\n");
           exit(-1);
       }
       for(int j=0;j<3*nvert;j++)
       {
           dIc0=dPyy[j]+dPzz[j]-dV[j]*(cm22+cm32)-tvol*(2*dcm2[j]*cm2+2*dcm3[j]*cm3);
           dIc1=-dPxy[j]-dV[j]*(-cm1*cm2)-tvol*(-dcm1[j]*cm2-cm1*dcm2[j]);
           dIc2=-dPxz[j]-dV[j]*(-cm1*cm3)-tvol*(-dcm1[j]*cm3-cm1*dcm3[j]);
           dIc4=dPxx[j]+dPzz[j]-dV[j]*(cm12+cm32)-tvol*(2*dcm1[j]*cm1+2*dcm3[j]*cm3);
           dIc5=-dPyz[j]-dV[j]*(-cm2*cm3)-tvol*(-dcm2[j]*cm3-cm2*dcm3[j]);
           dIc8=dPxx[j]+dPyy[j]-dV[j]*(cm12+cm22)-tvol*(2*dcm1[j]*cm1+2*dcm2[j]*cm2);
           dT=dIc0+dIc4+dIc8;
           //dP0[j]=dIc0[j]*Ic[4]+Ic[0]*dIc4[j]+dIc0[j]*Ic[8]+Ic[0]*dIc8[j]+dIc4[j]*Ic[8]+Ic[4]*dIc8[j]-2*Ic[1]*dIc1[j]-2*Ic[2]*dIc2[j]-2*Ic[5]*dIc5[j];
           dP0=dIc0*Ic[4]+Ic[0]*dIc4+dIc0*Ic[8]+Ic[0]*dIc8+dIc4*Ic[8]+Ic[4]*dIc8-2*Ic[1]*dIc1-2*Ic[2]*dIc2-2*Ic[5]*dIc5;
           dU=T/(3.0*sqrt(-3*P+pow(T,2)))*dT-1/(2*sqrt(-3*P+pow(T,2)))*dP0;
           dctheta=((-27*(2*pow(Ic1,4) + 2*pow(Ic2,4) + 6*Ic1*Ic2*Ic5*(-2*Ic0 + Ic4 + Ic8) + pow(Ic1,2)*(4*pow(Ic2,2) - 3*Ic0*Ic4 + pow(Ic4,2) - 2*pow(Ic5,2) + 3*Ic0*Ic8 + Ic4*Ic8 - 2*pow(Ic8,2)) + pow(Ic2,2)*(-2*pow(Ic4,2) - 2*pow(Ic5,2) + 3*Ic0*(Ic4 - Ic8) + Ic4*Ic8 + pow(Ic8,2)) + (pow(Ic4,2) + 4*pow(Ic5,2) - 2*Ic4*Ic8 + pow(Ic8,2))*(pow(Ic0,2) - pow(Ic5,2) + Ic4*Ic8 - Ic0*(Ic4 + Ic8))))/(4.*pow(pow(Ic0,2) + 3*pow(Ic1,2) + 3*pow(Ic2,2) + pow(Ic4,2) +3*pow(Ic5,2) - Ic4*Ic8 + pow(Ic8,2) - Ic0*(Ic4 + Ic8),2.5)))*dIc0+((27*(12*pow(Ic1,2)*Ic2*Ic5 + pow(Ic1,3)*(Ic4 - 2*Ic8) - 2*Ic2*Ic5*(3*pow(Ic2,2) + pow(Ic4,2) + 3*pow(Ic5,2) - Ic4*Ic8 + pow(Ic8,2)) + pow(Ic0,2)*(-2*Ic2*Ic5 + Ic1*(-Ic4 + Ic8)) + Ic1*(Ic4*pow(Ic5,2) + pow(Ic4,2)*Ic8 + 7*pow(Ic5,2)*Ic8 - 3*Ic4*pow(Ic8,2) + 2*pow(Ic8,3) + pow(Ic2,2)*(-8*Ic4 + 7*Ic8)) + Ic0*(pow(Ic1,3) + 2*Ic2*Ic5*(Ic4 + Ic8) + Ic1*(pow(Ic2,2) - pow(Ic4,2) - 8*pow(Ic5,2) + 4*Ic4*Ic8 - 3*pow(Ic8,2)))))/(2.*pow(pow(Ic0,2) + 3*pow(Ic1,2) + 3*pow(Ic2,2) + pow(Ic4,2) + 3*pow(Ic5,2) - Ic4*Ic8 + pow(Ic8,2) - Ic0*(Ic4 + Ic8),2.5)))*dIc1+((-18*(Ic0*Ic2 + 3*Ic1*Ic5 + Ic2*(-2*Ic4 + Ic8))*(pow(Ic0,2) + 3*pow(Ic1,2) + 3*pow(Ic2,2) + pow(Ic4,2) + 3*pow(Ic5,2) - Ic4*Ic8 + pow(Ic8,2) - Ic0*(Ic4 + Ic8)) -9*Ic2*(-2*pow(Ic0 + Ic4 + Ic8,3) + 9*(Ic0 + Ic4 + Ic8)*(-pow(Ic1,2) - pow(Ic2,2) + Ic0*Ic4 - pow(Ic5,2) + Ic0*Ic8 + Ic4*Ic8) + 27*(pow(Ic2,2)*Ic4 - 2*Ic1*Ic2*Ic5 + pow(Ic1,2)*Ic8 + Ic0*(pow(Ic5,2) - Ic4*Ic8))))/(2.*pow(pow(Ic0,2) + 3*pow(Ic1,2) + 3*pow(Ic2,2) + pow(Ic4,2) + 3*pow(Ic5,2) - Ic4*Ic8 + pow(Ic8,2) - Ic0*(Ic4 + Ic8),2.5)))*dIc2+((6*(pow(Ic0,2) - 3*pow(Ic1,2) + 6*pow(Ic2,2) - 2*pow(Ic4,2) - 3*pow(Ic5,2) + 2*Ic0*(Ic4 - 2*Ic8) + 2*Ic4*Ic8 + pow(Ic8,2))*(pow(Ic0,2) + 3*pow(Ic1,2) + 3*pow(Ic2,2) + pow(Ic4,2) + 3*pow(Ic5,2) - Ic4*Ic8 + pow(Ic8,2) - Ic0*(Ic4 + Ic8)) - 3*(-Ic0 + 2*Ic4 - Ic8)*(-2*pow(Ic0 + Ic4 + Ic8,3) + 9*(Ic0 + Ic4 + Ic8)*(-pow(Ic1,2) - pow(Ic2,2) + Ic0*Ic4 - pow(Ic5,2) + Ic0*Ic8 + Ic4*Ic8) + 27*(pow(Ic2,2)*Ic4 - 2*Ic1*Ic2*Ic5 + pow(Ic1,2)*Ic8 + Ic0*(pow(Ic5,2) - Ic4*Ic8))))/(4.*pow(pow(Ic0,2) + 3*pow(Ic1,2) + 3*pow(Ic2,2) + pow(Ic4,2) + 3*pow(Ic5,2) - Ic4*Ic8 + pow(Ic8,2) - Ic0*(Ic4 + Ic8),2.5)))*dIc4+((27*(-6*pow(Ic1,3)*Ic2 + 2*pow(Ic0,3)*Ic5 + pow(Ic1,2)*Ic5*(Ic4 - 8*Ic8) - 2*Ic1*Ic2*(3*pow(Ic2,2) + pow(Ic4,2) - 6*pow(Ic5,2) - Ic4*Ic8 + pow(Ic8,2)) + Ic5*(Ic4*pow(Ic5,2) - pow(Ic4,2)*Ic8 + pow(Ic5,2)*Ic8 - Ic4*pow(Ic8,2) + pow(Ic2,2)*(-8*Ic4 + Ic8)) - pow(Ic0,2)*(2*Ic1*Ic2 + 3*Ic5*(Ic4 + Ic8)) + Ic0*(7*pow(Ic1,2)*Ic5 + 2*Ic1*Ic2*(Ic4 + Ic8) + Ic5*(7*pow(Ic2,2) + pow(Ic4,2) - 2*pow(Ic5,2) + 4*Ic4*Ic8 + pow(Ic8,2)))))/(2.*pow(pow(Ic0,2) + 3*pow(Ic1,2) + 3*pow(Ic2,2) + pow(Ic4,2) + 3*pow(Ic5,2) - Ic4*Ic8 + pow(Ic8,2) - Ic0*(Ic4 + Ic8),2.5)))*dIc5+((6*(pow(Ic0,2) + 6*pow(Ic1,2) - 3*pow(Ic2,2) - 4*Ic0*Ic4 + pow(Ic4,2) - 3*pow(Ic5,2) + 2*Ic0*Ic8 + 2*Ic4*Ic8 - 2*pow(Ic8,2))*(pow(Ic0,2) + 3*pow(Ic1,2) + 3*pow(Ic2,2) + pow(Ic4,2) + 3*pow(Ic5,2) - Ic4*Ic8 + pow(Ic8,2) - Ic0*(Ic4 + Ic8)) - 3*(-Ic0 - Ic4 + 2*Ic8)*(-2*pow(Ic0 + Ic4 + Ic8,3) + 9*(Ic0 + Ic4 + Ic8)*(-pow(Ic1,2) - pow(Ic2,2) + Ic0*Ic4 - pow(Ic5,2) + Ic0*Ic8 + Ic4*Ic8) + 27*(pow(Ic2,2)*Ic4 - 2*Ic1*Ic2*Ic5 + pow(Ic1,2)*Ic8 + Ic0*(pow(Ic5,2) - Ic4*Ic8))))/(4.*pow(pow(Ic0,2) + 3*pow(Ic1,2) + 3*pow(Ic2,2) + pow(Ic4,2) + 3*pow(Ic5,2) - Ic4*Ic8 + pow(Ic8,2) - Ic0*(Ic4 + Ic8),2.5)))*dIc8;
           
           dtheta=-1.0/sqrt(1-pow(ctheta,2))*dctheta;
            if(isnormal(dtheta)!=1)
               dtheta=0.0;
           if(eindex==1)
               // T/3.0-2*U*cos(theta/3.0);
               deig=dT/3-2*dU*cos(theta/3.0)+2*U*sin(theta/3.0)*dtheta/3.0;
           else if(eindex==2)
               deig=dT/3-2*dU*cos(theta/3.0-2.0*PI/3.0)+2*U*sin(theta/3.0-2.0*PI/3.0)*dtheta/3.0;
           else if(eindex==3)
               deig=dT/3-2*dU*cos(theta/3.0+2.0*PI/3.0)+2*U*sin(theta/3.0+2.0*PI/3.0)*dtheta/3.0;
           dL0=dIc0-deig;
           dL1=dIc1;
           dL2=dIc2;
           dL3=dIc1;
           dL4=dIc4-deig;
           dL5=dIc5;
           dL6=dIc2;
           dL7=dIc5;
           dL8=dIc8-deig;
           
           if(vindex==1)
           {
               dvzdL0=((L2*L4 - L1*L5)*(pow(L1,2)*L2 + L2*pow(L4,2) - L1*(L0 + L4)*L5))/pow(pow(pow(L1,2) - L0*L4,2) + pow(L1*L2 - L0*L5,2) + pow(L2*L4 - L1*L5,2),1.5);
               
               dvzdL1=(3*pow(L1,2)*L2*(L0 + L4)*L5 + L0*L2*L4*(L0 + L4)*L5 - pow(L1,3)*(pow(L2,2) + pow(L5,2)) -L1*(2*pow(L2,2)*pow(L4,2) + 2*pow(L0,2)*pow(L5,2) + L0*L4*(pow(L2,2) + pow(L5,2))))/pow(pow(pow(L1,2) - L0*L4,2) +pow(L1*L2 - L0*L5,2) + pow(L2*L4 - L1*L5,2),1.5);
               
               dvzdL2=((pow(L1,2) - L0*L4)*(pow(L1,2)*L2 + L2*pow(L4,2) - L1*(L0 + L4)*L5))/pow(pow(pow(L1,2) - L0*L4,2) + pow(L1*L2 - L0*L5,2) + pow(L2*L4 - L1*L5,2),1.5);
               
               dvzdL4=((-(L1*L2) + L0*L5)*(-(L0*L1*L2) + pow(L0,2)*L5 + L1*(-(L2*L4) + L1*L5)))/pow(pow(pow(L1,2) - L0*L4,2) + pow(L1*L2 - L0*L5,2) + pow(L2*L4 - L1*L5,2),1.5);
               
               dvzdL5=((pow(L1,2) - L0*L4)*(-(L0*L1*L2) + pow(L0,2)*L5 + L1*(-(L2*L4) + L1*L5)))/pow(pow(pow(L1,2) - L0*L4,2) + pow(L1*L2 - L0*L5,2) + pow(L2*L4 - L1*L5,2),1.5);
               
               dvzdL8=0;
           }
           else if(vindex==2)
           {
               dvzdL0=((L2*L5 - L1*L8)*(pow(L2,3) - L1*L5*L8 + L2*(pow(L5,2) - L0*L8)))/pow(pow(L1*L2 - L0*L5,2) + pow(pow(L2,2) - L0*L8,2) + pow(L2*L5 - L1*L8,2),1.5);
               
               dvzdL1=-(((pow(L2,2) - L0*L8)*(pow(L2,3) - L1*L5*L8 + L2*(pow(L5,2) - L0*L8)))/pow(pow(L1*L2 - L0*L5,2) + pow(pow(L2,2) - L0*L8,2) + pow(L2*L5 - L1*L8,2),1.5));
               
               dvzdL2=(pow(L1,2)*L2*L5*L8 - pow(L1,3)*pow(L8,2) + L0*L2*L5*(-2*pow(L2,2) - pow(L5,2) + 2*L0*L8) + L1*(pow(L2,4) + L0*L8*(pow(L5,2) - L0*L8)))/pow(pow(L1*L2 - L0*L5,2) + pow(pow(L2,2) - L0*L8,2) + pow(L2*L5 - L1*L8,2),1.5);
               
               dvzdL4=0;
               
               dvzdL5=((pow(L2,2) - L0*L8)*(L0*pow(L2,2) - pow(L0,2)*L8 + L1*(L2*L5 - L1*L8)))/pow(pow(L1*L2 - L0*L5,2) + pow(pow(L2,2) - L0*L8,2) + pow(L2*L5 - L1*L8,2),1.5);
               
               dvzdL8=((L1*L2 - L0*L5)*(-(L0*pow(L2,2)) + pow(L0,2)*L8 + L1*(-(L2*L5) + L1*L8)))/pow(pow(L1*L2 - L0*L5,2) + pow(pow(L2,2) - L0*L8,2) + pow(L2*L5 - L1*L8,2),1.5);
           }
           else if(vindex=3)
           {
               dvzdL0=0;
               
               dvzdL1=((pow(L5,2) - L4*L8)*(pow(L2,2)*L5 + pow(L5,3) - L1*L2*L8 - L4*L5*L8))/pow(pow(L2*L4 - L1*L5,2) + pow(L2*L5 - L1*L8,2) + pow(pow(L5,2) - L4*L8,2),1.5);
               
               dvzdL2=((pow(L5,2) - L4*L8)*(-(L1*L2*L5) + pow(L1,2)*L8 + L4*(-pow(L5,2) + L4*L8)))/pow(pow(L2*L4 - L1*L5,2) + pow(L2*L5 - L1*L8,2) + pow(pow(L5,2) - L4*L8,2),1.5);
               
               dvzdL4=-(((L2*L5 - L1*L8)*(pow(L2,2)*L5 + pow(L5,3) - L1*L2*L8 - L4*L5*L8))/pow(pow(L2*L4 - L1*L5,2) + pow(L2*L5 - L1*L8,2) + pow(pow(L5,2) - L4*L8,2),1.5));
               
               dvzdL5=(pow(L2,3)*L4*L5 - L1*pow(L5,4) - L1*pow(L2,2)*L4*L8 + L1*(pow(L1,2) + pow(L4,2))*pow(L8,2) + L2*L5*(2*L4*pow(L5,2) - pow(L1,2)*L8 - 2*pow(L4,2)*L8))/pow(pow(L2*L4 - L1*L5,2) + pow(L2*L5 - L1*L8,2) + pow(pow(L5,2) - L4*L8,2),1.5);
               
               dvzdL8=-(((L2*L4 - L1*L5)*(L1*L2*L5 - pow(L1,2)*L8 + L4*(pow(L5,2) - L4*L8)))/pow(pow(L2*L4 - L1*L5,2) + pow(L2*L5 - L1*L8,2) + pow(pow(L5,2) - L4*L8,2),1.5));
           }
           
               
          dres2[j]=-2*vz*(dvzdL0*dL0+dvzdL1*dL1+dvzdL2*dL2+dvzdL4*dL4+dvzdL5*dL5+dvzdL8*dL8);    
    }

if(D1==NULL)
{
    memcpy(dres,dres2,3*nvert*sizeof(double));
    free(dres2);
}
else
{
   matrix_prod(dres2,1,nvert,D1,n,dres);
   matrix_prod(dres2+nvert,1,nvert,D1,n,dres+n);
   matrix_prod(dres2+2*nvert,1,nvert,D1,n,dres+2*n);
   free(dres2);
}
free(dPxx);
free(dPyy);
free(dPzz);
free(dPxy);
free(dPxz);
free(dPyz);
free(dV);
free(dcm1);
free(dcm2);
free(dcm3);
free(L);
free(Im);
free(Ic);
free(M);
}
       
 double Pderiv(int h,int k,double vol,double *D,double *E,double *F,double *dP)
 {
     double P;
     double x0,y0,z0,x1,y1,z1,x2,y2,z2;
     x0=D[0];
       y0=D[1];
       z0=D[2];
       x1=E[0];
       y1=E[1];
       z1=E[2];
       x2=F[0];
       y2=F[1];
       z2=F[2];
     P=vol/20.0*(2*D[h]*D[k]+2*E[h]*E[k]+2*F[h]*F[k]+D[h]*E[k]+D[k]*E[h]+D[h]*F[k]+D[k]*F[h]+E[h]*F[k]+E[k]*F[h]);
     if(h==0 && k==0)
     {
         dP[0]=(pow(x2,2)*(-(y1*z0) + y0*z1 - y2*z1 + y1*z2) + pow(x1,2)*(y2*z0 - y2*z1 - y0*z2 + y1*z2) + x1*x2*(-(y1*z0) + y2*z0 + y0*z1 - y2*z1 - y0*z2 + y1*z2) + pow(x0,2)*(-3*y2*z1 + 3*y1*z2)-2*x0*(x2*(y1*z0 - y0*z1 + y2*z1 - y1*z2) + x1*(-(y2*z0) + y2*z1 + y0*z2 - y1*z2)))/60.0;
         dP[1]=(3*pow(x1,2)*(y2*z0 - y0*z2) + pow(x2,2)*(-(y1*z0) + y2*z0 + y0*z1 - y0*z2) - 2*x1*x2*(y1*z0 - y2*z0 - y0*z1 + y0*z2) + pow(x0,2)*(y2*(z0 - z1) + (-y0 + y1)*z2) + x0*(2*x1*(y2*z0 - y2*z1 - y0*z2 + y1*z2) + x2*(-(y1*z0) + y2*z0 + y0*z1 - y2*z1 - y0*z2 + y1*z2)))/60.0;
         dP[2]=(3*pow(x2,2)*(-(y1*z0) + y0*z1) + pow(x1,2)*(-(y1*z0) + y2*z0 + y0*z1 - y0*z2) - 2*x1*x2*(y1*z0 - y2*z0 - y0*z1 + y0*z2) + pow(x0,2)*((y0 - y2)*z1 + y1*(-z0 + z2)) + x0*(2*x2*(-(y1*z0) + y0*z1 - y2*z1 + y1*z2) + x1*(-(y1*z0) + y2*z0 + y0*z1 - y2*z1 - y0*z2 + y1*z2)))/60.0;
     //wrt y
         dP[3]=-((pow(x0,2) + pow(x1,2) + x1*x2 + pow(x2,2) + x0*(x1 + x2))*(-(x2*z1) + x1*z2))/60.0;
         dP[4]=((pow(x0,2) + pow(x1,2) + x1*x2 + pow(x2,2) + x0*(x1 + x2))*(-(x2*z0) + x0*z2))/60.0;
         dP[5]=-((pow(x0,2) + pow(x1,2) + x1*x2 + pow(x2,2) + x0*(x1 + x2))*(-(x1*z0) + x0*z1))/60.0;
         
         dP[6]=((pow(x0,2) + pow(x1,2) + x1*x2 + pow(x2,2) + x0*(x1 + x2))*(-(x2*y1) + x1*y2))/60.0;
         dP[7]=-((pow(x0,2) + pow(x1,2) + x1*x2 + pow(x2,2) + x0*(x1 + x2))*(-(x2*y0) + x0*y2))/60.0;
         dP[8]=((pow(x0,2) + pow(x1,2) + x1*x2 + pow(x2,2) + x0*(x1 + x2))*(-(x1*y0) + x0*y1))/60.0;
         return P;
     }
     if(h==1 && k==1)
     {
         dP[0]=((pow(y0,2) + pow(y1,2) + y1*y2 + pow(y2,2) + y0*(y1 + y2))*(-(y2*z1) + y1*z2))/60.0;
         dP[1]=-((pow(y0,2) + pow(y1,2) + y1*y2 + pow(y2,2) + y0*(y1 + y2))*(-(y2*z0) + y0*z2))/60.0;
         dP[2]=((pow(y0,2) + pow(y1,2) + y1*y2 + pow(y2,2) + y0*(y1 + y2))*(-(y1*z0) + y0*z1))/60.0;
         
         dP[3]=(x2*(3*pow(y0,2)*z1 + pow(y2,2)*z1 + pow(y1,2)*(-z0 + z1) + y1*y2*(-z0 + z1) + y0*(-2*y1*z0 + 2*y1*z1 + 2*y2*z1)) + x0*(2*y0 + y1 + y2)*(-(y2*z1) + y1*z2) - x1*(3*pow(y0,2)*z2 + pow(y1,2)*z2 + y1*y2*(-z0 + z2) + pow(y2,2)*(-z0 + z2) + 2*y0*(y1*z2 + y2*(-z0 + z2))))/60.0;
         dP[4]=(-(x2*((3*pow(y1,2) + 2*y1*y2 + pow(y2,2))*z0 + pow(y0,2)*(z0 - z1) + y0*(2*y1 + y2)*(z0 - z1))) - x1*(y0 + 2*y1 + y2)*(-(y2*z0) + y0*z2) + x0*(pow(y0,2)*z2 + 3*pow(y1,2)*z2 + 2*y1*y2*(-z1 + z2) + pow(y2,2)*(-z1 + z2) + y0*(2*y1*z2 + y2*(-z1 + z2))))/60.0;
         dP[5]=(x2*(y0 + y1 + 2*y2)*(-(y1*z0) + y0*z1) + x1*((pow(y1,2) + 2*y1*y2 + 3*pow(y2,2))*z0 + pow(y0,2)*(z0 - z2) + y0*(y1 + 2*y2)*(z0 - z2)) - x0*(pow(y0,2)*z1 + 3*pow(y2,2)*z1 + y0*(2*y2*z1 + y1*(z1 - z2)) + pow(y1,2)*(z1 - z2) + 2*y1*y2*(z1 - z2)))/60.0;
         
         dP[6]=-((x2*y1 - x1*y2)*(pow(y0,2) + pow(y1,2) + y1*y2 + pow(y2,2) + y0*(y1 + y2)))/60.0;
         dP[7]=((x2*y0 - x0*y2)*(pow(y0,2) + pow(y1,2) + y1*y2 + pow(y2,2) + y0*(y1 + y2)))/60.0;
         dP[8]=((-(x1*y0) + x0*y1)*(pow(y0,2) + pow(y1,2) + y1*y2 + pow(y2,2) + y0*(y1 + y2)))/60.0;
         return P;
     }
     if(h==2 && k==2)
     {
         dP[0]=-((y2*z1 - y1*z2)*(pow(z0,2) + pow(z1,2) + z1*z2 + pow(z2,2) + z0*(z1 + z2)))/60.0;
         dP[1]=((y2*z0 - y0*z2)*(pow(z0,2) + pow(z1,2) + z1*z2 + pow(z2,2) + z0*(z1 + z2)))/60.0;
         dP[2]=((-(y1*z0) + y0*z1)*(pow(z0,2) + pow(z1,2) + z1*z2 + pow(z2,2) + z0*(z1 + z2)))/60.0;
         
         dP[3]=((x2*z1 - x1*z2)*(pow(z0,2) + pow(z1,2) + z1*z2 + pow(z2,2) + z0*(z1 + z2)))/60.0;
         dP[4]=-((x2*z0 - x0*z2)*(pow(z0,2) + pow(z1,2) + z1*z2 + pow(z2,2) + z0*(z1 + z2)))/60.0;
         dP[5]=((x1*z0 - x0*z1)*(pow(z0,2) + pow(z1,2) + z1*z2 + pow(z2,2) + z0*(z1 + z2)))/60.0;
         
         dP[6]=(-(x0*(2*z0 + z1 + z2)*(y2*z1 - y1*z2)) - x2*(-(y0*z1*(2*z0 + z1 + z2)) + y1*(3*pow(z0,2) + pow(z1,2) + z1*z2 + pow(z2,2) + 2*z0*(z1 + z2))) + x1*(-(y0*z2*(2*z0 + z1 + z2)) + y2*(3*pow(z0,2) + pow(z1,2) + z1*z2 + pow(z2,2) + 2*z0*(z1 + z2))))/60.0;
         dP[7]=(x1*(z0 + 2*z1 + z2)*(y2*z0 - y0*z2) + x2*(-(y1*z0*(z0 + 2*z1 + z2)) + y0*(pow(z0,2) + 2*z0*z1 + 3*pow(z1,2) + z0*z2 + 2*z1*z2 + pow(z2,2))) - x0*(-(y1*z2*(z0 + 2*z1 + z2)) + y2*(pow(z0,2) + 2*z0*z1 + 3*pow(z1,2) + z0*z2 + 2*z1*z2 + pow(z2,2))))/60.0;
         dP[8]=(-(x2*(y1*z0 - y0*z1)*(z0 + z1 + 2*z2)) - x1*(-(y2*z0*(z0 + z1 + 2*z2)) + y0*(pow(z0,2) + z0*z1 + pow(z1,2) + 2*z0*z2 + 2*z1*z2 + 3*pow(z2,2))) + x0*(-(y2*z1*(z0 + z1 + 2*z2)) + y1*(pow(z0,2) + z0*z1 + pow(z1,2) + 2*z0*z2 + 2*z1*z2 + 3*pow(z2,2))))/60.0;
         return P;
     }
     if(h==0 && k==1)
     {
         dP[0]=(2*x0*(2*y0 + y1 + y2)*(-(y2*z1) + y1*z2) + x1*(pow(y2,2)*(z0 - z1) + y0*y2*(2*z0 - z1 - z2) - 2*pow(y0,2)*z2 + 2*pow(y1,2)*z2 + y1*y2*(z0 - 2*z1 + z2)) + x2*(2*pow(y0,2)*z1 - 2*pow(y2,2)*z1 - y1*y2*(z0 + z1 - 2*z2) + pow(y1,2)*(-z0 + z2) + y0*y1*(-2*z0 + z1 + z2)))/120.0;
         dP[1]=(-2*x1*(y0 + 2*y1 + y2)*(-(y2*z0) + y0*z2) + x0*(pow(y2,2)*(z0 - z1) + y0*y2*(2*z0 - z1 - z2) - 2*pow(y0,2)*z2 + 2*pow(y1,2)*z2 + y1*y2*(z0 - 2*z1 + z2)) +x2*(2*(-pow(y1,2) + pow(y2,2))*z0 + pow(y0,2)*(z1 - z2) + y0*(y2*(z0 + z1 - 2*z2) - y1*(z0 - 2*z1 + z2))))/120.0;
         dP[2]=(2*x2*(y0 + y1 + 2*y2)*(-(y1*z0) + y0*z1) + x0*(2*pow(y0,2)*z1 - 2*pow(y2,2)*z1 - y1*y2*(z0 + z1 - 2*z2) + pow(y1,2)*(-z0 + z2) + y0*y1*(-2*z0 + z1 + z2)) + x1*(2*(-pow(y1,2) + pow(y2,2))*z0 + pow(y0,2)*(z1 - z2) + y0*(y2*(z0 + z1 - 2*z2) - y1*(z0 - 2*z1 + z2))))/120.0;
         
         dP[3]=(pow(x2,2)*(2*(y0 + y2)*z1 + y1*(-z0 + z1)) + pow(x0,2)*(-2*y2*z1 + 2*y1*z2) + pow(x1,2)*(y2*(z0 - z2) - 2*(y0 + y1)*z2) + x1*x2*(y2*(z0 + z1 - 2*z2) + 2*y0*(z1 - z2) -y1*(z0 - 2*z1 + z2)) + x0*(x1*(y2*(2*z0 - z1 - z2) - 4*y0*z2) + x2*(4*y0*z1 + y1*(-2*z0 + z1 + z2))))/120.0;
         dP[4]=(-(pow(x2,2)*(2*(y1 + y2)*z0 + y0*(z0 - z1))) + 2*pow(x1,2)*(y2*z0 - y0*z2) - x1*x2*(4*y1*z0 + y0*(z0 - 2*z1 + z2)) + pow(x0,2)*(2*(y0 + y1)*z2 + y2*(-z1 + z2)) + x0*(x1*(4*y1*z2 + y2*(z0 - 2*z1 + z2)) + x2*(-(y2*(z0 + z1 - 2*z2)) - 2*y1*(z0 - z2) + y0*(-2*z0 + z1 + z2))))/120.0;
         dP[5]=(2*pow(x2,2)*(-(y1*z0) + y0*z1) + x1*x2*(4*y2*z0 + y0*(z0 + z1 - 2*z2)) + pow(x1,2)*(2*(y1 + y2)*z0 + y0*(z0 - z2)) - pow(x0,2)*(2*y0*z1 + 2*y2*z1 + y1*(z1 - z2)) + x0*(-(x2*(4*y2*z1 + y1*(z0 + z1 - 2*z2))) + x1*(2*y2*(z0 - z1) + y0*(2*z0 - z1 - z2) + y1*(z0 - 2*z1 + z2))))/120.0;
         
         dP[6]=-((x2*y1 - x1*y2)*(x0*(2*y0 + y1 + y2) + x1*(y0 + 2*y1 + y2) + x2*(y0 + y1 + 2*y2)))/120.0;
         dP[7]=((x2*y0 - x0*y2)*(x0*(2*y0 + y1 + y2) + x1*(y0 + 2*y1 + y2) + x2*(y0 + y1 + 2*y2)))/120.0;
         dP[8]=-((x1*y0 - x0*y1)*(x0*(2*y0 + y1 + y2) + x1*(y0 + 2*y1 + y2) + x2*(y0 + y1 + 2*y2)))/120.0;
         return P;
     }
     if(h==0 && k==2)
     {
         dP[0]=(-2*x0*(2*z0 + z1 + z2)*(y2*z1 - y1*z2) + x1*(y2*(z0 - z1)*(2*z0 + 2*z1 + z2) + z2*(-(y0*(2*z0 + z1 + z2)) + y1*(z0 + 2*z1 + z2))) + x2*(-(y1*(z0 - z2)*(2*z0 + z1 + 2*z2)) + z1*(y0*(2*z0 + z1 + z2) - y2*(z0 + z1 + 2*z2))))/120.0;
         dP[1]=(2*x1*(z0 + 2*z1 + z2)*(y2*z0 - y0*z2) + x2*(-(y1*z0*(z0 + 2*z1 + z2)) + y2*z0*(z0 + z1 + 2*z2) + y0*(z1 - z2)*(z0 + 2*(z1 + z2))) + x0*(y2*(z0 - z1)*(2*z0 + 2*z1 + z2) + z2*(-(y0*(2*z0 + z1 + z2)) + y1*(z0 + 2*z1 + z2))))/120.0;
         dP[2]=(-2*x2*(y1*z0 - y0*z1)*(z0 + z1 + 2*z2) + x1*(-(y1*z0*(z0 + 2*z1 + z2)) + y2*z0*(z0 + z1 + 2*z2) + y0*(z1 - z2)*(z0 + 2*(z1 + z2))) + x0*(-(y1*(z0 - z2)*(2*z0 + z1 + 2*z2)) + z1*(y0*(2*z0 + z1 + z2) - y2*(z0 + z1 + 2*z2))))/120.0;
         
         dP[3]=((x2*z1 - x1*z2)*(x0*(2*z0 + z1 + z2) + x1*(z0 + 2*z1 + z2) + x2*(z0 + z1 + 2*z2)))/120.0;
         dP[4]=-((x2*z0 - x0*z2)*(x0*(2*z0 + z1 + z2) + x1*(z0 + 2*z1 + z2) + x2*(z0 + z1 + 2*z2)))/120.0;
         dP[5]=((x1*z0 - x0*z1)*(x0*(2*z0 + z1 + z2) + x1*(z0 + 2*z1 + z2) + x2*(z0 + z1 + 2*z2)))/120.0;
         
         dP[6]=(pow(x0,2)*(-2*y2*z1 + 2*y1*z2) + pow(x1,2)*(-(y0*z2) + y2*(2*z0 + 2*z1 + z2)) - pow(x2,2)*(-(y0*z1) + y1*(2*z0 + z1 + 2*z2)) + x1*x2*(y0*(z1 - z2) - y1*(2*z0 + 2*z1 + z2) + y2*(2*z0 + z1 + 2*z2)) + x0*(-(x2*((-2*y0 + y2)*z1 + y1*(4*z0 + z1))) + x1*((-2*y0 + y1)*z2 + y2*(4*z0 + z2))))/120.0;
         dP[7]=(x1*x2*((-2*y1 + y2)*z0 + y0*(z0 + 4*z1)) + 2*pow(x1,2)*(y2*z0 - y0*z2) - pow(x0,2)*(-(y1*z2) + y2*(2*z0 + 2*z1 + z2)) + pow(x2,2)*(-(y1*z0) + y0*(z0 + 2*(z1 + z2))) + x0*(-(x1*((y0 - 2*y1)*z2 + y2*(4*z1 + z2))) + x2*(y1*(-z0 + z2) + y0*(2*z0 + 2*z1 + z2) - y2*(z0 + 2*(z1 + z2)))))/120.0;
         dP[8]=(2*pow(x2,2)*(-(y1*z0) + y0*z1) + pow(x0,2)*(-(y2*z1) + y1*(2*z0 + z1 + 2*z2)) - x1*x2*((y1 - 2*y2)*z0 + y0*(z0 + 4*z2)) - pow(x1,2)*(-(y2*z0) + y0*(z0 + 2*(z1 + z2))) + x0*(x2*(y0*z1 - 2*y2*z1 + y1*(z1 + 4*z2)) + x1*(y2*(z0 - z1) - y0*(2*z0 + z1 + 2*z2) + y1*(z0 + 2*(z1 + z2)))))/120.0;
         return P;
     }
     if(h==1 && k==2)
     {
         dP[0]=-((y2*z1 - y1*z2)*(y0*(2*z0 + z1 + z2) + y1*(z0 + 2*z1 + z2) + y2*(z0 + z1 + 2*z2)))/120.0;
         dP[1]=((y2*z0 - y0*z2)*(y0*(2*z0 + z1 + z2) + y1*(z0 + 2*z1 + z2) + y2*(z0 + z1 + 2*z2)))/120.0;
         dP[2]=-((y1*z0 - y0*z1)*(y0*(2*z0 + z1 + z2) + y1*(z0 + 2*z1 + z2) + y2*(z0 + z1 + 2*z2)))/120.0;
         
         dP[3]=(-(x0*(2*z0 + z1 + z2)*(y2*z1 - y1*z2)) + x1*(y2*(z0 - z2)*(2*z0 + z1 + 2*z2) - z2*(2*y0*(2*z0 + z1 + z2) + y1*(z0 + 2*z1 + z2))) + x2*(-(y1*(z0 - z1)*(2*z0 + 2*z1 + z2)) + z1*(2*y0*(2*z0 + z1 + z2) + y2*(z0 + z1 + 2*z2))))/120.0;
         dP[4]=(x1*(z0 + 2*z1 + z2)*(y2*z0 - y0*z2) + x0*(-(y2*(z1 - z2)*(z0 + 2*(z1 + z2))) + z2*(y0*(2*z0 + z1 + z2) + 2*y1*(z0 + 2*z1 + z2))) - x2*(y0*(z0 - z1)*(2*z0 + 2*z1 + z2) + z0*(2*y1*(z0 + 2*z1 + z2) + y2*(z0 + z1 + 2*z2))))/120.0;
         dP[5]=(-(x2*(y1*z0 - y0*z1)*(z0 + z1 + 2*z2)) - x0*(y0*z1*(2*z0 + z1 + z2) + 2*y2*z1*(z0 + z1 + 2*z2) + y1*(z1 - z2)*(z0 + 2*(z1 + z2))) + x1*(y0*(z0 - z2)*(2*z0 + z1 + 2*z2) + z0*(y1*(z0 + 2*z1 + z2) + 2*y2*(z0 + z1 + 2*z2))))/120.0;
         
         dP[6]=(x0*(2*y0 + y1 + y2)*(-(y2*z1) + y1*z2) + x2*(2*pow(y0,2)*z1 + y0*(y2*z1 - y1*(4*z0 + z2)) - y1*(y1*(2*z0 + 2*z1 + z2) + y2*(2*z0 + z1 + 2*z2))) + x1*(-2*pow(y0,2)*z2 + y0*(y2*(4*z0 + z1) - y1*z2) + y2*(y1*(2*z0 + 2*z1 + z2) + y2*(2*z0 + z1 + 2*z2))))/120.0;
         dP[7]=(-(x1*(y0 + 2*y1 + y2)*(-(y2*z0) + y0*z2)) - x0*(y1*y2*(z0 + 4*z1) - 2*pow(y1,2)*z2 + pow(y2,2)*(z0 + 2*(z1 + z2)) + y0*(-(y1*z2) + y2*(2*z0 + 2*z1 + z2))) + x2*(-(y1*(2*y1 + y2)*z0) + pow(y0,2)*(2*z0 + 2*z1 + z2) + y0*(y1*(4*z1 + z2) + y2*(z0 + 2*(z1 + z2)))))/120.0;
         dP[8]=(x2*(y0 + y1 + 2*y2)*(-(y1*z0) + y0*z1) + x0*(-2*pow(y2,2)*z1 + y1*y2*(z0 + 4*z2) + pow(y1,2)*(z0 + 2*(z1 + z2)) + y0*(-(y2*z1) + y1*(2*z0 + z1 + 2*z2))) - x1*(-(y2*(y1 + 2*y2)*z0) + pow(y0,2)*(2*z0 + z1 + 2*z2) + y0*(y2*(z1 + 4*z2) + y1*(z0 + 2*(z1 + z2)))))/120.0;
         return P;
     }
     
    
 }
