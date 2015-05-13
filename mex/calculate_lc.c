#include"prepare.h"

void mult_vector(double A[3][3],double *y,double *x)
{
  x[0]=A[0][0]*y[0]+A[0][1]*y[1]+A[0][2]*y[2];
  x[1]=A[1][0]*y[0]+A[1][1]*y[1]+A[1][2]*y[2];
  x[2]=A[2][0]*y[0]+A[2][1]*y[1]+A[2][2]*y[2];
}
double sum_vector(double *vec1,int n)
{
  double res=0;
  for(int j=0;j<n;j++)
    res=res+vec1[j];
  return res;
}
    
  

void calculate_lcurve(int nE,int numfac,int numvert,int **tlist,double **vlist,double *angles,double **Eo,double **E0o,double *TIME,double *bright,double *dbrightx,double *dbrighty,double *dbrightz,double *dbrightb,double *dbrightl,double *dbrighto)
{
  /*First calculate normal,centroid,visibility*/
  /*REMEMBER ROTATION*/
  //   double normal[numfac][3],centroid[numfac][3],E[nE][3],E0[nE][3];
  double **normal,**centroid,**E,**E0;
 double **dSx,**dSy,**dSz;
  normal=mxCalloc(numfac,sizeof(double*));
  centroid=mxCalloc(numfac,sizeof(double*));
  E=mxCalloc(nE,sizeof(double*));
  E0=mxCalloc(nE,sizeof(double*));
  
  
  for(int i=0;i<numfac;i++)
  {
    normal[i]=mxCalloc(3,sizeof(double));
    centroid[i]=mxCalloc(3,sizeof(double));
  }
  for(int i=0;i<nE;i++)
  {
    E[i]=mxCalloc(3,sizeof(double));
    E0[i]=mxCalloc(3,sizeof(double));
  }
  double M[3][3],dMb[3][3],dMo[3][3],dMl[3][3];
  double *vect;
  vect=mxCalloc(3,sizeof(double));
  int *nbl;
  int *ibll;
  int *visiblel;
  
  double beta,lambda,omega;
  double mu,mu0;
  double *tb;
  
  double c=0.1;
  double In,dIx1,dIx2,dIx3,dIy1,dIy2,dIy3,dIz1,dIz2,dIz3;
  double x1,y1,z1,x2,y2,z2,x3,y3,z3;
  double v1[3],v2[3],v3[3],cv[3],normc;
  double dir[3],dir0[3],n[3];
  double dcx1[3],dcx2[3],dcx3[3],dcy1[3],dcy2[3],dcy3[3],dcz1[3],dcz2[3],dcz3[3];
  double dmux1,dmux2,dmux3,dmuy1,dmuy2,dmuy3,dmuz1,dmuz2,dmuz3;
  double dmu0x1,dmu0x2,dmu0x3,dmu0y1,dmu0y2,dmu0y3,dmu0z1,dmu0z2,dmu0z3;
  double *area;
//   double dareax[numfac][3],dareay[numfac][3],dareaz[numfac][3];
   double **dareax,**dareay,**dareaz;
  dareax=mxCalloc(numfac,sizeof(double*));
  dareay=mxCalloc(numfac,sizeof(double*));
  dareaz=mxCalloc(numfac,sizeof(double*));
  for(int i=0;i<numfac;i++)
  {
    dareax[i]=mxCalloc(3,sizeof(double));
    dareay[i]=mxCalloc(3,sizeof(double));
    dareaz[i]=mxCalloc(3,sizeof(double));
  }
//   double dnormalx1[numfac][3],dnormalx2[numfac][3],dnormalx3[numfac][3];
//   double dnormaly1[numfac][3],dnormaly2[numfac][3],dnormaly3[numfac][3];
//   double dnormalz1[numfac][3],dnormalz2[numfac][3],dnormalz3[numfac][3];
  double dnx1[3],dnx2[3],dnx3[3],dny1[3],dny2[3],dny3[3],dnz1[3],dnz2[3],dnz3[3];
  double dAx1,dAx2,dAx3,dAy1,dAy2,dAy3,dAz1,dAz2,dAz3;
  double dt1,dt2;
 

  int vert1,vert2,vert3;
//   double dSa[nE][lindex],dSb[nE][lindex],dSc[nE][lindex];
  
  dSx=mxCalloc(nE,sizeof(double*));
  dSy=mxCalloc(nE,sizeof(double*));
  dSz=mxCalloc(nE,sizeof(double*));
  for(int i=0;i<nE;i++)
  {
    dSx[i]=mxCalloc(numvert,sizeof(double));
    dSy[i]=mxCalloc(numvert,sizeof(double));
    dSz[i]=mxCalloc(numvert,sizeof(double));
  }
  double dub,dul,duo,du0b,du0l,du0o;
  double *dBb,*dBl,*dBo;
  dBb=mxCalloc(nE,sizeof(double));
  dBl=mxCalloc(nE,sizeof(double));
  dBo=mxCalloc(nE,sizeof(double));
  
  double sumbright=0;
  double suma,sumb,sumc;
  
  
 
   double **dnormalx1,**dnormalx2,**dnormalx3;
  double **dnormaly1,**dnormaly2,**dnormaly3;
  double **dnormalz1,**dnormalz2,**dnormalz3;
  dnormalx1=mxCalloc(numfac,sizeof(double*));
  dnormalx2=mxCalloc(numfac,sizeof(double*));
  dnormalx3=mxCalloc(numfac,sizeof(double*));
  dnormaly1=mxCalloc(numfac,sizeof(double*));
  dnormaly2=mxCalloc(numfac,sizeof(double*));
  dnormaly3=mxCalloc(numfac,sizeof(double*));
  dnormalz1=mxCalloc(numfac,sizeof(double*));
  dnormalz2=mxCalloc(numfac,sizeof(double*));
  dnormalz3=mxCalloc(numfac,sizeof(double*));
  for(int i=0;i<numfac;i++)
  {
    dnormalx1[i]=mxCalloc(3,sizeof(double));
    dnormalx2[i]=mxCalloc(3,sizeof(double));
    dnormalx3[i]=mxCalloc(3,sizeof(double));
    dnormaly1[i]=mxCalloc(3,sizeof(double));
    dnormaly2[i]=mxCalloc(3,sizeof(double));
    dnormaly3[i]=mxCalloc(3,sizeof(double));
    dnormalz1[i]=mxCalloc(3,sizeof(double));
    dnormalz2[i]=mxCalloc(3,sizeof(double));
    dnormalz3[i]=mxCalloc(3,sizeof(double));
  }
 
 

  area=mxCalloc(numfac,sizeof(double));
  tb=mxCalloc(nE,sizeof(double));
  visiblel=mxCalloc(nE*numfac,sizeof(int));
  ibll=mxCalloc(numfac*numfac,sizeof(int));
  nbl=mxCalloc(numfac,sizeof(int));
  double omega0;
  beta=angles[0];
  lambda=angles[1];
  omega=angles[2];
  omega0=angles[3];
 
  /*Rotate view and ill vectors*/
  for(int j=0;j<nE;j++)
  {
    rotate(beta,lambda,omega,omega0,TIME[j],M,dMb,dMl,dMo);
    mult_vector(M,Eo[j],E[j]);
    mult_vector(M,E0o[j],E0[j]);
  }
  
  prepare(numfac,numvert,vlist,tlist,normal,centroid,nbl,ibll);
  findblockers(numfac,numvert,nE,E,E0,normal,centroid,nbl,ibll,vlist,tlist,visiblel);
  
  
  
  
  
  
    
    /*Calculate area and derivatives*/
    for(int j=0;j<numfac;j++)
    {
      x1=vlist[tlist[j][0]-1][0];
      y1=vlist[tlist[j][0]-1][1];
      z1=vlist[tlist[j][0]-1][2];
      
      x2=vlist[tlist[j][1]-1][0];
      y2=vlist[tlist[j][1]-1][1];
      z2=vlist[tlist[j][1]-1][2];
      
      x3=vlist[tlist[j][2]-1][0];
      y3=vlist[tlist[j][2]-1][1];
      z3=vlist[tlist[j][2]-1][2];
      
      v1[0]=x2-x1;
      v1[1]=y2-y1;
      v1[2]=z2-z1;
      
      v2[0]=x3-x1;
      v2[1]=y3-y1;
      v2[2]=z3-z1;
      
      cross(v1,v2,cv);
      normc=NORM(cv);
      
      
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
      
      dareax[j][0]=1/(2*normc)*(cv[0]*dcx1[0]+cv[1]*dcx1[1]+cv[2]*dcx1[2]);
      dareax[j][1]=1/(2*normc)*(cv[0]*dcx2[0]+cv[1]*dcx2[1]+cv[2]*dcx2[2]);
      dareax[j][2]=1/(2*normc)*(cv[0]*dcx3[0]+cv[1]*dcx3[1]+cv[2]*dcx3[2]);
      
      dareay[j][0]=1/(2*normc)*(cv[0]*dcy1[0]+cv[1]*dcy1[1]+cv[2]*dcy1[2]);
      dareay[j][1]=1/(2*normc)*(cv[0]*dcy2[0]+cv[1]*dcy2[1]+cv[2]*dcy2[2]);
      dareay[j][2]=1/(2*normc)*(cv[0]*dcy3[0]+cv[1]*dcy3[1]+cv[2]*dcy3[2]);
      
      dareaz[j][0]=1/(2*normc)*(cv[0]*dcz1[0]+cv[1]*dcz1[1]+cv[2]*dcz1[2]);
      dareaz[j][1]=1/(2*normc)*(cv[0]*dcz2[0]+cv[1]*dcz2[1]+cv[2]*dcz2[2]);
      dareaz[j][2]=1/(2*normc)*(cv[0]*dcz3[0]+cv[1]*dcz3[1]+cv[2]*dcz3[2]);
      
      area[j]=0.5*normc; 
      for(int jk=0;jk<3;jk++)
      {
	dnormalx1[j][jk]=(normc*dcx1[jk]-cv[jk]*2*dareax[j][0])/(normc*normc);
	dnormalx2[j][jk]=(normc*dcx2[jk]-cv[jk]*2*dareax[j][1])/(normc*normc);
	dnormalx3[j][jk]=(normc*dcx3[jk]-cv[jk]*2*dareax[j][2])/(normc*normc);
	
	dnormaly1[j][jk]=(normc*dcy1[jk]-cv[jk]*2*dareay[j][0])/(normc*normc);
	dnormaly2[j][jk]=(normc*dcy2[jk]-cv[jk]*2*dareay[j][1])/(normc*normc);
	dnormaly3[j][jk]=(normc*dcy3[jk]-cv[jk]*2*dareay[j][2])/(normc*normc);
	
	dnormalz1[j][jk]=(normc*dcz1[jk]-cv[jk]*2*dareaz[j][0])/(normc*normc);
	dnormalz2[j][jk]=(normc*dcz2[jk]-cv[jk]*2*dareaz[j][1])/(normc*normc);
	dnormalz3[j][jk]=(normc*dcz3[jk]-cv[jk]*2*dareaz[j][2])/(normc*normc);
      }
      
    }
    
    
    for(int e=0;e<nE;e++)
    {
     
      
      
      rotate(beta,lambda,omega,omega0,TIME[e],M,dMb,dMl,dMo);
      for(int j=0;j<numfac;j++)
      {
	if(visiblel[e*numfac+j]==0) /*Not visible, skip*/
	  continue;
	
	for(int i=0;i<3;i++)
	{
	  dir[i]=E[e][i];
	  dir0[i]=E0[e][i];
	  n[i]=normal[j][i];
	  dnx1[i]=dnormalx1[j][i];
	  dnx2[i]=dnormalx2[j][i];
	  dnx3[i]=dnormalx3[j][i];
	  dny1[i]=dnormaly1[j][i];
	  dny2[i]=dnormaly2[j][i];
	  dny3[i]=dnormaly3[j][i];
	  dnz1[i]=dnormalz1[j][i];
	  dnz2[i]=dnormalz2[j][i];
	  dnz3[i]=dnormalz3[j][i];
	}
	x1=vlist[tlist[j][0]-1][0];
	y1=vlist[tlist[j][0]-1][1];
	z1=vlist[tlist[j][0]-1][2];
	
	x2=vlist[tlist[j][1]-1][0];
	y2=vlist[tlist[j][1]-1][1];
	z2=vlist[tlist[j][1]-1][2];
	
	x3=vlist[tlist[j][2]-1][0];
	y3=vlist[tlist[j][2]-1][1];
	z3=vlist[tlist[j][2]-1][2];
	
	dAx1=dareax[j][0];
	dAx2=dareax[j][1];
	dAx3=dareax[j][2];
	
	dAy1=dareay[j][0];
	dAy2=dareay[j][1];
	dAy3=dareay[j][2];
	
	dAz1=dareaz[j][0];
	dAz2=dareaz[j][1];
	dAz3=dareaz[j][2];
	
	dmux1=DOT(dnx1,dir);
	dmux2=DOT(dnx2,dir);
	dmux3=DOT(dnx3,dir);
	
	dmuy1=DOT(dny1,dir);
	dmuy2=DOT(dny2,dir);
	dmuy3=DOT(dny3,dir);
	
	dmuz1=DOT(dnz1,dir);
	dmuz2=DOT(dnz2,dir);
	dmuz3=DOT(dnz3,dir);
	
	dmu0x1=DOT(dnx1,dir0);
	dmu0x2=DOT(dnx2,dir0);
	dmu0x3=DOT(dnx3,dir0);
	
	dmu0y1=DOT(dny1,dir0);
	dmu0y2=DOT(dny2,dir0);
	dmu0y3=DOT(dny3,dir0);
	
	dmu0z1=DOT(dnz1,dir0);
	dmu0z2=DOT(dnz2,dir0);
	dmu0z3=DOT(dnz3,dir0);
	vert1=tlist[j][0]-1;
	vert2=tlist[j][1]-1;
	vert3=tlist[j][2]-1;
	
	
	mu=DOT(dir,n);
	mu0=DOT(dir0,n);
	/*Brightness*/
	In=mu*mu0*(1/(mu+mu0)+c);
	tb[e]=tb[e]+In*area[j];
	
	dt1=pow(mu/(mu+mu0),2)+c*mu;
	dt2=pow(mu0/(mu+mu0),2)+c*mu0;
	
	dIx1=dt1*dmu0x1+dt2*dmux1;
	dIx2=dt1*dmu0x2+dt2*dmux2;
	dIx3=dt1*dmu0x3+dt2*dmux3;
	
	dIy1=dt1*dmu0y1+dt2*dmuy1;
	dIy2=dt1*dmu0y2+dt2*dmuy2;
	dIy3=dt1*dmu0y3+dt2*dmuy3;
	
	dIz1=dt1*dmu0z1+dt2*dmuz1;
	dIz2=dt1*dmu0z2+dt2*dmuz2;
	dIz3=dt1*dmu0z3+dt2*dmuz3;
	//Derivative wrt vertices
	
	dSx[e][vert1]=dSx[e][vert1]+In*dAx1+dIx1*area[j];
	dSx[e][vert2]=dSx[e][vert2]+In*dAx2+dIx2*area[j];
	dSx[e][vert3]=dSx[e][vert3]+In*dAx3+dIx3*area[j];
	
	dSy[e][vert1]=dSy[e][vert1]+In*dAy1+dIy1*area[j];
	dSy[e][vert2]=dSy[e][vert2]+In*dAy2+dIy2*area[j];
	dSy[e][vert3]=dSy[e][vert3]+In*dAy3+dIy3*area[j];
	
	dSz[e][vert1]=dSz[e][vert1]+In*dAz1+dIz1*area[j];
	dSz[e][vert2]=dSz[e][vert2]+In*dAz2+dIz2*area[j];
	dSz[e][vert3]=dSz[e][vert3]+In*dAz3+dIz3*area[j];
	
	
	mult_vector(dMb,Eo[e],vect);
	dub=DOT(vect,n);
	mult_vector(dMl,Eo[e],vect);
	dul=DOT(vect,n);
	mult_vector(dMo,Eo[e],vect);
	duo=DOT(vect,n);
	mult_vector(dMb,E0o[e],vect);
	du0b=DOT(vect,n);
	mult_vector(dMl,E0o[e],vect);
	du0l=DOT(vect,n);
	mult_vector(dMo,E0o[e],vect);
	du0o=DOT(vect,n);
	
	dBb[e]=dBb[e]+area[j]*(dt1*du0b+dt2*dub);
	dBl[e]=dBl[e]+area[j]*(dt1*du0l+dt2*dul);
	dBo[e]=dBo[e]+area[j]*(dt1*du0o+dt2*duo);
	
       } /*End of j, numfac*/
      } /*End of e, nE*/
    
     
      sumbright=sum_vector(tb,nE);
      double sbeta=0,slambda=0,somega=0;
      for(int j=0;j<nE;j++)
      {
	bright[j]=nE*tb[j]/sumbright;
	sbeta=sbeta+dBb[j];
	slambda=slambda+dBl[j];
	somega=somega+dBo[j];
      }
      
      
      for(int sh_index=0;sh_index<numvert;sh_index++)
      {
	suma=0;
	sumb=0;
	sumc=0;
	
	for(int k=0;k<nE;k++)
	{
	  suma=suma+dSx[k][sh_index];
	  sumb=sumb+dSy[k][sh_index];
	  sumc=sumc+dSz[k][sh_index];
	}
	for(int i=0;i<nE;i++)
	{
	  // 	    
	  dbrightx[i*(numvert)+sh_index]=nE*(sumbright*dSx[i][sh_index]-tb[i]*suma)/pow(sumbright,2);
	  dbrighty[i*(numvert)+sh_index]=nE*(sumbright*dSy[i][sh_index]-tb[i]*sumb)/pow(sumbright,2);
	  dbrightz[i*(numvert)+sh_index]=nE*(sumbright*dSz[i][sh_index]-tb[i]*sumc)/pow(sumbright,2);
	}
      }
      
      
      for(int i=0;i<nE;i++)
      {
	dbrightb[i]=nE*(sumbright*dBb[i]-tb[i]*sbeta)/pow(sumbright,2);
	dbrightl[i]=nE*(sumbright*dBl[i]-tb[i]*slambda)/pow(sumbright,2);
	dbrighto[i]=nE*(sumbright*dBo[i]-tb[i]*somega)/pow(sumbright,2);
      }
     
}

      
      
      
      
      
      
      
      