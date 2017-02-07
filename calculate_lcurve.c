#include"utils.h"
#include"globals.h"

double phase_function(double *E,double *E0,double *params,double *dpdp);  
void calculate_lcurve(int *tlist,double *vlist,int numfac,int numvert,double *angles,double *Eo,double *E0o,int nE,double *TIME,double *bright,double *dbrightx,double *dbrighty,double *dbrightz,double *dbrightb,double *dbrightl,double *dbrighto,double *dbrightp,double *A,double *Alimit,double *dA,int rel,double *params)
{
   
  double *normal,*centroid,*E,*E0;
  double phase,dphase[3];
 double la=0;
 double ha=1; /*These are the albedo limits*/
 double alb_term=0;
  double *dSx,*dSy,*dSz;
  int albedo=0;
  normal=malloc(numfac*3*sizeof(double));
  
  E=malloc(3*nE*sizeof(double));
  E0=malloc(3*nE*sizeof(double));
  
  double M[3][3],dMb[3][3],dMo[3][3],dMl[3][3];
  double vect[3];
  double alb=0;
  double *dSdalb=NULL;
 if(A!=NULL)
 {
     albedo=1;
     dSdalb=calloc(nE*numfac,sizeof(double));
     la=Alimit[0];
     ha=Alimit[1];
 }
 
     
  int *visiblel;
  
  double beta,lambda,omega;
  double mu,mu0;
  double *tb;
  
  double c=0.1;
  if(params!=NULL)
      c=params[3];

  double In,dIx1,dIx2,dIx3,dIy1,dIy2,dIy3,dIz1,dIz2,dIz3;
  double x1,y1,z1,x2,y2,z2,x3,y3,z3;
  double v1[3],v2[3],v3[3],cv[3],normc,Scatt;
  double dir[3],dir0[3],n[3];
  double dcx1[3],dcx2[3],dcx3[3],dcy1[3],dcy2[3],dcy3[3],dcz1[3],dcz2[3],dcz3[3];
  double dmux1,dmux2,dmux3,dmuy1,dmuy2,dmuy3,dmuz1,dmuz2,dmuz3;
  double dmu0x1,dmu0x2,dmu0x3,dmu0y1,dmu0y2,dmu0y3,dmu0z1,dmu0z2,dmu0z3;
  double *area;

   double *dareax,*dareay,*dareaz;
  dareax=calloc(numfac*3,sizeof(double));
  dareay=calloc(numfac*3,sizeof(double));
  dareaz=calloc(numfac*3,sizeof(double));
 
  double dnx1[3],dnx2[3],dnx3[3],dny1[3],dny2[3],dny3[3],dnz1[3],dnz2[3],dnz3[3];
  double dAx1,dAx2,dAx3,dAy1,dAy2,dAy3,dAz1,dAz2,dAz3;
  double dt1,dt2;
 

  int vert1,vert2,vert3;

  
  dSx=calloc(nE*numvert,sizeof(double));
  dSy=calloc(nE*numvert,sizeof(double));
  dSz=calloc(nE*numvert,sizeof(double));
  
  double dub,dul,duo,du0b,du0l,du0o;
  double *dBb,*dBl,*dBo;
  dBb=calloc(nE,sizeof(double));
  dBl=calloc(nE,sizeof(double));
  dBo=calloc(nE,sizeof(double));
  
  double sumbright=0;
  double suma,sumb,sumc;
  
  
 
   double *dnormalx1,*dnormalx2,*dnormalx3;
  double *dnormaly1,*dnormaly2,*dnormaly3;
  double *dnormalz1,*dnormalz2,*dnormalz3;
  dnormalx1=calloc(numfac*3,sizeof(double));
  dnormalx2=calloc(numfac*3,sizeof(double));
  dnormalx3=calloc(numfac*3,sizeof(double));
  dnormaly1=calloc(numfac*3,sizeof(double));
  dnormaly2=calloc(numfac*3,sizeof(double));
  dnormaly3=calloc(numfac*3,sizeof(double));
  dnormalz1=calloc(numfac*3,sizeof(double));
  dnormalz2=calloc(numfac*3,sizeof(double));
  dnormalz3=calloc(numfac*3,sizeof(double));
  

  area=calloc(numfac,sizeof(double));
  tb=calloc(nE,sizeof(double));
  visiblel=calloc(nE*numfac,sizeof(int));
 
  double omega0;
  beta=angles[0];
  lambda=angles[1];
  omega=angles[2];
  omega0=angles[3];
 
  /*Rotate view and ill vectors*/
  for(int j=0;j<nE;j++)
  {
    rotate(beta,lambda,omega,omega0,TIME[j],M,dMb,dMl,dMo);
    mult_vector(M,&(Eo[3*j]),&(E[3*j]));
    mult_vector(M,&(E0o[3*j]),&(E0[3*j]));
  }
  
  
  FindActualBlockers(tlist,vlist,numfac,numvert,E,E0,nE,visiblel);
  int vr1,vr2,vr3;
  

    
    /*Calculate area and derivatives*/
    for(int j=0;j<numfac;j++)
    {
        vr1=tlist[3*j]-1; 
        vr2=tlist[3*j+1]-1;
        vr3=tlist[3*j+2]-1;
      x1=vlist[3*vr1];
      y1=vlist[3*vr1+1];
      z1=vlist[3*vr1+2];
      
      x2=vlist[3*vr2];
      y2=vlist[3*vr2+1];
      z2=vlist[3*vr2+2];
      
      x3=vlist[3*vr3];
      y3=vlist[3*vr3+1];
      z3=vlist[3*vr3+2];
      
      v1[0]=x2-x1;
      v1[1]=y2-y1;
      v1[2]=z2-z1;
      
      v2[0]=x3-x1;
      v2[1]=y3-y1;
      v2[2]=z3-z1;
      
      cross(v1,v2,cv);
      normc=NORM(cv);
      normal[3*j+0]=cv[0]/normc;
      normal[3*j+1]=cv[1]/normc;
      normal[3*j+2]=cv[2]/normc;
      
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
      
      dareax[3*j+0]=1/(2*normc)*(cv[0]*dcx1[0]+cv[1]*dcx1[1]+cv[2]*dcx1[2]);
      dareax[3*j+1]=1/(2*normc)*(cv[0]*dcx2[0]+cv[1]*dcx2[1]+cv[2]*dcx2[2]);
      dareax[3*j+2]=1/(2*normc)*(cv[0]*dcx3[0]+cv[1]*dcx3[1]+cv[2]*dcx3[2]);
      
      dareay[3*j+0]=1/(2*normc)*(cv[0]*dcy1[0]+cv[1]*dcy1[1]+cv[2]*dcy1[2]);
      dareay[3*j+1]=1/(2*normc)*(cv[0]*dcy2[0]+cv[1]*dcy2[1]+cv[2]*dcy2[2]);
      dareay[3*j+2]=1/(2*normc)*(cv[0]*dcy3[0]+cv[1]*dcy3[1]+cv[2]*dcy3[2]);
      
      dareaz[3*j+0]=1/(2*normc)*(cv[0]*dcz1[0]+cv[1]*dcz1[1]+cv[2]*dcz1[2]);
      dareaz[3*j+1]=1/(2*normc)*(cv[0]*dcz2[0]+cv[1]*dcz2[1]+cv[2]*dcz2[2]);
      dareaz[3*j+2]=1/(2*normc)*(cv[0]*dcz3[0]+cv[1]*dcz3[1]+cv[2]*dcz3[2]);
      
      area[j]=0.5*normc; 
      
      for(int jk=0;jk<3;jk++)
      {
	dnormalx1[3*j+jk]=(normc*dcx1[jk]-cv[jk]*2*dareax[3*j+0])/(normc*normc);
	dnormalx2[3*j+jk]=(normc*dcx2[jk]-cv[jk]*2*dareax[3*j+1])/(normc*normc);
	dnormalx3[3*j+jk]=(normc*dcx3[jk]-cv[jk]*2*dareax[3*j+2])/(normc*normc);
	
	dnormaly1[3*j+jk]=(normc*dcy1[jk]-cv[jk]*2*dareay[3*j+0])/(normc*normc);
	dnormaly2[3*j+jk]=(normc*dcy2[jk]-cv[jk]*2*dareay[3*j+1])/(normc*normc);
	dnormaly3[3*j+jk]=(normc*dcy3[jk]-cv[jk]*2*dareay[3*j+2])/(normc*normc);
	
	dnormalz1[3*j+jk]=(normc*dcz1[jk]-cv[jk]*2*dareaz[3*j+0])/(normc*normc);
	dnormalz2[3*j+jk]=(normc*dcz2[jk]-cv[jk]*2*dareaz[3*j+1])/(normc*normc);
	dnormalz3[3*j+jk]=(normc*dcz3[jk]-cv[jk]*2*dareaz[3*j+2])/(normc*normc);
      }
      
    }
    
    
    for(int e=0;e<nE;e++)
    {
        for(int i=0;i<3;i++)
        {
        dir[i]=E[3*e+i];
	  dir0[i]=E0[3*e+i];
        }
       
      rotate(beta,lambda,omega,omega0,TIME[e],M,dMb,dMl,dMo);
      for(int j=0;j<numfac;j++)
      {
	if(visiblel[e*numfac+j]==0) /*Not visible, skip*/
	  continue;
	
	for(int i=0;i<3;i++)
	{
	  
	  n[i]=normal[3*j+i];
	  dnx1[i]=dnormalx1[3*j+i];
	  dnx2[i]=dnormalx2[3*j+i];
	  dnx3[i]=dnormalx3[3*j+i];
	  dny1[i]=dnormaly1[3*j+i];
	  dny2[i]=dnormaly2[3*j+i];
	  dny3[i]=dnormaly3[3*j+i];
	  dnz1[i]=dnormalz1[3*j+i];
	  dnz2[i]=dnormalz2[3*j+i];
	  dnz3[i]=dnormalz3[3*j+i];
	}
	vr1=tlist[3*j]-1; /*LAST MODIFICATION HERE*/
        vr2=tlist[3*j+1]-1;
        vr3=tlist[3*j+2]-1;
	x1=vlist[3*vr1+0];
	y1=vlist[3*vr1+1];
	z1=vlist[3*vr1+2];
	
	x2=vlist[3*vr2+0];
	y2=vlist[3*vr2+1];
	z2=vlist[3*vr2+2];
	
	x3=vlist[3*vr3+0];
	y3=vlist[3*vr3+1];
	z3=vlist[3*vr3+2];
	
	dAx1=dareax[3*j+0];
	dAx2=dareax[3*j+1];
	dAx3=dareax[3*j+2];
	
	dAy1=dareay[3*j+0];
	dAy2=dareay[3*j+1];
	dAy3=dareay[3*j+2];
	
	dAz1=dareaz[3*j+0];
	dAz2=dareaz[3*j+1];
	dAz3=dareaz[3*j+2];
	
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
	vert1=tlist[3*j]-1;
	vert2=tlist[3*j+1]-1;
	vert3=tlist[3*j+2]-1;
	
	
	mu=DOT(dir,n);
	mu0=DOT(dir0,n);
	/*Brightness*/
        //Phase parameters are not NULL, so use them
       if(params!=NULL)
        {
            phase=phase_function(dir,dir0,params,dphase);
            
        }
        else
            phase=1;
        //Scattering Law?
        if(INI_HAPKE!=NULL)
            dhapke_bright(dir,dir0,mu,mu0,INI_HAPKE,INI_HAPKE[4],&Scatt,&dt2,&dt1); //Use Hapke
        else 
            Scatt=mu*mu0*(1.0/(mu+mu0)+c); //Use LSL
      
        //Variable albedo?
        if(albedo==0)
        {
            alb=0;
            alb_term=1;
        }
        else
        {
            alb=A[j];
            dSdalb[e*numfac+j]=phase*(ha-la)*exp(alb)/pow(exp(alb)+1,2)*Scatt*area[j];
            
            alb_term=(la+(ha-la)*exp(alb)/(exp(alb)+1));
        }
	In=alb_term*Scatt;
	tb[e]=tb[e]+phase*In*area[j];
       //printf("Facet %d: alb_term: %f Area: %f  mu: %f mu0: %f phase: %f Total: %f\n",j+1,alb_term,area[j],mu,mu0,phase,phase*In*area[j]);
	//Derivatives wrt phase params
        if(params!=NULL)
        {
            dbrightp[e*4]=dbrightp[e*4]+dphase[0]*In*area[j];
            dbrightp[e*4+1]=dbrightp[e*4+1]+dphase[1]*In*area[j];
            dbrightp[e*4+2]=dbrightp[e*4+2]+dphase[2]*In*area[j];
            dbrightp[e*4+3]=dbrightp[e*4+3]+alb_term*mu*mu0*phase*area[j];
        }
        //If not HAPKE, then calculate derivates of LSL wrt mu,mu0
        if(INI_HAPKE==NULL)
        {
            dt1=phase*alb_term*(pow(mu/(mu+mu0),2)+c*mu);
            dt2=phase*alb_term*(pow(mu0/(mu+mu0),2)+c*mu0);
        }
	
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
	
	dSx[e*numvert+vert1]=dSx[e*numvert+vert1]+phase*In*dAx1+dIx1*area[j];
	dSx[e*numvert+vert2]=dSx[e*numvert+vert2]+phase*In*dAx2+dIx2*area[j];
	dSx[e*numvert+vert3]=dSx[e*numvert+vert3]+phase*In*dAx3+dIx3*area[j];
	
	dSy[e*numvert+vert1]=dSy[e*numvert+vert1]+phase*In*dAy1+dIy1*area[j];
	dSy[e*numvert+vert2]=dSy[e*numvert+vert2]+phase*In*dAy2+dIy2*area[j];
	dSy[e*numvert+vert3]=dSy[e*numvert+vert3]+phase*In*dAy3+dIy3*area[j];
	
	dSz[e*numvert+vert1]=dSz[e*numvert+vert1]+phase*In*dAz1+dIz1*area[j];
	dSz[e*numvert+vert2]=dSz[e*numvert+vert2]+phase*In*dAz2+dIz2*area[j];
	dSz[e*numvert+vert3]=dSz[e*numvert+vert3]+phase*In*dAz3+dIz3*area[j];
	
	
	mult_vector(dMb,&(Eo[3*e]),vect);
	dub=DOT(vect,n);
	mult_vector(dMl,&(Eo[3*e]),vect);
	dul=DOT(vect,n);
	mult_vector(dMo,&(Eo[3*e]),vect);
	duo=DOT(vect,n);
	mult_vector(dMb,&(E0o[3*e]),vect);
	du0b=DOT(vect,n);
	mult_vector(dMl,&(E0o[3*e]),vect);
	du0l=DOT(vect,n);
	mult_vector(dMo,&(E0o[3*e]),vect);
	du0o=DOT(vect,n);
	
	dBb[e]=dBb[e]+area[j]*(dt1*du0b+dt2*dub);
	dBl[e]=dBl[e]+area[j]*(dt1*du0l+dt2*dul);
	dBo[e]=dBo[e]+area[j]*(dt1*du0o+dt2*duo);
	
       } /*End of j, numfac*/
       
      } /*End of e, nE*/
   double sum_paramc=0;
     double sumbright2=0;
      sumbright=sum_vector(tb,nE);
      sumbright2=pow(sumbright,2);
    if(rel!=0)
    {
      double sbeta=0,slambda=0,somega=0;
      for(int j=0;j<nE;j++)
      {
	bright[j]=nE*tb[j]/sumbright;
	sbeta=sbeta+dBb[j];
	slambda=slambda+dBl[j];
	somega=somega+dBo[j];
        sum_paramc+=dbrightp[j*4+3];
      }
      
      
      for(int sh_index=0;sh_index<numvert;sh_index++)
      {
	suma=0;
	sumb=0;
	sumc=0;
	
	for(int k=0;k<nE;k++)
	{
	  suma=suma+dSx[k*numvert+sh_index];
	  sumb=sumb+dSy[k*numvert+sh_index];
	  sumc=sumc+dSz[k*numvert+sh_index];
	}
	for(int i=0;i<nE;i++)
	{
	  // 	    
	  dbrightx[i*(numvert)+sh_index]=nE*(sumbright*dSx[i*numvert+sh_index]-tb[i]*suma)/sumbright2;
	  dbrighty[i*(numvert)+sh_index]=nE*(sumbright*dSy[i*numvert+sh_index]-tb[i]*sumb)/sumbright2;
	  dbrightz[i*(numvert)+sh_index]=nE*(sumbright*dSz[i*numvert+sh_index]-tb[i]*sumc)/sumbright2;
	}
      }
      
      
      for(int i=0;i<nE;i++)
      {
	dbrightb[i]=nE*(sumbright*dBb[i]-tb[i]*sbeta)/sumbright2;
	dbrightl[i]=nE*(sumbright*dBl[i]-tb[i]*slambda)/sumbright2;
	dbrighto[i]=nE*(sumbright*dBo[i]-tb[i]*somega)/sumbright2;
        
      }
      if(params!=NULL)
          for(int j=0;j<nE;j++)
        {
            
            dbrightp[j*4+3]=nE*(dbrightp[j*4+3]*sumbright-tb[j]*sum_paramc)/pow(sumbright,2);
      }
      if(albedo==1)
      {
        for(int i=0;i<numfac;i++)
        {
            suma=0;
            for(int k=0;k<nE;k++)
                suma=suma+dSdalb[k*numfac+i];
            for(int h=0;h<nE;h++)
                dA[h*numfac+i]=nE*(sumbright*dSdalb[h*numfac+i]-tb[h]*suma)/sumbright2;
        }
      }
       
    }
    else
    {
        memcpy(bright,tb,sizeof(double)*nE);
        memcpy(dbrightx,dSx,sizeof(double)*nE*numvert);
        memcpy(dbrighty,dSy,sizeof(double)*nE*numvert);
        memcpy(dbrightz,dSz,sizeof(double)*nE*numvert);
        memcpy(dbrightb,dBb,sizeof(double)*nE);
        memcpy(dbrightl,dBl,sizeof(double)*nE);
        memcpy(dbrighto,dBo,sizeof(double)*nE);
        if(albedo==1)
            memcpy(dA,dSdalb,sizeof(double)*numfac*nE);
    }
    
      free(dSx);
      free(dSy);
      free(dSz);
      free(dBb);
      free(dBl);
      free(dBo);
     free(normal);
     free(E);
     free(E0);
     free(dareax);
     free(dareay);
     free(dareaz);
     free(area);
     free(tb);
     free(dnormalx1);
     free(dnormalx2);
     free(dnormalx3);
     free(dnormaly1);
     free(dnormaly2);
     free(dnormaly3);
     free(dnormalz1);
     free(dnormalz2);
     free(dnormalz3);
     free(visiblel);
     free(dSdalb);
}

      
      
      
      
      
      
      
      
