#include"utils.h"
#define sec(x) (1.0/cos(x))
void dshadow(double mu,double mu0,double tth,double cal,double *,double *,double *,double *,double *,double *,double *,double *,double *);
void dhapke(double *p,double mu,double mu0,double alpha,double *rfh,double *rdfhdmu,double *rdfhdmu0);

void dhapke_bright(double *E,double *E0,double mu,double mu0,double *p,double th,double *rss,double *rdssdmu,double *rdssdmu0)
{
/* 
 * Calculate Hapke brightness function and derivatives wrt mu,mu0.
 */
  
  double tth,cth,cal,alpha;
  tth=tan(PI/180.0*th);
  cth=1.0/sqrt(1+PI*tth);
  cal=DOT(E,E0);
  alpha=acos(cal);
  
  double sh,mueiefi,mu0eiefi,dshdmu,dshdmu0,dmueiefidmu,dmueiefidmu0,dmu0eiefidmu,dmu0eiefidmu0;
  dshadow(mu,mu0,tth,cal,&sh,&mueiefi,&mu0eiefi,&dshdmu,&dshdmu0,&dmueiefidmu,&dmueiefidmu0,&dmu0eiefidmu,&dmu0eiefidmu0);
 
 
  
  double mu0new,dmu0newdmu,dmu0newdmu0,munew,dmunewdmu,dmunewdmu0,dnom,dnomdmu,dnomdmu0,sls,dslsdmu,dslsdmu0;
  double fh,dfhdmu,dfhdmu0,dfhdmu_temp,ss,dssdmu,dssdmu0;
  mu0new=cth*mu0eiefi;
  dmu0newdmu=cth*dmu0eiefidmu;
  dmu0newdmu0=cth*dmu0eiefidmu0;
  munew=cth*mueiefi;
  dmunewdmu=cth*dmueiefidmu;
  dmunewdmu0=cth*dmueiefidmu0;
  dnom=munew+mu0new;
  dnomdmu=dmunewdmu+dmu0newdmu;
  dnomdmu0=dmunewdmu0+dmu0newdmu0;
  sls=mu*mu0new/dnom;
  dslsdmu=((mu0new+mu*dmu0newdmu)*dnom-mu*mu0new*dnomdmu)/pow(dnom,2);
  dslsdmu0=(mu*dmu0newdmu0*dnom-mu*mu0new*dnomdmu0)/pow(dnom,2);
  
  dhapke(p,munew,mu0new,alpha,&fh,&dfhdmu,&dfhdmu0);
 
  
  dfhdmu_temp=dfhdmu;
  dfhdmu=dfhdmu*dmunewdmu+dfhdmu0*dmu0newdmu;
  dfhdmu0=dfhdmu0*dmu0newdmu0+dfhdmu_temp*dmunewdmu0;
  
  ss=fh*sls*sh;
  dssdmu=dfhdmu*sls*sh+fh*dslsdmu*sh+fh*sls*dshdmu;
  dssdmu0=dfhdmu0*sls*sh+fh*dslsdmu0*sh+fh*sls*dshdmu0;
  *rss=ss;
  *rdssdmu=dssdmu;
  *rdssdmu0=dssdmu0;
  
}

void dshadow(double mu,double mu0,double tth,double cal,double *rsh,double *rmueiefi,double *rmu0eiefi,double *rdshdmu,double *rdshdmu0,double *rdmueiefidmu,double *rdmueiefidmu0,double *rdmu0eiefidmu,double *rdmu0eiefidmu0)
{
//[sh,mueiefi,mu0eiefi,dshdmu,dshdmu0,dmueiefidmu,dmueiefidmu0,dmu0eiefidmu,dmu0eiefidmu0]
//Derivatives wrt mu,mu0
double EP2=1e-6;
double i,e;
i=acos(mu0);
e=acos(mu);
//ae and ae cannot be zero
if(i<EP2)
    i=EP2;

if(e<EP2)
    e=EP2;

double innerfi,fi,dedmu,didmu0,dfidmu,dfidmu0;

innerfi=(cal-mu*mu0)/(sin(i)*sin(e));
fi=acos(innerfi);
if(fi<EP2)
  fi=EP2;

if(isnan(fi))
  fi=PI-EP2;
if(innerfi>1-EP2)
  innerfi=1-EP2;
if(innerfi<-1+EP2)
  innerfi=-1+EP2;

dedmu=-1.0/sqrt(1.0-pow(mu,2));
didmu0=-1.0/sqrt(1.0-pow(mu0,2));
dfidmu=-1.0/sqrt(1.0-pow(innerfi,2))*(-mu0*sin(i)*sin(e)-(cal-mu*mu0)*sin(i)*cos(e)*(dedmu))/pow(sin(i)*sin(e),2);
dfidmu0=-1.0/sqrt(1.0-pow(innerfi,2))*(-mu*sin(i)*sin(e)-(cal-mu*mu0)*sin(e)*cos(i)*(didmu0))/pow(sin(i)*sin(e),2);

double f,dfdmu,dfdmu0,E1e,E1i,E2e,E2i;

f=exp(-2.0*tan(0.5*fi));
dfdmu=f*(-2.0*pow(sec(0.5*fi),2))*0.5*dfidmu;
dfdmu0=f*(-2.0*pow(sec(0.5*fi),2))*0.5*dfidmu0;
E1e=exp(-2.0/(tth*tan(e)*PI));
E1i=exp(-2.0/(tth*tan(i)*PI));
E2e=exp(-1.0/(pow(tth*tan(e),2)*PI));
E2i=exp(-1.0/(pow(tth*tan(i),2)*PI));

double dE1edmu,dE1idmu0,dE2edmu,dE2idmu0;
dE1edmu=E1e*(-2.0/(tth*PI)*(-pow(sec(e)/tan(e),2)))*dedmu;
dE1idmu0=E1i*(-2.0/(tth*PI)*(-pow(sec(i)/tan(i),2)))*didmu0;
dE2edmu=E2e*(-1.0/(pow(tth,2)*PI))*(-2*pow(1.0/tan(e),3)*pow(sec(e),2))*dedmu;
dE2idmu0=E2i*(-1.0/(pow(tth,2)*PI))*(-2*pow(1.0/tan(i),3)*pow(sec(i),2))*didmu0;
double top,dtopdmu,dtopdmu0,bot,dbotdmu,dbotdmu0;
double mu0ei0pi,dmu0ei0pidmu0;
double mue0e0,dmue0e0dmu,dmue0e0dmu0,mueiefi;
double dmueiefidmu,dmueiefidmu0,coef,dcoefdmu0,dcoefdmu,sh_bot,dsh_botdmu,dsh_botdmu0;
double sh,dshdmu,dshdmu0,mu0eiefi,dmu0eiefidmu,dmu0eiefidmu0;
double mu0ei00,dmu0ei00dmu0,mue0epi,dmue0epidmu;
if(i<e)
{
  
    top=pow(sin(0.5*fi),2)*E2i;
    dtopdmu=2*sin(0.5*fi)*cos(0.5*fi)*0.5*dfidmu*E2i;
    dtopdmu0=2*sin(0.5*fi)*cos(0.5*fi)*0.5*dfidmu0*E2i+pow(sin(0.5*fi),2)*dE2idmu0;
    bot=2-E1e-(fi/PI)*E1i;
    dbotdmu=-dE1edmu-(dfidmu/PI)*E1i;
    dbotdmu0=-(dfidmu0/PI)*E1i-(fi/PI)*dE1idmu0;
    
    
    mu0ei0pi=mu0+sin(i)*tth*E2i/(2-E1i);
    dmu0ei0pidmu0=1+((cos(i)*didmu0*tth*E2i+sin(i)*tth*dE2idmu0)*(2-E1i)-sin(i)*tth*E2i*(-dE1idmu0))/pow(2-E1i,2);
    
    
    mue0e0=mu+sin(e)*tth*E2e/(2-E1e);
    dmue0e0dmu=1+((cos(e)*dedmu*tth*E2e+sin(e)*tth*dE2edmu)*(2-E1e)-sin(e)*tth*E2e*(-dE1edmu))/pow(2-E1e,2);
    dmue0e0dmu0=0;
    mueiefi=mu+sin(e)*tth*(E2e-top)/bot;
    
    
    dmueiefidmu=1+((cos(e)*dedmu*tth*(E2e-top)+sin(e)*tth*(dE2edmu-dtopdmu))*bot-sin(e)*tth*(E2e-top)*dbotdmu)/pow(bot,2);
    dmueiefidmu0=((sin(e)*tth*(-dtopdmu0))*bot-sin(e)*tth*(E2e-top)*dbotdmu0)/pow(bot,2);
    coef=mu0/mu0ei0pi;
    
    dcoefdmu0=(mu0ei0pi-mu0*dmu0ei0pidmu0)/pow(mu0ei0pi,2);
    dcoefdmu=0;
    sh_bot=mue0e0*(1-f*(1-coef));
    
    dsh_botdmu=dmue0e0dmu*(1-f+f*coef)+mue0e0*(-dfdmu+dfdmu*coef);
    dsh_botdmu0=dmue0e0dmu0*(1-f+f*coef)+mue0e0*(-dfdmu0+dfdmu0*coef+f*dcoefdmu0);
    
    
    sh=mueiefi*coef/sh_bot;
    
    dshdmu=(dmueiefidmu*coef*sh_bot-mueiefi*coef*dsh_botdmu)/pow(sh_bot,2);
    dshdmu0=((dmueiefidmu0*coef+mueiefi*dcoefdmu0)*sh_bot-mueiefi*coef*dsh_botdmu0)/pow(sh_bot,2);
    
    mu0eiefi=mu0+sin(i)*tth*(cos(fi)*E2e+top)/bot;
    
    dmu0eiefidmu=(sin(i)*tth*(-sin(fi)*dfidmu*E2e+cos(fi)*dE2edmu+dtopdmu)*bot-dbotdmu*(sin(i)*tth*(cos(fi)*E2e+top)))/pow(bot,2);
    dmu0eiefidmu0=1+((cos(i)*didmu0*tth*(cos(fi)*E2e+top)+sin(i)*tth*(-sin(fi)*dfidmu0*E2e+dtopdmu0))*bot-dbotdmu0*(sin(i)*tth*(cos(fi)*E2e+top)))/pow(bot,2);
}
  
  
else
{
    top=pow(sin(0.5*fi),2)*E2e;
    dtopdmu=2*sin(0.5*fi)*cos(0.5*fi)*0.5*dfidmu*E2e+pow(sin(0.5*fi),2)*dE2edmu;
    dtopdmu0=2*sin(0.5*fi)*cos(0.5*fi)*0.5*dfidmu0*E2e;
    
    bot=2-E1i-(fi/PI)*E1e;
    dbotdmu=-(dfidmu/PI)*E1e-(fi/PI)*dE1edmu;
    dbotdmu0=-dE1idmu0-(dfidmu0/PI)*E1e;
    
   
    mu0ei00=mu0+sin(i)*tth*E2i/(2-E1i);
    
    dmu0ei00dmu0=1+((cos(i)*didmu0*tth*E2i+sin(i)*tth*dE2idmu0)*(2-E1i)-sin(i)*tth*E2i*(-dE1idmu0))/pow(2-E1i,2);
    
    mue0epi=mu+sin(e)*tth*E2e/(2-E1e);
    dmue0epidmu=1+((cos(e)*dedmu*tth*E2e+sin(e)*tth*dE2edmu)*(2-E1e)-sin(e)*tth*E2e*(-dE1edmu))/pow(2-E1e,2);
    
    
    mueiefi=mu+sin(e)*tth*(cos(fi)*E2i+top)/bot;
    
    dmueiefidmu=1+((cos(e)*dedmu*tth*(cos(fi)*E2i+top)+sin(e)*tth*(-sin(fi)*dfidmu*E2i+dtopdmu))*bot-sin(e)*tth*(cos(fi)*E2i+top)*dbotdmu)/pow(bot,2);
    dmueiefidmu0=(sin(e)*tth*(-sin(fi)*dfidmu0*E2i+cos(fi)*dE2idmu0+dtopdmu0)*bot-sin(e)*tth*(cos(fi)*E2i+top)*dbotdmu0)/pow(bot,2);
    
    coef=mu0/mu0ei00;
    dcoefdmu0=(mu0ei00-mu0*dmu0ei00dmu0)/pow(mu0ei00,2);
    
    
    
    sh_bot=mue0epi*(1-f*(1-coef));
    dsh_botdmu=dmue0epidmu*(1-f*(1-coef))+mue0epi*(-dfdmu+dfdmu*coef);
    dsh_botdmu0=mue0epi*(-dfdmu0+dfdmu0*coef+f*dcoefdmu0);
    sh=mueiefi*coef/sh_bot;
    
    dshdmu=(dmueiefidmu*coef*sh_bot-mueiefi*coef*dsh_botdmu)/pow(sh_bot,2);
    dshdmu0=((dmueiefidmu0*coef+mueiefi*dcoefdmu0)*sh_bot-mueiefi*coef*dsh_botdmu0)/pow(sh_bot,2);
    mu0eiefi=mu0+sin(i)*tth*(E2i-top)/bot;
    
    dmu0eiefidmu=(sin(i)*tth*(-dtopdmu)*bot-sin(i)*tth*(E2i-top)*dbotdmu)/pow(bot,2);
    dmu0eiefidmu0=1+((cos(i)*didmu0*tth*(E2i-top)+sin(i)*tth*(dE2idmu0-dtopdmu0))*bot-sin(i)*tth*(E2i-top)*dbotdmu0)/pow(bot,2);

  
}


*rsh=sh;

*rmueiefi=mueiefi;
*rmu0eiefi=mu0eiefi;
*rdshdmu=dshdmu;
*rdshdmu0=dshdmu0;
*rdmueiefidmu=dmueiefidmu;
*rdmueiefidmu0=dmueiefidmu0;
*rdmu0eiefidmu=dmu0eiefidmu;
*rdmu0eiefidmu0=dmu0eiefidmu0;

  
}
void dhapke(double *p,double mu,double mu0,double alpha,double *rfh,double *rdfhdmu,double *rdfhdmu0)
{
  double alb,h,S0,g,ta,ca,B0,bdnom,B,fhgdnom,fhg,sqalb,dnom,dnomdmu,dnom0,dnom0dmu0;
  double chmu,dchmudmu,chmu0,dchmu0dmu0,fm,dfmdmu,dfmdmu0,fh,dfhdmu,dfhdmu0;
  //function [fh,dfhdmu,dfhdmu0]=dhapke(p,mu,mu0,alpha)
  alb=p[0];
  h=p[1];
  S0=p[2];
  g=p[3];
  ta=tan(0.5*alpha);
  ca=cos(alpha);
  B0=S0*pow(1+g,2)/(alb*(1-g));
  bdnom=1+ta/h;
  B=B0/bdnom;
  fhgdnom=1+pow(g,2)+2*g*ca;
  fhg=(1-pow(g,2))/(pow(fhgdnom,1.5)); 
  //Particle Phase function
  sqalb=sqrt(1-alb);
  dnom=1+2*mu*sqalb;
  dnomdmu=2*sqalb;
  dnom0=1+2*mu0*sqalb;
  dnom0dmu0=2*sqalb;
  //Multiple scattering functions
  chmu=(1+2*mu)/dnom;
  dchmudmu=(2*dnom-(1+2*mu)*dnomdmu)/pow(dnom,2);
  chmu0=(1+2*mu0)/dnom0; 
  //H(alb,mu0)
  dchmu0dmu0=(2*dnom0-(1+2*mu0)*dnom0dmu0)/pow(dnom0,2);
  fm=chmu*chmu0-1;
  dfmdmu=dchmudmu*chmu0;
  dfmdmu0=dchmu0dmu0*chmu;
  fh=(1+B)*fhg+fm;
  dfhdmu=dfmdmu;
  dfhdmu0=dfmdmu0;
  *rfh=fh;
  *rdfhdmu=dfhdmu;
  *rdfhdmu0=dfhdmu0;
}
