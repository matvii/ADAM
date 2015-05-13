#include"prepare.h"

void multmatrix(double A[3][3],double B[3][3],double C[3][3])
{
  for(int i=0;i<3;i++)
  {
    for(int j=0;j<3;j++)
    {
  C[i][j]=A[i][0]*B[0][j]+A[i][1]*B[1][j]+A[i][2]*B[2][j];
    }
  }
}


void rotate(double beta,double lambda,double omega,double omega0,double t,double M[3][3],double dMb[3][3],double dMl[3][3],double dMo[3][3])
{
  double f,cf,sf,cb,sb,cl,sl;
  double fmat[3][3],dfm[3][3],blmat[3][3],Dblmat_bet[3][3],Dblmat_lam[3][3];
  f=omega*t+omega0;
  cf=cos(f);
  sf=sin(f);
  cb=cos(beta);
  sb=sin(beta);
  cl=cos(lambda);
  sl=sin(lambda);
  fmat[0][0]=cf;
  fmat[0][1]=sf;
  fmat[0][2]=0;
  fmat[1][0]=-sf;
  fmat[1][1]=cf;
  fmat[1][2]=0;
  fmat[2][0]=0;
  fmat[2][1]=0;
  fmat[2][2]=1;
  
  dfm[0][0]=-t*sf;
  dfm[0][1]=t*cf;
  dfm[0][2]=0;
  dfm[1][0]=-t*cf;
  dfm[1][1]=-t*sf;
  dfm[1][2]=0;
  dfm[2][0]=0;
  dfm[2][1]=0;
  dfm[2][2]=0;
  
  blmat[0][0]=cb*cl;
  blmat[0][1]=cb*sl;
  blmat[0][2]=-sb;
  blmat[1][0]=-sl;
  blmat[1][1]=cl;
  blmat[1][2]=0;
  blmat[2][0]=sb*cl;
  blmat[2][1]=sb*sl;
  blmat[2][2]=cb;
  
  Dblmat_bet[0][0]=-sb*cl;
  Dblmat_bet[0][1]=-sb*sl;
  Dblmat_bet[0][2]=-cb;
  Dblmat_bet[1][0]=0;
  Dblmat_bet[1][1]=0;
  Dblmat_bet[1][2]=0;
  Dblmat_bet[2][0]=cb*cl;
  Dblmat_bet[2][1]=cb*sl;
  Dblmat_bet[2][2]=-sb;
  
  Dblmat_lam[0][0]=-cb*sl;
  Dblmat_lam[0][1]=cb*cl;
  Dblmat_lam[0][2]=0;
  Dblmat_lam[1][0]=-cl;
  Dblmat_lam[1][1]=-sl;
  Dblmat_lam[1][2]=0;
  Dblmat_lam[2][0]=-sb*sl;
  Dblmat_lam[2][1]=sb*cl;
  Dblmat_lam[2][2]=0;
  
  multmatrix(fmat,blmat,M);
  multmatrix(fmat,Dblmat_bet,dMb);
  multmatrix(fmat,Dblmat_lam,dMl);
  multmatrix(dfm,blmat,dMo);
}
  
  