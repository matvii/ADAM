#include "prepare.h"
int is_in_triangle(double *point,double *dir,double *vert1,double *vert2,double *vert3)
{
  double s1[3];
  double s2[3];
  double b[3];
  double A[3][3];
  double det,gamma,beta;
  
  for(int i=0;i<3;i++) {
    s1[i]=vert1[i]-vert2[i];
    s2[i]=vert1[i]-vert3[i];
    b[i]=vert1[i]-point[i];
  }
  det=s1[0]*s2[1]*dir[2]-s1[0]*dir[1]*s2[2]+s2[0]*dir[1]*s1[2]-s2[0]*s1[1]*dir[2]+dir[0]*s1[1]*s2[2]-dir[0]*s2[1]*s1[2];
  if(fabs(det)<0.000001)
    return 0;
  gamma=(dir[2]*(s1[0]*b[1]-b[0]*s1[1])+dir[1]*(b[0]*s1[2]-s1[0]*b[2])+dir[0]*(s1[1]*b[2]-b[1]*s1[2]))/det;
  if(gamma<0||gamma>1)
    return 0;
  
  beta=(b[0]*(s2[1]*dir[2]-dir[1]*s2[2])+b[1]*(dir[0]*s2[2]-s2[0]*dir[2])+b[2]*(s2[0]*dir[1]-s2[1]*dir[0]))/det;
  if(beta<0||beta>1-gamma)
    return 0;
  return 1;
}
  