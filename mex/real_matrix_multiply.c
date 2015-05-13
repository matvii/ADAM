#include"prepare.h"
void mult_vector(double A[3][3],double *y,double *x)
{
  x[0]=A[0][0]*y[0]+A[0][1]*y[1]+A[0][2]*y[2];
  x[1]=A[1][0]*y[0]+A[1][1]*y[1]+A[1][2]*y[2];
  x[2]=A[2][0]*y[0]+A[2][1]*y[1]+A[2][2]*y[2];
}
void real_matrix_multiply(double *A,double *B,int k,int m,int n,double *C)
{
  double *v;
  double w;
  v=(double*)malloc(m*sizeof(double));
  for(int i=0;i<n;i++)
  {
    //initialize temp vector v
    for(int j=0;j<m;j++)
      v[j]=B[j*n+i];
     
    //calculate matrix element
    
    for(int j=0;j<k;j++)
    {
      w=0;
      for(int r=0;r<m;r++)
	w+=A[r+j*m]*v[r];
      C[j*n+i]=w;
    }
  }
}
void real_matrix_multiplyT(double *A,double *B,int k,int m,int n,double *C)
//A is kxm matrix, B is nxm matrix
//Calculates A*B'
{
  double *v;
  double w;
  for(int i=0;i<k;i++)
    for(int j=0;j<n;j++)
      for(int r=0;r<m;r++)
	C[i*n+j]+=A[i*m+r]*B[j*m+r];
}
void real_matrix_multiplyT_ele(double *A,int *B,int k,int m,double *C)
//A is kxm matrix, B is mxk matrix
//Calculates A.*B' (elementwise)
{
  for(int i=0;i<k;i++)
    for(int j=0;j<m;j++)
      C[i*m+j]=A[i*m+j]*B[j*k+i];
}
void complex_matrix_multiply(double complex *A,double complex *B,int k,int m,int n,double complex *C)
{
  double complex *v;
  double complex w;
  v=(double complex*)malloc(m*sizeof(double complex));
  for(int i=0;i<n;i++)
  {
    //initialize temp vector v
    for(int j=0;j<m;j++)
      v[j]=B[j*n+i];
     
    //calculate matrix element
    
    for(int j=0;j<k;j++)
    {
      w=0;
      for(int r=0;r<m;r++)
	w+=A[r+j*m]*v[r];
      C[j*n+i]=w;
    }
  }
}
void Calculate_Normals(int *tlist,double *vlist,int nfac,int nvert,double *normal)
{
  /*INPUT:
   * tlist 
   * vlist
   * NOTE: MATLAB uses columns first ordering, C uses rows first. We assume rows first because it is convenient
   * nfac number of facets
   * nvert number of vertices
   */
  /*OUTPUT:
   * Normal List of facet normals, rows first ordering since these will not leave C program
   
   */
 
 
  int j1,j2,j3;
  int vindex;
  double *v1,*v2,*v3;
  double*w;
  double side1[3],side2[3];
 
  double norm;
  double cnormal[3];
  
  
 
  //For each facet
  for(int j=0;j<nfac;j++)
  {
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
      side1[i]=v2[i]-v1[i];
      side2[i]=v3[i]-v1[i];
      
    }
    cross(side1,side2,cnormal);
    norm=NORM(cnormal);
    
    //Store centroids and normals
    for(int i=0;i<3;i++)
    {
      cnormal[i]=cnormal[i]/norm;
      normal[3*j+i]=cnormal[i];
    }
  }      
  
} 
 void Calculate_Normal_Derivative(double *w1,double *w2,double *w3,double *dndx1,double *dndx2,double *dndx3,double *dndy1,double *dndy2,double *dndy3,double *dndz1,double *dndz2,double *dndz3)
{
  //Calculate derivatives of normal of triangle (w1 w2 w3) wrt triangle vertex coordinates
  double x1,x2,x3,y1,y2,y3,z1,z2,z3;
  double v1[3],v2[3],cv[3];
  double normc;
  double dcx1[3],dcx2[3],dcx3[3],dcy1[3],dcy2[3],dcy3[3],dcz1[3],dcz2[3],dcz3[3];
  double dareax[3],dareay[3],dareaz[3];
  for(int i=0;i<3;i++)
  {
    v1[i]=w2[i]-w1[i];
    v2[i]=w3[i]-w1[i];
  }
  cross(v1,v2,cv);
  normc=NORM(cv);
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
void Convert_to_Matrix(double *xr,double *yr,double *zr,double R[3][3])
{
  R[0][0]=xr[0];
  R[0][1]=xr[1];
  R[0][2]=xr[2];
  R[1][0]=yr[0];
  R[1][1]=yr[1];
  R[1][2]=yr[2];
  R[2][0]=zr[0];
  R[2][1]=zr[1];
  R[2][2]=zr[2];
}
void Calculate_Area_and_Normal_Derivative(double *w1,double *w2,double *w3,double *n,double *dndx1,double *dndx2,double *dndx3,double *dndy1,double *dndy2,double *dndy3,double *dndz1,double *dndz2,double *dndz3,double *area,double *dAdx,double *dAdy,double *dAdz)
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
  area[0]=0.5*normc;
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
void mult_mat(double A[3][3],double B[3][3],double C[3][3])
{
  for(int i=0;i<3;i++)
  {
    for(int j=0;j<3;j++)
    {
  C[i][j]=A[i][0]*B[0][j]+A[i][1]*B[1][j]+A[i][2]*B[2][j];
    }
  }
} 
void transpose(double M[3][3],double Mt[3][3])
{

for(int i=0;i<3;i++)
  for(int j=0;j<3;j++)
    Mt[i][j]=M[j][i];
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