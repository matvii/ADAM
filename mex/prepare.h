#ifndef PREPARE_H
#define PREPARE_H
#define char16_t UINT16_T
#include<complex.h>
#include<stdlib.h>
#include<math.h>
#include <stdio.h>
#include <time.h>
#include "mex.h"
#define DOT(a,b) (a[0]*b[0]+a[1]*b[1]+a[2]*b[2])
#define NORM(a) sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2])
#define EP 1E-8
#define PI  3.141592653589793
#define OMEGA0 PI/180*0
#define LTS 299792.458
void Calculate_Area_and_Normal_Derivative(double *w1,double *w2,double *w3,double *n,double *dndx1,double *dndx2,double *dndx3,double *dndy1,double *dndy2,double *dndy3,double *dndz1,double *dndz2,double *dndz3,double *area,double *dAdx,double *dAdy,double *dAdz);
void Calculate_Radiance(int *tlist, double *vlist,int nfac,int nvert,double* angles,double *CE,double* CE0,double t0,double Gamma, double A,double R,double L,int N,double *Flux,double *Fldx,double *Fldy,double *Fldz,double *FldA,int deriv);
void mult_mat(double A[3][3],double B[3][3],double C[3][3]);
void transpose(double M[3][3],double Mt[3][3]);
void Convert_to_Matrix(double *xr,double *yr,double *zr,double R[3][3]);
void Calculate_Frame_Matrix(double *E,double *up,double R[3][3]);
void Calculate_Normals(int *tlist,double *vlist,int nfac,int nvert,double *normal);
void Calculate_Normal_Derivative(double *w1,double *w2,double *w3,double *dndx1,double *dndx2,double *dndx3,double *dndy1,double *dndy2,double *dndy3,double *dndz1,double *dndz2,double *dndz3);
void mult_vector(double A[3][3],double *y,double *x);
void Calculate_Temp(int *tlist, double *vlist,int nfac,int nvert,double* angles,double* E0,double t0,double Gamma, double A,double R,int N,double *Tres);
void Calculate_Temp_deriv(int *tlist, double *vlist,int nfac,int nvert,double* angles,double* E0,double t0,double Gamma, double A,double R,int N,double *Tres,double *Tdx,double *Tdy,double *Tdz,double *TdA);
void real_matrix_multiplyT(double *A,double *B,int k,int m,int n,double *C);
void real_matrix_multiplyT_ele(double *A,int *B,int k,int m,double *C);
void complex_matrix_multiply(double complex *A,double complex *B,int k,int m,int n,double complex *C);
void Calc_FTC_deriv(double *u1,double *v1,int nf,double a,double b,double c,double d,double g,double h,double complex *F,double complex *Fda,double complex *Fdb,double complex *Fdc,double complex *Fdd,double complex *Fdg,double complex *Fdh);
void FacetsOverHorizon(int *tlist,double *vlist,int nfac,int nvert,double *normal,double *centroid,int *NumofBlocks,int *IndexofBlocks);
void Calc_FTC(double *u1,double* v1,int nf,double a,double b,double c,double d,double g,double h,double complex *F);
void Calc_FTC_deriv(double *u1,double *v1,int nf,double a,double b,double c,double d,double g,double h,double complex *F,double complex *Fda,double complex *Fdb,double complex *Fdc,double complex *Fdd,double complex *Fdg,double complex *Fdh);
int is_in_triangle(double *point,double *dir,double *vert1,double *vert2,double *vert3);
void findblockers(int nfac,int nvert,int nE,double **E,double **E0,double **normal,double **centroid,int *nbl,int *ibll, double **vlist,int **tlist,int *visiblel);
void FindActualBlockers(int *tlist,double *vlist,int nfac,int nvert,double *E,double *E0,int nE,int *visible);
void prepare(int numfac,int numvert,double **vlist,int **tlist,double **normal,double  **centroid,int *nbl,int *ibll);
void prepare2(int numfac,int numvert,double **vlist,int **tlist,double **normal,double **centroid,int *nbl,int *ibll);
void findblockers2(int nfac,int nvert,double *E,double *E0,double **normal,double **centroid,int *nbl,int *ibll, double **vlist,int **tlist,int *visiblel);
void calc_lcurve_yorp(int nE,int numfac,int numvert,int LMAX,int **tlist,double **vlist,double **Eo,double **E0o,double *TIME,double *a,double *theta,double *phi,double *bright,double *dbright);
//void rotate(double beta,double lambda,double omega,double t,double M[3][3],double dMb[3][3],double dMl[3][3],double dMo[3][3]);
void rotate(double beta,double lambda,double omega,double omega0,double t,double M[3][3],double dMb[3][3],double dMl[3][3],double dMo[3][3]);
void rotate_yorp(double beta,double lambda,double omega,double nu,double t,double M[3][3],double dMb[3][3],double dMl[3][3],double dMo[3][3],double dMnu[3][3]);
double SH(int l, int m, double theta, double phi);
void cross(double A[3],double B[3],double C[3]);
void generate_points(int npoints,double **Q,double **cQ,double *v1,double *v2,double *v3);
void rotate2(double *W,double blmat[3][3]);
double fs(double x,int scandir);
double dfs(double x,int scandir);
void calc_lcurve_hapke(int nE,int numfac,int numvert,int LMAX,int **tlist,double **vlist,double **Eo,double **E0o,double *TIME,double *a,double *theta,double *phi,double *bright,double *dbright);
void dhapke_bright(double *E,double *E0,double mu,double mu0,double *p,double th,double *rss,double *rdssdmu,double *rdssdmu0);
void calc_lcurve(int nE,int numfac,int numvert,int LMAX,int **tlist,double **vlist,double **Eo,double **E0o,double *TIME,double *a,double *theta,double *phi,double *bright,double *dbright);
void calculate_lcurve(int nE,int numfac,int numvert,int **tlist,double **vlist,double *angles,double **Eo,double **E0o,double *TIME,double *bright,double *dbrightx,double *dbrighty,double *dbrightz,double *dbrightb,double *dbrightl,double *dbrighto);
#endif