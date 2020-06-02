#ifndef UTILS
#define UTILS
#define DOT(a,b) (a[0]*b[0]+a[1]*b[1]+a[2]*b[2])
#define NORM(a) sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2])
#define DET(a) (a[0]*a[4]*a[8]+a[1]*a[5]*a[6]+a[2]*a[3]*a[7]-a[2]*a[4]*a[6]-a[1]*a[3]*a[8]-a[0]*a[5]*a[7])
#define EP 1E-8
#define PI  3.141592653589793
#define LTS 299792.458
#include<complex.h>
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<omp.h>
#include<string.h>
#include"structs.h"
#ifdef USE_MKL
#include"mkl_lapacke.h"
#include "mkl_types.h"
#include "mkl_cblas.h"
#else
#include<cblas.h>
#include<lapacke.h>
#endif
#include"wcstools-3.9.2/libwcs/fitsfile.h"
void alb_smooth(double *ealb,int nAlbedo,double *albreg,double *dalbreg);
double Albedo_Term(int *tlist,double *vlist,int nfac,int nvert, double *Alimit,double *Alb,int index,double* dAlbv);
void albedo_smooth(int *tlist,double *vlist,int nfac,int nvert,double *ealb,double *Alim,double *res,double *drda);
void inertia(int *tlist,double* vlist,int nfac,int nvert,double *D,int m,int n,double *result,double *angle,double *dres);
void fit_vertex_model_to_LC_AO(LCstruct *LC,AOstruct *AO,OCstruct *OC,RDstruct *RD,CNTRstruct *CR);
void Generate_Normal_Deriv_Matrix_Pad(int *tlist,double *vlist,int nfac,int nvert,int *Nfacets,int *Facetlist,double *D,int padding);
void Scale_Matrix_with_Vector(double *vec,double *M,int m,int n,double *N);
void Calculate_Facet_Normals(int *tlist,double *vlist,int nfac,int nvert,int *Nfacets,int *Facetlist,double *Fnormals);
void Generate_Normal_Deriv_Matrix(int *tlist,double *vlist,int nfac,int nvert,int *Nfacets,int *Facetlist,double *D);
void Find_Facets(int *tlist,double *vlist,int nfac,int nvert,int **Nfacets,int **Facetlist);
int read_vector_fileI_alloc(char *filename,int **buffer,int bufsize);
void read_obj_file(char *file, int **tlist,double **vlist,int *nfac,int *nvert);
int read_state_fileI(char *filename,char *text,int *buffer,int n);
int read_state_file(char *filename,char *text,double *buffer,int n);
void write_obj_file(char *file,int *tlist,double *vlist,int nfac,int nvert);
void AdjFacet(int *tlist,double *vlist,int nfac,int nvert,int *A);
void localsmooth(int *tlist,double *vlist,int nfac,int nvert,double *ealb,double *Alim,double *res,double *drda);
int parse_vector(char *string,double *vec,int maxlength);
int parse_vectorI(char *string,int *vec,int maxlength);
void vector_regularization(double *V,int n,double *sV,double *dV);
void Calculate_RDs(int *tlist,double *vlist,int nfac,int nvert,double *angles,RDstruct  *RDs,double *offset,double *D,int dm,int dn,double *Weight,double *scale,double rexp,double *FT,double *FTdv,double *FTdoff,double *FTdsc,double *FTdxp,int deriv);
int read_vector_fileI(char *filename,int *buffer,int bufsize);
void mask_matrix(int m,int *mask,double **D,int *n);
int find_index(double *vect,int n,double x);
void readfits_rd(char* filename,double **buffer,int cx,int cy,int cxdim,int cydim,double *date,int *xsize,int *ysize,double *cdelt1,double *cdelt2);
RDstruct * process_rd_images(char **filenames,int nao,int *x0,int *y0,int *nx,int *ny, int *cx,int *cy,int *cxdim,int *cydim,double *dx,double *dy,double *dates,double *RadarFreq,double min_tim,double *E,double *TIME,int nephm,int* LowFreq);
void calc_image_fft_unnormed(double *M,int m,int n,double dx,double dy,double *zMr,double *zMi,double *Fx,double *Fy);
double sinc(double x);
int read_vector_file(char *filename,double *M,int n);
int calc_rot_frame(char *fitsfile,double *E,double *up);
void Fit_Occ(int *tlist,double *vlist,int nfac,int nvert,double *angles,double *up,double *E,double *V,double *TIME,double *offset,double *chords,int *type,int nchords,double *W,double *Chordoffset,double *dist,double *dx,double *dy,double *dz,double *dangles,double *dtox,double *dtoy,double *dChordoffset);
void find_chord(int *tlist,double *vlist2,int nfac,int nvert,double *offset,double *a,double *b,int *Adj,int *cledge,int *faedge,double *clpoint,double *fapoint,int *inters,double *dclpx,double *dclpy,double *dfapx,double *dfapy);
void calculate_OCs(int *tlist,double *vlist,int nfac,int nvert,double *angles,OCstruct* OC,double *offset,double *W,double *D,int dm,int dn,double *Chordoffset,double* OCdist,double* dOdv,double *dOdoff,double *dchordoffset);
void Calculate_Temp_deriv(int *tlist, double *vlist,int nfac,int nvert,double* angles,double* E0,double t0,double Gamma, double A,double R,int N,double *Tres,double *Tresdx,double *Tresdy,double *Tresdz, double *TresdA);
void Calculate_Temp(int *tlist, double *vlist,int nfac,int nvert,double* angles,double* E0,double t0,double Gamma, double A,double R,int N,double *Tres);
void Calculate_HF_deriv(int *tlist,double *vlist,int nfac,int nvert,double *angles,double *Eo,double *E0o,double *up,double TIME,double dist,double Gamma,double A,double Hdist,int N,double WL,double *freqx,double *freqy,int nfreq,double *offset,double *Fr,double *Fi,double *dFdxr,double *dFdxi,double  *dFdyr,double  *dFdyi,double *dFdzr,double *dFdzi,double *dFdAr,double *dFdAi,double *dFdoffr,double *dFdoffi);
void Calculate_Normal_Derivative(double *w1,double *w2,double *w3,double *dndx1,double *dndx2,double *dndx3,double *dndy1,double *dndy2,double *dndy3,double *dndz1,double *dndz2,double *dndz3);
void Calculate_Normals(int *tlist,double *vlist,int nfac,int nvert,double *normal);
void Calculate_HF(int *tlist,double *vlist,int nfac,int nvert,double *angles,double *Eo,double *E0o,double *up,double TIME,double dist,double Gamma,double A,double Hdist,int N,double WL,double *freqx,double *freqy,int nfreq,double *offset,double *Fr,double *Fi);
void Calculate_Radiance(int *tlist, double *vlist,int nfac,int nvert,double* angles,double *CE,double* CE0,double t0,double Gamma, double A,double R,double L,int N,double *Flux,double *Fldx,double *Fldy,double *Fldz,double *FldA,int deriv);
OCstruct  *read_occ(char* filename,double min_tim,double *up);
void find_chord(int *tlist,double *vlist2,int nfac,int nvert,double *offset,double *a,double *b,int *Adj,int *cledge,int *faedge,double *clpoint,double *fapoint,int *inters,double *dclpx,double *dclpy,double *dfapx,double *dfapy);
void triangulate_sphere(int nrows,double *t,double *f,int *ifp);
void Calculate_AO_deriv(int *tlist,double *vlist,int nfac,int nvert,double *angles,double *Eo,double *E0o,double *up,double TIME,double dist,double *freqx,double *freqy,int nfreq,double *offset,double *Fr,double *Fi,double *dFdxr,double *dFdxi,double *dFdyr,double *dFdyi,double *dFdzr,double *dFdzi,double *dFdAr,double *dFdAi,double *dFdoffr,double *dFdoffi,double *A,double *Alimit,double *dAr,double *dAi);
double Calculate_AO(int *tlist,double *vlist,int nfac,int nvert,double *angles,double *Eo,double *E0o,double *up,double TIME,double dist,double *freqx,double *freqy,int nfreq,double *offset,double *Fr,double *Fi,double *A,double *Alimit);
void Calculate_AOs(int *tlist,double *vlist,int nfac,int nvert,double *angles,AOstruct *AOs,double *offset,double *D,int dm,int dn,double *Weight,double *scale,double *FT,double *FTdv,double* FTdS,double *Albedo,double *Alimit,double *dA,int deriv);
AOstruct * process_ao_images(char **filenames,char **psfnames,int nao,int *x0,int *y0,int *nx,int *ny,int *xp0,int *yp0,int *npx,int *npy, double *dx,double *dy,double *dates,double min_tim,double *E,double *E0,double *up,double *TIME,int nephm,int* LowFreq);
void mul_cols(double *A,int m,int n,double *V);
void mul_rows(double *A,int m,int n,double *V);
void calc_image_fft(double *M,int m,int n,double dx,double dy,double *zMr,double *zMi,double *Fx,double *Fy);
int read_ephm_data(char *filename,double *TIME,double *E,double *E0);
void readfits(char* filename,double **buffer,int x0,int y0,int nx,int ny,double *date,int *xsize,int *ysize);
void write_shape_file(char *str,int *tlist,double *vlist,int nfac,int nvert);
struct LC  *read_lcurve(char* filename,double min_tim);
void add_vector_to_vlist(double *vlist,double *X,double *Y,int nvert);
void free_lc_struct(LCstruct *LC);
void replace_row(double *M,int m,int n,int k,double *N);
void replace_col(double *M,int m,int n,int k,double *N);
void print_matrix(double *M,int m,int n);
void print_matrixI(int *M,int m,int n);
void print_matrixC(complex double * M,int m,int n);
void print_submatrix(double *M,int m,int n,int k1,int l1,int k2,int l2);
void combine_matrices(double** A,double* B,int m,int n,int k);
double * join_matrices(double *A,double *B,int m,int n,int k);
void  set_el(double *A,int m,int n,double a,int i,int j);
void  set_elI(int *A,int m,int n,int a,int i,int j);
void  set_elS(short *A,int m,int n,int a,int i,int j);
double  get_el(double *A,int m,int n,int i,int j);
int  get_elI(int *A,int m,int n,int i,int j);
int  get_elS(short *A,int m,int n,int i,int j);
int ind2vec(int *A,int m,int **V,int i);
double sum_matelR(double *A,int m,int n,int *V,int l,int k);
double sum_matelC(double *A,int m,int n,int *V,int l,int k);
void read_shape(char* filename,int **facets2,double** vertices2,int *nfac2,int *nvert2,int type3);
void set_submatrix(double *A,int m0,int n0,double *B,int m1,int n1,int k,int l);
void Convert_to_Matrix(double *xr,double *yr,double *zr,double R[3][3]);
void cross(double *A,double *B,double *C);
void Calculate_Area_and_Normal_Derivative(double *w1,double *w2,double *w3,double *n,double *dndx1,double *dndx2,double *dndx3,double *dndy1,double *dndy2,double *dndy3,double *dndz1,double *dndz2,double *dndz3,double *area,double *dAdx,double *dAdy,double *dAdz);
void transpose(double M[3][3],double Mt[3][3]);
void transpose2(double M[3][3],double *Mt);
void Calculate_Frame_Matrix(double *E,double *up,double R[3][3]);
void find_neighborhood(int *tlist,double *vlist,int nfac,int nvert,int *E,int *N,int *E2,int *A);
void cross(double *A,double *B,double *C);
void FacetsOverHorizon(int *tlist,double *vlist,int nfac,int nvert,double *normal,double *centroid,int *NumofBlocks,int *IndexofBlocks);
void Calc_FTC(double *u1,double *v1,int nf,double a,double b,double c,double d,double g,double h,double complex *F);
int is_in_triangle(double *point,double *dir,double *vert1,double *vert2,double *vert3);
void FindActualBlockers(int *tlist,double *vlist,int nfac,int nvert,double *E,double *E0,int nE,int *visible);
void Calc_FTC(double *u1,double *v1,int nf,double a,double b,double c,double d,double g,double h,double complex *F);
void Calc_FTC_deriv(double *u1,double *v1,int nf,double a,double b,double c,double d,double g,double h,double complex *F,double complex *Fda,double complex *Fdb,double complex *Fdc,double complex *Fdd,double complex *Fdg,double complex *Fdh);
void rotate(double beta,double lambda,double omega,double omega0,double t,double M[3][3],double dMb[3][3],double dMl[3][3],double dMo[3][3]);
void mult_mat(double A[3][3],double B[3][3],double C[3][3]);
void mult_vector(double A[3][3],double *y,double *x);
double sum_vector(double *vec1,int n);
void dihedral_angle(int* tlist,double* vlist,int nfac,int nvert,double *res,int *EV,double *dresdx,double *dresdy,double *dresdz);
void convex_reg(int* tlist,double* vlist,int nfac,int nvert,double *D,int dm,int dn,double *res,double *drdv);
void zero_array(double *array,int count);
void area_reg(int *tlist,double *vlist,int nfac,int nvert,double *D,int dm,int dn,double *area,double *dArdv);
void dihedral_angle_reg(int *tlist,double *vlist,int nfac,int nvert,double *D,int dm,int dn,double *result,double *drsdv);
void mult_with_cons(double *A,int m,int n,double C);
void Sqrt3_Subdiv(int *tlist,double* vlist,int nfac,int nvert,int **tlistn,double **vlistn,int *nfacn,int *nvertn,double **D,int sdstep);
void butterfly_subdiv(int *tlist,double *vlist,int nfac,int nvert,int **tlistn,double **vlistn,int *nfacn,int *nvertn,double **D,int sdstep);
void calculate_lcs(int *tlist,double *vlist,int nfac,int nvert,double *angles,LCstruct *LC,double *D,int dm,int dn,double *LCout,double *dLCdv,double *Albedo,double *Alimit,double *dAlb,double *params,double *dparams,int deriv);
void write_matrix_file(char * str,double *M,int m,int n);
void write_matrix_fileI(char * str,int *M,int m,int n);
double SH(int l, int m, double theta, double phi);
void generate_sphere(int nrows,int *tlist,double *vlist);
void generate_ellipsoid(int nrows,double a,double b,double c,int *tlist,double *vlist);
int parse_ini(char *filename);
void fit_subdiv_model_to_LC_AO(LCstruct *LC,AOstruct *AO,OCstruct *OC,RDstruct *RD,CNTRstruct *CR);
void fit_oct_model_to_LC_AO(LCstruct *LC,AOstruct *AO,OCstruct *OC,RDstruct *RD,CNTRstruct *CR);
void octantoid_to_trimesh(double *a,int LMAX,int nrows,int *tlist,double *vlist,double *dvda,int padding);
void octantoid_reg(double *a,int LMAX,double *oreg,double *doreg);
void dhapke_bright(double *E,double *E0,double mu,double mu0,double *p,double th,double *rss,double *rdssdmu,double *rdssdmu0);
void calc_cam_angle(double *E,double angle,double *up,double *upr);
double calc_vol(int *tlist,double *vlist,int nfac,int nvert);
double vol_eq_dia(int *tlist,double *vlist,int nfac,int nvert);
void soft_maxdimz(int *tlist,double *vlist,int nfacn,int nvertn,double *D,int nvert,double zmax,double c,double *zm,double *dz);
int read_weight_file(char *filename,double *W,int max_size);
struct CNTR *read_contour(char *filename,double min_tim,int type,int rotate);
void Calculate_Contours(int *tlist,double *vlist,int nfac,int nvert,double *angles,CNTRstruct* CR,double *offset,double *D,int dm,int dn,double* CRdist,double* dOdv,double *dOdoff);
double minv(double *vlist,int nvert,int index);
double maxv(double *vlist,int nvert,int index);
int find_closest(double *p,double *plistx,double *plisty,int n);
int read_values_from_file(char *filename,double **fbuffer);
double center_of_mass(int *tlist,double *vlist,int nfac,int nvert,double *D,int dm,int dn,double *dv,int deriv);
#endif
