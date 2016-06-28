
#include"utils.h"

#include<lapacke.h>
#include<cblas.h>
#include<stdlib.h>
#include<stdio.h>
void matrix_transpose(double *A,int m,int n)
{
    /*Matrix transpose, overwrites A*/
    double *M=calloc(m*n,sizeof(double));
    for(int i=0;i<m;i++)
        for(int j=0;j<n;j++)
        {
            M[j*m+i]=A[i*n+j];
        }
    memcpy(A,M,m*n*sizeof(double));
    free(M);
}
void flip_dim(double *A,int m,int n,int dim)
{
    /*flips columns (dim==1) or rows (dim==2) of matrix A*/
    double temp;
    if(dim==1)
    {
        for(int i=0;i<n;i++)
            for(int j=0;j<0.5*m;j++)
            {
                temp=A[j*n+i];
                A[j*n+i]=A[(m-1-j)*n+i];
                A[(m-1-j)*n+i]=temp;
            }
    }
    else if(dim==2)
    {
        for(int i=0;i<m;i++)
            for(int j=0;j<0.5*n;j++)
            {
                temp=A[i*n+j];
                A[i*n+j]=A[i*n+n-1-j];
                A[i*n+n-1-j]=temp;
            }
    }
    else
    {
        fprintf(stderr,"In flip_dim, parameter dim must be either 1 or 2\n");
        exit(1);
    }
    
}
void real_matrix_multiplyT_ele(double *A,int *B,int k,int m,double *C)
//A is kxm matrix, B is mxk matrix
//Calculates A.*B' (elementwise)
{
  for(int i=0;i<k;i++)
    for(int j=0;j<m;j++)
      C[i*m+j]=A[i*m+j]*B[j*k+i];
}
void real_matrix_multiplyT(double *A,double *B,int m,int k,int n,double *C)
{
/*
 * A is mxk matrix, B is nxk matrix
 * Calculates A*B'
 */
   
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,m,n,k,1.0,A,k,B,k,0.0,C,n);
}
void matrix_prod(double  *A,int m,int k,double  *B,int n,double  *C)
{
    /*Calculate product C=A*B for  matrices
     * A is mxk matrix, B is kxn matrix
     * and the result C is mxn matrix
     */
   
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,1.0,A,k,B,n,0.0,C,n);
}
void solve_matrix_eqS(double *A,int m,double *B,double *X)
{
     /*A is symmetric mxm matrix, B is mx1 matrix, solution will be placed to X
     * Solve AX=B*/
   int info;
   int ipiv[m];
    memcpy(X,B,m*sizeof(double));
    
info = LAPACKE_dsysv( LAPACK_ROW_MAJOR, 'L',m, 1, A,m, ipiv,
                        X, 1 );
        /* Check for the exact singularity */
        if( info > 0 ) {
                fprintf(stderr, "The element of the diagonal factor " );
                fprintf(stderr, "D(%i,%i) is zero, so that D is singular;\n", info, info );
                fprintf(stderr, "the solution could not be computed.\n" );
                exit( 1 );
        }
}
void solve_matrix_eq(double *A,int m,double *B,double *X)
{
    /*A is mxm matrix, B is mx1 matrix, solution will be placed to X
     * Solve AX=B*/
   
    /*Contents of A are destroyed*/
    int info;
     int ipiv[m];
     memcpy(X,B,m*sizeof(double));
    
     
     info=LAPACKE_dgesv(LAPACK_ROW_MAJOR,m,1,A,m,ipiv,X,1);
     if (info != 0) fprintf(stderr, "failure in solve_matrix_eq with error %d\n", info);
}
void matrix_minus(double *A,int m,int n,double *B)
{
    /*Calculate A=A-B*/
    for(int j=0;j<m*n;j++)
        A[j]=A[j]-B[j];
}
void matrix_minus2(double *A,int m,int n,double *B,double **C)
{
    /*Calculate C=A-B*/
    *C=malloc(m*n*sizeof(double));
    for(int j=0;j<m*n;j++)
        (*C)[j]=A[j]-B[j];
}
void matrix_plus(double *A,int m,int n,double *B)
{
    /*Calculate A=A+B*/
    for(int j=0;j<m*n;j++)
        A[j]=A[j]+B[j];
}
void matrix_plus2(double *A,int m,int n,double *B,double *C)
{
    /*Calculate C=A+B*/
    
    for(int j=0;j<m*n;j++)
        C[j]=A[j]+B[j];
}
void matrix_transprod(double *A,int m,int n,double *B)
{
    /*Calculate A^T*A for mxn matrix. Assume memory is allocated for B*/
    cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,n,n,m,1.0,A,n,A,n,0.0,B,n);
}
void matrix_adddiag(double *A,double *B,int n,double lambda)
{
    memcpy(B,A,sizeof(double)*n*n);
   
    /*Calculate B=A+lambda*diag(A)*/
   for(int j=0;j<n;j++)
    B[j+j*n]+=lambda*A[j+j*n];
}
    

void matrix_prodplusdiag(double *J,int m,int n,double lambda,double **JTJ)
{
    /*Calculate res=J^T*J+lambda*diag(J^TJ), where J is mxn matrix */
    *JTJ=malloc(n*n*sizeof(double));
    cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,n,n,m,1.0,J,n,J,n,0.0,*JTJ,n);
    for(int j=0;j<n;j++)
        (*JTJ)[j+j*n]+=lambda*(*JTJ)[j+j*n];
   // print_matrix(JTJ,n,n);
}
void matrix_vectorprod(double *A,int m,int n,double *V,double *B,int trans)
{
    /*Calculate B=A*V (if trans=0) or B=A^T*V 
     * ASSUME MEMORY FOR B IS ALLOCATED BEFOREHAND */
    if(trans)
        cblas_dgemv(CblasRowMajor,CblasTrans,m,n,1.0,A,n,V,1,0.0,B,1);
    else
        cblas_dgemv(CblasRowMajor,CblasNoTrans,m,n,1.0,A,n,V,1,0.0,B,1);
}
