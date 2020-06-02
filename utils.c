#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<complex.h>
#include"structs.h"
#include"utils.h"
#include"globals.h"
double Albedo_Term(int *tlist,double *vlist,int nfac,int nvert, double *Alimit,double *Alb,int index,double* dAlbv)
/*Calculate albedo values for the each facet
 * Here albedo value is mean of vertex values
 */
{
    int v1,v2,v3;
    int talb;
    double la,ha;
    double Albedo;
    la=Alimit[0];
    ha=Alimit[1];
    v1=tlist[3*index]-1;
    v2=tlist[3*index+1]-1;
    v3=tlist[3*index+2]-1;
    Albedo=(la+ha)/2+(ha-la)/2*(tanh(Alb[v1])+tanh(Alb[v2])+tanh(Alb[v3]))/3;
    if(dAlbv!=NULL)
    {
        dAlbv[0]=(ha-la)/2*(1-pow(tanh(Alb[v1]),2))/3;
        dAlbv[1]=(ha-la)/2*(1-pow(tanh(Alb[v2]),2))/3;
        dAlbv[2]=(ha-la)/2*(1-pow(tanh(Alb[v3]),2))/3;
    }
    return Albedo;
}
void AdjFacet(int *tlist,double *vlist,int nfac,int nvert,int *A)
{
    /*A is nvertxnvert matrix, facets with common edge (i1,i2) are A(i1,i2), A(i2,i1)*/
    int i1,i2,i3;
    for(int j=0;j<nfac;j++)
    {
        i1=tlist[3*j]-1;
        i2=tlist[3*j+1]-1;
        i3=tlist[3*j+2]-1;
        set_elI(A,nvert,nvert,j+1,i1,i2);
        set_elI(A,nvert,nvert,j+1,i2,i3);
        set_elI(A,nvert,nvert,j+1,i3,i1);
    }
}
void alb_smooth(double *ealb,int nAlbedo,double *albreg,double *dalbreg)
{
    double sgn;
    double alb;
    for(int i=0;i<nAlbedo;i++)
    {
        alb=ealb[i];
        sgn=(alb > 0) ? 1 : ((alb < 0) ? -1 : 0);
        dalbreg[i*nAlbedo+i]=sgn;
        albreg[i]=abs(alb);
    }
}

void localsmooth(int *tlist,double *vlist,int nfac,int nvert,double *ealb,double *Alim,double *res,double *drda)
{
    /*
     * Facet albedo should be close to neighboring facets
     * INPUT:
     * alb - facet  of LOG albedo values
     * Alim - minimum and maximum albedo
     * OUTPUT:
     * res - nfac array
     * drda - nfac x nfac array
     */
    double *alb=calloc(nfac,sizeof(double));
    double albderiv=0;
    for(int i=0;i<nfac;i++)
        alb[i]=(Alim[0]+Alim[1])/2+(Alim[1]-Alim[0])/2*tanh(ealb[i]);
        
    int *Adj=calloc(nvert*nvert,sizeof(int));
    int i1,i2,i3,n1,n2,n3;
    AdjFacet(tlist,vlist,nfac,nvert,Adj);
    zero_array(drda,nfac*nvert);
    for(int j=0;j<nfac;j++)
    {
        i1=tlist[3*j]-1;
        i2=tlist[3*j+1]-1;
        i3=tlist[3*j+2]-1;
        albderiv=(Alim[1]-Alim[0])/2*(1-pow(tanh(ealb[j]),2));
        n1=get_elI(Adj,nvert,nvert,i2,i1)-1;
        n2=get_elI(Adj,nvert,nvert,i3,i2)-1;
        n3=get_elI(Adj,nvert,nvert,i1,i3)-1;
        res[j]=alb[j]-(alb[n1]+alb[n2]+alb[n3])/3.0;
        set_el(drda,nfac,nfac,albderiv,j,j);
        set_el(drda,nfac,nfac,-1.0/3.0*albderiv,j,n1);
        set_el(drda,nfac,nfac,-1.0/3.0*albderiv,j,n2);
        set_el(drda,nfac,nfac,-1.0/3.0*albderiv,j,n3);
    }
    free(alb);
    free(Adj);
}
int find_index(double *vect,int n,double x)
{
    double TOL=1e-2;
    int index=0;
    double dist=0;
    for(int i=0;i<n;i++)
    {
        dist=fabs(vect[i]-x);
        index=i;
      //  printf("dist: %f\n",dist);
        if(dist<TOL)
            return index;
    }
    fprintf(stderr,"NO VALID DATE WITHIN TOLERANCE LIMITS\n");
    {
        exit(1);
            return -1;
    }
}
double sinc(double x)
{
double limit=sqrt(sqrt(2.22e-16));
if(fabs(x*PI)>limit)
    return sin(x*PI)/(PI*x);
else
    return 1-pow(x*PI,2)/6.0+pow(x*PI,4)/120.0;
}
int calc_rot_frame(char *fitsfile,double *E,double *up)
{
    double upr[]={0,0.397748474527011,0.917494496447491};
    double M[3][3];
    double R[3][3];
    double Rz[3][3];
    double M1[3][3];
    double Mt[3][3];
    double Rot[3][3];
    char *head=NULL;
    char *image;
    float CD11,CD12,CD21,CD22;
    double angle1,angle2;
    int NOCD=0;
    double *buf;
    double det;
    int sgn=0;
    int imsize=0;
    int nlog;
  
    int lhead,nbhead;
    head=fitsrhead(fitsfile,&lhead,&nbhead);
    
    char M11[]="CD1_1";
    
    if(!hgetr4(head,"CD1_1",&CD11))
        NOCD=1;
    if(!hgetr4(head,"CD1_2",&CD12))
        NOCD=1;
    if(!hgetr4(head,"CD2_1",&CD21))
        NOCD=1;
    if(!hgetr4(head,"CD2_2",&CD22))
        NOCD=1;
    if(NOCD==1)
    {
                up[0]=upr[0];
                up[1]=upr[1];
                up[2]=upr[2];
                printf("CD cannot be read, defaulting to equ\n");
                return 0;
    }
    det=CD11*CD22-CD12*CD21;
    sgn=(det > 0) ? 1 : ((det < 0) ? -1 : 0);
    if(sgn>0)
        fprintf(stderr,"Warning: right-handed coordinate system in fits");
    angle1=atan2(sgn*CD12,CD22);
    angle2=atan2(-CD21,sgn*CD11);
   // printf("angle: %f %f\n",angle1*180/PI,angle2*180/PI);
    if(fabs(angle1-angle2)>0.01)
    {
        fprintf(stderr,"Rotation angle information in %s is inconsistent, set rotation manually\n",fitsfile);
        exit(-1);
    }
    Calculate_Frame_Matrix(E,upr,M);
    Rz[0][0]=cos(-angle1);
    Rz[0][1]=-sin(-angle1);
    Rz[0][2]=0;
    Rz[1][0]=sin(-angle1);
    Rz[1][1]=cos(-angle1);
    Rz[1][2]=0;
    Rz[2][0]=0;
    Rz[2][1]=0;
    Rz[2][2]=1;
    mult_mat(Rz,M,M1);
    transpose(M,Mt);
    mult_mat(Mt,M1,Rot);
    mult_vector(Rot,upr,up);
    free(head);
    
   return 1;
}
void mul_cols(double *A,int m,int n,double *V)
{
    /*Multiply each column of A with corresponding element of the vector V (=length n)*/
    for(int j=0;j<m*n;j++)
        A[j]*=V[j%n];
}
 void mul_rows(double *A,int m,int n,double *V)
{
    /*Multiply each row of A with corresponding element of the vector V (=length m)*/
    for(int i=0;i<m;i++)
    for(int j=0;j<n;j++)
        A[i*n+j]*=V[i];
}   
void write_shape_file(char *str,int *tlist,double *vlist,int nfac,int nvert)
{
    FILE *fp;
    fp=fopen(str,"w");
    if(fp==NULL)
    {
        fprintf(stderr,"Cannot open file %s for writing the shape\n",str);
        exit(-1);
    }
    fprintf(fp,"%d %d\n",nvert,nfac);
    for(int j=0;j<nvert;j++)
        fprintf(fp,"%.5f %.5f %.5f\n",vlist[3*j],vlist[3*j+1],vlist[3*j+2]);
    for(int j=0;j<nfac;j++)
        fprintf(fp,"%d %d %d\n",tlist[3*j],tlist[3*j+1],tlist[3*j+2]);
    fclose(fp);
}
void free_lc_struct(LCstruct *LC)
{
    for(int j=0;j<LC->nlc;j++)
    {
        free(LC->E[j]);
        free(LC->E0[j]);
        free(LC->lcs[j]);
        free(LC->TIME[j]);
    }
    free(LC->E);
    free(LC->E0);
    free(LC->lcs);
    free(LC->TIME);
    free(LC->nobs);
    free(LC);
}
void add_vector_to_vlist(double *vlist,double *X,double *Y,int nvert)
{
   for(int j=0;j<nvert;j++)
   {
       Y[3*j]=vlist[3*j]+X[j];
       Y[3*j+1]=vlist[3*j+1]+X[j+nvert];
       Y[3*j+2]=vlist[3*j+2]+X[j+2*nvert];
   }
}

void write_matrix_file(char * str,double *M,int m,int n)
{
    FILE *fp;
    fp=fopen(str,"w");
    if(fp==NULL)
    {
        fprintf(stderr,"Cannot open file %s for writing the matrix\n",str);
        exit(-1);
    }
    for(int j=0;j<m;j++)
    {
        for(int k=0;k<n;k++)
        {
            fprintf(fp,"%.9f ",M[j*n+k]);
        }
        if(m>1)
        fprintf(fp,"\n");
    }
    fclose(fp);
}
void write_matrix_fileI(char * str,int *M,int m,int n)
{
    FILE *fp;
    fp=fopen(str,"w");
    if(fp==NULL)
    {
        fprintf(stderr,"Cannot open file %s for writing the matrix\n",str);
        exit(-1);
    }
    for(int j=0;j<m;j++)
    {
        for(int k=0;k<n;k++)
        {
            fprintf(fp,"%d ",M[j*n+k]);
        }
        if(m>1)
        fprintf(fp,"\n");
    }
    fclose(fp);
}
void mult_with_cons(double *A,int m,int n,double C)
{
    for(int j=0;j<m*n;j++)
        A[j]=C*A[j];
}

double sum_vector(double *vec1,int n)
{
  double res=0;
  for(int j=0;j<n;j++)
    res=res+vec1[j];
  return res;
}

  
struct LC  *read_lcurve(char* filename,double min_tim)
{
  struct LC *tlc;
  char *buffer=calloc(2048,sizeof(char));
  tlc=malloc(sizeof(struct LC));
  FILE *fid;
  fid=fopen(filename,"r");
  if(fid==NULL)
  {
      perror("Cannot open lcurve file!");
      exit(-1);
  }
  int cal,nlc,nobs;
  int ncalib=0;
  double *br,cumsum;
  double E[3],E0[3];
  double lE,lE0;
  double time;
  int total=0;
  fgets(buffer,2048,fid);
  
  if(sscanf(buffer,"%d",&nlc)==0) /*Total number of lightcurves*/
  {
    fprintf(stderr,"Error reading lcurve file (total number of lightcurves)!");
      exit(-1);
  }
  
  (*tlc).nlc=nlc;
  (*tlc).nobs=malloc(nlc*sizeof(int));
  (*tlc).lcs=malloc(nlc*sizeof(double*));
  (*tlc).E=malloc(nlc*sizeof(double*));
  (*tlc).E0=malloc(nlc*sizeof(double*));
  (*tlc).TIME=malloc(nlc*sizeof(double*));
  (*tlc).rel=calloc(nlc,sizeof(int));
  
  tlc->calib=0;
  for(int i=0;i<nlc;i++)
  {
      fgets(buffer,2048,fid);
     
      if(sscanf(buffer,"%d %d",&nobs,&cal)!=2)
          {
    fprintf(stderr,"Error reading lcurve file!(nobs: %d,cal: %d), lc: %d\n",nobs,cal,i+1);
    fprintf(stderr,"buffer was: %s\n",buffer);
      exit(-1);
  }
 
          /*Only relative for now, remember to fix this*/
      /*TODO: Also remember to remove comments*/
      (*tlc).nobs[i]=nobs;
      (*tlc).lcs[i]=malloc(nobs*sizeof(double));
      (*tlc).E[i]=malloc(3*nobs*sizeof(double));
      (*tlc).E0[i]=malloc(3*nobs*sizeof(double));
      (*tlc).TIME[i]=malloc(nobs*sizeof(double));
      br=malloc(nobs*sizeof(double));
      cumsum=0;
      if(INI_LC_ARE_RELATIVE==1)
      {
          cal=0;
      }
      if(cal==0)
          (*tlc).rel[i]=1;
      for(int j=0;j<nobs;j++)
      {
          fgets(buffer,2048,fid);
         
         if(sscanf(buffer,"%lf %lf %lf %lf %lf %lf %lf %lf",&time,&(br[j]),E0,E0+1,E0+2,E,E+1,E+2)!=8)
             {
    perror("Error reading lcurve file!");
      exit(-1);
             }
          cumsum+=br[j];
          lE=NORM(E);
          lE0=NORM(E0);
          (*tlc).E0[i][3*j]=E0[0]/lE0;
          (*tlc).E0[i][3*j+1]=E0[1]/lE0;
          (*tlc).E0[i][3*j+2]=E0[2]/lE0;
          
          (*tlc).E[i][3*j]=E[0]/lE;
          (*tlc).E[i][3*j+1]=E[1]/lE;
          (*tlc).E[i][3*j+2]=E[2]/lE;
          
          (*tlc).TIME[i][j]=time-min_tim;
      }
      /*Copy relative brightness values to array*/
      /////////////////////////////////////////////////////////7
      //DO NOT USE CALIBRATED
     // cal=0;
      //////////////////////////////////////////////////////////
      
      
      if(cal==0)
      {
          for(int k=0;k<nobs;k++)
              (*tlc).lcs[i][k]=br[k]*nobs/cumsum;
      }
      else
      {
          for(int k=0;k<nobs;k++)
              (*tlc).lcs[i][k]=br[k];
          tlc->calib=1;
          printf("Lcurve %d is  calibrated\n",i+1);
          ncalib++;
      }
              
          
      free(br);
      total+=nobs;
      
  }
  free(buffer);
  fclose(fid);
  tlc->ntotal=total;
  tlc->ncalib=ncalib;
  return tlc;
}

void set_submatrix(double *A,int m0,int n0,double *B,int m1,int n1,int k,int l)
{
    /*Replace submatrix of A with matrix B. Upper left corner of B is placed at (k,l)
     * Indexing is from 0*/
    if(k+m1-1>=m0)
    {
        printf("Matrix size %d %d, submatrix %d %d, at %d %d\n",m0,n0,m1,n1,k,l);
        puts("Error: index too large in set_submatrix\n");
        exit(1);
    }
    if(l+n1-1>=n0)
    {
        printf("Matrix size %d %d, submatrix %d %d, at %d %d\n",m0,n0,m1,n1,k,l);
        puts("Error: index too large in set_submatrix\n");
        exit(1);
    }
    for(int i=0;i<m1;i++)
        for(int j=0;j<n1;j++)
            A[k*n0+l+j+i*n0]=B[j+i*n1];
}
void replace_row(double *M,int m,int n,int k,double *N)
{
    /*Replace kth row of the matrix mxn M  with N
     * Indexing starts from 0 */
    for(int j=0;j<n;j++)
        M[j+k*n]=N[j];
}
void replace_col(double *M,int m,int n,int k,double *N)
{
    /*Replace kth column of the matrix mxn M with N
     * Indexing start from 0*/
    for(int j=0;j<m;j++)
        M[j*n+k]=N[j];
}
int read_state_file(char *filename,char *text,double *buffer,int n)
{
    /*
     * In a file named filename, find a line starting with text,
     * then read the next line into buffer
     */ 
    if(n==0)
        return 0;
     FILE *fid;
    fid=fopen(filename,"r");
     char delims[]=" \t\r\n\f\v,";
     int count=0;
     int cmp=0;
     char *buff;
     buff=malloc(10000*sizeof(double));
     if(fid==NULL)
    {
        perror("Error opening file in read_state_file");
        exit(-1);
    }
    
    while(fgets(buff,10000*sizeof(double),fid)!=NULL)
    {
       // if(buff[0]=='#' || buff[0]==';')
        //    continue;
        
        
        if(strncmp(text,buff,strlen(text))==0)
        {
            fgets(buff,10000*sizeof(double),fid);
            count=parse_vector(buff,buffer,n);
            break;
            
    }
    }
    return count;
    free(buff);
    fclose(fid);
}
int read_state_fileI(char *filename,char *text,int *buffer,int n)
{
    /*
     * In a file named filename, find a line starting with text,
     * then read the next line into buffer
     */ 
    
    if(n==0)
        return 0;
     FILE *fid;
    fid=fopen(filename,"r");
     char delims[]=" \t\r\n\f\v,";
     int count=0;
     int cmp=0;
     char *buff;
     buff=malloc(10000*sizeof(int));
     if(fid==NULL)
    {
        perror("Error opening file in read_state_file");
        exit(-1);
    }
  
    while(fgets(buff,10000*sizeof(int),fid)!=NULL)
    {
       // if(buff[0]=='#' || buff[0]==';')
        //    continue;
        
        if(strncmp(text,buff,strlen(text))==0)
        {
            
            fgets(buff,10000*sizeof(int),fid);
            
            count=parse_vectorI(buff,buffer,n);
           
            break;
            
    }
    }
    return count;
    free(buff);
    fclose(fid);
}
void print_submatrix(double *M,int m,int n,int k1,int l1,int k2,int l2)
{
    /*Print submatrix of mxn matrix*/
    for(int j=0;j<k2;j++)
    {
        printf("\n");
        for(int i=0;i<l2;i++)
            printf(" %.4f ",M[(j+k1)*n+i+l1]);
    }
    printf("\n");
}
void print_matrix(double *M,int m,int n)
{
    /*Print mxn matrix*/
    for(int j=0;j<m;j++)
    {
        printf("\n");
        for(int i=0;i<n;i++)
            printf(" %.4f ",M[j*n+i]);
    }
    printf("\n");
}
void print_matrixI(int *M,int m,int n)
{
    /*Print mxn matrix*/
    for(int j=0;j<m;j++)
    {
        printf("\n");
        for(int i=0;i<n;i++)
            printf(" %d ",M[j*n+i]);
    }
    printf("\n");
}
void print_matrixC(double complex *M,int m,int n)
{
    /*Print mxn matrix*/
    for(int j=0;j<m;j++)
    {
        printf("\n");
        for(int i=0;i<n;i++)
            printf(" %.2f+%.2fI ",creal(M[j*n+i]),cimag(M[j*n+i]));
    }
    printf("\n");
}
void combine_matrices(double** A,double* B,int m,int n,int k)
{
    /*A is mxn matrix,B is kxn matrix, combine matrices to form
     *to form a matrix with B appended below A.
     * A IS DESTROYED IN THE PROCESS
     * USE ONLY IF A IS ORIGINALLY ALLOCATED WITH MALLOC/CALLOC */
    double *C;
    C=realloc(*A,(m+k)*n*sizeof(double));
    if(C==NULL)
    {
        puts("Error: Cannot reallocate memory!");
        exit(1);
    }
    for(int j=m*n;j<(m+k)*n;j++)
        C[j]=B[j-m*n];
    *A=C;
}
double * join_matrices(double *A,double *B,int m,int n,int k)
{
    /*A is mxn matrix, B is mxk, function returns pointer to mx(n+k) matrix
     * [A B] */
    double *C;
    C=malloc(m*(n+k)*sizeof(double));
    if(C==NULL)
        {
        puts("Error: Cannot allocate memory!");
        exit(1);
    }
    /*First copy A to C*/
    for(int i=0;i<m;i++)
        for(int j=0;j<n;j++)
            C[i*(n+k)+j]=A[i*n+j];
    /*Then copy B to C*/
    for(int i=0;i<m;i++)
        for(int j=0;j<k;j++)
            C[i*(n+k)+j+n]=B[i*k+j];
   return C;
}
void  set_el(double *A,int m,int n,double a,int i,int j)
{
    /*Given mxn matrix A, set element (i,j) to a 
     * indexing starts from 0*/
    if(i>m-1 || j>n-1)
    {
        puts("Error: error in set_element, index too large. Bailing out!");
        exit(1);
    }
    A[i*n+j]=a;
}
double  get_el(double *A,int m,int n,int i,int j)
{
    /*Given mxn matrix A, get element at (i,j) 
     * indexing starts from 0*/
    if(i>m-1 || j>n-1)
    {
        puts("Error: error in get_element, index too large. Bailing out!");
        exit(1);
    }
    return A[i*n+j];
}
void  set_elI(int *A,int m,int n,int a,int i,int j)
{
    /*Given mxn matrix A, set element (i,j) to a 
     * indexing starts from 0
     * This version is for ints*/
    if(i>m-1 || j>n-1)
    {
        puts("Error: error in set_element, index too large. Bailing out!");
        exit(1);
    }
    A[i*n+j]=a;
}
void  set_elS(short *A,int m,int n,int a,int i,int j)
{
    /*Given mxn matrix A, set element (i,j) to a 
     * indexing starts from 0
     * This version is for short ints*/
    if(i>m-1 || j>n-1)
    {
        puts("Error: error in set_element, index too large. Bailing out!");
        exit(1);
    }
    A[i*n+j]=a;
}
int  get_elS(short *A,int m,int n,int i,int j)
{
    /*Given mxn matrix A, get element at (i,j) 
     * indexing starts from 0*/
    if(i>m-1 || j>n-1)
    {
        puts("Error: error in get_element, index too large. Bailing out!");
        exit(1);
    }
    return A[i*n+j];
}
int  get_elI(int *A,int m,int n,int i,int j)
{
    /*Given mxn matrix A, get element at (i,j) 
     * indexing starts from 0*/
    if(i>m-1 || j>n-1)
    {
        puts("Error: error in get_element, index too large. Bailing out!");
        exit(1);
    }
    return A[i*n+j];
}
int ind2vec(int *A,int m,int **V,int i)
{
    /*Assume A is a binary mxm binary matrix 
     * Routine converts ith row to vector. eg
     * [0,1,0,1]->[1,3] and
     * returns the length of V*/
    int *W=malloc(m*sizeof(int));
    int count=0;
    for(int j=0;j<m;j++)
    {
        
        if(get_elI(A,m,m,i,j))
            W[count++]=j;
    }
    *V=W;
    return count;
}
double sum_matelR(double *A,int m,int n,int *V,int l,int k)
{
    /* For mxn matrix A, return sum of those elements in the kth row
     * that are indexed by the vector V of length l*/
    double res=0;
    for(int j=0;j<l;j++)
        res+=get_el(A,m,n,k,V[j]);
    return res;
}
double sum_matelC(double *A,int m,int n,int *V,int l,int k)
{
    /* For mxn matrix A, return sum of those elements in the kth column
     * that are indexed by the vector V of length l*/
    double res=0;
    for(int j=0;j<l;j++)
        res+=get_el(A,m,n,V[j],k);
    return res;
}
void read_shape(char* filename,int **facets2,double** vertices2,int *nfac2,int *nvert2,int type3)
{
    /* Read shape and return vertices and facets. If type3 is nonzero, then file is assumed to contain '3' between
     * every facet */
    int nvert,nfac;
    FILE *fid;
    fid=fopen(filename,"r");
    if(fid==0)
    {
        fprintf(stderr,"Cannot open shape file %s\n",filename);
        exit(-1);
    }
    if(fscanf(fid,"%d %d",&nvert,&nfac)!=2)
    {
        perror("Error while reading shape file");
        exit(-1);
    }
    
    int *facets;
    double *vertices;
    int temp;
    facets=malloc(3*nfac*sizeof(int));
    vertices=malloc(3*nvert*sizeof(double));
    for(int i=0;i<nvert;i++)
        if(fscanf(fid,"%lf %lf %lf",&vertices[i*3],&vertices[i*3+1],&vertices[i*3+2])!=3)
        {
            perror("Error while reading shape file");
            exit(-1);
        }
        if(type3)
            for(int i=0;i<nfac;i++)
            {
                fscanf(fid,"%d",&temp);
                if(fscanf(fid,"%d %d %d",&facets[i*3],&facets[i*3+1],&facets[i*3+2])!=3)
                {
                    perror("Error while reading shape file");
                    exit(-1);
                }
            }
            else
                for(int i=0;i<nfac;i++)
                {
                   
                    if(fscanf(fid,"%d %d %d",&facets[i*3],&facets[i*3+1],&facets[i*3+2])!=3)
                    {
                        perror("Error while reading shape file");
                        exit(-1);
                    }
                    
                }
                    fclose(fid);
                *nvert2=nvert;
            *nfac2=nfac;
            *facets2=facets;
            *vertices2=vertices;
}
void write_obj_file(char *file,int *tlist,double *vlist,int nfac,int nvert)
{
    /*
     * Write shape file in obj (wavefront) format
     */
    FILE *fp;
    fp=fopen(file,"w");
    if(fp==NULL)
    {
        fprintf(stderr,"Cannot open file %s for writing the shape in obj format\n",file);
        exit(-1);
    }
    
    for(int j=0;j<nvert;j++)
        fprintf(fp,"v %.4f %.4f %.4f\n",vlist[3*j],vlist[3*j+1],vlist[3*j+2]);
    for(int j=0;j<nfac;j++)
        fprintf(fp,"f %d %d %d\n",tlist[3*j],tlist[3*j+1],tlist[3*j+2]);
    fclose(fp);
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
void transpose2(double M[3][3],double *Mt)
{

for(int i=0;i<3;i++)
  for(int j=0;j<3;j++)
    Mt[i*3+j]=M[j][i];
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
void find_neighborhood(int *tlist,double *vlist,int nfac,int nvert,int *E,int *N,int *E2,int *A)
{
    /*E(i,j)=1 if there is and edge from vertex i to j, nvertxnvert matrix
     * N(i,j)=1, if vertex i belongs to facet j nvertxnfac matrix
     * E2(i,j)=k if edge(i,j) is in the triangle k, nvert x nvert matrix
     * A(i,j)=1 if facet i and j are adjacent, nfacxnfac matrix
     * ALL INDEXING IS FROM ZERO
     * NB: THIS WORKS ONLY FOR SHAPES WITHOUT BOUNDARY. IF BOUNDARY EXISTS, WE SHOULD TEST IF E2[i2*nvert+i1]>0 BEFORE SETTING A*/
    /*NB2: E2 is a bit problematic, maybe we should index from 1?*/
    int i1,i2,i3;
    
    for(int j=0;j<nfac;j++)
    {
        i1=tlist[3*j]-1;
        i2=tlist[3*j+1]-1;
        i3=tlist[3*j+2]-1;
        E[i1*nvert+i2]=1;
        E[i2*nvert+i1]=1;
        
        E[i2*nvert+i3]=1;
        E[i3*nvert+i2]=1;
        
        E[i1*nvert+i3]=1;
        E[i3*nvert+i1]=1;
        
        N[i1*nfac+j]=1;
        N[i2*nfac+j]=1;
        N[i3*nfac+j]=1;
        set_elI(E2,nvert,nvert,j,i1,i2);
        set_elI(E2,nvert,nvert,j,i2,i3);
        set_elI(E2,nvert,nvert,j,i3,i1);
        
    }
    for(int j=0;j<nfac;j++)
    {
        i1=tlist[3*j]-1;
        i2=tlist[3*j+1]-1;
        i3=tlist[3*j+2]-1;
        A[j*nfac+E2[i2*nvert+i1]]=1;
        A[j*nfac+E2[i3*nvert+i2]]=1;
        A[j*nfac+E2[i1*nvert+i3]]=1;
    }
}

void cross(double *A,double *B,double *C)
{
   /*Calculate the cross product of two vectors A and B, place the result
 * in the vector C
 */ 
  C[0]=A[1]*B[2]-A[2]*B[1];
  C[1]=A[2]*B[0]-A[0]*B[2];
  C[2]=A[0]*B[1]-A[1]*B[0];
}
void mult_vector(double A[3][3],double *y,double *x)
{
  x[0]=A[0][0]*y[0]+A[0][1]*y[1]+A[0][2]*y[2];
  x[1]=A[1][0]*y[0]+A[1][1]*y[1]+A[1][2]*y[2];
  x[2]=A[2][0]*y[0]+A[2][1]*y[1]+A[2][2]*y[2];
}

void zero_array(double *array,int count)
{
    memset(array,0,sizeof(double)*count);
}
void mask_matrix(int m,int *mask,double **D,int *n)
{
    /*
     * Generate mxm identity matrix, remove columns
     * corresponding to positions of ones in the array mask
     * Final matrix is mxn matrix
     * OUTPUT:
     * D and n
     */
    int count=0;
    for(int j=0;j<m;j++)
        if(mask[j]==0)
            count++;
    *n=count;
   
    double *M=calloc(m*count,sizeof(double));
    int i=0;
    for(int j=0;j<m;j++)
    {
        if(mask[j]==0)
        {
            M[i+j*count]=1;
            i++;
        }
    }
    *D=M;
}
int read_vector_file(char *filename,double *buffer,int bufsize)
{
    /*
     * Read <=bufsize vector of double values from a file. Lines starting with # or ;
     *are ignored. Any delimeter( \t\r\n\f\v,) can be used to separate values.
     */ 
    FILE *fid;
    fid=fopen(filename,"r");
     char delims[]=" \t\r\n\f\v,";
     char *token;
     char *buff,*filebuff;
     buff=malloc(100000);
     filebuff=malloc(100000);
     if(fid==NULL)
    {
        perror("Error opening file in read_vector_file");
        exit(-1);
    }
    int count=0;
    while(fgets(buff,100000,fid)!=NULL)
    {
        if(buff[0]=='#' || buff[0]==';')
            continue;
        memcpy(filebuff,buff,100000);
    token=strtok(filebuff,delims);
    while(token!=NULL && count<bufsize)
    {
        buffer[count]=atof(token);
        token=strtok(NULL,delims);
        count++;
    }
    }
    free(buff);
    free(filebuff);
    return count;
}
int read_vector_fileI_alloc(char *filename,int **buffer,int bufsize)
{
    /*
     * Read <=bufsize vector of int values from a file. Lines starting with # or ;
     *are ignored. Any delimeter( \t\r\n\f\v,) can be used to separate values.
     * Buffer is allocated HERE
     */ 
    FILE *fid;
    fid=fopen(filename,"r");
     char delims[]=" \t\r\n\f\v,";
     char *token;
     char *buff,*filebuff;
     buff=malloc(40000);
     filebuff=malloc(40000);
     int *buffer2=calloc(10000,sizeof(int));
     if(fid==NULL)
    {
        perror("Error opening file in read_vector_file");
        return -1;
    }
    int count=0;
    while(fgets(buff,40000,fid)!=NULL)
    {
        if(buff[0]=='#' || buff[0]==';')
            continue;
        memcpy(filebuff,buff,40000);
    token=strtok(filebuff,delims);
    while(token!=NULL && count<bufsize)
    {
        buffer2[count]=atoi(token);
        token=strtok(NULL,delims);
        count++;
    }
    }
    *buffer=calloc(count,sizeof(int));
    memcpy(*buffer,buffer2,sizeof(int)*count);
    free(buffer2);
    free(buff);
    free(filebuff);
    return count;
}
int parse_vector(char *string,double *vec,int maxlength)
{
    char delims[]=" \t\r\n\f\v,";
     char *token;
    token=strtok(string,delims);
    int count=0;
    while(token!=NULL && count<maxlength)
    {
        vec[count]=atof(token);
        token=strtok(NULL,delims);
        count++;
    }
    return count;
}
int parse_vectorI(char *string,int *vec,int maxlength)
{
    char delims[]=" \t\r\n\f\v,";
   
     char *token;
    token=strtok(string,delims);
    int count=0;
   
    while(token!=NULL && count<maxlength)
    {
       
        vec[count]=atoi(token);
        token=strtok(NULL,delims);
       
        count++;
    }
    return count;
}
int read_vector_fileI(char *filename,int *buffer,int bufsize)
{
    /*
     * Read <=bufsize vector of int values from a file. Lines starting with # or ;
     *are ignored. Any delimeter( \t\r\n\f\v,) can be used to separate values.
     */ 
     FILE *fid;
    fid=fopen(filename,"r");
     char delims[]=" \t\r\n\f\v,";
     char *token;
     char *buff,*filebuff;
     buff=malloc(100000);
     filebuff=malloc(100000);
     if(fid==NULL)
    {
        perror("Error opening file in read_vector_file");
        exit(-1);
    }
    int count=0;
    while(fgets(buff,100000,fid)!=NULL)
    {
        if(buff[0]=='#' || buff[0]==';')
            continue;
        memcpy(filebuff,buff,100000);
    token=strtok(filebuff,delims);
    while(token!=NULL && count<bufsize)
    {
        buffer[count]=atoi(token);
        token=strtok(NULL,delims);
        count++;
    }
    }
    free(buff);
    free(filebuff);
    return count;
}
 void vector_regularization(double *V,int n,double *sV,double *dV)
{
    /*
     * Calculate Sum(V[i])^2
     * and its derivatives wrt V[i]
     * OUTPUT:
     * n element array dV, where dV[i] is derivative wrt V[i]
     * */
    *sV=0.0;
    for(int j=0;j<n;j++)
    {
        (*sV)+=pow(V[j],2);
        dV[j]=2*V[j];
    }
}
void calc_cam_angle(double *E,double angle,double *up,double *upr)
{
    /*
     * Calculate the up vector that results from rotation the image plane
     *by angle (degrees). 
     *INPUT:
     *E observer direction as seen from the world frame 
     *up Camera up direction
     *OUTPUT:
     *upr: final camera up direction
     */ 
    double angle1=angle*PI/180.0;
    double M[3][3];
    double Mt[3][3];
    double Rz[3][3];
    double Rot[3][3];
    double M1[3][3];
    Calculate_Frame_Matrix(E,up,M);

    Rz[0][0]=cos(-angle1);
    Rz[0][1]=-sin(-angle1);
    Rz[0][2]=0;
    Rz[1][0]=sin(-angle1);
    Rz[1][1]=cos(-angle1);
    Rz[1][2]=0;
    Rz[2][0]=0;
    Rz[2][1]=0;
    Rz[2][2]=1;
    mult_mat(Rz,M,M1);
    transpose(M,Mt);
    mult_mat(Mt,M1,Rot);
    mult_vector(Rot,up,upr);
}

double calc_vol(int *tlist,double *vlist,int nfac,int nvert)
{
    /*Calculate the volume of triangular shape
     */
    double vol=0;
    int i1,i2,i3;
    double *v1,*v2,*v3;
    double tV[3];
    for(int j=0;j<nfac;j++)
    {
        i1=tlist[3*j]-1;
        i2=tlist[3*j+1]-1;
        i3=tlist[3*j+2]-1;
        v1=vlist+3*i1;
        v2=vlist+3*i2;
        v3=vlist+3*i3;
        cross(v2,v3,tV);
        vol+=1.0/6.0*(v1[0]*tV[0]+v1[1]*tV[1]+v1[2]*tV[2]);
    }
    return vol;
}
double vol_eq_dia(int *tlist,double *vlist,int nfac,int nvert)
{
    double vol;
    vol=calc_vol(tlist,vlist,nfac,nvert);
    return 2*cbrt(3.0/(4.0*PI)*vol);
}
int read_weight_file(char *filename,double *W,int max_size)
{
    /*Read information from  a file s
     * Example:
     * 1 0.5
     * 2 0.3
     * 5 0.2
     * Sets W[0]=0.5, W[1]=0.3 and W[5]=0.2 and leaves everything else unchanged
     * Maximum allowable index is max_size
     */
    char *buffer=calloc(2048,sizeof(char));
 
  FILE *fid;
  fid=fopen(filename,"r");
 
  int index=0;
  double w=0.1;
  int count=0;
  if(fid==NULL)
  {
      perror("Cannot open lcurve weights file!");
      exit(-1);
  }
  while(fgets(buffer,2000,fid)!=NULL)
    {
      
        if(buffer[0]=='#' || buffer[0]==';' ||buffer[0]=='\n' ||buffer[0]=='\0')
            continue;
        if(sscanf(buffer,"%d %lf",&index,&w)!=2)
        {
            fprintf(stderr,"Error parsing a line from the lcurve weight file (was: %s), ignoring\n",buffer);
            continue;
        }
        if(index>max_size)
        {
            fprintf(stderr,"Index %d in the lcurve weight file is larger than the number of lightcurves, ignoring\n",index);
            continue;
        }
       count++;
        W[index-1]=w;
    }
    return count;
}
struct CNTR *read_contour(char *filename,double min_tim,int type,int rotate)
{
    /*
     * Read contours from a file
     * min_tim is substracted from observation times
     * type=0 if contours are given as x-y coordinates
     * type=1 if contours are radius-angle
     * rotate=0: Ignore tilt, assume oriented north up in eq frame
     * rotate=1: Use tilt, assume ecliptic frame
     * ARE TIMES LT-CORRECTED??
     */
    struct CNTR *C;
    C=malloc(sizeof(struct CNTR));
    FILE *fid;
    int ncont;
    fid=fopen(filename,"r");
    double deg2rad=PI/180.0;
    double angle,radius;
     char delims[]=" \t\r\n\f\v,";
     char *token;
     char *buffer,*filebuff;
     double time,tilt,pixscale;
     double scale;
     double lE,lE0;
     double x,y;
     int nobs;
     buffer=malloc(10000);
     double E[3],E0[3];
     int count=0;
     if(fid==NULL)
    {
        perror("Error opening the contour file\n");
        exit(-1);
    }
    fgets(buffer,2048,fid);
     if(sscanf(buffer,"%d",&ncont)==0) /*Total number of contours*/
  {
    fprintf(stderr,"Error reading contour file (total number of contours)!");
      exit(-1);
  }
  (*C).ncont=ncont;
  (*C).nobs=calloc(ncont,sizeof(int));
  (*C).datax=calloc(ncont,sizeof(double*));
  (*C).datay=calloc(ncont,sizeof(double*));
  (*C).TIME=calloc(ncont,sizeof(double));
  (*C).E=calloc(3*ncont,sizeof(double));
  (*C).E0=calloc(3*ncont,sizeof(double));
  (*C).up=calloc(3*ncont,sizeof(double));
  (*C).distance=calloc(ncont,sizeof(double));
  //Loop over contours
  for(int j=0;j<ncont;j++)
  {
    fgets(buffer,2048,fid);
     
      if(sscanf(buffer,"%lf %lf %lf",E0,E0+1,E0+2)!=3)
          {
            fprintf(stderr,"Error reading E0 in the contour file!\n");
            exit(-1);
          }
      fgets(buffer,2048,fid);    
    if(sscanf(buffer,"%lf %lf %lf",E,E+1,E+2)!=3)
          {
            fprintf(stderr,"Error reading E in the contour file!\n");
            exit(-1);
        }
       fgets(buffer,2048,fid);
    
    if(sscanf(buffer,"%lf %lf",&time,&tilt)!=2)
          {
            fprintf(stderr,"Error reading (time,tilt) in the contour file!\n");
            exit(-1);
        }
         fgets(buffer,2048,fid);
        if(sscanf(buffer,"%lf",&pixscale)!=1)
          {
            fprintf(stderr,"Error reading pixscale in the contour file!\n");
            exit(-1);
        }
         fgets(buffer,2048,fid);
        if(sscanf(buffer,"%d",&nobs)!=1)
          {
            fprintf(stderr,"Error reading the number of points in the contour file!\n");
            exit(-1);
        }
        lE=NORM(E);
        lE0=NORM(E0);
        (*C).datax[j]=calloc(nobs,sizeof(double));
        (*C).datay[j]=calloc(nobs,sizeof(double));
        (*C).nobs[j]=nobs;
        (*C).E[3*j]=E[0]/lE;
        (*C).E[3*j+1]=E[1]/lE;
        (*C).E[3*j+2]=E[2]/lE;
        (*C).E0[3*j]=E0[0]/lE0;
        (*C).E0[3*j+1]=E0[1]/lE0;
        (*C).E0[3*j+2]=E0[2]/lE0;
        if(rotate==1)
        {
         (*C).up[3*j]=0;
         (*C).up[3*j+1]=0;
         (*C).up[3*j+2]=1;
        }
        else if(rotate==0)
        {
         (*C).up[3*j+0]=0;
                 (*C).up[3*j+1]=0.397748474527011;
                 (*C).up[3*j+2]=0.917494496447491;
                 tilt=0;
        }
        else
        {
            fprintf(stderr,"Unknown rotate value in read_contour call\n");
            exit(-1);
        }
        (*C).TIME[j]=time-min_tim;
        (*C).distance[j]=lE;
        count=count+nobs;
        //Now we will read the actual contour points
        
        //Scale to km
        scale=1.0/(lE*pixscale*149597871)*180/PI*3600;
        
        if(type==0) //We are reading x-y coordinate pairs
        {
        for(int k=0;k<nobs;k++)
        {
            fgets(buffer,2048,fid);
            
           if(sscanf(buffer,"%lf %lf",&x,&y)!=2)
          {
            fprintf(stderr,"Error parsing contour point  in the contour file!(contour %d, index %d)\n",j,k);
            exit(-1);
          }
          (*C).datax[j][k]=(cos(tilt*deg2rad)*x-sin(tilt*deg2rad)*y)/scale;
          (*C).datay[j][k]=(cos(tilt*deg2rad)*y+sin(tilt*deg2rad)*x)/scale;
        }
        }
        else
        {
            for(int k=0;k<nobs;k++)
        {
            fgets(buffer,2048,fid);
           if(sscanf(buffer,"%lf %lf",&angle,&radius)!=2)
          {
            fprintf(stderr,"Error parsing (angle,radius) point  in the contour file!(contour %d, index %d)\n",j,k);
            exit(-1);
          }
          
          x=radius*cos(angle); //IN RADIANS?
          y=radius*sin(angle);
          (*C).datax[j][k]=(cos(tilt*deg2rad)*x-sin(tilt*deg2rad)*y)/scale;
          (*C).datay[j][k]=(cos(tilt*deg2rad)*y+sin(tilt*deg2rad)*x)/scale;
        }
        }
            
  }
   free(buffer);
  fclose(fid);
  (*C).ntotal=count;
  return C;
}
  double minv(double *vlist,int nvert,int index)
{
    double min=1E15;
    for(int j=0;j<nvert;j++)
    {
        if(vlist[j*3+index]<min)
            min=vlist[j*3+index];
    }
    return min;
}
double maxv(double *vlist,int nvert,int index)
{
    double max=-1E15;
    for(int j=0;j<nvert;j++)
    {
        if(vlist[j*3+index]>max)
            max=vlist[j*3+index];
    }
    return max;
}
int find_closest(double *p,double *plistx,double *plisty,int n)
{
    int index=-1;
    double dist=1e9,dist2;
    for(int j=0;j<n;j++)
    {
        dist2=pow(plistx[j]-p[0],2)+pow(plisty[j]-p[1],2);
        if(dist2<dist)
        {
            dist=dist2;
            index=j;
        }
    }
    return index;
}
void read_obj_file(char *file, int **tlist,double **vlist,int *nfac,int *nvert)
{
    FILE *fid;
    char *buffer;
    int *ttlist=calloc(20000,sizeof(int));
    double *tvlist=calloc(20000,sizeof(double));
    int vcount=0;
    int tcount=0;
    char t;
    buffer=calloc(2048,sizeof(char));
   fid=fopen(file,"r");
     if(fid==NULL)
    {
        fprintf(stderr,"Cannot open file %s for reading the shape in obj format\n",file);
        exit(-1);
    }
    fgets(buffer,2048,fid);
    while(buffer[0]=='v')
    {
        if(sscanf(buffer,"%c %lf %lf %lf",&t,tvlist+3*vcount,tvlist+3*vcount+1,tvlist+3*vcount+2)!=4)
        {
            fprintf(stderr,"Error processing file %s (line %s)\n",file,buffer);
            exit(-1);
        }
        vcount++;
        if(3*vcount>20000-1)
        {
           fprintf(stderr,"Vertex count too large in file %s\n",file);
            exit(-1); 
        }
        fgets(buffer,2048,fid);
    }
     while(buffer[0]=='f')
    {
        if(sscanf(buffer,"%c %d %d %d",&t,ttlist+3*tcount,ttlist+3*tcount+1,ttlist+3*tcount+2)!=4)
        {
            fprintf(stderr,"Error processing file %s (line %s)\n",file,buffer);
            exit(-1);
        }
        tcount++;
         if(3*tcount>20000-1)
        {
           fprintf(stderr,"facet count too large in file %s\n",file);
            exit(-1); 
        }
        if(fgets(buffer,2048,fid)==NULL)
            break;
    }
    int *tlist2=calloc(3*tcount,sizeof(int));
    double *vlist2=calloc(3*vcount,sizeof(double));
    memcpy(tlist2,ttlist,3*tcount*sizeof(int));
    memcpy(vlist2,tvlist,3*vcount*sizeof(double));
    *nfac=tcount;
    *nvert=vcount;
    *vlist=vlist2;
    *tlist=tlist2;
    free(ttlist);
    free(tvlist);
    free(buffer);
    fclose(fid);
}
int read_values_from_file(char *filename,double **fbuffer)
{ 
/*
     * Read a vector of double values from a file. Lines starting with # or ;
     *are ignored. Any delimeter( \t\r\n\f\v,) can be used to separate values.
     */ 
/*Difference to read_vector_file: fbuffer is allocated in this routine
 */

 FILE *fid;
    fid=fopen(filename,"r");
     char delims[]=" \t\r\n\f\v,";
     char *token;
     char *buff,*filebuff;
     buff=malloc(10000*sizeof(double));
     filebuff=malloc(10000*sizeof(double));
     double *buffer=malloc(10000*sizeof(double));
     if(fid==NULL)
    {
        perror("Error opening file in read_values_from_file");
        exit(-1);
    }
    int count=0;
    while(fgets(buff,10000,fid)!=NULL)
    {
        if(buff[0]=='#' || buff[0]==';')
            continue;
        memcpy(filebuff,buff,10000);
        
    token=strtok(filebuff,delims);
    
    while(token!=NULL && count<10000)
    {
        buffer[count]=atof(token);
        token=strtok(NULL,delims);
        count++;
    }
    }
    *fbuffer=calloc(count,sizeof(double));
    memcpy(*fbuffer,buffer,count*sizeof(double));
    
    free(buff);
    free(buffer);
    free(filebuff);
    return count;
}
NOaARstruct* init_NOaARstruct(int nfac)
{
    /*
     * Alloc struct for normals, areas and their derivatives
     */
    NOaARstruct *noaar=calloc(1,sizeof(NOaARstruct));
    noaar->normal=calloc(3*nfac,sizeof(double));
    noaar->area=calloc(3*nfac,sizeof(double));
    noaar->dndx1=calloc(3*nfac,sizeof(double));
    noaar->dndx2=calloc(3*nfac,sizeof(double));
    noaar->dndx3=calloc(3*nfac,sizeof(double));
    
    noaar->dndy1=calloc(3*nfac,sizeof(double));
    noaar->dndy2=calloc(3*nfac,sizeof(double));
    noaar->dndy3=calloc(3*nfac,sizeof(double));
    
    noaar->dndz1=calloc(3*nfac,sizeof(double));
    noaar->dndz2=calloc(3*nfac,sizeof(double));
    noaar->dndz3=calloc(3*nfac,sizeof(double));
    
    noaar->dadx=calloc(3*nfac,sizeof(double));
    noaar->dady=calloc(3*nfac,sizeof(double));
    noaar->dadz=calloc(3*nfac,sizeof(double));
    return noaar;
}
void free_NOaARstruct(NOaARstruct* noaar)
{
    free(noaar->normal);
    free(noaar->area);
    free(noaar->dndx1);
    free(noaar->dndy1);
    free(noaar->dndz1);
    free(noaar->dndx2);
    free(noaar->dndy2);
    free(noaar->dndz2);
    free(noaar->dndx3);
    free(noaar->dndy3);
    free(noaar->dndz3);
    free(noaar->dadx);
    free(noaar->dady);
    free(noaar->dadz);
    free(noaar);
}

void Find_Facets(int *tlist,double *vlist,int nfac,int nvert,int **Nfacets,int **Facetlist)
{
    /*
     * For each vertex, find facets that contain the vertex
     * OUTPUT:
     * Nfacets: 1xnvert vector, number of facets corresponding to the vertex
     * Facetlist: nvertx6 matrix, contains list of facets. CHECK IF 6 IS ENOUGH!!!!
     */
    int *Facets=calloc(nvert,sizeof(int));
    int *Facetl=calloc(nvert*6,sizeof(int));
    int *Temp=calloc(nvert*nfac,sizeof(int));
    int v1,v2,v3;
    for(int j=0;j<nfac;j++)
    {
        v1=tlist[3*j]-1;
        v2=tlist[3*j+1]-1;
        v3=tlist[3*j+2]-1;
        Temp[v1*nfac+j]=1;
        Temp[v2*nfac+j]=1;
        Temp[v3*nfac+j]=1;
    }
    for(int j=0;j<nvert;j++)
        for(int k=0;k<nfac;k++)
            if(Temp[j*nfac+k]==1)
            {
                Facetl[6*j+Facets[j]]=k;
                Facets[j]++;
            }
    free(Temp);
    *Nfacets=Facets;
    *Facetlist=Facetl;
}
void Calculate_Facet_Normals(int *tlist,double *vlist,int nfac,int nvert,int *Nfacets,int *Facetlist,double *Fnormals)
{
    /*
     * Calculate facet normals given the model and neighborhood information
     * Fnormals nvertx3 matrix, allocated before
     */

    double *normals=calloc(nfac*3,sizeof(double));
    int count=0;
    Calculate_Normals(tlist,vlist,nfac,nvert,normals);
    double *temp,norm;
    for(int j=0;j<nvert;j++)
    {
        for(int k=0;k<Nfacets[j];k++)
        {
            Fnormals[3*j]+=normals[3*Facetlist[6*j+k]];
            Fnormals[3*j+1]+=normals[3*Facetlist[6*j+k]+1];
            Fnormals[3*j+2]+=normals[3*Facetlist[6*j+k]+2];
            
        }
        temp=Fnormals+3*j;
        norm=NORM(temp);
        Fnormals[3*j]=Fnormals[3*j]/norm;
        Fnormals[3*j+1]=Fnormals[3*j+1]/norm;
        Fnormals[3*j+2]=Fnormals[3*j+2]/norm;
    }
    free(normals);
}
void Generate_Normal_Deriv_Matrix(int *tlist,double *vlist,int nfac,int nvert,int *Nfacets,int *Facetlist,double *D)
{
    /*
     * Generate a matrix where columns consist of vertex normals
     * D is 3nvertxnvert matrix, preallocated
     */
    memset(D,0,sizeof(double)*3*nvert*nvert);
    double *Fnormals=calloc(nvert*3,sizeof(double));
    Calculate_Facet_Normals(tlist,vlist,nfac,nvert,Nfacets,Facetlist,Fnormals);
    for(int j=0;j<nvert;j++)
    {
        D[j*nvert+j]=Fnormals[3*j];
        D[nvert*(nvert+j)+j]=Fnormals[3*j+1];
        D[nvert*(2*nvert+j)+j]=Fnormals[3*j+2];
    }
    free(Fnormals);
}
void Generate_Normal_Deriv_Matrix_Pad(int *tlist,double *vlist,int nfac,int nvert,int *Nfacets,int *Facetlist,double *D,int padding)
{
    /*
     * Generate a matrix where columns consist of vertex normals
     * D is (3nvert+padding)xnvert+padding matrix, preallocated
     * D=[D 0
     *    0 I]
     */
    memset(D,0,sizeof(double)*3*nvert*nvert);
    double *Fnormals=calloc(nvert*3,sizeof(double));
    Calculate_Facet_Normals(tlist,vlist,nfac,nvert,Nfacets,Facetlist,Fnormals);
    for(int j=0;j<nvert;j++)
    {
        D[j*(nvert+3)+j]=Fnormals[3*j];
        D[(nvert+3)*(nvert+j)+j]=Fnormals[3*j+1];
        D[(nvert+3)*(2*nvert+j)+j]=Fnormals[3*j+2];
    }
    if(padding>0)
        for(int j=0;j<padding;j++)
            set_el(D,3*nvert+padding,nvert+padding,1.0,3*nvert+j,nvert+j);
    free(Fnormals);
}
void Scale_Matrix_with_Vector(double *vec,double *M,int m,int n,double *N)
{
    /*
     * N(i,:)=vec(i)*M(i,:)
     */
    for(int i=0;i<m;i++)
        for(int j=0;j<n;j++)
            N[i*n+j]=vec[j]*M[i*n+j];
}
