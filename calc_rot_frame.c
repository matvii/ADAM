#include"utils.h"

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
    if(abs(angle1-angle2)>0.01)
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
   return 1;
}
 /*   
void main(int argc,char** argv)
{
    double E[]={1.0,2.0,-1.0};
    double no=NORM(E);
    E[0]=E[0]/no;
    E[1]=E[1]/no;
    E[2]=E[2]/no;
    double up[3];
    calc_rot_frame(argv[1],E,up);
    print_matrix(up,1,3);
}
*/
