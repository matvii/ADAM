#include"utils.h"
#include"matrix_ops.h"
#include"globals.h"
void dihedral_angle_reg(int *tlist,double *vlist,int nfac,int nvert,double *D,int dm,int dn,double *result,double *drsdv)
{
    /*OUTPUT:
     * result a double
     * dresdx 1x3*nvert+3 array
     * regularization term: 1/(n1.n2)-1 for adjacent facets with normals n1 and n2
     */
    double mul=2;
    double cos_max_angle=0;
    double amul=1;
    double p=1;
    double sgn=0;
    double tres=0;
    if(INI_DIA_PARAMS)
    {
        mul=INI_DIA_PARAMS[2];
        cos_max_angle=cos(INI_DIA_PARAMS[1]*2*PI/360.0);
        amul=INI_DIA_PARAMS[0];
    }
    *result=0;
   int nedge=nfac+nvert-2;
   double *res,*drdx,*drdy,*drdz;
   int *EV;
   res=calloc(nedge,sizeof(double));
   drdx=calloc(nedge*nvert,sizeof(double));
   drdy=calloc(nedge*nvert,sizeof(double));
   drdz=calloc(nedge*nvert,sizeof(double));
//   
   if(res==NULL || drdx==NULL || drdy==NULL || drdz==NULL)
   {
       fprintf(stderr,"Calloc returned null\n");
       exit(1);
   }
   double *drsdx,*drsdy,*drsdz;
   
   double *dresdx,*dresdy,*dresdz;
   EV=malloc(nvert*nvert*sizeof(int));
   dihedral_angle(tlist,vlist,nfac,nvert,res,EV,drdx,drdy,drdz);
   if(D!=NULL && dm!=nvert)
    {
        puts("Error: Number of vertex coordinates is not equal to the number of rows in D.");
        exit(1);
    }
    if(D==NULL)
        dn=nvert;
    drsdx=calloc(nedge*dn,sizeof(double));
   drsdy=calloc(nedge*dn,sizeof(double));
   drsdz=calloc(nedge*dn,sizeof(double));
    if(D==NULL)
   for(int j=0;j<nedge;j++)
   {
       
       
       sgn=1.0;
       if(res[j]<0)
           sgn=-1.0;
       
       (*result)+=(sgn*1/res[j]-1);
       for(int k=0;k<nvert;k++)
       {
           
           drsdx[k]+=-sgn*1/pow(res[j],2)*drdx[j*nvert+k];
           drsdy[k]+=-sgn*1/pow(res[j],2)*drdy[j*nvert+k];
           drsdz[k]+=-sgn*1/pow(res[j],2)*drdz[j*nvert+k];
       }
   }
   else
   {
       dresdx=calloc(nvert,sizeof(double));
       dresdy=calloc(nvert,sizeof(double));
       dresdz=calloc(nvert,sizeof(double));
    for(int j=0;j<nedge;j++)
   {
       if(res[j]<0)
           sgn=-1.0;
       
       (*result)+=sgn*1/res[j]-1;
       for(int k=0;k<nvert;k++)
       {
          
           dresdx[k]+=-sgn*1/pow(res[j],2)*drdx[j*nvert+k];
           dresdy[k]+=-sgn*1/pow(res[j],2)*drdy[j*nvert+k];
           dresdz[k]+=-sgn*1/pow(res[j],2)*drdz[j*nvert+k];
       }
   }
   matrix_prod(dresdx,1,nvert,D,dn,drsdx);
   matrix_prod(dresdy,1,nvert,D,dn,drsdy);
   matrix_prod(dresdz,1,nvert,D,dn,drsdz);
   free(dresdx);
   free(dresdy);
   free(dresdz);
   }
   set_submatrix(drsdv,1,3*dn+3,drsdx,1,dn,0,0);
   set_submatrix(drsdv,1,3*dn+3,drsdy,1,dn,0,dn);
   set_submatrix(drsdv,1,3*dn+3,drsdz,1,dn,0,2*dn);
   free(drsdx);
   free(drsdy);
   free(drsdz);
   free(res);
   free(drdx);
   free(drdy);
   free(drdz);
   free(EV);
}
 /* 
  
int main()
 {
     //int tlist[]={1,2,3,1,3,4};//2x3
     //double vlist[]={-1,0,1,0,-1,0,0,1,0,0,10,10}; //4x3
//      int tlist[]={1,2,3,
//          1,3,4,
//          1,4,5,
//          1,5,2,
//          6,3,2,
//          6,4,3,
//          6,5,4,
//          6,2,5}; //8x3
//      double vlist[]={10,0,1,
//          1,0,0,
//          0,1,0,
//          -1,0,0,
//          0,-1,0,
//          0,0,-1}; //6x3
      int nfac=8,nfacn;
      int nvert=6,nvertn;
     int *tlist;
     double *vlist;
     double *vlistn;
     int *tlistn;
     double *D;
     char shapefile[]="shape2.txt";
     read_shape(shapefile,&tlist,&vlist,&nfac,&nvert,0);
      double ini_dims[]={90,90,70};
      mul_cols(vlist,nvert,3,ini_dims);
     Sqrt3_Subdiv(tlist,vlist,nfac,nvert,&tlistn,&vlistn,&nfacn,&nvertn,&D);
     int nedgen=nfacn+nvertn-2;
     int nedge=nfac+nvert-2;
     double res,*dresdv;
   
     int *EV=calloc(nvertn*nvertn,sizeof(int));
     
     double *res2=calloc(nedgen,sizeof(double));
     double *dresdx=calloc(nedgen*nvertn,sizeof(double));
      double *dresdy=calloc(nedgen*nvertn,sizeof(double));
       double *dresdz=calloc(nedgen*nvertn,sizeof(double));
    
  //  D=calloc(nvert*nvert,sizeof(double));
  //  for(int j=0;j<6;j++)
  //      D[j+nvert*j]=1;
     dresdv=calloc(3*nvertn+3,sizeof(double));
    //dihedral_angle(tlistn,vlistn,nfacn, nvertn,res2,EV,dresdx,dresdy,dresdz);
     dihedral_angle_reg(tlistn,vlistn,nfacn,nvertn,D,nvertn,nvert,&res,dresdv);
    printf("res: %.10f\n",res);
     //print_matrix(dresdv,1,3*nvert);
  //print_matrixI(EV,nvert,nvert);
//     write_matrix_file("/tmp/res.txt",res2,1,nedgen);
     write_matrix_file("/tmp/dresdv.txt",dresdv,1,3*nvert);
//      write_matrix_file("/tmp/dresdy.txt",dresdy,nedgen,nvertn);
//      write_matrix_file("/tmp/dresdz.txt",dresdz,nedgen,nvertn);
     
   free(dresdv);
   free(D);
    
 }
  
*/
