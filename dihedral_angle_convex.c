#include"utils.h"
void dihedral_angle_convex(int* tlist,double* vlist,int nfac,int nvert,double *res,int *angletype,int *EV,double *dresdx,double *dresdy,double *dresdz)
{
    /*For each edge calculate the (cosine of) dihedral angle, ie the angle between the normals of
     * adjacent facets
     * EV is the matrix indexing the fedges, ie if EV(i,j)=k>-1 then res(k)=
     * dihedral angle of edge between vertices (i,j) (and (j,i))
     * OUTPUT:
     * res nfac+nvert-2 array
     * angletype int nfac+nvert-2 array, 1 if angle is convex, 0 otherwise
     * EV nvertxnvert matrix
     * dresdx nedgexnvert matrix */
    
    int nedge;
    int *E,*N,*E2,*A;
    double *v1,*v2,*v3,n1[3],n2[3];
    double *w1,*w2,*w3;
    double dn1dx1[3],dn1dx2[3],dn1dx3[3];
    double dn1dy1[3],dn1dy2[3],dn1dy3[3];
    double dn1dz1[3],dn1dz2[3],dn1dz3[3];
    
    double dn2dx1[3],dn2dx2[3],dn2dx3[3];
    double dn2dy1[3],dn2dy2[3],dn2dy3[3];
    double dn2dz1[3],dn2dz2[3],dn2dz3[3];
    double area;
    double dAdx[3],dAdy[3],dAdz[3];
    double c1[3],c2[3],c[3];
    double anglesign;
    E=calloc(nvert*nvert,sizeof(int));
     N=calloc(nvert*nfac,sizeof(int));
     E2=calloc(nvert*nvert,sizeof(int));
     A=calloc(nfac*nfac,sizeof(int));
    find_neighborhood(tlist,vlist,nfac,nvert,E,N,E2,A);
   
    nedge=nfac+nvert-2;
    int ind;
    int f1,f2;
    int i1,i2,i3;
    int j1,j2,j3;
    int count=0;
    //Initialize EV
    for(int j=0;j<nvert*nvert;j++)
        EV[j]=-1;
    //Matrix to index edges
    for(int j=0;j<nvert;j++)
        for(int k=j;k<nvert;k++)
        {
            if(get_elI(E,nvert,nvert,j,k)==0)
                continue;
            f1=get_elI(E2,nvert,nvert,j,k);
            f2=get_elI(E2,nvert,nvert,k,j);
           
         //   printf("Count is %d\n",count);
            set_elI(EV,nvert,nvert,count,j,k);
            set_elI(EV,nvert,nvert,count,k,j);
            
            //Calculate normals and their derivatives
            i1=tlist[3*f1]-1;
            i2=tlist[3*f1+1]-1;
            i3=tlist[3*f1+2]-1;
            j1=tlist[3*f2]-1;
            j2=tlist[3*f2+1]-1;
            j3=tlist[3*f2+2]-1;
            v1=vlist+3*i1;
            v2=vlist+3*i2;
            v3=vlist+3*i3;
            w1=vlist+3*j1;
            w2=vlist+3*j2;
            w3=vlist+3*j3;
            Calculate_Area_and_Normal_Derivative(v1,v2,v3,n1,dn1dx1,dn1dx2,dn1dx3,dn1dy1,dn1dy2,dn1dy3,dn1dz1,dn1dz2,dn1dz3,&area,dAdx,dAdy,dAdz);
            Calculate_Area_and_Normal_Derivative(w1,w2,w3,n2,dn2dx1,dn2dx2,dn2dx3,dn2dy1,dn2dy2,dn2dy3,dn2dz1,dn2dz2,dn2dz3,&area,dAdx,dAdy,dAdz);
            //We don't care about areas here
            res[count]=DOT(n1,n2);
           c1[0]=(v1[0]+v2[0]+v3[0]);
           c1[1]=(v1[1]+v2[1]+v3[1]);
           c1[2]=(v1[2]+v2[2]+v3[2]);
           c2[0]=(w1[0]+w2[0]+w3[0]);
           c2[1]=(w1[1]+w2[1]+w3[1]);
           c2[2]=(w1[2]+w2[2]+w3[2]);
           c[0]=c2[0]-c1[0];
           c[1]=c2[1]-c1[1];
           c[2]=c2[2]-c1[2];
           anglesign=n1[0]*c[0]+n1[1]*c[1]+n1[2]*c[2];
           angletype[count]=0;
           if(anglesign<0)
               angletype[count]=1;
            dresdx[count*nvert+i1]=DOT(dn1dx1,n2);
            dresdy[count*nvert+i1]=DOT(dn1dy1,n2);
            dresdz[count*nvert+i1]=DOT(dn1dz1,n2);
            
            dresdx[count*nvert+i2]=DOT(dn1dx2,n2);
            dresdy[count*nvert+i2]=DOT(dn1dy2,n2);
            dresdz[count*nvert+i2]=DOT(dn1dz2,n2);
            
            dresdx[count*nvert+i3]=DOT(dn1dx3,n2);
            dresdy[count*nvert+i3]=DOT(dn1dy3,n2);
            dresdz[count*nvert+i3]=DOT(dn1dz3,n2);
            /**********************/
            dresdx[count*nvert+j1]+=DOT(dn2dx1,n1);
            dresdy[count*nvert+j1]+=DOT(dn2dy1,n1);
            dresdz[count*nvert+j1]+=DOT(dn2dz1,n1);
            
            dresdx[count*nvert+j2]+=DOT(dn2dx2,n1);
            dresdy[count*nvert+j2]+=DOT(dn2dy2,n1);
            dresdz[count*nvert+j2]+=DOT(dn2dz2,n1);
            
            dresdx[count*nvert+j3]+=DOT(dn2dx3,n1);
            dresdy[count*nvert+j3]+=DOT(dn2dy3,n1);
            dresdz[count*nvert+j3]+=DOT(dn2dz3,n1);
            count++;
        }
        free(E);
        free(E2);
        free(N);
        free(A);
}
/*
int main()
 {
     //int tlist[]={1,2,3,1,3,4};//2x3
     //double vlist[]={-1,0,1,0,-1,0,0,1,0,0,10,10}; //4x3
     int tlist[]={1,2,3,1,3,4,1,4,5,1,5,2,6,3,2,6,4,3,6,5,4,6,2,5}; //8x3
     double vlist[]={0,0,1,1,0,0,0,1,0,-1,0,0,0,-1,0,0,0,-1}; //6x3
     int nfac=8;
     int nvert=6;
     int nedge=nfac+nvert-2;
     double *res,*dresdx,*dresdy,*dresdz;
     int *EV;
     res=calloc(nedge,sizeof(double));
     dresdx=calloc(nedge*nvert,sizeof(double));
     dresdy=calloc(nedge*nvert,sizeof(double));
     dresdz=calloc(nedge*nvert,sizeof(double));
     EV=calloc(nvert*nvert,sizeof(int));
     dihedral_angle(tlist,vlist,nfac,nvert,res,EV,dresdx,dresdy,dresdz);
     print_matrix(res,1,nedge);
     print_matrix(dresdx,nedge,nvert);
    
 }
    
*/
