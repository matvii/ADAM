#include"utils.h"
#include"matrix_ops.h"
double center_of_mass(int *tlist,double *vlist,int nfac,int nvert,double *D,int dm,int dn,double *dv,int deriv)
{
    /*
     * dCdx 3xnvert matrix
     * dv 1x3*nvert+3 matrix if D is null
     * otherwise 1x3*dn+3 matrix
     * returns the distance of the center of mass from zero, truncated to 0 if <1e-2
     * 
     */
    double dndx1[3],dndx2[3],dndx3[3];
    double dndy1[3],dndy2[3],dndy3[3];
    double dndz1[3],dndz2[3],dndz3[3];
    double com[3];
    double normal[3],A=0,V=0;
    double cx=0.0,cy=0.0,cz=0.0;
    double dAdx[3],dAdy[3],dAdz[3];
    double *dVdx,*dVdy,*dVdz;
    double dpxdx1,dpxdx2,dpxdx3;
    double dpydy1,dpydy2,dpydy3;
    double dpzdz1,dpzdz2,dpzdz3;
    dVdx=calloc(nvert,sizeof(double));
    dVdy=calloc(nvert,sizeof(double));
    dVdz=calloc(nvert,sizeof(double));
    double *dcdx,*dcdy,*dcdz;
    dcdx=calloc(3*nvert,sizeof(double));
    dcdy=calloc(3*nvert,sizeof(double));
    dcdz=calloc(3*nvert,sizeof(double));
    double *v1,*v2,*v3;
    double tpx,tpy,tpz;
    int i1,i2,i3;
    for(int j=0;j<nfac;j++)
    {
        i1=tlist[3*j]-1;
        i2=tlist[3*j+1]-1;
        i3=tlist[3*j+2]-1;
        v1=vlist+3*i1;
        v2=vlist+3*i2;
        v3=vlist+3*i3;
        Calculate_Area_and_Normal_Derivative(v1,v2,v3,normal,dndx1,dndx2,dndx3,dndy1,dndy2,dndy3,dndz1,dndz2,dndz3,&A,dAdx,dAdy,dAdz);
        V=V+1.0/3.0*(v1[0]*normal[0]+v1[1]*normal[1]+v1[2]*normal[2])*A;
        tpx=(pow(v1[0]+v2[0],2)+pow(v2[0]+v3[0],2)+pow(v1[0]+v3[0],2));
        tpy=(pow(v1[1]+v2[1],2)+pow(v2[1]+v3[1],2)+pow(v1[1]+v3[1],2));
        tpz=(pow(v1[2]+v2[2],2)+pow(v2[2]+v3[2],2)+pow(v1[2]+v3[2],2));
        cx=cx+2*A*normal[0]*tpx;
        cy=cy+2*A*normal[1]*tpy;
        cz=cz+2*A*normal[2]*tpz;
        if(deriv>0)
        {
        dVdx[i1]+=1.0/3*((normal[0]+v1[0]*dndx1[0]+v1[1]*dndx1[1]+v1[2]*dndx1[2])*A+(v1[0]*normal[0]+v1[1]*normal[1]+v1[2]*normal[2])*dAdx[0]);
        dVdx[i2]+=1.0/3*((v1[0]*dndx2[0]+v1[1]*dndx2[1]+v1[2]*dndx2[2])*A+(v1[0]*normal[0]+v1[1]*normal[1]+v1[2]*normal[2])*dAdx[1]);
        dVdx[i3]+=1.0/3*((v1[0]*dndx3[0]+v1[1]*dndx3[1]+v1[2]*dndx3[2])*A+(v1[0]*normal[0]+v1[1]*normal[1]+v1[2]*normal[2])*dAdx[2]);
        
        dVdy[i1]+=1.0/3*((normal[1]+v1[0]*dndy1[0]+v1[1]*dndy1[1]+v1[2]*dndy1[2])*A+(v1[0]*normal[0]+v1[1]*normal[1]+v1[2]*normal[2])*dAdy[0]);
        dVdy[i2]+=1.0/3*((v1[0]*dndy2[0]+v1[1]*dndy2[1]+v1[2]*dndy2[2])*A+(v1[0]*normal[0]+v1[1]*normal[1]+v1[2]*normal[2])*dAdy[1]);
        dVdy[i3]+=1.0/3*((v1[0]*dndy3[0]+v1[1]*dndy3[1]+v1[2]*dndy3[2])*A+(v1[0]*normal[0]+v1[1]*normal[1]+v1[2]*normal[2])*dAdy[2]);
        
       
        dVdz[i1]+=1.0/3*((normal[2]+v1[0]*dndz1[0]+v1[1]*dndz1[1]+v1[2]*dndz1[2])*A+(v1[0]*normal[0]+v1[1]*normal[1]+v1[2]*normal[2])*dAdz[0]);
        dVdz[i2]+=1.0/3*((v1[0]*dndz2[0]+v1[1]*dndz2[1]+v1[2]*dndz2[2])*A+(v1[0]*normal[0]+v1[1]*normal[1]+v1[2]*normal[2])*dAdz[1]);
        dVdz[i3]+=1.0/3*((v1[0]*dndz3[0]+v1[1]*dndz3[1]+v1[2]*dndz3[2])*A+(v1[0]*normal[0]+v1[1]*normal[1]+v1[2]*normal[2])*dAdz[2]);
        
        
        
        dpxdx1=2*(v1[0]+v2[0])+2*(v1[0]+v3[0]);
        dpxdx2=2*(v1[0]+v2[0])+2*(v2[0]+v3[0]);
        dpxdx3=2*(v2[0]+v3[0])+2*(v1[0]+v3[0]);
        
        dpydy1=2*(v1[1]+v2[1])+2*(v1[1]+v3[1]);
        dpydy2=2*(v1[1]+v2[1])+2*(v2[1]+v3[1]);
        dpydy3=2*(v2[1]+v3[1])+2*(v1[1]+v3[1]);
        
        dpzdz1=2*(v1[2]+v2[2])+2*(v1[2]+v3[2]);
        dpzdz2=2*(v1[2]+v2[2])+2*(v2[2]+v3[2]);
        dpzdz3=2*(v2[2]+v3[2])+2*(v1[2]+v3[2]);
        
        
        //Column of dcdx is [dcxdx1,dcydx1,dczdx1]'
        dcdx[0*nvert+i1]+=(2*dAdx[0]*normal[0]+2*A*dndx1[0])*tpx+2*A*normal[0]*dpxdx1;
        dcdx[0*nvert+i2]+=(2*dAdx[1]*normal[0]+2*A*dndx2[0])*tpx+2*A*normal[0]*dpxdx2;
        dcdx[0*nvert+i3]+=(2*dAdx[2]*normal[0]+2*A*dndx3[0])*tpx+2*A*normal[0]*dpxdx3;
        
        dcdx[1*nvert+i1]+=(2*dAdx[0]*normal[1]+2*A*dndx1[1])*tpy;
        dcdx[1*nvert+i2]+=(2*dAdx[1]*normal[1]+2*A*dndx2[1])*tpy;
        dcdx[1*nvert+i3]+=(2*dAdx[2]*normal[1]+2*A*dndx3[1])*tpy;

        dcdx[2*nvert+i1]+=(2*dAdx[0]*normal[2]+2*A*dndx1[2])*tpz;
        dcdx[2*nvert+i2]+=(2*dAdx[1]*normal[2]+2*A*dndx2[2])*tpz;
        dcdx[2*nvert+i3]+=(2*dAdx[2]*normal[2]+2*A*dndx3[2])*tpz;
        
        dcdy[0*nvert+i1]+=(2*dAdy[0]*normal[0]+2*A*dndy1[0])*tpx;
        dcdy[0*nvert+i2]+=(2*dAdy[1]*normal[0]+2*A*dndy2[0])*tpx;
        dcdy[0*nvert+i3]+=(2*dAdy[2]*normal[0]+2*A*dndy3[0])*tpx;
        
        dcdy[1*nvert+i1]+=(2*dAdy[0]*normal[1]+2*A*dndy1[1])*tpy+2*A*normal[1]*dpydy1;
        dcdy[1*nvert+i2]+=(2*dAdy[1]*normal[1]+2*A*dndy2[1])*tpy+2*A*normal[1]*dpydy2;
        dcdy[1*nvert+i3]+=(2*dAdy[2]*normal[1]+2*A*dndy3[1])*tpy+2*A*normal[1]*dpydy3;
        
        dcdy[2*nvert+i1]+=(2*dAdy[0]*normal[2]+2*A*dndy1[2])*tpz;
        dcdy[2*nvert+i2]+=(2*dAdy[1]*normal[2]+2*A*dndy2[2])*tpz;
        dcdy[2*nvert+i3]+=(2*dAdy[2]*normal[2]+2*A*dndy3[2])*tpz;
        
       
        dcdz[0*nvert+i1]+=(2*dAdz[0]*normal[0]+2*A*dndz1[0])*tpx;
        dcdz[0*nvert+i2]+=(2*dAdz[1]*normal[0]+2*A*dndz2[0])*tpx;
        dcdz[0*nvert+i3]+=(2*dAdz[2]*normal[0]+2*A*dndz3[0])*tpx;
        
        dcdz[1*nvert+i1]+=(2*dAdz[0]*normal[1]+2*A*dndz1[1])*tpy;
        dcdz[1*nvert+i2]+=(2*dAdz[1]*normal[1]+2*A*dndz2[1])*tpy;
        dcdz[1*nvert+i3]+=(2*dAdz[2]*normal[1]+2*A*dndz3[1])*tpy;
        
        dcdz[2*nvert+i1]+=(2*dAdz[0]*normal[2]+2*A*dndz1[2])*tpz+2*A*normal[2]*dpzdz1;
        dcdz[2*nvert+i2]+=(2*dAdz[1]*normal[2]+2*A*dndz2[2])*tpz+2*A*normal[2]*dpzdz2;
        dcdz[2*nvert+i3]+=(2*dAdz[2]*normal[2]+2*A*dndz3[2])*tpz+2*A*normal[2]*dpzdz3;
        }
    }
    
    com[0]=1.0/48.0*1/V*cx;
    com[1]=1.0/48.0*1/V*cy;
    com[2]=1.0/48.0*1/V*cz;
    double V2=pow(V,2),C;
    C=NORM(com);
    if(deriv>0)
        if(D==NULL)
        {
            for(int j=0;j<nvert;j++)
            {
                dcdx[0*nvert+j]=1.0/48.0*(V*dcdx[0*nvert+j]-dVdx[j]*cx)/V2;
                dcdx[1*nvert+j]=1.0/48.0*(V*dcdx[1*nvert+j]-dVdx[j]*cy)/V2;
                dcdx[2*nvert+j]=1.0/48.0*(V*dcdx[2*nvert+j]-dVdx[j]*cz)/V2;
                
                dcdy[0*nvert+j]=1.0/48.0*(V*dcdy[0*nvert+j]-dVdy[j]*cx)/V2;
                dcdy[1*nvert+j]=1.0/48.0*(V*dcdy[1*nvert+j]-dVdy[j]*cy)/V2;
                dcdy[2*nvert+j]=1.0/48.0*(V*dcdy[2*nvert+j]-dVdy[j]*cz)/V2;
                
                dcdz[0*nvert+j]=1.0/48.0*(V*dcdz[0*nvert+j]-dVdz[j]*cx)/V2;
                dcdz[1*nvert+j]=1.0/48.0*(V*dcdz[1*nvert+j]-dVdz[j]*cy)/V2;
                dcdz[2*nvert+j]=1.0/48.0*(V*dcdz[2*nvert+j]-dVdz[j]*cz)/V2;
                
                dv[j]=1/C*(com[0]*dcdx[j]+com[1]*dcdx[nvert+j]+com[2]*dcdx[2*nvert+j]);
                dv[nvert+j]=1/C*(com[0]*dcdy[j]+com[1]*dcdy[nvert+j]+com[2]*dcdy[2*nvert+j]);
                dv[2*nvert+j]=1/C*(com[0]*dcdz[j]+com[1]*dcdz[nvert+j]+com[2]*dcdz[2*nvert+j]);
            }
            if(C<1e-2)
            {
                C=0;
                zero_array(dv,3*nvert+3);
                
            }
        }
        else
        {
            double *dvdx=calloc(nvert,sizeof(double));
            double *dvdy=calloc(nvert,sizeof(double));
            double *dvdz=calloc(nvert,sizeof(double));
            double *dvdt=calloc(dn,sizeof(double));
            for(int j=0;j<nvert;j++)
            {
                dcdx[0*nvert+j]=1.0/48.0*(V*dcdx[0*nvert+j]-dVdx[j]*cx)/V2;
                dcdx[1*nvert+j]=1.0/48.0*(V*dcdx[1*nvert+j]-dVdx[j]*cy)/V2;
                dcdx[2*nvert+j]=1.0/48.0*(V*dcdx[2*nvert+j]-dVdx[j]*cz)/V2;
                
                dcdy[0*nvert+j]=1.0/48.0*(V*dcdy[0*nvert+j]-dVdy[j]*cx)/V2;
                dcdy[1*nvert+j]=1.0/48.0*(V*dcdy[1*nvert+j]-dVdy[j]*cy)/V2;
                dcdy[2*nvert+j]=1.0/48.0*(V*dcdy[2*nvert+j]-dVdy[j]*cz)/V2;
                
                dcdz[0*nvert+j]=1.0/48.0*(V*dcdz[0*nvert+j]-dVdz[j]*cx)/V2;
                dcdz[1*nvert+j]=1.0/48.0*(V*dcdz[1*nvert+j]-dVdz[j]*cy)/V2;
                dcdz[2*nvert+j]=1.0/48.0*(V*dcdz[2*nvert+j]-dVdz[j]*cz)/V2;
                
                dvdx[j]=1/C*(com[0]*dcdx[j]+com[1]*dcdx[nvert+j]+com[2]*dcdx[2*nvert+j]);
                dvdy[j]=1/C*(com[0]*dcdy[j]+com[1]*dcdy[nvert+j]+com[2]*dcdy[2*nvert+j]);
                dvdz[j]=1/C*(com[0]*dcdz[j]+com[1]*dcdz[nvert+j]+com[2]*dcdz[2*nvert+j]);
            }
            matrix_prod(dvdx,1,nvert,D,dn,dvdt);
            set_submatrix(dv,1,3*dn+3,dvdt,1,dn,0,0);
            matrix_prod(dvdy,1,nvert,D,dn,dvdt);
            set_submatrix(dv,1,3*dn+3,dvdt,1,dn,0,dn);
            matrix_prod(dvdz,1,nvert,D,dn,dvdt);
            set_submatrix(dv,1,3*dn+3,dvdt,1,dn,0,2*dn);
            free(dvdx);
            free(dvdy);
            free(dvdz);
            free(dvdt);
            if(C<1e-2)
            {
                C=0;
                zero_array(dv,3*dn+3);
                
            }
        }
        return C;
        free(dcdx);
        free(dcdy);
        free(dcdz);
        free(dVdy);
        free(dVdx);
        free(dVdz);
        
}
