#include"utils.h"
#include"matrix_ops.h"
#include"num_of_threads.h"
#include"globals.h"
void Fit_Contour(int *tlist,double *vlist,int nfac,int nvert,double *angles,double *up,double *E,double *E0,double TIME,double *offset,double *datax,double *datay,int npoints,double *dist,double *dx,double *dy,double *dz,double *dangles,double *dtox,double *dtoy);
void Fit_Contour2(int *tlist,double *vlist,int nfac,int nvert,double *angles,double *up,double *E,double *E0,double TIME,double *offset,double *datax,double *datay,int npoints,double *dist,double *dx,double *dy,double *dz,double *dangles,double *dtox,double *dtoy);
void Calculate_Contours(int *tlist,double *vlist,int nfac,int nvert,double *angles,CNTRstruct* CR,double *offset,double *D,int dm,int dn,double* CRdist,double* dOdv,double *dOdoff)
{
    /*Construct derivative matrix wrt shape parameters corresponding to chord interserctions*/ 
    /*offset is 2*noc vector containing offsets
     * Chordoffset is ntotal vector containing offset of each chord in seconds (optional)
     * D optional derivative matrix by which the original derivative matrix is multiplied
     * dOdv is 4*ntotal x (3*nvert+3) (or (3*dn+3) if D!=NULL) matrix containing derivatives
     * dOdoff 4*ntotal x 2*noc matrix for derivatives wrt offset terms
     * dChordoff 4*notalxntotal matrix for derivatives wrt chord offsets
     */
    int noc=CR->ncont;
    int *nobs=CR->nobs;
    
    int *cumcount=calloc(noc+1,sizeof(int));
    if(D!=NULL && dm!=nvert)
    {
        fprintf(stderr,"nvert and dm must be equal if D is non-null\n");
        exit(-1);
    }
    if(D==NULL)
        dn=nvert;
    
    cumcount[0]=0;
    int ntotal;
    if(INI_CNTR_IS_SPARSE)
    {
        for(int j=1;j<=noc;j++)
            cumcount[j]=cumcount[j-1]+nobs[j-1];
        ntotal=CR->ntotal;
    }
     
        else
        {
          for(int j=1;j<=noc;j++)
            cumcount[j]=cumcount[j-1]+nobs[j-1]+nvert;  
          ntotal=CR->ntotal+nvert*CR->ncont;
        }
    
    zero_array(dOdv,ntotal*(3*dn+3));
    zero_array(dOdoff,ntotal*2*noc);
     
    for(int j=0;j<noc;j++)
    {
       
        double *dx=calloc((nobs[j]+nvert)*nvert,sizeof(double));
        double *dy=calloc((nobs[j]+nvert)*nvert,sizeof(double));
        double *dz=calloc((nobs[j]+nvert)*nvert,sizeof(double));
        double *dangles=calloc((nobs[j]+nvert)*3,sizeof(double));
        double *dtox=calloc((nobs[j]+nvert),sizeof(double));
        double *dtoy=calloc((nobs[j]+nvert),sizeof(double));
        double *dist=calloc(nobs[j]+nvert,sizeof(double));
        if(INI_CNTR_RAD)
       Fit_Contour2(tlist,vlist,nfac,nvert,angles,CR->up+3*j,CR->E+3*j,CR->E0+3*j,CR->TIME[j],offset+2*j,CR->datax[j],CR->datay[j],nobs[j],dist,dx,dy,dz,dangles,dtox,dtoy);
       else
        Fit_Contour(tlist,vlist,nfac,nvert,angles,CR->up+3*j,CR->E+3*j,CR->E0+3*j,CR->TIME[j],offset+2*j,CR->datax[j],CR->datay[j],nobs[j],dist,dx,dy,dz,dangles,dtox,dtoy);
        
        
        if(D!=NULL)
        {
            double *dx2=calloc((nobs[j]+nvert)*dn,sizeof(double));
            double *dy2=calloc((nobs[j]+nvert)*dn,sizeof(double));
            double *dz2=calloc((nobs[j]+nvert)*dn,sizeof(double));
            
            matrix_prod(dx,(nobs[j]+nvert),nvert,D,dn,dx2);
            matrix_prod(dy,(nobs[j]+nvert),nvert,D,dn,dy2);
            matrix_prod(dz,(nobs[j]+nvert),nvert,D,dn,dz2);
            if(INI_CNTR_IS_SPARSE)
            {
                memcpy(CRdist+cumcount[j],dist,nobs[j]*sizeof(double));
                set_submatrix(dOdv,ntotal,3*dn+3,dx2,nobs[j],dn,cumcount[j],0);
             
            set_submatrix(dOdv,ntotal,3*dn+3,dy2,nobs[j],dn,cumcount[j],dn);
            set_submatrix(dOdv,ntotal,3*dn+3,dz2,nobs[j],dn,cumcount[j],2*dn);
           
            set_submatrix(dOdv,ntotal,3*dn+3,dangles,nobs[j],3,cumcount[j],3*dn);
        
            set_submatrix(dOdoff,ntotal,2*noc,dtox,nobs[j],1,cumcount[j],2*j);
            set_submatrix(dOdoff,ntotal,2*noc,dtoy,nobs[j],1,cumcount[j],2*j+1);
            }
            else
            {
              memcpy(CRdist+cumcount[j],dist,(nobs[j]+nvert)*sizeof(double));  
            set_submatrix(dOdv,ntotal,3*dn+3,dx2,nvert+nobs[j],dn,cumcount[j],0);
             
            set_submatrix(dOdv,ntotal,3*dn+3,dy2,nvert+nobs[j],dn,cumcount[j],dn);
            set_submatrix(dOdv,ntotal,3*dn+3,dz2,nvert+nobs[j],dn,cumcount[j],2*dn);
           
            set_submatrix(dOdv,ntotal,3*dn+3,dangles,nvert+nobs[j],3,cumcount[j],3*dn);
        
            set_submatrix(dOdoff,ntotal,2*noc,dtox,nvert+nobs[j],1,cumcount[j],2*j);
            set_submatrix(dOdoff,ntotal,2*noc,dtoy,nvert+nobs[j],1,cumcount[j],2*j+1);
            }
            
            free(dx);
            free(dy);
            free(dz);
            free(dangles);
            free(dtox);
            free(dtoy);
            free(dx2);
            free(dy2);
            free(dz2);
            free(dist);
        }
        else
        {
            if(INI_CNTR_IS_SPARSE)
            {
                memcpy(CRdist+cumcount[j],dist,nobs[j]*sizeof(double));
                set_submatrix(dOdv,ntotal,3*nvert+3,dx,nobs[j],nvert,cumcount[j],0);
            set_submatrix(dOdv,ntotal,3*nvert+3,dy,nobs[j],nvert,cumcount[j],nvert);
            set_submatrix(dOdv,ntotal,3*nvert+3,dz,nobs[j],nvert,cumcount[j],2*nvert);
            
            set_submatrix(dOdv,ntotal,3*nvert+3,dangles,nobs[j],3,cumcount[j],3*nvert);
            
            set_submatrix(dOdoff,ntotal,2*noc,dtox,nobs[j],1,cumcount[j],2*j);
            set_submatrix(dOdoff,ntotal,2*noc,dtoy,nobs[j],1,cumcount[j],2*j+1);
            }
            else
            {
                memcpy(CRdist+cumcount[j],dist,(nobs[j]+nvert)*sizeof(double));  
           set_submatrix(dOdv,ntotal,3*nvert+3,dx,nobs[j],nvert,cumcount[j],0);
            set_submatrix(dOdv,ntotal,3*nvert+3,dy,nobs[j],nvert,cumcount[j],nvert);
            set_submatrix(dOdv,ntotal,3*nvert+3,dz,nobs[j],nvert,cumcount[j],2*nvert);
            
            set_submatrix(dOdv,ntotal,3*nvert+3,dangles,nobs[j],3,cumcount[j],3*nvert);
            
            set_submatrix(dOdoff,ntotal,2*noc,dtox,nobs[j],1,cumcount[j],2*j);
            set_submatrix(dOdoff,ntotal,2*noc,dtoy,nobs[j],1,cumcount[j],2*j+1);
            }
            free(dx);
            free(dy);
            free(dz);
            free(dangles);
            free(dtox);
            free(dtoy); 
            free(dist);
        }
        
    }
    free(cumcount);
}
/*
void main()
{
        int nfac,nfacn;
        int nvert,nvertn;
        OCstruct *OC;
        double *vlist,*vlistn;
        int *tlist,*tlistn;
        //read_shape("/tmp/Inishape.txt",&tlist,&vlist,&nfac,&nvert,0);
        read_shape("/tmp/testshape.txt",&tlist,&vlist,&nfac,&nvert,0);
        double *D;
        Sqrt3_Subdiv(tlist,vlist,nfac,nvert,&tlistn,&vlistn,&nfacn,&nvertn,&D,2);
    
        double angles[]={(90-52)*PI/180,272*PI/180,24*2*PI*1.0/8.40060,0};
        double min_tim=2443846.0;
        double offset[]={0,0};
        char filename[]="/tmp/Hertha.occ";
        double up[]={0,0.397748474527011,0.917494496447491};
        OC=read_occ(filename,min_tim,up);
      
        int ntotal=OC->ntotal;
        double *OCdist=calloc(4*ntotal,sizeof(double));
        double *dOdv=calloc(4*ntotal*(3*nvertn+3),sizeof(double));
        double *dOdoff=calloc(4*ntotal*2,sizeof(double));
         calculate_OCs(tlistn,vlistn,nfacn,nvertn,angles,OC,offset,NULL,nvert,nvert,OCdist,dOdv,dOdoff);
       write_matrix_file("/tmp/OCdist.txt",OCdist,4*ntotal,1);
        write_matrix_file("/tmp/dOdv.txt",dOdv,4*ntotal,3*nvertn+3);
        write_matrix_file("/tmp/dOdoff.txt",dOdoff,4*ntotal,2);
        write_matrix_file("/tmp/D.txt",D,nvertn,nvert);
         double *dOdv2=calloc(4*ntotal*(3*nvert+3),sizeof(double));
        calculate_OCs(tlistn,vlistn,nfacn,nvertn,angles,OC,offset,D,nvertn,nvert,OCdist,dOdv2,dOdoff);
         write_matrix_file("/tmp/dOdv2.txt",dOdv2,4*ntotal,3*nvert+3);
         write_matrix_file("/tmp/chords1.txt",OC->data[0],OC->nobs[0],4);
         write_matrix_file("/tmp/OCdist2.txt",OCdist,4*ntotal,1);
         write_matrix_file("/tmp/times1.txt",OC->TIME[0],OC->nobs[0],2);
        
        print_matrix(OC->E,1,3);
       
        print_matrix(angles,1,4);
        print_matrix(OC->up,1,3);
}
*/
