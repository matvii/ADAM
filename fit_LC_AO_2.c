#include<stdio.h>
#include<stdlib.h>
#include"utils.h"
#include"matrix_ops.h"
#include"structs.h"
int main()
{
    char shapefile[]="shape2.txt";
    char lcfile[]="Metis_lc";
    double ini_dims[]={90,90,70};
    double angles[]={66*PI/180,185*PI/180,24*2*PI*1.0/5.079177,0};
    double min_tim=2449830.78281;
    int *tlist,*tlistn;
    double *vlist,*vlist2,*vlistn,*vlistn2,*D;
    double *angles2;
    angles2=calloc(4,sizeof(double));
    int nvert,nfac,nfacn,nvertn;
    LCstruct *LC;
    read_shape(shapefile,&tlist,&vlist,&nfac,&nvert,0); //Read an initial shape
    LC=read_lcurve(lcfile,min_tim);
    nfacn=9*nfac;
    nvertn=2+9*nfac/2; //Here we assume sdstep=2
    mul_cols(vlist,nvert,3,ini_dims); //scale the initial shape model
    
    //Process AO data
    char file1[]="Metis/2.fits";
    char file2[]="Metis/1.fits";
    int nx[]={150,150};
    int ny[]={150,150};
    char **AOfiles=calloc(2,sizeof(char*));
    AOfiles[0]=file1;
    AOfiles[1]=file2;
    int zerovec[]={0,0};
    int fiftyvec[]={51,51};
    double dates[2];
    dates[0]=NAN;
    dates[1]=NAN;
    double pixscale[]={0.009942,0.009942};
    char **psfnames=calloc(2,sizeof(char*));
    char ephmfile[]="Metis/ephm.dat";
    int nao=2;
    AOstruct *AO;
    int nAOobs=2,reads;
    double *E=calloc(nAOobs*3,sizeof(double));
    double *E0=calloc(nAOobs*3,sizeof(double));
    double *TIME=calloc(nAOobs,sizeof(double));
    
    double up[]={0,0.397748474527011,0.917494496447491,0,0.397748474527011,0.917494496447491};
    reads=read_ephm_data(ephmfile,TIME,E,E0);
    AO=process_ao_images(AOfiles,psfnames,nao,fiftyvec,fiftyvec,nx,ny,zerovec,zerovec,zerovec,zerovec,pixscale,pixscale,dates,min_tim,E,E0,up,TIME);
    ////////////////DEBUG
   // write_matrix_file("/tmp/data0r.txt",AO->datar[0],1,AO->nobs[0]);
   // printf("%d\n",AO->nobs[0]);
   // printf("last datar: %f\n",AO->datar[0][AO->nobs[0]-1]);
   // printf("last datai: %f\n",AO->datai[0][AO->nobs[0]-1]);
   // write_matrix_file("/tmp/data0i.txt",AO->datai[0],1,AO->nobs[0]);
   // write_matrix_file("/tmp/data1r.txt",AO->datar[1],1,AO->nobs[1]);
   // write_matrix_file("/tmp/data1i.txt",AO->datai[1],1,AO->nobs[1]);
    //write_matrix_file("/tmp/freqx.txt",AO->freqx[0],1,AO->nobs[0]);
   // write_matrix_file("/tmp/freqy.txt",AO->freqy[0],1,AO->nobs[0]);
    /////////////////DEBUG
    free(E);
    free(E0);
    free(TIME);
   
    int nAOtotal=2*(AO->ntotal);
    int nAOcols=3*nvert+3+2*(AO->nao);
    int nAOoffsets=2*(AO->nao);
    double *offset=calloc(2*(AO->nao),sizeof(double));
    ///////////
    //double offset[]={0.0029,-0.0130,-0.0049,-0.0149};
    ///////////////
    double *offset2=calloc(2*(AO->nao),sizeof(double));
    double *AOout=calloc(nAOtotal,sizeof(double));
    double *AOdv=calloc(nAOtotal*nAOcols,sizeof(double));
    
    //
    //Variables for regularization terms
    double CRres; //Convex reg
    double *dCRdv;
    double *y,*y2;
    
    vlist2=calloc(3*nvert,sizeof(double));
   
    dCRdv=calloc(3*nvert+3,sizeof(double));
    
    double ANGres; //Dihedral angle regularization
    double *dANGdv;
   
    
    dANGdv=calloc(3*nvert+3,sizeof(double));
    
    double *Ares; //Area regularization
    double *dAdv;
   
    Ares=calloc(nfacn,sizeof(double));
    
    dAdv=calloc(nfacn*(3*nvert+3),sizeof(double));
    
    
    double *S,*S2;
    double *J,*JTJ,*JTJpd;
    
    int tp=0; 
    int Slength;
    int nLCtotal=LC->ntotal;
    tp=nLCtotal+nAOtotal+1+1+nfacn; //Number of rows in J
    Slength=nLCtotal+nAOtotal+1+1+nfacn; //LC points+AO points+convex reg+dihedral reg+angle reg
    int nJcols=3*nvert+3+nAOoffsets;
    S=calloc(Slength,sizeof(double)); 
    J=calloc(Slength*(nJcols),sizeof(double));
   
    JTJ=calloc((nJcols)*(nJcols),sizeof(double));
    JTJpd=calloc((nJcols)*(nJcols),sizeof(double));
    double *LCout,*dLCdv,*rhs,*X;
    LCout=calloc(nLCtotal,sizeof(double));
    dLCdv=calloc((nLCtotal)*(3*nvert+3),sizeof(double));
    rhs=calloc(nJcols,sizeof(double));
    X=calloc(nJcols,sizeof(double));
    //Weights
    double AOfit=0;
    double AOW=2;
    double lcW=2;
    double cW=50;
    double aW=50;
    double angW=2;
    double lambda=1;
    int decreased=1;
    double chisq=1e9;
    double chisq2;
    int NumofRounds=50;
   
    for(int k=0;k<NumofRounds;k++)
    {
        if(decreased==1)
        {
        //REMEMBER TO ZERO MATRICES
        Sqrt3_Subdiv(tlist,vlist,nfac,nvert,&tlistn,&vlistn,&nfacn,&nvertn,&D); //Do 2x subdivision. Remember to free allocated mem
        //Calculate LCs
       
        calculate_lcs(tlistn,vlistn,nfacn,nvertn,angles,LC,D,nvertn,nvert,LCout,dLCdv,NULL,NULL,NULL,1);
      //AO data
        Calculate_AOs(tlistn,vlistn,nfacn,nvertn,angles,AO,offset,D,nvertn,nvert,NULL,AOout,AOdv,1);
        //
        ///////////////////////////DEBUG
       // write_matrix_file("/tmp/AOout.txt",AOout,1,nAOtotal);
       // write_matrix_file("/tmp/AOdv.txt",AOdv,nAOtotal,nAOcols);
        ///////////////////////////DEBUG
        convex_reg(tlistn,vlistn,nfacn,nvertn,D,nvertn,nvert,&CRres,dCRdv);
        
        dihedral_angle_reg(tlistn,vlistn,nfacn,nvertn,D,nvertn,nvert,&ANGres,dANGdv);
        
        area_reg(tlistn,vlistn,nfacn,nvertn,D,nvertn,nvert,Ares,dAdv); 
       //REMEMBER TO CHECK CORRECTNESS OF AREA_REG
      
        //Multiply with weights;
        mult_with_cons(AOout,1,nAOtotal,AOW);
        mult_with_cons(AOdv,nAOtotal,nAOcols,AOW);
        mult_with_cons(LCout,1,nLCtotal,lcW);
        mult_with_cons(dLCdv,nLCtotal,3*nvert+3,lcW);
        CRres*=cW;
        
        mult_with_cons(dCRdv,1,3*nvert+3,cW);
        ANGres*=angW;
        mult_with_cons(dANGdv,1,3*nvert+3,angW);
        mult_with_cons(Ares,1,nfacn,aW);
        mult_with_cons(dAdv,nfacn,3*nvert+3,aW);
        
        //Remember to free and realloc matrices
        //Build the res vector and matrix
        set_submatrix(S,1,Slength,LCout,1,nLCtotal,0,0);
        set_submatrix(S,1,Slength,AOout,1,nAOtotal,0,nLCtotal);
        S[nLCtotal+nAOtotal]=CRres; //Convex reg terms
        S[nLCtotal+nAOtotal+1]=ANGres;
        set_submatrix(S,1,Slength,Ares,1,nfacn,0,nLCtotal+nAOtotal+2); //Area regularization
       
        matrix_transprod(S,Slength,1,&chisq);
      
        printf("chisq: %f\n",chisq); 
       
        //Construct Jacobian matrix
        
        set_submatrix(J,Slength,nJcols,dLCdv,nLCtotal,3*nvert+3,0,0);
        set_submatrix(J,Slength,nJcols,AOdv,nAOtotal,nAOcols,nLCtotal,0);
        set_submatrix(J,Slength,nJcols,dCRdv,1,3*nvert+3,nLCtotal+nAOtotal,0);
        set_submatrix(J,Slength,nJcols,dANGdv,1,3*nvert+3,nLCtotal+1+nAOtotal,0);
        set_submatrix(J,Slength,nJcols,dAdv,nfacn,3*nvert+3,nLCtotal+2+nAOtotal,0);
        
        
        //Calculate J^TJ+lambda*J
       // matrix_transprod(J,Slength,nJcols,JTJ);
         free(tlistn);
        free(vlistn);
        free(D);
        
        
        matrix_transprod(J,Slength,nJcols,JTJ);
         matrix_vectorprod(J,Slength,nJcols,S,rhs,1); //rhs=J^T*S;
        }
       matrix_adddiag(JTJ,JTJpd,nJcols,lambda); //JTJpd=JTJ+lambda*diag(JTJ)
      
       
        /////////////////////DEBUG
       // write_matrix_file("/tmp/S.txt",S,Slength,1);
       //  write_matrix_file("/tmp/J.txt",J,Slength,nJcols);
        //  write_matrix_file("/tmp/JTJpd.txt",JTJpd,nJcols,nJcols);
        //  write_matrix_file("/tmp/rhs.txt",rhs,nJcols,1);
          ///////////////////DEBUG
        solve_matrix_eq(JTJpd,nJcols,rhs,X); //Solve the LM 
        //////////////////DEBUG
       // write_matrix_file("/tmp/X.txt",X,1,nJcols);
        ////////////////DEBUG
    
        add_vector_to_vlist(vlist,X,vlist2,nvert);
        matrix_plus2(offset,1,nAOoffsets,&X[3*nvert+3],offset2);
       
        angles2[0]=angles[0]+X[3*nvert];
        angles2[1]=angles[1]+X[3*nvert+1];
        angles2[2]=angles[2]+X[3*nvert+2];
    
      
        //New shape is vlist2. Check for fit
       
        
        Sqrt3_Subdiv(tlist,vlist2,nfac,nvert,&tlistn,&vlistn,&nfacn,&nvertn,&D);
        
        calculate_lcs(tlistn,vlistn,nfacn,nvertn,angles2,LC,D,nvertn,nvert,LCout,dLCdv,NULL,NULL,NULL,0);
        Calculate_AOs(tlistn,vlistn,nfacn,nvertn,angles2,AO,offset2,D,nvertn,nvert,NULL,AOout,AOdv,0);
        convex_reg(tlistn,vlistn,nfacn,nvertn,D,nvertn,nvert,&CRres,dCRdv);
        dihedral_angle_reg(tlistn,vlistn,nfacn,nvertn,D,nvertn,nvert,&ANGres,dANGdv);
        area_reg(tlistn,vlistn,nfacn,nvertn,D,nvertn,nvert,Ares,dAdv); 
        mult_with_cons(LCout,1,nLCtotal,lcW);
        mult_with_cons(AOout,1,nAOtotal,AOW);
        CRres*=cW;
        
        ANGres*=angW;
        
        mult_with_cons(Ares,1,nfacn,aW);
        
       set_submatrix(S,1,Slength,LCout,1,nLCtotal,0,0);
        set_submatrix(S,1,Slength,AOout,1,nAOtotal,0,nLCtotal);
        S[nLCtotal+nAOtotal]=CRres; //Convex reg terms
        S[nLCtotal+nAOtotal+1]=ANGres;
        set_submatrix(S,1,Slength,Ares,1,nfacn,0,nLCtotal+nAOtotal+2); //Area regularization
       
        matrix_transprod(S,Slength,1,&chisq2);
         matrix_transprod(AOout,nAOtotal,1,&AOfit);
        printf("chisq2: %f\n",chisq2);
      printf("k=%d\n",k);
       if(chisq2<chisq)
       {
           printf("decreased!\n");
           memcpy(vlist,vlist2,sizeof(double)*(3*nvert));
           memcpy(offset,offset2,sizeof(double)*(nAOoffsets));
           printf("Offsets: %f %f %f %f\n",offset[0],offset[1],offset[2],offset[3]);
           decreased=1;
           angles[0]=angles2[0];
           angles[1]=angles2[1];
           angles[2]=angles2[2];
           lambda=0.1*lambda;
           zero_array(J,Slength*nJcols);
           zero_array(S,Slength);
       }
       else
       {
           lambda=10*lambda;
           decreased=0;
          
    }
    if(k==NumofRounds-1)
    {
        write_shape_file("/tmp/shape.txt",tlistn,vlistn,nfacn,nvertn);
        printf("Angles: %.4f %.4f %.8f\n",angles[0],angles[1],angles[2]);
    }
    free(vlistn);
    free(tlistn);
    free(D);
//     zero_array(dLCdv,nLCtotal*(3*nvert+3));
//     zero_array(AOdv,nAOcols*nAOtotal);
//     zero_array(dCRdv,3*nvert+3);
//     zero_array(dANGdv,3*nvert+3);
//     zero_array(dAdv,nfacn*(3*nvert+3));
    
}
free(tlist);
free(vlist);
free(dLCdv);
free(dCRdv);
free(dANGdv);
free(AOdv);
free(dAdv);
free(vlist2);
free(angles2);
free(Ares);
free(S);
free(J);
free(JTJ);
free(JTJpd);
free(LCout);
free(rhs);
free(X);
free_lc_struct(LC);
}
