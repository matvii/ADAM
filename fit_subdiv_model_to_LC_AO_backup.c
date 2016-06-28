#include<stdio.h>
#include<stdlib.h>
#include"utils.h"
#include"matrix_ops.h"
#include"structs.h"
#include"globals.h"
void fit_subdiv_model_to_LC_AO(LCstruct *LC,AOstruct *AO)
{
    //First initialize the initial shape
    int *tlist,*tlistn,nfac,nvert,nfacn,nvertn;
    double *vlist,*vlist2,*vlistn,*D;
    if(INI_SHAPE_FILE!=NULL)
        read_shape(INI_SHAPE_FILE,&tlist,&vlist,&nfac,&nvert,0);
    else
    {
        nfac=8*pow(INI_NROWS,2);
        nvert=4*pow(INI_NROWS,2)+2;
        vlist=calloc(3*nvert,sizeof(double));
        tlist=calloc(3*nfac,sizeof(int));
        generate_ellipsoid(INI_NROWS,ELLIPSOID_SEMI_A,ELLIPSOID_SEMI_B,ELLIPSOID_SEMI_C,tlist,vlist);
    }
   
    vlist2=calloc(3*nvert,sizeof(double));
    //Setup angles
    double angles[]={(90-INI_ANGLE_B)*PI/180,INI_ANGLE_L*PI/180,24*2*PI*1.0/INI_ANGLE_P,INI_ANGLE_PHI0*PI/180};
    double angles2[4]={0,0,0,0};
     Sqrt3_Subdiv(tlist,vlist,nfac,nvert,&tlistn,&vlistn,&nfacn,&nvertn,&D,INI_SD_LEVEL);
   
    //This is for calibrated lightcurves
    double *params,*params2;
   if(INI_PARAMS!=NULL)
   {
       params=calloc(4,sizeof(double));
       params2=calloc(4,sizeof(double));
       params[0]=INI_PARAMS[0];
       params[1]=INI_PARAMS[1];
       params[2]=INI_PARAMS[2];
       params[3]=INI_PARAMS[3];
   }
    int USE_CALIB=LC->calib;
    int ncalib=0;
    if(USE_CALIB==1)
        ncalib=3;
    //This is for calibrated lightcurves end
    int nAOtotal=0;
    int nAOcols=0;
    int nAO=0;
    int nAOoffsets=0;
    int nAOscale=0;
    double *AOoffset,*AOoffset2,*AOout,*AOdv,AOfit,*AOscale,*AOscale2,*AOds;
    if(INI_HAVE_AO)
    {
        nAOtotal=2*(AO->ntotal);
        nAOcols=3*nvert+3+2*(AO->nao);
        nAOoffsets=2*(AO->nao);
        nAO=AO->nao;
        AOoffset=calloc(2*(AO->nao),sizeof(double));
        
        AOoffset2=calloc(2*(AO->nao),sizeof(double));
        AOout=calloc(nAOtotal,sizeof(double));
        AOdv=calloc(nAOtotal*nAOcols,sizeof(double));
        AOfit=0;
        AOscale=NULL;
        AOscale2=NULL;
        
        
        if(INI_AO_SCALING)
        {
            AOscale=calloc(1*nAO,sizeof(double));
            AOscale2=calloc(1*nAO,sizeof(double));
            AOds=calloc(nAOtotal*nAO,sizeof(double));
            nAOscale=nAO;
        }
    }
    
    //Variables for regularization terms
    double CRres; //Convex reg
    double *dCRdv;
    double *y,*y2;
//      write_matrix_file("/tmp/data0r.txt",AO->datar[0],1,AO->nobs[0]);
//      write_matrix_file("/tmp/data0i.txt",AO->datai[0],1,AO->nobs[0]);
//      write_matrix_file("/tmp/data1r.txt",AO->datar[1],1,AO->nobs[1]);
//      write_matrix_file("/tmp/data1i.txt",AO->datai[1],1,AO->nobs[1]);
//      write_matrix_file("/tmp/freqx.txt",AO->freqx[0],1,AO->nobs[0]);
//      write_matrix_file("/tmp/freqy.txt",AO->freqy[0],1,AO->nobs[0]);
//     write_matrix_file("/tmp/psf0r.txt",AO->psfr[0],1,AO->nobs[0]);
//      write_matrix_file("/tmp/psf0i.txt",AO->psfi[0],1,AO->nobs[0]);
//      write_matrix_file("/tmp/psf1r.txt",AO->psfr[1],1,AO->nobs[1]);
//      write_matrix_file("/tmp/psf1i.txt",AO->psfi[1],1,AO->nobs[1]);
   
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
    int nJcols=3*nvert+3+ncalib+nAOoffsets+nAOscale;
    S=calloc(Slength,sizeof(double)); 
    J=calloc(Slength*(nJcols),sizeof(double));
   
    JTJ=calloc((nJcols)*(nJcols),sizeof(double));
    JTJpd=calloc((nJcols)*(nJcols),sizeof(double));
    double *LCout,*dLCdv,*rhs,*X;
    double LCfit;
    LCout=calloc(nLCtotal,sizeof(double));
    dLCdv=calloc((nLCtotal)*(3*nvert+ncalib+3),sizeof(double));
    rhs=calloc(nJcols,sizeof(double));
    X=calloc(nJcols,sizeof(double));
    //phase parameters
    double *dLCdp=NULL;
    if(USE_CALIB==1)
        dLCdp=calloc(nLCtotal*3,sizeof(double));
    
    //Weights
   
    double AOW=INI_AOW;
    double lcW=INI_LCW;
    double cW=INI_CW;
    double aW=INI_AW;
    double angW=INI_DW;
    double lambda=INI_LAMBDA;
    int decreased=1;
    double chisq=1e9;
    double chisq2;
    double Aresfit;
   // NUM_OF_ROUNDS=1;
  
    for(int k=0;k<NUM_OF_ROUNDS;k++)
    {
        if(decreased==1)
        {
        
        Sqrt3_Subdiv(tlist,vlist,nfac,nvert,&tlistn,&vlistn,&nfacn,&nvertn,&D,INI_SD_LEVEL); //Do  subdivision. Remember to free allocated mem
        //Calculate LCs
//         write_shape_file("/tmp/Inishape.txt",tlistn,vlistn,nfacn,nvertn);
        calculate_lcs(tlistn,vlistn,nfacn,nvertn,angles,LC,D,nvertn,nvert,LCout,dLCdv,NULL,NULL,NULL,params,dLCdp,1);
//        write_matrix_file("/tmp/LCout.txt",LCout,1,nLCtotal);
//         write_matrix_file("/tmp/dLCdv.txt",dLCdv,nLCtotal,3*nvert+3);
        //AO data
        if(INI_HAVE_AO)
        {
            Calculate_AOs(tlistn,vlistn,nfacn,nvertn,angles,AO,AOoffset,D,nvertn,nvert,NULL,AOscale,AOout,AOdv,AOds,1);
            mult_with_cons(AOout,1,nAOtotal,AOW);
            mult_with_cons(AOdv,nAOtotal,nAOcols,AOW);
            set_submatrix(S,1,Slength,AOout,1,nAOtotal,0,nLCtotal); //Construct Jacobian matrix
            set_submatrix(J,Slength,nJcols,AOdv,nAOtotal,nAOcols,nLCtotal,0);
            matrix_transprod(AOout,nAOtotal,1,&AOfit);
            if(INI_AO_SCALING)
            {
                mult_with_cons(AOds,nAOtotal,nAO,AOW);
                set_submatrix(J,Slength,nJcols,AOds,nAOtotal,nAO,nLCtotal,nAOcols);
            }
        }
//           write_matrix_file("/tmp/AOout.txt",AOout,1,nAOtotal);
//           write_matrix_file("/tmp/AOdv.txt",AOdv,nAOtotal,nAOcols);
        convex_reg(tlistn,vlistn,nfacn,nvertn,D,nvertn,nvert,&CRres,dCRdv);
    
        dihedral_angle_reg(tlistn,vlistn,nfacn,nvertn,D,nvertn,nvert,&ANGres,dANGdv);
      
        area_reg(tlistn,vlistn,nfacn,nvertn,D,nvertn,nvert,Ares,dAdv); 

        
        mult_with_cons(LCout,1,nLCtotal,lcW);
        mult_with_cons(dLCdv,nLCtotal,3*nvert+3,lcW);
        CRres*=cW;
        
        mult_with_cons(dCRdv,1,3*nvert+3,-cW);
        ANGres*=angW;
        mult_with_cons(dANGdv,1,3*nvert+3,-angW);
        mult_with_cons(Ares,1,nfacn,aW);
        mult_with_cons(dAdv,nfacn,3*nvert+3,-aW);
        
        
        //Build the res vector and matrix
        set_submatrix(S,1,Slength,LCout,1,nLCtotal,0,0);
       
        
       
        S[nLCtotal+nAOtotal]=CRres; //Convex reg terms
        S[nLCtotal+nAOtotal+1]=ANGres;
       
        set_submatrix(S,1,Slength,Ares,1,nfacn,0,nLCtotal+nAOtotal+2); //Area regularization
        
       
       matrix_transprod(LCout,nLCtotal,1,&LCfit);
        matrix_transprod(S,Slength,1,&chisq);
       matrix_transprod(Ares,nfacn,1,&Aresfit);
        printf("chisq: %f LCfit: %f AOfit: %f Convex reg: %f Dihedral Angle reg: %f Area reg:%f\n",chisq,LCfit,AOfit,pow(CRres,2),pow(ANGres,2),Aresfit); 
       
        //Construct Jacobian matrix
        
        set_submatrix(J,Slength,nJcols,dLCdv,nLCtotal,3*nvert+3,0,0);
        
        set_submatrix(J,Slength,nJcols,dCRdv,1,3*nvert+3,nLCtotal+nAOtotal,0);
        set_submatrix(J,Slength,nJcols,dANGdv,1,3*nvert+3,nLCtotal+1+nAOtotal,0);
        set_submatrix(J,Slength,nJcols,dAdv,nfacn,3*nvert+3,nLCtotal+2+nAOtotal,0);
        
        //include lc parameter derivatives
        if(USE_CALIB==1)
        set_submatrix(J,Slength,nJcols,dLCdp,nLCtotal,3,0,nJcols-3);
        
        //Calculate J^TJ+lambda*J
       // matrix_transprod(J,Slength,nJcols,JTJ);
         free(tlistn);
        free(vlistn);
        free(D);
//         write_matrix_file("/tmp/S.txt",S,Slength,1);
//          write_matrix_file("/tmp/J.txt",J,Slength,nJcols);
        matrix_transprod(J,Slength,nJcols,JTJ);
         matrix_vectorprod(J,Slength,nJcols,S,rhs,1); //rhs=J^T*S;
        }
       matrix_adddiag(JTJ,JTJpd,nJcols,lambda); //JTJpd=JTJ+lambda*diag(JTJ)
      
       solve_matrix_eq(JTJpd,nJcols,rhs,X); //Solve the LM 
 
    
        add_vector_to_vlist(vlist,X,vlist2,nvert);
        if(INI_HAVE_AO)
        {
        matrix_plus2(AOoffset,1,nAOoffsets,&X[3*nvert+3],AOoffset2);
        if(INI_AO_SCALING)
           matrix_plus2(AOscale,1,nAO,&X[3*nvert+3+nAOoffsets],AOscale2);
        }
        angles2[0]=angles[0]+X[3*nvert];
        angles2[1]=angles[1]+X[3*nvert+1];
        angles2[2]=angles[2]+X[3*nvert+2];
        if(USE_CALIB==1)
        {
            params2[0]=params[0]+X[nJcols-3];
            params2[1]=params[1]+X[nJcols-2];
            params2[2]=params[2]+X[nJcols-1];
            params2[3]=params[3];
        }
    
      
        //New shape is vlist2. Check for fit
       
        
        Sqrt3_Subdiv(tlist,vlist2,nfac,nvert,&tlistn,&vlistn,&nfacn,&nvertn,&D,INI_SD_LEVEL);
        
        calculate_lcs(tlistn,vlistn,nfacn,nvertn,angles2,LC,D,nvertn,nvert,LCout,dLCdv,NULL,NULL,NULL,params2,dLCdp,0);
        if(INI_HAVE_AO)
        {
        Calculate_AOs(tlistn,vlistn,nfacn,nvertn,angles2,AO,AOoffset2,D,nvertn,nvert,NULL,AOscale2,AOout,AOdv,NULL,0);
        mult_with_cons(AOout,1,nAOtotal,AOW);
        set_submatrix(S,1,Slength,AOout,1,nAOtotal,0,nLCtotal);
        }
         
        convex_reg(tlistn,vlistn,nfacn,nvertn,D,nvertn,nvert,&CRres,dCRdv);
        dihedral_angle_reg(tlistn,vlistn,nfacn,nvertn,D,nvertn,nvert,&ANGres,dANGdv);
        area_reg(tlistn,vlistn,nfacn,nvertn,D,nvertn,nvert,Ares,dAdv); 
        mult_with_cons(LCout,1,nLCtotal,lcW);
        
        CRres*=cW;
        
        ANGres*=angW;
        
        mult_with_cons(Ares,1,nfacn,aW);
        
       set_submatrix(S,1,Slength,LCout,1,nLCtotal,0,0);
       
        S[nLCtotal+nAOtotal]=CRres; //Convex reg terms
        S[nLCtotal+nAOtotal+1]=ANGres;
        set_submatrix(S,1,Slength,Ares,1,nfacn,0,nLCtotal+nAOtotal+2); //Area regularization
       
        matrix_transprod(S,Slength,1,&chisq2);
         matrix_transprod(AOout,nAOtotal,1,&AOfit);
        printf("Round: %d chisq2: %f \n",k+1,chisq2);
      //printf("k=%d\n",k);
       if(chisq2<chisq)
       {
           printf("decreased!\n");
        //   printf("scales: %f %f\n",AOscale2[0],AOscale2[1]);
           memcpy(vlist,vlist2,sizeof(double)*(3*nvert));
           
          // printf("Offsets: %f %f %f %f\n",AOoffset[0],AOoffset[1],AOoffset[2],AOoffset[3]);
           decreased=1;
           angles[0]=angles2[0];
           angles[1]=angles2[1];
           angles[2]=angles2[2];
           if(INI_HAVE_AO)
           {
               memcpy(AOoffset,AOoffset2,sizeof(double)*(nAOoffsets));
           if(INI_AO_SCALING)
              memcpy(AOscale,AOscale2,sizeof(double)*nAO);
           }
           lambda=0.1*lambda;
           zero_array(J,Slength*nJcols);
           zero_array(S,Slength);
           if(USE_CALIB==1)
               memcpy(params,params2,3*sizeof(double));
       }
       else
       {
           lambda=10*lambda;
           decreased=0;
          
    }
    if(k==NUM_OF_ROUNDS-1)
    {
        printf("Angles: %.4f %.4f %.8f\n",angles[0],angles[1],angles[2]);
       // printf("Offsets: %f %f %f %f\n",AOoffset[0],AOoffset[1],AOoffset[2],AOoffset[3]);
        write_shape_file(OUT_SHAPE_FILE,tlistn,vlistn,nfacn,nvertn);
        FILE* fid=fopen(OUT_PARAM_FILE,"w");
        FILE* fidlc=fopen(OUT_LC_FILE,"w");
        if(fid!=NULL)
        {
            fprintf(fid,"%.4f %.4f %.8f\n",90-angles[0]*180/PI,angles[1]*180/PI,24*2*PI/angles[2]);
            fclose(fid);
        }
        else
        {
            fprintf(stderr,"Cannot open %s\n",OUT_PARAM_FILE);
            printf("Angles: %.4f %.4f %.8f\n",90-angles[0]*180/PI,angles[1]*180/PI,24*2*PI/angles[2]);
        }
        if(fidlc!=NULL)
        {
            for(int j=0;j<LC->nlc;j++)
                zero_array(LC->lcs[j],LC->nobs[j]);
            Sqrt3_Subdiv(tlist,vlist,nfac,nvert,&tlistn,&vlistn,&nfacn,&nvertn,&D,INI_SD_LEVEL);
            calculate_lcs(tlistn,vlistn,nfacn,nvertn,angles,LC,D,nvertn,nvert,LCout,dLCdv,NULL,NULL,NULL,params,dLCdp,0);
            for(int j=0;j<LC->ntotal;j++)
                fprintf(fidlc,"%f ",LCout[j]);
            fclose(fidlc);
        }
            
    }
    free(vlistn);
    free(tlistn);
    if(INI_SD_LEVEL!=0)
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
//free(AOdv);
free(dAdv);
free(vlist2);

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
