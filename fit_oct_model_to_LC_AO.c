#include<stdio.h>
#include<stdlib.h>
#include"utils.h"
#include"matrix_ops.h"
#include"structs.h"
#include"globals.h"
void fit_oct_model_to_LC_AO(LCstruct *LC,AOstruct *AO,OCstruct *OC,RDstruct *RD)
{
    //First initialize the initial shape
    int *tlist;
    double *vlist,*vlist2,*D1,*D2;
   
    int nfac=8*pow(INI_NROWS,2);
    int nvert=4*pow(INI_NROWS,2)+2;
        vlist=calloc(3*nvert,sizeof(double));
        vlist2=calloc(3*nvert,sizeof(double));
        tlist=calloc(3*nfac,sizeof(int));
    double *a,*a2;
    int alength=3*pow(INI_LMAX+1,2);
    int al=pow(INI_LMAX+1,2);
    
    a=calloc(alength,sizeof(double));
    a2=calloc(alength,sizeof(double));
    a[0]=log(ELLIPSOID_SEMI_A)/0.2821;
    a[1*al]=log(ELLIPSOID_SEMI_B)/0.2821-a[0];
    a[2*al]=log(ELLIPSOID_SEMI_C)/0.2821-a[0];
    
    //Setup angles
    double angles[]={(90-INI_ANGLE_B)*PI/180,INI_ANGLE_L*PI/180,24*2*PI*1.0/INI_ANGLE_P,INI_ANGLE_PHI0*PI/180};
    double angles2[4]={0,0,0,0};
   angles2[3]=angles[3];
    //This is for calibrated lightcurves
    double *params=NULL,*params2=NULL;
    if(INI_PHASE_PARAMS!=NULL)
    {
        params=calloc(4,sizeof(double));
        params2=calloc(4,sizeof(double));
        params[0]=INI_PHASE_PARAMS[0];
        params[1]=INI_PHASE_PARAMS[1];
        params[2]=INI_PHASE_PARAMS[2];
        params[3]=INI_PHASE_PARAMS[3];
        memcpy(params2,params,sizeof(double)*4);
    }
     int nLCtotal=LC->ntotal;
    int USE_CALIB=LC->calib;
    int ncalib=0;
    if(USE_CALIB==1)
        ncalib=4;
    //This is for calibrated lightcurves end
    //LC albedo variegation
    int nAlbedo=0;
    double *eAlbedo=NULL,*eAlbedo2=NULL;
    double *dLCdalb=NULL;
    double *Alblimits=NULL;
    double *Albreg=NULL,*dAlbreg=NULL;
    int nAlbreg=0;
    if(INI_FIT_ALBEDO==1)
    {
        nAlbedo=nfac;
        eAlbedo=calloc(nfac,sizeof(double));
        eAlbedo2=calloc(nfac,sizeof(double));
        dLCdalb=calloc(nLCtotal*nfac,sizeof(double));
        Alblimits=calloc(2,sizeof(double));
        Alblimits[0]=INI_ALBEDO_MIN;
        Alblimits[1]=INI_ALBEDO_MAX;
        Albreg=calloc(nfac,sizeof(double));
        dAlbreg=calloc(nfac*nfac,sizeof(double));
        nAlbreg=nfac;
    }
    
   double ANGres;
   double *dANGdv=calloc(3*nvert+3,sizeof(double));
   double *dANGda=calloc(alength+3,sizeof(double));
    
   
   int nAOtotal=0;
    int nAOcols=0;
    int nAOcolst=0;
    int nAO=0;
    int nAOoffsets=0;
    int nAOscale=0;
    double *AOoffset,*AOoffset2,*AOout,*AOda,*AOdv,AOfit,*AOscale,*AOscale2,*AOds;
    if(INI_HAVE_AO)
    {
        nAOtotal=2*(AO->ntotal);
        nAOcolst=3*nvert+3+2*(AO->nao);
        nAOoffsets=2*(AO->nao);
        nAOcols=alength+3+2*(AO->nao);
        nAO=AO->nao;
        AOoffset=calloc(2*(AO->nao),sizeof(double));
        if(INI_INPUT_AO_OFFSET!=NULL)
        {
            read_vector_file(INI_INPUT_AO_OFFSET,AOoffset,2*(AO->nao));
        }
        AOoffset2=calloc(2*(AO->nao),sizeof(double));
        AOout=calloc(nAOtotal,sizeof(double));
        AOdv=calloc(nAOtotal*nAOcolst,sizeof(double));
        AOda=calloc(nAOtotal*nAOcols,sizeof(double));
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
    //Occultation data
    int nOCtotal=0;
    int nOCcols=0;
    int nOC=0;
    int nOCoffsets=0;
    int nOCcolst=0;
    int nChordoffsets=0;
    int nvectorreg=0;
    double *OCoffset,*OCoffset2,*OCout,*OCdv,*OCda,OCfit=0,*OCdoff;
    double *Chordoffset,*Chordoffset2;
    double *dChordoffset; //Derivative matrix for chord offsets
    double vectorreg; //regularization for free chord offsets
    double *dvectorreg; //derivative matrix for free chord offsets, 1xnChordoffsets 
    if(INI_HAVE_OC)
    {
        nOCtotal=4*(OC->ntotal);
        nOCcolst=3*nvert+3;
        nOCcols=alength+3;
        nOCoffsets=2*(OC->noc);
        nOC=OC->noc;
        OCoffset=calloc(nOCoffsets,sizeof(double));
        
        if(INI_OC_OFFSET!=NULL)
            memcpy(OCoffset,INI_OC_OFFSET,sizeof(double)*nOCoffsets);
            
        OCoffset2=calloc(nOCoffsets,sizeof(double));
        OCdv=calloc(nOCtotal*nOCcolst,sizeof(double));
        OCda=calloc(nOCtotal*(alength+3),sizeof(double));
        OCout=calloc(nOCtotal,sizeof(double));
        OCdoff=calloc(nOCtotal*nOCoffsets,sizeof(double));
        Chordoffset=calloc(OC->ntotal,sizeof(double));
        Chordoffset2=calloc(OC->ntotal,sizeof(double));
        if(INI_CHORD_OFFSET!=NULL)
        {
            memcpy(Chordoffset,INI_CHORD_OFFSET,sizeof(double)*OC->ntotal);
            memcpy(Chordoffset2,INI_CHORD_OFFSET,sizeof(double)*OC->ntotal);
        }
        if(INI_FREE_CHORD_NMR>0)
        {
            nChordoffsets=OC->ntotal;
           nvectorreg=1;
           dChordoffset=calloc(nOCtotal*nChordoffsets,sizeof(double));
           dvectorreg=calloc(nChordoffsets,sizeof(double));
        }
    }
    int nRDtotal=0;
    int nRDcols=0,nRDcolst=0;
    int nRD=0;
    int nRDoffsets=0;
    int nRDscale=0;
    int nRDexp=0;
    double RDexp=INI_RDEXP;
    double RDexp2=RDexp;
    double RDfit=0;
    double *RDda;
    double *RDoffset,*RDoffset2,*RDscale,*RDscale2,*RDout,*RDdv,*RDdoff,*RDdscale,*RDdexp;
    if(INI_HAVE_RD)
    {
        nRDtotal=2*(RD->ntotal);
        nRDcolst=3*nvert+3;
        nRDcols=alength+3;
        nRDoffsets=2*(RD->nRD);
        nRD=RD->nRD;
        nRDscale=nRD;
        nRDexp=0;
        RDoffset=calloc(nRDoffsets,sizeof(double));
        RDoffset2=calloc(nRDoffsets,sizeof(double));
        RDscale=calloc(nRDscale,sizeof(double));
        RDscale2=calloc(nRDscale,sizeof(double));
        RDout=calloc(nRDtotal,sizeof(double));
        RDdv=calloc(nRDtotal*nRDcolst,sizeof(double));
        RDda=calloc(nRDtotal*nRDcols,sizeof(double));
        RDdoff=calloc(nRDtotal*nRDoffsets,sizeof(double));
        RDdscale=calloc(nRDtotal*nRD,sizeof(double));
        RDdexp=calloc(nRDtotal,sizeof(double));
        
    }  
    //Variables for regularization terms
    double CRres; //Convex reg
    double *dCRdv,*dCRda;
    double *y,*y2;
    
   
    dCRdv=calloc(3*nvert+3,sizeof(double));
    dCRda=calloc(alength+3,sizeof(double));
    //Octantoid regularization
    double Ores;
    double *dOda;
    dOda=calloc(alength,sizeof(double));
    
    
    double *S,*S2;
    double *J,*JTJ,*JTJpd;
    
    ////////////////////////////////////
//     write_matrix_file("/tmp/data0r.txt",AO->datar[0],1,AO->nobs[0]);
//     write_matrix_file("/tmp/data0i.txt",AO->datai[0],1,AO->nobs[0]);
//     write_matrix_file("/tmp/data1r.txt",AO->datar[1],1,AO->nobs[1]);
//     write_matrix_file("/tmp/data1i.txt",AO->datai[1],1,AO->nobs[1]);
//     write_matrix_file("/tmp/freqx.txt",AO->freqx[0],1,AO->nobs[0]);
//     write_matrix_file("/tmp/freqy.txt",AO->freqy[0],1,AO->nobs[0]);
    /////////////////////////////////////////////////////
    int Slength;
   
    //Number of rows in J
    Slength=nLCtotal+nAOtotal+nOCtotal+nRDtotal+1+1+1+nvectorreg+nAlbreg; //LC points+AO points+OC points+RD points convex reg+oct reg+dihedral+[vector reg for free chords]
    int nJcols=alength+3+ncalib+nAlbedo+nAOoffsets+nAOscale+nOCoffsets+nChordoffsets+nRDoffsets+nRDscale+nRDexp;
    S=calloc(Slength,sizeof(double)); 
    J=calloc(Slength*(nJcols),sizeof(double));
    int padding=3+nAOoffsets; 
    JTJ=calloc((nJcols)*(nJcols),sizeof(double));
    JTJpd=calloc((nJcols)*(nJcols),sizeof(double));
    int OCoffsetcolpos=alength+3+nAOoffsets+nAOscale;
    int Chordoffsetcolpos=alength+3+nAOoffsets+nAOscale+nOCoffsets;
    int OCrowpos=nLCtotal+nAOtotal;
    int RDoffsetcolpos=alength+3+nAOoffsets+nAOscale+nOCoffsets+nChordoffsets;
    int RDscalecolpos=alength+3+nAOoffsets+nAOscale+nOCoffsets+nChordoffsets+nRDoffsets;
    int RDrowpos=nLCtotal+nAOtotal+nOCtotal;
    int Albcolpos=alength+3+nAOoffsets+nAOscale+nOCoffsets+nChordoffsets+nRDoffsets+nRDscale+nRDexp;
    int regpos=nLCtotal+nAOtotal+nOCtotal+nRDtotal;
     int Albregpos=regpos+1+1+1+nvectorreg;
    double *LCout,*dLCdv,*rhs,*X,*dLCda,*dLCdp=NULL;
    LCout=calloc(nLCtotal,sizeof(double));
    dLCdv=calloc((nLCtotal)*(3*nvert+3),sizeof(double));
    dLCda=calloc((nLCtotal)*(alength+3),sizeof(double));
    rhs=calloc(nJcols,sizeof(double));
    X=calloc(nJcols,sizeof(double));
   //Check if some parameters are fixed:
    double *Mask_Matrix;
    int *Mask;
    int nMask=0;
    double *MJTJ;
    double *MJTJpd;
    double *Mrhs;
    double *MX;
    if(INI_MASK_SET==1)
    {
        Mask=calloc(nJcols,sizeof(int));
        if(INI_FIX_SHAPE==1)
            for(int k=0;k<alength;k++)
                Mask[k]=1;
        if(INI_FIX_ANGLES==1)
            for(int k=alength;k<alength+3;k++)
                Mask[k]=1;
        if(INI_FREE_CHORD_NMR>0)
        {
            for(int k=0;k<nChordoffsets;k++)
                Mask[Chordoffsetcolpos+k]=1;
            for(int k=0;k<INI_FREE_CHORD_NMR;k++)
                Mask[Chordoffsetcolpos+INI_FREE_CHORD_LIST[k]-1]=0;
        }
       
        if(INI_PHASE_MASK!=NULL)
        {
            Mask[nJcols-4]=INI_PHASE_MASK[0];
            Mask[nJcols-3]=INI_PHASE_MASK[1];
            Mask[nJcols-2]=INI_PHASE_MASK[2];
            Mask[nJcols-1]=INI_PHASE_MASK[3];
            
        }
        mask_matrix(nJcols,Mask,&Mask_Matrix,&nMask);
        MJTJ=calloc(nMask*nMask,sizeof(double));
        MJTJpd=calloc(nMask*nMask,sizeof(double));
        Mrhs=calloc(nMask,sizeof(double));
        MX=calloc(nMask,sizeof(double));
        
    }
    //Weights
    int nvertn=0;
    double angW=INI_DW;
    double AOW=INI_AOW;
    double lcW=INI_LCW;
    double cW=INI_CW;
    double oW=INI_OW;
    double ocW=INI_OCW;
    double RDW=INI_RW;
    double lambda=INI_LAMBDA;
    int decreased=1;
    double chisq=1e9;
    double chisq2;
    double LCfit=0;
     double dec=10;
    double prevdec=dec;
    double threshold=INI_MINDEC;
    int count=0;
    int DONE=0;

   //Allocate for D matrix
    D1=calloc((3*nvert+3)*(alength+3),sizeof(double));
    D2=calloc((3*nvert+3+nAOoffsets)*(alength+3+nAOoffsets),sizeof(double));
  if(USE_CALIB==1)
        dLCdp=calloc(nLCtotal*4,sizeof(double));
  //Restore previous state if the state file is available
  if(INI_RESTORE_STATE)
  {
      //Restore angles
      read_state_file(INI_RESTORE_STATE,"#Angles",angles,4);
      angles[0]=(90-angles[0])*PI/180;
      angles[1]=angles[1]*PI/180;
      angles[2]=24*2*PI*1.0/angles[2];
      angles[3]=angles[3]*PI/180;
      //Restore shape
      read_state_file(INI_RESTORE_STATE,"#Params",a,alength);
      //Restore AO
      read_state_file(INI_RESTORE_STATE,"#AOoffset",AOoffset,nAOoffsets);
      read_state_file(INI_RESTORE_STATE,"#AOscale",AOscale,nAOscale);
      //Restore OC
      read_state_file(INI_RESTORE_STATE,"#OCCoffset",OCoffset,nOCoffsets);
      read_state_file(INI_RESTORE_STATE,"#Chordoffset",Chordoffset,nChordoffsets);
      //Restore albedo
      read_state_file(INI_RESTORE_STATE,"#Albedolog",eAlbedo,nfac);
  }
      
  //This is the main optimization loop
    for(int k=0;k<NUM_OF_ROUNDS;k++)
    {
        if(decreased==1)
        {
        octantoid_to_trimesh(a,INI_LMAX,INI_NROWS,tlist,vlist,D1,3); //D is 3nvert+3 x alength matrix
        octantoid_to_trimesh(a,INI_LMAX,INI_NROWS,tlist,vlist,D2,3+nAOoffsets); //TBD: Fix this ugly thing
        //We must pad D before using it to preserve angle and offset derivatives
        //Calculate LCs
//          write_shape_file("/tmp/Inishape.txt",tlist,vlist,nfac,nvert);
//          write_matrix_file("/tmp/a.txt",a,1,alength);
//           write_matrix_file("/tmp/D1.txt",D1,3*nvert+3,alength+3);
//           write_matrix_file("/tmp/D2.txt",D2,3*nvert+3+nAOoffsets,alength+3+nAOoffsets);
        
        calculate_lcs(tlist,vlist,nfac,nvert,angles,LC,NULL,nvert,nvert,LCout,dLCdv,eAlbedo,Alblimits,dLCdalb,params,dLCdp,1);
//          write_matrix_file("/tmp/LCout.txt",LCout,1,nLCtotal);
//          write_matrix_file("/tmp/dLCdv.txt",dLCdv,nLCtotal,3*nvert+3);
       //Convert to derivative matrix wrt a
       
         if(INI_FIT_ALBEDO==1)
            {
                mult_with_cons(dLCdalb,nLCtotal,nfac,lcW);
                localsmooth(tlist,vlist,nfac,nvert,eAlbedo,Alblimits,Albreg,dAlbreg);
                mult_with_cons(dAlbreg,nfac,nfac,-INI_ALBREGW);
                mult_with_cons(Albreg,1,nfac,INI_ALBREGW);
               
                set_submatrix(S,1,Slength,Albreg,1,nfac,0,Albregpos);
               
                set_submatrix(J,Slength,nJcols,dLCdalb,nLCtotal,nfac,0,Albcolpos);
               
                set_submatrix(J,Slength,nJcols,dAlbreg,nfac,nfac,Albregpos,Albcolpos);
            }
           
        matrix_prod(dLCdv,nLCtotal,3*nvert+3,D1,alength+3,dLCda); //dLCda is nLCtotal x alength+3 matrix
//         write_matrix_file("/tmp/dLCda.txt",dLCda,nLCtotal,alength+3);
      //  printf("Angles: %f %f %f %f\n",angles[0],angles[1],angles[2],angles[3]);
      //AO data
       // printf("before AOs\n");
        if(INI_HAVE_AO)
        {
        Calculate_AOs(tlist,vlist,nfac,nvert,angles,AO,AOoffset,NULL,nvert,nvert,NULL,AOscale,AOout,AOdv,AOds,1);
        //Convert to derivative wrt params a
        matrix_prod(AOdv,nAOtotal,nAOcolst,D2,nAOcols,AOda); //AOda is nAOtotal x (nAOcols)
//         write_matrix_file("/tmp/AOda.txt",AOda,nAOtotal,nAOcols);
//          write_matrix_file("/tmp/AOdv.txt",AOdv,nAOtotal,nAOcolst);
//         write_matrix_file("/tmp/AOout.txt",AOout,1,nAOtotal);
        mult_with_cons(AOout,1,nAOtotal,AOW);
        mult_with_cons(AOda,nAOtotal,nAOcols,AOW);
        set_submatrix(S,1,Slength,AOout,1,nAOtotal,0,nLCtotal);
        matrix_transprod(AOout,nAOtotal,1,&AOfit);
         set_submatrix(J,Slength,nJcols,AOda,nAOtotal,nAOcols,nLCtotal,0);
         if(INI_AO_SCALING)
        {
            mult_with_cons(AOds,nAOtotal,nAO,AOW);
            set_submatrix(J,Slength,nJcols,AOds,nAOtotal,nAO,nLCtotal,nAOcols);
        }
        }
        if(INI_HAVE_OC)
        {
            if(INI_FREE_CHORD_NMR>0)
                {
                    calculate_OCs(tlist,vlist,nfac,nvert,angles,OC,OCoffset,INI_OC_WEIGHT,NULL,nvert,nvert,Chordoffset,OCout,OCdv,OCdoff,dChordoffset);
                    vector_regularization(Chordoffset,nChordoffsets,&vectorreg,dvectorreg);
                    mult_with_cons(dChordoffset,nOCtotal,nChordoffsets,-ocW);
                    mult_with_cons(dvectorreg,1,nChordoffsets,-INI_CHRDW);
                    vectorreg*=INI_CHRDW;
                    S[regpos+3]=vectorreg; //NOTE ABSOLUTE ADDRESS HERE. FIX!
                    set_submatrix(J,Slength,nJcols,dChordoffset,nOCtotal,nChordoffsets,OCrowpos,Chordoffsetcolpos); //Chord offsets
                    set_submatrix(J,Slength,nJcols,dvectorreg,1,nChordoffsets,regpos+3,Chordoffsetcolpos);
                  
                }
                    
            else
                calculate_OCs(tlist,vlist,nfac,nvert,angles,OC,OCoffset,INI_OC_WEIGHT,NULL,nvert,nvert,Chordoffset,OCout,OCdv,OCdoff,NULL);
            matrix_prod(OCdv,nOCtotal,nOCcolst,D1,alength+3,OCda);
            mult_with_cons(OCout,1,nOCtotal,ocW);
            mult_with_cons(OCda,nOCtotal,nOCcols,-ocW); //Note the minus here.
            mult_with_cons(OCdoff,nOCtotal,nOCoffsets,-ocW);
            set_submatrix(S,1,Slength,OCout,1,nOCtotal,0,OCrowpos);
            set_submatrix(J,Slength,nJcols,OCda,nOCtotal,nOCcols,OCrowpos,0);
            set_submatrix(J,Slength,nJcols,OCdoff,nOCtotal,nOCoffsets,OCrowpos,OCoffsetcolpos);
            matrix_transprod(OCout,nOCtotal,1,&OCfit);
        } 
        if(INI_HAVE_RD)
            {
                Calculate_RDs(tlist,vlist,nfac,nvert,angles,RD,RDoffset,NULL,nvert,nvert,INI_RD_WEIGHT,RDscale,RDexp,RDout,RDdv,RDdoff,RDdscale,RDdexp,1);
                matrix_prod(RDdv,nRDtotal,nRDcolst,D1,nRDcols,RDda);
                
                mult_with_cons(RDout,1,nRDtotal,RDW);
                
                mult_with_cons(RDda,nRDtotal,nRDcols,RDW); //NO minus here
                mult_with_cons(RDdoff,nRDtotal,nRDoffsets,RDW);
                
                mult_with_cons(RDdscale,nRDtotal,nRDscale,RDW);
               
                if(nRDexp>0)
                {
                    mult_with_cons(RDdexp,nRDtotal,1,RDW);
                    set_submatrix(J,Slength,nJcols,RDdexp,nRDtotal,1,RDrowpos,RDscalecolpos+nRDscale);
                }
                
                set_submatrix(S,1,Slength,RDout,1,nRDtotal,0,RDrowpos);
              
                set_submatrix(J,Slength,nJcols,RDda,nRDtotal,nRDcols,RDrowpos,0);
                
                set_submatrix(J,Slength,nJcols,RDdoff,nRDtotal,nRDoffsets,RDrowpos,RDoffsetcolpos); //[RDda 0 RDdoff RDdscale RDdexp]
               
                set_submatrix(J,Slength,nJcols,RDdscale,nRDtotal,nRDscale,RDrowpos,RDscalecolpos);
                
                matrix_transprod(RDout,nRDtotal,1,&RDfit);
            }
        convex_reg(tlist,vlist,nfac,nvert,NULL,nvertn,nvert,&CRres,dCRdv);
        dihedral_angle_reg(tlist,vlist,nfac,nvert,NULL,nvert,nvert,&ANGres,dANGdv);
        matrix_prod(dANGdv,1,3*nvert+3,D1,alength+3,dANGda);
        //Convert to derivative wrt a
        matrix_prod(dCRdv,1,3*nvert+3,D1,alength+3,dCRda);
        //Octantoid regularization
        octantoid_reg(a,INI_LMAX,&Ores,dOda);

        
        mult_with_cons(LCout,1,nLCtotal,lcW);
        mult_with_cons(dLCda,nLCtotal,alength+3,lcW);
        CRres*=cW;
        mult_with_cons(dANGda,1,alength+3,-angW);
        ANGres*=angW;
        mult_with_cons(dCRda,1,alength+3,-cW);
       
        Ores*=oW;
        mult_with_cons(dOda,1,alength,-oW);
        
        
        //Build the res vector and matrix
        set_submatrix(S,1,Slength,LCout,1,nLCtotal,0,0);
        
       
        S[regpos]=CRres; //Convex reg terms
        S[regpos+1]=Ores; //Octantoid reg
       S[regpos+2]=ANGres;
        
       
        matrix_transprod(S,Slength,1,&chisq);
       matrix_transprod(LCout,nLCtotal,1,&LCfit);
        printf("Round: %d chisq: %4.2f LCfit: %4.2f AOfit: %4.2f OCfit: %4.2f RDfit: %4.2f Convex reg: %4.2f Octantoid reg: %4.2f Dihedral reg: %4.2f\n",k,chisq,LCfit,AOfit,OCfit,RDfit,pow(CRres,2),pow(Ores,2),pow(ANGres,2)); 
       //printf("scale: %f %f\n",AOscale[0],AOscale[1]);
        //Construct Jacobian matrix
        
        set_submatrix(J,Slength,nJcols,dLCda,nLCtotal,alength+3,0,0);
       
        set_submatrix(J,Slength,nJcols,dCRda,1,alength+3,regpos,0);
        set_submatrix(J,Slength,nJcols,dOda,1,alength,regpos+1,0);
        set_submatrix(J,Slength,nJcols,dANGda,1,alength,regpos+2,0);
        
           if(USE_CALIB==1)
            {
                mult_with_cons(dLCdp,nLCtotal,4,lcW);
                set_submatrix(J,Slength,nJcols,dLCdp,nLCtotal,4,0,nJcols-4);
            } 
       
         // write_matrix_file("/tmp/S.txt",S,1,Slength);
         // write_matrix_file("/tmp/J.txt",J,Slength,nJcols);
        
        //Calculate J^TJ+lambda*J
       // matrix_transprod(J,Slength,nJcols,JTJ);
         
        
        
        matrix_transprod(J,Slength,nJcols,JTJ);
         matrix_vectorprod(J,Slength,nJcols,S,rhs,1); //rhs=J^T*S;
         //////////////////

//          write_matrix_file("/tmp/a.txt",a,alength,1);
//          write_matrix_file("/tmp/dOda.txt",dOda,1,alength);
//          write_matrix_file("/tmp/dCRda.txt",dCRda,1,alength);
         if(INI_MASK_SET==1)
            {
                matrix_prod_ATBA(Mask_Matrix,nJcols,nMask,JTJ,MJTJ);
                matrix_prod_ATB(Mask_Matrix,nJcols,nMask,rhs,1,Mrhs);
            }
        }
        
        if(INI_MASK_SET==1)
        {
            matrix_adddiag(MJTJ,MJTJpd,nMask,lambda);
            solve_matrix_eq(MJTJpd,nMask,Mrhs,MX);
            matrix_prod(Mask_Matrix,nJcols,nMask,MX,1,X);
        }
        else
        {
        matrix_adddiag(JTJ,JTJpd,nJcols,lambda); //JTJpd=JTJ+lambda*diag(JTJ)
        
        solve_matrix_eq(JTJpd,nJcols,rhs,X); //Solve the LM 
        }
      
    
       // add_vector_to_vlist(vlist,X,vlist2,nvert);
       matrix_plus2(a,1,alength,X,a2);
       if(INI_HAVE_AO)
       {
        matrix_plus2(AOoffset,1,nAOoffsets,&X[alength+3],AOoffset2);
       if(INI_AO_SCALING)
           matrix_plus2(AOscale,1,nAO,&X[alength+3+nAOoffsets],AOscale2);
       }
       if(INI_HAVE_OC)
       {
           matrix_plus2(OCoffset,1,nOCoffsets,&X[OCoffsetcolpos],OCoffset2);
           if(INI_FREE_CHORD_NMR>0)
                matrix_plus2(Chordoffset,1,nChordoffsets,&X[Chordoffsetcolpos],Chordoffset2);
       }
       if(INI_FIT_ALBEDO==1)
        {
            for(int j=0;j<nfac;j++)
                eAlbedo2[j]=eAlbedo[j]+X[Albcolpos+j];
        }
        angles2[0]=angles[0]+X[alength];
        angles2[1]=angles[1]+X[alength+1];
        angles2[2]=angles[2]+X[alength+2];
        angles2[3]=angles[3];
       
      if(USE_CALIB==1)
        {
            params2[0]=params[0]+X[nJcols-4];
            params2[1]=params[1]+X[nJcols-3];
            params2[2]=params[2]+X[nJcols-2];
            params2[3]=params[3]+X[nJcols-1];
        }
         if(INI_HAVE_RD)
        {
            matrix_plus2(RDoffset,1,nRDoffsets,&X[RDoffsetcolpos],RDoffset2);
            matrix_plus2(RDscale,1,nRDscale,&X[RDscalecolpos],RDscale2);
            if(nRDexp>0)
                RDexp2=X[RDscalecolpos+nRDscale]+RDexp;
        }
        //New shape is a2. Check for fit
       
        octantoid_to_trimesh(a2,INI_LMAX,INI_NROWS,tlist,vlist2,D1,3);
       
        
        calculate_lcs(tlist,vlist2,nfac,nvert,angles2,LC,NULL,nvertn,nvert,LCout,dLCdv,eAlbedo2,Alblimits,dLCdalb,params2,dLCdp,0);
         if(INI_FIT_ALBEDO==1)
        {
                localsmooth(tlist,vlist2,nfac,nvert,eAlbedo2,Alblimits,Albreg,dAlbreg);
                mult_with_cons(Albreg,1,nfac,INI_ALBREGW);
                set_submatrix(S,1,Slength,Albreg,1,nfac,0,Albregpos);
        }
        if(INI_HAVE_AO)
        {
            Calculate_AOs(tlist,vlist2,nfac,nvert,angles2,AO,AOoffset2,NULL,nvertn,nvert,NULL,AOscale2,AOout,AOdv,NULL,0);
            mult_with_cons(AOout,1,nAOtotal,AOW);
            set_submatrix(S,1,Slength,AOout,1,nAOtotal,0,nLCtotal);
        }
        if(INI_HAVE_OC)
        {
            calculate_OCs(tlist,vlist2,nfac,nvert,angles2,OC,OCoffset2,INI_OC_WEIGHT,NULL,nvert,nvert,Chordoffset2,OCout,OCdv,OCdoff,NULL);
            mult_with_cons(OCout,1,nOCtotal,ocW);
            set_submatrix(S,1,Slength,OCout,1,nOCtotal,0,OCrowpos);
            if(INI_FREE_CHORD_NMR>0)
            {
                vector_regularization(Chordoffset2,nChordoffsets,&vectorreg,dvectorreg);
                S[regpos+3]=INI_CHRDW*vectorreg;
            }
        }
        if(INI_HAVE_RD)
        {
          Calculate_RDs(tlist,vlist2,nfac,nvert,angles2,RD,RDoffset2,NULL,nvert,nvert,INI_RD_WEIGHT,RDscale2,RDexp2,RDout,RDdv,RDdoff,RDdscale,RDdexp,0);  
          mult_with_cons(RDout,1,nRDtotal,RDW);
          set_submatrix(S,1,Slength,RDout,1,nRDtotal,0,RDrowpos);
        }
        
        dihedral_angle_reg(tlist,vlist2,nfac,nvert,NULL,nvert,nvert,&ANGres,dANGdv);
        convex_reg(tlist,vlist2,nfac,nvert,NULL,nvertn,nvert,&CRres,dCRdv);
        octantoid_reg(a2,INI_LMAX,&Ores,dOda);
        mult_with_cons(LCout,1,nLCtotal,lcW);
        
        CRres*=cW;
        Ores*=oW;
        ANGres*=angW;
        
       set_submatrix(S,1,Slength,LCout,1,nLCtotal,0,0);
        
        S[regpos]=CRres; //Convex reg terms
        S[regpos+1]=Ores;
        S[regpos+2]=ANGres;
        
       
        matrix_transprod(S,Slength,1,&chisq2);
         matrix_transprod(AOout,nAOtotal,1,&AOfit);
        printf("Round: %d chisq2: %7.2f\n",k+1,chisq2);
      
       if(chisq2<chisq)
       {
           prevdec=dec;
            dec=chisq-chisq2;
            count++;
           printf("decreased!\n");
           memcpy(a,a2,sizeof(double)*(alength));
           if(INI_HAVE_AO)
           {
           memcpy(AOoffset,AOoffset2,sizeof(double)*(nAOoffsets));
          if(INI_AO_SCALING)
              memcpy(AOscale,AOscale2,sizeof(double)*nAO);
           }
           if(INI_HAVE_OC)
           {
               memcpy(OCoffset,OCoffset2,sizeof(double)*(nOCoffsets));
               if(INI_FREE_CHORD_NMR>0)
                    memcpy(Chordoffset,Chordoffset2,sizeof(double)*(nChordoffsets));
           }
           if(USE_CALIB==1)
                memcpy(params,params2,4*sizeof(double));
           if(INI_HAVE_RD)
            {
                memcpy(RDoffset,RDoffset2,sizeof(double)*(nRDoffsets));
                memcpy(RDscale,RDscale2,sizeof(double)*nRDscale);
                RDexp=RDexp2;
            }
           decreased=1;
           angles[0]=angles2[0];
           angles[1]=angles2[1];
           angles[2]=angles2[2];
            if(INI_FIT_ALBEDO==1)
                memcpy(eAlbedo,eAlbedo2,sizeof(double)*nfac);
           lambda=lambda/INI_LAMBDADEC;
           zero_array(J,Slength*nJcols);
           zero_array(S,Slength);
           write_shape_file(OUT_SHAPE_FILE,tlist,vlist2,nfac,nvert);
       }
       else
       {
           lambda=INI_LAMBDAINC*lambda;
           decreased=0;
          
    }
    if(prevdec<threshold && dec<threshold)
        {
            printf("Two previous optimization steps below threshold, stopping optimization loop.\n");
            DONE=1;
        }
        if(lambda>INI_LAMBDAMAX)
        {
            printf("lambda is larger than LambdaMax, halting optimization\n");
            DONE=1;
        }
    if(k==NUM_OF_ROUNDS-1 || DONE==1)
    {
        
        
        printf("Angles: %.4f %.4f %.8f\n",angles[0],angles[1],angles[2]);
        if(INI_HAVE_OC)
        {
            printf("Occultation offsets:\n");
            print_matrix(OCoffset,1,nOCoffsets);
            if(INI_FREE_CHORD_NMR>0)
                 {
                     printf("OCC chord offsets(seconds):\n");
                     print_matrix(Chordoffset,1,nChordoffsets);
                 }
        }
        
      //  printf("Offsets: %f %f %f %f\n",offset[0],offset[1],offset[2],offset[3]);
      printf("Writing shape information to %s\n",OUT_SHAPE_FILE);
        write_shape_file(OUT_SHAPE_FILE,tlist,vlist,nfac,nvert);
        if( OUT_OBJSHAPE_FILE!=NULL)
            {
                printf("Writing obj shape to: %s\n",OUT_OBJSHAPE_FILE);
                write_obj_file(OUT_OBJSHAPE_FILE,tlist,vlist,nfac,nvert);
            }
        FILE* fid=fopen(OUT_PARAM_FILE,"w");
         FILE* fidlc=fopen(OUT_LC_FILE,"w");
         if(fidlc!=NULL)
            {
                for(int j=0;j<LC->nlc;j++)
                    zero_array(LC->lcs[j],LC->nobs[j]);
                 octantoid_to_trimesh(a,INI_LMAX,INI_NROWS,tlist,vlist,D1,3);
                calculate_lcs(tlist,vlist,nfac,nvert,angles,LC,NULL,nvert,nvert,LCout,dLCdv,NULL,NULL,NULL,NULL,NULL,0);
                for(int j=0;j<LC->ntotal;j++)
                    fprintf(fidlc,"%f\n",-LCout[j]);
                fclose(fidlc);
            }
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
        
            if(params!=NULL)
            {
            printf("Phase parameters: \n");
                print_matrix(params,1,4);
            }
            
        if(OUT_SHAPE_PARAM_FILE!=NULL)
            write_matrix_file(OUT_SHAPE_PARAM_FILE,a,1,alength);
        if(INI_WRITE_STATE_FILE!=NULL)
            {
                printf("Writing state file %s:\n",INI_WRITE_STATE_FILE);
                FILE *fp;
                fp=fopen(INI_WRITE_STATE_FILE,"w");
                if(fp==NULL)
                {
                    fprintf(stderr,"Cannot open file %s for writing the final state file\n",INI_WRITE_STATE_FILE);
                    exit(-1);
                }
                //write angles
                fprintf(fp,"#Angles:\n");
                fprintf(fp,"%.2f %.2f %.8f %.2f\n",90-angles[0]*180/PI,angles[1]*180/PI,24*2*PI/angles[2],angles[3]*180/PI);
                if(INI_PHASE_PARAMS!=NULL)
                {
                    fprintf(fp,"#Phaseparams:\n");
                    fprintf(fp,"%.4f %.4f %.4f %.4f\n",params[0],params[1],params[2],params[3]);
                }
                if(INI_FIT_ALBEDO==1)
                {
                    fprintf(fp,"#Albedo: %d\n",nfac);
                    for(int j=0;j<nfac;j++)
                        fprintf(fp,"%.2f ",Alblimits[0]+(Alblimits[1]-Alblimits[0])*exp(eAlbedo[j])/(exp(eAlbedo[j])+1.0));
                  fprintf(fp,"\n");
                  fprintf(fp,"#Albedolog: %d\n",nfac);
                    for(int j=0;j<nfac;j++)
                        fprintf(fp,"%.2f ",eAlbedo[j]);
                  fprintf(fp,"\n");
                }  
                //Shape parameters
                fprintf(fp,"#Shape:\n");
                 fprintf(fp,"#LMAX: %d\n",INI_LMAX);
                 fprintf(fp,"#Params: %d\n",alength);
                for(int j=0;j<alength;j++)
                    fprintf(fp,"%.4f ",a[j]);
                fprintf(fp,"\n");
               
                fprintf(fp,"#Polyhedron: %d %d\n",nfac,nvert);
                for(int j=0;j<3*nfac;j++)
                    fprintf(fp,"%d ",tlist[j]);
                fprintf(fp,"\n");
                for(int j=0;j<3*nvert;j++)
                    fprintf(fp,"%.4f ",vlist[j]);
                 fprintf(fp,"\n");
                if(INI_HAVE_AO)
                {
                fprintf(fp,"#AOoffset: %d\n",nAOoffsets);
                for(int j=0;j<nAOoffsets;j++)
                    fprintf(fp,"%.4f ",AOoffset[j]);
                fprintf(fp,"\n");
                fprintf(fp,"#AOscale: %d\n",nAOscale);
                for(int j=0;j<nAOscale;j++)
                    fprintf(fp,"%.4f ",AOscale[j]);
                fprintf(fp,"\n");
                }
                if(INI_HAVE_OC)
                {
                    fprintf(fp,"#OCCoffset: %d\n",nOCoffsets);
                    for(int j=0;j<nOCoffsets;j++)
                        fprintf(fp,"%.4f ",OCoffset[j]);
                    fprintf(fp,"\n");
                    fprintf(fp,"#Chordoffset: %d\n",nChordoffsets);
                    for(int j=0;j<nChordoffsets;j++)
                        fprintf(fp,"%.4f ",Chordoffset[j]);
                    fprintf(fp,"\n");
                }
                if(INI_HAVE_RD)
                {
                    fprintf(fp,"#RDoffset: %d\n",nRDoffsets);
                    for(int j=0;j<nRDoffsets;j++)
                        fprintf(fp,"%.4f ",RDoffset[j]);
                    fprintf(fp,"\n");
                    fprintf(fp,"#RDscale: %d\n",nRDscale);
                    for(int j=0;j<nRDscale;j++)
                        fprintf(fp,"%.4f ",RDscale[j]);
                }
                
               fclose(fp); 
            }
            if(INI_FIT_ALBEDO==1 && INI_ALBEDO_OUT_FILE!=NULL)
            {
                double *falbedo=calloc(nfac,sizeof(double));
                for(int j=0;j<nfac;j++)
                    falbedo[j]=Alblimits[0]+(Alblimits[1]-Alblimits[0])*exp(eAlbedo[j])/(exp(eAlbedo[j])+1.0);
                write_matrix_file(INI_ALBEDO_OUT_FILE,falbedo,1,nfac);
                free(falbedo);
            }
        return;
    }
    
//     zero_array(dLCdv,nLCtotal*(3*nvert+3));
//     zero_array(AOdv,nAOcols*nAOtotal);
//     zero_array(dCRdv,3*nvert+3);
//     zero_array(dANGdv,3*nvert+3);
//     zero_array(dAdv,nfacn*(3*nvert+3));
    
}

}
