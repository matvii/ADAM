#include<stdio.h>
#include<stdlib.h>
#include"utils.h"
#include"matrix_ops.h"
#include"structs.h"
#include"globals.h"
void Calculate_AOs_fit(int *tlist,double *vlist,int nfac,int nvert,double *angles,AOstruct *AOs,double *offset,double *D,int dm,int dn,double *Weight,double *scale,double *FT,double *FTdv,double* FTdS,double *Albedo,double *Alimit,double *dA,int deriv);
void check_model_fit(LCstruct *LC,AOstruct *AO,OCstruct *OC,RDstruct *RD,CNTRstruct *CR)
{
    //First initialize the initial shape
    int *tlist,*tlistn,nfac,nvert,nfacn,nvertn;
    double *vlist,*vlist2,*vlistn,*D;
    int nLCtotal=LC->ntotal;
    INI_SD_LEVEL=-1;
    if(INI_INPUTOBJSHAPE!=NULL)
        read_obj_file(INI_INPUTOBJSHAPE,&tlist,&vlist,&nfac,&nvert);
    else if(INI_INPUTSHAPE!=NULL)
        read_shape(INI_INPUTSHAPE,&tlist,&vlist,&nfac,&nvert,0);
        
    else
        read_shape(OUT_SHAPE_FILE,&tlist,&vlist,&nfac,&nvert,0);
    tlistn=tlist;
    vlistn=vlist;
    nfacn=nfac;
    nvertn=nvert;
    //Setup angles
    double angles[]={(90-INI_ANGLE_B)*PI/180,INI_ANGLE_L*PI/180,24*2*PI*1.0/INI_ANGLE_P,INI_ANGLE_PHI0*PI/180};
    double angles2[4]={0,0,0,0};
    angles2[3]=angles[3];
    
   INI_VERBOSE=1;
    //This is for calibrated lightcurves
    double *params=NULL,*params2=NULL;
     int ncalib=0;
    if(INI_PHASE_PARAMS!=NULL)
    {
        params=calloc(3,sizeof(double));
        params2=calloc(3,sizeof(double));
        params[0]=INI_PHASE_PARAMS[0];
        params[1]=INI_PHASE_PARAMS[1];
        params[2]=INI_PHASE_PARAMS[2];
       // params[3]=INI_PHASE_PARAMS[3];
        memcpy(params2,params,sizeof(double)*3);
        ncalib=3;
    }
    int USE_CALIB=LC->calib;
    
   
       
        
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
        nAlbedo=nfacn;
        eAlbedo=calloc(nfacn,sizeof(double));
        eAlbedo2=calloc(nfacn,sizeof(double));
        dLCdalb=calloc(nLCtotal*nfacn,sizeof(double));
        Alblimits=calloc(2,sizeof(double));
        Alblimits[0]=INI_ALBEDO_MIN;
        Alblimits[1]=INI_ALBEDO_MAX;
        Albreg=calloc(nfacn,sizeof(double));
        dAlbreg=calloc(nfacn*nfacn,sizeof(double));
        nAlbreg=nfacn;
    }
        
        
    int nAOtotal=0;
    int nAOcols=0;
    int nAO=0;
    int nAOoffsets=0;
    int nAOscale=0;
    double *AOoffset,*AOoffset2,*AOout,*AOdv,AOfit=0,*AOscale,*AOscale2,*AOds;
    if(INI_HAVE_AO)
    {
        nAOtotal=2*(AO->ntotal);
        nAOcols=3*nvert+3+2*(AO->nao);
        nAOoffsets=2*(AO->nao);
        nAO=AO->nao;
        AOoffset=calloc(2*(AO->nao),sizeof(double));
        if(INI_INPUT_AO_OFFSET!=NULL)
        {
            read_vector_file(INI_INPUT_AO_OFFSET,AOoffset,2*(AO->nao));
        }
        AOoffset2=calloc(2*(AO->nao),sizeof(double));
        AOout=calloc(nAOtotal,sizeof(double));
        AOdv=calloc(nAOtotal*nAOcols,sizeof(double));
        AOfit=0;
        AOscale=NULL;
        AOscale2=NULL;
        INI_AO_TOTAL_BRIGHT=calloc(nAO,sizeof(double));
        
        if(INI_AO_SCALING)
        {
            AOscale=calloc(1*nAO,sizeof(double));
            AOscale2=calloc(1*nAO,sizeof(double));
            AOds=calloc(nAOtotal*nAO,sizeof(double));
            nAOscale=nAO;
        }
    }
    int *ncumAOcount;
    double *verbose_ao_fit;
    if(INI_VERBOSE==1 && INI_HAVE_AO)
   {
       verbose_ao_fit=calloc(nAO,sizeof(double));
       ncumAOcount=calloc(nAO+1,sizeof(int));
       for(int j=1;j<=nAO;j++)
           ncumAOcount[j]=ncumAOcount[j-1]+AO->nobs[j-1];
   }
    int nCRtotal=0;
    int nCRcols=0;
    int nCRoffsets=0;
    int nCR=0;
    double* CRoffset,*CRoffset2;
    double *CRdv;
    double *CRout;
    double *CRdoff;
    double CRfit;
    if(INI_HAVE_CNTR)
    {
        if(INI_CNTR_IS_SPARSE)
            nCRtotal=(CR->ntotal);
        else
        nCRtotal=(CR->ntotal)+nvertn*CR->ncont;
        nCRcols=3*nvert+3; //Columns in derivative matrix 3*vertices+angles
        nCRoffsets=2*(CR->ncont);
        nCR=CR->ncont;
        CRout=calloc(nCRtotal,sizeof(double));
        CRoffset=calloc(2*nCR,sizeof(double));
        CRoffset2=calloc(2*nCR,sizeof(double));
        CRdv=calloc(nCRtotal*(3*nvert+3),sizeof(double)); //derivative matrix, vertices+angles
        CRdoff=calloc(nCRtotal*2*nCR,sizeof(double)); //derivative matrix, offsets
    }
   // print_matrix(AO->up,AO->nao,3);
    //Occultation data
    int nOCtotal=0;
    int nOCcols=0;
    int nOC=0;
    int nOCoffsets=0;
    int nChordoffsets=0;
    int nvectorreg=0;
    double *OCoffset,*OCoffset2,*OCout,*OCdv,OCfit=0,*OCdoff;
    double *Chordoffset,*Chordoffset2;
    double *dChordoffset; //Derivative matrix for chord offsets
    double vectorreg; //regularization for free chord offsets
    double *dvectorreg; //derivative matrix for free chord offsets, 1xnChordoffsets 
    if(INI_HAVE_OC)
    {
        nOCtotal=4*(OC->ntotal);
        nOCcols=3*nvert+3;
        nOCoffsets=2*(OC->noc);
        nOC=OC->noc;
        OCoffset=calloc(nOCoffsets,sizeof(double));
        if(INI_OC_OFFSET!=NULL)
            memcpy(OCoffset,INI_OC_OFFSET,sizeof(double)*nOCoffsets);
        OCoffset2=calloc(nOCoffsets,sizeof(double));
        OCdv=calloc(nOCtotal*nOCcols,sizeof(double));
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
    int nRDcols=0;
    int nRD=0;
    int nRDoffsets=0;
    int nRDscale=0;
    int nRDexp=0;
    double RDexp=INI_RDEXP;
    double RDexp2=RDexp;
    double *RDoffset,*RDoffset2,*RDscale,*RDscale2,*RDout,*RDdv,*RDdoff,*RDdscale,*RDdexp;
    if(INI_HAVE_RD)
    {
        nRDtotal=2*(RD->ntotal);
        nRDcols=3*nvert+3;
        nRDoffsets=2*(RD->nRD);
        nRD=RD->nRD;
        nRDscale=nRD;
        nRDexp=0;
        RDoffset=calloc(nRDoffsets,sizeof(double));
        RDoffset2=calloc(nRDoffsets,sizeof(double));
        RDscale=calloc(nRDscale,sizeof(double));
        RDscale2=calloc(nRDscale,sizeof(double));
        RDout=calloc(nRDtotal,sizeof(double));
        RDdv=calloc(nRDtotal*nRDcols,sizeof(double));
        RDdoff=calloc(nRDtotal*nRDoffsets,sizeof(double));
        RDdscale=calloc(nRDtotal*nRD,sizeof(double));
        RDdexp=calloc(nRDtotal,sizeof(double));
        
    }  
        
    //Variables for regularization terms
    double CRres; //Convex reg
    double *dCRdv;
    double *y,*y2;
// printf("AO nobs: %d %d total: %d\n",AO->nobs[0],AO->nobs[1],AO->ntotal);
//           write_matrix_file("/tmp/data0r.txt",AO->datar[0],1,AO->nobs[0]);
//           write_matrix_file("/tmp/data0i.txt",AO->datai[0],1,AO->nobs[0]);
//           write_matrix_file("/tmp/data1r.txt",AO->datar[1],1,AO->nobs[1]);
//           write_matrix_file("/tmp/data1i.txt",AO->datai[1],1,AO->nobs[1]);
//           write_matrix_file("/tmp/freqx.txt",AO->freqx[0],1,AO->nobs[0]);
//           write_matrix_file("/tmp/freqy.txt",AO->freqy[0],1,AO->nobs[0]);
//          write_matrix_file("/tmp/psf0r.txt",AO->psfr[0],1,AO->nobs[0]);
//           write_matrix_file("/tmp/psf0i.txt",AO->psfi[0],1,AO->nobs[0]);
//           write_matrix_file("/tmp/psf1r.txt",AO->psfr[1],1,AO->nobs[1]);
//           write_matrix_file("/tmp/psf1i.txt",AO->psfi[1],1,AO->nobs[1]);
//     write_matrix_file("/tmp/data1r.txt",RD->datar[0],1,RD->nobs[0]);
//            write_matrix_file("/tmp/data1i.txt",RD->datai[0],1,RD->nobs[0]);
//            write_matrix_file("/tmp/freqx.txt",RD->freqx[0],1,RD->nobs[0]);
//            write_matrix_file("/tmp/freqy.txt",RD->freqy[0],1,RD->nobs[0]); 
//     exit(1);
    dCRdv=calloc(3*nvert+3,sizeof(double)); 
    
    double ANGres; //Dihedral angle regularization
    double *dANGdv;
    dANGdv=calloc(3*nvert+3,sizeof(double));
    
    double *Ares; //Area regularization
    double *dAdv;
    
    Ares=calloc(nfacn,sizeof(double));
    
    dAdv=calloc(nfacn*(3*nvert+3),sizeof(double));
    
    
    double *S,*Sp;
    double *J,*JTJ,*JTJpd;
    
   
    int Slength;
    double zm=0;
    double *dzmdz;
    int softmaxz=0;
    if(INI_ZMAX_WEIGHT>0)
    {
        softmaxz=1;
        dzmdz=calloc(nvert,sizeof(double));
    }
    int areareglength=0;
    if(INI_AW>0) //Option to disable Area regularization
    {
        areareglength=nfacn;
        
    }
    Slength=nLCtotal+nAOtotal+nOCtotal+nRDtotal+1+1+areareglength+nvectorreg+nAlbreg+softmaxz+nCRtotal; //LC points+AO points+OC points+RD points+convex reg+dihedral reg+area reg+[vector reg for free chords]+albedo reg+restricted z reg+contour_points
   
    int nJcols=3*nvert+3+ncalib+nAlbedo+nAOoffsets+nAOscale+nOCoffsets+nChordoffsets+nRDoffsets+nRDscale+nRDexp+nCRoffsets; //shape params+offset params and scale
    int OCoffsetcolpos=3*nvert+3+nAOoffsets+nAOscale;
    int Chordoffsetcolpos=3*nvert+3+nAOoffsets+nAOscale+nOCoffsets;
    int OCrowpos=nLCtotal+nAOtotal;
    int AOscalepos=3*nvert+3+nAOoffsets;
    int RDoffsetcolpos=3*nvert+3+nAOoffsets+nAOscale+nOCoffsets+nChordoffsets;
    int RDscalecolpos=3*nvert+3+nAOoffsets+nAOscale+nOCoffsets+nChordoffsets+nRDoffsets;
    int RDrowpos=nLCtotal+nAOtotal+nOCtotal;
    int Albcolpos=3*nvert+3+nAOoffsets+nAOscale+nOCoffsets+nChordoffsets+nRDoffsets+nRDscale+nRDexp;
    int regpos=nLCtotal+nAOtotal+nOCtotal+nRDtotal;
    int phasecolpos=3*nvert+3+nAlbedo+nAOoffsets+nAOscale+nOCoffsets+nChordoffsets+nRDoffsets+nRDscale+nRDexp;
    int Convregpos=regpos;
    int Angregpos=regpos+1;
    int Arearegpos=regpos+2;
    int Zmaxregpos=Arearegpos+areareglength;
    int Vectorregpos=Zmaxregpos+softmaxz;
    int Albregpos=Vectorregpos+nvectorreg;
    int CRoffsetcolpos=3*nvert+3+ncalib+nAlbedo+nAOoffsets+nAOscale+nOCoffsets+nChordoffsets+nRDoffsets+nRDscale+nRDexp;
    int CRrowpos=nLCtotal+nAOtotal+nOCtotal+nRDtotal+1+1+areareglength+nvectorreg+nAlbreg+softmaxz;
    S=calloc(Slength+nJcols,sizeof(double)); 
    Sp=calloc(Slength,sizeof(double));
    J=calloc(Slength*(nJcols),sizeof(double));
    
    JTJ=calloc((nJcols)*(nJcols),sizeof(double));
    JTJpd=calloc((nJcols)*(nJcols),sizeof(double));
    //Check if some parameters are fixed:
    double *Mask_Matrix;
    int *Mask;
    int nMask=0;
    double *MJTJ;
    double *MJTJpd;
    
    double *MX;
    double *JM;
    double *lhs;
    
  double *d_old=calloc(nJcols,sizeof(double));
  double *d=calloc(nJcols,sizeof(double));
    if(INI_MASK_SET==1)
    {
        Mask=calloc(nJcols,sizeof(int));
        
        if(INI_FIX_SHAPE==1)
            for(int k=0;k<3*nvert;k++)
                Mask[k]=1;
        if(INI_FIX_ANGLES==1)
            for(int k=3*nvert;k<3*nvert+3;k++)
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
            Mask[phasecolpos]=INI_PHASE_MASK[0];
            Mask[phasecolpos+1]=INI_PHASE_MASK[1];
            Mask[phasecolpos+2]=INI_PHASE_MASK[2];
            
            
        }
        
        mask_matrix(nJcols,Mask,&Mask_Matrix,&nMask);
        MJTJ=calloc(nMask*nMask,sizeof(double));
        MJTJpd=calloc(nMask*nMask,sizeof(double));
       
        MX=calloc(nMask+Slength,sizeof(double));
        
        JM=calloc(Slength*nMask,sizeof(double));
        
    }
                
        
        
    double *LCout,*dLCdv,*rhs,*X;
    double LCfit;
    double RDfit=0;
    LCout=calloc(nLCtotal,sizeof(double));
    dLCdv=calloc((nLCtotal)*(3*nvert+3),sizeof(double));
    rhs=calloc(nJcols,sizeof(double));
    X=calloc(nJcols,sizeof(double));
    //phase parameters
    double *dLCdp=NULL;
   
    if(INI_PHASE_PARAMS!=NULL)
        dLCdp=calloc(nLCtotal*3,sizeof(double));
    
    //Weights
    
    double AOW=INI_AOW;
    double lcW=INI_LCW;
    double RDW=INI_RW;
    double cW=INI_CW;
    double aW=INI_AW;
    double ocW=INI_OCW;
    double angW=INI_DW;
    double lambda=INI_LAMBDA;
    int decreased=1;
    double chisq=1e9;
    double chisq2;
    double Aresfit;
    double prevdec=10;
    double dec=10;
    double threshold=INI_MINDEC;
    int count=0;
    int DONE=0;
    int alength=0,LMAX=0;
   if(INI_RESTORE_STATE)
  {
      //Restore angles
     
     
      read_state_file(INI_RESTORE_STATE,"#Angles",angles,4);
      angles[0]=(90-angles[0])*PI/180;
      angles[1]=angles[1]*PI/180;
      angles[2]=24*2*PI*1.0/angles[2];
      angles[3]=angles[3]*PI/180;
      //Restore shape
      
      //Restore AO
      read_state_file(INI_RESTORE_STATE,"#AOoffset",AOoffset,nAOoffsets);
      read_state_file(INI_RESTORE_STATE,"#AOscale",AOscale,nAOscale);
      //Restore OC
      read_state_file(INI_RESTORE_STATE,"#OCCoffset",OCoffset,nOCoffsets);
      read_state_file(INI_RESTORE_STATE,"#Chordoffset",Chordoffset,nChordoffsets);
      //Restore albedo
      read_state_file(INI_RESTORE_STATE,"#Albedolog",eAlbedo,nfacn);
    read_state_file(INI_RESTORE_STATE,"#AOtotalbright",INI_AO_TOTAL_BRIGHT,nAO);
  }
  //  print_matrix(INI_AO_TOTAL_BRIGHT,1,nAO);
       
            if(INI_DW_DEC!=1.0)
            angW=INI_DW*pow(INI_DW_DEC,count);
            Sqrt3_Subdiv(tlist,vlist,nfac,nvert,&tlistn,&vlistn,&nfacn,&nvertn,&D,INI_SD_LEVEL); //Do  subdivision. Remember to free allocated mem
            //Calculate LCs
            
            calculate_lcs(tlistn,vlistn,nfacn,nvertn,angles,LC,D,nvertn,nvert,LCout,dLCdv,eAlbedo,Alblimits,dLCdalb,params,dLCdp,0);
         
            /////////////DEBUG////////////////////////////////
//               write_shape_file("/tmp/Inishape.txt",tlistn,vlistn,nfacn,nvertn);
//               write_matrix_file("/tmp/D.txt",D,nvertn,nvert);
//                      write_matrix_file("/tmp/LCout.txt",LCout,1,nLCtotal);
 //                     write_matrix_file("/tmp/dLCdv.txt",dLCdv,nLCtotal,3*nvert+3);
//                      write_matrix_file("/tmp/dLCdp.txt",dLCdp,nLCtotal,4);
            ///////////DEBUG/////////////////////////////////////
            if(INI_FIT_ALBEDO==1)
            {
                mult_with_cons(dLCdalb,nLCtotal,nfacn,lcW);
                localsmooth(tlistn,vlistn,nfacn,nvertn,eAlbedo,Alblimits,Albreg,dAlbreg);
                mult_with_cons(dAlbreg,nfacn,nfacn,-INI_ALBREGW);
                mult_with_cons(Albreg,1,nfacn,INI_ALBREGW);
                set_submatrix(S,1,Slength,Albreg,1,nfacn,0,Albregpos);
             //   set_submatrix(J,Slength,nJcols,dLCdalb,nLCtotal,nfacn,0,Albcolpos);
              //  set_submatrix(J,Slength,nJcols,dAlbreg,nfacn,nfacn,Albregpos,Albcolpos);
            }
                
            //AO data
            if(INI_HAVE_AO)
            {
                ////////////DEBUG/////////////REMOVE THIS!!!!!!!!!!!!!!!
//                 AO->psfr[0]=NULL;
//                 AO->psfi[0]=NULL;
//                 AO->psfr[1]=NULL;
//                 AO->psfi[1]=NULL;
//                 zero_array(AO->datar[0],AO->nobs[0]);
//                 zero_array(AO->datai[0],AO->nobs[0]);
//                 zero_array(AO->datar[1],AO->nobs[1]);
//                 zero_array(AO->datai[1],AO->nobs[1]);
                ///////////DEBUG////////////////////////////////////////
              
                Calculate_AOs_fit(tlistn,vlistn,nfacn,nvertn,angles,AO,AOoffset,D,nvertn,nvert,INI_AO_WEIGHT,AOscale,AOout,AOdv,AOds,eAlbedo,Alblimits,NULL,0);
                /////////////DEBUG//////////////////
//                   write_matrix_file("/tmp/AOout.txt",AOout,1,nAOtotal);
//                        write_matrix_file("/tmp/AOdv.txt",AOdv,nAOtotal,nAOcols);
                /////////////DEBUG////////////////
                mult_with_cons(AOout,1,nAOtotal,AOW);
              //  mult_with_cons(AOdv,nAOtotal,nAOcols,AOW);
                set_submatrix(S,1,Slength,AOout,1,nAOtotal,0,nLCtotal); //Construct Jacobian matrix
               // set_submatrix(J,Slength,nJcols,AOdv,nAOtotal,nAOcols,nLCtotal,0);
                matrix_transprod(AOout,nAOtotal,1,&AOfit);
                         
                if(INI_AO_SCALING)
                {
                 //   mult_with_cons(AOds,nAOtotal,nAO,AOW);
                 //   set_submatrix(J,Slength,nJcols,AOds,nAOtotal,nAO,nLCtotal,nAOcols);
                }
                if(INI_VERBOSE)
                    for(int k=0;k<nAO;k++)
                        for(int k2=0;k2<AO->nobs[k];k2++)
                            verbose_ao_fit[k]+=pow(AOout[k2+ncumAOcount[k]],2)+pow(AOout[k2+ncumAOcount[k]+AO->ntotal],2);
                        
                    
            }
            if(INI_HAVE_CNTR)
            {
                Calculate_Contours(tlistn,vlistn,nfacn,nvertn,angles,CR,CRoffset,D,nvertn,nvert,CRout,CRdv,CRdoff);
               
                mult_with_cons(CRout,1,nCRtotal,INI_CNTR_WEIGHT);
                mult_with_cons(CRdv,nCRtotal,3*nvert+3,-INI_CNTR_WEIGHT);
                mult_with_cons(CRdoff,nCRtotal,nCRoffsets,-INI_CNTR_WEIGHT);
                set_submatrix(S,1,Slength,CRout,1,nCRtotal,0,CRrowpos);
               // set_submatrix(J,Slength,nJcols,CRdv,nCRtotal,3*nvert+3,CRrowpos,0);
               // set_submatrix(J,Slength,nJcols,CRdoff,nCRtotal,nCRoffsets,CRrowpos,CRoffsetcolpos);
                matrix_transprod(CRout,nCRtotal,1,&CRfit);
         
            }
            if(INI_HAVE_OC)
            {
                if(INI_FREE_CHORD_NMR>0)
                {
                    calculate_OCs(tlistn,vlistn,nfacn,nvertn,angles,OC,OCoffset,INI_OC_WEIGHT,D,nvertn,nvert,Chordoffset,OCout,OCdv,OCdoff,dChordoffset);
                    vector_regularization(Chordoffset,nChordoffsets,&vectorreg,dvectorreg);
                    mult_with_cons(dChordoffset,nOCtotal,nChordoffsets,-ocW);
                    mult_with_cons(dvectorreg,1,nChordoffsets,-INI_CHRDW);
                    vectorreg*=INI_CHRDW;
                    S[Vectorregpos]=vectorreg; //NOTE ABSOLUTE ADDRESS HERE. FIX!
                 //   set_submatrix(J,Slength,nJcols,dChordoffset,nOCtotal,nChordoffsets,OCrowpos,Chordoffsetcolpos); //Chord offsets
                 //   set_submatrix(J,Slength,nJcols,dvectorreg,1,nChordoffsets,Vectorregpos,Chordoffsetcolpos);
                }
                    
                else
                    calculate_OCs(tlistn,vlistn,nfacn,nvertn,angles,OC,OCoffset,INI_OC_WEIGHT,D,nvertn,nvert,Chordoffset,OCout,OCdv,OCdoff,NULL);
                 /////////////DEBUG////////////////////////////////
//                  write_matrix_file("/tmp/OCout.txt",OCout,1,nOCtotal);
//                  write_matrix_file("/tmp/OCdv.txt",OCdv,nOCtotal,nOCcols);
//                  write_matrix_file("/tmp/OCdoff.txt",OCdoff,nOCtotal,nOCoffsets);
//                 write_matrix_file("/tmp/dChordoffset.txt",dChordoffset,nOCtotal,nChordoffsets);
//                 exit(1);
                /////////////DEBUG////////////////////////////////
                mult_with_cons(OCout,1,nOCtotal,ocW);
                //mult_with_cons(OCdv,nOCtotal,nOCcols,-ocW); //Note the minus here.
               // mult_with_cons(OCdoff,nOCtotal,nOCoffsets,-ocW);
                set_submatrix(S,1,Slength,OCout,1,nOCtotal,0,OCrowpos);
                //set_submatrix(J,Slength,nJcols,OCdv,nOCtotal,nOCcols,OCrowpos,0);
                //set_submatrix(J,Slength,nJcols,OCdoff,nOCtotal,nOCoffsets,OCrowpos,OCoffsetcolpos);
                matrix_transprod(OCout,nOCtotal,1,&OCfit);
                
            }
            if(INI_HAVE_RD)
            {
                Calculate_RDs(tlistn,vlistn,nfacn,nvertn,angles,RD,RDoffset,D,nvertn,nvert,INI_RD_WEIGHT,RDscale,RDexp,RDout,RDdv,RDdoff,RDdscale,RDdexp,0);
                //printf("RDW: %f RDexp: %f\n",RDW,RDexp);
                ///DEBUG
//                 write_matrix_file("/tmp/RDout.txt",RDout,1,nRDtotal);
//                 write_matrix_file("/tmp/RDdv.txt",RDdv,nRDtotal,nRDcols);
//                 write_matrix_file("/tmp/RDdoff.txt",RDdoff,nRDtotal,nRDoffsets);
//                 write_matrix_file("/tmp/RDdscale.txt",RDdscale,nRDtotal,nRDscale);
//                 write_matrix_file("/tmp/RDdexp.txt",RDdexp,1,nRDtotal);
                /////DEBUG///////////////////////////////
                mult_with_cons(RDout,1,nRDtotal,RDW);
              //  mult_with_cons(RDdv,nRDtotal,nRDcols,RDW); //NO minus here
              //  mult_with_cons(RDdoff,nRDtotal,nRDoffsets,RDW);
               // mult_with_cons(RDdscale,nRDtotal,nRDscale,RDW);
                
                set_submatrix(S,1,Slength,RDout,1,nRDtotal,0,RDrowpos);
              //  set_submatrix(J,Slength,nJcols,RDdv,nRDtotal,nRDcols,RDrowpos,0);
              //  set_submatrix(J,Slength,nJcols,RDdoff,nRDtotal,nRDoffsets,RDrowpos,RDoffsetcolpos); //[RDdv 0 RDdoff RDdscale RDdexp]
               // set_submatrix(J,Slength,nJcols,RDdscale,nRDtotal,nRDscale,RDrowpos,RDscalecolpos);
                if(nRDexp>0)
                {
                    mult_with_cons(RDdexp,nRDtotal,1,RDW);
                    set_submatrix(J,Slength,nJcols,RDdexp,nRDtotal,1,RDrowpos,RDscalecolpos+nRDscale);
                }
                matrix_transprod(RDout,nRDtotal,1,&RDfit);
//                 printf("RDdoff pos: %d %d\n",RDrowpos,RDoffsetcolpos);
//                 printf("RDscale pos: %d %d\n",RDrowpos,RDscalecolpos);
//                 printf("RDdexp pos: %d %d\n",RDrowpos,RDscalecolpos+nRDscale);
//                 printf("RDfit: %f\n",RDfit);
                
            }
            convex_reg(tlistn,vlistn,nfacn,nvertn,D,nvertn,nvert,&CRres,dCRdv);
            ///////////DEBUG///////////////////////////
            // write_matrix_file("/tmp/dCRdv.txt",dCRdv,1,3*nvert+3);
            ///////////DEBUG///////////////////////////
            dihedral_angle_reg(tlistn,vlistn,nfacn,nvertn,D,nvertn,nvert,&ANGres,dANGdv);
            ///////////DEBUG///////////////////////////
            // write_matrix_file("/tmp/dANGdv.txt",dANGdv,1,3*nvert+3);
            ///////////DEBUG///////////////////////////
            if(INI_AW>0)
            {
                area_reg(tlistn,vlistn,nfacn,nvertn,D,nvertn,nvert,Ares,dAdv); 
                mult_with_cons(Ares,1,nfacn,aW);
              //  mult_with_cons(dAdv,nfacn,3*nvert+3,-aW);
                set_submatrix(S,1,Slength,Ares,1,nfacn,0,Arearegpos); //Area regularization
             //   set_submatrix(J,Slength,nJcols,dAdv,nfacn,3*nvert+3,Arearegpos,0);
                matrix_transprod(Ares,nfacn,1,&Aresfit);
            } 
            
            ///////////DEBUG///////////////////////////
//             write_matrix_file("/tmp/dAdv.txt",dAdv,nfacn,3*nvert+3);
            ///////////DEBUG///////////////////////////
           if(INI_ZMAX_WEIGHT>0)
           {
               soft_maxdimz(tlistn,vlistn,nfacn,nvertn,D,nvert,INI_ZMAX,1.0,&zm,dzmdz);
               zm*=INI_ZMAX_WEIGHT;
             //  mult_with_cons(dzmdz,1,nvert,-INI_ZMAX_WEIGHT);
               S[Zmaxregpos]=zm;
           //    set_submatrix(J,Slength,nJcols,dzmdz,1,nvert,Zmaxregpos,2*nvert);
           }
            mult_with_cons(LCout,1,nLCtotal,lcW);
          //  mult_with_cons(dLCdv,nLCtotal,3*nvert+3,lcW);
            CRres*=cW;
            
         //   mult_with_cons(dCRdv,1,3*nvert+3,-cW);
            ANGres*=angW;
         //   mult_with_cons(dANGdv,1,3*nvert+3,-angW);
           
             
            
            //Build the res vector and matrix
            set_submatrix(S,1,Slength,LCout,1,nLCtotal,0,0);
            
            
            
            S[Convregpos]=CRres; //Convex reg terms
            S[Angregpos]=ANGres;
             
            
            
           
            matrix_transprod(LCout,nLCtotal,1,&LCfit);
            matrix_transprod(S,Slength,1,&chisq);
            
            printf("chisq: %4.2f LCfit: %4.2f AOfit: %4.2f OCfit: %4.2f RDfit: %4.2f Convex reg: %4.2f Dihedral Angle reg: %4.2f",chisq,LCfit,AOfit,OCfit,RDfit,pow(CRres,2),pow(ANGres,2)); 
            if(INI_HAVE_CNTR>0)
                printf(" CNTRfit: %4.2f\n",CRfit);
            else
                printf("\n");
            //Construct Jacobian matrix
           if(INI_VERBOSE && INI_HAVE_AO)
           {
               printf("AOfits: ");
               for(int j=0;j<nAO;j++)
                   printf(" %7.4f",verbose_ao_fit[j]);
               printf("\n");
           }
                   
       
        if(INI_SD_LEVEL!=0)
            free(D);
        //     zero_array(dLCdv,nLCtotal*(3*nvert+3));
        //     zero_array(AOdv,nAOcols*nAOtotal);
        //     zero_array(dCRdv,3*nvert+3);
        //     zero_array(dANGdv,3*nvert+3);
        //     zero_array(dAdv,nfacn*(3*nvert+3));
        
    
    free(tlist);
    free(vlist);
    free(dLCdv);
    free(dCRdv);
    free(dANGdv);
    //free(AOdv);
    free(dAdv);
   
    
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
