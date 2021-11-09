#include<stdio.h>
#include<stdlib.h>
#include"utils.h"
#include"matrix_ops.h"
#include"structs.h"
#include"globals.h"

void fit_model_albedo(LCstruct *LC,AOstruct *AO,double** feAlbedo,int* alblength)
{
    double *INI_AO_TOTAL_BRIGHT;
    //First initialize the initial shape
    int *tlist,*tlistn,nfac,nvert,nfacn,nvertn;
    double *vlist,*vlist2,*vlistn,*D;
    int nLCtotal=0;
    void (*SubdivFunc)(int*,double*,int,int,int **,double **,int *,int*,double **,int);
    if(INI_SUBDIV_TYPE_BUTTERFLY==1)
        SubdivFunc=&butterfly_subdiv;
    else
        SubdivFunc=&Sqrt3_Subdiv;
            
    if(INI_SHAPE_FILE!=NULL)
        read_shape(INI_SHAPE_FILE,&tlist,&vlist,&nfac,&nvert,0);
    else if(INI_RESTORE_STATE)
    {
         read_state_fileI(INI_RESTORE_STATE,"#NFacets",&nfac,1);
         read_state_fileI(INI_RESTORE_STATE,"#NVertices",&nvert,1);
         vlist=calloc(3*nvert,sizeof(double));
         tlist=calloc(3*nfac,sizeof(int));
         int rv=0;
         rv=read_state_fileI(INI_RESTORE_STATE,"#PolyhedronFacets",tlist,3*nfac);
         read_state_file(INI_RESTORE_STATE,"#PolyhedronVertices",vlist,3*nvert);
        nfacn=nfac;
        nvertn=nvert;
    }
   
    *alblength=nfac;
    
    if(INI_SHAPE_FILE==NULL && INI_RESTORE_STATE==0)
    {
        fprintf(stderr,"Error: Either INI_SHAPE_FILE or INI_RESTORE_STATE must be set for albedo fitting\n");
        exit(-1);
    }
   
    
    vlist2=calloc(3*nvert,sizeof(double));
    //Setup angles
    double angles[]={(90-INI_ANGLE_B)*PI/180,INI_ANGLE_L*PI/180,24*2*PI*1.0/INI_ANGLE_P,INI_ANGLE_PHI0*PI/180};
    double angles2[4]={0,0,0,0};
    angles2[3]=angles[3];
   
    
    //This is for calibrated lightcurves
    double *params=NULL,*params2=NULL;
    int ncalib=0;
    int USE_CALIB=0;
    int nAlbedo=0;
    double *eAlbedo=NULL,*eAlbedo2=NULL;
    double *dLCdalb=NULL;
    double *Alblimits=NULL;
    double *Albreg=NULL,*dAlbreg=NULL;
    int nAlbreg=0;   
    if(INI_HAVE_LC)
    {
        nLCtotal=LC->ntotal;
        
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
        USE_CALIB=LC->calib;
        
        
        
        //This is for calibrated lightcurves end
        //LC albedo variegation
        
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
    }
    
    int nAOtotal=0;
    int nAOcols=0;
    int nAO=0;
    int nAOoffsets=0;
    int nAOscale=0;
    double *AOoffset,*AOoffset2,*AOout,*AOdv,AOfit=0,*AOscale,*AOscale2,*AOds,*dAOdAlb,*ao_total_bright;
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
        ao_total_bright=calloc(nAO,sizeof(double));
        
        if(INI_AO_SCALING)
        {
            AOscale=calloc(1*nAO,sizeof(double));
            AOscale2=calloc(1*nAO,sizeof(double));
            AOds=calloc(nAOtotal*nAO,sizeof(double));
            nAOscale=nAO;
        }
        if(INI_FIT_AO_ALBEDO)
        {
            dAOdAlb=calloc(nAOtotal*nfacn,sizeof(double));
        }
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
    
    int nRDtotal=0;
    int nRDcols=0;
    int nRD=0;
    int nRDoffsets=0;
    int nRDscale=0;
    int nRDexp=0;
    double RDexp=INI_RDEXP;
    double RDexp2=RDexp;
    double *RDoffset,*RDoffset2,*RDscale,*RDscale2,*RDout,*RDdv,*RDdoff,*RDdscale,*RDdexp;
   
    
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
    int inerreglength=0;
    double inerres,inerangle;
    double *inerdvi=NULL,*inerdv=NULL;
    if(INI_INER_WEIGHT>0)
    {
        inerreglength=1;
        inerdv=calloc(3*nvert+3,sizeof(double));
        inerdvi=calloc(3*nvertn+3,sizeof(double));
    }
    int areareglength=0;
    if(INI_AW>0) //Option to disable Area regularization
    {
        areareglength=nfacn;
        
    }
    int nCOM=0;
    double CoM=0.0;
    double *CoMdv;
    int info=0;
    if(INI_COM_WEIGHT>0)
    {
        nCOM=1;
        CoMdv=calloc(3*nvert+3,sizeof(double));
        
    }
    Slength=nLCtotal+nAOtotal+nOCtotal+nRDtotal+1+1+areareglength+nvectorreg+nAlbreg+softmaxz+nCRtotal+nCOM+inerreglength; //LC points+AO points+OC points+RD points+convex reg+dihedral reg+area reg+[vector reg for free chords]+albedo reg+restricted z reg+contour_points
    //Albedo fit only, no shape regularization
    Slength=nLCtotal+nAOtotal+nAlbreg;
    int nJcols=3*nvert+3+ncalib+nAlbedo+nAOoffsets+nAOscale+nOCoffsets+nChordoffsets+nRDoffsets+nRDscale+nRDexp+nCRoffsets; //shape params+offset params and scale
    nJcols=ncalib+nAlbedo;
    int OCoffsetcolpos=3*nvert+3+nAOoffsets+nAOscale;
    int Chordoffsetcolpos=3*nvert+3+nAOoffsets+nAOscale+nOCoffsets;
    int OCrowpos=nLCtotal+nAOtotal;
    int AOscalepos=3*nvert+3+nAOoffsets;
    int RDoffsetcolpos=3*nvert+3+nAOoffsets+nAOscale+nOCoffsets+nChordoffsets;
    int RDscalecolpos=3*nvert+3+nAOoffsets+nAOscale+nOCoffsets+nChordoffsets+nRDoffsets;
    int RDrowpos=nLCtotal+nAOtotal+nOCtotal;
    int Albcolpos=3*nvert+3+nAOoffsets+nAOscale+nOCoffsets+nChordoffsets+nRDoffsets+nRDscale+nRDexp;
    //We are fitting albedo only,
    Albcolpos=0;
    int regpos=nLCtotal+nAOtotal+nOCtotal+nRDtotal;
    regpos=nLCtotal+nAOtotal;
    int phasecolpos=3*nvert+3+nAlbedo+nAOoffsets+nAOscale+nOCoffsets+nChordoffsets+nRDoffsets+nRDscale+nRDexp;
    int Convregpos=regpos;
    int Angregpos=regpos+1;
    int Arearegpos=regpos+2;
    int Zmaxregpos=Arearegpos+areareglength;
    int Vectorregpos=Zmaxregpos+softmaxz;
    int Albregpos=Vectorregpos+nvectorreg;
    Albregpos=regpos;
    int CRoffsetcolpos=3*nvert+3+ncalib+nAlbedo+nAOoffsets+nAOscale+nOCoffsets+nChordoffsets+nRDoffsets+nRDscale+nRDexp;
    int CRrowpos=nLCtotal+nAOtotal+nOCtotal+nRDtotal+1+1+areareglength+nvectorreg+nAlbreg+softmaxz;
    int CoMregpos=CRrowpos+nCRtotal;
    int inerregpos=CoMregpos+nCOM;
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
    int first_time=1;
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
                if(INI_FIX_VERTEX_NBR>0)
                {
                    for(int k=0;k<INI_FIX_VERTEX_NBR;k++)
                        if(INI_FIX_VERTEX_LIST[k]<=nvert)
                        {
                            Mask[INI_FIX_VERTEX_LIST[k]-1]=1;
                            Mask[INI_FIX_VERTEX_LIST[k]-1+nvert]=1;
                            Mask[INI_FIX_VERTEX_LIST[k]-1+2*nvert]=1;
                        }
                        
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
        printf("Restoring previous state from %s\n",INI_RESTORE_STATE);
        
        
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
      
        //Restore albedo
        int albcount=0;
        if(INI_HAVE_LC && INI_FIT_ALBEDO)
        {
            printf("Restoring albedo\n");
            albcount=read_state_file(INI_RESTORE_STATE,"#Albedolog",eAlbedo,nfacn);
            if(albcount!=nfacn)
            {
                printf("There are %d albedo values, but %d facets. Discarding loaded albedos\n",albcount,nfacn);
                zero_array(eAlbedo,nfacn);
            }
        }
    }
    vlistn=vlist;
    tlistn=tlist;
    nfacn=nfac;
    nvertn=nvert;
    for(int k=0;k<NUM_OF_ROUNDS;k++)
    {
        start:
        if(decreased==1)
        {
            if(INI_DW_DEC!=1.0)
                angW=INI_DW*pow(INI_DW_DEC,count);
          //  (*SubdivFunc)(tlist,vlist,nfac,nvert,&tlistn,&vlistn,&nfacn,&nvertn,&D,INI_SD_LEVEL); //Do  subdivision. Remember to free allocated mem
            
            //Calculate LCs
            if(INI_HAVE_LC)
            {
                calculate_lcs(tlistn,vlistn,nfacn,nvertn,angles,LC,NULL,nvertn,nvert,LCout,dLCdv,eAlbedo,Alblimits,dLCdalb,params,dLCdp,1);
                
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
                    set_submatrix(J,Slength,nJcols,dLCdalb,nLCtotal,nfacn,0,Albcolpos);
                    set_submatrix(J,Slength,nJcols,dAlbreg,nfacn,nfacn,Albregpos,Albcolpos);
                }
                mult_with_cons(LCout,1,nLCtotal,lcW);
                //mult_with_cons(dLCdv,nLCtotal,3*nvert+3,lcW);
                set_submatrix(S,1,Slength,LCout,1,nLCtotal,0,0);
                matrix_transprod(LCout,nLCtotal,1,&LCfit);
                //set_submatrix(J,Slength,nJcols,dLCdv,nLCtotal,3*nvert+3,0,0);
                
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
                
                Calculate_AOs(tlistn,vlistn,nfacn,nvertn,angles,AO,AOoffset,NULL,nvertn,nvert,INI_AO_WEIGHT,AOscale,AOout,AOdv,AOds,eAlbedo,Alblimits,dAOdAlb,1);
                /////////////DEBUG//////////////////
                //                   write_matrix_file("/tmp/AOout.txt",AOout,1,nAOtotal);
                //                        write_matrix_file("/tmp/AOdv.txt",AOdv,nAOtotal,nAOcols);
                /////////////DEBUG////////////////
                mult_with_cons(AOout,1,nAOtotal,AOW);
                mult_with_cons(AOdv,nAOtotal,nAOcols,AOW);
                set_submatrix(S,1,Slength,AOout,1,nAOtotal,0,nLCtotal); //Construct Jacobian matrix
               // set_submatrix(J,Slength,nJcols,AOdv,nAOtotal,nAOcols,nLCtotal,0);
                matrix_transprod(AOout,nAOtotal,1,&AOfit);
                
               
                if(INI_FIT_AO_ALBEDO)
                {
                    mult_with_cons(dAOdAlb,nAOtotal,nfacn,AOW);
                    set_submatrix(J,Slength,nJcols,dAOdAlb,nAOtotal,nfacn,nLCtotal,Albcolpos);
                }
                
            }
            
            
            
            
            
            matrix_transprod(S,Slength,1,&chisq);
            
            printf("chisq: %4.2f LCfit: %4.2f AOfit: %4.2f",chisq,LCfit,AOfit,OCfit); 
           
                printf("\n");
            //Construct Jacobian matrix
            
            
            
            
            
            //include lc parameter derivatives
            if(INI_HAVE_LC && INI_PHASE_PARAMS!=NULL)
            {
                mult_with_cons(dLCdp,nLCtotal,3,lcW);
                set_submatrix(J,Slength,nJcols,dLCdp,nLCtotal,3,0,phasecolpos);
            }
           
            //Calculate J^TJ+lambda*J
            // matrix_transprod(J,Slength,nJcols,JTJ);
           
            /////////////DEBUG////////////////////////////////
            //  write_matrix_file("/tmp/S.txt",S,Slength,1);
            // write_matrix_file("/tmp/J.txt",J,Slength,nJcols);
            //  printf("Convex pos: %d Angpos: %d CRpos: %d CRoffsetcolpos: %d\n",Convregpos,Angregpos,CRrowpos,CRoffsetcolpos);
            //                          exit(1);
            /////////////DEBUG////////////////////////////////     
            matrix_transprod(J,Slength,nJcols,JTJ);
            //  matrix_vectorprod(J,Slength,nJcols,S,rhs,1); //rhs=J^T*S;
            memcpy(Sp,S,sizeof(double)*Slength);
            if(INI_MASK_SET==1)
            {
                matrix_prod_ATBA(Mask_Matrix,nJcols,nMask,JTJ,MJTJ);
                
                matrix_prod(J,Slength,nJcols,Mask_Matrix,nMask,JM); //JM is SLength x nMask matrix
                //  matrix_prod_ATB(JM,Slength,nMask,S,1,Mrhs); //Mrhs is nMask x 1 matrix
            }
            
        }
        if(INI_MASK_SET==1)
        {
            // matrix_adddiag(MJTJ,MJTJpd,nMask,lambda);
            matrix_max_diag(d_old,nMask,MJTJ,d);
            memcpy(d_old,d,sizeof(double)*nMask);
            zero_array(MX,nMask+Slength);
            memcpy(MX,Sp,sizeof(double)*Slength);
            matrix_concat_special2(JM,Slength,nMask,d,lambda,&lhs);
            info=solve_matrix_eq_QR(lhs,Slength+nMask,nMask,MX);
            free(lhs);
            matrix_prod(Mask_Matrix,nJcols,nMask,MX,1,X);
            
            
            //solve_matrix_eqS(MJTJpd,nMask,Mrhs,MX);
            //matrix_prod(Mask_Matrix,nJcols,nMask,MX,1,X);
        }
        else
        {
            //  matrix_adddiag(JTJ,JTJpd,nJcols,lambda); //JTJpd=JTJ+lambda*diag(JTJ)
            
            matrix_max_diag(d_old,nJcols,JTJ,d);
            memcpy(d_old,d,sizeof(double)*nJcols);
            matrix_concat_special2(J,Slength,nJcols,d,lambda,&lhs);
            // matrix_concat_special(J,Slength,nJcols,JTJ,lambda,&lhs); //S has been allocated to Slength+nJcols
            zero_array(S,nJcols+Slength);
            memcpy(S,Sp,sizeof(double)*Slength);
            info=solve_matrix_eq_QR(lhs,Slength+nJcols,nJcols,S); //Solve the LM 
            free(lhs);
            memcpy(X,S,nJcols*sizeof(double));
        }
       
        
       // add_vector_to_vlist(vlist,X,vlist2,nvert);
        
        
        if(INI_HAVE_LC && INI_FIT_ALBEDO==1)
        {
            for(int j=0;j<nfacn;j++)
                eAlbedo2[j]=eAlbedo[j]+X[Albcolpos+j];
        }
        if(INI_HAVE_LC && INI_PHASE_PARAMS!=NULL)
        {
            params2[0]=params[0]+X[phasecolpos];
            params2[1]=params[1]+X[phasecolpos+1];
            params2[2]=params[2]+X[phasecolpos+2];
            // params2[3]=params[3]+X[phasecolpos+3];
        }
        
        if(INI_HAVE_LC)
        {
            
        calculate_lcs(tlistn,vlistn,nfacn,nvertn,angles,LC,D,nvertn,nvert,LCout,dLCdv,eAlbedo2,Alblimits,dLCdalb,params2,dLCdp,0);
        
        if(INI_FIT_ALBEDO==1)
        {
            localsmooth(tlistn,vlistn,nfacn,nvertn,eAlbedo2,Alblimits,Albreg,dAlbreg);
            mult_with_cons(Albreg,1,nfacn,INI_ALBREGW);
            set_submatrix(S,1,Slength,Albreg,1,nfacn,0,Albregpos);
        }
         mult_with_cons(LCout,1,nLCtotal,lcW);
         set_submatrix(S,1,Slength,LCout,1,nLCtotal,0,0);
        }
        if(INI_HAVE_AO)
        {
            Calculate_AOs(tlistn,vlistn,nfacn,nvertn,angles,AO,AOoffset,NULL,nvertn,nvert,INI_AO_WEIGHT,AOscale,AOout,AOdv,NULL,eAlbedo2,Alblimits,dAOdAlb,0);
            mult_with_cons(AOout,1,nAOtotal,AOW);
            set_submatrix(S,1,Slength,AOout,1,nAOtotal,0,nLCtotal);
        }
        
       
        
        
        matrix_transprod(S,Slength,1,&chisq2);
        //matrix_transprod(AOout,nAOtotal,1,&AOfit);
        printf("Round: %d chisq2: %4.2f lambda: %4.2f \n",k+1,chisq2,lambda);
        //printf("k=%d\n",k);
        if(chisq2<chisq)
        {
            prevdec=dec;
            dec=chisq-chisq2;
            count++;
            printf("decreased!\n");
            
            //   printf("scales: %f %f\n",AOscale2[0],AOscale2[1]);
            //memcpy(vlist,vlist2,sizeof(double)*(3*nvert));
            
            // printf("Offsets: %f %f %f %f\n",AOoffset[0],AOoffset[1],AOoffset[2],AOoffset[3]);
            decreased=1;
           
            if(INI_FIT_ALBEDO==1)
                memcpy(eAlbedo,eAlbedo2,sizeof(double)*nfacn);
            if(INI_PHASE_PARAMS!=NULL)
                memcpy(params,params2,sizeof(double)*3);
            
            lambda=lambda/INI_LAMBDADEC;
            zero_array(J,Slength*nJcols);
            zero_array(S,Slength);
            write_shape_file(OUT_SHAPE_FILE,tlistn,vlistn,nfacn,nvertn);
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
            printf("Lambda is larger than LambdaMax, stopping loop\n");
            DONE=1;
        }
        if(k==NUM_OF_ROUNDS-1 || DONE==1)
            break;
        
        
    }
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
    *feAlbedo=eAlbedo;
}
