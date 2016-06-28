#include"utils.h"
#include"matrix_ops.h"
void main()
{
//     int tlist[]={1,2,3,
//         1,3,4,
//         2,4,3,
//         1,2,4}; //4 facets
//     double vlist[]={0.0,-2.0,0.0
//         ,0.5,0.0,-1.0,
//         0.0,1.0,1.0,
//         -3,1,4};
//     
//     int nvert=4;
//     int nfac=4;
    int *tlist;
    double *vlist;
    int nfac,nvert;
    char file[]="mshape.txt";
    read_shape(file,&tlist,&vlist,&nfac,&nvert,0);
    printf("nfac: %d nvert: %d\n",nfac,nvert);
    int nobs[]={29,29,29};
    int nao=3;
    int ntpoints=3*29;
    //double E[]={1,0,0};
    double E2[]={1,0.1,0.1};
    double E[9];
    E[0]=1;
    E[1]=0;
    E[2]=0;
    E[6]=1;
    E[7]=0;
    E[8]=0;
    double norm=NORM(E2);
    //printf("norm: %f\n",norm);
    for(int j=0;j<3;j++)
        E[j+3]=E2[j]/norm;
    double E0[]={1,0,0,1,0,0,1,0,0};
    double TIME[]={0.1,0.2,-0.1};
    double distance[]={0.00137879506,0.00137879506,0.00137879506};
    double scale[]={1,1,1};
    double up[]={0,0,1,0,0,1,0,0,1};
    double *datar=calloc(29,sizeof(double));
    double *datai=calloc(29,sizeof(double));
    double freqx[]={-1.0000,   -0.9300,   -0.8600,   -0.7900,   -0.7200,   -0.6500,   -0.5800,   -0.5100,   -0.4400,   -0.3700,   -0.3000,
        -0.2300,   -0.1600,   -0.0900,   -0.0200,    0.0500,    0.1200,    0.1900,    0.2600,    0.3300,    0.4000,
        0.4700,    0.5400,    0.6100,    0.6800,    0.7500,    0.8200,    0.8900,    0.9600};
    double freqy[]={1.2900,    1.2200,    1.1500,    1.0800,    1.0100,    0.9400,    0.8700,    0.8000,    0.7300,    0.6600,
        0.5900,    0.5200,    0.4500,    0.3800,    0.3100,    0.2400,    0.1700,    0.1000,    0.0300,   -0.0400,   -0.1100,
        -0.1800,   -0.2500,   -0.3200,   -0.3900,   -0.4600,   -0.5300,   -0.6000,   -0.6700,
    };
    double *psf=calloc(nobs[0],sizeof(double));
    for(int j=0;j<nobs[0];j++)
        psf[j]=2;
    double freqy2[]={-0.3,0.05};
    double freqy3[]={-0.5,-0.1};
    double freqx2[]={0.1,0.15};
    double angles[]={0.1,0.3,30,0};
    double offset[]={0.1,0.2,0.5,-0.1,0,-0.3};
    double D[]={1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1};
    double *Weight;
    double *FT,*FTdv;
    double *AOscale=calloc(nao,sizeof(double));
    AOscale[0]=1;
    AOscale[1]=0.5;
    AOscale[2]=-0.2;
    double *AOds=calloc(nao*2*ntpoints,sizeof(double));
    FT=calloc(2*ntpoints,sizeof(double));
    FTdv=calloc(2*ntpoints*(3*nvert+2*nao+3),sizeof(double));
    AOstruct AO;
    AO.nao=nao;
    AO.nobs=nobs;
    AO.datar=calloc(nao,sizeof(double*));
    AO.datai=calloc(nao,sizeof(double*));
    AO.freqx=calloc(nao,sizeof(double*));
    AO.freqy=calloc(nao,sizeof(double*));
    AO.psfr=calloc(nao,sizeof(double*));
    AO.psfi=calloc(nao,sizeof(double*));
    printf("psf: %f %f\n",psf[0],psf[nobs[0]-1]);
    AO.psfr[0]=psf;
    AO.psfi[0]=psf;
    AO.datar[0]=datar;
    AO.datai[0]=datai;
    AO.freqx[0]=freqx;
    AO.freqy[0]=freqy;
    AO.datar[1]=datar;
    AO.datai[1]=datai;
    AO.freqx[1]=freqx;
    AO.freqy[1]=freqy;
    AO.datar[2]=datar;
    AO.datai[2]=datai;
    AO.freqx[2]=freqx;
    AO.freqy[2]=freqy;
    AO.E=E;
    AO.E0=E0;
    AO.TIME=TIME;
    AO.distance=distance;
    AO.scalex=scale;
    AO.scaley=scale;
    AO.up=up;
    
    Calculate_AOs(tlist,vlist,nfac,nvert,angles,&AO,offset,NULL,nvert,nvert,NULL,AOscale,FT,FTdv,AOds,1);

    //print_matrix(FT,1,2*ntpoints);
    //print_matrix(FTdv,2*ntpoints,3*nvert+2*nao+3);
    write_matrix_file("/tmp/FT.txt",FT,2*ntpoints,1);
  write_matrix_file("/tmp/FTdv.txt",FTdv,2*ntpoints,3*nvert+2*nao+3);
  write_matrix_file("/tmp/AOds.txt",AOds,2*ntpoints,nao);
    free(FT);
    free(FTdv);
    free(AO.datar);
    free(AO.datai);
    free(AO.freqx);
    free(AO.freqy);
}
