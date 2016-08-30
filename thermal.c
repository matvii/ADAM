#include"utils.h"
void Calculate_Temp(int *tlist, double *vlist,int nfac,int nvert,double* angles,double* E0,double t0,double Gamma, double A,double R,int N,double *Tres);
int main(int argc,char **argv)
{
    char arg;
    int count=1;
    int *tlist;
    double time;
    char *out_file=NULL;
    int VERBOSE=0;
    double *vlist;
    int nfac=0,nvert=0;
    double E0[]={0,0,0};
    double t0=0,Gamma=50,A=0.1,R=1.0;
    int N=1024;
    double *Tres;
    double angles[]={0,0,0};
    double Angles[]={0,0,0};
    char shape_file[100];
    shape_file[0]='\0';
    if(argc==1)
    {
        printf("This program calculates the temperature distribution of an asteroid. Shadowing effects are taken into account, but not self-heating due to reflected radiation.\n");
        printf("Usage:\n");
        printf("thermal -a beta lambda P -E e0 e1 e2 -S shapefile -A albedo -G Gamma -N nFouriercoeffs -t obstime [-o outfile -R Heliodist -V]\n");
        printf("Where beta lambda are the directions of rotation axis, P rotation period in hours, [e0,e1,e2] is the asteroid->Sun vector (in AU), shapefile is the file containing vertex and facet list,\n");
        printf("nFouriercoeffs is the length of FFT used, obstime is the observation time (JD), Gamma is the thermal inertia, Heliodist is the distance (AU) between the asteroid and the Sun (if not set, norm of the vector [e0,e1,e2] is used)\n");
        printf("-V prints values of input parameters\n");
        printf("if outfile is not set, results are printed to stdout\n");
        exit(0);
    }
    while(count<argc)
    {
        if(argv[count][0]=='-')
            arg=argv[count][1];
        else
        {
            fprintf(stderr,"Bad argument %s\n",argv[count]);
            exit(1);
        }
        
        switch(arg)
        {
            case 'a':
                Angles[0]=atof(argv[count+1]);
                Angles[1]=atof(argv[count+2]);
                Angles[2]=atof(argv[count+3]);
                angles[0]=(90-Angles[0])*PI/180;
                angles[1]=Angles[1]*PI/180;
                angles[2]=24.0*2*PI/Angles[2];
                count=count+4;
                break;
            case 'S':
                strcpy(shape_file,argv[count+1]);
                count=count+2;
                break;
            case 'E':
                E0[0]=atof(argv[count+1]);
                E0[1]=atof(argv[count+2]);
                E0[2]=atof(argv[count+3]);
                count=count+4;
                R=NORM(E0);
                E0[0]=E0[0]/R;
                E0[1]=E0[1]/R;
                E0[2]=E0[2]/R;
                break;
            case 'A':
                A=atof(argv[count+1]);
                count=count+2;
                break;
            case 'N':
                N=atoi(argv[count+1]);
                count=count+2;
                break;
            case 'G':
                Gamma=atof(argv[count+1]);
                count=count+2;
                break;
            case 'o':
                out_file=calloc(100,sizeof(char));
                strcpy(out_file,argv[count+1]);
                count=count+2;
                break;
            case 't':
                time=atof(argv[count+1]);
                count=count+2;
                break;
            case 'R':
                R=atof(argv[count+1]);
                count=count+2;
                break;
            case 'V':
                VERBOSE=1;
                count=count+2;
                break;
                
            default:
                printf("Unknown argument: %s\n",argv[count]);
                count=count+1;
        }
    }
    if(shape_file[0]=='\0')
    {
        fprintf(stderr,"Shape file is not set, aborting\n");
        exit(1);
    }
    if(E0[0]==0&&E0[1]==0&&E0[2]==0)
    {
        fprintf(stderr,"E is not set, aborting\n");
        exit(1);
    }
     read_shape(shape_file,&tlist,&vlist,&nfac,&nvert,0);
    if(VERBOSE==1)
    {
        printf("Angles: %f %f %f\n",Angles[0],Angles[1],Angles[2]);
        printf("Input Shape file: %s\n",shape_file);
        printf("Vertices: %d Facets: %d\n",nvert,nfac);
        printf("E: %f %f %f\n",E0[0],E0[1],E0[2]);
        printf("Range: %f Albedo: %f Gamma: %f\n",R,A,Gamma);
        printf("Number of Fourier coefficients: %d\n",N);
    }
    
    
    
     read_shape(shape_file,&tlist,&vlist,&nfac,&nvert,0);
     Tres=calloc(nfac,sizeof(double));
   Calculate_Temp(tlist,vlist,nfac,nvert,angles,E0,time,Gamma,A,R,N,Tres);
   if(out_file!=NULL)
   {
       FILE *fid=fopen(out_file,"w");
       for(int j=0;j<nfac;j++)
           fprintf(fid,"%4.2f ",Tres[j]);
       fprintf(fid,"\n");
       fclose(fid);
   }
   else
   {
       for(int j=0;j<nfac;j++)
           fprintf(stdout,"%4.2f ",Tres[j]);
       fprintf(stdout,"\n");
   }
}
