%This program demonstrates how to reconstruct a shape from disk-resolved
%observations. Simulated data is used.
high_res=0; %Set this to 1 if you want to use full resolution for reconstruction.
%otherwise effective pixel size used for reconstruction is 2dx, where dx is the data pixel size.
%In any case, there is not much difference, but the computational time is
%considerable less.
%Simulation

%%%%CLEAR VARIABLES
clear FT;
%%%%%%FIRST INIT VARIABLES AND LOAD DATA:
arcsec=0.01; %Data image resolution (arcsecs per pixel)
dist=1.25; %Distance in AU to asteroid
%Set ephemeris information: 
%E-direction of earth as seen from the asteroid
%E0-sun direction
%TIME observation time (in days)
k=1;
FT.E{k}=[1,0.1,0.1];
FT.E0{k}=[1.5,0.1,0.1];
FT.E{k}=FT.E{k}/norm(FT.E{k});
FT.E0{k}=FT.E0{k}/norm(FT.E0{k});
FT.TIME{k}=0;
k=2;
FT.E{k}=[0.1 ,1,0.1];
FT.E0{k}=[0.1,1.5,0.1];
FT.E{k}=FT.E{k}/norm(FT.E{k});
FT.E0{k}=FT.E0{k}/norm(FT.E0{k});
FT.TIME{k}=0.1;
k=3;
FT.E{k}=[0.1 ,0.1,1];
FT.E0{k}=[0.1,0.2,1.5];
FT.E{k}=FT.E{k}/norm(FT.E{k});
FT.E0{k}=FT.E0{k}/norm(FT.E0{k});
FT.TIME{k}=0.2;
k=4;
FT.E{k}=[-1,0.1,0.1];
FT.E0{k}=[-1.5,0.1,0.1];
FT.E{k}=FT.E{k}/norm(FT.E{k});
FT.E0{k}=FT.E0{k}/norm(FT.E0{k});
FT.TIME{k}=0.1;
k=5;
FT.E{k}=[1,0.1,0.1];
FT.E0{k}=[1.5,0.1,0.1];
FT.E{k}=FT.E{k}/norm(FT.E{k});
FT.E0{k}=FT.E0{k}/norm(FT.E0{k});
FT.TIME{k}=0.5;
k=6;
FT.E{k}=[0.1 ,1,0.1];
FT.E0{k}=[0.1,1.5,0.1];
FT.E{k}=FT.E{k}/norm(FT.E{k});
FT.E0{k}=FT.E0{k}/norm(FT.E0{k});
FT.TIME{k}=0.5;
k=7;
FT.E{k}=[0.1 ,0.1,1];
FT.E0{k}=[0.1,0.2,1.5];
FT.E{k}=FT.E{k}/norm(FT.E{k});
FT.E0{k}=FT.E0{k}/norm(FT.E0{k});
FT.TIME{k}=0.7;
k=8;
FT.E{k}=[-1,0.1,0.1];
FT.E0{k}=[-1.5,0.1,0.1];
FT.E{k}=FT.E{k}/norm(FT.E{k});
FT.E0{k}=FT.E0{k}/norm(FT.E0{k});
FT.TIME{k}=0.2;
scale=arcsec*ones(1,8); %pixel scale
figure;
for k=1:8
FT.distance{k}=dist;
FT.up{k}=[0,0,1]; %Camera orientation
 filename=strcat('AO_ex/aoim',num2str(k),'.mat');
 
 load(filename); %Load data
isiz=size(M); 
M=M'; %Original simulated data was transposed, transpose back

%M(M>0)=1; %Geometric scattering

subplot(2,4,k); imagesc(M);axis xy;axis equal
% pause
origo=[isiz(2)/2+1,isiz(1)/2+1];
X=[-isiz(2)/2*scale(k):scale(k):(isiz(2)/2-1)*scale(k)]; %image coordinates
Y=[-isiz(1)/2*scale(k):scale(k):(isiz(1)/2-1)*scale(k)];
[Z,Fx,Fy]=calc_image_fft(M,X,Y,origo,arcsec,arcsec); %Calculate 2D FFT and scale wrt. total brightness
Iy=Fy>=0; %Discard half of frequencies (not needed)
I0=(Fx==0) & (Fy==0); %
%NB: HERE WE IGNORE 0-frequency, this special case not implemented,
%Because not needed. Will cause NaN is used.
%%%Use all the data?
if high_res==0
    Ires=(abs(Fy))<=1/(4*scale(k))&(abs(Fx))<=1/(4*scale(k));
else
    Ires=true(size(Fy));
end
Index=Iy&~I0&Ires; %Note we assume here that image the is a square. More care is required otherwise.
FT.data{k}=Z(Index); 
FT.freq{k}=[Fx(Index) Fy(Index)];
end

%%%%%%%%%
%OPTIMIZATION PART
%%%%%%%%%%%%

%Here is spin vector direction, omega is rad/day
%This guess should be a good one, as this program is not meant to fit spin
%vector, only do some small refinitions
omega=24*2*pi*1/5.3852808;
beta=(90-21)*pi/180;
lambda=72*pi/180;
angles=[beta,lambda,omega];
%Initialize the shape
nrows=10; %Corresponds to the number of facets
[THETA,PHI,IFP,ADJ]=triangulate_sphere2(nrows);
LMAX=6; %Degree of spherical harmonics

a=1.05;
b=1;
c=0.95;
%%%%%%%

%%%%%%
%SCALING!!!
%NOTE: correct scale must chosen carefully to obtain the best convergence,
%it is better to choose too small a diameter than too large
escl=100;
a1=a*escl;
b1=b*escl;
c1=c*escl;

x1=a1*sin(THETA).*cos(PHI);
y1=b1*sin(THETA).*sin(PHI);
z1=c1*cos(THETA);
%Calculate Coefficient matrix
B=zeros(length(PHI),(LMAX+1)^2);
for j=0:LMAX
    for k=-j:j
        B(:,j*(j+1)+k+1)=SH(j,k,THETA,PHI)';
    end
end
r=sqrt(x1.^2+y1.^2+z1.^2);

%There are better ways...
    xc=B\log(a1*ones(1,size(x1,2)))';
    yc=B\(log(b1*ones(1,size(x1,2)))'-B*xc);
     zc=B\(log(c1*ones(1,size(x1,2)))'-B*xc);
    a=[xc' yc' zc' beta lambda omega];
    


    
    
 nft=size(FT.E,2); 
mask=false(1,length(a)+2*nft); %can be used to freeze some parameters


offset=zeros(1,2*nft); %Image offset wrt projected shape (two coordinates)
ftw=2; %Data weight
cw=5; %Convex weight
ow=5; %Octantoid weight
scale=0; %scale, not used here
alambda=0.001;
ichisq=Inf
decreased=1;
chisq=Inf;
 %LM Optimization
 figure;
for j=1:200
    disp('alambda')
    alambda
    if decreased==1
   [tlist,vlist,dvda]=Octantoid_to_Trimesh(a(1:3*(LMAX+1)^2),nrows); %Convert to the usual trimesh
   nvert=size(vlist,1);
   angles=a(3*(LMAX+1)^2+1:3*(LMAX+1)^2+3);
    [oreg,doreg]=Octantoid_Reg(a,LMAX); 
   [creg,dcregdx,dcregdy,dcregdz]=Convex_Reg(tlist,vlist);
   dcreg=[dcregdx dcregdy dcregdz]*dvda;  %chain rule
    
%[M,dMdv,dMdA,dMdoff,dMscale]=Generate_AOFT_Matrix(tlist,vlist,angles,offset,scale,FT,1); %Data
  [M,dMdx,dMdy,dMdz,dMdA,dMdoff]=Generate_AOFT_Matrix_mex(tlist,vlist,angles,offset,FT,1); 
  dMdv=[dMdx dMdy dMdz];
    da=[ftw*dMdv*dvda ftw*dMdA ftw*dMdoff ;-ow*doreg zeros(1,3+2*nft);-cw*dcreg zeros(1,3+2*nft)]; %Note the minus sign
    ynew=[ftw*M;ow*oreg;cw*creg];
    da(:,mask)=[]; 
  
    end
    disp('Initial chisq:')
    ichisq=ynew'*ynew
    disp('Reg terms:')
    ynew(end-1).^2
    ynew(end).^2
  B=da'*da;
  delta=zeros(size(mask));
  delta(~mask)=(B+alambda*diag(diag(B)))\(da'*(ynew)); %LM
 
  %Take a step
 delta1=delta(1:3*(LMAX+1)^2+3); 
  delta2=delta(3*(LMAX+1)^2+4:end);
  a2=a+delta1;
  offset2=offset+delta2;
  %Test if we get better fit
 [tlist2,vlist2,dvda]=Octantoid_to_Trimesh(a2(1:3*(LMAX+1)^2),nrows); %Convert to the usual trimesh
   angles2=a2(3*(LMAX+1)^2+1:3*(LMAX+1)^2+3);
    oreg2=Octantoid_Reg(a2,LMAX);
   creg2=Convex_Reg(tlist2,vlist2);

  % M2=Generate_AOFT_Matrix(tlist2,vlist2,angles2,offset2,scale,FT,0); %Data
    M2=Generate_AOFT_Matrix_mex(tlist2,vlist2,angles2,offset2,FT,0); %Data
    ynew2=[ftw*M2;ow*oreg2;cw*creg2];
  
  
  chisq2=ynew2'*ynew2;
  if chisq2<ichisq
      %Ok, new parameters are better
      a=a2;
      offset=offset2;
      alambda=0.1*alambda;
      disp('chisq decreased')
      decreased=1;
      draw_triangles_3spha_exp(a(1:3*(LMAX+1)^2),nrows); %draw shape
      drawnow;
  else 
      alambda=10*alambda;
      decreased=0;
  end
  if alambda>1e3
      break;
  end
end
