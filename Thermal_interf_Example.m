
%Setup observation Geometry
clear FT;
R=1;
A=0.1;
Gamma=100;
N=32;

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
FT.TIME{k}=-0.1;


k=3;
FT.E{k}=[0.1 ,0.1,1];
FT.E0{k}=[0.1,0.2,1.5];
FT.E{k}=FT.E{k}/norm(FT.E{k});
FT.E0{k}=FT.E0{k}/norm(FT.E0{k});
FT.TIME{k}=-0.2;

k=4;
FT.E{k}=[-1,0.1,0.1];
FT.E0{k}=[-1.5,0.1,0.1];
FT.E{k}=FT.E{k}/norm(FT.E{k});
FT.E0{k}=FT.E0{k}/norm(FT.E0{k});
FT.TIME{k}=-0.03;


k=5;
FT.E{k}=[1,0.1,0.1];
FT.E0{k}=[1.5,0.1,0.1];
FT.E{k}=FT.E{k}/norm(FT.E{k});
FT.E0{k}=FT.E0{k}/norm(FT.E0{k});
FT.TIME{k}=-0.15;


k=6;
FT.E{k}=[-1 ,0.5,0.5];
FT.E0{k}=[-1,1,0.5];
FT.E{k}=FT.E{k}/norm(FT.E{k});
FT.E0{k}=FT.E0{k}/norm(FT.E0{k});
FT.TIME{k}=-0.01;

k=7;
FT.E{k}=[0.1 ,0.1,1];
FT.E0{k}=[0.1,0.2,1.5];
FT.E{k}=FT.E{k}/norm(FT.E{k});
FT.E0{k}=FT.E0{k}/norm(FT.E0{k});
FT.TIME{k}=-0.17;

k=8;
FT.E{k}=[-1,0.1,0.1];
FT.E0{k}=[-1.5,0.1,0.1];
FT.E{k}=FT.E{k}/norm(FT.E{k});
FT.E0{k}=FT.E0{k}/norm(FT.E0{k});
FT.TIME{k}=-0.09;

k=9;
FT.E{k}=[1,0.1,0.1];
FT.E0{k}=[1.5,0.1,0.1];
FT.E{k}=FT.E{k}/norm(FT.E{k});
FT.E0{k}=FT.E0{k}/norm(FT.E0{k});
FT.TIME{k}=0.05;

k=10;
FT.E{k}=[0.1 ,1,0.1];
FT.E0{k}=[0.1,1.5,0.1];
FT.E{k}=FT.E{k}/norm(FT.E{k});
FT.E0{k}=FT.E0{k}/norm(FT.E0{k});
FT.TIME{k}=0.07;

k=11;
FT.E{k}=[0.1 ,1,0.1];
FT.E0{k}=[0.1,1.5,0.1];
FT.E{k}=FT.E{k}/norm(FT.E{k});
FT.E0{k}=FT.E0{k}/norm(FT.E0{k});
FT.TIME{k}=0.09;

arcsec=0.005;
nsiz=size(FT.E,2);
scale=arcsec*ones(1,nsiz);


F=350e9; %Observation frequency
L=299792458/F;
for j=1:nsiz
    FT.distance{j}=1.5;
    FT.Hdistance{j}=1;
    FT.up{j}=[0,0,1];
    FT.WL{j}=L;
end
   
%Read data
for k=1:11
      filename=strcat('Thermal_ex_data/','kleo_uv_',num2str(k),'.dat');
     M=dlmread(filename);
     
     FT.data{k}=M(4,:)'+1i*(M(5,:)');
     FT.freq{k}=[M(1,:)' -M(2,:)']; %Note the rotation here, as the original images were rotated
     
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Draw a sample image
figure;
k=11;
[X,Y]=meshgrid([-0.2:5e-3:0.2-5e-3],[-0.2:5e-3:0.2-5e-3]);
fx=FT.freq{k}(:,1);
fy=FT.freq{k}(:,2);
Z=FT.data{k};
res=zeros(size(X));
for j=1:length(Z);
res=res+Z(j)*exp(2*pi*1i*(fx(j)*X+fy(j)*Y));
end
imagesc([-0.2:5e-3:0.2-5e-3],[-0.2:5e-3:0.2-5e-3],real(res)); axis xy; axis off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Reconstruct shape
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    beta=69*pi/180;
    lambda=72*pi/180;
    omega=24*2*pi*1/5.3852808;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%Intial shape, Ellipsoid for octantoids,
%Subdivision surfaces can use arbitrary
nrows=10;
[THETA,PHI,IFP,ADJ]=triangulate_sphere2(nrows);
LMAX=4;



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

%%%%%This is for octantoids
%Calculate Coefficient matrix
B=zeros(length(PHI),(LMAX+1)^2);
for j=0:LMAX
    for k=-j:j
        B(:,j*(j+1)+k+1)=SH(j,k,THETA,PHI)';
    end
end
r=sqrt(x1.^2+y1.^2+z1.^2);


    xc=B\log(a1*ones(1,size(x1,2)))';
    yc=B\(log(b1*ones(1,size(x1,2)))'-B*xc);
     zc=B\(log(c1*ones(1,size(x1,2)))'-B*xc);
    a=[xc' yc' zc' beta lambda omega]; %Parameter vector
   

%Weights for data and regularization
ftw=0.15;
cw=1;
ow=1;
%%%%%%%%%%%%%%    
nft=size(FT.E,2);
offset=zeros(1,2*nft); %Offsets
mask=false(1,length(a)+2*nft+nft);
%mask(3*(LMAX+1)^2+[ 1 2 3])=true; %Do not fit rotation angles
scale=zeros(1,nft); %Exponent
alambda=0.001;
ichisq=Inf
decreased=1;
chisq=Inf;
weight=1;
N=512; %Number of FFT samples
Gamma=100; %Thermal inertia, it is useless to fit this, since the shape depends only weakly on thermal inertia
%Just the usual LM optimization loop
for j=1:100
    disp('alambda')
    alambda
    if decreased==1
   [tlist,vlist,dvda]=Octantoid_to_Trimesh(a(1:3*(LMAX+1)^2),nrows); %Convert to the usual trimesh
   
    angles=a(3*(LMAX+1)^2+[1:3]);
   [oreg,doreg]=Octantoid_Reg(a,LMAX); %Octantoid regularization
    [creg,dcregdx,dcregdy,dcregdz]=Convex_Reg(tlist,vlist); %Convex regularization
[M,dMdx,dMdy,dMdz,dMdA,dMdoff,dMdscale]=Generate_HF_Matrix_mex(tlist,vlist,angles,offset,scale,FT,Gamma,N,1);

  dMda=[dMdx dMdy dMdz]*dvda; %Chain rule
    dcregda=[dcregdx dcregdy dcregdz]*dvda;
    da=[ftw*dMda ftw*dMdA ftw*dMdoff ftw*dMdscale ;-ow*doreg zeros(1,3+2*nft+nft); -cw*dcregda zeros(1,3+2*nft+nft)];
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
  delta(~mask)=(B+alambda*diag(diag(B)))\(da'*(ynew));
  
  delta1=delta(1:3*(LMAX+1)^2+3)';
  
  delta2=delta(3*(LMAX+1)^2+4:end-nft)';
  delta3=delta(end-nft+1:end)';
  a2=a+delta1';
  
  offset2=offset+delta2';
  scale2=scale+delta3';
 [tlist2,vlist2,~]=Octantoid_to_Trimesh(a2(1:3*(LMAX+1)^2),nrows); %Convert to the usual trimesh
    nvert=size(vlist2,1);
    angles2=a2(3*(LMAX+1)^2+[1:3]);
   oreg2=Octantoid_Reg(a2,LMAX); %Octantoid regularization
    creg2=Convex_Reg(tlist2,vlist2); 
 

    [M2,~]=Generate_HF_Matrix_mex(tlist2,vlist2,angles2,offset2,scale2,FT,Gamma,N,0);
    
    ynew2=[ftw*M2;ow*oreg2;cw*creg2];
  
  
  chisq2=ynew2'*ynew2;
  if chisq2<ichisq
      a=a2;
      offset=offset2;
      scale=scale2;
      alambda=0.1*alambda;
      disp('chisq decreased')
      decreased=1;
      draw_triangles_3spha_exp(a(1:3*(LMAX+1)^2),10);
      drawnow;
  else 
      alambda=10*alambda;
      decreased=0;
  end
  if alambda>1e6
      break;
  end
end
