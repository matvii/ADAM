%See comments in the file AO_Inver_Octantoid_Example
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
scale=arcsec*ones(1,8);
figure;
for k=1:8
FT.distance{k}=dist;
FT.up{k}=[0,0,1]; %Instrument orientation
 filename=strcat('AO_ex/aoim',num2str(k),'.mat'); %Load data
 load(filename);
isiz=size(M);
M=M'; %Transpose, since the data was transposed

origo=[isiz(2)/2+1,isiz(1)/2+1];
subplot(2,4,k); imagesc(M);axis xy;axis equal
X=[-isiz(2)/2*scale(k):scale(k):(isiz(2)/2-1)*scale(k)];
Y=[-isiz(1)/2*scale(k):scale(k):(isiz(1)/2-1)*scale(k)];
[Z,Fx,Fy]=calc_image_fft(M,X,Y,origo,arcsec,arcsec); %Calculate FFT
Iy=Fy>=0;
I0=(Fx==0)&(Fy==0);
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
%%%%Here we init subdivision, vlist is is the set parameters
%Initial shape guess, a triaxial ellipsoid
nrows=5;
[tlist,vlist]=generate_sphere(nrows);
nvert=size(vlist,1);
vlist(:,1)=80*vlist(:,1);
vlist(:,2)=40*vlist(:,2);
vlist(:,3)=40*vlist(:,3);
sdstep=2; %Number of subdivision steps


nft=size(FT.E,2); %Number of data images
offset=zeros(1,2*nft);

ftw=3; %Weight for data
cw=10; %Convex regularization
angw=0.5; %dihedral angle reg
aw=90; %area regularization
scale=0;
alambda=0.1; 
ichisq=Inf;
decreased=1;
chisq=Inf;
weight=1;
%Here we have typical LM optimization
figure;
for j=1:100
    disp('alambda')
    alambda
    if decreased==1
    
  



 [tlistn,vlistn,D]=Sqrt3_Subdiv(tlist,vlist,sdstep);
 
 [Creg,dregdx,dregdy,dregdz]=Convex_Reg(tlistn,vlistn);

 [M,dMdx,dMdy,dMdz,dMdA,dMdoff]=Generate_AOFT_Matrix_mex(tlistn,vlistn,angles,offset,FT,1); %Data
 nvertn=size(vlistn,1);

 [ang1,dangdx,dangdy,dangdz]=dihedral_angle_reg(tlistn,vlistn);
[areg,dSdx,dSdy,dSdz]=area_reg(tlistn,vlistn);

ynew=[ftw*M;cw*Creg ;aw*areg';angw*ang1]; %Data matrix
%Jacobian matrix, chain rule
da=[ftw*dMdx*D ftw*dMdy*D ftw*dMdz*D ftw*dMdA ftw*dMdoff ;-cw*dregdx*D -cw*dregdy*D -cw*dregdz*D zeros(1,3+2*nft); -aw*dSdx*D -aw*dSdy*D -aw*dSdz*D zeros(size(dSdx*D,1),3+2*nft); -angw*dangdx*D -angw*dangdy*D -angw*dangdz*D zeros(1,3+2*nft)];
    ichisq=ynew'*ynew;


    ftfit=ftw^2*M'*M;
   
    end
    disp('Initial chisq:')
    ichisq
    disp('FT:')
    ftfit
   
    disp('Reg terms:')
    disp('convex')
    cw^2*Creg^2
    disp('area')
    aw^2*areg*areg'
    disp('dihedral')
    angw^2*ang1^2
   
  B=da'*da;
 
  delta=(B+alambda*diag(diag(B)))\(da'*(ynew)); %LM
  
  delta1=delta(1:3*nvert); %Shape parameters
  dangles=delta(3*nvert+1:3*nvert+3); %spin axis
  
  delta2=delta(3*nvert+4:end); %Offsets
  offset2=offset+delta2';
  angles2=angles+dangles';
  vlist2=vlist+reshape(delta1,[],3); %New values
  
  %Test if new values are better:
 [tlistn2,vlistn2]=Sqrt3_Subdiv(tlist,vlist2,sdstep);

    M2=Generate_AOFT_Matrix_mex(tlistn2,vlistn2,angles2,offset2,FT,0);
    
    res2=area_reg(tlistn2,vlistn2);
    ang2=dihedral_angle_reg(tlistn2,vlistn2);
 
 Creg2=Convex_Reg(tlistn2,vlistn2);
ynew2=[ftw*M2;cw*Creg2;aw*res2';angw*ang2];
 
  
  
  chisq2=ynew2'*ynew2;
  if chisq2<ichisq %Is the value better than old?
      vlist=vlist2;
      offset=offset2;
      angles=angles2;
      alambda=1/10*alambda;
      disp('chisq decreased')
      decreased=1;
       trisurf(tlistn2,vlistn2(:,1),vlistn2(:,2),vlistn2(:,3)); axis equal
  drawnow;
  else 
      alambda=10*alambda;
      decreased=0;
  end
  if alambda>1e3
      break;
  end
  
 
end
 
