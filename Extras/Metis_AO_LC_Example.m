LC=read_lcurve_struct('Metis_lc',1); %Read the lcurve data
metis_ao_prep; %Prepare AO data
T0=2433222.0; %zero time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:size(LC.TIME,2)
    LC.TIME{j}=LC.TIME{j}-T0;
end
for j=1:size(FT.TIME,2)
FT.TIME{j}=FT.TIME{j}-T0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


beta=(90-24)*pi/180;
    lambda=185*pi/180;
    omega=24*2*pi*1/5.079177;
nrows=10; %How many facets
[THETA,PHI,IFP,ADJ]=triangulate_sphere2(nrows);
LMAX=5; %Maximum spherical harmonics degree

a=1.05; %Initial sphere
b=1;
c=0.95;
%%%%%%%

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

%This is a mess, just initialize parameters to some values
    xc=B\log(a1*ones(1,size(x1,2)))';
    yc=B\(log(b1*ones(1,size(x1,2)))'-B*xc);
     zc=B\(log(c1*ones(1,size(x1,2)))'-B*xc);
    a=[xc' yc' zc' beta lambda omega];
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %Weights corresponding to data and regularization terms
  nft=size(FT.E,2); %Number of observations
offset=zeros(1,2*nft); %Offset information for each image
scale=0; %scaling term, scaling is exp(scale)
ftw=2; %Weight of FT data
cw=3; %Convex regularization weight
ow=6; %Octantoid regularization. Increase this is shape looks funny
lcw=2; %Weight of lcurves

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



mask=false(1,length(a)+2*nft);

alambda=1;
ichisq=Inf
decreased=1;
chisq=Inf;



figure;
%Just the usual LM optimization loop
for j=1:100
    disp('alambda')
    alambda
    if decreased==1
        [tlist,vlist,dvda]=Octantoid_to_Trimesh(a(1:3*(LMAX+1)^2),nrows); %Convert to the usual trimesh
        nvert=size(vlist,1);
        angles=a(3*(LMAX+1)^2+[1:3]);
        [ldata,dlcdx,dlcdy,dlcdz,dlcdA]=Generate_LC_Matrix(tlist,vlist,angles,LC) ; %lcurve part
        [oreg,doreg]=Octantoid_Reg(a,LMAX); %Octantoid regularization
        [creg,dcregdx,dcregdy,dcregdz]=Convex_Reg(tlist,vlist); %Convex regularization
        [M,dMdx,dMdy,dMdz,dMdA,dMdoff]=Generate_AOFT_Matrix_mex(tlist,vlist,angles,offset,FT,1);
        %[M,dMdv,dMdA,dMdoff]=Generate_AOFT_Matrix(tlist,vlist,angles,offset,scale,FT,1);
        
   %Derivative matrix for LM algorithm, note minus sign in the reg terms
        da=[lcw*[dlcdx dlcdy dlcdz]*dvda lcw*dlcdA zeros(size(dlcdx,1),2*nft); 
            ftw*[dMdx dMdy dMdz]*dvda ftw*dMdA ftw*dMdoff ;
            -ow*doreg zeros(1,3+2*nft); 
            -cw*[dcregdx dcregdy dcregdz]*dvda zeros(1,3+2*nft)];
   
  
        ynew=[lcw*ldata;ftw*M;ow*oreg; cw*creg]; 
    end
    disp('Initial chisq:')
    ichisq=ynew'*ynew
    disp('chisq FT:');
    ftw^2*M'*M
    disp('chisq LC:');
    lcw^2*ldata'*ldata
    disp('Reg terms:')
    ynew(end-1).^2
    ynew(end).^2
  B=da'*da;
  delta=zeros(size(mask));
  delta(~mask)=(B+alambda*diag(diag(B)))\(da'*(ynew)); %Invert matrix
  
  delta1=delta(1:3*(LMAX+1)^2+3); %Shape parameters
  delta2=delta(3*(LMAX+1)^2+4:end); %offsets
 % delta3=delta(end); %scaling not used here
  a2=a+delta1;
  offset2=offset+delta2;
  %scale2=scale+delta3;
  angles2=a2(3*(LMAX+1)^2+[1 2 3]);
  [tlist2,vlist2]=Octantoid_to_Trimesh(a2(1:3*(LMAX+1)^2),nrows); %Convert to the usual trimesh
  [ldata2,~]=Generate_LC_Matrix(tlist2,vlist2,angles2,LC); 
     oreg2=Octantoid_Reg(a2,LMAX);
   creg2=Convex_Reg(tlist2,vlist2);
 
    M2=Generate_AOFT_Matrix_mex(tlist2,vlist2,angles2,offset2,FT,0);
%M2=Generate_AOFT_Matrix(tlist2,vlist2,angles2,offset2,scale,FT,0);
  ynew2=[lcw*ldata2;ftw*M2;ow*oreg2; cw*creg2];
  
  chisq2=ynew2'*ynew2;
  if chisq2<ichisq %Decreased?
      a=a2;
      offset=offset2;
      alambda=0.1*alambda;
    %scale=scale2;
      disp('chisq decreased')
      decreased=1;
      draw_triangles_3spha_exp(a(1:3*(LMAX+1)^2),nrows); %Draw object
      drawnow;
  else 
      alambda=10*alambda;
      decreased=0;
  end
  if alambda>1e3
      break;
  end
end
