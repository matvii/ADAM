function [M,dMdx,dMdy,dMdz,dMdA,dMdoff,dMdscale]=Generate_HF_Matrix_mex(tlist,vlist,angles,offsets,scale,FTstruct,Gamma,N,deriv)
%INPUT:
%Generate derivative matrix from 
%tlist,vlist triangular mesh 
%angles: [beta,lambda,omega]
%offsets: offsets on the projection plane, 2*nobs vector, where nobs is number of observations
%FTstuct
%E
%E0
%TIME
%freq nx2
%data n complex values

%obsWavelength: Wavelength used
%Gamma: Thermal inertia
%N number of Fourier coefficients. 
%deriv =1 if partial derivatives are to be calculated, 0 otherwise
%%%%%
%Output:
%M data-model at the frequency points selected by FTstruct.freq
%dMdx,dMdy,dMdz: partial derivatives wrt x,y and z coordinates of vertices
%dMdA, nfreqx3 matrix, partial derivatives wrt angles
%dMdoff partial derivatives wrt offsets
%dMdscale  partial derivatives wrt scaling term
PSF=[];
nE=size(FTstruct.E,2);
nvert=size(vlist,1);

 ndata=zeros(nE+1,1);
for j=1:nE
   ndata(j+1)=size(FTstruct.data{j},1);
end
if isfield(FTstruct,'psf') && ~isempty(FTstruct.psf)
    PSF=cell2mat(FTstruct.psf');
end
cndata=cumsum(ndata);
 
  if length(scale)>1
     ms=1;
    % nscale=nE;
     mscale=zeros(cndata(end),1);
     for j=1:nE
        mscale(cndata(j)+1:cndata(j+1))=exp(scale(j));
     end
  elseif length(scale)==1;
     ms=0;
  else
      ms=-1;
  end
  if deriv==0
   % parfor j=1:nE
    %[FT{j},~,~,~,~,~]=calc_ft_projection_with_deriv2(IFP,vlist,angles,FTstruct.E{j},FTstruct.E0{j},FTstruct.TIME{j},FTstruct.distance{j},FTstruct.freq{j},offsets(2*j-1:2*j),0,EQ);
   %FT{j}=exp(scale(j))*Calc_Heat_Flux(tlist,vlist,angles,FTstruct.E{j},FTstruct.E0{j},FTstruct.up{j},FTstruct.TIME{j},FTstruct.distance{j},FTstruct.freq{j},Gamma,0.1,FTstruct.Hdistance{j},N,FTstruct.WL{j},offsets(2*j-1:2*j),0);
   FT=Calculate_HF(tlist,vlist,angles,FTstruct.E,FTstruct.E0,FTstruct.TIME,FTstruct.freq,FTstruct.up,FTstruct.distance,FTstruct.Hdistance,Gamma,0.1,N,FTstruct.WL,offsets,0);
  % M=[];
   if ~isempty(PSF)
    FT=PSF.*FT;
    end
if ms==0
    FT=exp(scale)*FT;
end
if ms==1
    FT=mscale.*FT;
end
   dat=cell2mat(FTstruct.data(:));
   M=[real(dat)-real(FT);imag(dat)-imag(FT)];

  
clear FT;
%dMdv=[];
dMdx=[];
dMdy=[];
dMdz=[];
dMdA=[];
dMdoff=[];
dMdscale=[];
return;
end  


%%%%Calculate derivatives

   % [FT{j},FTdx{j},FTdy{j},FTdz{j},FTdA{j},FTdoff{j}]=calc_ft_projection_with_deriv2(IFP,vlist,angles,FTstruct.E{j},FTstruct.E0{j},FTstruct.TIME{j},FTstruct.distance{j},FTstruct.freq{j},offsets(2*j-1:2*j),1,EQ);
[FT,FTdx,FTdy,FTdz,FTdA,FTdoff]=Calculate_HF(tlist,vlist,angles,FTstruct.E,FTstruct.E0,FTstruct.TIME,FTstruct.freq,FTstruct.up,FTstruct.distance,FTstruct.Hdistance,Gamma,0.1,N,FTstruct.WL,offsets,1);
%FTdA{j}=zeros(size(FT{j},1),3);
if ~isempty(PSF)
         FT=PSF.*FT;
       %  disp('psf added')
         FTdx=repmat(PSF,1,nvert).*FTdx;
         FTdy=repmat(PSF,1,nvert).*FTdy;
         FTdz=repmat(PSF,1,nvert).*FTdz;
         FTdA=repmat(PSF,1,3).*FTdA;
         FTdoff=repmat(PSF,1,2*nE).*FTdoff;
         
    end
    if ms==0
    FT=exp(scale)*FT;
    FTdx=exp(scale)*FTdx;
     FTdy=exp(scale)*FTdy;
      FTdz=exp(scale)*FTdz;
      FTdA=exp(scale).*FTdA;
      FTdoff=exp(scale).*FTdoff;
      dMdscale=[real(FT);imag(FT)];
    end
    if ms==1
        FT=mscale.*FT;
        FTdx=repmat(mscale,1,nvert).*FTdx;
         FTdy=repmat(mscale,1,nvert).*FTdy;
         FTdz=repmat(mscale,1,nvert).*FTdz;
         FTdA=repmat(mscale,1,3).*FTdA;
         FTdoff=repmat(mscale,1,2*nE).*FTdoff;
         dMdscale1=zeros(cndata(end),nE);
         dMdscale2=zeros(cndata(end),nE);
      % keyboard
        for j=1:nE
           dMdscale1(cndata(j)+1:cndata(j+1),j)=real(FT(cndata(j)+1:cndata(j+1))); 
           dMdscale2(cndata(j)+1:cndata(j+1),j)=imag(FT(cndata(j)+1:cndata(j+1)));
        end
        dMdscale=[dMdscale1;dMdscale2];
        clear dMdscale1 dMdscale2;
    end
    %%%%%Form derivative matrix
dat=cell2mat(FTstruct.data(:));
   M=[real(dat)-real(FT);imag(dat)-imag(FT)];
%M=[real(FT) ;imag(FT)];
dMdx=[real(FTdx);imag(FTdx)];
dMdy=[real(FTdy);imag(FTdy)];
dMdz=[real(FTdz);imag(FTdz)];
dMdA=[real(FTdA);imag(FTdA)];
dMdoff=[real(FTdoff);imag(FTdoff)];
 
     

clear FT;
 clear FTdx;
 clear FTdy;
 clear FTdz;
 clear FTdA;
end