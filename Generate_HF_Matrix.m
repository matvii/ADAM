function [M,dMdx,dMdy,dMdz,dMdA,dMdoff,dMdscale]=Generate_HF_Matrix(tlist,vlist,angles,offsets,scale,FTstruct,Gamma,N,deriv)
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
nE=size(FTstruct.E,2);
if length(scale)>1
    ms=true;
    nscale=nE;
else
    ms=false;
    scale=scale*ones(1,nE);
    nscale=nE;
end
nvert=size(vlist,1);

ndata=zeros(nE+1,1);
for j=1:nE
   ndata(j+1)=size(FTstruct.data{j},1);
end
cndata=cumsum(ndata);
 M=zeros(cndata(end)*2,1);
    
if deriv==0
    FT=cell(1,nE);
    parfor j=1:nE
    %[FT{j},~,~,~,~,~]=calc_ft_projection_with_deriv2(IFP,vlist,angles,FTstruct.E{j},FTstruct.E0{j},FTstruct.TIME{j},FTstruct.distance{j},FTstruct.freq{j},offsets(2*j-1:2*j),0,EQ);
    FT{j}=exp(scale(j))*Calc_Heat_Flux(tlist,vlist,angles,FTstruct.E{j},FTstruct.E0{j},FTstruct.up{j},FTstruct.TIME{j},FTstruct.distance{j},FTstruct.freq{j},Gamma,0.1,FTstruct.Hdistance{j},N,FTstruct.WL{j},offsets(2*j-1:2*j),0);
    end
    M=[];
 for j=1:nE
   
   M((2*cndata(j)+1):2*cndata(j+1),1)=[real(FTstruct.data{j})-real(FT{j}) ;imag(FTstruct.data{j})-imag(FT{j})];
   
   end

dMdx=[];
dMdy=[];
dMdz=[];
dMdA=[];
dMdoff=[];
dMdscale=[];
return;
end  


%%%%Calculate derivatives
FT=cell(1,nE);
FTdx=cell(1,nE);
FTdy=cell(1,nE);
FTdz=cell(1,nE);
FTdA=cell(1,nE);
freqsize=zeros(1,nE);
%keyboard
parfor j=1:nE
   % [FT{j},FTdx{j},FTdy{j},FTdz{j},FTdA{j},FTdoff{j}]=calc_ft_projection_with_deriv2(IFP,vlist,angles,FTstruct.E{j},FTstruct.E0{j},FTstruct.TIME{j},FTstruct.distance{j},FTstruct.freq{j},offsets(2*j-1:2*j),1,EQ);
[FT{j},FTdx{j},FTdy{j},FTdz{j},FTdA{j},FTdoff{j}]=Calc_Heat_Flux(tlist,vlist,angles,FTstruct.E{j},FTstruct.E0{j},FTstruct.up{j},FTstruct.TIME{j},FTstruct.distance{j},FTstruct.freq{j},Gamma,0.1,FTstruct.Hdistance{j},N,FTstruct.WL{j},offsets(2*j-1:2*j),1);
%FTdA{j}=zeros(size(FT{j},1),3);
FT{j}=exp(scale(j))*FT{j};
FTdoff{j}=exp(scale(j))*FTdoff{j};
FTdx{j}=exp(scale(j))*FTdx{j};
 FTdy{j}=exp(scale(j))*FTdy{j};
FTdz{j}=exp(scale(j))*FTdz{j};
 FTdA{j}=exp(scale(j))*FTdA{j};
end
%%%%%Form derivative matrix
%keyboard
dMdx=zeros(cndata(end)*2,nvert);
dMdy=zeros(cndata(end)*2,nvert);
dMdz=zeros(cndata(end)*2,nvert);
dMdA=zeros(cndata(end)*2,3);
dMdoff=zeros(cndata(end)*2,2*nE);
dMdscale=zeros(cndata(end)*2,nscale);
for j=1:nE
    dMdx((2*cndata(j)+1):2*cndata(j+1),1:nvert)=[real(FTdx{j});imag(FTdx{j})];
    dMdy((2*cndata(j)+1):2*cndata(j+1),1:nvert)=[real(FTdy{j});imag(FTdy{j})];
    dMdz((2*cndata(j)+1):2*cndata(j+1),1:nvert)=[real(FTdz{j});imag(FTdz{j})];
     
    dMdA((2*cndata(j)+1):2*cndata(j+1),1:3)=[real(FTdA{j});imag(FTdA{j})];
   M((2*cndata(j)+1):2*cndata(j+1))=[real(FTstruct.data{j})-real(FT{j}) ;imag(FTstruct.data{j})-imag(FT{j})];
    %substract this from the data
end
if ms==true
for j=1:nE
    dMdoff((2*cndata(j)+1):2*cndata(j+1),[2*j-1 2*j])=[real(FTdoff{j}(:,1)) real(FTdoff{j}(:,2));imag(FTdoff{j}(:,1)) imag(FTdoff{j}(:,2))];
    dMdscale((2*cndata(j)+1):2*cndata(j+1),j)=[real(FT{j}); imag(FT{j})];
end
else
   for j=1:nE
    dMdoff((2*cndata(j)+1):2*cndata(j+1),[2*j-1 2*j])=[real(FTdoff{j}(:,1)) real(FTdoff{j}(:,2));imag(FTdoff{j}(:,1)) imag(FTdoff{j}(:,2))];
    dMdscale((2*cndata(j)+1):2*cndata(j+1),1)=[real(FT{j}); imag(FT{j})];
    end 
end
end