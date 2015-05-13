function [M,dMdv,dMdA,dMdoff,dMdscale]=Generate_AOFT_Matrix(tlist,vlist,angles,offsets,scale,FTstruct,deriv)
%INPUT:
%tlist: nfacx3 facets
%vlist: nvertx3 vertices
%angles: beta,lambda,omega
%offsets: 1 x 2*nobs vector, offstes
%deriv ==1 is derivatives are to be calculated
%OUTPUT:
%M: vector, data values-model values
%dMdv partial derivatives wrt. vertex coordinates, [d/dx d/dy d/dz]
%dMdA partial derivatives wrt. angles
%dMdoff partial derivatives wrt offsets
%dMdscale partial derivatives wrt scaling term

nE=size(FTstruct.E,2);
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
    [FT{j},~,~,~,~,~]=Calc_AOimage_FT(tlist,vlist,angles,FTstruct.E{j},FTstruct.E0{j},FTstruct.up{j},FTstruct.TIME{j},FTstruct.distance{j},FTstruct.freq{j},offsets(2*j-1:2*j),0);
   if isfield(FTstruct,'psf') && ~isempty(FTstruct.psf{j})
            FT{j}=FTstruct.psf{j}.*FT{j};
   end
        FT{j}=exp(scale)*FT{j}; 
    end
   
  for j=1:nE
   
   M((2*cndata(j)+1):2*cndata(j+1))=[real(FTstruct.data{j})-real(FT{j}) ;imag(FTstruct.data{j})-imag(FT{j})];
   
   end

  
clear FT;
dM=[];
return;
end  


%%%%Calculate derivatives
FT=cell(1,nE);
FTdx=cell(1,nE);
FTdy=cell(1,nE);
FTdz=cell(1,nE);
FTdA=cell(1,nE);

for j=1:nE
   [FT{j},FTdx{j},FTdy{j},FTdz{j},FTdA{j},FTdoff{j}]=Calc_AOimage_FT(tlist,vlist,angles,FTstruct.E{j},FTstruct.E0{j},FTstruct.up{j},FTstruct.TIME{j},FTstruct.distance{j},FTstruct.freq{j},offsets(2*j-1:2*j),1);
  if isfield(FTstruct,'psf') && ~isempty(FTstruct.psf{j})
            FT{j}=FTstruct.psf{j}.*FT{j};
            FTdx{j}=repmat(FTstruct.psf{j},1,nvert).*FTdx{j};
            FTdy{j}=repmat(FTstruct.psf{j},1,nvert).*FTdy{j};
            FTdz{j}=repmat(FTstruct.psf{j},1,nvert).*FTdz{j};
            
            FTdA{j}=repmat(FTstruct.psf{j},1,3).*FTdA{j};
            FTdoff{j}=repmat(FTstruct.psf{j},1,2).*FTdoff{j};
end
        FT{j}=exp(scale)*FT{j};
        FTdx{j}=exp(scale)*FTdx{j};
        FTdy{j}=exp(scale)*FTdy{j};
        FTdz{j}=exp(scale)*FTdz{j};
        FTdA{j}=exp(scale)*FTdA{j};
        FTdoff{j}=exp(scale)*FTdoff{j};
end

%%%%%Form derivative matrix


dMdv=zeros(cndata(end)*2,3*nvert);
dMdA=zeros(cndata(end)*2,3);
dMdoff=zeros(cndata(end)*2,2*nE);
dMdscale=zeros(cndata(end)*2,1);
for j=1:nE
    dMdv((2*cndata(j)+1):2*cndata(j+1),1:nvert)=[real(FTdx{j});imag(FTdx{j})];
    dMdv((2*cndata(j)+1):2*cndata(j+1),nvert+1:2*nvert)=[real(FTdy{j});imag(FTdy{j})];
    dMdv((2*cndata(j)+1):2*cndata(j+1),2*nvert+1:3*nvert)=[real(FTdz{j});imag(FTdz{j})];
     
    dMdA((2*cndata(j)+1):2*cndata(j+1),1:3)=[real(FTdA{j});imag(FTdA{j})];
   M((2*cndata(j)+1):2*cndata(j+1))=[real(FTstruct.data{j})-real(FT{j}) ;imag(FTstruct.data{j})-imag(FT{j})];
    %substract this from the data
end
 



for j=1:nE
    dMdoff((2*cndata(j)+1):2*cndata(j+1),[2*j-1 2*j])=[real(FTdoff{j}(:,1)) real(FTdoff{j}(:,2));imag(FTdoff{j}(:,1)) imag(FTdoff{j}(:,2))];
    dMdscale((2*cndata(j)+1):2*cndata(j+1),1)=[real(FT{j}); imag(FT{j})];
end


clear FT;
 clear FTdx;
 clear FTdy;
 clear FTdz;
end