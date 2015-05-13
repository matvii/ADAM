function [M,dMdx,dMdy,dMdz,dMdA,dMdoff]=Generate_AOFT_Matrix_mex(tlist,vlist,angles,offsets,FTstruct,deriv)
%INPUT:
%tlist: nfacx3 facets
%vlist: nvertx3 vertices
%angles: beta,lambda,omega
%offsets: 1 x 2*nobs vector, offstes
%deriv ==1 is derivatives are to be calculated
%FTstruct.PSF is the FT of the PSF. if PSF is nonempty, it should be the same size as FTstruct.data
%OUTPUT:
%M: vector, data values-model values
%dMdx dMdy dMdz partial derivatives wrt. vertex coordinates
%dMdA partial derivatives wrt. angles
%dMdoff partial derivatives wrt offsets

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
 M=zeros(cndata(end)*2,1);
if deriv==0
    FT=cell(1,nE);
    
    FT=Calculate_AO(tlist,vlist,angles,FTstruct.E,FTstruct.E0,FTstruct.TIME,FTstruct.freq,FTstruct.up,FTstruct.distance,offsets,[],0);
    
    if ~isempty(PSF)
    FT=PSF.*FT;
end
   dat=cell2mat(FTstruct.data(:));
   M=[real(dat)-real(FT);imag(dat)-imag(FT)];

  
clear FT;
dMdv=[];
dMdA=[];
dMdoff=[];
return;
end  


%%%%Calculate derivatives
%%%%Calculate derivatives



    [FT,FTdx,FTdy,FTdz,FTdA,FTdoff]=Calculate_AO(tlist,vlist,angles,FTstruct.E,FTstruct.E0,FTstruct.TIME,FTstruct.freq,FTstruct.up,FTstruct.distance,offsets,[],1);
 
 
    if ~isempty(PSF)
         FT=PSF.*FT;
       %  disp('psf added')
         FTdx=repmat(PSF,1,nvert).*FTdx;
         FTdy=repmat(PSF,1,nvert).*FTdy;
         FTdz=repmat(PSF,1,nvert).*FTdz;
         FTdA=repmat(PSF,1,3).*FTdA;
         FTdoff=repmat(PSF,1,2*nE).*FTdoff;
         
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