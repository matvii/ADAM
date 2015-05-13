function [M,dMdx,dMdy,dMdz,dMdA,dMdoff,dMdscale,dMdMs]=Generate_AOFT_Matrix_mex_scaled_M(tlist,vlist,angles,offsets,scale,Modelscale,FTstruct,deriv)
%INPUT:
%tlist: nfacx3 facets
%vlist: nvertx3 vertices
%angles: beta,lambda,omega
%scale: either empty, scalar or vector
%offsets: 1 x 2*nobs vector, offstes
%deriv ==1 is derivatives are to be calculated
%Modelscale = model scale
%FTstruct.PSF is the FT of the PSF. if PSF is nonempty, it should be the same size as FTstruct.data
%OUTPUT:
%M: vector, data values-model values
%dMdx dMdy dMdz partial derivatives wrt. vertex coordinates
%dMdA partial derivatives wrt. angles
%dMdoff partial derivatives wrt offsets
%keyboard
vlist=exp(Modelscale)*vlist;
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
%    
%    % nscale=1;
% end
if deriv==0
    FT=Calculate_AO(tlist,vlist,angles,FTstruct.E,FTstruct.E0,FTstruct.TIME,FTstruct.freq,FTstruct.up,FTstruct.distance,offsets,[],0);
    
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
%dMdA=[];
%dMdoff=[];
%dMscale=[];
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
dMdx=[real(FTdx);imag(FTdx)]*exp(Modelscale);
dMdy=[real(FTdy);imag(FTdy)]*exp(Modelscale);
dMdz=[real(FTdz);imag(FTdz)]*exp(Modelscale);
dMdA=[real(FTdA);imag(FTdA)];
dMdoff=[real(FTdoff);imag(FTdoff)];
%keyboard 
dMdMsc=FTdx*vlist(:,1)+FTdy*vlist(:,2)+FTdz*vlist(:,3);
dMdMs=[real(dMdMsc);imag(dMdMsc)];

clear FT;
 clear FTdx;
 clear FTdy;
 clear FTdz;
 clear FTdA;
end