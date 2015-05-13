function [M,dMdx,dMdy,dMdz,dMdA,dMdoff,dMdscale]=Generate_RD_Matrix(tlist,vlist,angles,offsets,scale,FTstruct,rexp,deriv)
%Scale every image separately
%Generate derivative matrix from 
%FTstuct:
%E
%E0
%TIME
%freq nx2
%data n complex values
%radarfreq radar frequency
%offsets: 1x2nobs vector
%scale: 1xnobs vector
%rexp: radar parameter, exponent of cosine law. Note that
%exponent=exp(rexp) (to guarantee nonnegativity)
%Radar image orientation:
%Frequency first coordinate, increases from top to bottom
%Range second coordinate, increases from left to right
nE=size(FTstruct.E,2);

if deriv==0
   
    FT=Calculate_RD(tlist,vlist,angles,FTstruct.E,FTstruct.TIME,FTstruct.freq,FTstruct.radarfreq,offsets,scale,rexp,0);
   
   dat=cell2mat(FTstruct.data(:));
   M=[real(dat)-real(FT);imag(dat)-imag(FT)];
   
    
dMdx=[];
dMdy=[];
dMdz=[];
dMdA=[];
dMdoff=[];
dMdscale=[]
return;
end



%%%%Calculate derivatives

%keyboard
freqsize=zeros(1,nE);
for j=1:nE
     freqsize(j)=size(FTstruct.freq{j},1);
end
[FT,FTdx,FTdy,FTdz,FTdA,FTdoff]=Calculate_RD(tlist,vlist,angles,FTstruct.E,FTstruct.TIME,FTstruct.freq,FTstruct.radarfreq,offsets,scale,rexp,1);



%%%%%Form derivative matrix

dMdx=[real(FTdx);imag(FTdx)];
dMdy=[real(FTdy);imag(FTdy)];
dMdz=[real(FTdz);imag(FTdz)];
dMdA=[real(FTdA);imag(FTdA)];
dMdoff=[real(FTdoff);imag(FTdoff)];
dat=cell2mat(FTstruct.data(:));
   M=[real(dat)-real(FT);imag(dat)-imag(FT)];
   
%Offset and scale
dMscalere=zeros(sum(freqsize),1);
dMscaleim=zeros(sum(freqsize),1);

offindex=1;
for j=1:nE
  
    dMscalere(offindex:sum(freqsize(1:j)),j)=real(FT(offindex:sum(freqsize(1:j))));
    dMscaleim(offindex:sum(freqsize(1:j)),j)=imag(FT(offindex:sum(freqsize(1:j))));
    
    offindex=offindex+freqsize(j);
end    
   

dMdscale=[dMscalere;dMscaleim];



end