function [flux,dF,dFdA]=Calc_Radiance(tlist,vlist,angles,CE,CE0,t0,Gamma,A,R,N,obsWavelength)
%Calculate heat flux from the object to direction E

%Note that facet blocks another facet if facet centroid is eclipsed.
%Derivative matrix nfac x 3nvert dF=[M1 M2 M3], where M1=dH/dx M2=dH/dy
%M3=dH/dz, ie First row of M1 contains dH1/dx1 dH1/dx2 dH1/dx3 dH1/dx4...
%dH1/dxi, where i is the vertex index
%Probably should consider subfacets TODO

%UNITS HERE ARE Watts/(m^2*Hz*sr)
%Wavelength in m
L=10e-6; %in microns
%%%%%%%%%%%%%%%%%%%%%55%
%L=5.7e-4; %10 microns
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin==11
    L=obsWavelength;
end
    
%L=0.3e-3;

hc=6.62606957e-34; %Planck's
cc=299792458; %speed of light
kc=1.3806488e-23; %Boltzmann

nu=cc/L;

C1=2*hc*nu^3/cc^2;
C2=hc*nu/kc;
beta=angles(1);
lambda=angles(2);
omega=angles(3);
[res,dT,dA]=Calc_Temp(tlist,vlist,angles,CE0,t0,Gamma,A,R,N,1);
%keyboard
%[res,dT,~]=Calc_Poly_Temperature(tlist,vlist,angles,CE0,t0,Gamma,A,R,N,1);
%Rotate
[rotmatrix,~,~,~]=rot_matrix(beta,lambda,omega,t0,0);
E=(rotmatrix*CE')';
[normal,~,~,~,visible]=Calc_Vis(vlist,tlist,E,E); %Facets visible to CE
visible=double(visible);
nfac=size(tlist,1);
nvert=size(vlist,1);
mu=(E*normal');
flux=zeros(1,nfac);
dF=zeros(nfac,3*nvert);
dFdA=zeros(nfac,3);
for j=1:nfac
    if visible(j)==0 %Facet is not visible to the observer
        continue;
    end
  
    vi1=tlist(j,1);
    vi2=tlist(j,2);
    vi3=tlist(j,3);
   
  
    planck=C1/(exp(C2/res(j))-1);
    flux(j)=planck;
    if res(j)==0
        flux(j)=0;
    end
 
    dF(j,vi1)=(C1)*1/(exp(C2/(res(j)))-1)^2*exp(C2/(res(j)))*C2/1/res(j)^2*dT(j,vi1);
    dF(j,vi2)=(C1)*1/(exp(C2/(res(j)))-1)^2*exp(C2/(res(j)))*C2*1/res(j)^2*dT(j,vi2);
    dF(j,vi3)=(C1)*1/(exp(C2/(res(j)))-1)^2*exp(C2/(res(j)))*C2*1/res(j)^2*dT(j,vi3);
    
    dF(j,vi1+nvert)=(C1)*1/(exp(C2/(res(j)))-1)^2*exp(C2/(res(j)))*C2*1/res(j)^2*dT(j,vi1+nvert);
    dF(j,vi2+nvert)=(C1)*1/(exp(C2/(res(j)))-1)^2*exp(C2/(res(j)))*C2*1/res(j)^2*dT(j,vi2+nvert);
    dF(j,vi3+nvert)=(C1)*1/(exp(C2/(res(j)))-1)^2*exp(C2/(res(j)))*C2*1/res(j)^2*dT(j,vi3+nvert);
    
    dF(j,vi1+2*nvert)=(C1)*1/(exp(C2/(res(j)))-1)^2*exp(C2/(res(j)))*C2*1/res(j)^2*dT(j,vi1+2*nvert);
    dF(j,vi2+2*nvert)=(C1)*1/(exp(C2/(res(j)))-1)^2*exp(C2/(res(j)))*C2*1/res(j)^2*dT(j,vi2+2*nvert);
    dF(j,vi3+2*nvert)=(C1)*1/(exp(C2/(res(j)))-1)^2*exp(C2/(res(j)))*C2*1/res(j)^2*dT(j,vi3+2*nvert);
    dFdA(j,:)=(C1)*1/(exp(C2/(res(j)))-1)^2*exp(C2/(res(j)))*C2*1/res(j)^2*dA(j,:);
    if flux(j)==0
        dF(j,:)=0;
        dFdA(j,:)=0;
    end
end

end
    