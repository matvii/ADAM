function [tlist,vlist,dvda]=Octantoid_to_Trimesh(a,nrows)
%Convert octantoid representation a=[ax,by,cz], where
%ax by cz are 1 x (LMAX+1)^2 vectors
%To the standard tlist,vlist representation. Also outputs the derivative
%matrix dx/da
LMAX=sqrt(length(a)/3)-1;
assert(ceil(LMAX)==floor(LMAX),'The length of a must be 3*(LMAX+1)^2')
[THETA,PHI,tlist,ADJ]=triangulate_sphere2(nrows);

SH_table=zeros((LMAX+1)^2,length(PHI));
    for l=0:LMAX
        for k=-l:l
            SH_table(l*(l+1)+k+1,:)=SH(l,k,THETA,PHI);
        end
    end
    shx=a(1:(LMAX+1)^2);
    shy=a((LMAX+1)^2+1:2*(LMAX+1)^2);
    shz=a(2*(LMAX+1)^2+1:3*(LMAX+1)^2);
    B=zeros(length(PHI),(LMAX+1)^2);
    for j=0:LMAX
        for k=-j:j
            B(:,j*(j+1)+k+1)=SH(j,k,THETA,PHI)';
        end
    end
    vlist(:,1)=exp(B*shx').*sin(THETA)'.*cos(PHI)';
    vlist(:,2)=exp(B*(shx+shy)').*sin(THETA)'.*sin(PHI)';
    vlist(:,3)=exp(B*(shx+shz)').*cos(THETA)';
    nvert=size(vlist,1);
  

%keyboard
dxda=[kron(vlist(:,1),ones(1,(LMAX+1)^2)).*B zeros(nvert,(LMAX+1)^2) zeros(nvert,(LMAX+1)^2)];
dyda=[kron(vlist(:,2),ones(1,(LMAX+1)^2)).*B kron(vlist(:,2),ones(1,(LMAX+1)^2)).*B zeros(nvert,(LMAX+1)^2)];
dzda=[kron(vlist(:,3),ones(1,(LMAX+1)^2)).*B zeros(nvert,(LMAX+1)^2) kron(vlist(:,3),ones(1,(LMAX+1)^2)).*B];
dvda=[dxda;dyda;dzda];
end