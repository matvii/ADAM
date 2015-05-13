function draw_triangles_3spha_exp(a,nrows)
LMAX=sqrt(length(a)/3)-1;
[THETA,PHI,IFP]=triangulate_sphere(nrows);
shx=a(1:(LMAX+1)^2);
   shy=a(((LMAX+1)^2+1):2*(LMAX+1)^2);
   shz=a((2*(LMAX+1)^2+1):3*(LMAX+1)^2);

B=zeros(length(PHI),(LMAX+1)^2);
for j=0:LMAX
    for k=-j:j
        B(:,j*(j+1)+k+1)=SH(j,k,THETA,PHI)';
    end
end

vlist2=[];
vlist2(:,1)=exp(B*shx').*sin(THETA)'.*cos(PHI)';
vlist2(:,2)=exp(B*(shx'+shy')).*sin(THETA)'.*sin(PHI)';
vlist2(:,3)=exp(B*(shx'+shz')).*cos(THETA)';
trisurf(IFP,vlist2(:,1),vlist2(:,2),vlist2(:,3))
axis equal
end