function res=Calc_image_FT(M,X,Y,dx,dy,Fx,Fy)
%image M, with coordinates X Y
%Coordinate of image point M(i,j) at (i,j)
%is (Y(i,j),X(i,j))  Using MATLAB convention
%Calculated at frequencies Fx Fy
%Note that res(i,j) corresponds to frequency Fy(i,j) Fx(i,j) 
Ny=size(X,1);
Nx=size(X,2);
res=zeros(size(Fx));
for j=1:Nx*Ny
    res=res+M(j)*exp(-2*pi*1i*(X(j)*Fx+Y(j)*Fy)); %Each pixel corresponds to delta function
end
 res=res.*sinc(dx*Fx).*sinc(dy*Fy)*dx*dy;
 end