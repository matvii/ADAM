function [Z,Fx,Fy]=calc_image_fft(oimage,X,Y,origo,resx,resy)
%keyboard
dx=resx;
dy=resy;
Lx=X(end)-X(1)+dx;
Ly=Y(end)-Y(1)+dy;
Nx=length(X);
Ny=length(Y);
offx=Nx/2-origo(1)+1; 
offy=Ny/2-origo(2)+1;
fnx=[-1/(2*dx):1/Lx:1/(2*dx)-1/Lx];
fny=[-1/(2*dy):1/Ly:1/(2*dy)-1/Ly];
[Fx,Fy]=meshgrid(fnx,fny);
[phasx,phasy]=meshgrid(exp(-2*pi*1i*offx*[0:Nx-1]/Nx),exp(-2*pi*1i*offy*[0:Ny-1]/Ny));

%keyboard
Z=fft2(fftshift(oimage))*dx*dy;
Z=Z/Z(1,1); %Scale
%keyboard
Z=phasx.*phasy.*Z;
Z=fftshift(Z); %zero frequency shifted to the image center
end
