%Calculate Fourier transforms of images
%Initialize cell variable FT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FT:
%FT.E - Asteroid->obs. instrument, unit vector
%FT.E0 -Asteroid->Sun, unit vector
%FT.TIME - LT-corrected TIME
%FT.freq - nfreq x 2 vector, spatial frequencies used
%FT.data -nfreq x 1 vector, complex values corresponding to the spatial
%frequencies
%FT.pdf - nfreq x 1 vector, Fourier transform of the system PSF evaluated at the
%spatial frequencies (can be empty)
%FT.scale - pixel scale
%FT.up - camera orientation, unit vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear FT;
filelist=dir('Metis/*.fits');
TIME=[];
im=[];
pixscale=0.009942; %arcsecs per pixel
t0=2444771.793820;
scale=pixscale*ones(1,size(filelist,1));
for j=1:size(filelist,1)
    filename=strcat('Metis/',filelist(j).name);
    Info=fitsinfo(filename);
    isize=size(Info.PrimaryData.Keywords,1);
    for k=1:isize
 if strcmp(Info.PrimaryData.Keywords{k,1},'MJD-OBS')
    time= Info.PrimaryData.Keywords{k,2}+2400000.5;
    break;
 end
    end
      
   
TIME(j)=time;
nobs=size(TIME,2);
E=zeros(nobs,3);
E0=zeros(nobs,3);
Filename{j}=filename;
end
TIMElt=zeros(1,nobs);

M=dlmread('Metis/ephm.dat'); %Ephemeris information from the file
%%%%%
%File structure:
%MJD-OBS E0 E LT_corrected_time
%%%%%%
for j=1:nobs
    [I,J]=min(abs(M(:,1)-TIME(j)));
    if I>1e-5
        error('MJD-OBS time not found in ephm.dat')
    end
    E0(j,:)=M(J,2:4);
    E(j,:)=M(J,5:7);
    TIMElt(j)=M(J,8);
end
ep=23.4375*pi/180; %equatorial-ecliptic tilt
Req=[1 0 0;0 cos(ep) sin(ep);0 -sin(ep) cos(ep)];
 for k=1:nobs
FT.E{k}=E(k,:);
dist=norm(E(k,:));
FT.E{k}=FT.E{k}/dist;
FT.scale{k}=scale(k);
FT.E0{k}=E0(k,:)/norm(E0(k,:));
FT.TIME{k}=TIMElt(k);

FT.distance{k}=dist;
%up vector
up=(Req*[0,0,1]')'; %Observations are in EQ coordinates, we use ecliptic

M=fitsread(Filename{k});
M=M(51:200,51:200);  %Take only relevant part of the image, but take enough surrounding so that the FT by FFT is affected by aliasing
im{k}=M;
npix=150;
FT.up{k}=up;
isiz=size(M);

origo=[isiz(2)/2+1,isiz(1)/2+1]; %Assumed image origin
X=[-isiz(2)/2*scale(k):scale(k):(isiz(2)/2-1)*scale(k)]; %Or is it otherway round? CHECK THIS!!! better use square images
Y=[-isiz(1)/2*scale(k):scale(k):(isiz(1)/2-1)*scale(k)];
[Z,Fx,Fy]=calc_image_fft(M,X,Y,origo,scale(k),scale(k)); %Calculate Fourier transform (normalized) of the image with the FFT
%%%%%%%%%%

     [Xm,Ym]=meshgrid(X,Y); %Matrix form
   
    %PSF, just a guess since we don't know the true PSF
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    width=20e-3; %HWHM
     psf=1./(1+(Xm.^2+Ym.^2)/(width)^2); %We use Lorentzian distribution here
    %%%%%%%%%%%%%%%%%%%%
Fx(abs(Fx)<1e-14)=0; %Just to avoid issues with numerical precision
Fy(abs(Fy)<1e-14)=0;


Zpsf=calc_image_fft(psf,X,Y,origo,scale(k),scale(k)); %FT of the PSF, normalized



Iy=Fy>=0;
Ix0=Fx==0;
Iy0=Fy==0;
pInd=abs(Zpsf)>1e-8;
Index=Iy&~(Ix0&Iy0)&pInd;
FT.data{k}=Z(Index);
FT.freq{k}=[Fx(Index) Fy(Index)];
FT.psf{k}=Zpsf(Index);
 end
  

