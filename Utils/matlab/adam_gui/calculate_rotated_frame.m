function [upr,angle]=calculate_rotated_frame(fitsfile,E)
import matlab.io.*
fptr=fits.openFile(fitsfile);
NOCD=0;
try
CD11=fits.readKeyDbl(fptr,'CD1_1');
catch
    NOCD=1;
end
try
CD12=fits.readKeyDbl(fptr,'CD1_2');
catch
    NOCD=1;
end
try
CD21=fits.readKeyDbl(fptr,'CD2_1');
catch
    NOCD=1;
end
try
CD22=fits.readKeyDbl(fptr,'CD2_2');
catch
    NOCD=1;
end
if NOCD==1
    upr=[0    0.3977    0.9175];
    angle=0;
    fits.closeFile(fptr);
    return;
end
    
    
sgn=sign(CD11*CD22-CD12*CD21);
if sgn>0
    disp('sign is positive')
end
angle1=atan2d(sgn*CD12,CD22);
angle2=atan2d(-CD21,sgn*CD11);
%up=[0,0,1];
%keyboard
up=[0    0.3977    0.9175];
if(abs(angle1-angle2)>0.5)
    error('Rotation angle information in fits file do not agree, probably due to skewness');
end
    angle=deg2rad(angle1);
M=Ortho_Proj(E,up); %Project E to [0,0,1];
%angle=angle-23.4375*pi/180;

Rz=[cos(angle) -sin(angle) 0;sin(angle) cos(angle) 0; 0 0 1];
angle=rad2deg(angle);
Mrot=M'*Rz'*M;
upr=(Mrot*up')';
fits.closeFile(fptr);