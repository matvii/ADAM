function [im,date,cdelt1,cdelt2,crpix1,crpix2]=process_RD_files(fitsfile,x,y)
%return 2x x 2y image centered at (crpix1,crpix2)
import matlab.io.*
fptr=fits.openFile(fitsfile);
date=0;
try
date=fits.readKeyDbl(fptr,'JDMEAN');
catch
    date=NaN;
end
try 
    cdelt1=fits.readKeyDbl(fptr,'CDELT1');
catch
    cdelt1=NaN;
end
try 
    cdelt2=fits.readKeyDbl(fptr,'CDELT2');
catch
    cdelt2=NaN;
end
try 
    crpix1=fits.readKeyDbl(fptr,'CRPIX1');
catch
    crpix1=NaN;
end
try 
    crpix2=fits.readKeyDbl(fptr,'CRPIX2');
catch
    crpix2=NaN;
end
try 
    naxis1=fits.readKey(fptr,'NAXIS1');
catch
    naxis1=NaN;
end
try 
    naxis2=fits.readKey(fptr,'NAXIS2');
catch
    naxis2=NaN;
end


M=fits.readImg(fptr);
fits.closeFile(fptr);
if x>0 && y>0
cx=round(crpix1);
cy=round(crpix2);
nx=2*x;
ny=2*y;
x0=cx-x;
y0=cy-y;
if isnan(crpix1) || isnan(crpix2)
    error('Error reading CRPIX keyword');
end
im=M(y0:ny-1+y0,x0:nx-1+x0);
else
im=M;
end
end