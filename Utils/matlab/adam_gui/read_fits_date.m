function date=read_fits_date(fitsfile)
import matlab.io.*
fptr=fits.openFile(fitsfile);
date=0;
try
date=fits.readKeyDbl(fptr,'MJD-OBS');
catch
    date=NaN;
end
if isnan(date)
    try
        date=fits.readKeyDbl(fptr,'JDMEAN')-2400000.5;
    catch
        date=NaN;
    end
end
fits.closeFile(fptr);
end