function date=read_fits_date_RD(fitsfile)
import matlab.io.*
fptr=fits.openFile(fitsfile);
date=0;


    try
        date=fits.readKeyDbl(fptr,'JDMEAN');
    catch
        date=NaN;
    end
fits.closeFile(fptr);
end