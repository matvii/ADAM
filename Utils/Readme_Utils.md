#Utils

##MATLAB

####Displaying the shape
[tlist,vlist]=read_shape('path_to_shape_file',0);
figure; trisurf(tlist,vlist(:,1),vlist(:,2),vlist(:,3)); axis equal

####Displaying the projections corresponding to the AO files
Assuming that MATLAB and the ini file are in the same directory (and paths to AO files are correct relative to MATLAB working directory), projections of reconstructed shape can be displayed with function 
Display_ao_projections('inifile')

##Python

####Displaying the shape:
python ./Display_shape.py shapefile

####Plotting occultations:

python ./Plot_Occ.py -T JD0 -O OCfile -S shapefile -A " beta,lambda,P" -o " offsetx,offsety" -I occultationindex
Where JD0 is the zero time for the model, OCfile is the occultation file, shapefile is the Shapefile, pole direction and rotation period are beta, lambda and P. Note the extra space after ". Quotation marks and space are required if beta or offsetx are neqative. 
Offset values (in km) of plane projection are offsetx and offsety. Occultationindex is the index of occultation if there are several occultations in the file.
This program can be also used to plot an ellipsoid:
python ./Plot_Occ.py -T JD0 -O OCfile -E sizex,sizey,sizez -A " beta,lambda,P" -o " offsetx,offsety" -I occultationindex
It is useful for determining good initial values for OCCOffset in ini file.
Example:
./Plot_Occ.py -T 2443846.0 -O Hertha.occ -E 20,10,30 -A " 51.9582,275.2514,8.40059630"

####Ephm file generation:
Assuming MJD-OBS keyword is present in the fits files, ephemeris information can be retrieved from Horizons web service using command
python ./fetch_ephm_fits.py asteroidnumber '*.fits' ephm.dat
Similarly, if observation dates are listed in a file 'dates':
python ./fetch__ephm_time_file asteroidnumber dates ephm.dat    
