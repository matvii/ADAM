First build mex-files by running build_mex.m. Both gcc and icc are supported, but icc preferred since it produces considerably more effective code. 
Number of threads used for computing (routines Generate_AOFT_Matrix_mex and Generate_RD_Matrix) can be set in file num_of_threads.h
We cannot provide real data with the examples, since the most data sets are not publicly available and we do not have permission to distribute them.
Currently there are following example programs:
AO_invert_Octantoid_Example.m, which reconstructs a shape using octantoid representation from simulated AO data

AO_invert_Subdiv_Example.m, which uses subdivision surface for shape representation. 

Thermal_Interf_Example.m reconstructs shape from simulated thermal infrared interferometry observations.

Metis_LC_Subdiv_Example.m reconstructs a shape of asteroid Metis from lightcurves only (data provided by DAMIT) using subdivision surfaces
Metic_LC_Octantoid_Example.m reconstructs a shape using octantoid representation. It is instructive  to compare resulting shapes.
It is evident that we cannot reconstruct unique nonconvex shape from lightcurves only.

Finally, in the Extras folder there is a more complicated example (AO DATA IS NOT INCLUDED, SO IT CANNOT BE RUN) that reconstructs a shape asteroid Metis from lightcurves and two adaptive optics images  with a nontrivial PSF.

Suggestions, comments, bug reports and data are welcome!

Matti.viikinkoski@gmail.com
