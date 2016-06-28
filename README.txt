ADAM is a program for 3D shape reconstruction of asteroids from various observations. Current version supports (relative) lightcurves, adaptive optics, occultations, and range-Doppler radar images. Source code for thermal modelling is also included, but not currently integrated into main program due to lack of available observations (Juno is the only one currently).

This program uses several libraries: Kissfft, wcstools and iniparser are included. In addition, lapack, lapacke, cblas and blas are required. To compile, replace -lreflapacke -lreflapack -latlas -latlcblas in the Makefile with correct libraries. Command make adam will then produce the program. Number of concurrent threads can be set in num_of_threads.h. ADAM is tested in Linux with gcc 5.3.0.

For documentation, read the accompanying manual and see the example 135_1.ini (which can be run with  ./adam 135_1.ini). Available config options are also described in the example ini adam.ini.
