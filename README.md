#ADAM: All-Data Asteroid Modelling

ADAM is a program for 3D shape reconstruction of asteroids from disk-resolved observations. Current version supports (relative) lightcurves, adaptive optics, occultations, and range-Doppler radar images. Raw images are used directly, boundary extraction is not needed. Source code for thermal modelling is also included, but not currently integrated into main program due to lack of available observations (Juno is currently the only one).

##Required libraries
- KissFFT (https://sourceforge.net/projects/kissfft)
- Iniparser (https://github.com/ndevilla/iniparser)
- Wcstools (http://tdc-www.harvard.edu/wcstools)
- Lapacke
- Lapack
- Blas
- Cblas

KissFFT, Iniparser and Wcstools folders are included for convenience; they come with their own licenses.

##Building
- Build Iniparser
- Build Wcstools
- Change Makefile to reflect installed libraries, i.e. change reflapacke, reflapack, atlas and atlcbas.
- make adam

This program has been tested on Linux with gcc 5.3.

##Usage
ADAM uses ini files to process data. For an example of asteroid Hertha, 135_1.ini and 135_1_oct.ini (usage: ./adam 135_1.ini). For all available config options and short descriptions, see the file Adam.ini.
Adam.pdf contains some useful tips. For theoretical background, read [Shape reconstruction from generalized projections](http://urn.fi/URN:ISBN:978-952-15-3673-1).

##License

This software is licensed under [CC Attribution 4.0 international License](https://creativecommons.org/licenses/by/4.0/legalcode).

If you use ADAM in your research, please cite
Viikinkoski, M; Kaasalainen, M.; Durech, J.: *ADAM: a general method for for using various data types in asteroid reconstruction*, A&A Vol 576, 2015.

##Contact
Bug reports, data, feature suggestions and comments are welcome.

Matti Viikinkoski (matti.viikinkoski@gmail.com)

##Updates
####10.7.2016
- Added support for occultation chord offsets
- Added support for optimizing occultation chord offsets

####23.7.2016
- Added support for calibrated lightcurves

####TODO
- Documentation
- Albedo variegation
- Hapke scattering law
