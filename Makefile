#Set CC if icc if you want to use icc
CC=gcc
#Location of INTEL math kernel library (if it exists)
#MKLROOT=/opt/intel/compilers_and_libraries_2016.2.181/linux/mkl
ifeq ($CC,icc)
CFLAGS=-O2 -std=c99  -qopenmp  -xHost -pthread 
LD=xild
AR=xiar
else
CFLAGS=-O2  -std=c99 -fopenmp  -march=native -pthread 
endif
ifdef MKLROOT
LDLIBS=-L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core  -L$(MKLROOT)/../compiler/lib/intel64 -liomp5 -lpthread -lm -ldl -fopenmp  wcstools-3.9.2/libwcs/libwcs.a -I$(MKLROOT)/include 
CFLAGS:=$(CFLAGS) -I$(MKLROOT)/include 
DEFINE=USE_MKL
else
LDLIBS=-lm  -lreflapacke -lreflapack -latlas -latlcblas -fopenmp  wcstools-3.9.2/libwcs/libwcs.a 
#LDLIBS=-lm  -llapacke -llapack -lblas -lgslcblas -fopenmp -I/usr/include/lapacke wcstools-3.9.2/libwcs/libwcs.a
LDLIBS=-lm  -lopenblas -fopenmp -lpthread wcstools-3.9.2/libwcs/libwcs.a

DEFINE=NO_MKL
endif
CPPFLAGS=$(addprefix -D,$(DEFINE))
#LDLIBS=-L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core  -L$(MKLROOT)/../compiler/lib/intel64 -liomp5 -lpthread -lm -ldl -fopenmp  wcstools-3.9.2/libwcs/libwcs.a
LIBS=phase_function.o read_occ.o Fit_Occ.o calculate_OCs.o find_chord.o  SH.o generate_sphere.o fit_subdiv_model_to_LC_AO_QR.o parse_ini.o calculate_lcs.o Calc_FTC.o calculate_lcurve.o dihedral_angle.o dihedral_angle_reg.o utils.o convex_reg.o area_reg.o sqrt3_subdiv.o  triangulate_sphere.o matrix_ops.o rotate.o is_in_triangle.o FacetsOverHorizon.o FindActualBlockers.o Calculate_AOs.o Calculate_AO.o calc_image_fft_sinc.o read_ephm_data.o process_ao_images.o octantoid_reg.o octantoid_to_trimesh.o fit_oct_model_to_LC_AO.o process_rd_images.o calc_image_fft_unnormed.o Calculate_RDs.o Calculate_Range_Doppler.o readfits_rd.o kiss_fftndr.o kiss_fft.o kiss_fftr.o kiss_fftnd.o readfits.o hapke_brightness.o soft_maxdimz.o

adam:	adam.c $(LIBS)
	$(CC) $(CPPFLAGS) $(CFLAGS)    adam.c $(LIBS)  $(LDLIBS) -Iiniparser/src -Liniparser -liniparser -o adam
	
thermal:	thermal.c $(LIBS)
	$(CC) $(CPPFLAGS) $(CFLAGS) thermal.c Calculate_Temp.c $(LIBS) $(LDLIBS)  -Iiniparser/src -Liniparser -liniparser -o thermal
libs:	$(LIBS)

kiss_fftndr.o:	Kissfft/tools/kiss_fftndr.c
	$(CC) $(CFLAGS) -c $<

kiss_fft.o:	Kissfft/kiss_fft.c
	$(CC) $(CFLAGS) -c $<

kiss_fftr.o:	Kissfft/tools/kiss_fftr.c
	$(CC) $(CFLAGS) -c $<

kiss_fftnd.o:	Kissfft/tools/kiss_fftnd.c
	$(CC) $(CFLAGS) -c $<


%.o:	%.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $<

clean:
	rm *.o 

