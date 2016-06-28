CC=gcc

CFLAGS=-O3 -std=c99  -fopenmp  -march=native
LDLIBS=-lm  -lreflapacke -lreflapack -latlas -latlcblas -fopenmp  wcstools-3.9.2/libwcs/libwcs.a

adam:	adam.c 
	$(CC) $(CFLAGS)    adam.c  phase_function.c read_occ.c Fit_Occ.c calculate_OCs.c find_chord.c  SH.c generate_sphere.c fit_subdiv_model_to_LC_AO.c parse_ini.c calculate_lcs.c Calc_FTC.c calculate_lcurve.c dihedral_angle.c dihedral_angle_reg.c utils.c convex_reg.c area_reg.c sqrt3_subdiv.c  triangulate_sphere.c matrix_ops.c rotate.c is_in_triangle.c FacetsOverHorizon.c FindActualBlockers.c Calculate_AOs.c Calculate_AO.c calc_image_fft_sinc.c read_ephm_data.c process_ao_images.c octantoid_reg.c octantoid_to_trimesh.c fit_oct_model_to_LC_AO.c process_rd_images.c calc_image_fft_unnormed.c Calculate_RDs.c Calculate_Range_Doppler.c readfits_rd.c Kissfft/tools/kiss_fftndr.c Kissfft/kiss_fft.c Kissfft/tools/kiss_fftr.c Kissfft/tools/kiss_fftnd.c readfits.c  $(LDLIBS) -Iiniparser/src -Liniparser -liniparser -o adam

test_fit: test_fit.o  calculate_lcs.o calculate_lcurve.o dihedral_angle.o dihedral_angle_reg.o utils.o convex_reg.o area_reg.o sqrt3_subdiv.o  matrix_ops.o rotate.o is_in_triangle.o FacetsOverHorizon.o FindActualBlockers.o

test_fit.o: test_fit.c structs.h calculate_lcs.o dihedral_angle_reg.o  convex_reg.o area_reg.o sqrt3_subdiv.o 

fit_LC_AO: fit_LC_AO.c calculate_lcs.c Calc_FTC.c calculate_lcurve.c dihedral_angle.c dihedral_angle_reg.c utils.c convex_reg.c area_reg.c sqrt3_subdiv.c  matrix_ops.c rotate.c is_in_triangle.c FacetsOverHorizon.c FindActualBlockers.c Calculate_AOs.c Calculate_AO.c calc_image_fft.o read_ephm_data.o process_ao_images.o kiss_fftndr.o kiss_fft.o kiss_fftr.o kiss_fftnd.o readfits.o
fit_LC_AO.o: fit_LC_AO.c

calculate_lcs.o: calculate_lcs.c  calculate_lcurve.o  

dihedral_angle_reg.o: dihedral_angle_reg.c dihedral_angle.o  

dihedral_angle.o: dihedral_angle.c 

utils.o:	utils.c 

matrix_ops.o:	matrix_ops.c

convex_reg.o:	convex_reg.c   

FacetsOverHorizon.o: FacetsOverHorizon.c  

FindActualBlockers.o: FindActualBlockers.c 

is_in_triangle.o:	is_in_triangle.c

area_reg.o:	area_reg.c 

Calc_FTC.o: Calc_FTC.c

readfits.o: readfits.c

sqrt3_subdiv.o: sqrt3_subdiv.c 

calculate_lcurve.o: calculate_lcurve.c  

rotate.o:	rotate.c
is_in_triangle.o:	is_in_triangle.c

calc_image_fft: calc_image_fft.o utils.c Kissfft/tools/kiss_fftndr.c Kissfft/kiss_fft.c Kissfft/tools/kiss_fftr.c Kissfft/tools/kiss_fftnd.c
	$(CC) calc_image_fft.o utils.c Kissfft/tools/kiss_fftndr.c Kissfft/kiss_fft.c Kissfft/tools/kiss_fftr.c Kissfft/tools/kiss_fftnd.c  -lm  -o calc_image_fft -fsanitize=address

calc_image_fft.o: calc_image_fft.c
	$(CC) -c calc_image_fft.c

read_ephm_data.o: read_ephm_data.c

process_ao_images.o: process_ao_images.c
process_ao_images: process_ao_images.o utils.c read_ephm_data.c readfits.c calc_image_fft.c readfits.c Kissfft/tools/kiss_fftndr.c Kissfft/kiss_fft.c Kissfft/tools/kiss_fftr.c Kissfft/tools/kiss_fftnd.c

kiss_fftndr.o: Kissfft/tools/kiss_fftndr.c
	$(CC) -c Kissfft/tools/kiss_fftndr.c
kiss_fft.o:	Kissfft/kiss_fft.c
	$(CC) -c Kissfft/kiss_fft.c
kiss_fftr.o: Kissfft/tools/kiss_fftr.c
	$(CC) -c Kissfft/tools/kiss_fftr.c
kiss_fftnd.o: Kissfft/tools/kiss_fftnd.c
	$(CC) -c Kissfft/tools/kiss_fftnd.c

Calculate_AOs.o: Calculate_AOs.c
Calculate_AO.o: Calculate_AO.c

temp: process_ao_images.c rotate.c is_in_triangle.c FacetsOverHorizon.c FindActualBlockers.c Calculate_AOs.c Calculate_AO.c Calc_FTC.c utils.c read_ephm_data.c readfits.c calc_image_fft.c readfits.c Kissfft/tools/kiss_fftndr.c Kissfft/kiss_fft.c Kissfft/tools/kiss_fftr.c Kissfft/tools/kiss_fftnd.c

generate_sphere.o: generate_sphere.c

triangulate_sphere.o: triangulate_sphere.c
