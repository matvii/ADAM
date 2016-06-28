CC=gcc

CFLAGS=-O3 -std=c99  -fopenmp  -march=native -w
LDLIBS=-lm  -lreflapacke -lreflapack -latlas -latlcblas -fopenmp  wcstools-3.9.2/libwcs/libwcs.a

adam:	adam.c 
	$(CC) $(CFLAGS)    adam.c  phase_function.c read_occ.c Fit_Occ.c calculate_OCs.c find_chord.c  SH.c generate_sphere.c fit_subdiv_model_to_LC_AO.c parse_ini.c calculate_lcs.c Calc_FTC.c calculate_lcurve.c dihedral_angle.c dihedral_angle_reg.c utils.c convex_reg.c area_reg.c sqrt3_subdiv.c  triangulate_sphere.c matrix_ops.c rotate.c is_in_triangle.c FacetsOverHorizon.c FindActualBlockers.c Calculate_AOs.c Calculate_AO.c calc_image_fft_sinc.c read_ephm_data.c process_ao_images.c octantoid_reg.c octantoid_to_trimesh.c fit_oct_model_to_LC_AO.c process_rd_images.c calc_image_fft_unnormed.c Calculate_RDs.c Calculate_Range_Doppler.c readfits_rd.c Kissfft/tools/kiss_fftndr.c Kissfft/kiss_fft.c Kissfft/tools/kiss_fftr.c Kissfft/tools/kiss_fftnd.c readfits.c  $(LDLIBS) -Iiniparser/src -Liniparser -liniparser -o adam

