 %c99 or gnu99  flag must be included in compiler options
mex CFLAGS="\$CFLAGS -fopenmp -std=c99 -O3" LDFLAGS="\$LDFLAGS -fopenmp" mex/find_vis_mex.c mex/find_vis.c  mex/findblockers.c  mex/cross.c mex/is_in_triangle.c  -output Calc_Vis
mex CFLAGS="\$CFLAGS -fopenmp -std=c99 -O3" LDFLAGS="\$LDFLAGS -fopenmp" mex/calculate_lc_mex.c mex/calculate_lc.c mex/rotate.c mex/findblockers.c mex/find_vis.c mex/is_in_triangle.c mex/SH.c mex/cross.c -output calculate_lc_mex
mex CFLAGS="\$CFLAGS -fopenmp -std=c99 -O3" LDFLAGS="\$LDFLAGS -fopenmp" mex/FacetsOverHorizon_mex.c mex/FacetsOverHorizon.c mex/cross.c -output FacetsOverHorizon_mex
mex CFLAGS="\$CFLAGS -fopenmp -std=c99 -O3" LDFLAGS="\$LDFLAGS -fopenmp" mex/Calc_FTC2_mex.c -output calc_ftc2
mex CFLAGS="\$CFLAGS -fopenmp -std=c99 -O3" LDFLAGS="\$LDFLAGS -fopenmp" mex/Calculate_AO.c mex/Calculate_AO_mex2.c mex/cross.c mex/rotate.c mex/FindActualBlockers.c mex/FacetsOverHorizon.c mex/is_in_triangle.c mex/Calc_FTC2_mx.c mex/real_matrix_multiply.c -output Calculate_AO
mex  CFLAGS="\$CFLAGS  -std=c99 -O3 -Wall -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" -output Calculate_RD mex/Calculate_Range_Doppler_calloc.c mex/is_in_triangle.c mex/cross.c mex/FacetsOverHorizon.c mex/rotate.c mex/Calc_FTC2_mx.c mex/Calculate_RD_radar_mex_cell2.c
mex  CFLAGS="\$CFLAGS  -std=c99 -O3 -Wall -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" mex/Calculate_Temp_deriv.c mex/FindActualBlockers.c ...
mex/rotate.c mex/cross.c mex/FacetsOverHorizon.c mex/is_in_triangle.c mex/real_matrix_multiply.c  mex/Calculate_Radiance.c mex/Calculate_HF.c mex/Calc_FTC2_mx.c mex/Calculate_HF_mex.c mex/Kissfft/kiss_fftr.c mex/Kissfft/kiss_fft.c -output Calculate_HF
%
