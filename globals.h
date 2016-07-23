#ifndef GLOBALS
#define GLOBALS

extern int     INI_HAVE_LC;
extern int     INI_HAVE_AO;
extern int      INI_HAVE_OC;
extern int      INI_HAVE_HF;
extern int      INI_HAVE_RD;
extern double INI_ANGLE_B;
extern double INI_ANGLE_L;
extern double INI_ANGLE_P;
extern double INI_ANGLE_PHI0;
extern int INI_LMAX;
extern int INI_SD_LEVEL;
extern int INI_LC_ARE_RELATIVE;
extern char *INI_SHAPE_FILE;
extern char *OUT_SHAPE_FILE;
extern char *OUT_PARAM_FILE;
extern char *OUT_SHAPE_PARAM_FILE;
extern char *INI_INPUT_AO_OFFSET;
extern char *INI_OUTPUT_AO_OFFSET;
extern char *INI_INPUT_RD_OFFSET;
extern char *INI_OUTPUT_RD_OFFSET;
double *INI_OC_OFFSET;
extern char *OUT_LC_FILE;
extern int USE_ELLIPSOID;
extern double ELLIPSOID_SEMI_A;
extern double ELLIPSOID_SEMI_B;
extern double ELLIPSOID_SEMI_C;
extern double INI_MIN_TIM;
extern double INI_RDEXP;
extern int INI_NROWS;
extern int NUM_OF_ROUNDS;
extern double INI_LCW;
extern double INI_AOW;
extern double INI_CW;
extern double INI_AW;
extern double INI_DW;
extern double INI_OW;
extern double INI_OCW;
extern double INI_RW;
extern double INI_CHRDW;
extern double INI_LAMBDA;
extern double INI_LAMBDAINC;
extern double INI_LAMBDAMAX;
extern double INI_MINDEC;
extern double INI_CALIBLCW;
extern double *INI_PARAMS;
extern double *INI_AO_WEIGHT;
extern double *INI_OC_WEIGHT;
extern double *INI_RD_WEIGHT;
extern double *INI_PHASE_PARAMS;
extern int INI_AO_SCALING;
extern AOstruct *INI_AO;
extern LCstruct *INI_LC;
extern OCstruct *INI_OC;
extern RDstruct *INI_RD;


extern int INI_FIX_SHAPE;
extern int INI_FIX_ANGLES;
extern int INI_MASK_SET;
extern int *INI_PARAMETER_MASK;
extern double *INI_CHORD_OFFSET;
extern int *INI_FREE_CHORD_LIST;
extern int INI_FREE_CHORD_NMR;
extern int *INI_PHASE_MASK;
#endif
