#ifndef STRUCTS
#define STRUCTS
typedef struct LC
    {
        int nlc; /*number of lightcurves*/
        int *nobs; /*Array, number of points per lightcurve*/
        int *rel; /*Relative or calibrated?*/
        int ncalib; /*Number of calibrated lightcurves*/
        int calib;
        int ntotal; /*Total number of points*/
        double **lcs; /*Observations, lcs[n-1] points to observations of nth curve*/
        double **E; /*View directions*/
        double **E0; /*Illumination directions*/
        double **TIME; /*obs times*/
    } LCstruct;
typedef  struct AO
    {
        int nao; /*number of AO images*/
        int ntotal; /*total number of observation points*/
        int *nobs; /*Array, number of points per ao image*/
        double **datar; /*Observations, Fourier transform of images*/
        double **datai;
        double **psfr; /*Fourier transformed psf, optional*/
        double **psfi;
        double **freqx; /*Frequency points of FT, corresponding to data*/
        double **freqy;
        double *E; /*view directions, 3*nao array*/
        double *E0; /*obs directions 3*nao array*/
        double *TIME; /*obs time, nao array*/
        double *distance; /*distance in au, nao array*/
        double *scalex; /*pixel scale, nao array*/
        double *scaley; /*pixel scale, nao array*/
        double *up; /*orientation of instrument 3*nao array*/
    } AOstruct;

    typedef struct RD
    {
        int nRD;
        int ntotal;
        int *nobs;
        double **datar; /*Observations, Fourier transform of images*/
        double **datai;
        double **psfr; /*Fourier transformed psf, optional*/
        double **psfi;
        double **freqx; /*Frequency points of FT, corresponding to data*/
        double **freqy;
        double *E; /*view directions, 3*nao array*/
        double *TIME; /*obs time, nao array*/
        double *distance; /*distance in au, nao array*/
        double *scalex; /*pixel scale, nao array*/
        double *scaley; /*pixel scale, nao array*/
        double *rfreq; /*Radar frequency*/
    } RDstruct;
        
typedef  struct HF
    {
        int nhf; /*number of thermal images*/
        int ntotal; /*total number of observation points*/
        int *nobs; /*Array, number of points per thermal image*/
        double **datar; /*Observations, Fourier transform of images*/
        double **datai;
        
        double **freqx; /*Frequency points of FT, corresponding to data*/
        double **freqy;
        double *E; /*view directions, 3*nao array*/
        double *E0; /*obs directions 3*nao array*/
        double *TIME; /*obs time, nao array*/
        double *distance; /*distance in au, nao array*/
        double *Hdistance; /*distance to Sun*/
        double *scalex; /*pixel scale, nao array*/
        double *scaley; /*pixel scale, nao array*/
        double *up; /*orientation of instrument 3*nao array*/
        double *WL; /* obs wavelength*/
    } HFstruct;
typedef    struct OC
    {
        int noc; /*Number of occultations*/
        int *nobs; /*Array, number of chords per occultation*/
        int ntotal; /*Total number of chords*/
        double **data; /*Chords, points to nobsx4 array where coordinates of starting and ending points*/
        double **TIME; /*observation times*/
        double *E; /*view directions nocx3 array, normalized to unit*/
        double *up; /*camera up direction, nocx3 array, normalized to unit*/ 
        double *V; /*Velocity in km/s*/
        double *distance; /*distance in km*/
        double **etime; /*errors in observation times points to nobsx2 array*/
        int **type; /*chord types, points to noc array*/
    } OCstruct;
    
typedef struct CNTR
{
    int ncont; /*Number of contours*/
    int *nobs; /*Array, number of points per contour*/
    int ntotal; /*Total number of points*/
    double **datax; /*Points to an array with nobs[index] doubles, x-coordinates of contours*/
    double **datay; /*Points to an array with nobs[index] doubles, y-coordinates of contours*/
    double *TIME; /*observation times, ncont array*/
    double *E; /*view directions ncontx3 array, normalized to unit*/
    double *E0; /*sun directions ncontx3 array, normalized to unit*/
    double *up; /*camera up direction ncontx3 array, normalized to unit*/
    double *distance; /*distance in km*/
} CNTRstruct;
typedef struct NOaAR
{
    /*
     * Struct to hold normals and areas of facets and their derivatives wrt vertices
     */
    double *normal;
    double *dndx1; //3*nfac vector
    double *dndx2;
    double *dndx3;
    double *dndy1;
    double *dndy2;
    double *dndy3;
    double *dndz1;
    double *dndz2;
    double *dndz3;
    
    double *area;
    double *dadx; //3*nfac vector
    double *dady;
    double *dadz;
} NOaARstruct;
#endif
