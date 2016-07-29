/* Function prototypes for FFT routines */

#define GAUSS_FWHM_MM (8.0)  /* FWHM used for Gaussian filter (mm) in x,y,z */
#define GAUSS_STD_DEV (0.0)  /* Standard deviation used for Gaussian filter *
                              * Overrides FWHM if not zero (0.0)            */
#define GAUSS_MAX     (1.0)  /* Maximum parameter used for Gaussian filter  */
#define GAUSS_MEAN    (0.0)  /* Mean parameter used for Gaussian filter     */
#define NO_FILTER 0          /* Forward and inverse FFT only thus no filter */
#define GAUSSIAN_FILTER 1
#define LOW_PASS_FILTER 2
#define HIGH_PASS_FILTER 3
#define IDENTITY_FILTER 4

int isPowerOf2(unsigned long n);
int nextPowerOf2(unsigned long n);

void fourn(float data[], unsigned long nn[], int ndim, int isign);
void rlft3(float ***fftdata, float **speq, unsigned long nn1, unsigned long nn2,
                   unsigned long nn3, int isign);

void normaliseInverse(float ***fftdata, unsigned long nn1, unsigned long nn2,
                   unsigned long nn3);

void loadImageData3d(float ***imagedata,
                     long i_dim, long j_dim, long k_dim,
                     long i_pad, long j_pad, long k_pad,
                     float ***iFFTdata);

void unloadImageData3d(float ***imagedata,
                       long i_dim, long j_dim, long k_dim,
                       long i_pad, long j_pad, long k_pad,
                       float ***fftdata);

void loadFilterData3d(float ***filterdata,
                      long i_dim, long j_dim, long k_dim,
                      long i_pad, long j_pad, long k_pad,
                      float ***fFFTdata);

void displayImage(float ***imagedata, int nrows, int ncols);

void applyFFT(float ***imagedata, float*** filterdata,
              int nslices, int nrows, int ncols);


void applyFilter3d(float*** iFFTdata, float*** fFFTdata,
                   float** ispeq,  float** fspeq,
                   int nn1, int nn2, int nn3);

float*** allocateKernelMatrix(unsigned int n1, unsigned int n2, unsigned int n3);

float*** genGaussianKernel(float* step,
                           unsigned int n1, unsigned int n2, unsigned int n3);

float*** genIdentityKernel(unsigned int n1, unsigned int n2, unsigned int n3);

void freeKernelMatrix(float*** matrix, unsigned int n1, unsigned int n2);

void psfGaussian(float*** psf, int xdim, int ydim, int zdim,
                 double xstep, double ystep, double zstep,
                 float xfwhm, float yfwhm, float zfwhm,
                 float xstdev, float ystdev, float zstdev,
                 float xcntrd, float ycntrd, float zcntrd, int normalise);

void astroGaussian(float* xi, int xidim, float* gauss,
                   float max, float mean, float sigma);
