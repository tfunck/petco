int smooth_image(data* image, double fwhm);

int gaussian_filter( double* input, double fwhm, unsigned int* N, double* step);
int forward(double *input, fftw_complex* signal, fftw_complex* SIGNAL, int* N);
int inverse(double *output, fftw_complex* signal, fftw_complex* SIGNAL, int* N);
double* create_filter(double fwhm, int* N, double* step);
int remove_pad(double* img, double* img_pad, int* N, int* M);
double* pad(double* img, int* N, int* M, int* N_pad);
double* pad_like(double* img, int* N, int* M);
int filter(fftw_complex* im1, fftw_complex* im2, int n);
int fftshift(double* img, int* n);
void swap(double* img, int i1, int i2);

