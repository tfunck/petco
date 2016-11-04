int smooth(float* mask, float* image, int* count, double fwhm, double dt, double step, double lambda,  int nthreads, int VERBOSE);

struct grad_args {
    float* image;
    float* dx;
    float* dy;
    float* dz;
    int* count;
    int thread;
    int nthreads;
};

struct div_args {
    float* image;
    int* count;
    float* dx;
    float* dy;
    float* dz;
    int thread;
    int nthreads;
};

struct update_img_args {
    float* temp_dx;
    float* temp_dy;
    float* temp_dz;
    float* img_dz;
    float* img_dy;
    float* img_dx;
    float** big_tensor;
    int* count;
    int thread;
    int nthreads;
};

struct diffusion_args {
    float* image;
    float* temp_image;
    double scale_factor;
    int thread;
    int* count;
    int nthreads;
};


struct structure_tensor_args {
    float** big_tensor;
    float* dx;
    float* dy;
    float* dz;
    int* count;
    double lambda;
    int thread;
    int nthreads;
};
