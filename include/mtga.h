#ifndef TRUE
#  define TRUE 1
#  define FALSE 0
#endif

#define EXIT_SUCCESS	0

/******************************
*			STRUCTURES	 	   *
*******************************/

typedef struct {
    char name[30];
    double lambda;
    double k2mean;
    double k2max;
    double k2min;
} tracer;

typedef struct {
   mihandle_t* outputimage;
   char* outputimage_file;
   char* smoothedImageFile;
   mihandle_t* smoothedImage;
   misize_t* sizes;
   float* timeWidths;
   float* times;
   float* referencemask;
   float* brainmask;
   float* referenceArray;
   float** arterialTAC;
   tracer *tracerInfo;
} imageData;

typedef struct {
  char* pet_fn;
  char* ref_mask_fn;
  char* brain_mask_fn;
  char* time_widths_fn;
  char* arterial_fn;
  char* output_fn;
  char* analysis_type; 
  char* time_unit;
  char* output_unit;
  char* regression_type;
  char* analysis;
  double injected_dose;
  double body_weight;
  double k2;
  double conversion_factor;
  double LC;
  double PG;
  int start_frame;

} parameters;



/******************************
*   GLOBAL VARIABLES	      *
*******************************/


/******************************
* FUNCTIONS		      *
*******************************/
double** importArterialTAC (char* arterial_filename, int* points, parameters* userOptions);
extern void smoothFrames(mihandle_t *petimage, mihandle_t *tempimage, unsigned  int sizes[4], int filtertype, double step[4], parameters userOptions);
extern  void find_reference_values(char* petimage_file, misize_t *sizes, float **reference_array,  float* referencemask) ;
extern float* get_arterial_TAC(int nframes,float* times, char* arterial_filename, parameters* userOption);
extern double loadpetimage( int z, misize_t *sizes, misize_t *start, misize_t *count,  char* petimage_file, float *slab, float *brainmask,  float *referencemask, parameters userOptions)  ;
int getdimensions(char *image_file, mihandle_t *image, int num_dimensions, misize_t *sizes, midimhandle_t *dimensions, double* step, double* separations, double* starts);
void getTimeWidths(int nframes, int startFrame,  float* timeWidths, char* inputfile, parameters userOptions);
void getTimes(mihandle_t* volume, int nframes,float* times, float* timeWidths, char* inputfile, parameters userOptions);
void loadmask_hyperslab(char* mask_file, float** mask);
tracer* getTracerParameterValues(mihandle_t *petimage, parameters* userOptions);
void createImage(char* newfilename, int num_dim, misize_t *sizes, double *separations, double *starts);

