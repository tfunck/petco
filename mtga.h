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
  int analysis_type; /*0 = Logan Plot, 1=Patlak Plot, 2=SRTM, 3=Split Frames, 4=Gen. Brain Mask, 5=ROI Stats, 6=SUV*/
  int write;
  int min; /*converts from seconds to minutes*/
  int input3d;
  int input4d;
  int test;
  int integrate;
  int bp;
  int nosmooth;
  int arterial;
  int arterial_convert_to_seconds;
  float RAMthreshold;
  int RAMexceeded;
  int testStartTime;
  char *regression_type;
  double brainmask_threshold;
  int    autok2;
  double injected_dose;
  double patient_weight;
  double start_frame;
  double end_frame;
  double k2;
  double conversion_factor;
  double LC;
  double PG;
  int CMRg;

} parameters;


/******************************
*   GLOBAL VARIABLES	      *
*******************************/


/******************************
* FUNCTIONS		      *
*******************************/
extern float** importArterialTAC (char* arterial_filename, int* points, parameters* userOption ) ;
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

