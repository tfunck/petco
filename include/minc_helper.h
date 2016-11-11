#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <malloc.h>
#include <libgen.h>
#include <math.h>
#include <omp.h>
#include "minc2.h"
#define TRUE 1
#define FALSE 0


int VERBOSE;

typedef struct temp_point {
    int index;
	int id;
	int z;
	int y;
	int x;
	double wx;
	double wy;
	double wz;
	float value;

} point;

typedef struct image {
    char* filename;
    char* values_filename;
    //char* like_filename
    char** groups;
    double* data;
	point* points;
    misize_t* wcount;
    misize_t* wstarts; 
    unsigned int* count;
    double* start;
    double* step;
    double val;
    int n_roi;
    int ngroups;
    int ndim;
    int index;
    int n;
    int n3d;
    int label;
    int tmax;
    int zmax;
    int ymax;
    int xmax;
    double xstep;
    double ystep;
    double zstep;
    double zstart;
    double ystart;
    double xstart;
    int tstep;
    int smooth;
} data;


/*typedef struct temp_vertex {
  double z;
  double y;
  double x;
  double dist;

  int index;
} vertex;*/

typedef struct temp_vertex {
  double z;
  double y;
  double x;
  int i, j, k; //voxel indices
  double dist;
  double val;
  double previous_distance;
  //for voxels the pointer below should 
  //point to the nearest vertex. For vertices
  //the pointer should go to the ngh 
  //of that vertex
  struct temp_vertex *previous;
  //the info below only applies to vertices
  //on a mesh, not to voxels
  struct temp_vertex** ngh;
  struct temp_edge** incoming_edges;
  struct temp_edge** outgoing_edges;
  int n_incoming_edges;
  int n_outgoing_edges;
  int inside_brain;
  int inside;
  int border;
  int nngh;
  int index;
  int n;
} vertex;




int allocate_vol(data** ptr, int* n);
int allocate_atlas(data** ptr, int* n);
double* malloc_check(int dim1, int dim2, int dim3, int dim4, int ndim, int dtype);
int writeVolume(char* outputfile, void* array, double starts[3],  double separations[3],  misize_t count[3], int dtype  );
void* readVolume(data* volume, int ndim, int dtype);
int check_input_files(char *file_inputs[]);
int checkDimensions(data* data1, data* data2);
void pexit(char* string1, char* string2, int number);
int check_file(char *file_input);
int read_obj(data* surface);
int read_dimensions(mihandle_t* img, data* dataume, int ndim);
int get_dimensions(mihandle_t* img, data* dataume);
char* strip_path_ext(char* instr);
double nearest_neighbour(double* array, int  zmax, int sz, double oz, int  ymax,int  sy, double  oy, int xmax,int  sx, double  ox,int  xi,int yi,int zi);
int allocate_data(data** ptr, int* n);
void createVolume(char* newfilename, int num_dim, misize_t *sizes, double *separations, double *starts, mitype_t data_type);
int volume_to_obj(data* volume);
int read_frame(data* image, int t, void* array, mitype_t dsize);
double trilinear(double* array, int zmax, double sz, double oz, int ymax, double sy, double oy,int xmax, double sx, double ox, double xi, double yi, double zi);
int nearest_voxel(double vertex, double step, double start, int nvox, int direction /*-1=smallest, 1=largest*/);
double nearest_voxel_world(double vertex, double step, double start, int nvox, int direction /*-1=smallest, 1=largest*/);
float voxel2real(int location, float min, float step);
int real2voxel(float location, float min, float step);
int writeVolume4d(char* outputfile, void* array,  double starts[4],  double separations[4],  misize_t  count[4], mitype_t dtype  );
int writeVolumeSubset(char* outputfile, void* array, misize_t start[3], double starts[3],  double separations[3],   misize_t  count[3], mitype_t dtype  );
data* copyVolume(data* vol);
int extrema(double* array, int tmax, int zmax, int ymax, int xmax, double* min, double*max);
