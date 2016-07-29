
struct tempVector;
struct node;

typedef struct tempVector {
  double z;
  double y;
  double x;
  float observed_activity;
  float test_activity;
  float blurred_activity;
  float true_activity;
  float correctionFactor;
  float *distance;
  //for voxels the pointer below should 
  //point to the nearest vertex. For vertices
  //the pointer should go to the neighbours 
  //of that vertex
  struct tempVector* nearest;
  int inside;
  int index;
  //the info below only applies to vertices
  //on a mesh, not to voxels
  int layer;
  struct tempVector* pair;
  struct tempVector** neighbours;
  int nneighbours;
  int nvox;
  float **voxNgh;
} vector;

typedef struct {
    double starts[3], maxima[3], step[3];
    unsigned long count[3];
} dimensions;


typedef struct node {
 struct node*** array;
 double bounds[8][3];
 int nvertices;
 double center[3];
 vector** vertices;
} telescopingArray;

typedef struct {
    vector offset, offset2, offset3, offset4;
    vector line1;
    vector line2;
    vector line3;
    vector line4;
    double mag1; //store the magnitudes of the three line vectors 
    double mag2; //because this is useful info that is annoying to 
    double mag3; //calculate more than once.
} plane;
double** rref(double** matrix, int colSize, int rowSize);
vector* searchTelescopingArray( telescopingArray* matrix, vector* testVector, int scaleFactor, int previousDepth, int maxDepth );
telescopingArray* createTelescopingArray( telescopingArray* matrix, double currentBounds[8][3], int previousDepth, int maxDepth, int scaleFactor, vector** vertices, int nvertices);
void voxel_vertexPairs(vector** vertices, double maxima[3], double minima[3], double step[3], float*** voxelArray,  double xmin, double ymin, double zmin, vector vtx_min, vector vtx_max, vector**** local_voxels );
int pointinMesh(vector* point, vector** inner_mesh, vector** outer_mesh, int n_meshPts);
int interpolate_sphere( char * input_object_filename, int** n_ngh, int*** ngh ); 
int traceRay(vector* point, vector** inner_mesh, vector** outer_mesh, int n_meshPts, plane** local, double rotate );
int intersection( vector *x1, vector* u, vector* v, vector* x0, vector** inner_mesh, vector** outer_mesh,  int meshPoint, int n_meshPoints, double rotate, double** matrix, int rect2 ) ;
plane** findPlanes(vector **inner_mesh, vector** outer_mesh,  int n_meshPts );
plane* calcLines(vector* point1, vector* point2, vector* point3 );

int getdimensions(char *image_file, mihandle_t *image, int num_dimensions, misize_t *sizes, midimhandle_t *dimensions, double* step,  double starts[3]);

void blurImage( unsigned long *sizes, double* step, float*** voxeldata, float fwhm);
float voxel2real(int location, float min, float step);
int real2voxel(float location, float min, float step);
