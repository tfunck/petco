 typedef struct  {
    float** matrix;
    int start;
    int row;
    int col;
} vectorArg;   
               

void scalar_mult(float scalar, vectorArg *argsA, vectorArg *argsB);   
void matrix_mult(vectorArg *argsA, vectorArg *argsB, vectorArg *argsC); 
void matrix_add(vectorArg *A,  vectorArg *B, vectorArg *C, int operation);   
void gramm_schmidt(float **input, float **output,int  m, int n);      
void normalize(float **A,int row, int col, float *sum); 
void qrdecomp(vectorArg *input, vectorArg *Q, vectorArg *R, int col, int row);                
void proj(float **A, int startA, float **B, int startB, int n, float **C );  
void solve_qr(float** Q, float** R, float *Ct, float* theta1, float* theta2, int col, int row, vectorArg* argsB, vectorArg* Qt);     
void transpose(float **A, int col, int row, float** B);  
float srtmBasisFn(float *Ct, float* Cr, int nframes,  
                   float *function, float *ans, float *response, 
                   vectorArg *argsQ, vectorArg *argsR, vectorArg *argsA, float* weights, int ctStep, int ctStart ) ;
float dot(vectorArg *argsA, vectorArg *argsB);   
float magnitude(float *A, int row);        
float* scalar2vector_mult(float scalar, float* A, int row, int col);   




float srtmBasisResult(float *Ct, float *Cr,  int nframes, float* weights,  int basisMax,double  lambda, vectorArg *Q_array, vectorArg *R_array, float **basis, float *theta3Array, vectorArg *argsB, vectorArg *argsQt, int *numMax ); 

