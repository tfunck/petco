#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

void eigen(double *data, double* vector, double* value){
    int n=3;
    gsl_matrix_view m = gsl_matrix_view_array(data, n, n);
    gsl_vector *eval = gsl_vector_alloc(n);
    gsl_matrix *evec = gsl_matrix_alloc(n,n);

    gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(n);
    
    gsl_eigen_symmv(&m.matrix, eval, evec, w);
   
    gsl_eigen_symmv_free(w);

    gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_ASC);
    for(int i=0; i < 3; i++){
        double eval_i = gsl_vector_get (eval, i);
        printf("eigenvalue: %g\n", eval_i);
        gsl_vector_view evec_i= gsl_matrix_column (evec, i);
        //gsl_vector_fprintf (stdout,  &evec_i.vector, "%g");
    }
    value[0]=gsl_vector_get(eval, 0);
    value[1]=gsl_vector_get(eval, 1);
    value[2]=gsl_vector_get(eval, 2);
    vector[0]=gsl_matrix_get(evec, 0, 0);
    vector[1]=gsl_matrix_get(evec, 0, 1);
    vector[2]=gsl_matrix_get(evec, 0, 2);
    vector[3]=gsl_matrix_get(evec, 1, 0);
    vector[4]=gsl_matrix_get(evec, 1, 1);
    vector[5]=gsl_matrix_get(evec, 1, 2);
    vector[6]=gsl_matrix_get(evec, 2, 0);
    vector[7]=gsl_matrix_get(evec, 2, 1);
    vector[8]=gsl_matrix_get(evec, 2, 2);

}
void inverse(double* data, gsl_matrix * inverse, gsl_permutation * perm, gsl_matrix * m ){
    // Define the dimension n of the matrix
    // and the signum s (for LU decomposition)
    int s;
    // Define all the used matrices
    m=gsl_matrix_view_array(data, n, n);
    // Fill the matrix m
    // Make LU decomposition of matrix m
    gsl_linalg_LU_decomp (m.matrix, perm, &s);
    // Invert the matrix m
    gsl_linalg_LU_invert (m.matrix, perm, inverse);
    gsl_vector_view evec_i= gsl_matrix_column (evec, i);

}

int main(){
int n=3;
double data[]={2, 0, 0,  0, 3, 0, 0 , 0, 5};
double *vector=malloc(sizeof(double)*n);
double *matrix=malloc(sizeof(double)*n*n);
eigen(data,  matrix,vector );

gsl_matrix * m = gsl_matrix_alloc (3, 3);
gsl_matrix * inverse = gsl_matrix_alloc (3, 3);
gsl_permutation * perm = gsl_permutation_alloc (3);

 /*printf("Eigenvalue:\n%g\t%g\t%g\nEigenvector:\n", vector[0], vector[1], vector[2]);
    for(int i=0; i<3; i++) {
        for(int j=0; j<3; j++) {
            printf("(%d) %f\t",i*3+j, matrix[i*3+j]);
        }
        printf("\n");
    }
*/

return 0;
}
