#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include "minc2.h"
#include "hdf5.h"
#include "volume_io.h"
#include "minc_helper.h"
#include "anisotropic_diffusion.h"
#include "3x3-C/dsyevh3.h"
#include "3x3-C/dsyevd3.h"
#include "pthread.h"
#include <unistd.h>


//#include <fftw3.h>
//#include <rfftw.h>
#define TENSOR 1


void printMat(float *mat){
	for(int y=0; y<3; y++ ){
		for(int x=0; x< 3; x++){
			printf("%f\t", mat[y*3+x]);
		}
	printf("\n");
	}	
printf("\n");
}
void printMat_d(double *mat){
	for(int y=0; y<3; y++ ){
		for(int x=0; x< 3; x++){
			printf("%f\t", mat[y*3+x]);
		}
	printf("\n");
	}	
printf("\n");
}
int isin(int number, int* array, int dim){

    for(int i=0; i<dim; i++) if( number==array[i]) return 1;

return 0;

}

void print_stats(double* data, int n){
    float total=0, min=data[0], max=data[0], mean=0, sd, sum2=0;
    for(int  i=0; i < n; i++){
        total += data[i];
        sum2 +=  (data[i] * data[i]);
        if(min > data[i]) min=data[i];
        if(max < data[i]) max=data[i];
    }
    mean = total / n;
    sd=sqrt((sum2 - (total*total)/n)/n);
    
    printf("Min: %3.3f\tMax: %3.3f\tMean: %3.3f\tStdDev: %3.3f", min, max, mean, sd);
}

int mult(double *vol1, double *vol2, int o_rows, int i_cols, int i_rows, int o_cols, double* output ){
    int index1, index2, i, m, n;
    float sum=0;
        
    if(i_cols != i_rows) {
        printf("Error: Invalid matrix dimensions.\n"); 
        return(1);
    }

    for(int m=0; m<o_rows; m++){//for outer rows
        for(int n=0; n<o_cols; n++){
            sum=0;
            for(int i=0; i<i_rows; i++){
                index1=m*i_cols+i;
                index2=i*o_cols+n;
                //printf("%f(%d) x %f(%d) + ", vol1[index1], index1, vol2[index2], index2 );
                sum += vol1[index1]*vol2[index2];
            }
            //printf("=%f\t%d\n", sum, m*o_cols+n );
            output[m*o_cols+n]=sum;
        }
    }
    //printMat(output);
    return(0);
}

int  sort_eigenvalues(double* eigenvalues, int* l_eig, int* m_eig, int* s_eig){
    int E0gtE1=eigenvalues[0] >= eigenvalues[1];
    int E0gtE2=eigenvalues[0] >= eigenvalues[2];
    int E1gtE2=eigenvalues[1] >= eigenvalues[2];
    //printf("\t%d %d %d\n", E0gtE1,E0gtE2,E1gtE2);
    if( E0gtE1 && E0gtE2 && E1gtE2 ){//E0 > E1 > E2
        *l_eig=0;
        *m_eig=1;
        *s_eig=2; 
    } 
    else if( E0gtE1 && E0gtE2 && E1gtE2 == FALSE){//E0 > E2 > E1
        *l_eig=0;
        *m_eig=2;
        *s_eig=1;                
    }
    else if( E0gtE1 == FALSE && E0gtE2 ==FALSE &&  E1gtE2 == FALSE ){ //E2 > E1 > E0 
        *l_eig=2;
        *m_eig=1;
        *s_eig=0;                        
    }
    else if (E0gtE1 && E0gtE2 ==FALSE && E1gtE2==FALSE ){ //E2 > E0 > E1
        *l_eig=2;
        *m_eig=0;
        *s_eig=1;  
    }
    else if (E0gtE1 ==FALSE && E0gtE2 && E1gtE2 ) { //E1 > E0 > E2
        *l_eig=1;
        *m_eig=0;
        *s_eig=2;  
    }
    else if(E0gtE1==FALSE && E0gtE2==FALSE && E1gtE2){ //E1 > E2 > E0
        *l_eig=1;
        *m_eig=2;
        *s_eig=0;  
    } else pexit("Problem with logical operations when sorting eigenvalues\n", "", 1);
    return(0);
}
double diffusivity_function1(float alpha, float lambda2){
    //     1
    //-------------
    //  1 + alpha^2 / lambda^2
    return(1/( 1 + ((alpha*alpha) / lambda2)));
}

float inverse3x3(double* m, double* out){
    //[ 0, 1, 2 ]
    //[ 3, 4, 5 ]
    //[ 6, 7, 8 ]
    //Get minors
    double d0 = m[4] * m[8] - m[7] * m[5];
    double d1 = m[3] * m[8] - m[6] * m[5];
    double d2 = m[3] * m[7] - m[6] * m[4];
    double d3 = m[1] * m[8] - m[7] * m[2];
    double d4 = m[0] * m[8] - m[6] * m[2];
    double d5 = m[0] * m[7] - m[6] * m[1];
    double d6 = m[1] * m[5] - m[4] * m[2];
    double d7 = m[0] * m[5] - m[3] * m[2];
    double d8 = m[0] * m[4] - m[3] * m[1];
    double det = m[0] * d0 - m[3] * d3 + m[6] * d6;
    //printf("%f %f\n", d0, m[0]);
    //printf("%f = %f - %f + %f\n", det, m[0] * d0, m[3] * d3, m[6] * d6);
    if(det != 0){
        //[ 0, -3, 6 ]
        //[ -1, 4, -7 ]
        //[ 2, -5, 8 ]
        out[0]= d0 / det;
        out[1]=-d3 / det;
        out[2]= d6 / det;
        out[3]=-d1 / det;
        out[4]= d4 / det;
        out[5]=-d7 / det;
        out[6]= d2 / det;
        out[7]=-d5 / det;
        out[8]= d8 / det;

    }
    return(det);
}


int gradient_threaded(void* args){
    double* image= ((struct grad_args*) args)->image; 
    float* dx=((struct grad_args*) args)->dx; 
    float* dy=((struct grad_args*) args)->dy; 
    float* dz=((struct grad_args*) args)->dz; 
    int thread=((struct grad_args*) args)->thread;
    int nthreads=((struct grad_args*) args)->nthreads;
    int* count=((struct grad_args*) args)->count;
    
    int zmax=count[0];
    int ymax=count[1];
    int xmax=count[2];
    int indexp, indexm, index;
    int thread_step=nthreads;
    int m=1, p=1;
    for(int z=thread; z<zmax; z+=thread_step){
        for(int y=0; y < ymax; y++){
            for(int x=0; x < xmax; x++){
                index=z*ymax*xmax+y*xmax+x;
                if(z-1 < 0 || z+1 >= zmax || y-1 < 0 || y+1 >= ymax || x-1 < 0 || x+1 >= xmax ){
                    dx[index]=dy[index]=dz[index]=0;
                }
                else{
                    indexp=z*ymax*xmax+y*xmax+x+p;
                    indexm=z*ymax*xmax+y*xmax+x-m;
                    dx[index]=(float)image[indexp] - image[indexm]; //gradient in x direction
                    indexp=z*ymax*xmax+(y+p)*xmax+x;
                    indexm=z*ymax*xmax+(y-m)*xmax+x;
                    dy[index]=(float) image[indexp] - image[indexm]; //gradient in y direction
                    indexp=(z+p)*ymax*xmax+y*xmax+x;
                    indexm=(z-m)*ymax*xmax+y*xmax+x;
                    dz[index]=(float) image[indexp] - image[indexm]; //gradient in z direction
                }
            }
        }
    }
    return(0);

}
int gradient_multithread(double* image, float* dx, float* dy, float *dz, int* count, int nthreads){
    int rc;
    pthread_t threads[nthreads];
    struct grad_args thread_args[nthreads];
    for(int t=0; t< nthreads; t++){
        thread_args[t].image=image;
        thread_args[t].dx=dx;
        thread_args[t].dy=dy;
        thread_args[t].dz=dz;
        thread_args[t].thread=t;
        thread_args[t].nthreads=nthreads;
        thread_args[t].count=count;
        rc= pthread_create(&threads[t], NULL, gradient_threaded, (void *) &thread_args[t] ); //gradient(image, img_dx, img_dy, img_dz);
    }
    for(int t=0; t<nthreads;t++){
        rc = pthread_join(threads[t], NULL);
    }
    return(0);
}

int divergence_threaded(void* args){
    float* image= ((struct div_args*) args)->image; 
    int* count= ((struct div_args*) args)->count;
    float* dx=((struct div_args*) args)->dx; 
    float* dy=((struct div_args*) args)->dy; 
    float* dz=((struct div_args*) args)->dz; 
    int thread=((struct div_args*) args)->thread;
    int nthreads=((struct div_args*) args)->nthreads;

    int zmax=count[0];
    int ymax=count[1];
    int xmax=count[2];
    int indexp, indexm, index;
    float sum;
    int m=1, p=1;
    int thread_step=nthreads;
    //Div(I)= dI/dx + dI/dy + dI/dz
    for(int z=thread; z<zmax; z+=thread_step){
        for(int y=0; y < ymax; y++){
            for(int x=0; x < xmax; x++){
                index=z*ymax*xmax+y*xmax+x;
                if(z-1 < 0 || z+1 >= zmax || y-1 < 0 || y+1 >= ymax || x-1 < 0 || x+1 >= xmax ){
                    dx[index]=dy[index]=dz[index]=0;
                }
                else{
                    sum=0;
                    //finite difference in x dimension
                    indexp=z*ymax*xmax+y*xmax+x+p;
                    indexm=z*ymax*xmax+y*xmax+x-m;
                    sum += dx[indexp] - dx[indexm]; 
                    //finite difference in y dimension                   
                    indexp=z*ymax*xmax+(y+p)*xmax+x;
                    indexm=z*ymax*xmax+(y-m)*xmax+x;
                    sum += dy[indexp] - dy[indexm]; 
                    //finite difference in z dimension                   
                    indexp=(z+p)*ymax*xmax+y*xmax+x;
                    indexm=(z-m)*ymax*xmax+y*xmax+x;
                    sum += dz[indexp] - dz[indexm];
                    image[index]=sum;
                }
            }
        }
    }
    return(0);
}
int divergence_multithread(float* image,  int* count,float* dx, float* dy, float *dz, int nthreads){
    int rc;
    pthread_t threads[nthreads];
    struct div_args thread_args[nthreads];
    for(int t=0; t< nthreads; t++){
        thread_args[t].image=image;
        thread_args[t].count=count;
        thread_args[t].dx=dx;
        thread_args[t].dy=dy;
        thread_args[t].dz=dz;
        thread_args[t].thread=t;
        thread_args[t].nthreads=nthreads;
        rc= pthread_create(&threads[t], NULL, divergence_threaded, (void *) &thread_args[t] ); //gradient(image, img_dx, img_dy, img_dz);
        if(rc!=0) pexit("Error: creating thread","" , 1);
    }
    for(int t=0; t<nthreads;t++){
        rc = pthread_join(threads[t], NULL);
        if(rc!=0) pexit("Error: joining thread", "", 1);
    }
    return(0);
}

int update_img_threaded(void* args){
    float* temp_dx= ((struct update_img_args*) args)->temp_dx ; 
    float* temp_dy= ((struct update_img_args*) args)->temp_dy ; 
    float* temp_dz= ((struct update_img_args*) args)->temp_dz ;
    float* img_dx= ((struct update_img_args*) args)->img_dx;
    float* img_dy= ((struct update_img_args*) args)->img_dy;
    float* img_dz= ((struct update_img_args*) args)->img_dz;
    float** big_tensor= ((struct update_img_args*) args)->big_tensor ; 
    int* count= ((struct update_img_args*) args)->count ;
    int thread= ((struct update_img_args*) args)->thread;
    int nthreads= ((struct update_img_args*) args)->nthreads;
    int thread_step=nthreads;
    float* tensor_ptr;
    int zmax=count[0];
    int ymax=count[1];
    int xmax=count[2];
    for(int z=thread; z<zmax; z+=thread_step) 
        for(int y=0; y<ymax; y++) 
            for(int x=0; x<xmax; x++) {
                int index=z*ymax*xmax+y*xmax+x;
                tensor_ptr=big_tensor[index];

                //#################################
                //# EQ:  Tensor * Grad( Image_i ) #
                //# T[0, 1, 2]   [dI/dx]          #
                //#  [3, 4, 5] X [dI/dy]          #  
                //#  [6, 7, 8]   [dI/dz]          #
                //#################################
    #if TENSOR
                temp_dx[index]=(float) tensor_ptr[0]*img_dx[index]+tensor_ptr[1]*img_dy[index]+tensor_ptr[2]*img_dz[index];
                temp_dy[index]=(float) tensor_ptr[3]*img_dx[index]+tensor_ptr[4]*img_dy[index]+tensor_ptr[5]*img_dz[index];
                temp_dz[index]=(float) tensor_ptr[6]*img_dx[index]+tensor_ptr[7]*img_dy[index]+tensor_ptr[8]*img_dz[index];
    #else
                temp_dx[index]=(float) img_dx[index];
                temp_dy[index]=(float) img_dy[index];
                temp_dz[index]=(float) img_dz[index];

    #endif
    }
    return(0);
}

int update_img_multithread(float* temp_dx,float* temp_dy,float* temp_dz,float* img_dx, float* img_dy, float* img_dz, float** big_tensor,  int* count, int nthreads){
    int rc;
    pthread_t threads[nthreads];
    struct update_img_args thread_args[nthreads];
    for(int t=0; t< nthreads; t++){
        thread_args[t].temp_dz=temp_dz;
        thread_args[t].temp_dy=temp_dy;
        thread_args[t].temp_dx=temp_dx;
        thread_args[t].img_dz=img_dz;
        thread_args[t].img_dy=img_dy;
        thread_args[t].img_dx=img_dx;
        thread_args[t].big_tensor=big_tensor;
        thread_args[t].count=count;
        thread_args[t].thread=t;
        thread_args[t].nthreads=nthreads;
        rc= pthread_create(&threads[t], NULL, update_img_threaded, (void *) &thread_args[t] ); //gradient(image, img_dx, img_dy, img_dz);
        if(rc!=0) pexit("Error: creating thread","" , 1);
    }
    for(int t=0; t<nthreads;t++){
        rc = pthread_join(threads[t], NULL);
        if(rc!=0) pexit("Error: joining thread", "", 1);
    }
    return(0);
}

int diffuse_threaded(void* args){
    double* image= ((struct diffusion_args*) args)->image; 
    float* temp_image=((struct diffusion_args*) args)->temp_image;
    double scale_factor=((struct diffusion_args*) args)->scale_factor;
    int nthreads=((struct diffusion_args*) args)->nthreads;
    int thread=((struct diffusion_args*) args)->thread;
    int* count=((struct diffusion_args*) args)->count;
    int thread_step=nthreads;
    int zmax=count[0];
    int ymax=count[1];
    int xmax=count[2];
    
    for(int z=thread; z<zmax; z+=thread_step)
        for(int y=0; y<ymax; y++) 
            for(int x=0; x<xmax; x++) {
               int index=z*ymax*xmax+y*xmax+x;
                //Need to figure out a good number for delta t
                image[index] += (scale_factor * (double) temp_image[index]);
    }
    return(0);
}
int diffuse_multithread(double* image, float* temp_image, double scale_factor, int* count, int nthreads  ){
    int rc;
    pthread_t threads[nthreads];
    struct diffusion_args thread_args[nthreads];
    for(int t=0; t< nthreads; t++){
        thread_args[t].temp_image=temp_image;
        thread_args[t].image=image;
        thread_args[t].scale_factor=scale_factor;
        thread_args[t].thread=t;
        thread_args[t].nthreads=nthreads;
        thread_args[t].count=count;
        rc= pthread_create(&threads[t], NULL, diffuse_threaded, (void *) &thread_args[t] ); //gradient(image, img_dx, img_dy, img_dz);
        if(rc!=0) pexit("Error: creating thread","" , 1);
    }
    for(int t=0; t<nthreads;t++){
        rc = pthread_join(threads[t], NULL);
        if(rc!=0) pexit("Error: joining thread", "", 1);
    }
    return(0);
}




int structure_tensor_threaded(void* args){
    float** big_tensor=((struct structure_tensor_args*) args)->big_tensor;
    float* dx=((struct structure_tensor_args*) args)->dx;
    float* dy=((struct structure_tensor_args*) args)->dy;
    float* dz=((struct structure_tensor_args*) args)->dz;
    int* count=((struct structure_tensor_args*) args)->count;
    double lambda=((struct structure_tensor_args*) args)->lambda;
    int thread=((struct structure_tensor_args*) args)->thread;
    int nthreads=((struct structure_tensor_args*) args)->nthreads;
    int thread_step=nthreads;
	double eigenvectors[3][3];
    double *eigenvectors_ptr=&(eigenvectors[0][0]);
    double eigenvalues[3]; //*eigenvalues=malloc(sizeof(double) * 3);
    double tensor[3][3];
    float *tensor_ptr;
    double* tensor_dptr=&(tensor[0][0]);
    double *D=calloc(9, sizeof(*D));//Array for storing diagonalized eigen vectors
	double *eigenvectors_inverse=malloc(sizeof(*eigenvectors)*9);
    double *temp=malloc(sizeof(*temp)*9);
    char tempoutputfilename[100];
    double lambda2=lambda*lambda;
    int zmax=count[0];
    int ymax=count[1];
    int xmax=count[2];
    float dxdx; 
    float dxdy; 
    float dxdz;  
    float dydz; 
    float dydy;  
    float dzdz;
    int trigger=0, l_eig, m_eig, s_eig; 



    for(int z=thread; z<zmax; z+=thread_step) for(int y=0; y<ymax; y++) for(int x=0; x<xmax; x++) {
        int index=z*ymax*xmax+y*xmax+x;
        tensor_ptr=big_tensor[index];
        //Calculate structural tensor.
        dxdx=dx[index]* dx[index];  
        dydy=dy[index]* dy[index]; 
        dzdz=dz[index]* dz[index];
        dxdy=dx[index]* dy[index];  
        dxdz=dx[index]* dz[index]; 
        dydz=dy[index]* dz[index]; 
        //######################
        //# Tensor             # 
        //# [dxdx, dxdy, dxdz] #
        //# [dydx, dydy, dydz] #
        //# [dzdx, dzdy, dzdz] #
        //######################
        tensor[0][0]=(double) dxdx; 
        tensor[0][1]=(double) dxdy;
        tensor[0][2]=(double) dxdz;
        tensor[1][0]=(double) dxdy;
        tensor[1][1]=(double) dydy;
        tensor[1][2]=(double) dydz;
        tensor[2][0]=(double) dxdz;
        tensor[2][1]=(double) dydz;
        tensor[2][2]=(double) dzdz;
        float tensor_sum=0;
        for(int i=0; i < 3; i++) for(int j=0; j<3; j++) {
            tensor_sum += (float) tensor[i][j];
        }
        //Find eigenvalues and eigenvectors for structure tensor 
        if( tensor_sum == 0){
            for(int i=0; i < 9; i++){ 
                if(i==0 || i==4 || i==8) tensor_ptr[i]=1;
                else tensor_ptr[i]=0;
                //printMat(tensor_ptr);
            } 
        }
        else{
            dsyevd3(tensor, eigenvectors, eigenvalues);
            //dsyevh3(tensor, eigenvectors, eigenvalues); 
            //find largest eigenvalue
            sort_eigenvalues(eigenvalues, &l_eig,  &m_eig, &s_eig );
            eigenvalues[l_eig]=diffusivity_function1(eigenvalues[l_eig], lambda2);
            eigenvalues[m_eig]=diffusivity_function1(eigenvalues[m_eig], lambda2);
            eigenvalues[s_eig]=diffusivity_function1(eigenvalues[s_eig], lambda2);
            D[0]=eigenvalues[0]; //D[ e1, 0, 0]
            D[4]=eigenvalues[1]; // [ 0, e2, 0]
            D[8]=eigenvalues[2]; // [ 0, 0, e3]
            //########################################################################
            //# Reconstruct structure tensor with modified eigen values              #
            //# New Tensor = Eigenvectors x New Eigenvalues x  Inverse(Eigenvectors) #
            //########################################################################
            float det=inverse3x3( eigenvectors_ptr, eigenvectors_inverse);
                //printf("%f\n", det );
                //printMat_d(tensor);
                //printMat_d(eigenvectors_ptr);
            if( det != 0 ){
                mult( eigenvectors_ptr, D, 3, 3, 3, 3, temp);
                mult(temp, eigenvectors_inverse, 3, 3, 3, 3, tensor_dptr);
                for(int i=0; i<9; i++) tensor_ptr[i] = (float) tensor_dptr[i];
            }
        }
    }
    free(eigenvectors_inverse);
    return(0);
}
int structure_tensor_multithread(float** big_tensor, float* dx, float* dy, float* dz,int* count, double lambda, int nthreads ){ 
    int rc;
    pthread_t threads[nthreads];
    struct structure_tensor_args thread_args[nthreads];
    for(int t=0; t< nthreads; t++){
        thread_args[t].big_tensor=big_tensor;
        thread_args[t].dx=dx;
        thread_args[t].dy=dy;
        thread_args[t].dz=dz;
        thread_args[t].count=count;
        thread_args[t].lambda=lambda;
        thread_args[t].thread=t;
        thread_args[t].nthreads=nthreads;
        rc= pthread_create(&threads[t], NULL, structure_tensor_threaded, (void *) &thread_args[t] ); //gradient(image, img_dx, img_dy, img_dz);
        if(rc!=0) pexit("Error: creating thread","" , 1);
    }
    for(int t=0; t<nthreads;t++){
        rc = pthread_join(threads[t], NULL);
        if(rc!=0) pexit("Error: joining thread", "", 1);
    }
    return(0);

    return(0);
}
    

/***********************************
Name:   smooth
Inputs: mask - mask of discrete regions
        image - image to be smoothed
        fwhm - maximum number of iterations, determines FWHM


***********************************/
int smooth(double* mask, double* image, int* count, double fwhm, double dt, double step, double lambda,  int nthreads, int VERBOSE){

    double totalTime=fwhm*fwhm / (16*log(2));
    int maxIterations= (int) round( totalTime/ dt);
    dt=totalTime/maxIterations;
    int zmax=count[0];
    int ymax=count[1];
    int xmax=count[2];
    int nvox=zmax*ymax*xmax;
    int index, row;
    int n[3]={zmax, ymax, xmax};
    char tempoutputfilename[100];
    float *dx=calloc(nvox, sizeof(float));//Gradient for mask 
    float *dy=calloc(nvox, sizeof(float)); 
    float *dz=calloc(nvox, sizeof(float));
    float *img_dx=calloc(nvox, sizeof(float));//Gradient for image at iteration i
    float *img_dy=calloc(nvox, sizeof(float)); 
    float *img_dz=calloc(nvox, sizeof(float)); 
    float *temp_dx=calloc(nvox, sizeof(float));//Temporary gradient, before divergence 
    float *temp_dy=calloc(nvox, sizeof(float)); 
    float *temp_dz=calloc(nvox, sizeof(float)); 
    float *temp_image=calloc(nvox, sizeof(float));
    float** big_tensor=malloc(nvox*sizeof(*big_tensor));
    
    for(int x=0; x < nvox; x++) big_tensor[x]=calloc(9,sizeof(**big_tensor));

    double dx2 =step*step;
    //Calculate gradient of structural image
    for(int iterations=1; iterations <= maxIterations; iterations++){
        if(iterations==1){
            gradient_multithread(mask, dx, dy, dz, n, nthreads);
            //writeVolume("dz.mnc", dz, image->start,  image->step, image->wcount, MI_TYPE_FLOAT  );
            structure_tensor_multithread(big_tensor, dx, dy, dz, n, lambda, nthreads );
        }
        //EQ: Image_i+1 =  Image_i + (delta_t/(4*delta_x^2)) * Div(Tensor * Grad( Image_i ) )
        gradient_multithread(image, img_dx, img_dy, img_dz, n, nthreads);

        //sprintf(tempoutputfilename,"img_dx_%d.mnc", iterations );
        //writeVolume(tempoutputfilename, img_dx, mask->start,  mask->step, mask->wcount, MI_TYPE_FLOAT  );
        update_img_multithread( temp_dx, temp_dy, temp_dz, img_dx, img_dy, img_dz, big_tensor, n, nthreads);
        divergence_multithread(temp_image, n, temp_dx, temp_dy, temp_dz, nthreads);
        diffuse_multithread(image, temp_image, 0.25*dt/dx2, n , nthreads  );

        //sprintf(tempoutputfilename,"outputfile_%d.mnc", iterations );
        //writeVolume(tempoutputfilename, image->data, mask->start,  mask->step, mask->wcount, MI_TYPE_DOUBLE  );

        if(VERBOSE){
            printf("%d:\t", iterations);
            print_stats(image, zmax*ymax*xmax );
            printf("\n");
        }
   }
    
    free(dx);
    free(dy);
    free(dz);
    free(img_dx);
    free(img_dy);
    free(img_dz);
    free(temp_dx);
    free(temp_dy);
    free(temp_dz);
    free(temp_image);
    for(int x=0; x < nvox; x++) free(big_tensor[x]);
    free(big_tensor);
    return(0);
}

