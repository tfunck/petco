#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <fftw3.h>
#include "minc2.h"
#include "hdf5.h"
#include "minc_helper.h"


void swap(double* img, int i1, int i2){
    double tmp;
    tmp=img[i1];
    img[i1]=img[i2];
    img[i2]=tmp;
}


int fftshift(double* img, int* n){
    double xmid=n[2]/2;
    double ymid=n[1]/2;
    double zmid=n[0]/2;
    int xc =ceil(xmid);
    int xf=floor(xmid); 
    int yc =ceil(ymid);
    int yf=floor(ymid); 
    int zc =ceil(zmid);
    int zf=floor(zmid); 

    for(int z=0; z< n[0]; z++)
        for(int y=0; y< n[1]; y++)
            for(int x=0; x< xf; x++){
                int i=z*n[1]*n[2] + y*n[2] + x;
                int x2=z*n[1]*n[2] + y*n[2] + (xc+x);    
                swap(img, i, x2);
            }

    for(int z=0; z< n[0]; z++)
        for(int y=0; y< yf; y++)
            for(int x=0; x< n[2]; x++){
                int i=z*n[1]*n[2] + y*n[2] + x;
                int y2=z*n[1]*n[2] + (yc+y)*n[2] + x;  
                swap(img, i, y2);
            }
    for(int z=0; z< zf; z++)
        for(int y=0; y< n[1]; y++)
            for(int x=0; x< n[2]; x++){
                int i=z*n[1]*n[2] + y*n[2] + x;
                int z2=(zc+z)*n[1]*n[2] + y*n[2] + x;         
                swap(img, i, z2);
            }

    return(0);
}

int filter(fftw_complex* im1, fftw_complex* im2, int n){
    for(int i=0; i < n; i++){
        im1[i][0]=im1[i][0] * im2[i][0] - im1[i][1] * im2[i][1];
        im1[i][1]=im1[i][0] * im2[i][1] + im1[i][1] * im2[i][0];
    }
    return(0);
}


double* pad_like(double* img, int* N, int* M){
    for(int i=0; i<3; i++) if(N[i] > M[i]) pexit("Image to pad like is smaller than input image", "",1);
    double* img_pad=calloc(M[0] * M[1] * M[2], sizeof(double*));
    int mid[3]={ (int) round(M[0]/2),   (int) round(M[1]/2),   (int) round(M[2]/2)  };
    int pad_dim[3] = {(N[0]-1)/2, (N[1]-1)/2, (N[2]-1)/2};
    int start[3]= {mid[0]-pad_dim[0], mid[1]-pad_dim[1],mid[2]-pad_dim[2]};
    int end[3]= {mid[0]+pad_dim[0], mid[1]+pad_dim[1],mid[2]+pad_dim[2]};

    for(int z=start[0]; z < end[0] ; z++ )
        for(int y=start[1]; y < end[1]; y++ )
            for(int x=start[2]; x < end[2]; x++ ){
                int i2=z*M[1]*M[2] + y * M[2] + x;
                int i1=(z-start[0])*N[1]*N[2] + (y-start[1]) * N[2] + (x-start[2]);
                img_pad[i2]=img[i1];
            }
    
    return(img_pad);
}


double* pad(double* img, int* N, int* M, int* N_pad){
    int pad_dim[3] = {(M[0]-1)/2, (M[1]-1)/2, (M[2]-1)/2};
    for(int i=0; i<3; i++) N_pad[i]= N[i]+2*pad_dim[i];
   
    double* img_pad=calloc(N_pad[0] * N_pad[1] * N_pad[2], sizeof(*img_pad));
    //printf("%d %d %d\n", M[0], M[1], M[2] );
    //printf("%d %d %d\n", N[0], N[1], N[2] );
    //printf("%d %d %d\n", N_pad[0], N_pad[1], N_pad[2]);
    for(int z=pad_dim[0]; z < pad_dim[0]+N[0]; z++ )
        for(int y=pad_dim[1]; y < pad_dim[1]+N[1]; y++ )
            for(int x=pad_dim[2]; x < pad_dim[2]+N[2]; x++ ){
                int i2=z*N_pad[1]*N_pad[2] +  y * N_pad[2] + x;
                int i1=(z-pad_dim[0])*N[2]*N[1] + (y-pad_dim[1]) * N[2] + (x-pad_dim[2]);
                img_pad[i2]=img[i1];
            }
    
    return(img_pad);
}


int remove_pad(double* img, double* img_pad, int* N, int* M){
    int pad_dim[3] = {(M[0]-1)/2, (M[1]-1)/2, (M[2]-1)/2};
    int N_pad[3];
    for(int i=0; i<3; i++) N_pad[i]= N[i]+2*pad_dim[i];

    for(int z=pad_dim[0]; z < pad_dim[0]+N[0]; z++ )
        for(int y=pad_dim[1]; y < pad_dim[1]+N[1]; y++ )
            for(int x=pad_dim[2]; x < pad_dim[2]+N[2]; x++ ){
                int i2=z*N_pad[1]*N_pad[2] +  y * N_pad[2] + x;
                int i1=(z-pad_dim[0])*N[2]*N[1] + (y-pad_dim[1]) * N[2] + (x-pad_dim[2]);
                img[i1]=img_pad[i2];
            }
    
    return(0);
}

double gaussian(double a, double b, double c, double x){ return(a * exp(- b * (x-c)*(x-c) ));}

double* create_filter(double fwhm, int* N, double* step){
    /*Create Gaussian kernel*/
    double pi=3.14159265359;
    double t=0.000001; //tolerance level for the precision of the approximation of the gaussian
    double sd=fwhm/(2*sqrt(2*log(2)));
    double a=1/(sd*sqrt(2*pi));
    double b=1/(2*sd*sd);
    double min = sd * sqrt( -2* log(t/a)  );
    int minima[3]= { ceil(min/fabsf(step[0])),   ceil(min/fabsf(step[1])),   ceil(min/fabsf(step[2])) }; //minimal distance in voxels needed to achieve desired precision 
    int length[3]={minima[0]*2+1, minima[1]*2+1, minima[2]*2+1 };
    double start[3]={};
    double** gauss1d=malloc(sizeof(*gauss1d)* 3);
    double* out = calloc(length[0]*length[1]*length[2], sizeof(*out) );
    for(int i=0; i<3; i++){ 
        N[i]=length[i];//Pass length of kernel so that it will be accessible outside of this function. 
        gauss1d[i]=malloc(length[i]*sizeof(**gauss1d));
        for(int k=0; k< minima[i]; k++){
            int k1=minima[i]-k-1;
            int k2=minima[i]+1+k;
            //gaussian is symmetric so only need to calculate one side
            gauss1d[i][k1]=gauss1d[i][k2]= gaussian(a, b, 0, (k+1)*step[i] );

        }
        gauss1d[i][ minima[i]  ]=gaussian(a, b, 0, 0 );
        //for(int k=0; k< length[i]; k++) printf("%d %f\n", i, gauss1d[i][k]);
    }

    for(int z=0; z<length[0]; z++ )
        for(int y=0; y<length[1]; y++ )
            for(int x=0; x<length[2]; x++ ){
                int index=z*length[1]*length[2]+y*length[2]+x;
                out[index]=gauss1d[0][z]*gauss1d[1][y]*gauss1d[2][x];
            } 
    
    for(int i=0; i<3; i++) free(gauss1d[i]);
    free(gauss1d);
    return(out);
}


int inverse(double *output, fftw_complex* signal, fftw_complex* SIGNAL, int* N){
    fftw_plan p;
        
    p = fftw_plan_dft_3d(N[0], N[1], N[2], SIGNAL, signal, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(p); /* repeat as needed */

    for(int i=0; i<N[0]*N[1]*N[2]; i++) output[i]=signal[i][0]/(N[0]*N[1]*N[2]);
     
    fftw_destroy_plan(p);
    return(0);
}

int forward(double *input, fftw_complex* signal, fftw_complex* SIGNAL, int* N){
    fftw_plan p;
        
    for(int i=0; i<N[0]*N[1]*N[2]; i++) signal[i][0]=input[i];
    
    p = fftw_plan_dft_3d(N[0], N[1], N[2], signal, SIGNAL, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p); /* repeat as needed */
 
    fftw_destroy_plan(p);
    return(0); 
}



int gaussian_filter( double* input, double fwhm, unsigned int* N, double* step){
    int M[3];
    int N_pad[3];
    int Mmax;
    fftw_complex *signal, *SIGNAL, *f, *F;
    double* fgaussian=create_filter(fwhm, M, step); 
    double* input_pad=pad(input, N, M, N_pad);
    double* fgaussian_pad=pad_like(fgaussian, M, N_pad);

    int Nmax=N_pad[0] * N_pad[1] * N_pad[2];
    Mmax=M[0]*M[1]*M[2];
   
    signal = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nmax);
    SIGNAL = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nmax);
    f = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nmax);
    F = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nmax);
    forward(input_pad, signal, SIGNAL, N_pad);  //Forward FFT signal
    forward(fgaussian_pad, f, F, N_pad);        //Forward FFT filter
    filter(SIGNAL, F, Nmax );
    inverse(input_pad, signal, SIGNAL, N_pad);  //Inverse filter

    fftshift(input_pad, N_pad);
    remove_pad(input, input_pad, N, M );

    free(input_pad);
    free(fgaussian_pad);
    free(fgaussian);
    fftw_free(f);
    fftw_free(F);
    fftw_free(signal);
    fftw_free(SIGNAL);
    return(0);
}


int smooth_image(data* image, double fwhm){
    int i=gaussian_filter(image->data,  fwhm, image->count, image->step);
    return(i);
}

/*int main(int argc, char** argv){
    if(strcmp(argv[1], "-help") == 0 ){
        useage();
    }

    data image;
    unsigned char* output_fn=argv[ 2]; //printf("%s\n", output_fn); exit(1);
    mihandle_t img;


    image.filename=argv[1];
    image.data = (double*) readVolume( &image , 1, MI_TYPE_DOUBLE);

    gaussian_filter(&image,  2.5 );
    writeVolume(output_fn, image.data, image.start, image.step, image.wcount, MI_TYPE_DOUBLE);
  
 
    return(0);
}


void useage(){
    printf("Useage: mask1.mnc mask2.mnc ... maskn.mnc out.mnc\n");
    printf("Purpose: Concatenate a series of masks and give them unique labels.\n");
    exit(1);
}*/
