#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <omp.h>
#include "minc2.h"
#include "hdf5.h"
#include "minc_helper.h"

void useage();

int main(int argc, char** argv){
    if(argc != 4 || strcmp(argv[1], "-help") == 0 ){
        useage();
    }

    data image, mask;

    //readVolume(data* volume, int read_hyperslab /*1 for read, o.w. don't read*/, mitype_t dtype)
    unsigned char* image_vol;
    image.filename=argv[1];
    double new_step_size=atof(argv[2]);
    char* output_filename=argv[3];
    float* array, *new, min, max;
    int tmax, zmax, ymax, xmax;
    double wmax_x, wmax_y, wmax_z;
    double* new_step;
    double* new_starts;
    misize_t *wstarts, *wcount;
    mihandle_t img;
    float dxdydz; 

    array = (float*) readVolume(&image, 1, MI_TYPE_FLOAT);
    
    dxdydz=image.zstep*image.ystep*image.xstep;
    
    new_step=malloc(image.ndim*sizeof(*new_step));
    new_starts=malloc(image.ndim*sizeof(*new_starts));
    wstarts=malloc(image.ndim*sizeof(*wstarts));
    wcount=malloc(image.ndim*sizeof(*wcount));


    zmax=image.zmax;
    ymax=image.ymax;
    xmax=image.xmax;
    if(image.tmax == 0) tmax=1;
    else tmax=image.tmax;

    wmax_z=zmax*image.zstep + image.zstart;
    wmax_y=ymax*image.ystep + image.ystart;
    wmax_x=xmax*image.xstep + image.xstart;

    int s=4-image.ndim;
    int zi=1-s;
    int yi=2-s;
    int xi=3-s; 
    if( image.ndim==4){ 
        wcount[0]=image.tmax; 
        new_step[0]=1;
        new_starts[0]=0;
        wstarts[0];
    }

    wstarts[zi]=wstarts[yi]=wstarts[xi]=0;
    wcount[zi]=(misize_t) ceil((wmax_z - image.zstart)/new_step_size) ; 
    wcount[yi]=(misize_t) ceil((wmax_y - image.ystart)/new_step_size) ; 
    wcount[xi]=(misize_t) ceil((wmax_x - image.xstart)/new_step_size) ; 
    new_step[zi]=new_step_size; 
    new_step[yi]=new_step_size; 
    new_step[xi]=new_step_size;
    new_starts[zi]=image.zstart;
    new_starts[yi]=image.ystart;
    new_starts[xi]=image.xstart;


    new=calloc( wcount[zi] * wcount[yi] * wcount[xi], sizeof(*new));
    createVolume(output_filename, image.ndim, wcount, new_step, new_starts);
    miopen_volume(output_filename, MI2_OPEN_RDWR, &img);
    if(image.ndim==4) wcount[0]=1;
    min=max=0; //new[0];
    for( int t=0; t<tmax; t++){
        wstarts[0]=t;
        read_frame(&image, t, array, MI_TYPE_FLOAT);
        for(unsigned int i=0; i< wcount[zi] *wcount[yi]*wcount[xi]; i++) new[i]=0;
        //#pragma omp parallel  
        for(unsigned int z=0; z< zmax; z++){
            for(unsigned int y=0; y< ymax; y++){
                for(unsigned int x=0; x < xmax; x++){
                    int i = (int) floor((z*image.zstep)/new_step[zi]);
                    int j = (int) floor((y*image.ystep)/new_step[yi]); 
                    int k = (int) floor((x*image.xstep)/new_step[xi]); 
                    unsigned int index2=i*wcount[yi]*wcount[xi]+j*wcount[xi]+k;
                    unsigned int index=z*ymax*xmax+y*xmax+x;
                    //if(array[index] < -10) printf("%f\n", array[index]);
                    new[index2] += array[index]*dxdydz;
                    if(new[index2] > max) max=new[index2];
                    if(new[index2] < min) min=new[index2];
                    //printf("%f\n", new[index2]);
                }
            } 
        }    
        printf("Max: %f\tMin: %f\n", max, min);
	    int status=miset_real_value_hyperslab(img, MI_TYPE_FLOAT, wstarts, wcount, new);
    }
    /**/

    miset_volume_range ( img, max,  min);
    miclose_volume(img);
    free(new);
    free(array);
    return(0);
}


void useage(){
    printf("Useage: image.mnc step downsampled.mnc\n");
    printf("Purpose: Downsample by integrating over voxels.\n");
    exit(1);
}
