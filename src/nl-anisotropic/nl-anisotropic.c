#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "math.h"
#include "minc2.h"
#include "hdf5.h"
#include "volume_io.h"
#include "minc_helper.h"
#include "nl-anisotropic.h"
#include "dsyevh3.h" //"3x3-C/dsyevh3.h"
#include "dsyevd3.h" //"3x3-C/dsyevd3.h"
#include "pthread.h"

void useage(){
    printf("\nName: nl-anisotropic\n");
    printf("Useage: <-verbose> fwhm dt lambda  mask.mnc image.mnc  output.mnc\n");
    printf("Purpose: Performs smoothing based on heat-diffusion equation\n");
    printf("\tImage_i+1 =  Image_i + (dt/(4*dx^2)) * Div( u(Grad(Mask) * Grad( Image_i ))\n");
    printf("\tThe mask image is used to define the structure tensor along which which determines\n ");
    printf("\tthe direction of diffusion. The parameter \"lambda\" is used to control the ammount\n");
    printf("\tof diffusion along the maximum gradient.\n");
    printf("\tlambda << 1 = diffusion across gradient goes to 0, i.e. more anisotropic\n");
    printf("\tlambda >> 1 = diffusion across gradient goes to 1, i.e. more isotropic\n\n");
    exit(1);
}

int main(int argc, char** argv){
    int argi=1;
    int expected_args=7;
    int VERBOSE=FALSE;
    if(argc>1){
        if(strcmp(argv[argi], "-verbose") == 0){ 
            VERBOSE=TRUE;
            argi++;
            expected_args++;
        }
    }
    if(argc != expected_args){ 
        useage();
    }
    
    data mask, image;
    double fwhm=atof(argv[argi++]);
    double dt=atof(argv[argi++]);
    double totalTime=fwhm*fwhm / (16*log(2));
    int maxIterations= (int) round( totalTime/ dt);
    dt=totalTime/maxIterations;
    double lambda=atof(argv[argi++]);
    mask.filename=argv[argi++];
    image.filename=argv[argi++];
    char* outputfilename=argv[argi++];
    int max_cpu=1; //sysconf(_SC_NPROCESSORS_ONLN);
    float *mask_data, *image_data;

    printf("dt = %f\n", dt);
    if(dt > (double) 1/7) {printf("Error: dt=%f is greater than 1/7 and hence will lead to unstable diffusion."); exit(1);}
    printf("Maximum number of iterations %d\n", maxIterations);

    mask_data =(float*) readVolume(&mask, 1, MI_TYPE_FLOAT);
    image_data=(float*) readVolume(&image, 1, MI_TYPE_FLOAT);
    if (checkDimensions(&mask, &image) == TRUE ) pexit("Mask and image dimensions are not the same.", "", 1);

    if( smooth(mask_data, image_data, mask.count, fwhm, dt,  mask.step[0], lambda, max_cpu,  VERBOSE)==0){
        printf("Completed.\n");
    }
    else{
        printf("Failed\n");
        return(1);
    }
    writeVolume(outputfilename, image_data, image.start, image.step, image.wcount, MI_TYPE_DOUBLE );
    printf("Output written to: %s\n", outputfilename);
    return(0);
}
                                   