#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <omp.h>
#include "minc2.h"
#include "hdf5.h"
#include "minc_helper.h"

void useage();

int min(int a, int b){
    if( a < b){ 
        return a;
    }
    else return b;
}

int max(int a, int b){
    if( a > b) return a;
    else return b;
}

int main(int argc, char** argv){
    if(argc != 4 || strcmp(argv[1], "-help") == 0 ){
        useage();
    }

    int nimages=0;
    data image, mask;
    misize_t wstart[3];
    int maxima[3];
    //readVolume(data* volume, int read_hyperslab /*1 for read, o.w. don't read*/, mitype_t dtype)
    int* mask_vol;
    float* image_vol;
    float* new;
    image.filename=argv[1];
    mask.filename=argv[2];
    char* output_filename=argv[3];
    int xmin, xmax=0;
    int ymin, ymax=0; 
    int zmin, zmax=0;
    mihandle_t img;
    misize_t wstarts[3], wcount[3];
    wstarts[1]=wstarts[2]=0;

    mask_vol=(int*)readVolume(&mask, 1, MI_TYPE_INT);
    xmin=mask.xmax;
    ymin=mask.ymax;
    zmin=mask.zmax;
    
    //#pragma omp parallel  
    printf("z,y,x: %d %d %d\n", mask.zmax, mask.ymax, mask.xmax);
    for(unsigned long z=0; z<mask.zmax; z++){
        for(unsigned long y=0; y< mask.ymax; y++){
            for(unsigned long x=0; x < mask.xmax; x++){
                unsigned long index= z*(unsigned long) mask.ymax*mask.xmax+y*(unsigned long)mask.xmax + x;

                if( mask_vol[index] > 0){
                    zmin=min(z, zmin);
                    zmax=max(z, zmax);
                    ymin=min(y, ymin);
                    ymax=max(y, ymax);
                    xmin=min(x, xmin);
                    xmax=max(x, xmax);
                } 
            }
        }
    }
    float wz_min= zmin*mask.step[0]+mask.start[0];
    float wz_max= zmax*mask.step[0]+mask.start[0];
    float wy_min= ymin*mask.step[1]+mask.start[1];
    float wy_max= ymax*mask.step[1]+mask.start[1];
    float wx_min= xmin*mask.step[2]+mask.start[2];
    float wx_max= xmax*mask.step[2]+mask.start[2];

    free(mask_vol);
    readVolume(&image, 0, MI_TYPE_DOUBLE);

    image.start[0] = wz_min;
    image.start[1] = wy_min;
    image.start[2] = wx_min;
    wstart[0]=zmin;
    wstart[1]=ymin;
    wstart[2]=xmin;   

    new = malloc(sizeof(float) *image.count[0]*image.wcount[1]*image.wcount[2]);
    image_vol =  malloc(sizeof(float) *image.wcount[1]*image.wcount[2]);
    image.wcount[0] = (zmax - zmin) ;
    image.wcount[1] = (ymax - ymin) ;
    image.wcount[2] = (xmax - xmin) ;
    wcount[0]=1; 
    wcount[1]=mask.ymax;
    wcount[2]=mask.xmax;
    printf("Cropped Dimensions:\n");
    printf("\tVoxel Min: %d %d %d\n", zmin, ymin,xmin);
    printf("\tVoxel Max: %d %d %d\n", zmax, ymax, xmax);
    printf("\tMinima: %3.3f %3.3f %3.3f\n", wz_min, wy_min, wx_min);
    printf("\tMaxima: %3.3f %3.3f %3.3f\n", wz_max, wy_max, wx_max);
    printf("\tSizes: %d %d %d\n", image.wcount[0], image.wcount[1], image.wcount[2]);
    miopen_volume(image.filename, MI2_OPEN_READ, &img);
    printf("Opened MINC file: %s\n", image.filename);
    for(int z=wstart[0]; z< wstart[0]+image.wcount[0]; z++){
        wstarts[0]=z; 
        if(miget_real_value_hyperslab(img, MI_TYPE_FLOAT, wstarts, wcount, image_vol) != MI_NOERROR) pexit(image.filename, "Could not read values.", 1);
        //#pragma omp parallel  
        for(int y=wstart[1]; y<wstart[1]+image.wcount[1]; y++){
            for(int x=wstart[2]; x < wstart[2]+image.wcount[2]; x++){
                //unsigned long index=/*z*mask.ymax*mask.xmax+*/y*mask.xmax + x;
                unsigned long index=/*z*mask.ymax*mask.xmax+*/y*mask.xmax + x;
                unsigned long index2=(z-wstart[0])*image.wcount[1]*image.wcount[2]+(y-wstart[1])*image.wcount[2] + x-wstart[2];
                new[index2]=image_vol[index];
            }
        }
    }
    miclose_volume(img);
    free(image_vol);
    writeVolume(output_filename, new, image.start, image.step, image.wcount, MI_TYPE_FLOAT  );
    free(new);
    return(0);
}


void useage(){
    printf("Useage: image.mnc mask.mnc cropped_image.mnc\n");
    printf("Purpose: Crop image according to a binary mask\n");
    exit(1);
}
