#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include "minc2.h"
#include "hdf5.h"
#include "minc_helper.h"

void useage();

int main(int argc, char** argv){
    if(argc != 3 || strcmp(argv[1], "-help") == 0 ){
        useage();
    }
    mihandle_t img;
    int nimages=0;
    data image;
    image.filename=argv[1];
    char* output_filename=argv[2];
    float* out;
    image.data=readVolume(&image, 1, MI_TYPE_DOUBLE);

    out=calloc(image.n3d, sizeof(*out));
    for(int t=0; t<image.tmax; t++){
        for(int z=0; z<image.zmax; z++){
            for(int y=0; y< image.ymax; y++){
                for(int x=0; x < image.xmax; x++){
                    int index3d=z*image.ymax*image.xmax+y*image.xmax+x;
                    int index4d=t*image.zmax*image.ymax*image.xmax + z*image.ymax*image.xmax + y*image.xmax + x;
                    out[index3d] += image.data[index4d];
                }
            }
        }
    }


    writeVolume(output_filename, out, &(image.start[1]), &image.step[1], &image.wcount[1], MI_TYPE_FLOAT  );

    return(0);
}


void useage(){
    printf("Useage: mincsum input_4d.mnc output_3d.mnc\n");
    printf("Purpose: sum over time frames\n");
    exit(1);
}
