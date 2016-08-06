#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include "minc2.h"
//#include "hdf5.h"
#include "minc_helper.h"

void useage();

int main(int argc, char** argv){
    if(argc <= 3 || strcmp(argv[1], "-help") == 0 ){
        useage();
    }
    mihandle_t img;
    int nimages=0;
    data image, *out1, *out2, *left, *right;
    //readVolume(data* volume, int read_hyperslab /*1 for read, o.w. don't read*/, mitype_t dtype)
    image.filename=argv[1];
    char* output_filename1=argv[2];
    char* output_filename2=argv[3];

    image.data=readVolume(&image, 1, MI_TYPE_DOUBLE);
    out1=copyVolume(&image);
    out2=copyVolume(&image);
    out1->data=calloc(out1->n, sizeof(*out1->data));
    out2->data=calloc(out2->n, sizeof(*out2->data));
    if(image.xstep < 0 ) {
        left=out2;
        right=out1;
    } else {
        left=out1;
        right=out2;
    }  


    for(int z=0; z<image.zmax; z++){
        for(int y=0; y< image.ymax; y++){
            for(int x=0; x < image.xmax; x++){
                int index=z*image.ymax*image.xmax+y*image.xmax+x;
                if(x < image.xmax/2) {
                    left->data[index] = image.data[index];
                }
                else{ 
                    right->data[index]=image.data[index];
                }
            }
        }
    }

    left->filename = output_filename1;
    right->filename = output_filename2; 

    printf("%f %f %f\n", image.start[0], image.start[1], image.start[2]);
    printf("%f %f %f\n", image.step[0], image.step[1], image.step[2]);
    printf("%d %d %d\n", image.wcount[0], image.wcount[1], image.wcount[2]);

    printf("%f %f %f\n", left->start[0], left->start[1], left->start[2]);
    printf("%f %f %f\n", left->step[0], left->step[1], left->step[2]);
    printf("%d %d %d\n", left->wcount[0], left->wcount[1], left->wcount[2]);

    writeVolume(left->filename, left->data, left->start, left->step, left->wcount, MI_TYPE_DOUBLE  );
    writeVolume(right->filename, right->data, right->start, right->step, right->wcount, MI_TYPE_DOUBLE  );

    return(0);
}


void useage(){
    printf("Useage: mincsplit input.mnc left.mnc right.mnc\n");
    printf("Purpose: Split the image down the x axis\n");
    exit(1);
}
