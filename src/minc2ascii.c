#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include "minc2.h"
#include "hdf5.h"
#include "minc_helper.h"

void useage();

int main(int argc, char** argv){
    if(argc <= 3 || strcmp(argv[1], "-help") == 0 ){
        useage();
    }
    mihandle_t img;
    int nimages=0;
    data image;
    char* dtype_str = argv[1];
    int dtype;
    int mi_dtype;
    void* array;
    //readVolume(data* volume, int read_hyperslab /*1 for read, o.w. don't read*/, mitype_t dtype)
    image.filename=argv[2];
    FILE* outputfile=fopen(argv[3], "wt+");
    image.data=readVolume(&image, 1, MI_TYPE_DOUBLE);

    
    /*if( strcmp(dtype_str, "-int")){ 
        dtype = sizeof(int);
        mi_dtype = MI_TYPE_INT;
    } else if (strcmp(dtype_str, "-float")){
        dtype = sizeof(float);
        mi_dtype = MI_TYPE_FLOAT;
    } else if (strcmp(dtype_str, "-double")){
        dtype = sizeof(double);
        mi_dtype = MI_TYPE_DOUBLE;
    }*/
    fprintf(outputfile, "%d %d %d\n", image.count[0], image.count[1], image.count[2]);
    fprintf(outputfile, "%f %f %f\n", image.step[0], image.step[1], image.step[2]);
    for(int z=0; z<image.zmax; z++){
        for(int y=0; y< image.ymax; y++){
            for(int x=0; x < image.xmax; x++){
                int index=z*image.ymax*image.xmax+y*image.xmax+x;
                    if( strcmp(dtype_str, "-int")==0) fprintf(outputfile, "%d ",(int) image.data[index]);
                    else if (strcmp(dtype_str, "-float")==0)   fprintf(outputfile, "%f ", image.data[index]);
                    else if (strcmp(dtype_str, "-double")==0)  fprintf(outputfile, "%f ", image.data[index]);
                
            }
            fprintf(outputfile, "\n");
        }
        //fprintf(outputfile, "\n");
    }



    return(0);
}


void useage(){
    printf("Useage: minc2ascii <data type> input.mnc left.mnc right.mnc\n");
    printf("Purpose: Convert 3D minc to ascii\n");
    exit(1);
}
