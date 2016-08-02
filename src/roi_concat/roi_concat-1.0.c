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
    if(strcmp(argv[1], "-help") == 0 ){
        useage();
    }


    //readVolume(data* volume, int read_hyperslab /*1 for read, o.w. don't read*/, mitype_t dtype)
    int nmasks=argc-2;
    printf("N Masks: %d\n", nmasks);
    data image, mask[nmasks];
    unsigned char* output_fn=argv[ argc-1]; //printf("%s\n", output_fn); exit(1);
    unsigned char** mask_array=malloc(sizeof(*mask_array) * nmasks);
    unsigned int* out;
    mihandle_t img;
    int old_zmax, old_ymax, old_xmax;


    for(int i=0; i<nmasks; i++){
        mask[i].filename = argv[i+1];
        mask_array[i] = (unsigned char*) readVolume( &(mask[i]) , 1, MI_TYPE_UBYTE);
        if(i != 0){
            if( old_zmax !=  mask[i].zmax || old_ymax !=  mask[i].ymax || old_xmax !=  mask[i].xmax ){ 
                    pexit("Dimensions do not match for", mask[i].filename, 1 );
            }
        }
        old_zmax = mask[i].zmax;
        old_ymax = mask[i].ymax;
        old_xmax = mask[i].xmax;
    }
//createVolume(char* newfilename, int num_dim, misize_t *sizes, double *separations, double *starts, mitype_t data_type)
    createVolume(output_fn, 3, mask[0].wcount, mask[0].step, mask[0].start, MI_TYPE_UINT);
    out=calloc(mask[0].n,sizeof(*out) );
    for(unsigned int z=0; z< mask[0].zmax; z++){
        for(unsigned int y=0; y< mask[0].ymax; y++){
            for(unsigned int x=0; x < mask[0].xmax; x++){
                int index=z*mask[0].ymax*mask[0].xmax+y*mask[0].xmax +x;
                for(int i=0; i<nmasks; i++){
                    if(mask_array[i][index] != 0){ 
                        out[index]= i+1;
                        //printf("%d\n", out[index]);
                    }
                }
            }
        } 
    }    

    miopen_volume(output_fn, MI2_OPEN_RDWR, &img);

    miset_volume_range ( img, nmasks,  0);
    int status=miset_real_value_hyperslab(img, MI_TYPE_UINT, mask[0].wstarts, mask[0].wcount, out);
    miclose_volume(img);
    return(0);
}


void useage(){
    printf("Useage: mask1.mnc mask2.mnc ... maskn.mnc out.mnc\n");
    printf("Purpose: Concatenate a series of masks and give them unique labels.\n");
    exit(1);
}
