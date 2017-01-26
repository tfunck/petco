#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <omp.h>
#include "minc2.h"
#include "hdf5.h"
#include "minc_helper.h"

void useage();



int check_dim(float wi, float cur_big_i, float old_step, float new_step, float l[2], int idx[2] ){
    //Check if the current world coordinate of the small voxel (wi) is within the current big voxel
    idx[0]= (int) floor(cur_big_i/new_step);
    if( wi > cur_big_i +new_step ){
        //Current small voxel is past the current big voxel
        /*
    cur_big_i
         0      1 wi    2
         |   .  | .   . |
         |   .  | .   . |
         |___.__|_.___._|
              \/ \/
              l0 l1    
        */
        l[1]=wi-(cur_big_i+new_step); //segment of small voxel that is past the edge of the big voxel
        //if( l[1] > old_step ){ printf("%f %f %f\n", wi, cur_big_i, new_step); } 
        l[0]=old_step - l[1]; //segment of the small voxel that is before the edge of the big voxel
        /*if( l[0] < 0 || l[1] < 0){
            printf("cur_big_i: %f\n", cur_big_i);
            printf("wi: %f\n", wi);
            printf("old step: %1.8f\n", old_step);
            printf("new step: %1.8f\n", new_step);
            printf("l[0] = %f\n", l[0]);
            printf("l[1] = %f\n", l[1]);
            printf("%f + %f = %f\n", l[0], l[1], l[0]+l[1] );
        }*/
        idx[1]= (int) floor( (cur_big_i+new_step) / new_step); //Update index 
        return(2);
    } else{
        l[0]=old_step; //
        return(1);
    }
}

int main(int argc, char** argv){
    if(argc != 4 || strcmp(argv[1], "-help") == 0 ){
        useage();
    }

    data image, mask;

    //readVolume(data* volume, int read_hyperslab /*1 for read, o.w. don't read*/, mitype_t dtype)
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
    new_step[zi]=new_step_size; 
    new_step[yi]=new_step_size; 
    new_step[xi]=new_step_size;
    wcount[zi]=(misize_t) ceil((wmax_z - image.zstart)/new_step[zi]); 
    wcount[yi]=(misize_t) ceil((wmax_y - image.ystart)/new_step[yi]); 
    wcount[xi]=(misize_t) ceil((wmax_x - image.xstart)/new_step[xi]);

    float new_vox_size=new_step[zi] * new_step[yi] * new_step[xi];
    float old_vox_size=image.zstep * image.ystep * image.xstep;
    new_starts[zi]=image.zstart;
    new_starts[yi]=image.ystart;
    new_starts[xi]=image.xstart;
    int z_big_idx[2], y_big_idx[2], x_big_idx[2];
    int nz, ny, nx;
    float lx[2], ly[2], lz[2];


    int n_big_vox=wcount[zi]*wcount[yi]*wcount[xi];
    new=calloc( n_big_vox, sizeof(*new));
    createVolume(output_filename, image.ndim, wcount, new_step, new_starts, MI_TYPE_FLOAT);
    miopen_volume(output_filename, MI2_OPEN_RDWR, &img);
    if(image.ndim==4) wcount[0]=1;
    min=max=0; //new[0];

    for( int t=0; t<tmax; t++){
        wstarts[0]=t;
        read_frame(&image, t, array, MI_TYPE_FLOAT);
        for(unsigned int i=0; i< wcount[zi] * wcount[yi] * wcount[xi]; i++) new[i]=0;
        float cur_big_z=0;//new_step[zi];
        for(unsigned int z=0; z< (int) zmax; z++){
            float cur_big_y=0;//new_step[yi];
            float wz=z*image.zstep;
            if( (nz=check_dim(wz, cur_big_z, image.zstep, new_step[zi], lz, z_big_idx )) == 2) cur_big_z += new_step[zi];
            for(unsigned int y=0; y< (int) ymax; y++){
                
                float wy=y*image.ystep;
                if( (ny=check_dim(wy, cur_big_y, image.ystep, new_step[yi], ly, y_big_idx )) == 2) cur_big_y += new_step[yi];
                float cur_big_x=0;//new_step[xi];
                for(unsigned int x=0; x < (int) xmax; x++){
                    float wx=x*image.xstep;
                    int small_vox_idx=z*xmax*ymax+y*xmax+x;
                    if( (nx=check_dim(wx, cur_big_x, image.xstep, new_step[xi], lx, x_big_idx )) == 2) cur_big_x += new_step[xi];
                    //
                    for(int k=0; k < nz; k++){
                        int bz=z_big_idx[k];
                        for( int j=0; j<ny; j++){
                            int by=y_big_idx[j];
                            for(int i=0; i<nx; i++){
                                int bx=x_big_idx[i];
                                int big_vox_idx=bz*(int) wcount[yi]*(int) wcount[xi]+by*(int)wcount[xi]+bx;
                                //factor is the volume of the small (old) voxels.
                                float factor=(lx[i]*ly[j]*lz[k]) / old_vox_size;
                                //if (factor < 0) {printf("lx: %f\tly: %f\tlz %f\t%f\n", factor); exit(1);}
                                //if(big_vox_idx >   wcount[zi]* wcount[yi]* wcount[xi]) printf("z:%d,%d  y:%d,%d  x:%d,%d\n", bz, wcount[zi],by, wcount[yi], bx, wcount[xi] );
                                /*if(x == xmax -1 && y==ymax-1){ 
                                    printf("Okay! %f\t %f %f %f\n", factor, lx[i], ly[j], lz[k] );
                                    printf("small vox idx %d %d %d, %d\n",z,y,x, small_vox_idx);
                                    printf("small world coord %f %f %f\n", wz, wy, wx);
                                    printf("big vox idx %d %d %d, %d\n", bz, by, bx, big_vox_idx);
                                    printf("big world coord ");
                                    if(nz==2) printf("%f ", cur_big_z-new_step[zi]);
                                    else printf("%f ", cur_big_z);
                                    if(ny==2) printf("%f ", cur_big_y-new_step[yi]);
                                    else printf("%f ", cur_big_y);
                                    if(nx==2) printf("%f ", cur_big_x-new_step[xi]);
                                    else printf("%f ", cur_big_x);
                                    printf("\n");
                                    printf("%d %d %d\n\n", i, j, k);
                                }*/
                                if(big_vox_idx < n_big_vox){
                                    new[big_vox_idx] +=  array[small_vox_idx] *  factor  ; 
                                } 
                                /*if( big_vox_idx == 14759  ) { // if(x == wcount[xi] -1 && y == wcount[yi] -1 ){ 
                                    printf("Okay! %f\t %f %f %f\n", factor, lx[i], ly[j], lz[k] );
                                    printf("small vox idx %d %d %d, %d\n",z,y,x, small_vox_idx);
                                    printf("small world coord %f %f %f\n", wz, wy, wx);
                                    printf("big vox idx %d %d %d, %d\n", bz, by, bx, big_vox_idx);
                                    printf("big world coord ");
                                    if(nz==2) printf("%f ", cur_big_z-new_step[zi]);
                                    else printf("%f ", cur_big_z);
                                    if(ny==2) printf("%f ", cur_big_y-new_step[yi]);
                                    else printf("%f ", cur_big_y);
                                    if(nx==2) printf("%f ", cur_big_x-new_step[xi]);
                                    else printf("%f ", cur_big_x);
                                    printf("\n");
                                    printf("new = %f\n\n", new[big_vox_idx] );
                                    exit(1);
                                }/**/
                           }
                        }
                    }
                }
            }
        }
        max=min=new[0];
        for(int i=0; i<n_big_vox ; i++){
            if(new[i] > max) max=new[i];
            if(new[i] < min) min=new[i];
        }
        //printf("Max: %f\tMin: %f\n", max, min);
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
