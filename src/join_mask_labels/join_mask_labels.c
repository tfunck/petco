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


int get_start(int x, int y){
    int diff=x-y;

    if( diff < 0 ) return(0);
    
    return(diff);

}

int fz(int i /*y*/, int j /*x*/, int k /*z*/, int* dim){ return(k*dim[1]*dim[0] + i*dim[0] + j);}
int fy(int i /*z*/, int j /*x*/, int k /*y*/, int* dim){ return(i*dim[2]*dim[1] + k*dim[1] + j);}
int fx(int i /*z*/, int j /*y*/, int k /*x*/, int* dim){ return(i*dim[1]*dim[2] + j*dim[2] + k);}

void write_z(float i /*y*/, float j /*x*/, float k /*z*/, FILE* fn){ fprintf(fn,"%f,%f,%f", k, i, j  ); }
void write_y(float i /*z*/, float j /*x*/, float k /*y*/, FILE* fn){ fprintf(fn,"%f,%f,%f", i, k, j  ); }
void write_x(float i /*z*/, float j /*y*/, float k /*x*/, FILE* fn){ fprintf(fn,"%f,%f,%f", i, j, k  ); }

int find_min(int* slab, int kernel_offset, int voxel_id,int dim[3], int cvox[3], int k, float* step, float* start, int (*f)(), float*** knots, int vox, int* nknots){
    int *cur_i=&(cvox[0]), *cur_j=&(cvox[1]),*cur_k=&(cvox[2]);
    int i_start=get_start(*cur_i, kernel_offset);
    int j_start=get_start(*cur_j, kernel_offset);
    int i_end=i_start+kernel_offset+1;
    int j_end=j_start+kernel_offset+1;
    int dim1=dim[0], dim2=dim[1], dim3=dim[2];
    int new_min=0;
    float min=dim[0]*step[0], d;
    
    for(int i=i_start; i < i_end; i++){
        for(int j=j_start; j<j_end; j++){
            int index=(*f)(i,j,k,dim);
            int value=slab[index];
            if( value == voxel_id  ){
                //Calculates distance of test voxel from projection of current voxel onto current plane
                d=step[0] * step[0] * (*cur_i-i) * (*cur_i-i) + step[1] * step[1] * (*cur_j-j) * (*cur_j-j) ; 
                if( d < min){ 
                    min=d; 
                }
                
                if( d <= min ){
                    d = min;
                    if(i_start == 0) *cur_i=i;
                    else *cur_i= i ;
                    if(j_start == 0) *cur_j=j;
                    else *cur_j= j ;
                    *cur_k=k;
                    new_min=1;
                }
            }
        }
    }
    if(new_min==1){
        //printf("Min: %f\t%d %d %d\n", min, *cur_i, *cur_j, *cur_k);
        *nknots += 1;
        knots[vox]=realloc(knots[vox], sizeof(**knots) * *nknots );
        knots[vox][*nknots-1]=NULL;
        knots[vox][*nknots-1]=realloc(knots[vox][*nknots-1], sizeof(***knots) * 3 );
        knots[vox][*nknots-1][0]=*cur_i * step[0]+start[0];
        knots[vox][*nknots-1][1]=*cur_j * step[1]+start[1];
        knots[vox][*nknots-1][2]=*cur_k * step[2]+start[2];
    }
    return(0);
}

int search_slab(int* slab,int kernel_offset, int cvox[3], int k, int dim[3], int voxel_id, float* step, float* start, int (*f)(int, int, int, int* ), float*** knots, int vox, int* nknots ){
    float iw, jw, kw;
    if(k==dim[2] ) return(0);
    find_min(slab, kernel_offset, voxel_id, dim, cvox, k, step, start, f, knots, vox, nknots);
    search_slab(slab,  kernel_offset,  cvox, k+1, dim, voxel_id, step, start, f, knots,vox, nknots);

    return(0);
}

int write_knots(FILE* knots_file, float*** knots,int n, int voxel_id, int nknots, void* (*write)(float, float, float, FILE*) ){
    fprintf(knots_file, "%d:",voxel_id);
    for(int i=0; i<nknots; i++){
        //fprintf(knots_file, "%3.3f,%3.3f,%3.3f", knots[n-1][i][0], knots[n-1][i][1], knots[n-1][i][2]);
        write(knots[n-1][i][0], knots[n-1][i][1], knots[n-1][i][2], knots_file);
        if(i<nknots-1) fprintf(knots_file, ";");
    }
    fprintf(knots_file, "\n");

    return(0);
}

int set_parameters(int* dim, int amax, int bmax, int cmax, float*start, float astart, float bstart, float cstart, float* step, float astep, float bstep, float cstep){
    dim[0]=amax;
    dim[1]=bmax;
    dim[2]=cmax;

    start[0]=astart;
    start[1]=bstart;
    start[2]=cstart;

    step[0]=astep;
    step[1]=bstep;
    step[2]=cstep;
    return(0);
}

int main(int argc, char** argv){
    if(argc != 7 || strcmp(argv[1], "-help") == 0 ){
        useage();
    }

    data mask;
    int* mask_vol;
    mask.filename=argv[1];
    char* search_dim=argv[2];
    int k=atoi(argv[3]);
    int search_radius=atoi(argv[4]);
    char* knot_fn=argv[5];
    char* output_fn=argv[6];
    FILE* knots_file=fopen(knot_fn, "w+");
    int kernel_offset;
    int i,j;
    int n=0;
    int max;
    int dim[3];
    int cvox[3];
    float start[3];
    float step[3];
    int* voxel_ids=NULL;
    float*** knots=NULL;
    int (*f)(int,int,int,int*);
    void (*write)(float, float, float, FILE*);
    int nknots=0;


    mask_vol=(int*)readVolume(&mask, 1, MI_TYPE_INT);
    //Figure out dimension along which we will use the splines
    if( strcmp(search_dim, "z") == 0){
        set_parameters(dim, mask.ymax, mask.xmax, mask.zmax, start, mask.ystart, mask.xstart, mask.zstart, step, mask.ystep, mask.xstep, mask.zstep);
        f=&fz;
        write=&write_z;
    }
    else if ( strcmp(search_dim, "y") == 0){ 
        set_parameters(dim, mask.zmax, mask.xmax, mask.ymax, start, mask.zstart, mask.xstart, mask.ystart, step, mask.zstep, mask.xstep, mask.ystep);
        f=&fy;
        write=&write_y;
    }
    else if ( strcmp(search_dim, "x") == 0){
        set_parameters(dim, mask.zmax, mask.ymax, mask.xmax, start, mask.zstart, mask.ystart, mask.xstart, step, mask.zstep, mask.ystep, mask.xstep);
        f=&fx;
        write=&write_x;
    }
    else { 
        printf("Error could not recognize dimension: %s. Please run with <z/y/x>.\n", search_dim);
        exit(0);
    }

    kernel_offset=round((search_radius + 1) / 2);

    int max_knots=100; 
    for(i=0; i < dim[0]; i++){
        for(j=0; j < dim[1]; j++){
            int index = (*f)(i,j,k,dim);
            int voxel_id=mask_vol[index];
            if( voxel_id == 12 /* > 0.1*/  ){
                n++;

                voxel_ids=realloc( voxel_ids, sizeof(*voxel_ids)*n );
                voxel_ids[n-1]=voxel_id;
                
                
                nknots=1;
                knots=realloc(knots, sizeof(*knots) * n  );
                knots[n-1]=malloc(sizeof(**knots) * nknots );
                knots[n-1][0]=malloc(sizeof(***knots) * 3);
                knots[n-1][0][0]=i*step[0]+start[0];
                knots[n-1][0][1]=j*step[1]+start[1];
                knots[n-1][0][2]=k*step[2]+start[2];

                cvox[0]=i;
                cvox[1]=j;
                cvox[2]=k;
                search_slab( mask_vol,  kernel_offset, cvox, k+1, dim, voxel_id, step, start, f, knots, n-1, &nknots);
                write_knots(knots_file, knots,n, voxel_id, nknots, write);
                printf("n: %d\tnknots: %d\n",n, nknots );
                //exit(0);
                if( n == max_knots) exit(0);
            }
        }
    }
    printf("Done.\n"); 

    return(0);
}


void useage(){
    printf("Useage: mask.mnc <axis: z, y, x> <axis start value> <voxel search radius> knots.txt output.mnc  \n");
    printf("Purpose: Join mask labels across slices. The axis start value refers to the \n");
    printf("Example: join_mask_labels mymask.mnc z 0 10 mymask_interpol.mnc\n");
    exit(1);
}
