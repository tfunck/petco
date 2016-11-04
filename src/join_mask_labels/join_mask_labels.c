#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <omp.h>
#include "minc2.h"
#include "hdf5.h"
#include "minc_helper.h"

const int I=369; 
const int J=604; //1312; 
const int K=11; 
const int Kp=12;//12; 

void useage();
int set_parameters(int* dim, int amax, int bmax, int cmax, float*start, float astart, float bstart, float cstart, float* step, float astep, float bstep, float cstep);
int get_lin_domain(double i0, double j0, double k0, double i1, double j1, double k1, double t, double* cur_i, double* cur_j, double* cur_k);
int set_dim_parameters(data* mask, char* search_dim, int* dim, float* start, float* step  );
float inverse3x3(double* m);
int min(int a, int b);
int max(int a, int b);
int get_start(int x, int y);
int fz(int i /*y*/, int j /*x*/, int k /*z*/, int* dim){ return(k*dim[1]*dim[0] + i*dim[0] + j);}
int fy(int i /*z*/, int j /*x*/, int k /*y*/, int* dim){ return(i*dim[2]*dim[1] + k*dim[1] + j);}
int fx(int i /*z*/, int j /*y*/, int k /*x*/, int* dim){ return(i*dim[1]*dim[2] + j*dim[2] + k);}

void write_z(float i /*y*/, float j /*x*/, float k /*z*/, FILE* fn){ fprintf(fn,"%f,%f,%f", k, i, j  ); }
void write_y(float i /*z*/, float j /*x*/, float k /*y*/, FILE* fn){ fprintf(fn,"%f,%f,%f", i, k, j  ); }
void write_x(float i /*z*/, float j /*y*/, float k /*x*/, FILE* fn){ fprintf(fn,"%f,%f,%f", i, j, k  ); }


int search(int i, int j, int k, int* found_next_voxel, int voxel_id, int* slab, int* cvox, int* nvox, float* step, int* dim, double* min,  int (*f)()){
    int index=(*f)(i,j, k,dim);
    int value=slab[index];
    double d;
    
    if( value == voxel_id  ){
        d=(double) step[0] * step[0] * (cvox[0]-i) * (cvox[0]-i) + step[1] * step[1] * (cvox[1]-j) * (cvox[1]-j); 

        if( d < *min ){
            *min=d;
            nvox[0]=i;
            nvox[1]=j;
            nvox[2]=k;
            *found_next_voxel=1;
        }
    }
    return(0);
}

int find_min(int* slab, int kernel_offset, int voxel_id,int* dim, int* cvox, int* nvox, int k, float* step, float* start, int (*f)()){
    int *cur_i=&(cvox[0]), *cur_j=&(cvox[1]),*cur_k=&(cvox[2]);
    int index, index_min=0;
    double min=(double) dim[0]*step[0], d;
    int found_next_voxel=0;
    int counter=0;
    int i0, i1, j0, j1, i, j, i_start, i_end, j_start, j_end;

    //1 Find closest next voxel
    while( found_next_voxel == 0 && counter < kernel_offset){
        int i_start = cvox[0]-counter;
        int j_start = cvox[1]-counter;
        int i_end = cvox[0] + counter + 1;
        int j_end = cvox[1] + counter + 1;
        if(i_start < 0) i_start=0;
        if(j_start < 0) j_start=0;
        if(i_end > dim[0]) i_end=dim[0];
        if(j_end > dim[1]) j_end=dim[1];
        
        i0=i_start;
        i1=i_end-1;
        for(j=j_start; j<j_end; j++){ 
            search(i0, j, k, &found_next_voxel, voxel_id, slab, cvox, nvox, step, dim, &min, f);
            search(i1, j, k, &found_next_voxel, voxel_id, slab, cvox, nvox, step, dim, &min, f);
        }

        j0=j_start;
        j1=j_start-1;
        for(i=i_start; i < i_end; i++){  
            search(i, j0, k, &found_next_voxel, voxel_id, slab, cvox, nvox, step, dim, &min, f);
            search(i, j1, k, &found_next_voxel, voxel_id, slab, cvox, nvox, step, dim, &min, f);
        }

        if(found_next_voxel == 1) break;
        counter++;
    }
    //printf("Counter: %d\n", counter);

    return(found_next_voxel);
}




int fill_voxels(int* mask_vol, int* cvox, int* nvox, float* step, float* start, int* dim, float* distances,int voxel_id,int (*findex)(int, int, int, int*)){
    //If closest next voxel is found, fill in voxel between if the line between is closest
    double t_step=1/(ceil(fabs(cvox[2]-nvox[2])));  
    double i0=voxel2real(cvox[0],start[0], step[0]); 
    double j0=voxel2real(cvox[1],start[1], step[1]);
    double k0=voxel2real(cvox[2],start[2], step[2]);
    double i1=voxel2real(nvox[0],start[0],step[0]);
    double j1=voxel2real(nvox[1],start[1],step[1]);
    double k1=voxel2real(nvox[2],start[2],step[2]);
    double i, j, k;
    double d;

    for(double t=0; t < 1; t += t_step ){
        get_lin_domain(i0, j0, k0, i1, j1, k1, t, &i, &j, &k); 

        int i_vox=real2voxel((float) i, start[0], step[0]);
        int j_vox=real2voxel((float) j, start[1], step[1]);
        int k_vox=real2voxel((float) k, start[2], step[2]);

        double i_vox_world= ((double) i_vox)*step[0]+start[0];
        double j_vox_world= ((double) j_vox)*step[1]+start[1];
        double k_vox_world= ((double) k_vox)*step[2]+start[2];
        d= ( i_vox_world - i )*( i_vox_world - i );
        d+=( j_vox_world - j )*( j_vox_world - j );
        d+=( k_vox_world - k )*( k_vox_world - k );

        int index=(*findex)(i_vox,j_vox,k_vox,dim);
        if( d < distances[index] && mask_vol[index] == 0 ){
            distances[index]= d;
            mask_vol[index]= -1 * voxel_id; 
        }
    }
    //Update current voxel

    return(0);
}

int search_slab(int* mask_vol,int kernel_offset, int* cvox, int k, int* dim, int voxel_id, float* step, float* start, int (*f)(int, int, int, int* ), int search_depth, float* distances, int* last_k, int sign ){
    int max_slice=cvox[2]+search_depth;
    int found_next_voxel=0;
    int nvox[3];
    if( k==dim[2] || ( k >  max_slice && sign==1) || (k<0 && sign == -1)  ) return(0);
    found_next_voxel = find_min(mask_vol, kernel_offset, voxel_id, dim, cvox, nvox, k, step, start, f);
    if( found_next_voxel == 1) {
        fill_voxels(mask_vol, cvox, nvox, step, start,dim, distances, voxel_id, f);
        cvox[0]=nvox[0];
        cvox[1]=nvox[1];
        cvox[2]=nvox[2];
        if( k > *last_k) *last_k=k;
    }

    search_slab(mask_vol,kernel_offset, cvox, k+sign*search_depth, dim, voxel_id, step, start, f, search_depth, distances, last_k, sign);
    return(0);
}


int get_knots(int* dim, int* mask_vol, int k, int* last_k, int* n, int search_depth, float* step, float* start, int search_radius, int (*f)(int, int, int, int* ), float* distances, int sign ){
    int kernel_offset=round((search_radius + 1) / 2);
    int m;
    int cvox[3];

    for(int i=0; i < dim[0]; i++){
        for(int j=0; j < dim[1]; j++){
            int index = (*f)(i,j,k,dim);
            int voxel_id=mask_vol[index];
            if( voxel_id > 0 /*== label*/  ){
                *n += 1;
                cvox[0]=i;
                cvox[1]=j;
                cvox[2]=k;
                
                search_slab(mask_vol, kernel_offset, cvox, k+sign*search_depth, dim, voxel_id, step, start, f, search_depth, distances, last_k, sign);
            }
        }
    }


    return(0);
}


int iterate_over_knots(int* dim, int* mask_vol, int* n, int search_depth, float* step, float* start, int search_radius, int (*f)(int, int, int, int* ), float* distances,  int* last_k,  int k_start, int k_end, int sign){
    for(int k=k_start; k != k_end; k += sign * 1){
        //printf("\r%3.1f\%", sign, (double) 100* k/dim[2]);
        printf("%d\n", k);
        fflush(stdout);
        get_knots(dim, mask_vol, k, last_k,  n, search_depth, step, start, search_radius, f, distances,sign);
        if( *last_k > k && sign==1){ 
            k=*last_k;
            //printf("\t%d\n", *last_k);
        }
        else if( *last_k < k && sign==-1){
            k=*last_k;
        }
        //if( k > 50 ) break;
    }
    printf("\n");
    return(0);
}
int main(int argc, char** argv){
    if(argc != 5 || strcmp(argv[1], "-help") == 0 ){
        useage();
    }

    data mask;
    int* mask_vol;
    int argi=1;
    mask.filename=argv[argi++];
    int search_radius=atoi(argv[argi++]);
    int search_depth=atoi(argv[argi++]);
    char* output_fn=argv[argi++];
    char* search_dim;
    int kernel_offset;
    int i,j;
    int n;
    int max;
    int dim[3];
    int last_k;
    float start[3];
    float step[3];
    float* distances;
    char* dim_list[]={"y", "z", "x"};
    int (*f)(int,int,int,int*);
    void (*write)(float, float, float, FILE*);


    mask_vol=(int*)readVolume(&mask, 1, MI_TYPE_INT);
    distances=malloc( mask.n3d * sizeof(*distances));
    for(int i=0; i < mask.n3d; i++){
        distances[i]=(float) mask.n3d ;//set a big number for default distance
    }

    for(int d=0; d<1; d++){
        search_dim=dim_list[d];
        printf("Filling axis: %s\n", search_dim);
        if( strcmp(search_dim, "z") == 0) f=&fz;
        else if ( strcmp(search_dim, "y") == 0) f=&fy;
        else if ( strcmp(search_dim, "x") == 0) f=&fx;
        else useage();

        set_dim_parameters(&mask, search_dim, dim, start, step);
        printf("Positive dimension step\n");
        iterate_over_knots(dim, mask_vol, &n, search_depth, step, start, search_radius, f, distances,  &last_k,  0, dim[2],  1);
        printf("Negative dimension step\n");
        iterate_over_knots(dim, mask_vol, &n, search_depth, step, start, search_radius, f, distances,  &last_k, last_k, 0,  -1);
    }
    for(int i=0; i < mask.n3d; i++){
        if(mask_vol[i] < 0) mask_vol[i] *= -1;
    }

    writeVolume(output_fn, mask_vol, mask.start, mask.step, mask.wcount, MI_TYPE_INT  );
    printf("Done.\n"); 

    return(0);
}


void useage(){
    printf("Useage: mask.mnc <axis: z, y, x> <axis start value> <voxel  search radius> output.mnc  \n");
    printf("Purpose: Join mask labels across slices. The axis start value refers to the \n");
    printf("Example: join_mask_labels mymask.mnc z 0 10 mymask_interpol.mnc\n");
    exit(1);
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



int set_dim_parameters(data* mask, char* search_dim, int* dim, float* start, float* step  ){
    //Figure out dimension along which we will use the splines
    if( strcmp(search_dim, "z") == 0){
        set_parameters(dim, mask->ymax, mask->xmax, mask->zmax, start, mask->ystart, mask->xstart, mask->zstart, step, mask->ystep, mask->xstep, mask->zstep);
    }
    else if ( strcmp(search_dim, "y") == 0){
        set_parameters(dim, mask->zmax, mask->xmax, mask->ymax, start, mask->zstart, mask->xstart, mask->ystart, step, mask->zstep, mask->xstep, mask->ystep);
    }
    else if ( strcmp(search_dim, "x") == 0){
        set_parameters(dim, mask->zmax, mask->ymax, mask->xmax, start, mask->zstart, mask->ystart, mask->xstart, step, mask->zstep, mask->ystep, mask->xstep);
    }
    else { 
        printf("Error could not recognize dimension: %s. Please run with <z/y/x>.\n", search_dim);
        exit(0);
    }

    return(0);
}
float inverse3x3(double* m){
    //[ 0, 1, 2 ]
    //[ 3, 4, 5 ]
    //[ 6, 7, 8 ]
    //Get minors
    double d0 = m[4] * m[8] - m[7] * m[5];
    double d1 = m[3] * m[8] - m[6] * m[5];
    double d2 = m[3] * m[7] - m[6] * m[4];
    double d3 = m[1] * m[8] - m[7] * m[2];
    double d4 = m[0] * m[8] - m[6] * m[2];
    double d5 = m[0] * m[7] - m[6] * m[1];
    double d6 = m[1] * m[5] - m[4] * m[2];
    double d7 = m[0] * m[5] - m[3] * m[2];
    double d8 = m[0] * m[4] - m[3] * m[1];
    double det = m[0] * d0 - m[3] * d3 + m[6] * d6;
    //printf("%f %f\n", d0, m[0]);
    //printf("%f = %f - %f + %f\n", det, m[0] * d0, m[3] * d3, m[6] * d6);
    if(det != 0){
        //[ 0, -3, 6 ]
        //[ -1, 4, -7 ]
        //[ 2, -5, 8 ]
        m[0]= d0 / det;
        m[1]=-d3 / det;
        m[2]= d6 / det;
        m[3]=-d1 / det;
        m[4]= d4 / det;
        m[5]=-d7 / det;
        m[6]= d2 / det;
        m[7]=-d5 / det;
        m[8]= d8 / det;

    }
    return(det);
}
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

int get_lin_domain(double i0, double j0, double k0, double i1, double j1, double k1, double t, double* cur_i, double* cur_j, double* cur_k){
    //printf("y0: %f, y1: %f, t: %f\n", y0, y1, t);
    *cur_i = i0 + (i1 - i0) * t;
    *cur_j = j0 + (j1 - j0) * t;
    *cur_k = k0 + (k1 - k0) * t;
   return(0); 
}





