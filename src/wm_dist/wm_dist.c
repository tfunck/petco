#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include "minc2.h"
#include "hdf5.h"
#include "volume_io.h"
#include "minc_helper.h"
#include "pthread.h"
#include <unistd.h>
void useage();

struct wm_vol_args{
    data* img;
    int* img_vol;
    float* mat; 
    int WM;
    int** border;
    int thread;
    int nthreads;
    int n;
};

int check_for_wm(int* img_vol, int z, int y, int x, int zmax, int ymax, int xmax, int WM){
    for( int i=-1; i <= 1; i++ )
        for( int j=-1; j <= 1; j++ ) 
            for( int k=-1; k <= 1; k++ ){
                int zi=z+i;
                int yi=y+j;
                int xi=x+k; 
                if( img_vol[zi*ymax*xmax+yi*xmax+xi ] == WM) return(1);
            }
    return(0);
}

int** wm_gm_border(data* img, int GM, int WM, int* img_vol, float* mat, int* n){
    int zmax=img->count[0];
    int ymax=img->count[1];
    int xmax=img->count[2];
    int** border=NULL;
    *n=0;
    /*****************
    *   Find border  * 
    ******************/
    for(int z=0; z<zmax; z++ ){
        for(int y=0; y<ymax; y++ ){
            for(int x=0; x<xmax; x++ ){
                int index=z*ymax*xmax+y*xmax+x;
                int val=img_vol[index];
                if(val==GM){
                    if( check_for_wm(img_vol, z, y, x, zmax, ymax, xmax, WM) == 1){
                        *n += 1;
                        border=realloc(border, sizeof(*border) * *n);
                        border[*n-1]=malloc( sizeof(**border) * 3);
                        border[*n-1][0]=z;
                        border[*n-1][1]=y;
                        border[*n-1][2]=x;
                        //printf("%d: %d %d %d\n", *n-1, z, y, x );
                    }
                } 
            }
        }
    }

    return(border);
}
int isin_border(int z,int y, int x, int**  border,int n){
    for(int i=0; i<n; i++) 
        if( z==border[i][0] && y==border[i][1] && x==border[i][2]  ) 
            return(1);
    return(0);
}

int find_border(int* img_vol, int* inner, int z, int y, int x, int zmax, int ymax, int xmax, int*** border, int* nborder, float* distances, int WM, int *chunk_size){
    int index0=z*ymax*xmax+y*xmax+x;
    int xymax=ymax*xmax;
    int** tborder=*border;
    int z0, z1, y1, y0, x0,x1;

    if(z+1 >= zmax) z1=z; else z1=z+1;
    if(z-1 < 0 ) z0=0 ; else z0=z-1;
    if(y+1 >= ymax ) y1=z ; else y1=y+1;
    if(y-1 < 0 ) y0=0 ; else y0=y-1;
    if(x+1 >= xmax ) x1=z; else x1=x+1;
    if(x-1 < 0 ) x0=0; else x0=x-1;
    for(int i=z0; i <= z1; i++ )
        for( int j=y0; j <= y1; j++ ) 
            for( int k=x0; k <= x1; k++ ){
                    int index=i*xymax+j*xmax+k;
                    //IF DISTANCE == MIN POSSIBLE DISTANCE, SET VOXEL TO INSIDE?
                    if( img_vol[index] == WM){
                        if( inner[index] == 0 ){
                            float d=sqrt( (k-x)*(k-x) + (j-y)*(j-y) + (i-z)*(i-z) )+distances[index0];
                            if(*nborder >= *chunk_size){
                                *chunk_size += *chunk_size;
                                tborder=realloc(tborder, sizeof(*tborder) * *chunk_size);
                            }
                            tborder[*nborder]=malloc(sizeof(**tborder)*3);
                            tborder[*nborder][0]=i; 
                            tborder[*nborder][1]=j; 
                            tborder[*nborder][2]=k;
                            if( d < distances[index] || distances[index]==0 ){
                                 distances[index]=d;
                            }
                            *nborder += 1;
                        }
                    }
                }
             
    *border=tborder;
    return(0);
}


void wm_dist_threaded(void* args){
    data* img= ((struct wm_vol_args*) args)->img;
    int WM= ((struct wm_vol_args*) args)->WM;
    int* img_vol=((struct wm_vol_args*) args)->img_vol;
    float* mat=((struct wm_vol_args*) args)->mat;
    int** gm_border=((struct wm_vol_args*) args)->border;
    int n=((struct wm_vol_args*) args)->n;
    int thread= ((struct wm_vol_args*) args)->thread;
    int nthreads= ((struct wm_vol_args*) args)->nthreads;
    int thread_step=nthreads;
    //Calculate inner and outer border voxels
    //Calculate distances to outer border voxels
    //Add border voxels to inner region
    int zmax=img->zmax;
    int ymax=img->ymax;
    int xmax=img->xmax;
    int max =img->n3d;
    int* inner=calloc(max, sizeof(*inner));
    float* distances=calloc(max, sizeof(distances));
    int** border=NULL, **old_border=NULL;
    int nborder=0, old_nborder=0;
    int counter=0;
    char testfn[100];
    int sum=0, wm_total;
  
    for(int i =0;i<max; i++) if(img_vol[i]==WM) wm_total++;
    int chunk_size = wm_total/5, base_chunk_size;
    base_chunk_size=chunk_size;

    old_border=malloc(sizeof(*old_border) * chunk_size);


    for(int i=thread; i<n; i += thread_step){
        int z=gm_border[i][0];
        int y=gm_border[i][1];
        int x=gm_border[i][2];
        int index=z*ymax*xmax+y*xmax+x;
            find_border(img_vol, inner, z, y, x, zmax, ymax, xmax, &old_border, &old_nborder, distances, WM, &chunk_size);
        printf("WM,NBorder\n");
        counter=0;
        while (old_nborder > 0){
            nborder=0;
            //Allocate memory for a new WM border
            border=malloc(sizeof(*border) * chunk_size);

            for(int i=0; i< old_nborder; i++ ){
                if(counter % 10000000 ) printf("%3.3f\%\r", (float) 100. * i / old_nborder ); 
                int cz=old_border[i][0];
                int cy=old_border[i][1];
                int cx=old_border[i][2];
                index=cz*ymax*xmax+cy*xmax+cx;
                    find_border(img_vol, inner, cz, cy, cx, zmax, ymax, xmax, &border, &nborder, distances, WM, &chunk_size);
            }
           
            //set border voxels to inside
            for(int i=0; i<nborder; i++){
                int cz=border[i][0];
                int cy=border[i][1];
                int cx=border[i][2];
                int index=cz*ymax*xmax+cy*xmax+cx;
                //int index=/*z*/border[i][0]*ymax*xmax+/*y*/border[i][1]*xmax+/*x*/border[i][2];
                if(inner[index]==0) sum++; 
                inner[index]=1;
                
            }
            //printf("Completed: %3.3f\% =  %d / %d, old_nborder=%d\n", (float) 100* sum / wm_total, sum, wm_total, old_nborder);
            printf("%3.3f,%d\n", (float) 100* sum / wm_total, old_nborder);
            counter++;

            free(old_border);
            old_nborder=nborder;
            old_border=border;
            border=NULL;

            if(counter==10) break;
        }
        sprintf(testfn,"test_%d.mnc", counter );
        writeVolume(testfn, distances, img->start, img->step, img->wcount, MI_TYPE_FLOAT );


    }



}

int wm_dist_multithreaded(data* img, int* img_vol, int** border, int WM, float* mat, int n, int nthreads){
    int rc;
    pthread_t threads[nthreads];
    struct wm_vol_args thread_args[nthreads];
    for(int t=0; t< nthreads; t++){
        thread_args[t].img_vol=img_vol;
        thread_args[t].img=img;
        thread_args[t].WM=WM;
        thread_args[t].mat=mat;
        thread_args[t].n=n;
        thread_args[t].border=border;
        thread_args[t].thread=t;
        thread_args[t].nthreads=nthreads;
        rc= pthread_create(&threads[t], NULL, wm_dist_threaded, (void *) &thread_args[t] ); 
        if(rc!=0) pexit("Error: creating thread","" , 1);
    }
    for(int t=0; t<nthreads;t++){
        rc = pthread_join(threads[t], NULL);
        if(rc!=0) pexit("Error: joining thread", "", 1);
    }
    return(0);
}


int main(int argc, char** argv){
    int argi=1;
    int expected_args=5;
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
    
    data img;
    img.filename=argv[argi++];
    int GM=atoi(argv[argi++]);
    int WM=atoi(argv[argi++]);
    char* outputfilename=argv[argi++];
    int max_cpu=1; //sysconf(_SC_NPROCESSORS_ONLN);
    int *img_vol;
    int n=0;
    float* mat=NULL;

    img_vol=(int*) readVolume(&img, 1, MI_TYPE_INT);
    int** gm_border=wm_gm_border(&img, GM, WM, img_vol,  mat, &n);
    printf("Number of GM-WM border voxels: %d\n", n);
    wm_dist_multithreaded(&img, img_vol, gm_border, WM, mat, n, max_cpu);
    

    printf("Output written to: %s\n", outputfilename);
    return(0);
}
                            
void useage(){
    printf("Name: wm_dist\n");
    printf("Purpose: calculate distances between GM voxels via WM\n");
    printf("Useage: wm_dist <classified.mnc> <GM value> <WM value> <matrix.csv>\n");
    exit(1);
}

    /*int i0=z*xymax+y*xmax+x*x0; //oo0
    int i1=z*xymax+y0*xmax+x*x0; //o00
    int i2=z0*xymax+y0*xmax+x*x0; //000

    int i3=z0*xymax+y0*xmax+x*x; //00o
    int i4=z0*xymax+y*xmax+x*x; //0oo
    int i5=z*xymax+y0*xmax+x*x; //o0o

    int i6=z0*xymax+y*xmax+x*x0; //0o0
    int i7=z0*xymax+y*xmax+x*x; //0oo

    int i8=z*xymax+y*xmax+x*x1; //oo1
    int i9=z*xymax+y1*xmax+x*x1; //o11
    int i10=z1*xymax+y1*xmax+x*x1; //111

    int i11=z1*xymax+y1*xmax+x*x; //11o
    int i12=z1*xymax+y*xmax+x*x; //1oo
    int i13=z*xymax+y1*xmax+x*x; //o1o

    int i14=z1*xymax+y*xmax+x*x1; //1o1
    int i15=z1*xymax+y*xmax+x*x; //1oo

    int i16=z1*xymax+y0*xmax+x*x0; //100
    int i17=z1*xymax+y1*xmax+x*x0; //110
    int i18=z1*xymax+y1*xmax+x*x1; //111

    int i19=z1*xymax+y*xmax+x*x1; //101
    int i20=z0*xymax+y0*xmax+x*x1; //001
    int i21=z0*xymax+y1*xmax+x*x1; //011

    int i22=z1*xymax+y0*xmax+x*x; //10o
    int i23=z*xymax+y0*xmax+x*x1; //o01
    int i24=z*xymax+y1*xmax+x*x1; //o11

    int i22=z1*xymax+y*xmax+x*x; //1oo
    int i23=z0*xymax+y*xmax+x*x; //0oo
    int i24=z*xymax+y*xmax+x*x1; //oo1
    int i25=z*xymax+y*xmax+x*x0; //oo0*/
