#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <malloc.h>
#include <math.h>
#include <pthread.h>
#include <unistd.h>
#include "minc2.h"
#include "minc_helper.h"
#define TRUE 1
#define FALSE 0
#define MAX 99999 //FIXME: Only makes sense for brains at the mm scale, may have to be adapted for other applications
int VERBOSE=FALSE;


struct temp_node;
struct temp_node;
struct temp_edge;
typedef struct temp_node {
  double z;
  double y;
  double x;
  int i, j, k; //voxel indices
  double dist;
  double val;
  //for voxels the pointer below should 
  //point to the nearest node. For vertices
  //the pointer should go to the ngh 
  //of that node
  //the info below only applies to vertices
  //on a mesh, not to voxels
  struct temp_node** ngh;
  struct temp_edge** outgoing_edges;
  int n_outgoing_edges;
  int inside;
  int border;
  int nngh;
  int index;
} node;


void* wm_dist_threaded(void*);
double dist(double a, double b);
int check_if_border(node* c_vtx);
void useage();



struct wm_vol_args{
    data* img;
    node** mesh;
    node** gm_border;
    int* img_vol;
    int* valid_indices;
    int nvalid;
    float* mat; 
    int WM;
    int thread;
    int nthreads;
    int n;
};




typedef struct temp_edge {
    node* start;
    node* end;
    float length;
    int travelled;

} edge;


node** read_vertices_from_volume(data* anatomy,int* img_vol, int* nvertices, int WM, int GM){
    int zmax=anatomy->zmax;
    int ymax=anatomy->ymax;
    int xmax=anatomy->xmax;
    node** mesh=malloc(zmax*ymax*xmax*sizeof(node*));
    if(mesh ==NULL) pexit("Error allocating mesh array", "", 1);
    
    for(int z=0; z < zmax; z++ ){
        for(int y=0; y < ymax; y++ ){
            for(int x=0; x < xmax; x++ ){
                unsigned long i=z*ymax*xmax+xmax*y+x;
                if( img_vol[i] > 0 ){
                    mesh[i]=malloc(sizeof(**mesh));
                    mesh[i]->i=x;
                    mesh[i]->j=y;
                    mesh[i]->k=z;
                    mesh[i]->x=voxel2real(x,anatomy->xstart, anatomy->xstep); 
                    mesh[i]->y=voxel2real(y,anatomy->ystart, anatomy->ystep); 
                    mesh[i]->z=voxel2real(z,anatomy->zstart, anatomy->zstep); 
                    mesh[i]->index=i;
                    mesh[i]->val=img_vol[i];
                    mesh[i]->dist=0;
                    mesh[i]->n_outgoing_edges=0;
                    mesh[i]->outgoing_edges=NULL;
                    mesh[i]->border=FALSE;
                    mesh[i]->ngh=NULL;
                    mesh[i]->nngh=0;
                    //else 
                    mesh[i]->inside=FALSE;
                }
            else mesh[i]=NULL;
            }
        }
    }
    *nvertices=zmax*ymax*xmax;
    return(mesh);
}


node*** copy_node(node **mesh, int nmesh, int nnode){
    node*** list=malloc(sizeof(**list)*nmesh);//Allocate list of list of nodes
    list[0]=mesh; //Point first element of list to existing node
    for(int i=1; i<nmesh; i++){ //Starting from the second element, create new lists of nodes
        list[i]=malloc(sizeof(**list)*nnode); //Allocate memory for the new list of nodes
        for(int j=0; j<nnode; j++){
            if(mesh[j] != NULL){ 
                //for all the node points that exist(i.e., are not NULL), copy elements from existing node list (i.e., mesh) to a new node list. 
                list[i][j]=malloc(sizeof(***list));
                list[i][j]->z=mesh[j]->z;
                list[i][j]->y=mesh[j]->y;
                list[i][j]->x=mesh[j]->x;
                list[i][j]->i=mesh[j]->i;
                list[i][j]->j=mesh[j]->j;
                list[i][j]->k=mesh[j]->k; 
                list[i][j]->dist=mesh[j]->dist;
                list[i][j]->val=mesh[j]->val; 
                list[i][j]->ngh=mesh[j]->ngh; 
                list[i][j]->outgoing_edges=mesh[j]->outgoing_edges; 
                list[i][j]->n_outgoing_edges=mesh[j]->n_outgoing_edges; 
                list[i][j]->inside=mesh[j]->inside; 
                list[i][j]->border=mesh[j]->border; 
                list[i][j]->nngh=mesh[j]->nngh; 
                list[i][j]->index=mesh[j]->index; 
            }
        }

    }
   return(list); 

}

/*Name: Connect_nodes
 *Purpose: For a given node/node (c_vtx), look around the neighbouring voxels for those that are inside the brain (inside_brain==TRUE) 
 *          When a neigbouring voxel is found, add a neighbour to (ngh) and point to the neighbour. This way we can keep track of all
 *          the current node's neighbours via pointers stored in ngh.
 */

int connect_nodes(int zmax, int ymax, int xmax, node* c_vtx, node** mesh, int WM ){
    node* t_vtx;
    int x=c_vtx->i;
    int y=c_vtx->j;
    int z=c_vtx->k;
    for(int c=-1; c<=1; c++){
        for(int b=-1; b<=1; b++){
            for(int a=-1; a<=1; a++){
                if( (a==0 && b==0 && c==0) ==FALSE) {
                    int zi=z+c;
                    int yi=y+b;
                    int xi=x+a;
                    if( zi >= 0 && zi < zmax && yi >= 0 && yi < ymax && zi >= 0 && zi < zmax){
                        int index=zi*ymax*xmax+yi*xmax+xi;
                        t_vtx=mesh[index];
                        if(t_vtx != NULL){ 
                            if( t_vtx->val == WM){
                                //if( c_vtx->val == 2) printf("Neigbour: %d\n", t_vtx->i);
                                c_vtx->nngh += 1;
                                c_vtx->ngh=realloc( c_vtx->ngh, (c_vtx->nngh+1) * sizeof(node*));
                                c_vtx->ngh[ c_vtx->nngh]=t_vtx;

                            }
                        }
                    }
                }
            }
        }
    }
    return(0);
}
int create_edges(node* c_vtx, node** mesh){
    for(int k=1; k< c_vtx->nngh; k++ ){
            node* t_vtx=mesh[ c_vtx->ngh[k]->index ] ;
            //if we find a neighbouring node that is outside

            double length=sqrt( dist(c_vtx->x, t_vtx->x) + dist(c_vtx->y, t_vtx->y) + dist(c_vtx->z, t_vtx->z)  );
            edge* e=malloc(sizeof(edge));
            e->start=c_vtx;
            e->end=t_vtx;
            e->travelled=FALSE;
            e->length=(float) length;
            c_vtx->outgoing_edges=realloc(c_vtx->outgoing_edges, sizeof(edge*) * (1+c_vtx->n_outgoing_edges));
            c_vtx->outgoing_edges[c_vtx->n_outgoing_edges]=e;
            c_vtx->n_outgoing_edges += 1;

            

            //increment number of outgoing and incoming edges
        }
    if(c_vtx->n_outgoing_edges > 27 ) { printf("Problem\n"); exit(1);}

    return(0);
}

int linkVerticesToNeighbours(data* anatomy, node** mesh, int nvertices,int GM, int WM){
    int i, k;
    int ngh_count;
    edge* e;
    node* c_vtx, *t_vtx;
    int j=0;
    int zmax=anatomy->zmax;
    int ymax=anatomy->ymax;
    int xmax=anatomy->xmax;
    for (i=0; i< nvertices; i++){
        if(mesh[i] != NULL){
            if( mesh[i]->val > 0 ){
                //printf("%d\n", i % 100);
                //if( i % 100000 == 0 ) printf("i: %2.2f\n", 100*((float)nvertices - i)/nvertices);
                c_vtx=mesh[i];

                //the first neighbour points to itself. this may seem strange but its easier to work with an array that contains all
                //the points we are interested in.
                c_vtx->ngh=malloc( (c_vtx->nngh+1) * sizeof(node*));
                c_vtx->ngh[0]=mesh[i];
               
                connect_nodes(zmax, ymax,xmax, c_vtx, mesh, WM );
                create_edges(c_vtx, mesh);
                //Because GM voxels can only have edges with GM voxels, if they have any outgoing edges then they
                //must be border points
                if(c_vtx->val == GM && c_vtx->n_outgoing_edges >= 1) c_vtx->border=TRUE;
                
            }
        }
    }



    return(0);
}

int length(char *array[]) {
    int i;
    for(i=0; array[i] != NULL; i++);
    return(i);
}


int check_input_files(char *file_inputs[]){
    int n_file_inputs=length(file_inputs);
    int i;
    if(VERBOSE) printf("Number of file inputs: %d\n", n_file_inputs);
    for(i=0; i < n_file_inputs; i++){
        if( access( file_inputs[i], R_OK ) != -1 ) {
            if(VERBOSE) printf("Can read file: %s\n", file_inputs[i]);
        } else {
            printf("Error: could not access %s\n", file_inputs[i]);
            return 1;
        }
    }
    return 0 ;
}

double dist(double a, double b){
    return ((a-b)*(a-b));
}


void find_border(node** inner_border, int ninner,  node** outer_border, int* nouter, float* distances, int* inside, int* temp_outer_border  ){
    *nouter=0;

    for(int i=0; i < ninner; i++){ //For inner border nodes
        node* c_vtx=inner_border[i];
        if(c_vtx->n_outgoing_edges > 26 ) printf("Too many edges!\n");
        for(int j=0; j<c_vtx->n_outgoing_edges; j++){
            edge* e=c_vtx->outgoing_edges[j];
            node* t_vtx=e->end;
            float cur_d=distances[c_vtx->index];
             
            if( inside[t_vtx->index] == FALSE){

                float d=cur_d + 1; //e->length;
                float dt=distances[t_vtx->index];
                if( d < dt || dt == 0){
                    distances[t_vtx->index]=d;
                }
                if( temp_outer_border[t_vtx->index] == FALSE ){
                    //Add target node to outer border
                    *nouter += 1;
                    outer_border[*nouter-1]=t_vtx;
                    temp_outer_border[t_vtx->index]=TRUE;
                }
            }
        }
    }
    for(int i=0; i < *nouter; i++){ 
        inside[outer_border[i]->index]=TRUE;
        temp_outer_border[outer_border[i]->index]=FALSE;
    }

}

void wm_dist(data* img, node** mesh, int* img_vol,  node** gm_border, int n, int WM, int* valid_indices, int nvalid, int start, int step){
    //Calculate inner and outer border voxels
    //Calculate distances to outer border voxels
    //Add border voxels to inner region
    int zmax=img->zmax;
    int ymax=img->ymax;
    int xmax=img->xmax;
    int max =img->n3d;
    int wm_total=0;
    double min_step_size=img->zstep+img->ystep+img->xstep;
    for(int i=0; i<max; i++) if(img_vol[i]==WM) wm_total++;
    int chunk_size = wm_total/1;
    node** inner_border=malloc(sizeof(*inner_border)*wm_total);
    node** outer_border=malloc(sizeof(*outer_border)*wm_total);
    node** tmp=NULL;
    int* inside=calloc(max, sizeof(*inside));
    int* temp_outer_border=calloc(max, sizeof(temp_outer_border));
    float* distances=calloc(max, sizeof(*distances));
    int nouter=0, ninner=0, n_chunks=1;
    int sum=0;
    int counter=0;
    char testfn[100];

    for(int i=start; i < n; i += step){ //Iterate over nodes on WM-GM border
        if(i % 100 == 0) printf("Thread %d:\t%3.3f\n", start, (float) 100. * i/(n/step));
        inner_border[0]=gm_border[i]; //Set initial inner border
        inside[inner_border[0]->index]=TRUE;
        nouter=inner_border[0]->n_outgoing_edges; //Set initial number of outer border nodes
        ninner=1; //Set starting number of inner border to one, because we are starting with a single node

        //printf("%d\n", nouter);
        /**/
        while( nouter > 0){
            find_border(inner_border, ninner, outer_border, &nouter, distances, inside, temp_outer_border);
            tmp=inner_border;
            //The boundary of nodes expands outwards to the outer border.
            inner_border=outer_border;
            //The outer border points to the list in what was previously inner border. This is so that we don't have to reallocate memory for the new outer border points. 
            outer_border=tmp; 
            ninner=nouter;
        }
        for(int j=0; j < nvalid; j++){ 
            int index=valid_indices[j];
            mesh[index]->dist=0; 
            inside[index]=FALSE;
            for(int k=0; k<mesh[index]->n_outgoing_edges; k++){
                edge* e=mesh[index]->outgoing_edges[k];
                e->travelled=FALSE;
            }
        }
        
        //for(int j=0; j < max; j++) if(mesh[j] != NULL) if(mesh[j]->border ==TRUE) distances[mesh[j]->index] += 1;
        if(i==6){ 
            sprintf(testfn,"test_wm.mnc", i );
            writeVolume(testfn, distances, img->start, img->step, img->wcount, MI_TYPE_FLOAT );
            exit(0);
        }
        
    }

}

void* wm_dist_threaded(void* args){
    data* img= ((struct wm_vol_args*) args)->img;
    int WM= ((struct wm_vol_args*) args)->WM;
    int* img_vol=((struct wm_vol_args*) args)->img_vol;
    float* mat=((struct wm_vol_args*) args)->mat;
    node** gm_border=((struct wm_vol_args*) args)->gm_border;
    node** mesh=((struct wm_vol_args*) args)->mesh;
    int n=((struct wm_vol_args*) args)->n;
    int nvalid=((struct wm_vol_args*) args)->nvalid;
    int* valid_indices=((struct wm_vol_args*) args)->valid_indices;
    int thread= ((struct wm_vol_args*) args)->thread;
    int nthreads= ((struct wm_vol_args*) args)->nthreads;
    wm_dist(img,  mesh, img_vol,  gm_border, n,  WM, valid_indices, nvalid, thread, nthreads);
}

int wm_dist_multithreaded(data* img, node** mesh, int* img_vol, node** gm_border, int WM, float* mat, int n, int* valid_indices, int nvalid, int nthreads){
    int rc;
    pthread_t threads[nthreads];
    struct wm_vol_args thread_args[nthreads];
    for(int t=0; t< nthreads; t++){
        thread_args[t].img_vol=img_vol;
        thread_args[t].img=img;
        thread_args[t].WM=WM;
        thread_args[t].mat=mat;
        thread_args[t].n=n;
        thread_args[t].nvalid=nvalid;
        thread_args[t].valid_indices=valid_indices;
        thread_args[t].gm_border=gm_border;
        thread_args[t].mesh=mesh;
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






/***************************************************************************************************
 *Author: Thomas Funck 
 *Contact: thomas.funck@mail.mcgill.ca
 *Date: January 15, 2016
 *Location: Montreal Neurological Institute, Montreal QC
 *
 *Name: surf_dist  
 *Synopsis: - program to calculate minimum distance from a region to any point within anatomically constrained space
 *
 *Inputs: 
 *          1) Mesh composed of vertices in .obj format
 *          2) Binary mask file (.txt) containing "1" for vertices inside region
 *
 *Outputs: 
 *          1) File (.txt) containing distance values at each node
 * 
 *
 *Description:
 *  Calculates distances from 3D region within anatomically constrained region. 
 *
 *
 *Version History:
 *
 *January, 2016
 *   Version-1.0: 
 *
 *
 * ****************************************************************/
int main(int argc, char** argv){
    if (argc < 4 || strcmp(argv[0],"-help") ==0  ) useage();
    
    int i=1;
    float dt=0.1;
    if(strcmp(argv[i++],"-dt") == 0){
        dt=atof(argv[i++]);
    }
    else i=1;
    data anatomy;
    double* dist_vol;
    char* anatomy_filename=argv[i++];
    int GM=atoi(argv[i++]);
    int WM=atoi(argv[i++]);
    int nvertices;
    int* img_vol;
    char *file_inputs[]={anatomy_filename,  NULL}; //={mesh_filename, node_values_filename};
    int nthreads=1; //sysconf(_SC_NPROCESSORS_ONLN);
    float* mat=NULL;
    VERBOSE=FALSE;
    //Check input files to make sure they exist
    if (check_input_files(file_inputs) != 0) exit(1);
    
    anatomy.filename=anatomy_filename;
    img_vol=(int*) readVolume(&anatomy, 1, MI_TYPE_INT);
    dist_vol=malloc(sizeof(*dist_vol) * anatomy.n3d);
    //Read in node locations and values

    if(VERBOSE) printf("Reading in vertices.\n");
    node** mesh=read_vertices_from_volume(&anatomy, img_vol, &nvertices, WM, GM  );
    if(VERBOSE) printf("Number of vertices: %d\n", nvertices);

    if(VERBOSE) printf("Linking to neighbours.\n");
    linkVerticesToNeighbours(&anatomy,  mesh, nvertices, GM, WM);

    if(VERBOSE) printf("GM label: %d\nWM label: %d\n", GM, WM); 
    /*Find GM-WM Border*/
    int n=0, k=0;
    int nvalid=0;
    int* valid_indices=malloc(sizeof(*valid_indices)*anatomy.n3d);
    for(int j=0; j < nvertices; j++) if(mesh[j] != NULL) if(mesh[j]->border == TRUE) n++;
    node** gm_border=malloc(sizeof(*gm_border) * n);
    for(int j=0; j < nvertices; j++){ 
        if(mesh[j] != NULL){
            valid_indices[nvalid++]=j; //Find list of non-NULL indices
            if(mesh[j]->border == TRUE){
                gm_border[k++]=mesh[j];
            }
        }
    }
    int *tmp=realloc(valid_indices, sizeof(*valid_indices) * nvalid );
    valid_indices=tmp;
    /*Find voxel to voxel distances through the WM*/
    //wm_dist(&anatomy, mesh, img_vol, gm_border, n, WM, 0, 1);
    wm_dist_multithreaded(&anatomy, mesh, img_vol, gm_border, WM, mat, n, valid_indices, nvalid, nthreads);

    return 0;
}

void useage(){
    printf("Name:\nwm_dist ~ calculate distances in wm.\n");
    printf("Description:\n");
    printf("Useage:\n<-dt increment> input.mnc GM WM \n");
    exit(1);
}
