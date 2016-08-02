#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <malloc.h>
#include <math.h>
#include "minc2.h"
#define MAX 99999

#define TRUE 1
#define FALSE 0
int VERBOSE=FALSE;
int interpolate_sphere( char * input_object_filename, int** n_ngh, int*** ngh ); 
double dist(double a, double b);

struct temp_vertex;
struct temp_edge;

typedef struct temp_vertex {
  double z;
  double y;
  double x;
  //for voxels the pointer below should 
  //point to the nearest vertex. For vertices
  //the pointer should go to the ngh 
  //of that vertex
  struct temp_vertex *previous;
  //the info below only applies to vertices
  //on a mesh, not to voxels
  struct temp_vertex** ngh;
  struct temp_edge** incoming_edges;
  struct temp_edge** outgoing_edges;
  int n_incoming_edges;
  int n_outgoing_edges;
  float label;
  int nngh;
  int index;
  int visited;
} vertex;

struct group {
    vertex **set;
    int n;
} ;

struct freq_array {
    float value;
    int freq;
};

typedef struct temp_edge {
    vertex* start;
    vertex* end;
    double length;
    int travelled;

} edge;

void pexit(char* string, char* string2, int number) {
    printf("%s %s\n", string, string2); 
    exit(number);
}

double dist(double a, double b){
    return ((a-b)*(a-b));
}

vertex** read_vertices(char* mesh_filename, char* vertex_values_filename, int* nvertices, float inside_label){
    char buffer1[100],  buffer2[100]; //should probably improve buffer length so that it is not a fixed size
    int i; 
    FILE *mesh_file=fopen(mesh_filename, "rt");
    FILE *vertex_values_file=fopen(vertex_values_filename, "rt");
    char dlm[]=" "; 
    vertex* mesh;
    vertex** mesh_ptr;


    if(vertex_values_file==NULL) {
        printf("Error opening %s\n", mesh_filename);
        exit(1);
    }
    if(mesh_file==NULL) {
        printf("Error opening %s\n", mesh_filename);
        exit(1);
    }

    fgets(buffer1, sizeof(buffer1), mesh_file) ;
    //read nvertices from file
    strtok(buffer1, dlm);  
    strtok(NULL, dlm); 
    strtok(NULL, dlm);  
    strtok(NULL, dlm); 
    strtok(NULL, dlm); 
    strtok(NULL, dlm); 
    *nvertices=atoi(strtok(NULL, dlm));

    mesh=malloc(sizeof(vertex) * *nvertices);
    mesh_ptr=malloc(sizeof(vertex*) * *nvertices);
    for(i=0; i< *nvertices; i++){
        fgets(buffer1, sizeof(buffer1), mesh_file);
        fgets(buffer2, sizeof(buffer2), vertex_values_file); 
        mesh[i].x=atof(strtok(buffer1, dlm)); 
        mesh[i].y=atof(strtok(NULL, dlm)); 
        mesh[i].z=atof(strtok(NULL, dlm)); 
        mesh[i].index=i;
        mesh[i].label=atof(buffer2);
        mesh[i].n_incoming_edges=mesh[i].n_incoming_edges=0;
        mesh[i].incoming_edges=mesh[i].outgoing_edges=NULL;
        mesh[i].visited=FALSE;
        mesh_ptr[i]=&mesh[i];
    }

    fclose(mesh_file);
    fclose(vertex_values_file);
    return(mesh_ptr);
}


int linkVerticesToNeighbours(vertex** mesh, int* n_ngh, int** ngh, int nvertices, int max_vertex_length){
    int i, k;
    int ngh_count;
    edge* e;
    vertex* c_vtx, *t_vtx;
    
    for (i=0; i< nvertices; i++){
        c_vtx=mesh[i];
 	    c_vtx->ngh=malloc(sizeof(vertex*) * (n_ngh[i]) );
	    c_vtx->nngh=n_ngh[i];
	    //the first neighbour points to itself. this may seem strange but its easier to work with an array that contains all
	    //the points we are interested in.

	    for(k=0; k< n_ngh[i]; k++){
		    //because k starts at one we have to subtract k by 1 in ngh, otherwise we spill over into
		    //adjascent arrays
            //first check that the neighbour actually exists
            if(ngh[i][k] < nvertices  ){
		        c_vtx->ngh[k]=mesh[ ngh[i][k] ];
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


int search_neighbours(vertex* c_vtx, struct group* vertex_group, int inside_label){
    vertex* t_vtx; 
     c_vtx->visited=TRUE;//take vertex out of exterior points
    //make sure outgoing edges are pointing to EXTERIOR vertices
    //printf("%d\n",c_vtx->n_outgoing_edges );
    for(int i=0; i< c_vtx->nngh; i++){
        t_vtx=c_vtx->ngh[i];
        if( t_vtx->visited != TRUE && t_vtx->label == inside_label ){
            vertex_group->n += 1;
            vertex_group->set=realloc(vertex_group->set, sizeof(*vertex_group->set) * vertex_group->n);
            if(vertex_group->set == NULL) pexit("Error:", "could not allocate enough memory in <search_neighbours>", 1);
            vertex_group->set[vertex_group->n -1]=t_vtx;
            search_neighbours(t_vtx, vertex_group, inside_label);
        }
    } 

    return(0);
}



struct group* find_holes(vertex** mesh, int n_vtx, int* n_groups, int inside_label ){
    struct group *vertex_groups=NULL;
    vertex* c_vtx;
    *n_groups=0;

    //for all vertices in the mesh
    for(int i=0; i< n_vtx; i++){
        c_vtx=mesh[i]; //set current vertex to mesh[i]
        if (c_vtx->label == inside_label && c_vtx->visited != TRUE) {
            //We have found an inner vertex, so we will create a new
            //set of pointers to this vertex and its surrounding neighbours
            *n_groups += 1;
            vertex_groups=realloc(vertex_groups, sizeof(*vertex_groups) * *n_groups ); 
            vertex_groups[*n_groups-1].set=malloc(sizeof(*(vertex_groups[*n_groups-1].set)));
            vertex_groups[*n_groups-1].set[0]=c_vtx;
            vertex_groups[*n_groups-1].n=1;
            //search for the current vertex's surounding neighbours that
            //are also outside. 
            search_neighbours(c_vtx, &(vertex_groups[*n_groups-1]), inside_label );
            //if(VERBOSE) printf(" %d\n", vertex_groups[*n_groups-1].n);
        }
    }
    return(vertex_groups);
}



int value_in_array(vertex* t_vtx, struct freq_array* array, int narray){

    for (int i=0; i< narray; i++){   
        if(t_vtx->label == array[i].value){
            array[i].freq++;
            return(0);
            }
        }   

return(1);
}


float most_freq_label(vertex* c_vtx, float inside_label, float* output){
    vertex* t_vtx;
    struct freq_array* array=malloc(sizeof(*array)), max;
    int narray=0;
    max.freq=0;

    for(int k=0; k< c_vtx->nngh; k++){
        t_vtx=c_vtx->ngh[k];
        if( t_vtx->label != inside_label ){
            if( value_in_array(t_vtx, array, narray) != 0){
                narray++;
                array=realloc(array, sizeof(*array)* narray);
                array[narray-1].freq=1;
                array[narray-1].value=t_vtx->label;
            }
        }
    }

    if(narray > 0){
        for( int i=0; i<narray; i++){
            if(max.freq < array[i].freq){
                max.freq=array[i].freq;
                max.value=array[i].value;
            }
        }
    *output=max.value;
    free(array);
    return(0);
    }
    else return(1);
}



int fill_holes(struct group* vertex_groups, float threshold, int n_groups, int inside_label){
    int max_vertices=0; 
    float vertex_threshold;
    vertex* c_vtx, *t_vtx;
    float new_label;
    int found_inside_vertex;

    for(int i=0; i < n_groups; i++){
        if( max_vertices < vertex_groups[i].n) {
            max_vertices=vertex_groups[i].n;
        }
    }
    vertex_threshold = max_vertices * threshold;
    for(int i=0; i< n_groups; i++){
        if ( vertex_groups[i].n < vertex_threshold ){
            printf("Filling in group of size: %d\n", vertex_groups[i].n );
            //While some vertices are not yet
            found_inside_vertex=TRUE;
            while( found_inside_vertex==TRUE ){
                found_inside_vertex=FALSE;
                for(int j=0; j<vertex_groups[i].n; j++){
                    c_vtx=vertex_groups[i].set[j];
                    if( c_vtx->label == inside_label){
                        if(most_freq_label(c_vtx, inside_label, &new_label) ==0){
                            
                            c_vtx->label=new_label;
                        }
                        found_inside_vertex=TRUE;
                    }
                        
                }
            }
        }
    }

    return(0);
}

int main(int argc, char** argv){
    if (argc != 6 ){
        printf("surf_defrag - defragment object file based\n");
        printf("Useage: input_mesh.obj  mesh_values.txt  inside_label threshold  output.txt\n");
        exit(1);
    }
    int i=1;
    char* mesh_filename=argv[i++];
    char* vertex_values_filename=argv[i++]; 
    float inside_label=atoi(argv[i++]); //value for label we want to close
    float threshold=atof(argv[i++]);
    FILE* outputfile=fopen(argv[i++], "w+");
    int nvertices, n_border_vtx,n_inner_vtx, n_outer_vtx, n_frontier_vtx=-1;
    int *n_ngh, **ngh;
    int* output_labels;
    int max_vertex_length=10;
    char *file_inputs[]={mesh_filename, vertex_values_filename, NULL}; //={mesh_filename, vertex_values_filename};
    int n_groups=0;
    struct group *vertex_groups=NULL;
    
    //Check input files to make sure they exist
    if (check_input_files(file_inputs) != 0) exit(1);
    
    //Read in vertex locations and values
    if(VERBOSE) printf("Reading in vertices.\n");
    vertex** mesh=read_vertices(mesh_filename,  vertex_values_filename, &nvertices, inside_label  );
    if(VERBOSE) printf("Number of vertices: %d\n", nvertices);
    output_labels=malloc(sizeof(*output_labels) * nvertices);
    for(int i=0; i<nvertices; i++) output_labels[i]=0;
    //Find vertex ngh
    if(VERBOSE) printf("Interpolate sphere.\n");
    interpolate_sphere( mesh_filename , &n_ngh, &ngh);
    
    //Link vertex ngh together
    if(VERBOSE) printf("Linking to neighbours.\n");
    linkVerticesToNeighbours(mesh,  n_ngh, ngh,  nvertices, max_vertex_length);

    
    vertex_groups=find_holes( mesh, nvertices, &n_groups, inside_label);
    fill_holes(vertex_groups, threshold, n_groups, inside_label);

    float out;
    for(int i=0; i< nvertices; i++) {
        fprintf( outputfile, "%d\n",(int) round(mesh[i]->label));
    }
    return 0;
}


