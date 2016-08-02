#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <malloc.h>
#include <math.h>
#include "minc2.h"
#include "minc_helper.h"
#define TRUE 1
#define FALSE 0
#define MAX 99999 //FIXME: Only makes sense for brains at the mm scale, may have to be adapted for other applications
int VERBOSE=FALSE;



struct temp_vertex;
struct temp_edge;

/*typedef struct temp_vertex {
  double z;
  double y;
  double x;
  int i, j, k; //voxel indices
  double dist;
  double val;
  double previous_distance;
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
  int inside;
  int border;
  int nngh;
  int index;
  int n;
} vertex;*/


typedef struct temp_edge {
    vertex* start;
    vertex* end;
    double length;
    int travelled;

} edge;

double dist(double a, double b);
int check_if_border(vertex* c_vtx);
void useage();



int length(char *array[]);
int interpolate_sphere( char * input_object_filename, int** n_ngh, int*** ngh ); 

vertex* read_vertices_from_volume(data* anatomy, data* roi, int* nvertices){
    unsigned long i=0; 
    int zmax=anatomy->zmax;
    int ymax=anatomy->ymax;
    int xmax=anatomy->xmax;
    vertex* mesh=malloc(zmax*ymax*xmax*sizeof(vertex));
    if(mesh ==NULL) pexit("Error allocating mesh array", "", 1);
    
    for(int z=0; z < zmax; z++ ){
        for(int y=0; y < ymax; y++ ){
            for(int x=0; x < xmax; x++ ){
                unsigned long index=z*ymax*xmax+xmax*y+x;
                i=index;
                mesh[i].i=x;
                mesh[i].j=y;
                mesh[i].k=z;
                mesh[i].x=voxel2real(x,anatomy->xstart, anatomy->xstep); 
                mesh[i].y=voxel2real(y,anatomy->ystart, anatomy->ystep); 
                mesh[i].z=voxel2real(z,anatomy->zstart, anatomy->zstep); 
                mesh[i].index=i;
                mesh[i].val=roi->data[index];
                mesh[i].dist=0;
                mesh[i].n_incoming_edges=mesh[i].n_incoming_edges=0;
                mesh[i].incoming_edges=mesh[i].outgoing_edges=NULL;
                mesh[i].border=FALSE;
                mesh[i].n=0;
                mesh[i].ngh=NULL;
                mesh[i].nngh=0;
                if(anatomy->data[index] > 0 ) mesh[i].inside_brain=TRUE;
                else mesh[i].inside_brain=FALSE;
                if(mesh[i].val != 0){
                    mesh[i].inside=TRUE;
                    mesh[i].dist=0;
                }
                else {
                    mesh[i].dist=0;
                    mesh[i].inside=FALSE;
                }
            }
        }
    }
    *nvertices=zmax*ymax*xmax;
    return(mesh);
}


int linkVerticesToNeighbours(data* anatomy, vertex* mesh, int nvertices, int max_vertex_length){
    int i, k;
    int ngh_count;
    edge* e;
    vertex* c_vtx, *t_vtx;
    int j=0;
    int zmax=anatomy->zmax;
    int ymax=anatomy->ymax;
    int xmax=anatomy->xmax;
    for (i=0; i< nvertices; i++){
        if( mesh[i].inside_brain == TRUE ){
            //printf("%d\n", i % 100);
            if( i % 100000 == 0 ) printf("i: %2.2f\n", 100*((float)nvertices - i)/nvertices);
            c_vtx=&(mesh[i]);

            //the first neighbour points to itself. this may seem strange but its easier to work with an array that contains all
            //the points we are interested in.
            int x=c_vtx->i;
            int y=c_vtx->j;
            int z=c_vtx->k;
            c_vtx->ngh=malloc( (c_vtx->nngh+1) * sizeof(vertex*));
            c_vtx->ngh[0]=&(mesh[i]);
            
            for(int c=-1; c<=1; c++){
                for(int b=-1; b<=1; b++){
                    for(int a=-1; a<=1; a++){
                        if( (a==0 && b==0 && c==0) ==FALSE) {
                            int zi=z+c;
                            int yi=y+b;
                            int xi=x+a;
                            if( zi >= 0 && zi < zmax && yi >= 0 && yi < ymax && zi >= 0 && zi < zmax){
                                int index=zi*ymax*xmax+yi*xmax+xi;
                                t_vtx=&(mesh[index]); 
                                if( t_vtx->inside_brain == TRUE){
                                    c_vtx->nngh += 1;
                                    c_vtx->ngh=realloc( c_vtx->ngh, (c_vtx->nngh+1) * sizeof(vertex*));
                                    //find vertex that corresponds to this neighbouring vertex
                                   
                                    //because k starts at one we have to subtract k by 1 in ngh, otherwise we spill over into
                                    //adjascent arrays
                                    //first check that the neighbour actually exists
                                    c_vtx->ngh[ c_vtx->nngh]=t_vtx;
                                }
                            }
                        }
                    }
                }
            }
            //printf("Okay: %d\n", c_vtx->nngh);
            for(int k=0; k< c_vtx->nngh; k++ ){
                    t_vtx=&( mesh[ c_vtx->ngh[k]->index ] );
                    //if we find a neighbouring vertex that is outside
                    if(t_vtx->inside != TRUE){
                        double length=sqrt( dist(c_vtx->x, t_vtx->x) + dist(c_vtx->y, t_vtx->y) + dist(c_vtx->z, t_vtx->z)  );
                        if (length < max_vertex_length) {
                            //if the vertex mesh[i] is connected to a neigbour that is outside
                            //then create an edge connecting the two
                            e=malloc(sizeof(edge));
                            e->start=c_vtx;
                            e->end=t_vtx;
                            e->travelled=FALSE;
                            e->length=length;
                            c_vtx->outgoing_edges=realloc(c_vtx->outgoing_edges, sizeof(edge*) * (1+c_vtx->n_outgoing_edges));
                            c_vtx->outgoing_edges[c_vtx->n_outgoing_edges]=e;
                            c_vtx->n_outgoing_edges += 1;
                            t_vtx->incoming_edges=realloc(t_vtx->incoming_edges, sizeof(edge*) * (1+t_vtx->n_incoming_edges));
                            t_vtx->incoming_edges[t_vtx->n_incoming_edges]=e;
                            t_vtx->n_incoming_edges += 1;

                        } 
                        //set current vertex's border flag to true
                        if(c_vtx->inside == TRUE){ 
                            c_vtx->border=TRUE;
                            //printf("In: %d\tBorder:\t%d\n", c_vtx->inside, c_vtx->border);
                        }
                    }
                    //increment number of outgoing and incoming edges
                }
            }
        }

    for (i=0; i< nvertices; i++) if(mesh[i].inside_brain==TRUE) mesh[i].border=check_if_border(&(mesh[i]));


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


int find_euclidian_distance(vertex** mesh, int nvertices, int inside_label){
    double x_sum=0, y_sum=0, z_sum=0;
    double euclidian_distance;
    int n=0;
    vertex* c_vtx; 
    for(int i =0; i<nvertices; i++){
        c_vtx=mesh[i];
        if(c_vtx->val == inside_label ){
            x_sum += c_vtx->x;
            y_sum += c_vtx->y;
            z_sum += c_vtx->z;
            n+=1;
        }
    }
    x_sum /= n;
    y_sum /= n;
    z_sum /= n;

    for(int i =0; i<nvertices; i++){
        c_vtx=mesh[i];
        if(c_vtx->val != inside_label ){
            euclidian_distance=sqrt(dist(c_vtx->x, x_sum) + dist(c_vtx->y, y_sum) + dist(c_vtx->z, z_sum));
            c_vtx->dist = c_vtx->dist- euclidian_distance;
        }
    }

return(0);
}


int check_if_border(vertex* c_vtx){
    vertex* t_vtx;

    //update border voxels
    for(int j=0; j<c_vtx->nngh; j++){
        t_vtx=c_vtx->ngh[j];
            if( t_vtx->inside == TRUE){ 
                if(c_vtx->inside != TRUE )
                    //printf("target is inside, current is outside\n");
                    //Current is outside, target is inside, therefore current is border point.
                    return(TRUE);
            }
            else {
                if(c_vtx->inside ==TRUE){
                //Current vertex is inside, target is outside, therefore current is border
                    return(TRUE);
                }   
            }  
    }
    return(FALSE);
}

int add_to_list(void* list, int list_size, int* n, void* element ){
    *n += 1; //increase the total number of elements in the list by 1
    if(list_size==sizeof(double) ){ //If we are dealing with a double.
        double** list_d=list; //list d is a double pointer because "list" is the address of a pointer, so you have to dereference it twice to get to what its point to.
        double* element_d= element; //a double pointer that is used to get to what "element" is pointing to
        *list_d=realloc(*list_d, sizeof(**list_d) * *n); 
        if( *list_d ==NULL) return(1);
        (*list_d)[*n-1]=*element_d;
    } 

    return(0);
}

double projection(vertex* first_vertex, vertex* second_vertex,vertex* third_vertex, double* x, double* y, double* z){
    //Derivation of shortest distance from first_vertex to edge between second and third vertex:
    //
    //1) D=(x2-x0-t(x1-x0))^2 + (y2-y0-t(y1-y0))^2 + (z2-z0-t(z1-z0))^2
    //2) dD/dt = 2(x2-x0-t(x1-x0))(x1-x0) + 2(y2-y0-t(y1-y0))(y1-y0) + 2(z2-z0-t(z1-z0))(z1-z0) = 0
    //3) 2(x2-x0)(x1-x0)-2*t(x1-x0)(x1-x0) + 2(y2-y0)(y1-y0)-2*t(y1-y0)(y1-y0) + 2(z2-z0)(z1-z0)-2*t(z1-z0)(z1-z0) =0
    //4) (x2-x0)(x1-x0) + (y2-y0)(y1-y0)+ (z2-z0)(z1-z0) =t(x1-x0)(x1-x0) + t(y1-y0)(y1-y0) + t(z1-z0)(z1-z0)
    //
    //      (x2-x0)(x1-x0) + (y2-y0)(y1-y0)+ (z2-z0)(z1-z0)
    //5) t=---------------------------------------------------
    //          (x1-x0)^2 + (y1-y0)^2 + (z1-z0)^2
    //       
    double x0=second_vertex->x;//inside
    double y0=second_vertex->y;
    double z0=second_vertex->z;
    double x1=third_vertex->x;//inside/
    double y1=third_vertex->y;
    double z1=third_vertex->z;   
    double x2=first_vertex->x;//outside
    double y2=first_vertex->y;
    double z2=first_vertex->z;   
    double num= (x2-x0)*(x1-x0)+(y2-y0)*(y1-y0)+(z2-z0)*(z1-z0);
    double den= (x1-x0)*(x1-x0) + (y1-y0)*(y1-y0) + (z1-z0)*(z1-z0);
    double t=num/den;
    if( t< 0) t=0;
    else if (t >1) t=1; 
    *x=x0+t*(x1-x0);
    *y=y0+t*(y1-y0);
    *z=z0+t*(z1-z0);
    //printf("\tSecond\tThird\tFirst\n");
    //printf("x %f %f, %f --> %f\n", x0, x1, t , *x);
    //printf("y %f %f, %f --> %f\n", y0, y1, t , *y);
    //printf("z %f %f, %f --> %f\n", z0, z1,t , *z);
    return(0);
}

double interpolate_distance(vertex* vertex1, vertex* vertex2, double x, double y, double z){
    double x0=vertex1->x;
    double y0=vertex1->y;
    double z0=vertex1->z;
    double x1=vertex2->x;
    double y1=vertex2->y;
    double z1=vertex2->z; 
    double distance1=vertex1->dist;
    double distance2=vertex2->dist;
    double total_distance=sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)+(z1-z0)*(z1-z0));
    double part_distance=sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0)+(z-z0)*(z-z0));
    //when p==p0 -> part_distance==0 --> second_distance==0
    double out=distance1*(1-part_distance/total_distance )+distance2*(part_distance/total_distance);

    return(out);
}

int is_neighbour(vertex* current, vertex* target){
    for(int i=0; i< current->nngh; i++){
       vertex* ngh=current->ngh[i];
       if(ngh == target){
           return(TRUE);
       }
    }
    return(FALSE);
}

double find_smallest_distance(vertex* first_vertex){
    vertex* second_vertex, *third_vertex;
    double* distance_list=NULL;
    double distance;
    double min_distance;
    int ndistances;
    double x, y, z;
    int n_virtual;
    ndistances=0;

    //printf("First: %f %f %f\n", first_vertex->x, first_vertex->y, first_vertex->z);
    //check the neighbouring vertices around the first vertex
    for(int j=0; j < first_vertex->nngh; j++){
        second_vertex=first_vertex->ngh[j];
        //if any of the surroudning vertices are inside
        if(second_vertex->inside == TRUE){
            //printf("Second: %f %f %f\n", second_vertex->x, second_vertex->y, second_vertex->z);
            //then check the neighbours around the second vertex to find another vertex in common with the first one
            for(int i=0; i < second_vertex->nngh; i++){
                third_vertex=second_vertex->ngh[i];
                if ( third_vertex->inside == TRUE && is_neighbour(first_vertex, third_vertex)==TRUE  ){
                    //printf("Third: %f %f %f\n", third_vertex->x, third_vertex->y, third_vertex->z);
                    projection(first_vertex, second_vertex, third_vertex, &x, &y, &z);
                    distance=interpolate_distance(second_vertex, third_vertex, x, y, z);
                    distance += sqrt(dist(x, first_vertex->x)+dist(y, first_vertex->y)+dist(z, first_vertex->z));
                    if ( add_to_list(&distance_list, sizeof(distance), &ndistances, &distance ) != 0 ) 
                        //pexit("Error: could not allocate memory in <add_to_list>.",  __FILE__, __LINE__, 1 );
                        pexit("Error: could not allocate memory in <add_to_list>.", "", 1 );
                    //printf("Projection: %d, %f\n", ndistances, distance );
                }
            }
            
            distance=sqrt(dist(second_vertex->x, first_vertex->x)+dist(second_vertex->y, first_vertex->y)+dist(second_vertex->z, first_vertex->z));

            distance += second_vertex->dist;
            if( add_to_list(&distance_list, sizeof(distance), &ndistances, &distance) != 0 ) 
                pexit("Error: could not allocate memory in <add_to_list>.", "", 1 );
            //printf("Edge: %d, %f\n", ndistances, distance );
        }
    }
    
    min_distance=distance_list[0];
    for(int i=1; i< ndistances; i++){
        //printf("%d: %f\n", i, distance_list[i]);
        if( min_distance > distance_list[i]) min_distance=distance_list[i];
    }
    //printf("Minimum: %f\n", min_distance);
    free(distance_list);
    return(min_distance);
}


//need to initialize edges

int find_distances(vertex* mesh,double* dist_vol, int n_vtx,  float dt){
    //For time "t" with increment "dt"
    float t=0;
    int n_outer_vtx=-1;
    vertex *c_vtx, *t_vtx, *p_vtx;
    edge* e;
    double distance, cp_dist;
    int var;
    int skip_alt=0;
    int n_border_vtx=-1;
    int counter=0;
    //while(n_outer_vtx != 0 ){
    //while the number of border points is not equal to 0 (

    while(n_border_vtx != 0 ){
        if (counter % 100 == 0) printf("%3.2f\t%d\t%d\n",t, n_outer_vtx, n_border_vtx);
        //Set the number of outer vertices to 0
        n_outer_vtx=n_border_vtx=0;
        //increment t by dt (keeps track of total distance travelled along edges)
        t += dt;
        //for all vertiex points
        for(int i=0; i< n_vtx; i++){
            if( mesh[i].inside_brain == TRUE){
                c_vtx=&(mesh[i]); //current vertex "c_vtx"
                //Count the number of vertices that are not within the region of interest (ROI) (i.e., they are exterior) and
                //are connected by edges to other vertices 
                if (c_vtx->inside != TRUE && c_vtx->n_incoming_edges > 0) n_outer_vtx += 1;
                //Check if c_vtx is an interior, border point
                if( c_vtx->inside == TRUE && c_vtx->border==TRUE){
                    var=0;
                    //for edges leaving the current vertex
                    for(int j=0; j < c_vtx->n_outgoing_edges; j++){
                        e=c_vtx->outgoing_edges[j];
                        t_vtx=e->end; //the target vertex is the one at the end of the current edge
                        //if this edge leads to an outer vertex, check if we have travelled
                        //to this outer vertex
                        if(e->travelled != TRUE && e->end->inside != TRUE ){
                            var += 1;
                            //the total distance to this outer vertex is the
                            //distance travelled to the current vertex plus the edge 
                            //length
                            
                            distance =  c_vtx->dist + e->length; //find_smallest_distance(t_vtx);
                            //printf("%f %f\n", distance, t );
                            if (distance <= t){
                                //if t is greater than or equal to distance,
                                //then the edge has been travelled
                                e->travelled=TRUE; //mark edge as travelled
                                if(skip_alt == 1){
                                    t_vtx->previous=c_vtx;
                                }
                                    //if there are no more edges connecting to the target vertex
                                    //then it goes from being outside to inside. 
                                    //Its exterior neighbours are also added to the list of frontier points
                                t_vtx->inside=TRUE;
                                t_vtx->border=TRUE;
                                for(int k=1; k < t_vtx->n_incoming_edges; k++){  
                                    t_vtx->incoming_edges[k]->travelled=TRUE;
                                }
                                dist_vol[t_vtx->index]=distance; 
                                t_vtx->dist=distance; 
                            }
                        }
                    }
                    //There are no more untravelled edges for this vertex
                    //so it is no longer a border point
                    if(var == 0) c_vtx->border=FALSE;
                    else n_border_vtx++;
                }
            }
        }
    counter++;
    }
    //for(int i=0; i< n_vtx; i++) if(mesh[i]->dist == MAX ){
    //mesh[i]->dist=0;
    //}
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
 *          1) File (.txt) containing distance values at each vertex
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
 *
 *
 *
 *
 *
 *
 * ****************************************************************/
int main(int argc, char** argv){
    if (argc < 4 || strcmp(argv[0],"-help") ==0  ) useage();
    
    int i=1;
    float dt=0.1;
    if(strcmp(argv[i++],"-dt") == 0){
        dt=atof(argv[i++]);
        printf("%f\n");
    }
    else i=1;
    data anatomy, roi;
    double* dist_vol;
    char* anatomy_filename=argv[i++];
    char* roi_filename=argv[i++]; 
    char* dist_vol_filename=argv[i++];
    int nvertices, n_border_vtx,n_inner_vtx, n_outer_vtx, n_frontier_vtx=-1;
    int *n_ngh, **ngh;
    int max_vertex_length=4;
    char *file_inputs[]={anatomy_filename, roi_filename, NULL}; //={mesh_filename, vertex_values_filename};
    float inside_label=1, outside_label=0;
    vertex  **border, **inside, **outside, **frontier;
    VERBOSE=TRUE;
    //Check input files to make sure they exist
    if (check_input_files(file_inputs) != 0) exit(1);
    
    anatomy.filename=anatomy_filename;
    roi.filename=roi_filename;
    anatomy.data=(double*) readVolume(&anatomy, 1, MI_TYPE_DOUBLE);
    roi.data=(double*)readVolume(&roi, 1, MI_TYPE_DOUBLE);
    dist_vol=malloc(sizeof(*dist_vol) *anatomy.xmax*anatomy.ymax*anatomy.zmax);
    //Read in vertex locations and values
    if(VERBOSE) printf("Reading in vertices.\n");
    vertex* mesh=read_vertices_from_volume(&anatomy, &roi, &nvertices  );
    if(VERBOSE) printf("Number of vertices: %d\n", nvertices);

    //Link vertex ngh together
    if(VERBOSE) printf("Linking to neighbours.\n");
    linkVerticesToNeighbours(&anatomy,   mesh, nvertices, max_vertex_length);
    printf("Find distances\n");
    find_distances( mesh, dist_vol, nvertices,  dt);


    writeVolume(dist_vol_filename, dist_vol, anatomy.start,  anatomy.step, anatomy.wcount, MI_TYPE_DOUBLE  );

    return 0;
}

void useage(){
    printf("Name:\nsurf_dist ~ calculate distances along mesh surface.\n");
    printf("Description:\nFor a given mask defined on a surface mesh, calculate the minimum distance from thie region to every vertex exterior ");
    printf("to this region on the mesh.\n");
    printf("Useage:\n<-dt> increment> input_mesh.obj vertex_values.txt output.txt\n");
    exit(1);
}
