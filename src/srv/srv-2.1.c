#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include "minc2.h"
#include "meshPVC_for_srv.h"
#include "minc_helper.h"


void linkVerticesToNeighbours(  vector* mesh, int* n_ngh, int** ngh, int nvertices){
  int i, k;

 for (i=0; i< nvertices; i++){
 	mesh[i].neighbours=malloc(sizeof(vector*) * (n_ngh[i]+1) );
	mesh[i].nneighbours=n_ngh[i]+1;
	mesh[i].neighbours[0]=&mesh[i];
    mesh[i].nvox=0;
    mesh[i].voxNgh=NULL;
	//the first neighbour points to itself. this may seem strange but its easier to work with an array that contains all
	//the points we are interested in. 
	for(k=1; k< n_ngh[i]+1; k++){
		//because k starts at one we have to subtract k by 1 in ngh, otherwise we spill over into
		//adjascent arrays
		mesh[i].neighbours[k]=&mesh[ ngh[i][k-1] ];
	} 
 }

};

int readData(char* outerMeshfilename, char* innerMeshfilename, int* nvertices, vector** inner_mesh,  vector** outer_mesh  ){
char buffer1[100], buffer2[100], buffer3[100], buffer4[100], buffer5[100], buffer6[100]; //should probably improve buffer length so that it is not a fixed size
int i; 
int vertices;
FILE *outerMeshFile;
FILE *innerMeshFile;
FILE *outer_activityFile;
FILE *inner_activityFile;
char *token;
char dlm[]=" "; 
int xmin=0, ymin=0, zmin=0;
int nodeCounter=0;

outerMeshFile=fopen(outerMeshfilename, "rt");
innerMeshFile=fopen(innerMeshfilename, "rt"); 
fgets(buffer1, sizeof(buffer1), outerMeshFile) ;
fgets(buffer2, sizeof(buffer2), innerMeshFile);

//read nvertices from file
strtok(buffer1, dlm);  
strtok(NULL, dlm); 
strtok(NULL, dlm);  
strtok(NULL, dlm); 
strtok(NULL, dlm); 
strtok(NULL, dlm); 
*nvertices=atoi(strtok(NULL, dlm));
*inner_mesh=malloc(sizeof(vector) * *nvertices);
*outer_mesh=malloc(sizeof(vector) * *nvertices);

for(i=0; i< *nvertices; i++){
	fgets(buffer1, sizeof(buffer1), outerMeshFile);
	fgets(buffer2, sizeof(buffer2), innerMeshFile);
       		 
        /*OUTER*/
        (*outer_mesh)[i].x=atof(strtok(buffer1, dlm)); 
        (*outer_mesh)[i].y=atof(strtok(NULL, dlm)); 
        (*outer_mesh)[i].z=atof(strtok(NULL, dlm)); 
        (*outer_mesh)[i].layer=2; //2 == outer layer
        (*outer_mesh)[i].index=i;
        /*INNER*/
        (*inner_mesh)[i].x=atof(strtok(buffer2, dlm)); 
        (*inner_mesh)[i].y=atof(strtok(NULL, dlm)); 
        (*inner_mesh)[i].z=atof(strtok(NULL, dlm)); 
        (*inner_mesh)[i].layer=0; //0 == inner layer
        (*inner_mesh)[i].index=i;
        /*Mid*/
	(*outer_mesh)[i].pair=&(*inner_mesh)[i];
	(*inner_mesh)[i].pair=&(*outer_mesh)[i];

}
 

fclose(outerMeshFile);
fclose(innerMeshFile);
return(0);
}

void voxel_vertexPairs(vector** vertices , double maxima[3], double minima[3], double step[3], float *voxelArray,  double global_xmin, double global_ymin, double global_zmin, vector vtx_min, vector vtx_max, misize_t* count ){
    float xmin=minima[2], ymin=minima[1], zmin=minima[0];
    float xmax=maxima[2], ymax=maxima[1], zmax=maxima[0];
    float xstep = fabsf(step[2]);
    float ystep = fabsf(step[1]);
    float zstep = fabsf(step[0]);
    float x, y, z;
    int u, v, w;
    int i, j, invoxels=0;
    vector testVoxel;
    vector** inner_mesh, **outer_mesh;
    float min; 
    int index;
    int signal=0;

    //printf("%f\t%f\t%f\n", xmin, xstep, xmax);
    //printf("%f\t%f\t%f\n", ymin, ystep, ymax);
    //printf("%f\t%f\t%f\n", zmin, zstep, zmax);
    for(z=zmin; z <= zmax /*+ zstep*/; z += zstep){
	    for(y=ymin; y <=  ymax/* + ystep*/ ; y += ystep){
            for(x=xmin; x <=  xmax /*+ xstep*/; x += xstep){
                testVoxel.x=x;
                testVoxel.y=y;
                testVoxel.z=z;                
                testVoxel.nearest=vertices[0]; 
               	if(testVoxel.nearest->layer == 0){
                    //we need to know if the vertex we found is on the inner
                    //or outer layer (actually i'm not sure that this is true) and
                    //then get its corresponding pair on the opposite layer
                    inner_mesh=testVoxel.nearest->neighbours;
                    outer_mesh=testVoxel.nearest->pair->neighbours;
                }
                else{
                    outer_mesh=testVoxel.nearest->neighbours;
                    inner_mesh=testVoxel.nearest->pair->neighbours;
                }
                //printf("Vox: %f %f %f \n", z, y, x );
                for(i=0; i < vertices[i]->nneighbours ; i++ )
                    //printf("Vtx: %f %f %f \n", vertices[i]->z, vertices[i]->y, vertices[i]->x );

		        testVoxel.inside=pointinMesh( &testVoxel, inner_mesh, outer_mesh,  testVoxel.nearest->nneighbours );
			    if(testVoxel.inside == 1){
			        //Convert real voxel location to index values
                    u=(int) round((testVoxel.x-global_xmin)/step[2]);  
			        v=(int) round((testVoxel.y-global_ymin)/step[1]);  
			        w=(int) round((testVoxel.z-global_zmin)/step[0]); 
                    index=w*count[1]*count[2]+v*count[2]+u; 
			        voxelArray[index]= (unsigned short) 1;
			        invoxels++;
                }
                //printf("Not inside\n");

		    }
        }
    }

    //if(invoxels == 0 ){ 
    //    printf("Number of Voxels in Mesh: %d\n", invoxels);
    //}   

}

int main(int argc, char** argv){         
    if(argc != 6) {
        printf("Usage: <MRI file> <Outer surface vertices> <Inner surface vertices> <Scaling Factor>  <Output MINC file> \n");
        exit(1);
    }        
    int argi=1, i=0, x, y, z, w;          
    int nvertices=0;
    int *n_ngh_inner, **ngh_inner;
    int *n_ngh_mid, **ngh_mid;
    int *n_ngh_outer, **ngh_outer;
    vector* inner_mesh, *mid_mesh, *outer_mesh;
    vector min, max;
    vector vox_min, vox_max;
    int maxDepth=2, counter, surf_nvertices, loc_nvertices;
    int currentDepth=0;
    vector** vertices=NULL;
    float bounds[8][3];
    double maxima[3];
    double minima[3], start[3];	
    misize_t native_count[3];
    double xmin, ymin, zmin;
    double native_step[3];
    double zmax, ymax, xmax;
    char* mrifilename=argv[argi++];
    char* outer_surface_file=argv[argi++];
    char* inner_surface_file=argv[argi++];
    unsigned long scaleFactor=(unsigned long) atoi(argv[argi++]); //must be an integer to keep the count an integer
    char* outputfilename=argv[argi++];
    double step[3];
    misize_t count[3];
    data vol;
   

    vol.filename=mrifilename;
    readVolume( &vol, 0, MI_TYPE_DOUBLE );
    
    for(i=0; i<3; i++){
        start[i]=vol.start[i];
        native_count[i]=vol.count[i];
        native_step[i]=vol.step[i];
    }
    

    step[0]=(double)native_step[0]/scaleFactor;
    step[1]=(double)native_step[1]/scaleFactor; 
    step[2]=(double)native_step[2]/scaleFactor;
    zmin=start[0];
    ymin=start[1];
    xmin=start[2];
    zmax=minima[0]+native_count[0]*native_step[0]; 
    ymax=minima[1]+native_count[1]*native_step[1]; 
    xmax=minima[2]+native_count[2]*native_step[2];
    count[0]=native_count[0]*scaleFactor; 
    count[1]=native_count[1]*scaleFactor; 
    count[2]=native_count[2]*scaleFactor;
    int nvox=native_count[0]*native_count[1]*native_count[2];
   //printf("Counts: %d %d %d = %d\n", count[0], count[1], count[2], count[0]*count[1]*count[2]*sizeof(float)); 
   // printf("Starts: %f %f %f \n", zmin, ymin, xmin); 
    float* array=calloc(nvox, sizeof(*array));
    
	readData(outer_surface_file , inner_surface_file ,   &nvertices, &inner_mesh,   &outer_mesh  );

    interpolate_sphere( outer_surface_file , &n_ngh_outer, &ngh_outer);
    interpolate_sphere( inner_surface_file , &n_ngh_inner, &ngh_inner);
	//interpolate_sphere( argv[4], &n_ngh_mid, &ngh_mid);

	linkVerticesToNeighbours( inner_mesh, n_ngh_inner, ngh_inner, nvertices);
	linkVerticesToNeighbours( outer_mesh, n_ngh_outer, ngh_outer, nvertices);
	//linkVerticesToNeighbours( mid_mesh, n_ngh_mid, ngh_mid, nvertices);

for(counter=0; counter < nvertices; counter++){
	surf_nvertices=outer_mesh[counter].nneighbours;
	loc_nvertices=2*surf_nvertices;
	vertices=malloc( loc_nvertices * sizeof(vector*));
	if(counter % 1000 == 0) printf("Super Resolution: %f\n", 100*(float) counter / nvertices );

	for(i=0; i < surf_nvertices; i++){
	  vertices[i]=inner_mesh[counter].neighbours[i];
	  vertices[i+surf_nvertices]=outer_mesh[counter].neighbours[i];
	}


	//Find minimum and maximum for a set of neighbouring vertices (both inner and outer vertices)
    for(i=0; i<loc_nvertices; i++){
		if( i ==0){
		min.x=vertices[0]->x;
		min.y=vertices[0]->y;
		min.z=vertices[0]->z;
		max.x=vertices[0]->x;
		max.y=vertices[0]->y;
		max.z=vertices[0]->z;
		}
		else if(vertices[i]->x < min.x) min.x=vertices[i]->x;  
		else if(vertices[i]->y < min.y) min.y=vertices[i]->y;  
		else if(vertices[i]->z < min.z) min.z=vertices[i]->z;  
		if(vertices[i]->x > max.x) max.x=vertices[i]->x;  
		if(vertices[i]->y > max.y) max.y=vertices[i]->y;  
		if(vertices[i]->z > max.z) max.z=vertices[i]->z;  
	}

  	vox_min.y=ymin; //dimensions used to define CIVET vertices
    vox_min.z=zmin;
  	vox_max.y=ymin; //dimensions used to define CIVET vertices. start at the minima
    vox_max.z=zmin; //to cut down on time
  

    if( step[2] > 0  ){
    //Assuming that the step size if positive, while the minimum voxel value is less than the 
    //minimum vertex value, increase the voxel value. This way, the min voxel value will approach 
    //the min vertex value from the left.
    //If step size for x is positive:
    //x globel min
    //      |-------|-------|-------|-------|-------|-------|-------|-------|-------|
    //  vox_min.x--->------->------->------->------>o  |
    //                                                 min vertex value
    //
    
	    vox_min.x=xmin;  
	    vox_max.x=xmin;  
        while(vox_min.x /*voxel*/ < min.x /*vertex*/) vox_min.x += step[2];
    	while(vox_max.x < max.x) vox_max.x += step[2];
        vox_min.x -= step[2];//reduce the voxel minima so that they are smaller than the vertex minima

    }
    else{
	    vox_min.x=xmax;  
	    vox_max.x=xmax;  
        while(vox_min.x > min.x) vox_min.x += step[2];
    	while(vox_max.x > max.x) vox_max.x += step[2];
        vox_max.x -= step[2];
    }
    //printf("%f\n", vox_min.x);
    if(  step[1] > 0  ){
	    vox_min.y=ymin;  
	    vox_max.y=ymin;  
        while(vox_min.y < min.y) vox_min.y += step[1];
    	while(vox_max.y < max.y) vox_max.y += step[1];
        vox_min.y -= step[1];
    }
    else{
	    vox_min.y=ymax;  
	    vox_max.y=ymax;  
        while(vox_min.y > min.y) vox_min.y += step[1];
    	while(vox_max.y > max.y) vox_max.y += step[1];
        vox_max.y -= step[1];
    }

    if(  step[0] > 0  ){
	    vox_min.z=zmin;  
	    vox_max.z=zmin;  
        while(vox_min.z < min.z) vox_min.z += step[0];
    	while(vox_max.z < max.z) vox_max.z += step[0];
        vox_min.z -= step[0];
    }
    else{
	    vox_min.z=zmax;  
	    vox_max.z=zmax;  
        while(vox_min.z > min.z) vox_min.z += step[0];
    	while(vox_max.z > max.z) vox_max.z += step[0];
        vox_max.z -= step[0];
    }

    //if step is negative, the maximum value will be too low because we will hav added an extra step before exiting
    //the while loop. So we havce to subtract this exra step from the maximum value if the step is negative.

	/*Store the local max and min in maxima and minima*/
    minima[2]=vox_min.x; maxima[2]=vox_max.x;
	minima[1]=vox_min.y; maxima[1]=vox_max.y;
	minima[0]=vox_min.z; maxima[0]=vox_max.z;
	voxel_vertexPairs(vertices, maxima, minima, step, array,  xmin, ymin, zmin, min, max,  count);
    free(vertices);
    //if (counter > 100) break;
}
//for(z=0; z<native_count[0]; z++) for(y=0; y<native_count[1]; y++) for(x=0; x<native_count[2]; x++){
//    int index=   z*native_count[1]*native_count[2] + y*native_count[2] + x; 
//    if( array[index] >= 1 )  printf("%d %d %d\t%d\n", z, y, x, array[index]); 
//}

printf("Writing Array to MINC file %s.\n", outputfilename);

if( writeVolume( outputfilename, array, start,  native_step,  native_count, MI_TYPE_FLOAT ) != 0) exit(1) ;


return(0);
}
