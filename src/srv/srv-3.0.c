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

int readData(char* outerMeshfilename, char* innerMeshfilename, char* maskfilename, int* nvertices, vector** inner_mesh,  vector** outer_mesh, int** mask  ){
char buffer1[100], buffer2[100], buffer3[100], buffer4[100], buffer5[100], buffer6[100]; //should probably improve buffer length so that it is not a fixed size
int i; 
int vertices;
FILE *outerMeshFile;
FILE *innerMeshFile;
FILE *maskfile;
char *token;
char dlm[]=" "; 
int xmin=0, ymin=0, zmin=0;
int nodeCounter=0;

outerMeshFile=fopen(outerMeshfilename, "rt");
innerMeshFile=fopen(innerMeshfilename, "rt"); 
maskfile=fopen(maskfilename, "rt");

fgets(buffer1, sizeof(buffer1), outerMeshFile) ;
fgets(buffer2, sizeof(buffer2), innerMeshFile);
fgets(buffer3, sizeof(buffer3), maskfile);

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
*mask=calloc(*nvertices, sizeof(*mask));

for(i=0; i< *nvertices; i++){
	fgets(buffer1, sizeof(buffer1), outerMeshFile);
	fgets(buffer2, sizeof(buffer2), innerMeshFile);
    fgets(buffer3, sizeof(buffer3), maskfile);

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
	(*mask)[i]=atoi(buffer3);
}
 
fclose(maskfile);
fclose(outerMeshFile);
fclose(innerMeshFile);
return(0);
}

int voxel_vertexPairs(vector** vertices , double maxima[3], double minima[3], double step[3], float *voxelArray,  double global_xmin, double global_ymin, double global_zmin, vector vtx_min, vector vtx_max, misize_t* count ){
    float xmin=minima[2], ymin=minima[1], zmin=minima[0];
    float xmax=maxima[2], ymax=maxima[1], zmax=maxima[0];
    float xstep = step[2];
    float ystep = step[1];
    float zstep = step[0];
    float x=xmin, y=ymin, z=zmin;
    int u, v, w;
    int i=0, j, invoxels=0;
    vector testVoxel;
    vector** inner_mesh, **outer_mesh;
    float min; 
    int index;
    int signal=0;

    int zi_max= fabs((zmax-zmin)/zstep)+1;
    int yi_max= fabs((ymax-ymin)/ystep)+1;
    int xi_max= fabs((xmax-xmin)/xstep)+1;
    //printf("New\n%f\t%f\t%f, %d\n", xmin, xstep, xmax, xi_max);
    //printf("%f\t%f\t%f, %d\n", ymin, ystep, ymax, yi_max);
    //printf("%f\t%f\t%f, %d\n", zmin, zstep, zmax, zi_max);


    for(int zi=0; zi <= zi_max ; zi++ ){
        y=ymin;
	    for(int yi=0; yi <= yi_max; yi++){
            x=xmin;
            for(int xi=0; xi <= xi_max; xi++ ){
    //for(z=zmin; z <= zmax ; z){
	    //for(y=ymin; y <=  ymax ; y ){
            //for(x=xmin; x <=  xmax ; x ){
                //printf("2 . %f %f %f\n", z, y, x);
                //for(z=zmin; z <= zmax /*+ zstep*/; z += zstep){
	            //for(y=ymin; y <=  ymax/* + ystep*/ ; y += ystep){
                //for(x=xmin; x <=  xmax /*+ xstep*/; x += xstep){
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
                x+=xstep;
		    }
            y+=ystep;
        }
        z+=zstep;
    }

    //printf("%d %d\n", i, zi_max);
    //if(invoxels == 0 ){ 
        //printf("Number of Voxels in Mesh: %d\n", invoxels);
    //}   
    return(invoxels);
}

int main(int argc, char** argv){         
    if(argc != 6 && argc != 7 ) {
        printf("Usage:<Optional mask> <MRI file> <Outer surface vertices> <Inner surface vertices> <Scaling Factor>  <Output MINC file> \n");
        exit(1);
    }        
    int argi=1, i=0, x, y, z, w;          
    int nvertices=0;
    int *n_ngh_inner, **ngh_inner;
    int *n_ngh_mid, **ngh_mid;
    int *n_ngh_outer, **ngh_outer;
	int* mask;
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
	char* maskfilename;
	if(argc==7) maskfilename=argv[argi++];
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
    printf("Counts: %d %d %d\n", native_count[0], native_count[1], native_count[2]); 
    printf("Starts: %f %f %f\n", zmin, ymin, xmin); 
    printf("Step: %f %f %f\n", native_step[0], native_step[1], native_step[2]);
    float* array=calloc(nvox, sizeof(*array));
    int nvoxels_in_mesh=0;
	readData(outer_surface_file , inner_surface_file , maskfilename,   &nvertices, &inner_mesh,   &outer_mesh, &mask  );

    interpolate_sphere( outer_surface_file , &n_ngh_outer, &ngh_outer);
    interpolate_sphere( inner_surface_file , &n_ngh_inner, &ngh_inner);
	//interpolate_sphere( argv[4], &n_ngh_mid, &ngh_mid);

	linkVerticesToNeighbours( inner_mesh, n_ngh_inner, ngh_inner, nvertices);
	linkVerticesToNeighbours( outer_mesh, n_ngh_outer, ngh_outer, nvertices);
	//linkVerticesToNeighbours( mid_mesh, n_ngh_mid, ngh_mid, nvertices);

for(counter=0; counter < nvertices; counter++){

  	if(maskfilename != NULL) {
		if( mask[counter] == 0){ 
			continue;
	  	}
	}
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

    vox_min.z = nearest_voxel_world(min.z, step[0], start[0], count[0], -1);
    vox_min.y = nearest_voxel_world(min.y, step[1], start[1], count[1], -1);
    vox_min.x = nearest_voxel_world(min.x, step[2], start[2], count[2], -1);
    vox_max.z = nearest_voxel_world(max.z, step[0], start[0], count[0], 1);
    vox_max.y = nearest_voxel_world(max.y, step[1], start[1], count[1], 1);
    vox_max.x = nearest_voxel_world(max.x, step[2], start[2], count[2], 1);
    //printf("\n%f %f %f, %f %f %f\n", vox_min.z, vox_min.y, vox_min.x, vox_max.z, vox_max.y, vox_max.x ); 
    //printf("%f %f %f, %f %f %f\n", vox_min.z, vox_min.y, vox_min.x, vox_max.z, vox_max.y, vox_max.x ); 
    //printf("y %f %f %f %d --> %f\n", max.y, step[1], start[1], count[1], vox_max.y );

	/*Store the local max and min in maxima and minima*/
    if(step[2] > 0 ) { 
        minima[2]=vox_min.x; 
        maxima[2]=vox_max.x;
    } else {
        minima[2]=vox_max.x; 
        maxima[2]=vox_min.x;
    } 
    if(step[1] > 0 ) {
	    minima[1]=vox_min.y; 
        maxima[1]=vox_max.y;
	} else{
	    minima[1]=vox_max.y; 
        maxima[1]=vox_min.y;
    }
    if(step[0] > 0) {
        minima[0]=vox_min.z; 
        maxima[0]=vox_max.z;
    } else {
        minima[0]=vox_max.z; 
        maxima[0]=vox_min.z;
    }
    //printf("Vertex\n");
    //printf("z %f  %f < %f\t %f < %f\n", step[0], vox_min.z, min.z , max.z, vox_max.z);
    //printf("y %f  %f < %f\t %f < %f\n", step[1], vox_min.y, min.y , max.y, vox_max.y);
    //printf("x %f  %f < %f\t %f < %f\n", step[2], vox_min.x, min.x , max.x, vox_max.x);
    nvoxels_in_mesh += voxel_vertexPairs(vertices, maxima, minima, step, array,  xmin, ymin, zmin, min, max,  count);
    free(vertices);
    //if (counter > 10000) break;
}
printf("Total of %d voxels identified in mesh\n", nvoxels_in_mesh);
//for(z=0; z<native_count[0]; z++) for(y=0; y<native_count[1]; y++) for(x=0; x<native_count[2]; x++){
//    int index=   z*native_count[1]*native_count[2] + y*native_count[2] + x; 
//    if( array[index] >= 1 )  printf("%d %d %d\t%d\n", z, y, x, array[index]); 
//}

printf("Writing Array to MINC file %s.\n", outputfilename);
printf("Out %f %f %f\n", start[0], start[1], start[2]);
printf("Out %f %f %f\n", step[0], step[1], step[2]);
printf("Out %d %d %d\n", count[0], count[1], count[2]);
if( writeVolume( outputfilename, array, start,  step,  count, MI_TYPE_FLOAT ) != 0) exit(1) ;


return(0);
}
