#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <assert.h>
#include <string.h>
#include <minc2.h>
#include <hdf5.h>
#include "mtga.h"
#include <volume_io.h>

float* get_arterial_TAC(int nframes,float* times, char* arterial_filename, parameters* userOptions  ){
    float* referenceTAC=malloc(nframes * sizeof(*referenceTAC) );
    int npoints=0;
    float** arterialTAC=importArterialTAC( arterial_filename, &npoints, userOptions); //read the arterial input from text file
        
    //interpolate arterial values to mid-frame time points  
    //printf("npoints: %d\n", npoints);  
    for(int i=0; i<nframes; i++){
        //printf("times[%d]=%f\n", i, times[i]);
        for(int j=0; j<npoints; j++){
            //printf("arterial time: %f\n", arterialTAC[j][0]);
            if(  times[i] <= arterialTAC[j][0]) {
                int max=j;
                int min=j-1;
                //Bilinear interpolation of the arterial input
                float d=arterialTAC[max][0] - arterialTAC[min][0];
                float d0 = 1-(times[i] - arterialTAC[min][0]) / d;
                float d1= 1-(arterialTAC[max][0] - times[i]) / d;
                //printf("max: %d min: %d\n", max, min);
                referenceTAC[i] = d0*arterialTAC[min][1] + d1 * arterialTAC[max][1];
                break;
            }
        }
        if( times[i] > arterialTAC[npoints-1][0]) referenceTAC[i]=arterialTAC[npoints-1][1];
    } //exit(0);
    free(arterialTAC);
    return(referenceTAC);
}


float** importArterialTAC (char* arterial_filename, int* points, parameters* userOptions){
    char buffer[100];
    char *dlm="	";
    int counter=0;
    float value, time;
    FILE* arterial_file=fopen(arterial_filename, "rt");
    float** arterialTAC=NULL; // 2D array for time [0] and radioactivity [1]
    //printf("Time\tValue\n");
    while(fgets(buffer, sizeof(buffer), arterial_file) != NULL){
        counter++;
        arterialTAC=realloc(arterialTAC, counter * sizeof(*arterialTAC));
        arterialTAC[counter-1]=malloc(sizeof(**arterialTAC)*2);

        arterialTAC[counter-1][0]=atof(strtok(buffer, dlm));
        if( userOptions->arterial_convert_to_seconds ==TRUE)  arterialTAC[counter-1][0] *= 60;

        arterialTAC[counter-1][1]=atof(strtok(NULL, dlm));

        printf("%s\t:\t%f\t%f\n", buffer, arterialTAC[counter-1][0], arterialTAC[counter-1][1]);
    }
    *points=counter;
    printf("Imported arterial TAC file.\n");

    fclose(arterial_file);
    return(arterialTAC);
}


void getTimeWidths(int nframes, int startFrame, float* timeWidths, char* inputfile, parameters userOptions)
{
char command[200];
FILE *input_f=fopen(inputfile, "rt");
char buffer[100], dlm[]=" ", *token;
int t;
int unit;

if(userOptions.min==TRUE) unit=60.0;
else unit=1.0;

for(t=0; t<nframes; t++){
            fgets(buffer, sizeof(buffer), input_f);
            timeWidths[t]=atof(buffer);
    }

//for(t=0; t<nframes; t++) printf("Frame: %d, Width: %f\n", t, timeWidths[t]);

fclose(input_f);
}


int getdimensions(char *image_file, mihandle_t *image, int num_dimensions, misize_t *sizes, midimhandle_t *dimensions, double* step, double* separations, double* starts) {
int counter;
int i, j, num_offset=0;
float temp_step[4];

  /* check for error on opening */
  if ( miopen_volume(image_file, MI2_OPEN_READ, image) != MI_NOERROR ) {
    fprintf(stderr, "Error opening input file %s.\n", image_file);
    return(1);
  }

  miget_volume_dimensions(*image, MI_DIMCLASS_ANY, MI_DIMATTR_ALL, MI_DIMORDER_FILE, num_dimensions, dimensions);
   
		if (miget_dimension_sizes(dimensions, num_dimensions, sizes)	== 0 ){
		    printf("File %s has dimensions", image_file);
			for(counter=0; counter < num_dimensions; counter++) printf(" %u", sizes[counter]);
		    printf(".\n");
		}
		else {
		    printf("Error: Could not get dimensions for %s\n", image_file); 
		    return 1;
		}	
 
    if(num_dimensions==4)/*Get offset values only for 4d pet image, not 3d masks*/
    {
    miget_dimension_starts(dimensions, MI_DIMORDER_FILE, 4, starts);
    miget_dimension_separations(dimensions, MI_DIMORDER_FILE, 4, separations);

    /*Now we change the order of the step array so that it can be passed to the gaussian kernel and FFT functions in the    same order as the 3d slice volume is allocated. This is a little confusing but it makes sense when you consider that the image_array is allocated according to x, y, t and not x, y, z*/
    step[0]=(float)separations[1]; /*set t-step to first position in step array*/
    step[1]=(float)separations[2]; /*set y-step to first position in step array*/
    step[2]=(float)separations[3]; /*set x-step to first position in step array*/
    printf("Average separations: t=%f, z=%f y=%f x=%f\nStarts: %f %f %f %f\n", separations[0], separations[1], separations[2], separations[3], starts[0], starts[1], starts[2], starts[3] );
    }


return 0;
}

void find_reference_values(char *petimage_file, misize_t *sizes,  float **reference_array, float* referencemask){
misize_t location[4], count[4], start[4];
int x, y, z, t, result, counter;
int tmax, xmax, ymax, zmax;
int maskindex, slabindex;
int *reference_count;
float voxelvalue;
float* slab;

tmax=(int) sizes[0];/*save sizes in individual variables to make it easier to keep track of*/
zmax=(int) sizes[1];
ymax=(int) sizes[2];
xmax=(int) sizes[3];

start[0]=start[1]=start[2]=start[3]=(misize_t) 0; /*update start z (slice)*/
count[0]=(misize_t) 1; /* slab end time */
count[1]=(misize_t) zmax; /* slab end z*/
count[2]=(misize_t) ymax; /* slab end y */
count[3]=(misize_t) xmax; /* slab end x */

slab=malloc( ymax * xmax * zmax * sizeof(float) );
if(slab == NULL) {
	printf("Error allocating memory for slab\n");
	exit(0);
}
*reference_array=malloc( tmax * sizeof(float) ); /*  */
reference_count=malloc( tmax * sizeof(int));/*array for keeping track of number of reference voxels in each time frame*/

mihandle_t image2;
  if ( miopen_volume( petimage_file, MI2_OPEN_READ, &image2) != MI_NOERROR ) {
    fprintf(stderr, "Error opening input file %s.\n", "tst.mnc");
    exit(1);
  }

for(t=0; t<tmax; t++) { reference_count[t]=0; (*reference_array)[t]=0;}
	for(t=0; t < tmax; t++) /*find ref values by finding petimage values that are in ref region*/
	{
		start[0]=t;
		result=miget_real_value_hyperslab( image2, MI_TYPE_FLOAT, start, count, slab );
		if ( result == MI_ERROR ) /* Get hyperslab according to start and count */
		{ 
			fprintf(stderr, "Error getting hyperslab for 4D slab\n"); 
			exit(1); 
		}	

		for(z=0; z < zmax; z++) 
		{
			for(y=0; y < ymax; y++)
			{
				for(x=0; x < xmax; x++)
				{
					maskindex=( z * ( ymax * xmax) + (y * xmax) + x );
					//slabindex=( t * (ymax * xmax * zmax) + z *( ymax * xmax) + (y * xmax) + x );
				  	slabindex=( z *( ymax * xmax) + (y * xmax) + x );

                    if ( referencemask[maskindex]  == 1){
					     (*reference_array)[t] += slab[ slabindex];
					     reference_count[t] += 1;     
					}
				}
			}
		}
	}

for (t=0; t<tmax; t++) { (*reference_array)[t] /= (float) reference_count[t]; printf("REF[%d]=%f\t%d\n", t , (*reference_array)[t], reference_count[t]); }
free(reference_count);
free(slab);
miclose_volume(image2);

}




double loadpetimage( int z, misize_t *sizes, misize_t *start, misize_t *count, char* petimage_file, float *slab, float *brainmask, float  *referencemask, parameters userOptions)
{
int x, y, t, counter, i=0.0;
int slabindex, maskindex;
int tmax, zmax, ymax, xmax;
misize_t *start_ptr; /*Points to first or second element of start depending on if input is 4d or 3d*/
misize_t *count_ptr; /*Points to first or second element of count depending on if input is 4d or 3d*/
mihandle_t image;

miopen_volume(petimage_file, MI2_OPEN_READ, &image);

tmax=sizes[0]; zmax=sizes[1]; ymax=sizes[2]; xmax=sizes[3];
	printf("loading: %d %d %d %d\n%d %d %d %d\n", start[0], start[1], start[2], start[3], count[0], count[1], count[2], count[3]);
if ( userOptions.input3d == TRUE ) start_ptr=&(start[1]); 
else start_ptr=start;

if( (start_ptr == &(start[1]) || start_ptr == start ) == FALSE ) {fprintf(stderr, "Error setting start pointer\n"); return 0;}

if ( userOptions.input3d == TRUE ) count_ptr=&(count[1]);
else count_ptr=count;

if( (count_ptr == &(count[1]) || count_ptr == count ) == FALSE ) {fprintf(stderr, "Error setting count pointer\n"); return 0;}

if ( miget_real_value_hyperslab(image, MI_TYPE_FLOAT, start_ptr, count_ptr, slab) == MI_ERROR )
	 /* Read in smoothed PET data according to start and count for current slice , i.e. z	*/
 	{ 
		fprintf(stderr, "Error getting hyperslab for slice %d\n", z); 
		return(1); 
	}
miclose_volume(image);

return 0;
}


void loadmask_hyperslab(char* mask_file, float** mask)
{
misize_t start[3], count[3];
misize_t sizes[3];
int index, i, t;
int xmax, ymax, zmax, z, y, x;
midimhandle_t  dimensions[3];
mihandle_t image;

getdimensions(mask_file, &image, 3, sizes, dimensions, NULL, NULL, NULL);
zmax=sizes[0];ymax=sizes[1];  xmax=sizes[2]; 

*mask=malloc( zmax * ymax * xmax * sizeof(float));
if(mask == NULL) {fprintf(stderr, "Error allocating memory for mask %s", mask_file ); exit(1); }


start[0]=start[1]=start[2]=(misize_t) 0; 

count[0]=(misize_t) zmax; /* slab end z */
count[1]=(misize_t) ymax; /* slab end y */
count[2]=(misize_t) xmax; /* slab end x */
	if (miget_real_value_hyperslab(image, MI_TYPE_FLOAT, start, count, *mask) == MI_ERROR)
	{ 
		fprintf(stderr, "Error getting hyperslab for mask %s", mask_file); 
		exit(1); 
	}
//for(i=0; i<xmax*ymax*zmax; i++) if( (*mask)[i] == 1.0) printf("%s %f\n", mask_file, (*mask)[i]  );
miclose_volume(image);

}

void createImage( char* newfilename, int num_dim, misize_t *sizes, double *separations, double *starts)
{
    mihandle_t hvol;
    midimhandle_t hdim[num_dim];
    long count[num_dim];
    double altSeparations[num_dim];
	double  altStarts[num_dim];
	int CZ, CY, CX, CT, r, NDIMS = num_dim;

    if(num_dim==4)
	{
	 CT=sizes[0]; CZ=sizes[1]; CY=sizes[2]; CX=sizes[3];
    r = micreate_dimension("time", MI_DIMCLASS_TIME,  MI_DIMATTR_REGULARLY_SAMPLED, CT, &hdim[0]);
	r = micreate_dimension("zspace", MI_DIMCLASS_SPATIAL, MI_DIMATTR_REGULARLY_SAMPLED, CZ, &hdim[1]);  
    r = micreate_dimension("yspace", MI_DIMCLASS_SPATIAL, MI_DIMATTR_REGULARLY_SAMPLED, CY, &hdim[2]);
    r = micreate_dimension("xspace", MI_DIMCLASS_SPATIAL, MI_DIMATTR_REGULARLY_SAMPLED, CX, &hdim[3]);
    if (r == MI_NOERROR) printf("Created Dimensions.\n");
    else printf("Could not create dimensions.\n");

    r = miset_dimension_starts(hdim, NDIMS, starts);
    
    if (r == MI_NOERROR) printf("Created starts.\n");
    else printf("Could not create starts.\n");

    r = miset_dimension_separations(hdim, NDIMS, separations); 
   	if (r == MI_NOERROR) printf("Created separations.\n");
	else printf("Could not create seperations.\n");

    r = micreate_volume(newfilename, NDIMS, hdim, MI_TYPE_FLOAT, MI_CLASS_REAL, NULL, &hvol);
    
    r = miset_slice_scaling_flag(hvol, FALSE); /* indicate that we will be using slice scaling */
    	if (r == MI_NOERROR) printf("Set Scaling.\n");
	else printf("Could not set scaling.\n");

    r = micreate_volume_image(hvol);
   	 if (r == MI_NOERROR) printf("Created output volume.\n");
	 	else printf("Error creating output volume.\n");        
	printf("\nDimensions for %s:\n\t%d %d %d %d\n\t%f %f %f %f\n\t%f %f %f %f\n\n", newfilename, sizes[0], sizes[1], sizes[2], sizes[3], separations[0], separations[1], separations[2], separations[3],  starts[0], starts[1], starts[2], starts[3] );                                                    
	}
	else if(num_dim==3)
	{
	altSeparations[0]=separations[1];
	altSeparations[1]=separations[2];
	altSeparations[2]=separations[3];
	altStarts[0]=starts[1];
	altStarts[1]=starts[2];
	altStarts[2]=starts[3];

	 CZ= abs(sizes[1]); CY=abs(sizes[2]); CX=abs(sizes[3]);
	r = micreate_dimension("zspace", MI_DIMCLASS_SPATIAL, MI_DIMATTR_REGULARLY_SAMPLED, CZ, &hdim[0]);  
    	r = micreate_dimension("yspace", MI_DIMCLASS_SPATIAL, MI_DIMATTR_REGULARLY_SAMPLED, CY, &hdim[1]);
    	r = micreate_dimension("xspace", MI_DIMCLASS_SPATIAL, MI_DIMATTR_REGULARLY_SAMPLED, CX, &hdim[2]);
		if (r == MI_NOERROR) printf("Created Dimensions.\n");
		else printf("Could not create dimensions.\n");
	
	
    	r = miset_dimension_starts(hdim, NDIMS, altStarts);
    		if (r == MI_NOERROR) printf("Created starts.\n");
		else printf("Could not create starts.\n");
	r = miset_dimension_separations(hdim, NDIMS, altSeparations); 
		if (r == MI_NOERROR) printf("Created separations.\n");
		else printf("Could not create seperations.\n");
    	r = micreate_volume(newfilename, NDIMS, hdim, MI_TYPE_FLOAT, MI_CLASS_REAL, NULL, &hvol);
	r = miset_slice_scaling_flag(hvol, FALSE); /* indicate that we will be using slice scaling */
		if (r == MI_NOERROR) printf("Set Scaling.\n");
		else printf("Could not set scaling.\n");
    	r = micreate_volume_image(hvol); 
    	if (r == MI_NOERROR) printf("Created output volume.\n");
		else printf("Error creating output volume.\n");

	printf("Dimensions for %s:\n\t %d %d %d\n\t %f %f %f\n\t %f %f %f\n\n", newfilename, sizes[1], sizes[2], sizes[3], separations[1], separations[2], separations[3], altStarts[0], altStarts[1], altStarts[2] );

	}
    miclose_volume(hvol);
}

