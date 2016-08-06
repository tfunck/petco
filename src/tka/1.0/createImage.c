#include <stdio.h>
#include <stdlib.h>`
#include <math.h>
#include <assert.h>
#include <string.h>
#include <minc2.h>
#include <hdf5.h>
#include "mtga.h"
#include <volume_io.h>


void main()
{
	mihandle_t hvol, *hvol_ptr;
	char* newfilename="testImage.mnc"; 
	int num_dim=3;
	unsigned int sizes[3]={10, 10, 10};
	double separations[3]={1, 1, 1};
	double starts[3]={0, 0, 0};
    midimhandle_t hdim[num_dim];
    long count[3]={10, 10, 10};
	
	 int CZ, CY, CX, CT, r, NDIMS = num_dim;
	hvol_ptr=&hvol;
	printf("Dimensions for %s:\n\t %d %d %d\n\t %f %f %f\n\t %f %f %f\n\n", newfilename, sizes[0], sizes[1], sizes[2],
													  separations[0], separations[1], separations[2],
													   starts[0], starts[1], starts[2] );

	 CZ= abs(sizes[0]); CY=abs(sizes[1]); CX=abs(sizes[2]);
	r = micreate_dimension("zspace", MI_DIMCLASS_SPATIAL, MI_DIMATTR_REGULARLY_SAMPLED, CZ, &hdim[0]);  
    r = micreate_dimension("yspace", MI_DIMCLASS_SPATIAL, MI_DIMATTR_REGULARLY_SAMPLED, CY, &hdim[1]);
    r = micreate_dimension("xspace", MI_DIMCLASS_SPATIAL, MI_DIMATTR_REGULARLY_SAMPLED, CX, &hdim[2]);
		if (r == MI_NOERROR) printf("Created Dimensions.\n");
		else printf("Could not create dimensions.\n");
	
	
    r = miset_dimension_starts(hdim, NDIMS, starts);
    		if (r == MI_NOERROR) printf("Created starts.\n");
		else printf("Could not create starts.\n");
	r = miset_dimension_separations(hdim, NDIMS, separations); 
			if (r == MI_NOERROR) printf("Created separations.\n");
		else printf("Could not create seperations.\n");
    r = micreate_volume(newfilename, NDIMS, hdim, MI_TYPE_FLOAT, MI_CLASS_REAL, NULL, hvol_ptr);
 
	r = miset_slice_scaling_flag(*hvol_ptr, FALSE); /* indicate that we will be using slice scaling */
		if (r == MI_NOERROR) printf("Set Scaling.\n");
		else printf("Could not set scaling.\n");
    r = micreate_volume_image(*hvol_ptr); 
    	if (r == MI_NOERROR) printf("Created output volume.\n");
		else printf("Error creating output volume.\n");

float matrix[CZ*CY*CX];
	for(r=0; r<CZ*CY*CX; r++) matrix[r]=r;

miset_real_value_hyperslab(*hvol_ptr, MI_TYPE_FLOAT, starts, count, matrix);

}


