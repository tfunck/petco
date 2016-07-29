#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <minc2.h>
#include <hdf5.h>
#include "volume_io.h"
#include "nrutil.h"
#include "fft.h"


int getdimensions(char *image_file, mihandle_t *image, int num_dimensions, misize_t sizes[3], midimhandle_t *dimensions, double step[3], double starts[3]  )
{
 int counter;
 int i, j, num_offset=0;
 float temp_step[4];

  /* check for error on opening */
  if ( miopen_volume(image_file, MI2_OPEN_READ, image) != MI_NOERROR ) 
  {
    fprintf(stderr, "Error opening input file %s.\n", image_file);
    return(1);
  }

  miget_volume_dimensions(*image, MI_DIMCLASS_ANY, MI_DIMATTR_ALL, MI_DIMORDER_FILE, num_dimensions, dimensions);
     
		if (miget_dimension_sizes(dimensions, num_dimensions, sizes)	== 0 )
		{
		//printf("File %s has dimensions", image_file);
		//	for(counter=0; counter < num_dimensions; counter++) printf(" %u", sizes[counter]);
		//printf(".\n");
		}
		else 
		{
		printf("Error: Could not get dimensions for %s\n", image_file); 
		return 1;
		}	


miget_dimension_starts(dimensions, MI_DIMORDER_FILE, num_dimensions, starts);
miget_dimension_separations(dimensions, MI_DIMORDER_FILE, num_dimensions, step);

return 0;
}

/*void blurImage(unsigned long *sizes, double* step, float*** voxeldata, float fwhm)
{
float*** filterdata;
int  pad_nslices, pad_nrows, pad_ncols;
    float steps[3]={step[0], step[1], step[2]}; 
     pad_nslices = nextPowerOf2(sizes[0]);
     pad_nrows   = nextPowerOf2(sizes[1]);
     pad_ncols   = nextPowerOf2(sizes[2]);
     printf("The padded image dimensions are (%d, %d, %d)\n", pad_nslices, pad_nrows, pad_ncols);
     filterdata=NULL;//genGaussianKernel(step, pad_nslices, pad_nrows, pad_ncols, fwhm);
     applyFFT( voxeldata, filterdata, sizes[0], sizes[1], sizes[2] ); //Perform the 3d filtering to this hyperslab 
   
     //probably need to free filter data? or move up its data allocation...      
}*/






