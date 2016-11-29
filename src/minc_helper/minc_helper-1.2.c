#include "minc_helper.h"

int real2voxel(float location, float min, float step){
//simple function to convert real to voxel coordinates
    //step=1;
    //if ( fmod((location-min), step) == 0 ) printf("Problem! %f\n",(location-min) / step );
    int out=(int) round((location-min)/step);
    //printf("(%f - %f) / %f --> %d\n", location,min, step, out);
    return(out );
}

float voxel2real(int location, float min, float step){
//simple function to convert voxel to real coordinates
    return (float) min+location*step;
}


double* malloc_check(int dim1, int dim2, int dim3, int dim4, int dtype, int ndim){
   size_t n=0;
   if (ndim==4){
       n=dim1*dim2*dim3*dim4;
    }
    else{
	    n=(unsigned long) dim2*dim3*dim4;
	}
   	void *pointer=malloc( dtype*n);
   	if( pointer == NULL){
       printf("Error: Insufficient memory. %3.3f GB required\n", dtype*n*0.0000000009313226 );
       exit(0);
   }

   return pointer;
}

char* strip_path_ext(char* instr) {
    char *basestr=basename(instr);
    char *retstr;
    char *lastdot;
    if (basestr == NULL) 
        return NULL;
    if ((retstr = malloc (strlen (basestr) + 1)) == NULL) 
        return NULL;
    strcpy (retstr, basestr);
    lastdot = strrchr (retstr, '.');
    if (lastdot != NULL)
        *lastdot = '\0';
    return retstr;
}

int allocate_atlas(data** ptr, int* n){
    *n += 1;
    data* temp=realloc(*ptr, sizeof(**ptr) * *n);
    
    if(temp == NULL) pexit("Error:", "could not allocate memory for realloc in <allocate_data>", 1);
    else *ptr=temp;
    
    (*ptr)[*n-1].label=*n-1; 
    (*ptr)[*n-1].ngroups=0;
    (*ptr)[*n-1].tmax=(*ptr)[*n-1].zmax=(*ptr)[*n-1].ymax=(*ptr)[*n-1].xmax=0;
    (*ptr)[*n-1].n=0;
    //Initialize pointers in vol to NULL.
    (*ptr)[*n-1].groups=NULL;
    (*ptr)[*n-1].filename=NULL;
    //(*ptr)[*n-1].like_filename=NULL;
    (*ptr)[*n-1].count=NULL;
    (*ptr)[*n-1].start=NULL;
    (*ptr)[*n-1].step=NULL;
    (*ptr)[*n-1].data=NULL;
    (*ptr)[*n-1].wcount=NULL;
    (*ptr)[*n-1].wstarts=NULL;
    (*ptr)[*n-1].values_filename=NULL;

    return(0);
}

int allocate_vol(data** ptr, int* n){
    *n += 1;
    void* temp=realloc(*ptr, sizeof(**ptr) * *n);
    if(temp == NULL) pexit("Error:", "could not allocate memory for realloc in <allocate_data>", 1);
    else *ptr=temp;
    (*ptr)[*n-1].filename=NULL;
    (*ptr)[*n-1].values_filename=NULL;
    (*ptr)[*n-1].label=*n-1; 
    (*ptr)[*n-1].ngroups=0;
    (*ptr)[*n-1].tmax=(*ptr)[*n-1].zmax=(*ptr)[*n-1].ymax=(*ptr)[*n-1].xmax=0;
    (*ptr)[*n-1].n=0;
    (*ptr)[*n-1].ngroups=0;
    //Initialize pointers in vol to NULL.
    (*ptr)[*n-1].groups=NULL;
    (*ptr)[*n-1].count=NULL;
    (*ptr)[*n-1].start=NULL;
    (*ptr)[*n-1].step=NULL;
    (*ptr)[*n-1].data=NULL;
    (*ptr)[*n-1].wcount=NULL;
    (*ptr)[*n-1].wstarts=NULL;

    return(0);
}


double nearest_neighbour(double* array, int  zmax, int sz, double oz, int  ymax,int  sy, double  oy, int xmax,int  sx, double  ox,int  xi,int yi,int zi){
    int zlo, zhi, ylo, yhi, xlo, xhi;
    double  x1, x0, y1, y0, z0, z1;
    int xn, yn, zn; 
    int i=0, j=0;
	/*
       for (i=0; i* step[2]+min[2] <=  xi ; i++);
       x0=min[2] + (i-1)*step[2] ; 
	   xlo=i-1;
       x1=min[2] + i*step[2];
       xhi=i;			
            
       if( pow(x0 - xi, 2) - pow(x1 - xi, 2) <= 0.0 ) xn=xlo;
       else xn=xhi;
            
       for (i=0; i* step[1]+min[1] <=  yi ; i++);
       y0=min[1] + (i-1)*step[1];
	   y1=min[1] + i*step[1];
       ylo= i-1;
       yhi=i;

       if( (pow(y0 - yi, 2) - pow(y1 - yi, 2)) <= 0.0) yn=ylo;
       else yn=yhi;
       
       for (i=0; i* step[0]+min[0] <=  zi ; i++);
	   z0=min[0] + (i-1)*step[0];
	   z1=min[0]+ i*step[0];
       zlo=i-1;
       zhi=i;

       if( pow(z0 - zi, 2) - pow(z1 - zi, 2) >= 0) zn=zlo;
       else zn=zhi;

       //if( isnan(array[zn * image->xmax*image->ymax + yn*image->xmax + xn])){
       //    printf("%f = array[%d][%d][%d]\n",array[zn * image->xmax*image->ymax + yn*image->xmax + xn], zn, yn ,xn ); 
       //    exit(0);
       // }
       return array[zn * image->xmax*image->ymax + yn*image->xmax + xn];
	   */
	return(0);

}

int extrema(double* array, int tmax, int zmax, int ymax, int xmax, double* min, double*max){
    *min=*max=array[0];
    for(int t=0; t < tmax; t++)
        for(int z=0; z < zmax; z++)
            for(int y=0; y < ymax; y++)
                for(int x=0; x < xmax; x++){
                    int i=t*zmax*ymax*xmax + z*ymax*xmax+y*xmax+x;
                    if( array[i] > *max  ) *max=array[i]; 
                    else if( array[i] < *min) *min=array[i];   
                }

 return(0);
}

int sign(double x){
    if (x > 0) return 1;
    if (x < 0) return -1;
    return 0;
}

int nearest_voxel(double vertex, double step, double start, int nvox, int direction /*-1=smallest, 1=largest*/){
    int output;
    if( sign(step) == 1){ /*positive step*/
        if(direction == -1) output=ceil((vertex-step-start)/step);      
        else output=floor((vertex+step-start)/step); 
    } else{  
        if(direction == -1) output=floor((vertex+step-start)/step);         
        else output=ceil((vertex-step-start)/step);
    }
    return( output);
}

double nearest_voxel_world(double vertex, double step, double start, int nvox, int direction /*-1=smallest, 1=largest*/){
    int output;
    if( sign(step) == 1){ /*positive step*/
        if(direction == -1) output=ceil((vertex-step-start)/step);      
        else output=floor((vertex+step-start)/step); 
    } else{  
        if(direction == -1) output=floor((vertex+step-start)/step);         
        else output=ceil((vertex-step-start)/step);
    }
    return( (double) output*step+start);
}


double trilinear(double* array, int zmax, double sz, double oz, int ymax, double sy, double oy,int xmax, double sx, double ox, double xi, double yi, double zi)
{	int i=0, j=0;
    int xymax=xmax*ymax;
    int zlo = nearest_voxel(zi, sz, oz, zmax, -1); 
    int ylo = nearest_voxel(yi, sy, oy, ymax, -1); 
    int xlo = nearest_voxel(xi, sx, ox, xmax, -1); 
	int zhi=zlo+1;
	int yhi=ylo+1;
	int xhi=xlo+1;
    double  z0=oz+sz*zlo;
	double  y0=oy+sy*ylo;
    double  x0=ox+sx*xlo;
    double  z1=oz+sz*zhi;
	double  y1=oy+sy*yhi;
    double  x1=ox+sx*xhi;
	double c00, c01, c10, c11, c0, c1, c;
    double xd=(xi-x0)/(x1-x0);
    double yd=(yi-y0)/(y1-y0);
    double zd=(zi-z0)/(z1-z0);
	/*printf("origin: %f %f %f\n", oz, oy, ox);
	printf("point: %f %f %f\n", zi, yi, xi);
	printf("z: %d %f %d\n", zlo, zi, zhi);
	printf("y: %d %f %d\n", ylo, yi, yhi);
	printf("x: %d %f %d\n", xlo, xi, xhi);
	printf("%d %d %d\n", zmax, ymax, xmax);
	/**/
	  //find nearest voxels for real location x, y, z

	/*  LINEAR INTERPOLATION    */
	c00=array[zlo*xymax + ylo*xmax + xlo]*(1-xd) + array[zlo*xymax+  ylo*xmax + xhi]*xd; 
	c10=array[zlo*xymax + yhi*xmax + xlo]*(1-xd) + array[zlo*xymax + yhi*xmax + xhi]*xd;
	c01=array[zhi*xymax + ylo*xmax + xlo]*(1-xd) + array[zhi*xymax + ylo*xmax + xhi]*xd;
	c11=array[zhi*xymax + yhi*xmax + xlo]*(1-xd) + array[zhi*xymax + yhi*xmax + xhi]*xd;
	c0=c00*(1-yd) + c10*yd;
	c1=c01*(1-yd)+c11*yd;
	c = c0*(1-zd)+c1*zd; // interpolation value is the blurred value derived from the test image
    return c;
}




int writeVolume(char* outputfile, void* array, double starts[3],  double separations[3],  misize_t  count[3], mitype_t dtype  ){
    misize_t start[3]={0,0,0};

    int r=writeVolumeSubset(outputfile, array, start, starts,  separations,count, dtype);
    return(r);
}


int writeVolumeSubset(char* outputfile, void* array, misize_t start[3], double starts[3],  double separations[3],  misize_t  count[3], mitype_t dtype  ) {
mihandle_t outputimage;
midimhandle_t hdim[3];
misize_t wstart[3]={0,0,0};
int counter;
int r, z, y, x, index;
double max=0.0, min=0.0, value;
int xmax=(int) count[2]; //this can cause problems if maxima are not whole numbers
int ymax=(int) count[1];
int zmax=(int) count[0];

//printf("%f %f %f\n", starts[0], starts[1], starts[2]);
//printf("%f %f %f\n", separations[0], separations[1], separations[2]);
//printf("%d %d %d\n", count[0], count[1], count[2]);

if(dtype == MI_TYPE_INT) min=max=(double) ((int*) array)[0];
else if(dtype == MI_TYPE_USHORT)  min=max=(double) ((unsigned short*) array)[0];
else if(dtype == MI_TYPE_UINT) min=max=(double) ((unsigned int*) array)[0];
else if(dtype == MI_TYPE_UBYTE) min=max=(double) ((char*) array)[0];
else if(dtype == MI_TYPE_DOUBLE) min=max=((double*) array)[0];
else if(dtype == MI_TYPE_FLOAT){ min=max=(double) ((float*) array)[0];}
else {printf("Error: Could not recognize data type for writing MINC file.\n"); exit(1);}
#pragma omp parallel  
for (z=start[0]; z < start[0]+count[0]; z++) {
    for (y=start[1]; y < start[1]+count[1]; y++) {
        for (x=start[2]; x < start[2]+count[2]; x++) {
            index=z*ymax*xmax+y*xmax+x;
            if(dtype == MI_TYPE_INT) value=(double)  ((int*) array)[index];
            else if(dtype == MI_TYPE_USHORT) value= (double)  ((unsigned short*) array)[index];
            else if(dtype == MI_TYPE_UINT)  value=(double) ((unsigned int*) array)[index];
            else if(dtype == MI_TYPE_UBYTE) value=(double) ((char*) array)[index];
            else if(dtype == MI_TYPE_DOUBLE)value=((double*) array)[index];
            else if(dtype == MI_TYPE_FLOAT) value=(double) ((float*) array)[index];
            if (max < value) max= value;
            if (min > value) min= value;
            //if( ((double*) array)[x] > 0 ) printf("%f\t%f\t%f\n", max, value, ((float*) array)[x] );
            //if( value >= 1) printf("%f\n", value); 
        }
    }
}
//    printf("Min: %f  Max: %f\n", min, max);
    r = micreate_dimension("zspace", MI_DIMCLASS_SPATIAL, MI_DIMATTR_REGULARLY_SAMPLED, count[0], &hdim[0]);
    r = micreate_dimension("yspace", MI_DIMCLASS_SPATIAL, MI_DIMATTR_REGULARLY_SAMPLED, count[1], &hdim[1]);
    r = micreate_dimension("xspace", MI_DIMCLASS_SPATIAL, MI_DIMATTR_REGULARLY_SAMPLED, count[2], &hdim[2]);
    r = miset_dimension_starts(hdim, 3, starts);
    r = miset_dimension_separations(hdim, 3, separations);
    r = micreate_volume(outputfile, 3, hdim, dtype, MI_CLASS_REAL, NULL, &outputimage);
    r = miset_slice_scaling_flag(outputimage, FALSE ); /* indicate that we will be using slice scaling */
    r = micreate_volume_image(outputimage);
    if (r == MI_NOERROR); //printf("Created output volume.\n");
    else printf("Error creating output volume.\n");
    miset_volume_range ( outputimage, max,  min);
	int status=miset_real_value_hyperslab(outputimage, dtype, start, count, array);
	if ( status != MI_NOERROR ) fprintf(stderr, "\tError setting hyperslab\n");
    miclose_volume(outputimage);
return 0;
}

int writeVolume4d(char* outputfile, void* array, double starts[4],  double separations[4],  misize_t  count[4], mitype_t dtype  ) {
    mihandle_t outputimage;
    midimhandle_t hdim[4];
    misize_t wstart[4]={0,0,0,0};
    int counter;
    int r, z, y, x,t, index;
    double max=0.0, min=0.0, value;
    int tmax=(int) count[0];
    int zmax=(int) count[1];
    int ymax=(int) count[2];
    int xmax=(int) count[3]; //this can cause problems if maxima are not whole numbers

    //printf("%f %f %f\n", starts[0], starts[1], starts[2]);
    //printf("%f %f %f\n", separations[0], separations[1], separations[2]);
    //printf("%d %d %d\n", count[0], count[1], count[2]);

    if(dtype == MI_TYPE_INT) min=max=(double) ((int*) array)[0];
    else if(dtype == MI_TYPE_USHORT)  min=max=(double) ((unsigned short*) array)[0];
    else if(dtype == MI_TYPE_UINT) min=max=(double) ((unsigned int*) array)[0];
    else if(dtype == MI_TYPE_UBYTE) min=max=(double) ((char*) array)[0];
    else if(dtype == MI_TYPE_DOUBLE) min=max=((double*) array)[0];
    else if(dtype == MI_TYPE_FLOAT){ min=max=(double) ((float*) array)[0];}
    else {printf("Error: Could not recognize data type for writing MINC file.\n"); exit(1);}
    //#pragma omp parallel  
    for (t=0; t < count[0]; t++) {
        for (z=1; z < count[1]; z++) {
            for (y=2; y < count[2]; y++) {
                for (x=3; x < count[3]; x++) {
                    index=t*zmax*ymax*xmax+z*ymax*xmax+y*xmax+x;
                    if(dtype == MI_TYPE_INT) value=(double)  ((int*) array)[index];
                    else if(dtype == MI_TYPE_USHORT) value= (double)  ((unsigned short*) array)[index];
                    else if(dtype == MI_TYPE_UINT)  value=(double) ((unsigned int*) array)[index];
                    else if(dtype == MI_TYPE_UBYTE) value=(double) ((char*) array)[index];
                    else if(dtype == MI_TYPE_DOUBLE)value=((double*) array)[index];
                    else if(dtype == MI_TYPE_FLOAT) value=(double) ((float*) array)[index];
                    if (max < value) max= value;
                    if (min > value) min= value;
                    //if( ((double*) array)[x] > 0 ) printf("%f\t%f\t%f\n", max, value, ((float*) array)[x] );
                    //if( value >= 1) printf("%f\n", value); 
                }
            }
        }
    }
    //printf("Min: %f  Max: %f\n", min, max);
    r = micreate_dimension("time", MI_DIMCLASS_TIME,  MI_DIMATTR_REGULARLY_SAMPLED, count[0], &hdim[0]);
    r = micreate_dimension("zspace", MI_DIMCLASS_SPATIAL, MI_DIMATTR_REGULARLY_SAMPLED, count[1], &hdim[1]);
    r = micreate_dimension("yspace", MI_DIMCLASS_SPATIAL, MI_DIMATTR_REGULARLY_SAMPLED, count[2], &hdim[2]);
    r = micreate_dimension("xspace", MI_DIMCLASS_SPATIAL, MI_DIMATTR_REGULARLY_SAMPLED, count[3], &hdim[3]);
    r = miset_dimension_starts(hdim, 4, starts);
    r = miset_dimension_separations(hdim, 4, separations);
    r = micreate_volume(outputfile, 4, hdim, dtype, MI_CLASS_REAL, NULL, &outputimage);
    r = miset_slice_scaling_flag(outputimage, FALSE ); /* indicate that we will be using slice scaling */
    r = micreate_volume_image(outputimage);
    if (r == MI_NOERROR); //printf("Created output volume.\n");
    else printf("Error creating output volume.\n");
    miset_volume_range ( outputimage, max,  min);
	int status=miset_real_value_hyperslab(outputimage, dtype, wstart, count, array);
	if ( status != MI_NOERROR ) fprintf(stderr, "\tError setting hyperslab\n");
    miclose_volume(outputimage);
return 0;
}




void createVolume(char* newfilename, int num_dim, misize_t *sizes, double *separations, double *starts, mitype_t data_type){
    mihandle_t hvol;
    midimhandle_t hdim[num_dim];
    misize_t count[num_dim];
    double altSeparations[num_dim];
	double  altStarts[num_dim];
	int CZ, CY, CX, CT, r, NDIMS = num_dim;

    if(num_dim==4){
	    CT=sizes[0]; CZ=sizes[1]; CY=sizes[2]; CX=sizes[3];
        r = micreate_dimension("time", MI_DIMCLASS_TIME,  MI_DIMATTR_REGULARLY_SAMPLED, CT, &hdim[0]);
        r = micreate_dimension("zspace", MI_DIMCLASS_SPATIAL, MI_DIMATTR_REGULARLY_SAMPLED, CZ, &hdim[1]);  
        r = micreate_dimension("yspace", MI_DIMCLASS_SPATIAL, MI_DIMATTR_REGULARLY_SAMPLED, CY, &hdim[2]);
        r = micreate_dimension("xspace", MI_DIMCLASS_SPATIAL, MI_DIMATTR_REGULARLY_SAMPLED, CX, &hdim[3]);
        if (r == MI_NOERROR);// printf("Created Dimensions.\n");
        else printf("Could not create dimensions.\n");

        r = miset_dimension_starts(hdim, NDIMS, starts);
        
        if (r == MI_NOERROR);// printf("Created starts.\n");
        else printf("Could not create starts.\n");

        r = miset_dimension_separations(hdim, NDIMS, separations); 
        if (r == MI_NOERROR);// printf("Created separations.\n");
        else printf("Could not create seperations.\n");

        r = micreate_volume(newfilename, NDIMS, hdim, data_type, MI_CLASS_REAL, NULL, &hvol);
        
        r = miset_slice_scaling_flag(hvol, FALSE); /* indicate that we will be using slice scaling */
        if (r == MI_NOERROR);// printf("Set Scaling.\n");
        else printf("Could not set scaling.\n");

        r = micreate_volume_image(hvol);
         if (r == MI_NOERROR) ;// printf("Created output volume.\n");
            else printf("Error creating output volume.\n");        
        printf("\nDimensions for %s:\n\t%d %d %d %d\n\t%f %f %f %f\n\t%f %f %f %f\n\n", newfilename, sizes[0], sizes[1], sizes[2], sizes[3], separations[0], separations[1], separations[2], separations[3],  starts[0], starts[1], starts[2], starts[3] );                                                    
	}
	else if(num_dim==3){
        altSeparations[0]=separations[0];
        altSeparations[1]=separations[1];
        altSeparations[2]=separations[2];
        altStarts[0]=starts[0];
        altStarts[1]=starts[1];
        altStarts[2]=starts[2];

        CZ= abs(sizes[0]); CY=abs(sizes[1]); CX=abs(sizes[2]);
        r = micreate_dimension("zspace", MI_DIMCLASS_SPATIAL, MI_DIMATTR_REGULARLY_SAMPLED, CZ, &hdim[0]);  
        r = micreate_dimension("yspace", MI_DIMCLASS_SPATIAL, MI_DIMATTR_REGULARLY_SAMPLED, CY, &hdim[1]);
        r = micreate_dimension("xspace", MI_DIMCLASS_SPATIAL, MI_DIMATTR_REGULARLY_SAMPLED, CX, &hdim[2]);
        if (r == MI_NOERROR);// printf("Created Dimensions.\n");
        else printf("Could not create dimensions.\n");
        
        r = miset_dimension_starts(hdim, NDIMS, altStarts);
        if (r == MI_NOERROR);// printf("Created starts.\n");
        else printf("Could not create starts.\n");
        r = miset_dimension_separations(hdim, NDIMS, altSeparations); 
        if (r == MI_NOERROR);// printf("Created separations.\n");
        else printf("Could not create seperations.\n");
        r = micreate_volume(newfilename, NDIMS, hdim, data_type, MI_CLASS_REAL, NULL, &hvol);
        r = miset_slice_scaling_flag(hvol, FALSE); /* indicate that we will be using slice scaling */
        if (r == MI_NOERROR);// printf("Set Scaling.\n");
        else printf("Could not set scaling.\n");
        r = micreate_volume_image(hvol); 
        if (r == MI_NOERROR);// printf("Created output volume.\n");
        else printf("Error creating output volume.\n");

        //printf("Dimensions for %s:\n\t %d %d %d\n\t %f %f %f\n\t %f %f %f\n\n", newfilename, sizes[1], sizes[2], sizes[3], separations[1], separations[2], separations[3], altStarts[0], altStarts[1], altStarts[2] );
	}
    miclose_volume(hvol);
}


int read_frame(data* image, int t, void* array, mitype_t  dtype ){
    mihandle_t img;
    if(image->ndim==4){
        image->wstarts[0]=t;
        image->wcount[0]=1;
    }
    miopen_volume(image->filename, MI2_OPEN_READ, &img);
    if(miget_real_value_hyperslab(img, dtype, image->wstarts, image->wcount, array) != MI_NOERROR) pexit(image->filename, "Could not read values.", 1);
    miclose_volume(img);
}

int read_dimensions(mihandle_t* img, data* volume, int ndim){
    int spatial_dim;
    midimhandle_t dim[ndim];
    volume->ndim=ndim;
    volume->wcount=malloc(ndim * sizeof(*volume->wcount));
    volume->wstarts=malloc(ndim * sizeof(*volume->wstarts)); 
    volume->count=malloc(ndim * sizeof(*volume->count));
    misize_t* temp_count=malloc(ndim*sizeof(*temp_count));
    volume->start=malloc(ndim * sizeof(*volume->start));
    volume->step=malloc(ndim * sizeof(*volume->step));
    miget_volume_dimensions( *img, MI_DIMCLASS_ANY, MI_DIMATTR_ALL, MI_DIMORDER_FILE, ndim, dim);
    if(miget_dimension_sizes(dim, ndim, temp_count /*volume->count*/) != MI_NOERROR) pexit(volume->filename, "Could not get sizes.", 1);
    for(int i=0; i<ndim; i++) volume->count[i]=temp_count[i]; 
	free(temp_count);
    if(miget_dimension_starts(dim, MI_DIMORDER_FILE, ndim, volume->start) != MI_NOERROR) pexit(volume->filename, "Could not get starts.", 1);
    if(miget_dimension_separations(dim, MI_DIMORDER_FILE, ndim, volume->step) != MI_NOERROR) pexit(volume->filename, "Could not get step sizes.", 1);
    for(int i=0; i<ndim; i++){ 
        volume->wstarts[i]=0;
        volume->wcount[i]= volume->count[i];
    } 

    //FIXME assumes dimensions must be t, z, y, x
    spatial_dim=ndim-3;
    if(ndim==4){ 
        volume->tmax=volume->count[0];
        volume->tstep=volume->step[0];
    }
    else  volume->tstep=volume->tmax=1;
    volume->zstart=volume->start[spatial_dim];
    volume->ystart=volume->start[spatial_dim+1];
    volume->xstart=volume->start[spatial_dim+2];
    volume->zstep=volume->step[spatial_dim];
    volume->ystep=volume->step[spatial_dim+1];
    volume->xstep=volume->step[spatial_dim+2];
    volume->zmax=volume->count[spatial_dim];
    volume->ymax=volume->count[spatial_dim+1];
    volume->xmax=volume->count[spatial_dim+2];
    if (ndim==4) volume->n= volume->tmax * volume->zmax * volume->ymax * volume->xmax;
    else  volume->n= volume->zmax * volume->ymax * volume->xmax;
    volume->n3d=volume->zmax * volume->ymax * volume->xmax;
    //volume->data=malloc(sizeof(*volume->data)*volume->n);
    
    //mifree_dimension_handle(dim); 
    return(0);
}

int get_dimensions(mihandle_t* img, data* volume){
    int ndim, spatial_dim;
    midimhandle_t *dim;
    if( miget_volume_dimension_count( *img, MI_DIMCLASS_ANY, MI_DIMATTR_ALL, &ndim ) != MI_NOERROR) pexit(volume->filename, "Could not read number of dimensions.", 1);
    read_dimensions(img, volume, ndim);
}

int set_pointer_start(int n1, int n2, misize_t* array1, misize_t* array2){

   return(0);
}

data* copyVolume(data* vol){
    data* out=malloc(sizeof(data));
    out->filename=vol->filename;
    out->values_filename=vol->values_filename;
    out->groups=vol->groups;
    out->data=vol->data;
	out->points=vol->points;
    out->wcount=vol->wcount;
    out->wstarts=vol->wstarts; 
    out->count=vol->count;
    out->start=vol->start;
    out->step=vol->step;
	out->val=vol->val;
    out->n_roi=vol->n_roi;
    out->ngroups=vol->ngroups;
    out->ndim=vol->ndim;
	out->index=vol->index;
    out->n=vol->n;
    out->label=vol->label;
    out->tmax=vol->tmax;
    out->zmax=vol->zmax;
    out->ymax=vol->ymax;
    out->xmax=vol->xmax;
    out->xstep=vol->xstep;
    out->ystep=vol->ystep;
    out->zstep=vol->zstep;
    out->tstep=vol->tstep;
    return(out);
}

void* readVolume(data* volume, int read_hyperslab /*1 for read, o.w. don't read*/, mitype_t dtype){
    mihandle_t img;
    void* ptr=NULL;
    int status=MI_NOERROR;
    int offset;
    int sizeof_dtype;
    
    if(dtype == MI_TYPE_INT ) sizeof_dtype=sizeof(int);
    else if(dtype == MI_TYPE_UINT) sizeof_dtype=sizeof(unsigned int);
    else if(dtype == MI_TYPE_UBYTE) sizeof_dtype=sizeof(unsigned char);
    else if(dtype == MI_TYPE_FLOAT) sizeof_dtype=sizeof(float);
    else if (dtype == MI_TYPE_DOUBLE) sizeof_dtype=sizeof(double);

    if( check_file(volume->filename ) == 0 ){
        if( miopen_volume(volume->filename, MI2_OPEN_READ, &img) != MI_NOERROR) pexit("Could not open", volume->filename, 1);
        get_dimensions(&img, volume);
        if(read_hyperslab==1) {

            ptr=malloc_check(volume->tmax, volume->zmax, volume->ymax, volume->xmax, sizeof_dtype, volume->ndim);
            status=miget_real_value_hyperslab(img, dtype, volume->wstarts, volume->wcount, ptr); 
	}
        if(status != MI_NOERROR) {
            miclose_volume(img);
            pexit(volume->filename, "Could not read values.", 1);
        }
        miclose_volume(img);
        return(ptr);
    }
    else{
        pexit(volume->filename, "Could not find file.", 1);
    }
}


int volume_to_obj(data* volume){
    int i, n; 
	int zmax=volume->zmax;
	int ymax=volume->ymax;
	int xmax=volume->xmax;
	volume->points=malloc(sizeof(*volume->points) * volume->n);
    #pragma omp parallel  
    for(int z=0; z < zmax; z++){ 
	    for(int y=0; y < volume->ymax; y++){
	    	for(int x=0; x < volume->xmax; x++){
			  i=z*ymax*xmax+y*xmax+x;
			  volume->points[i].x=x; 
			  volume->points[i].y=y; 
			  volume->points[i].z=z; 
			  volume->points[i].wx=x*volume->step[2]+volume->start[2]; 
			  volume->points[i].wy=y*volume->step[1]+volume->start[1]; 
			  volume->points[i].wz=z*volume->step[0]+volume->start[0]; 
			  //printf("a %f %f %f\t  %f %f %f\n",volume->points[i].wz, volume->points[i].wy,volume->points[i].wx,  volume->start[2], volume->step[2],  fmod(volume->points[i].wx-volume->start[2], volume->step[2]));
			  volume->points[i].index=i;
			  volume->points[i].value=volume->data[i];
			}
		}
	} 
    
    return(0);
}

int read_obj(data* surface){
    char buffer1[100],  buffer2[100]; //should probably improve buffer length so that it is not a fixed size
    int i, nvertices; 
    char* mesh_fn=surface->filename; 
    char* values_fn=surface->values_filename;
    FILE *mesh_file=fopen(mesh_fn, "rt");
    FILE *value_file=fopen(values_fn, "rt");
    char dlm[]=" "; 
    mihandle_t img;
    //Get world coordinate dimensions for mesh
    /*if( miopen_volume(surface->like_filename, MI2_OPEN_READ, &img) != MI_NOERROR) pexit("Could not open", volume->filename, 1);
    get_dimensions(&img, surface);
    miclose(img);*/
    if(mesh_file==NULL) {
        printf("Error opening %s\n", mesh_fn);
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

    nvertices=surface->n=atoi(strtok(NULL, dlm));
    surface->points = malloc(sizeof(*surface->points) * nvertices );
    
    int n=0;
    for(i=0; i< nvertices; i++){
    //while(fgets(buffer1, sizeof(buffer1), mesh_file) != NULL && fgets(buffer2, sizeof(buffer2), value_file) != NULL){
	    fgets(buffer1, sizeof(buffer1), mesh_file);
      	fgets(buffer2, sizeof(buffer2), value_file);  
        surface->points[i].x=atof(strtok(buffer1, dlm)); 
        surface->points[i].y=atof(strtok(NULL, dlm)); 
        surface->points[i].z=atof(strtok(NULL, dlm));
 	    surface->points[i].wx=surface->points[i].x; //x*surface->step[2]+surface->start[2];
  	    surface->points[i].wy=surface->points[i].y; //y*surface->step[1]+surface->start[1]; 
 	    surface->points[i].wz=surface->points[i].z; //z*surface->step[0]+surface->start[0]; 
        surface->points[i].index=i;
        surface->points[i].value=atof(buffer2);
        //printf("%f\n", surface->points[i].value);
        n++;
        //if( n > 5) exit(0);
    }
    surface->n=nvertices;
    fclose(mesh_file);
    fclose(value_file);
    return(0);
}


int check_file(char *file_input){
    if( access( file_input, R_OK ) != -1 ) {
       //printf("Can read file: %s\n", file_input);
    } else {
        printf("Error: could not access %s\n", file_input);
        return 1;
    }
    return 0 ;
}

void pexit(char* string1, char* string2, int number) {
    fprintf(stderr, "Error: %s\t%s\n", string1, string2); 
    exit(number);
}

int checkDimensions(data* vol1, data* vol2){
    int o1=0, o2=0;
    if (vol1->ndim ==4) o1=1;
    if (vol2->ndim ==4) o2=1;

    for(int i=0; i<3; i++) if( vol1->count[i+o1] != vol2->count[i+o2] ) return 1;
    for(int i=0; i<3; i++) if( vol1->start[i+o1] != vol2->start[i+o2] ) return 1;
    for(int i=0; i<3; i++) if( vol1->step[i+o1] != vol2->step[i+o2] ) return 1;
    return 0;
}


