#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include "minc2.h"
#include "hdf5.h"
#include "minc_helper.h"
//#include "volume_io.h"
/*----------------
*|VERSION HISTORY|
*-----------------
* Version 2.3: 
*   Bug: Does not give correct values when using multiple atlases, only the first atlas seems to work
*   Solution: This was fixed by correctly reading the labels stored in the images
* Version 2.4:  
*   Bug: The program does not seem to read in surfaces correctly.
*   Solution:   
*/


struct roi_aggregate {
    int n;
    int *indices;
    int **labels;
    int *nlabels;
};




struct result {
    int n;
	int roi;
    double avg;
    double max;
    double min;
    double sd;
    double vol;
};

int stats_3D(data** images,   int nimages, data* atlas, int natlas, point*** roi_list, int* roi_list_n , int n_roi, int t, int max_i, int max_k, int** nat_id, int** unique_id, int* n_nat_id, char* interpolation_method, struct result ***output );
double interpolate(double* array, int zmax, double sz, double oz, int ymax, double sy, double oy,int xmax, double sx, double ox, double xi, double yi, double zi, char* interpolation_method);
char* concatenate_groups(char** groups, int ngroups);
void useage();
int mapping(int val, int* list1, int* list2, int max);



double  interpolate(double* array, int zmax, double sz, double oz, int ymax, double sy, double oy,int xmax, double sx, double ox, double xi, double yi, double zi, char* interpolation_method){
  	if( fmod(xi-ox , sx) == 0 && fmod(yi-oy, sy) == 0 && fmod(zi-oz, sz)== 0 ) {
	  int z=(zi-oz)/sz;
	  int y=(yi-oy)/sy;
	  int x=(xi-ox)/sx;
		//printf("%d %f\n",(int) (z*ymax*xmax + y*xmax + x) , array[(int) (z*ymax*xmax + y*xmax + x) ] ); 
	  return(array[ (int) (z*ymax*xmax + y*xmax + x) ]);
	}
	else if( strcmp(interpolation_method, "-nearest") == 0) return(nearest_neighbour(array, zmax, sz, oz, ymax, sy, oy, xmax, sx, ox, xi, yi, zi));
    else if( strcmp(interpolation_method, "-linear")==0 ){ 
       return(trilinear(array, zmax, sz, oz, ymax, sy, oy, xmax, sx, ox, xi, yi, zi));
    }
    else pexit("Error: Undefined interpolation method: ", interpolation_method, 1);

    return(0.0);
}

char* concatenate_groups(char** groups, int ngroups){
    char* all_groups=NULL;

    if( ngroups > 0){
        int l=strlen(groups[0]);
        all_groups=realloc(all_groups, sizeof(*all_groups) * l);
        strcpy(all_groups, groups[0]);
        for(int i=1; i<ngroups; i++){
            l=strlen(groups[i]);
            all_groups=realloc(all_groups, sizeof(*all_groups) * l);
            strcat(all_groups, groups[i] );
        }
    }
    return(all_groups);

}


int stats_3D(data** images,   int nimages, data* atlas, int natlas, point*** roi_list, int* roi_list_n , int n_roi, int t, int max_i, int max_k, int** nat_id, int** unique_id, int* n_nat_id, char* interpolation_method, struct result ***output ){
    double sum[nimages][n_roi], sumsqrd[nimages][n_roi], avg[nimages][n_roi], sd[nimages][n_roi], max[nimages][n_roi], min[nimages][n_roi];
    int n[nimages][n_roi];
    int i, j, k, x, total;
    int l;
    char* image_names[nimages];
    char* mask_names[n_roi];
    double voxel_volume=images[0]->step[0]*images[0]->step[1]*images[0]->step[2];
    double value;
    data* cimage;
    for(j=0; j<nimages; j++){
        for(i=0; i< n_roi; i++) {
            avg[j][i]=sum[j][i]=sumsqrd[j][0]=sd[j][i]=0;
            n[j][i]=0;
            max[j][i]=min[j][i]=images[j]->data[0];
        }

    }
    int voxel, vertex;

    for(j=0; j<nimages; j++){
        for( i=0; i< images[j]->n ;i++){
            if(isnan(images[j]->data[i]) || isinf(images[j]->data[i]) ) images[j]->data[i]=0;
        }
    }
	
    for(j=0; j<nimages; j++){
	  	int d=images[j]->ndim-3;//dimension offset
        /*Iterate through index values*/
	  	for(int roi_index=0; roi_index < n_roi; roi_index++){
        	for(i=0; i < roi_list_n[roi_index]  ; i++){
			  	int index=roi_list[roi_index][i]->index;
                cimage=images[j];
				
                value=interpolate(images[j]->data, images[j]->zmax, images[j]->step[d+0], images[j]->start[d+0], images[j]->ymax, images[j]->step[d+1], images[j]->start[d+1],  images[j]->xmax, images[j]->step[d+2], images[j]->start[d+2], roi_list[roi_index][i]->wx, roi_list[roi_index][i]->wy, roi_list[roi_index][i]->wz, interpolation_method); //cimage->data[index];
                sum[j][roi_index] += value;
                sumsqrd[j][roi_index] += value*value;
                n[j][roi_index] += 1;
                if(max[j][roi_index] > value) max[j][roi_index]= value;
                if(min[j][roi_index] < value) min[j][roi_index]= value;
            }

		//printf("%d %f %f", n[j][roi_index], sum[j][roi_index], sum[j][roi_index]/n[j][roi_index] );
        }
    }
    
    
    	
    for(k=0; k<nimages; k++){
        j=images[k]->label;
        for(i=0; i < n_roi; i++){
            if(n[k][i] > 0){ 
                avg[k][i]=sum[k][i] / n[k][i];
                double variance= (sumsqrd[k][i]- pow((sum[k][i]),2)/n[k][i])/n[k][i];
                if (variance > 0.000001) { //FIXME: this shouldn't be hardcoded
                    sd[k][i]= sqrt(variance );
                } else {
                    sd[k][i]=0;
                }
            }
            else{ 
                avg[k][i]=0;
                sd[k][i]=0;
            }
            output[j][i][t].avg=avg[k][i];
            output[j][i][t].sd=sd[k][i];
            output[j][i][t].n=n[k][i];
            output[j][i][t].vol=n[k][i];
            output[j][i][t].min=min[k][i];
            output[j][i][t].max=max[k][i];
            total=0;
            //Find the mask number for roi i
            for(x=0; x < natlas; x++){ 
                total+=n_nat_id[x]; 
                if(i < total) break;
            }
			output[j][i][t].roi=mapping(i, unique_id[x], nat_id[x],  n_nat_id[x]);
        }
    }
    
    return(0);
}

void print_output(FILE* output_f, data* images, data* atlas, int natlas, struct result ***output, int a /*atlas number*/, int i/*roi number*/, int j/*image number*/, int t, int mean_only ){
  	int total=0;
  	if(output_f == NULL){
	  	/*for(int g=0; g<images[j].ngroups; g++) printf("%s,", images[j].groups[g]); 
		for(int a=0; a<natlas; a++){
			total+=atlas[a].ngroups;
		 	if( j <= total){
				for( int g=0; g<atlas[a].ngroups; g++){
					printf("%s,", atlas[a].groups[g]);
				}
		 		break;
		 	}
		}*/
		if(mean_only != TRUE){ 
            for( int g=0; g<images[j].ngroups; g++) printf("%s,", images[j].groups[g]); 
            for( int g=0; g<atlas[a].ngroups; g++) printf("%s,", atlas[a].groups[g]);
            printf( "%d,%d,%d,",images[j].ndim, output[j][i][t].roi, t );
        }
        printf( "%f",output[j][i][t].avg);
        if( mean_only != TRUE) printf( ",%f,%f,%f,%f", output[j][i][t].sd, output[j][i][t].min, output[j][i][t].max, output[j][i][t].vol );
        printf("\n");

	}
	else{
		//for(int g=0; g<images[j].ngroups; g++) fprintf(output_f, "%s,", images[j].groups[g]); 
	    //for(int k=0; k < natlas; k++) for( int g=0; g<atlas[k].ngroups; g++) fprintf(output_f, "%s,", atlas[k].groups[g]); 
	  	//fprintf(output_f, "%d,%d,%d,%f,%f,%f,%f,%f\n",images[j].ndim, output[j][i][t].roi,   t, output[j][i][t].avg, output[j][i][t].sd, output[j][i][t].min, output[j][i][t].max, output[j][i][t].vol );
		for( int g=0; g<images[j].ngroups; g++) fprintf(output_f,"%s,", images[j].groups[g]); 
        for( int g=0; g<atlas[a].ngroups; g++) fprintf(output_f,"%s,", atlas[a].groups[g]);
        fprintf(output_f, "%d,%d,%d,%f,%f,%f,%f,%f\n",images[j].ndim, output[j][i][t].roi,   t, output[j][i][t].avg, output[j][i][t].sd, output[j][i][t].min, output[j][i][t].max, output[j][i][t].vol );
    }
}

int get_descriptive_stats(data* images, int nimages, data* atlas /**/, int natlas, int nframes, point*** roi_list, int* roi_list_n , int n_roi, int** nat_id, int** unique_id,  int* n_nat_id,  char* interpolation_method, char* output_fn, int mean_only){
    int i, j, t, index, n_temp_images;
    unsigned int total=0;
    data** temp_images=NULL;
    struct result ***output=malloc(sizeof(*output) * nimages );
    int max_i, max_k;
	FILE* output_f=NULL;
	if(output_fn != NULL){
	  output_f = fopen(output_fn, "wb");
	  if (output_f==NULL) pexit("Error: could not write to", output_fn , 1);
	}

    for (i=0; i< nimages; i++) {
        output[i]=malloc(sizeof(**output) * (n_roi));
        for (j=0; j< n_roi; j++) output[i][j]=NULL;
    }
	if (nframes <= 0) nframes=1;

    //Initialize output: output[nimages][nroi][nframes in image]
    for(i=0; i<nimages; i++){
        output[i]=malloc(sizeof(**output)*n_roi);
        for (j=0; j< n_roi; j++){
            output[i][j]=malloc(sizeof(***output)*images[i].tmax);
        }
    }

	//printf("Frame: ");	
	for(t=0; t< nframes; t++){
        //printf("%d\n", t);
        n_temp_images=0;
        //Find images that still need to be analyzed
        for(i=0; i<nimages; i++){
            //image must have more frames than t
            if ( images[i].tmax >= t) {
                read_frame( &(images[i]), t, images[i].data, MI_TYPE_DOUBLE );
                n_temp_images++; //Have to dynamically allocate number of temp images because they may not all have the same number of time frames
                temp_images=realloc(temp_images, sizeof(*temp_images) * n_temp_images);
                temp_images[n_temp_images-1] = &(images[i]);
                //Have to allocate memory for output
                
            }
        }
        stats_3D(temp_images,  n_temp_images, atlas, natlas, roi_list, roi_list_n , n_roi, t, max_i, max_k, nat_id, unique_id, n_nat_id, interpolation_method, output);
        if(output_f != NULL) {printf("\rCompleted: %3.1f\% ",(float) t/nframes * 100); fflush(stdout);}
        
    }

	for(t=0; t< nframes; t++){
        for(j=0; j<nimages; j++){
            if( images[j].tmax>=t){/*Iterate over images*/
                for(int a=0; a < natlas; a++){
                    for(i=0; i<atlas[a].n_roi; i++){ /*Iterate over ROI*/
                        print_output(output_f, images, atlas, natlas, output, a, unique_id[a][i], j, t, mean_only);
                    }
                }
            }
        }
    }


	printf("\n");
    char* short_fn;
 
    free(temp_images);
    for (i=0; i< nimages; i++) {
        for (j=0; j< n_roi; j++) free(output[i][j]);
        free(output[i]); 
    }
    free(output);
    return(0);
}



int isin(int val, int* unique_list, int nunique){
	for(int i=0; i< nunique; i++) if( unique_list[i] == val) return(1);
	return(0);
}

int append_to_list(int value, int* list, int length){

	return(0);
}

int mapping(int val, int* list1, int* list2, int max){


    for(int i=0; i<max; i++)
        if(list1[i] == val) return list2[i];


  printf("Error: Should have found equivalent value for %d\n.", val);
  exit(1);

}


int make_unique(data* atlas,  int natlas, int* n_roi, int*** nat_id /*Native id values in atlases*/, int **nunique /*Number of unique id values in each atlas*/, int*** unique_id){
  	///EXAMPLE:
  	///Masks may be loaded as atlases that contain multiple labeled regions of interest
	///that are not unique across volumes. For example, the label "1" may be used for the
	///motor cortex in one mask but also be used to identify the amygdala in another mask.
	///			Old Labels		New and Unique Labels
  	///ROI A: 	1,5				1,2
	///ROI B:   1,2				3,4
	///ROI C:   3,4,5			5,6,7
	int total_unique=0;
	int unq_val=0;
    int atlas_n_roi=0;
	int counter=0, new_label, val;
    int** nat_id_ptr, *nunique_ptr, **unique_id_ptr;
    
	*nat_id=malloc(natlas*sizeof(**nat_id));
	*unique_id=malloc(natlas*sizeof(**unique_id));
    nat_id_ptr=*nat_id;
    unique_id_ptr=*unique_id;
    for(int i=0; i<natlas; i++){ 
       unique_id_ptr[i]=nat_id_ptr[i]=NULL;
    }
	*nunique=calloc(natlas, sizeof(**nunique));
    nunique_ptr=*nunique;
	//Create list of ROI values
    //for(int i=0; i< natlas; i++) for(int j=0; j < atlas[i].n; j++) printf("", atlas[i].points[j].value);

  	for(int i=0; i < natlas; i++){
        atlas_n_roi=0;
		for(int v =0; v < atlas[i].n ; v++){
		    //get value of roi i at location v
		  	val= (int) round(atlas[i].points[v].value);
            if( val != 0 ){
			  	//Initialize for first value
				if(isin(val, nat_id_ptr[i], nunique_ptr[i]) == 0 ){
			  		nunique_ptr[i] += 1; //increase the number of unique values by one
					nat_id_ptr[i]=realloc(nat_id_ptr[i], sizeof(**nat_id_ptr) * nunique_ptr[i] );
					unique_id_ptr[i]=realloc(unique_id_ptr[i], sizeof(**unique_id_ptr) *  nunique_ptr[i] );
                  //  printf("Val: %d, nunique[%d]=%d\n", val, i, nunique_ptr[i]);
                    nat_id_ptr[i][nunique_ptr[i]-1]=val;
                    unique_id_ptr[i][nunique_ptr[i]-1]=unq_val;
					unq_val++;
                    atlas_n_roi++;
				}
			}
		}
        atlas[i].n_roi=atlas_n_roi;
	}  
	
    /*for(int i=0; i < natlas; i++){
        for(int j=0; j< nunique_ptr[i] ; j++){
		    printf("%d -->  %d\n", nat_id_ptr[i][j],  unique_id_ptr[i][j] );
	    }
    }/**/ 
	*n_roi=unq_val;

	/* Example:
	 *------------------------------------------
	 *|Atlas	|	Nat ID	|	Unique ID		|
	 *|-----------------------------------------|
	 *|0		|	1,3,5	|	1,2,3,4,5,6,7,8	|
	 *|1		|	1,2		|					|
	 *|2		|	2,3,4,5	|					|
	 *-------------------------------------------
	 */
	return(0); 

}

int aggregate_roi(data* atlas, int* n_roi, int natlas,  point**** roi_list /*list of points for each roi*/, int** roi_list_n /*number of points in each roi*/,   
					int*** nat_id /*list of ids in atlases*/, int** n_nat_id /*number of ids in each atlas*/, int*** unique_id ){
    /*Aggrate the ROIs in volumetric or surface-based masks into a single array*/
    /*This is done by creating a list of the label values that uniquely species id values in each of the atlases. For each unique label
    *value, a list of points is created that contains all the points specified by that label.*/
    int n=0;
    int l=0;
    int** nat_id_ptr, *n_nat_id_ptr, **unique_id_ptr, *roi_list_n_ptr;
	
    make_unique(atlas, natlas, n_roi, nat_id, n_nat_id, unique_id);
	unique_id_ptr=*unique_id;
    nat_id_ptr=*nat_id;
    n_nat_id_ptr=*n_nat_id; 
    *roi_list=malloc(sizeof(**roi_list) * *n_roi);
	*roi_list_n=calloc(*n_roi, sizeof(**roi_list_n));
    point*** roi_list_ptr=*roi_list;
	for(int i=0; i < *n_roi; i++) roi_list_ptr[i]=NULL;
    roi_list_n_ptr=*roi_list_n;
	/***********************************************
    *EXAMPLE:
    *   roi->indices    roi->labels     roi->nlabels
    *  -----------------------------------------------
    *   4               1               1
    *   32              2, 3            2
    *   444             4               1
    *   1234            5, 6, 7         3
    *   1555            8, 9            2
    ************************************************/

	/*For atlas in list of atlases */
	for(int atlas_i=0; atlas_i < natlas; atlas_i++){
		/*For i in atlas */
		for(int index=0; index < atlas[atlas_i].n; index++){
			/*If value for atlas at i does not equal 0...*/
			int id=round(atlas[atlas_i].points[index].value);	
			if( id != 0){
				/*find the unique id for i in atlas*/
				int uid = mapping(id, nat_id_ptr[atlas_i], unique_id_ptr[atlas_i], n_nat_id_ptr[atlas_i] );
				/*and add it to list of roi*/
				roi_list_n_ptr[uid] += 1;
				int ni=roi_list_n_ptr[uid]; //ni is the number of points covered by the unique id number, uid.
                //printf("%d %d\n", uid,ni );
				roi_list_ptr[uid] = realloc( roi_list_ptr[uid], sizeof(*roi_list_ptr)*  ni );
				roi_list_ptr[uid][ni-1]=&(atlas[atlas_i].points[index]); //set pointer to point in atlas i 

            }
		}
	}

    return(0);
}


int make_ngroups_equal(data** vol,int nvol){
    int max_ngroups=0;

    for(int i=0; i < nvol; i++){ 
        if( (*vol)[i].ngroups > max_ngroups) max_ngroups=(*vol)[i].ngroups;
    }
    for(int i=0; i < nvol; i++){
        char*** groups=&((*vol)[i].groups);
        int* ngroups=&( (*vol)[i].ngroups);
        if( *ngroups < max_ngroups){
            *groups=realloc(*groups, sizeof(**groups) * max_ngroups);
            if(groups == NULL) pexit("Could not allocated memory for groups associated with file:", (*vol)[i].filename, 1 );
            for(int j=*ngroups; j < max_ngroups; j++){
                (*groups)[j]="empty";
            }
         *ngroups=max_ngroups;    
        }
    }
    

    return(0);
}



int isgroup(char** argv, int* i, char*** groups, int* ngroups){
    if(strcmp(argv[*i], "-g") == 0 ){
        *i += 1;
        *ngroups += 1;
        *groups=realloc(*groups, sizeof(**groups) * *ngroups);
        if(groups == NULL) pexit("Could not allocated memory for groups associated with file:", argv[*i], 1 );
        (*groups)[*ngroups-1]=argv[*i];
        return(0);
    }
    else return(1);

}

int add_groups(char** argv, int* i, char*** groups, int* ngroups){
    while( strcmp(argv[*i], "-g") == 0 ){
        *i += 1;
        *ngroups += 1;
        *groups=realloc(*groups, sizeof(**groups) * *ngroups);
        if(groups == NULL) pexit("Could not allocated memory for groups associated with file:", argv[*i], 1 );
        (*groups)[*ngroups-1]=argv[*i];
        *i += 1;
    }
    
    return(0);
}



int parse_sub_input(int *i, char** argv, data* ptr){
  	//function to go through sub-inputs for images, vol. masks and surf. masks, given by user
    int out=1;
    char* image_type=argv[*i];
    *i += 1;
    //if out equals 0, then we found additional parameters
    //if out equals 1, then we did not find additional parameters
    add_groups(argv, i, &(ptr->groups), &(ptr->ngroups));
    if( strcmp(image_type, "-v")==0  || strcmp(image_type, "-i")==0  ){
        ptr->filename=argv[*i];
        if(check_file(ptr->filename) != MI_NOERROR ) pexit("Error: could not read", ptr->filename , 1 );
        return(1);
    }
    else if (strcmp(image_type, "-s")==0){
        ptr->filename=argv[*i];
        *i += 1;
        //ptr->like_filename=argv[*i];
        //*i += 1;
        ptr->values_filename=argv[*i];
        if(check_file(ptr->filename) != MI_NOERROR) pexit("Error: could not read", ptr->filename, 1 );
        if(check_file(ptr->values_filename) != MI_NOERROR ) pexit("Error: could not read", ptr->values_filename, 1 );
        return(1);
    }
    *i+=1;
    return(0);
}


int parse_input(int argc, char** argv, data** images, int* nimages, data** atlas, int *natlas, char** interpolation_method, char** output_fn, int* mean_only ){
    int i;
    char* linear="-linear";
    char* nearest="-nearest";
    char* mean_only_str="-mean_only";
    *interpolation_method=linear;
    *atlas=NULL;
    data* tmp=NULL;
    for( i=1; i < argc; i++ ){
        if(strcmp(argv[i], "-i" ) == 0 ){
            allocate_vol(&(*images), nimages);
            parse_sub_input(&i, argv, &((*images)[*nimages-1]) );
        }
        else if(strcmp(argv[i], "-v")  == 0 || strcmp(argv[i], "-s" ) ==0 ){
            allocate_atlas( &(*atlas), natlas);
            parse_sub_input(&i, argv,  &((*atlas)[*natlas-1]) );
        }
       else if(strcmp(argv[i], "-o" )==0 ){
        	i++;
			*output_fn=argv[i];
		}
        else if( strcmp(argv[i], nearest) == 0  ){
            *interpolation_method=nearest;
        }
        else if( strcmp(argv[i], mean_only_str) == 0  ){
            *mean_only=TRUE;
        }
    }
    make_ngroups_equal(images, *nimages);
    make_ngroups_equal(atlas, *natlas);

return(0);
}

int read_atlas( data* atlas ){
  	//If there is no values filename, then we are dealing with a volume
  	if( atlas->values_filename == NULL ){
		atlas->data=(double*)readVolume(atlas, 1, MI_TYPE_DOUBLE);
  	} else{
  	//if there is a values filename, then we have a surface
		read_obj(atlas);
  	}

	return(MI_NOERROR);
}

//Get statistics for 3D or 4D images with an arbitrary number of masks
int main(int argc, char** argv){
    VERBOSE=FALSE;
    if(argc <= 4 || strcmp(argv[1], "-help") == 0 ){
	  	useage();
	}
    int n_aux_inputs=2;
    char* interpolation_method=NULL;
	char* output_fn=NULL;
    int nimages=0;
    int natlas=0;
    int nframes=0;
    int spatial_dim=0;
    int nmask_voxels;
	int n_roi=0;
    int mean_only=FALSE;
	int* n_nat_id=NULL;
    int** nat_id=NULL;
    int** unique_id=NULL;
    data* images=NULL;
    data* atlas=NULL;
    struct roi_aggregate atlas_voxels;
	point*** roi_list=NULL; /*list of points for each roi*/
   	int* roi_list_n=NULL; /*number of points in each roi*/
    parse_input(argc, argv, &images, &nimages, &atlas, &natlas, &interpolation_method, &output_fn, &mean_only );
    //printf("Received Inputs:\nImages\tVolumes\tSurfs\n%d\t%d\n", nimages, natlas);
    if(nimages == 0 || natlas==0 ) useage();
    //load volumetric masks 
	//printf("Reading volume\n");
	for(int i=0; i< natlas; i++){
    	if( read_atlas(&atlas[i]) == MI_ERROR  ) pexit("Could not read data from" , atlas[i].filename,  1);
		if(atlas[i].values_filename == NULL) {
            volume_to_obj( &(atlas[i]));
           
		}
	}
	//printf("Aggregating ROI for volume\n");
    aggregate_roi(atlas, &n_roi, natlas, &roi_list, &roi_list_n, &nat_id, &n_nat_id, &unique_id );
    //FIXME might not have same number of vertices for all meshes
    //printf("Nimages= %d, Nmasks=%d\n", nimages, natlas);
    //printf("Output: %s\n", output_fn);
    //Load the dimensions for the images
     
    
    for(int i=0; i<nimages; i++){
        mihandle_t img;
		//printf("Opening image: %s\n", images[i].filename);
        miopen_volume(images[i].filename, MI2_OPEN_READ, &img);
        get_dimensions(&img, &images[i]);
        miclose_volume(img);
        //FIXME assumes that the first dimension is time, but this doesn't have to be the case
        if(images[i].tmax > nframes ) nframes=images[i].tmax; //update maximum number of frames.
		images[i].data=malloc_check(nframes, images[i].zmax, images[i].ymax, images[i].xmax, sizeof(double), images[i].ndim);
    }
	//printf("Getting descriptive stats"); fflush(stdout);
    get_descriptive_stats(images, nimages, atlas, natlas, nframes, roi_list, roi_list_n, n_roi, nat_id, unique_id, n_nat_id, interpolation_method, output_fn, mean_only); 
    return(0);
}



void useage(){
        printf("\nName:\n\tmincgroupstats -- get descriptive statistics from multiple images with multiple regions of interest\n");
        printf("Description:\n");
        printf("\tGet descriptive statistics (mean, std. dev., min, max, volume) from volumetric images using mask and atlases (defined both volumetriclly and on surfaces). All volumetric files must have the same dimensions and the surface .obj files must be defined in the same space as the volumetric images.\n");
        printf("Useage:\n");
        printf("\t-i <.mnc>\t\tTarget volumetric MINC2 file.\n");
        printf("\t-v <.mnc>\t\tVolume or surface atlas or mask.\n");
        printf("\t\t\t\tVolumetric atlas/mask: -v <Volume.mnc>\n");
        printf("\t\t\t\tSurface atlas/mask: -v <surface.obj>  <roi.txt>\n");
        printf("\t-g <label>\t\tAdd label to image or atlas.\n");
        printf("\t\t\t\tAdd \'-g <label>\' after \'-v\' but before the mask information.\n");
        printf("Options:\n");
        printf("\t-nearest\t\tUse nearest neighbour interpolation.\n");
        printf("\t-linear\t\t\tUse trilinear interpolation (default).\n");
        printf("\nExample:\n\tmincgroupstats -i -g PET pet.mnc -i -g fMRI fmri.mnc -v -g Classify .mnc -s -g BA3 midsurface.obj surface_roi.txt\n\n");
        printf("\tOutput:\n");
        printf("\tLabel1\tLabel2\t\tnFrames\tROI#\tFrame\tMean\tStd. Dev.\tMin\tMax\tVolume\n");
        printf("\tPET\tClassify\t3\t3\t0\t...\n");
        printf("\tPET\tClassify\t3\t2\t0\t...\n");
        printf("\tPET\tClassify\t3\t1\t0\t...\n");
        printf("\tPET\tBA3\t\t3\t1\t0\t...\n");
        printf("\tfMRI\tClassify\t3\t3\t0\t...\n");
        printf("\tfMRI\tClassify\t3\t2\t0\t...\n");
        printf("\tfMRI\tClassify\t3\t1\t0\t...\n");
        printf("\tfMRI\tBA3\t\t3\t1\t0\t...\n");

        exit(1);
}


