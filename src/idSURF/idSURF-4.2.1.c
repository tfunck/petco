#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include "pthread.h"
#include "minc2.h"
#include "minc_helper.h"
#include "gaussianiir3d.h"

void useage();
data** parse_input(int argc, char** argv, int* nthreads, int* nimages, int *nmasks, int* max_iterations, int* max_avg, float* fwhm, float* tolerance, char** outputfilename, int* smooth_only, int* cubic_averaging_regiong );
int local_allocate_vol(data* ptr, int* n, char* argv_next);
void* idSURF_threaded(void* args);

void* idSURF(data* image, mihandle_t* outputimage, char* outputfilename, int*** regions, int *n_avg_region,  int* mask_label,  int max_iterations, double fwhm, int p, double tolerance, data* masks, int max_avg,   int smooth_only, data* first_guess, int thread, int nthreads);


int** copy_int_2darray(int** input, int dim1, int dim2 ){
    int** output=malloc(sizeof(*input) * dim1);
    for(int i=0; i<dim1; i++){ 
        output[i]=malloc(sizeof(**input)*dim2 );
        for(int j=0; j<dim2; j++) output[i][j]=input[i][j];
    }
    return(output);
}

int isin(int border_index_ngh, int** region, int nregion){
    for(int i=0; i< nregion; i++){ 
        if( border_index_ngh == region[i][0]){ 
            return(1);
        }
    }
    return(0);
}

int isin2(int x, int* list, int n){
    for(int i=0; i< n; i++) if( x == list[i] ) return(1);
    return(0);
}
int swap_border(int** border, int i, int j){
    int temp, tempd;
    temp=border[i][0];
    tempd=border[i][1];
    border[i][0]=border[j][0];
    border[i][1]=border[j][1];
    border[j][0]=temp;
    border[j][1]=tempd;
    return(0);
}

int sort_by_distance(int** array, int n){
    int nswaps=1;

    while(nswaps > 0){
        nswaps=0; //set number of swaps to 0
        for(int i=0; i<n-1; i++){
            if(array[i][1] > array[i+1][1]){
                //if element i+1 is larger than element i, swap them and their distances
               swap_border(array, i, i+1);
               nswaps+=1;
            }
        }
    }
    return(0);
}

int append_int_array(int*** region, int voxel, int nregion,int** border, int nborder){
    int** tmp=realloc(region[voxel], sizeof(**region) * (nregion+nborder));
    if(tmp == NULL) return(1);
    region[voxel]=(int**) tmp;
    for(int i=nregion; i<nregion+nborder; i++){ 
        region[voxel][i]=malloc(sizeof(***region) * 2);
        if(region[voxel][i] ==NULL) return(1);
        for(int j=0; j<2; j++){ 
            region[voxel][i][j]=border[i-nregion][j];
        }
    }

    return(0);
}
/********************************************************************************************
*Name: create_border 
*Purpose:   
*Method:     
*Inputs:    
*Output:    
*********************************************************************************************/
int** create_border(int* ngh, int nngh, int depth){
    int** border=malloc(sizeof(*border) * nngh );
    for(int j=0; j< nngh; j++) border[j]=malloc(sizeof(**border)*2);
    for(int j=0; j< nngh; j++) {
        border[j][0]=ngh[j];
        border[j][1]=depth;
    }
    return(border);
}


/********************************************************************************************
*Name: initialize_region 
*Purpose:   Allocate memory for a 2D array of integers, "regions[i]", where i represents the 
            index value for a given image voxel. It will serve to keep track of which voxels
            will be averaged over when smoothing voxel i. The first dimension is allocated 
            based on the number of voxels that are within the averaging region of voxel i.
            Initially, the averaging region contains the immeadiately neighbouring voxels to
            voxel i. The second dimension of the array contains two entries, the first stores
            the index value of the jth voxel in the averaging region. The second stores the 
            "depth" (i.e., how many iterations had to be performed to reach voxel j) associated
            with voxel j. The depth starts at 1, because we are on the first iteration.  
*Method:     
*Inputs:    
*Output:    
*********************************************************************************************/

int initialize_region(int*** region, int** ngh, int i, int nregion, int* depth){
    region[i]=malloc(sizeof(**region)*nregion);
    if(region[i] == NULL) return(1);
    for(int j=0; j<nregion; j++)  {
        region[i][j]=malloc(sizeof(***region)*2);
        if( region[i][j]==NULL) return(1);
        if(j==0){
            region[i][j][0]=i;
            region[i][j][1]=*depth;
            *depth += 1; 
        } else{
            region[i][j][0]=ngh[i][j-1];
            region[i][j][1]=*depth; 
        }
    }

    return(0);
}

int get_box(int x, int m, int max){
    int mid=  (m+1)/2;
    int value=0;
    if(max==0){ /*find minimum*/
        value= x-mid+1;
        if(value < 0) return(0); /*bounds of box go below 0*/
        else return(value);
    } else if (max > 1) {
        value =x+mid;
        if(value < max) return(m); 
        else return(max-x-1); /*bounds of box go above max*/

    } else pexit("Error: get_box does not have <bound> set to positive integer", "", 1);

}

/********************************************************************************************
*Name: cubic_averaging_region
*Purpose:   For each voxel within a mask, find neighbouring voxels over which we will average 
*Method:     
*Inputs:    
*Output:    
*********************************************************************************************/
int find_cubic_averaging_region(data* img, int*** region, int* n_region,  double* mask, int max_ngh){
    int nvoxels=img->n;
    int nngh;
    printf("Cubic averaging region\n");
    for(int z=0; z< img->zmax; z++){
        //printf("%d\n", z);
        for(int y=0; y< img->ymax; y++){
            for(int x=0; x< img->xmax; x++){
                int xlmin=get_box(x, max_ngh, 0);
                int ylmin=get_box(y, max_ngh, 0);
                int zlmin=get_box(z, max_ngh, 0);
                int xlmax=get_box(x, max_ngh, img->zmax);
                int ylmax=get_box(y, max_ngh, img->ymax);
                int zlmax=get_box(z, max_ngh, img->xmax);
                int index=z*img->ymax*img->xmax + y*img->xmax + x;
                nngh=0;
                //printf("\nx: %d %d\n", xlmin, xlmax);
                //printf("y: %d %d\n", ylmin, ylmax);
                //printf("z: %d %d\n", zlmin, zlmax);
                if( mask[index] != 0 ){
                    for(int k=0 ; k < zlmax ; k++){
                        for(int j=0 ; j < ylmax ; j++){
                            for(int i=0 ; i <  xlmax; i++){
                                int t=(zlmin+k)*img->ymax*img->xmax+(ylmin+j)*img->xmax + (xlmin+i);
                                if( mask[t] != 0 ){
                                    nngh++;
                                    region[index]=realloc(region[index], sizeof(*region[index])*nngh);
                                    region[index][nngh-1]=malloc(2*sizeof(**region[index]));
                                    if(region[index] == NULL ) pexit("Could not allocate memory in function <find_cubic_averaging_region", "", 1);
                                    if(mask[t] == 0 ) printf("mask[%d]=%d\n", index, t); 
                                    region[index][nngh-1][0]=t;
                                    region[index][nngh-1][1]=1;
                                    n_region[index]=nngh;
                                }     
                            }
                        }
                    }/**/
                //printf("nngh[%d] = %d\n", index, n_region[index]);
                }
            }
        }
    }
    return(0);
}



/********************************************************************************************
*Name: find_averaging_region
*Purpose:   For each voxel within a mask, find neighbouring voxels over which we will average 
*Method:    Uses while loop to expand outwards from central voxel and find neighbouring voxels.
            This is done by expanding a border outwards and adding these values to the list of 
            points region. This continues until more than the maximum number of neighbouring 
            voxels is found or there are no more border points.
            Only the closest neighbouring voxels are kept such that we end up with max_ngh 
            neighbouring voxels to average over. 
*Inputs:    
*Output:   region:
                dim1=voxels in mask
                dim2=neighbouring voxels over which we will average 
                dim3=voxel index, voxel depth (distance from central voxel)
*********************************************************************************************/
int find_averaging_region(int*** region, int* n_region, int* mask_label, int **ngh, int* nngh, int nvoxels, double* mask, int max_ngh){
    int i,z,y,x, border_index, border_index_ngh;
    int nborder=0, old_nborder, depth, nregion;
    int **border=NULL, **old_border;
    
    //for all voxels in the image
    for(i=0; i<nvoxels; i++){
        if(mask[i] != 0 && nngh[i] > 0){
            //printf("Voxel: %d\n", i);
            //initial number of points on border
            nborder=nngh[i];
            //initial number of points in averaging region
            nregion=1+nngh[i];
            //depth of search starts at 0
            depth=1;
            //initialize region
            //initialize_region(int*** region, int** ngh, int i, int n_init, int* depth)
            if(  initialize_region( region, ngh, i, nregion, &depth) != 0) printf("Could not allocate memory to initialize region\n");
            //the first border is simply the neighbours of the current voxel
            border=create_border(ngh[i], nngh[i], depth); 
            //append these first border points to averaging region
            if( append_int_array(region, i, nregion, border, nborder) != 0) printf("Could not allocate memory to append array.\n");

            //while the number of voxels in averaging region is less than desired number
            //of maximum number of averaging neighbours.
            while( nregion < max_ngh){
                //printf("\tDepth: %d\n", depth);
                //copy new border to old border
                old_border=copy_int_2darray(border, nborder, 2);
                old_nborder=nborder;
                if(old_nborder == 0) {
                    //for(int j=0; j < nregion; j++) free(region[i][j]);
                    //free(region[i]);
                    //region[i]=NULL;
                    mask[i]=-1; 
                    break;
                }
                //printf("\told_nborder: %d\n", old_nborder);
                //free border so that we can create a new border 
                for(int j=0; j < nborder; j++) free(border[j]);
                free(border);
                //set border to NULL so that we can realloc 
                border=NULL;
                //set number of border points to 0
                nborder=0;
                //increase the depth of the search
                depth++;
                //for the old border points...
                for(int j=0; j < old_nborder; j++){
                    //index value of border point
                    border_index=old_border[j][0]; 
                    //for neighbours of border point
                    for(int k=0; k < nngh[border_index]; k++ ){
                        //index value for neighbour of border voxel
                        border_index_ngh= ngh[border_index][k];
                        //if the current border point's neigbour is not yet in region
                        if( isin( border_index_ngh, region[i], nregion) == 0 && isin( border_index_ngh, border, nborder) == 0    ){
                            nborder++;
                            //printf("\t\t\tnborder=%d\t", nborder);
                            int** tmp=realloc(border, sizeof(*border) * nborder);
                            if(tmp==NULL) printf("Error allocating memory for border");
                            else border=tmp;
                            int* tmp2=malloc( sizeof(**border) * 2);
                            if(tmp2==NULL) printf("Error allocating memory for border");
                            else border[nborder-1]=tmp2;
                            //printf("index=%d\n", border_index_ngh  );
                            border[nborder-1][0]=border_index_ngh;
                            border[nborder-1][1]=depth; //assign distance (i.e. depth of search)
                        }
                    }
                }
                //printf("\tNew border points: %d\n", nborder);
                for(int j=0; j < old_nborder; j++) free(old_border[j]);
                free(old_border);
                old_border=NULL;
                
                //we now have a new border points that we can add to the averaging region
                //for(int v=0; v< nregion; v++) printf("\t1 region -> %d\n", region[i][v][0]);
                if( append_int_array(region, i, nregion, border, nborder)!= 0) printf("Could not allocate memory to appent array.\n");
                nregion += nborder;
            }
            if(mask[i] != -1){
                sort_by_distance(region[i], nregion);
                region[i]=realloc(region[i], sizeof(*region)*max_ngh);
                n_region[i]=max_ngh;
            }
            else{
                //printf("Deviant region: %d\n", nregion);
                n_region[i]=nregion;
                mask_label[i]  *= -1;
            }
        for(int j=0; j < nborder; j++) free(border[j]);
        free(border);
        }
    }

    return(0);
}

/********************************************************************************************
*Name: find_neighbours
*Purpose:   1) Find immeadiate neighbours surrounding voxel points
            2) At the same time, consolidate information from multiple masks into a single 
            array, <mask_label>    
*Method:    Neighbouring points are stored in <ngh>, with <nngh> keeping track of how many
            neighbouring voxels there are. 
*Inputs:    
*Output:    
*********************************************************************************************/
int find_neighbours(int max_ngh, int** ngh, int* nngh, int* mask_label, int label, data* mask){
    //variables
    int kernelx[]={1, -1, 0, 0, 0, 0};
    int kernely[]={0, 0, 1, -1, 0, 0};
    int kernelz[]={0, 0, 0, 0, 1, -1};
    int zmax=mask->zmax;
    int ymax=mask->ymax;
    int xmax=mask->xmax;
    int index;
    int voxel;
    for(int z=0; z < zmax; z++){
        for(int y=0; y < ymax; y++){
            for(int x=0; x < xmax; x++){
                voxel=z*ymax*xmax + y*xmax + x;
                if ( mask->data[voxel] != 0 ) {
                    mask_label[voxel]=label+1; //Set mask_label to consolidate info from multiple masks into a single array
                    for(int i=0; i< 6; i++){ //Iterate through 6 neighbouring voxels surrounding z,y,x
                        int xoffset=kernelx[i];
                        int yoffset=kernely[i];
                        int zoffset=kernelz[i];
                        if( z + zoffset >= 0 && y+yoffset >= 0 && x+xoffset >= 0 ){
                            index=(z+zoffset)*ymax*xmax + (y+yoffset)*xmax + (x+xoffset);
                            if (mask->data[index] != 0){
                                nngh[voxel]++;
                                //add neighbour 
                                int size=sizeof(*(ngh[voxel]));
                                ngh[voxel]=realloc(ngh[voxel], size * nngh[voxel]);
                                ngh[voxel][nngh[voxel]-1]=index;
                            }
                        }
                    }
                }
                if(nngh[voxel]==0){
                    free(ngh[voxel]);
                    mask_label[voxel]=0;
                    nngh[voxel]=0;
                }
            }
        }
    }


    return(0);
}



int get_total_weight(double* total_weight_list, int*** region, int nvox, int p, int* n_avg_region){
    for(int i=0; i < nvox; i++){
        if( region[i] != NULL){
            total_weight_list[i]=0;
            for(int j=0; j< n_avg_region[i]; j++) total_weight_list[i] += 1/pow(region[i][j][1], p);
        }
    }
    return(0);
}


double idw(int **region, double* data, double total_weight, int max_ngh, int p){
    double value=0, weight;
    int index;
    for(int j=0; j< max_ngh; j++){
        index=region[j][0];
        weight=1/pow(region[j][1], p);
        value += data[index]  * weight;
    }
    return(value/total_weight);
}

/********************************************************************************************
*Name: avg
*Inputs:    
*Purpose:   Averages over image values in set region 
*Method:    Locations are stored as index values in <region> array, which are then used to 
*           access the corresponding image values. 
*Output:    Average, stored as double
*********************************************************************************************/

double avg(int **region, double* data,  int max_ngh ){
    double value=0;
    int index;
    for(int j=0; j< max_ngh; j++){
        index=region[j][0];
        value += data[index]  ;
    }
    //printf("%f %d\n", value, max_ngh);
    return(value/max_ngh);
}

int copy_with_threshold(double* test, double*  observed, int* mask_label, int  nvox, double* energy, double* max, int threshold){
    *energy = 0;
    *max = 0;

    for(int i=0; i<nvox; i++) {//for all voxels
        if(observed[i] > *max){ //calculate the maximum observed value
            *max=observed[i]; 
        }
        if(mask_label[i] > 0 ){
            test[i]=observed[i];
            *energy += test[i]*test[i];
        }
        else{
            test[i]=0;
        }
    }

    return(0);
}


int SURF2(double* data, int***  region, int* n_avg_region,   int nvox, int* mask_label, int label){
    double total_weight;
    double value;
    double* temp=calloc(nvox, sizeof(*temp));

    for(int i=0; i < nvox; i++){
        if(n_avg_region[i] > 0 && mask_label[i]==label ) {
            value=avg(region[i], data, n_avg_region[i]);
        } else value=data[i];
        temp[i]=value;
    }
    for(int i=0; i<nvox; i++) data[i]=temp[i];
   
    free(temp);
    return(0);
}

int SURF(double* data, int***  region, int* n_avg_region, int nvox, int p){
    double total_weight;
    double value;
    double* temp=calloc(nvox, sizeof(*temp));

    for(int i=0; i < nvox; i++){
        if(n_avg_region[i] > 0) {
            //printf("%d\t%d\n", i, n_avg_region[i]);
            value=avg(region[i], data, n_avg_region[i]);
            //printf("%f\n", value); 
            //value=idw(region[i], data, total_weight_list[i], n_avg_region[i], p);
            //if(value < 0) printf("%f\t%f\n", data[i], value);
        } else value=0; //data[i];
        temp[i]=value;
    }
    for(int i=0; i<nvox; i++) data[i]=temp[i];
   
    free(temp);
    return(0);
}

int is_unique(double* data, int**  region,  int max_ngh, int voxel, double* norm_value){
    double mean=0;
    double sum2=0;
    double sd;
    int k=2;

    for(int i=0; i < max_ngh; i++){
        int index=region[i][0];
        sum2 += data[index] * data[index];
        mean += data[index];
    }     
    mean /= max_ngh;
    sd = sqrt((sum2 - (max_ngh * mean * mean) ) / (max_ngh-1));
    *norm_value= (data[voxel] - mean)/sd;
    //printf("%f %f %f\n",data[voxel], mean, sd); 
    if( data[voxel] > mean + 1*sd || data[voxel] < mean+1*sd) return(1);
    return(0);
}

int local_smoothing(double *data, int** region, int* n_avg_region, int* mask_label, int max_ngh){
    double d, new_mean;
    int n_pos, index, n;
    double value, old_mean;
   
    value=old_mean=new_mean=n_pos=0;
    //calculate mean for region and new value to replace negative one
    for(int j=0; j< max_ngh; j++){
        index=region[j][0];
        old_mean += data[index];
        if(data[index] > 0){ 
            value += data[index];
            n_pos++;
        }
    }
    if(n_pos > 0) { 
        value /= n_pos;
        data[region[0][0]]=value;
        n_pos++; //new positive value
        old_mean /= max_ngh;
        //replace negative value with positive one
        //Get new mean
        for(int j=0; j< max_ngh; j++){ 
            new_mean+=data[region[j][0]];
        }
        new_mean /= max_ngh;
        //calculate difference between old mean and new mean, 
        //divide by number of positive voxels in region 
        d=(old_mean - new_mean)*((float)max_ngh/n_pos);
        //lower positive values so that local mean stays the same
        new_mean=n=0;
        for(int j=0; j< max_ngh; j++){
            index=region[j][0];
            if(data[index] > 0){ 
                data[index] += d;  
                n++ ; 
            } //else printf("N  %d\t%f\t%d\n", j, data[index], index);
            new_mean+=data[index];
            //printf("%f,", data[index]);
        }
       new_mean /= max_ngh;
         //if(  new_mean / old_mean > 2) {printf("%f\t%f\t%f;\t%d %d %d\n", old_mean, d, new_mean, n_pos, n, max_ngh); exit(0);}
    }
    return(0);
}

int eliminate_zeros(double* data, int***  regions, int* n_avg_region, int* mask_label, int nvox, int p){ 
    //, double* noise){
    double total_weight;
    int iterations=0, max_iterations=10;
    int negative_values=1;
    int nlabel=0;
    int found_negative, n;
    double global_average;
    int negatives=1; 
    int iter=0;
    double norm_value;


        printf("Negatives: ");
        while(negatives > 0 && iter < 10){
            for(int i=0; i < nvox; i++){
                if( mask_label[i] > 0){
                    if( data[i] < 0){
                        local_smoothing(data, regions[i], n_avg_region, mask_label,n_avg_region[i] );
                    }
                }   
            }   
            negatives=0;
            for(int i=0; i < nvox; i++) if( data[i] < 0 && mask_label[i] > 0) negatives++;
            printf("%d ", negatives);
            iter++;
        }
        printf("\n");
        global_average=n=0;
        for(int i=0; i < nvox; i++){
            if( mask_label[i] > 0){
                global_average+=data[i];
                n++;
            }
        }
    n=0;
    for(int i=0; i < nvox; i++){
        if(mask_label[i] > 0){ 
            if( data[i] < 0){ 
                n++;    
                data[i]=0;
            }   
        }
    }
    printf("Thesholded %d voxels to 0\n", n  );
    return(nlabel);
}

double calc_error(double* observed, double* blurred, int* mask_labels, int nvox, double* last_error, double* error, double tolerance, int iteration){
    double error_change;

    *error=0;
    for(int i=0; i<nvox; i++){
        if(mask_labels[i] > 0){
            *error += (observed[i]-blurred[i])*(observed[i]-blurred[i]);
        }
    }

    if( iteration > 1){
        //printf("old error: %f ", *last_error); 
        error_change = fabs(*error / *last_error - 1)* 100;
    } else error_change=tolerance+1;

    *last_error = *error;
    return(error_change);
}

struct idargs {
    data* image; 
    char* outputfilename; 
    int*** regions;
    int *n_avg_region; 
    int* mask_label;  
    int max_iterations; 
    double fwhm;
    int p; 
    double tolerance; 
    data* masks; 
    int nmasks; 
    int max_avg; 
    int smooth_only; 
    data* first_guess;
    int nthreads;
    mihandle_t* outputimage;
    float* temp;
    int t;
};


int idSURF_multithread(data* image,mihandle_t* outputimage,  int*** regions, int *n_avg_region, int* mask_label,  int max_iterations, double fwhm, int p, double tolerance, data* masks, int nmasks, int max_avg, int smooth_only, data* first_guess, int nthreads){
    int rc;
    pthread_t threads[nthreads];
    struct idargs thread_args[nthreads];
    for(int t=0; t< nthreads; t++){
        thread_args[t].image=image;
        thread_args[t].regions=regions;
        thread_args[t].n_avg_region=n_avg_region;
        thread_args[t].mask_label=mask_label;
        thread_args[t].max_iterations;
        thread_args[t].fwhm=fwhm;
        thread_args[t].p=p;
        thread_args[t].tolerance=tolerance;
        thread_args[t].masks=masks;
        thread_args[t].nmasks=nmasks;
        thread_args[t].max_avg;
        thread_args[t].smooth_only=smooth_only;
        thread_args[t].first_guess=first_guess;
        thread_args[t].nthreads=nthreads;
        thread_args[t].outputimage=outputimage;


        thread_args[t].t=t;
        rc= pthread_create(&threads[t], NULL, idSURF_threaded, (void *) &thread_args[t] ); //gradient(image, img_dx, img_dy, img_dz);
    }
    for(int t=0; t<nthreads;t++){
        rc = pthread_join(threads[t], NULL);
    }
    return(0);

}

void*  idSURF_threaded(void* args){

    data* image=((struct idargs*) args)->image ; 
    char* outputfilename=((struct idargs*) args)->outputfilename; 
    int*** regions=((struct idargs*) args)->regions;
    int *n_avg_region=((struct idargs*) args)->n_avg_region; 
    int* mask_label=((struct idargs*) args)->mask_label;  
    int max_iterations=((struct idargs*) args)->max_iterations; 
    double fwhm=((struct idargs*) args)->fwhm;
    int p=((struct idargs*) args)->p; 
    double tolerance=((struct idargs*) args)->tolerance; 
    data* masks=((struct idargs*) args)->masks; 
    int max_avg=((struct idargs*) args)->max_avg; 
    int smooth_only=((struct idargs*) args)->smooth_only; 
    data* first_guess=((struct idargs*) args)->first_guess;
    int nthreads=((struct idargs*) args)->nthreads;
    int thread=((struct idargs*) args)->t;
    float *temp=((struct idargs*) args)->temp;
    mihandle_t* outputimage = ((struct idargs*) args)->outputimage;
    idSURF(image, outputimage, outputfilename, regions, n_avg_region, mask_label,  max_iterations, fwhm,  p, tolerance, masks, max_avg,  smooth_only, first_guess, thread, nthreads);
}

void* idSURF(data* image, mihandle_t* outputimage, char* outputfilename, int*** regions, int *n_avg_region,  int* mask_label,  int max_iterations, double fwhm, int p, double tolerance, data* masks, int max_avg,   int smooth_only, data* first_guess, int thread, int nthreads){
    int nvox=image->zmax*image->ymax*image->xmax; //WARNING: Becareful not use image->n (which is the total number of voxels in 4d) as opposed to nvox (which is the number of voxels in 3d)

    char command[300];
    char residual_str[300];
    char test_str[300];
    char outputFilename[50];
    mihandle_t testimage;

    mihandle_t output_observed_image;
    misize_t starts[image->ndim];
    misize_t sizes[image->ndim];
    misize_t* wcount_3d;
    misize_t wstarts[]={0 , 0, 0};
    double* start_3d;
    double* step_3d;
    double max=0.0, min=0.0, obs_min=0.0, obs_max=0.0;
    double* test;
    double* residuals;
    double* noise;
    double* blurred;
    double* observed;
    float* temp;
    double residual;
    int voxel, max_frames, n;
    double energy=0;
    double last_error, error;
    double* image_data=malloc(nvox*sizeof(*image_data));

    float sum=0;
    int nlabels=0;
    double sum_obs=0, sum_blr=0, sum_res=0, sum_ts=0, sum_tr=0, sum_r2=0;
    VERBOSE=1;
    //Allocate memory
    test=malloc(nvox*sizeof(*test));
    residuals=calloc(nvox,sizeof(*residuals));
    noise=calloc(nvox,sizeof(*noise));
    blurred=malloc(nvox*sizeof(*blurred));
    temp=malloc(nvox*sizeof(*temp));


    if(image->tmax==0){
        wcount_3d=image->wcount;
        start_3d=image->start;
        step_3d=image->step;
    }
    else {
        wcount_3d=&(image->wcount[1]);
        start_3d=&(image->start[1]);
        step_3d=&(image->step[1]);
    }
    if (image->tmax == 0) max_frames=1;
    else max_frames= image->tmax;

    for(int i=0; i < nvox; i++) if( nlabels < mask_label[i]) nlabels=mask_label[i];
    
    //writeVolume(  "test.mnc", test, start_3d, step_3d,  wcount_3d, MI_TYPE_DOUBLE );
    //writeVolume("noise.mnc", mask_label , start_3d, step_3d,  wcount_3d, MI_TYPE_INT  );
    for(int t=thread; t < max_frames; t+=nthreads){
        printf("Frame: %d\n", t);
        read_frame(image, t, image_data, MI_TYPE_DOUBLE);

        observed=image_data;
        for(int label=0; label <= nlabels; label++){
            if( masks[label].smooth == TRUE  ){
                printf("Smoothing label:%d\n", label+1);
                for(int k=0; k< 3; k++) SURF2( observed,  regions, n_avg_region,  nvox, mask_label, label+1);
            }
        }
        eliminate_zeros(observed, regions, n_avg_region, mask_label,  nvox, p); //, noise);
        
        double initial_energy=0, old_max=0;
        if(first_guess->filename != NULL) read_frame(first_guess, t,  image_data, MI_TYPE_DOUBLE); 
        else first_guess->data=observed;
        copy_with_threshold(test,first_guess->data, mask_label, nvox, &initial_energy, &old_max, 0); 

        
        /********
        *ITERATE*
        *********/
        for( int iteration=1; iteration <= max_iterations; iteration++  ){
            printf("%d:\t", iteration);
            if(smooth_only != 1){
                
                sum=n=0; 
                for(int t=0; t < nvox; t++) {
                    if(mask_label[t] == 3) { 
                        sum += test[t];
                        n++;
                    }
                }
                /*************
                *PSF Modeling*
                **************/
                float sigma=fwhm/2.35482005;
                int numsteps=10;
                for(int k=0; k<nvox; k++) temp[k]=(float) test[k]; 
                //sprintf(test_str, "test_%d_%d.mnc", t, iteration);
                //writeVolume(test_str, test, start_3d, step_3d,  wcount_3d, MI_TYPE_DOUBLE  );
                gaussianiir3d(temp, image->xmax, image->ymax, image->zmax, sigma, numsteps);
                //sprintf(test_str, "test_blur_%d_%d.mnc", t, iteration);
                //writeVolume(test_str, temp, start_3d, step_3d,  wcount_3d, MI_TYPE_FLOAT  );
                //sprintf(command, "mincblur -clobber -quiet -3dfwhm %1.2f %1.2f %1.2f %s test", fwhm, fwhm, fwhm, "test.mnc");
                //printf("%s\n", command);
                //system(command);

                //miopen_volume("test_blur.mnc", MI2_OPEN_READ, &testimage);
                //miget_real_value_hyperslab(testimage, MI_TYPE_FLOAT, tstarts, wcount_3d, temp);
                //miclose_volume(testimage);

                for(int k=0; k<nvox; k++) blurred[k]=(double) temp[k];
                /********************
                *Residual and Update* 
                *********************/
                int disp_count=0;
                sum_r2=0;
                float residual_avg =0; int n=0;
                for( voxel=0; voxel < nvox; voxel++){
                    residuals[voxel]=0;
                    if(mask_label[voxel] > 0){
                        if(blurred[voxel] == 0) residual=1; 
                        else residual=observed[voxel]/blurred[voxel];

                        //if(residual > 3) residual=3; //make sure residual is not too large
                        if(residual <= 0) residual=1; //{ printf("zero: mask_label %d\tresidual %f\tobserved %f\tblurred %f\ttest %f\n",mask_label[voxel], residual, observed[voxel], blurred[voxel], test[voxel]); residual=1;}
                        test[voxel] *= residual; //update image
                        residuals[voxel]=residual;
                    }
                }
                //SURF(residuals, regions, n_avg_region, region_t_weights, nvox, p);
                //sprintf(residual_str, "residuals_%d_%d.mnc", t, iteration);
                //writeVolume(residual_str, residuals, start_3d, step_3d,  wcount_3d, MI_TYPE_DOUBLE);
            } else { printf("Smoothing only\n");}

            n=0; 
            sum_tr=0;
            for(int t=0; t < nvox; t++) {
                if(mask_label[t] == 3) { 
                    sum_tr += test[t];
                    n++;
                }
            }

            /******************
            *Noise attenuation*
            *******************/
            SURF( test,  regions, n_avg_region, nvox, p);
            double error_change=calc_error(observed, blurred, mask_label, nvox, &last_error, &error,tolerance, iteration);
 
            n=0; 
            sum_ts=0;
            sum_obs=0;
            sum_blr=0;
            sum_res=0;
            for(int t=0; t < nvox; t++) {
                if(mask_label[t] == 3) { 
                    sum_ts += test[t];
                    sum_obs += observed[t];
                    sum_blr += blurred[t];
                    sum_res += residuals[t];
                    n++;
                }
            }
           //printf("t: %3.2f, tr: %3.2f, ts: %3.2f o: %3.2f, b: %3.2f, r:%3.2f, n: %d\t", sum/n, sum_tr/n, sum_ts/n, sum_obs/n, sum_blr/n, sum_res/n, n );
           printf("\n");
            if( error_change < tolerance && smooth_only != 1) break;
        }  
        /*********************
        *Write output to file*
        **********************/
        if(image->ndim == 4){
            starts[0]=t;
            starts[1]=starts[2]=starts[3]=0;
            sizes[0]=1;
            sizes[1]=image->zmax;
            sizes[2]=image->ymax;
            sizes[3]=image->xmax;
            printf("%llu %llu %llu %llu\n", starts[0], starts[1], starts[2], starts[3]);
            printf("%llu %llu %llu %llu\n", sizes[0], sizes[1], sizes[2], sizes[3]);
        } else {
            starts[0]=starts[1]=starts[2]=0;
            sizes[0]=image->zmax;
            sizes[1]=image->ymax;
            sizes[2]=image->xmax;
            printf("%llu %llu %llu\n", starts[0], starts[1], starts[2]);
            printf("%llu %llu %llu\n", sizes[0], sizes[1], sizes[2]);
        }
        energy=n=0;
        for (voxel=0; voxel < nvox; voxel++) {
            if(mask_label[voxel] > 0 ){
                if (max < test[voxel]) max=test[voxel];
                if (min > test[voxel]) min=test[voxel];
                if (obs_max < test[voxel]) obs_max=test[voxel];
                if(obs_min > test[voxel]) obs_min=test[voxel];
                if(observed[voxel] > 0 && test[voxel] > 0){
                    energy += test[voxel] / observed[voxel];
                    n++;
                }
            }
        }

        //printf(" Change: Energy = %3.3f, Max = %3.3f\n", energy/n, max/old_max);
        fflush(stdout);
        printf("%f %f\n", max, min);
        miset_volume_range(*outputimage, max, min);
  
        miset_real_value_hyperslab(*outputimage, MI_TYPE_DOUBLE, starts, sizes, test);
        miset_volume_range(output_observed_image, obs_max, obs_min);

    }
    free(image_data);
    free(test);
    free(residuals);
    free(noise);
    free(blurred);
    free(temp);
}

int*  get_nlabels(data* mask, int* mask_array, int* nlabels){
    int* labels=NULL; //labels stores the integer values of atlas
    for( int z=0;  z < mask->zmax; z++ ) {
        for( int y=0; y < mask->ymax; y++) {
            for( int x=0; x < mask->xmax; x++){
                int index=z*mask->ymax*mask->xmax + y*mask->xmax+x;
                if( mask_array[index] != 0 ){
                    int label = mask_array[index];
                    if ( isin2(label, labels, *nlabels) != TRUE){ //label has already been found
                        (*nlabels) += 1;
                        labels=realloc(labels, *nlabels * sizeof(*labels));
                        labels[*nlabels-1]=label;
                    }  //first time coming across this label
                    //int idx=find(label, labels, *nlabels);
                }
            }
        }
    }
    return(labels);
}




int main(int argc, char** argv){      
    if(argc == 1 || strcmp(argv[1], "-help") == 0 ) useage();
    int argi=1;
    char* outputfilename;
    int nimages=0;
    int nmasks=0;
    int nsurfaces=0;
    int nframes=0;
    int spatial_dim=0;
    int nvox3d;
    int p=1;
    double* avg_region_weights;
    float fwhm, tolerance;
    int** near_ngh; //nearest neighbours for each voxel in the masks
    int* n_near_ngh; //number of nearest neighbours for each voxel in the masks 
    int*** avg_region; //region of voxels over which to average 
    int *n_avg_region; //number of voxels in averaging region. maximum value, i.e. max_avg, is set by user
    int *mask_label;//keeps track of which mask to use for the nearest neighbours and averaging region
    int max_iterations; //maximum number of iterations for iterative deconvolution with idSURF
    int max_avg; //maximum number of voxels over which to average.
    int nthreads=1; //sysconf(_SC_NPROCESSORS_ONLN);
    int smooth_only=0;
    int cubic_averaging_region=FALSE;
    data** temp=parse_input(argc, argv, &nthreads, &nimages, &nmasks, &max_iterations, &max_avg, &fwhm, &tolerance,  &outputfilename, &smooth_only, &cubic_averaging_region);
    data* image=temp[0];
    data* masks=temp[1];
    data* first_guess=temp[2];
    mihandle_t img;
    mihandle_t outputimage;

    readVolume(image, 0, MI_TYPE_DOUBLE );

    if(first_guess->filename != NULL) first_guess->data=(double*) readVolume(first_guess, 1, MI_TYPE_DOUBLE);
    nvox3d=image->zmax*image->ymax*image->xmax;


    //Create output volume
    createVolume(outputfilename, image->ndim, image->wcount, image->step, image->start, MI_TYPE_DOUBLE);
    miopen_volume(outputfilename, MI2_OPEN_RDWR, &outputimage);

    //Allocate memory for neighbours and averaging region
    near_ngh=malloc(sizeof(*near_ngh)* nvox3d);
    n_near_ngh=calloc(nvox3d, sizeof(*n_near_ngh));
    avg_region=malloc(sizeof(*avg_region)*nvox3d);
    n_avg_region=calloc(nvox3d, sizeof(*n_avg_region));
    mask_label=malloc(sizeof(*mask_label) * nvox3d );
    avg_region_weights=malloc(sizeof(*avg_region_weights) * nvox3d);
    for(int i=0; i<nvox3d;i++) mask_label[i]=-1;
    for(int i=0; i<nvox3d; i++) avg_region[i]=NULL; //Initialize values of region to NULL so that we can use realloc on them




    for(int i=0; i<nmasks; i++){
        //read in mask volume 
        masks[i].data=(double*) readVolume(&masks[i], 1, MI_TYPE_DOUBLE );
        if( masks[i].data == NULL  ) pexit("Could not read data from" , masks[i].filename,  1);
        printf("Finding neighbours...");
        find_neighbours(max_avg, near_ngh, n_near_ngh, mask_label, i, &masks[i]);
        printf("Done.\n");
        printf("Cubic Averaging Regiong: %d\n", cubic_averaging_region);
        if( cubic_averaging_region == TRUE){ 
            printf("Finding cubic averaging region...");
            find_cubic_averaging_region(image, avg_region, n_avg_region, masks[i].data,  max_avg);
            printf("Done.\n");
        }
        else{
            //For voxels in masks, find nearest (i.e., adjacent) neighbours
            //find averaging region for each voxel that is within a mask
            printf("Finding averaging region...");fflush(stdout);
            find_averaging_region(avg_region, n_avg_region, mask_label, near_ngh, n_near_ngh, nvox3d, masks[i].data, max_avg );
            printf("Done.\n");
            //find distance weights for averaging region
            //free(masks[i].data);
        }
        printf("Getting region weights..."); fflush(stdout);
        get_total_weight(avg_region_weights, avg_region, nvox3d,  p, n_avg_region);
        printf("Done.\n");
    }

    printf("Performing iterative deconvolution.\n");
    if( nthreads > 1) idSURF_multithread(image,  &outputimage, avg_region, n_avg_region, mask_label, max_iterations, fwhm,  p, tolerance, masks, nmasks, max_avg, smooth_only, first_guess, nthreads);
    else{
        idSURF(image, &outputimage, outputfilename, avg_region, n_avg_region, mask_label,  max_iterations, fwhm,  p, tolerance, masks, max_avg,  smooth_only, first_guess, 0, 1);
    }

    miclose_volume(outputimage);
    return 0;
}
  




void useage(){
        printf("\nName:\n\tidSURF-4.2.1 -- partial-volume correct PET image using masks of anatomical regions\n");
        printf("Description:\n");
        printf("\tPerform partial-volume correction PET image using labeled mask of anatomical regions\n");
        printf("Useage:\n");
        printf("\t-pet <.mnc>\t\tTarget volumetric MINC2 file.\n");
        printf("\t-mask <.mnc>\t\tVolumetric MINC2 file with ROI labels.\n");
        printf("\t-o <.mnc>\t\tPVC output file.\n");
        printf("\t-max-iterations \tMaximum number of iterations per frame (default, set to 10).\n");
        printf("\t-nvoxel-to-average \tNumber of voxels over which to smooth (default, set to 64).\n");
        printf("\t-tolerance \t\tTolerance for deconvolution (default, set to 0.01).\n");
        printf("\t-fwhm\t\t\tfull-width at half-maximum of PET scanner\n");
        printf("Optional:\n");
        printf("\t-cubic_averaging_region\tUse a cubic averaging region\n");
        printf("\t-smooth_only\t\tNo PVC, only non-linear smoothing\n");
        printf("\t-first_guess\t\tManually choose an image to serve as first guess for deconvolution\n");
        printf("Example:\n\tidSURF -fwhm 2.5 -tolerance 0.01 -nvoxel-to-average 64 -pet pet_image.mnc -mask classify.mnc -m -o pvc.mnc\n");
        exit(1);
}





data** parse_input(int argc, char** argv, int* nthreads,  int* nimages, int *nmasks, int* max_iterations, int* max_avg, float* fwhm, float* tolerance, char** outputfilename, int* smooth_only, int *cubic_averaging_region ){
    int i, temp_i=0;
    char* linear="-linear";
    char* nearest="-nearest";
    int set_pet=0, set_mask=0,set_max_iterations=0, set_nvoxel_avg=0, set_fwhm=0, set_tolerance=0, set_outputfilename=0, set_smooth_only=0, set_first_guess=0;
    data** output=malloc(sizeof(*output) * 3);
    for( i=1; i < argc; i++ )
        if(strcmp(argv[i], "-mask" ))
            *nmasks += 1;

    data *masks=malloc(sizeof(data) * *nmasks);
    data *images=malloc(sizeof(data));
    data *first_guess=malloc(sizeof(data)); 
    output[0]=images;
    output[1]=masks;
    output[2]=first_guess;
    first_guess->filename=NULL;
    *nmasks=0;
    
    for( i=1; i < argc; i++ ){
        if(strcmp(argv[i], "-pet" ) == 0 ){
            if( local_allocate_vol(images, nimages, argv[i+1]) != 0 ) ; //FIXME, shouldn't have a useless conditional  
            i++;
            images[0].filename=argv[i];
            printf("PET: %s\n", images[0].filename);
            set_pet=1;
        }
        else if(strcmp(argv[i], "-mask" ) ==0 ){
            if( local_allocate_vol(masks, nmasks, argv[i+1]) != 0 ) ;//FIXME, shouldn't have a useless conditional 
            if( masks[*nmasks-1].smooth == TRUE) i++; 
            i++;

            masks[*nmasks-1].filename=argv[i];
            printf("Mask: %s\n",  masks[*nmasks-1].filename );
            set_mask=1; 
        }
        else if(strcmp(argv[i], "-fwhm" )==0 ){
            i++;
            *fwhm=atof(argv[i]); 
            set_fwhm=1; 
        }
        else if(strcmp(argv[i], "-nthreads" )==0 ){
            i++;
            *nthreads=atoi(argv[i]); 
            //Not sure if this is needed:
            //set_nthreads=1; 
        }
        else if(strcmp(argv[i], "-max-iterations" )==0 ){
            i++;
            *max_iterations=atoi(argv[i]); 
            set_max_iterations=1;
        }
        else if(strcmp(argv[i], "-nvoxel-to-average" )==0 ){
            i++;
            *max_avg=atoi(argv[i]); 
            set_nvoxel_avg=1; 
        }
        else if(strcmp(argv[i], "-tolerance" )==0 ) {
            i++;
            *tolerance=atof(argv[i]);
            set_tolerance=1;
        }
        else if(strcmp(argv[i], "-smooth_only" )==0 ) {
            *smooth_only=1;
            set_smooth_only=1;
        }
        else if(strcmp(argv[i], "-cubic_averaging_region" )==0 ) {
            *cubic_averaging_region=1;
        }
        else if(strcmp(argv[i], "-first_guess" )==0 ) {
            first_guess=malloc(sizeof(data));
            if( local_allocate_vol(first_guess, &temp_i, argv[i+1]) != 0 ) ;
            i++;
            first_guess->filename=argv[i];
            output[2]=first_guess;
        }
        else if(strcmp(argv[i], "-o" )==0 ) {
            i++;
            *outputfilename=argv[i];
            set_outputfilename=1;
        }
    }
    if( set_smooth_only != 1) {
        if(set_pet != 1 ||set_mask != 1 || set_fwhm != 1 || set_max_iterations != 1 || set_nvoxel_avg != 1 || set_fwhm != 1 || set_tolerance != 1 || set_outputfilename != 1 ){
            useage();
        }
    } else if( set_smooth_only == 1){
        if(set_pet != 1 ||set_mask != 1 || set_max_iterations != 1 || set_nvoxel_avg != 1 || set_outputfilename != 1){          printf("Got here!\n");
            useage();
        }
    }
    return(output);
}

int local_allocate_vol(data* ptr, int* n, char* argv_next){
    *n += 1;
   
    ptr[*n-1].label=*n-1; 
    ptr[*n-1].ngroups=0;
    ptr[*n-1].tmax=ptr[*n-1].zmax=ptr[*n-1].ymax=ptr[*n-1].xmax=0;
    ptr[*n-1].n=0;
    ptr[*n-1].ngroups=0;
    //Initialize pointers in vol to NULL.
    ptr[*n-1].groups=NULL;
    ptr[*n-1].filename=NULL;
    ptr[*n-1].count=NULL;
    ptr[*n-1].start=NULL;
    ptr[*n-1].step=NULL;
    ptr[*n-1].data=NULL;
    ptr[*n-1].wcount=NULL;
    ptr[*n-1].wstarts=NULL;
    if(strcmp(argv_next, "-smooth") == 0) {
        ptr[*n-1].smooth=TRUE;
    }
    else ptr[*n-1].smooth=FALSE; 

    return(0);
}



