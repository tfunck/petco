#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>
#include "minc2.h"
//#include "gaussianiir3d.h"
#include "minc_helper.h"
#include "minc_fftw-1.0.h"

void gaussj(float **a, int n, float **b, int m);
int local_allocate_vol(data* ptr, int* n);
void useage();
data** parse_input(int argc, char** argv,  int* nimages, int *nmasks, double* fwhm, char** outputfilename );
int local_allocate_vol(data* ptr, int* n);
struct linfo {
    int nlabels;
    int *labels;
    int *nLocalVoxels;
};

int find(int x, int* list, int n){
    //find first element in list that equals x
	//and return the index number of that element
    for(int i=0; i< n; i++){ 
        if( x == list[i] ){ 
            return(i);
        }
    }
    return(-1);
}



int isin(int x, int* list, int n){
    for(int i=0; i< n; i++){ 
        if( x == list[i] ){ 
            return(1);
        }
    }
    return(0);
}


struct linfo* get_nlabels(data* mask, int* mask_array){
    int nlabels=0;
    struct linfo* label_info=malloc(sizeof(*label_info));
    int* labels=NULL; //labels stores the integer values of atlas
    int* nLocalVoxels=NULL;
    for( int z=0;  z < mask->zmax; z++ ) {
        for( int y=0; y < mask->ymax; y++) {
            for( int x=0; x < mask->xmax; x++){
                int index=z*mask->ymax*mask->xmax + y*mask->xmax+x;
                if( round(mask_array[index]) > 0 ){
                    int label = round(mask_array[index]);
                    if ( isin(label, labels, nlabels) != TRUE){ //label has already been found
                        nlabels++;
                        labels=realloc(labels, nlabels * sizeof(*labels));
                        nLocalVoxels=realloc(nLocalVoxels, nlabels * sizeof(*nLocalVoxels));
                        labels[nlabels-1]=label;
                        nLocalVoxels[nlabels-1]=0;
                    }  //first time coming across this label
                    int idx=find(label, labels, nlabels);
                    nLocalVoxels[idx]++;
                }
            }
        }
    }
    for(int i=0; i<nlabels; i++) printf("Label: %d %d\n", i, labels[i]);
    label_info->labels=labels;
    label_info->nlabels=nlabels;
    label_info->nLocalVoxels=nLocalVoxels;
    return(label_info);
}


int gtm_solve(data* mask, int* mask_array, int nlabels, float** gtm, double fwhm[3], int* labels, int* nLocalVoxels){
    double* rsfTotal;
    int i, j, k, n=0;
    mihandle_t image;
    int* mask1, *mask2;
    float* rsf;
    float** zero_array; //Used for inverting GTM with gaussj 
    int nvox = mask->n;
    int label;
    double* temp=malloc(nvox*sizeof(*temp)); //
    int x;
    
    //Allocate array for RSF 
    rsf=calloc( nvox, sizeof(*rsf));
    //Allocate zero_array that keeps track of GTM for gaussj
    zero_array=malloc(sizeof(*zero_array) * nlabels);
    for(i=0; i< nlabels; i++){ 
        zero_array[i]=calloc(nlabels,  sizeof(**zero_array));
    }
    
    //Calculate the regional response function (RSF)
    printf("Label\tNVoxels\tRSF\n");
    for(i=0; i < nlabels; i++){
        label=labels[i];
        for(int index=0;  index < nvox; index++ ){
            if( mask_array[index] == label ){
                temp[index]=1;
            }
        }
		anisotropic_gaussian_filter(temp, fwhm, mask->count, mask->step); 
        //gaussian_filter(temp, fwhm[0], mask->count, mask->step); 
        for(int k=0; k<nvox; k++){ 
            rsf[k]=temp[k];
            temp[k]=0;
        }
        //find the proportion of each blurred mask i in mask j
        //cycle through each mask j...
        for(j=0; j < nlabels; j++){
            int label2 = labels[j];
            //gtm[j][i]=0; //initialize to 0, just in case it was automatically set to another value
            for(k=0; k < nvox; k++){
                //Find the amount of rsf of mask i in mask j
                if( mask_array[k] == label2){ 
                    //for this k we are in our mask i, so we add up the values of mask j at k
                    //printf("%d %f\n", mask_array[k], rsf[k]);
                    gtm[j][i] += (float) rsf[k]; //this voxel is in region, so we add rsf j to sum
                }
            }
        }
    }

    for(i=0; i < nlabels; i++){
        for(j=0; j < nlabels; j++){
            if(nLocalVoxels[j] > 0) gtm[j][i] /=  nLocalVoxels[j]; //Divide the total signal in mask j 
            else gtm[j][i]=0;
            printf("%lf\t", gtm[j][i]);
        } printf("\n");
    }
    printf("\n");
    //for(j=0; j < nlabels; j++)  printf("%d\t%d\n", labels[j], nLocalVoxels[j] );
    //invert GTM
    gaussj(gtm, nlabels, zero_array, 1);
    for(i=0; i < nlabels; i++){
        for(j=0; j < nlabels; j++){
            printf("%lf\t", gtm[i][j]);
        } printf("\n");
    } printf("\n");


    for(j=0; j < nlabels; j++){ 
        free(zero_array[j]);
    }
    free(rsf);
    free(zero_array);
    free(temp); 
return 0;
}

int gtm_apply(data* image, data* masks, int* mask_array, int* labels, int* nLocalVoxels, int nlabels, char* outputfilename, float** gtm){
    int i,j,  k, index;
    int nVoxels=masks[0].n;
    int max_frames;
    float*** gtmPVC_image;
    double max=0.0, min=0.0;
    float obsValues[nlabels];
    float trueValues[nlabels];
    float* array=calloc(nVoxels, sizeof(*array));
    misize_t *sizes=calloc(image->ndim, sizeof(*sizes));
    misize_t *starts=calloc(image->ndim , sizeof(*starts));
    mihandle_t outputimage;
    miopen_volume(outputfilename, MI2_OPEN_RDWR, &outputimage);


    if (image->tmax == 0) max_frames=1;
    else max_frames= image->tmax;
    float *temp=malloc(nVoxels*sizeof(*temp));
    if(image->data == NULL) image->data=malloc(nVoxels*sizeof(double));
    for(int t=0; t < max_frames; t++){
        read_frame(image, t, image->data, MI_TYPE_DOUBLE);
        for(i=0; i< nlabels; i++) obsValues[i]=0;
        //Initialize the average value in this region and the number of voxels found in this mask region to 0

        //Cycle through the entire image
        for(k=0; k < nVoxels; k++){
            //if the point k is in mask i then...
            if(  mask_array[k] > 0){
                int val = mask_array[k];
                int i = find(val, labels, nlabels);
                if( i < 0 ) pexit("Error: label list does not contain value", "", 1);
                //add the value of the pet image at k to the observed value
                //printf("%d %d %f\n", i, val, image->data[k]); 
                obsValues[i] += image->data[k];
                //and increment the number of voxels by 1
            }
        }

        //In order to prevent dividing by 0 in the case where a mask is empty
        //check that there are voxels in the mask before taking the average of the 
        //observed values within the mask.
        for(i=0; i< nlabels; i++){ 
            if (nLocalVoxels > 1) obsValues[i] /= nLocalVoxels[i];
            else obsValues[i]=0;
        }
        printf("%d\tObserved\n",t );
        for(int x=0; x <nlabels; x++) printf("%d\t%d\t%f\n", labels[x], x, obsValues[x]);
        //Multiply observed values by inverted GTM
        for(int x=0; x< nlabels; x++){ 
            trueValues[x]=0; 
            for(int y=0; y < nlabels; y++ ){ 
                trueValues[x] += gtm[x][y] * obsValues[y];
            }
        }  
        printf("PVC\n");
        for(int x=0; x <nlabels; x++) printf("%d\t%d\t%f\n", labels[x], x, trueValues[x]);
        printf("\n");
        for(index=0; index < nVoxels; index++){
            int val = mask_array[index];
            if( val > 0){
                int idx=find(val, labels , nlabels);
                array[index]=(float) trueValues[idx];     
            }
        }
        if(image->ndim == 4){
            starts[0]=t;
            starts[1]=starts[2]=starts[3]=0;
            sizes[0]=1;
            sizes[1]=image->zmax;
            sizes[2]=image->ymax;
            sizes[3]=image->xmax;
        } else {
            starts[0]=starts[1]=starts[2]=0;
            sizes[0]=image->zmax;
            sizes[1]=image->ymax;
            sizes[2]=image->xmax;
        }
        if(t==0){min=max=(double) array[0];}
        for (int l=0; l < nlabels; l++) {
            if (max <(double) trueValues[l]) max=(double) trueValues[l];
            if (min >(double) trueValues[l]) min=(double) trueValues[l];
        }
        //for(int a=0; a < 4; a++) printf("%d ", sizes[a] ); printf("\n");
        //for(int a=0; a < 4; a++) printf("%d ", starts[a] ); printf("\n");

        miset_real_value_hyperslab(outputimage, MI_TYPE_FLOAT, starts, sizes, array);
    }
    //printf("Min: %f\tMax: %f\n", min, max);
    miset_volume_range(outputimage, max, min);
    miclose_volume(outputimage);

    free(sizes);
    free(starts);
    free(array);
    return 0;
}


void useage(){
    printf("gtm -z <z fwhm> -y <y fwhm> -x <x fwhm> -pet pet.mnc -mask mask1.mnc -o output.mnc\n");
    exit(1);
}


/*
   The equation for the weights is: 

   w_i,j = (1 / n_vox) integral_ROI_j {  RSF_i(r)  } dr
   
   where, n_pix is equal to the number of voxels in ROI_j. The weighting factors w_i,j represent the 
   contribution of each domain Di to any ROI_j, i.e.,  how much of the signal from ROI_i ends up in ROI_j. 

   | obs_1 |     | w_1,1   w_2,1   ... w_N,1  |   | True_1 |
   | obs_2 |     | w_1,2   w_2,2   ... w_N,2  |   | True_2 |
   | .     |  =  | .                   .      | X | .      |
   | .     |     | .                   .      |   | .      |
   | .     |     | .                   .      |   | .      |
   | obs_n |     | w_1,N   w_2,N   ... w_N,N  |   | True_N |

*/

/*
 * NOTES:
 * 3.2 - 	Doesn't compile against old minc library. Must be libminc from the github
 * 			repository.
 * */

int main(int argc, char** argv){
    if(argc == 1) useage();
    else if( argc > 1 && strcmp(argv[1], "-help")==0) useage();
    int max, inc;
    char* gtmFilename;
    char* petfilename;
    char* outputfilename;
    int i,j,k, n, index;
    int nVoxels;
    int nimages=0;
    int nlabels=0;
    double fwhm[3];
    float **gtm;
    int* mask_array;
    float* temp_array;
    int* nLocalVoxels; 
    int* labels;    
    int nmasks=0;
    data** temp=parse_input(argc, argv, &nimages, &nmasks, fwhm, &outputfilename );
    data* image=temp[0];
    data* masks=temp[1];

    image->data=(double*) readVolume(image, 1, MI_TYPE_DOUBLE ); //Load PET image
    temp_array=(float*) readVolume(masks, 1, MI_TYPE_FLOAT ); //Load mask image
    mask_array=malloc(sizeof(*mask_array)*image->n);
    for(int i=0; i<image->n; i++) mask_array[i]=round(temp_array[i]);
    free(temp_array);

    createVolume(outputfilename, image->ndim, image->wcount, image->step, image->start, MI_TYPE_FLOAT); //Create output MINC file

    struct linfo* label_info=get_nlabels(masks, mask_array); //Figure out how many and what labels are in the mask image
    nlabels=label_info->nlabels; //Number of labels
    nLocalVoxels=label_info->nLocalVoxels; //Number of voxels with this label
    labels=label_info->labels; //List of the labels in the mask file
   
    //Allocate memory for GTM matrix
    gtm=malloc(nlabels * sizeof(*gtm));
    for(i=0; i < nlabels; i++){ 
        gtm[i]=calloc( nlabels, sizeof(**gtm));
    }

    //Get rsf from each label and solve the GTM
    gtm_solve(masks, mask_array, nlabels,  gtm, fwhm, labels, nLocalVoxels);
    //Apply gtm to each frame
    gtm_apply(image, masks, mask_array, labels, nLocalVoxels, nlabels, outputfilename, gtm);

    for(i=0; i < nlabels; i++){ 
        free(gtm[i]);
    }
    free(mask_array);
    free(labels);
    free(nLocalVoxels);
    free(image->data);  

    return 0;
}


data** parse_input(int argc, char** argv,  int* nimages, int *nmasks, double fwhm[3], char** outputfilename ){
    int i;
    char* linear="-linear";
    char* nearest="-nearest";
    int set_pet=0, set_mask=0,set_max_iterations=0, set_nvoxel_avg=0, set_fwhm=0, set_tolerance=0, set_outputfilename=0;
    data** output=malloc(sizeof(*output) * 2);
    for( i=1; i < argc; i++ )
        if(strcmp(argv[i], "-mask" ))
            *nmasks += 1;

    data *masks=malloc(sizeof(data) * *nmasks);
    data *images=malloc(sizeof(data));
    output[0]=images;
    output[1]=masks;
    *nmasks=0;
    
    for( i=1; i < argc; i++ ){
        if(strcmp(argv[i], "-pet" ) == 0 ){
            if( local_allocate_vol(images, nimages) != 0 ) ; //FIXME, shouldn't have a useless conditional  
            i++;
            images[0].filename=argv[i];
            printf("PET: %s\n", images[0].filename);
            set_pet=1;
        }
        else if(strcmp(argv[i], "-mask" ) ==0 ){
            if( local_allocate_vol(masks, nmasks) != 0 ) ;//FIXME, shouldn't have a useless conditional  
            i++;
            masks[*nmasks-1].filename=argv[i];
            printf("Mask: %s\n",  masks[*nmasks-1].filename );
            set_mask=1; 
        }
        else if(strcmp(argv[i], "-z" )==0 ){
            i++;
            fwhm[0]=atof(argv[i]); 
            set_fwhm+=1; 
        } 
 		else if(strcmp(argv[i], "-y" )==0 ){
            i++;
            fwhm[1]=atof(argv[i]); 
            set_fwhm+=1; 
        }
        else if(strcmp(argv[i], "-x" )==0 ){
            i++;
            fwhm[2]=atof(argv[i]); 
            set_fwhm+=1; 
        }
        else if(strcmp(argv[i], "-o" )==0 ){
            i++;
            *outputfilename=argv[i]; 
            set_outputfilename=1; 
        }
    }
    if(set_pet != 1 || set_mask != 1 || set_fwhm != 3 || set_outputfilename != 1 ){
        useage();
    }
return(output);
}

int local_allocate_vol(data* ptr, int* n){
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

    return(0);
}

