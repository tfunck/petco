#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include "minc2.h"
#include "hdf5.h"
#include "volume_io.h"
#include "minc_helper.h"
#include "3x3-C/dsyevh3.h"
#include "3x3-C/dsyevd3.h"

void printMat(float *mat){
	for(int y=0; y<3; y++ ){
		for(int x=0; x< 3; x++){
			printf("%f\t", mat[y*3+x]);
		}
	printf("\n");
	}	
printf("\n");
}
void printMat_d(double *mat){
	for(int y=0; y<3; y++ ){
		for(int x=0; x< 3; x++){
			printf("%f\t", mat[y*3+x]);
		}
	printf("\n");
	}	
printf("\n");
}
int isin(int number, int* array, int dim){

    for(int i=0; i<dim; i++) if( number==array[i]) return 1;

return 0;

}


double jacobi_iteration(data* input, data* mask){
    int zmax=input->count[0];
    int ymax=input->count[1];
    int xmax=input->count[2];
    double dz=input->step[0];
    double dy=input->step[1];
    double dx=input->step[2];
    double update=0;
    int indexp, indexm, index;
    int p=1, m=1;
    double average_residual=0;
    int n=0;
    if(dz!=dy ||dz!=dx || dx!=dy) printf("Error: Assumes all voxel dimensions are of equal size.\n");
    for(int z=1; z<zmax-1; z++){
        for(int y=1; y < ymax-1; y++){
            for(int x=1; x < xmax-1; x++){
                index=z*ymax*xmax+y*xmax+x;
                if( mask->data[index]==1){
                    update=0;
                    //if(z-1 < 0 || z+1 >= zmax || y-1 < 0 || y+1 >= ymax || x-1 < 0 || x+1 >= xmax ){
                    //    continue;     
                    // }
                    //else{
                        //dx
                        indexp=z*ymax*xmax+y*xmax+x+p;
                        indexm=z*ymax*xmax+y*xmax+x-m;
                        update =  input->data[indexp] + input->data[indexm];
                        //dy
                        indexp=z*ymax*xmax+(y+p)*xmax+x;
                        indexm=z*ymax*xmax+(y-m)*xmax+x;
                        update += input->data[indexp] + input->data[indexm]; //gradient in y direction
                        //dz
                        indexp=(z+p)*ymax*xmax+y*xmax+x;
                        indexm=(z-m)*ymax*xmax+y*xmax+x;
                        update += (float) input->data[indexp] + input->data[indexm]; //gradient in z direction
                        //finalize update
                        update = update /6;

                        average_residual += fabsf(input->data[index]- update);
                        input->data[index] = update;
                        n++;
                   // }
                }
            }
        }
    }
    return( (double )average_residual/n );
}




int solve(data* mask, data* image, double tolerance, char* outputfilename){
    double residual=1;
    int iterations=0; 
    while(residual>tolerance){
        //Calculate gradient of temporary image
        //EQ: Grad( Image_i )
        //Can turn this into function, nicer to read.
        residual = jacobi_iteration(image, mask);
        iterations++;
        printf("%d:\t%f\n", iterations, residual);
    }
    
    writeVolume(outputfilename, image->data, image->start,  image->step, image->wcount, MI_TYPE_DOUBLE );
    
    return(0);
}

int main(int argc, char** argv){
    if(argc != 5){ 
        printf("\nUseage: tolerance mask.mnc image.mnc output.mnc\n");
        printf("Details This program performs surface-based aniosotropic diffusion filtering on a target image using surface normals.\n");
        exit(1);
    }
    
    int argi=1;
    data mask, image;
    double tolerance=atof(argv[argi++]);
    mask.filename=argv[argi++];
    image.filename=argv[argi++];
    char* outputfilename=argv[argi++];

    mask.data=(double*) readVolume(&mask, 1, MI_TYPE_DOUBLE);
    image.data=(double*)readVolume(&image, 1, MI_TYPE_DOUBLE);
    if (checkDimensions(&mask, &image) == TRUE ) pexit("Mask and image dimensions are not the same.", "", 1);
    if( solve(&mask, &image,  tolerance, outputfilename) == 0) {
        printf("Completed.\n"); 
        return(0);
    }
    else{
        printf("Failed.\n");
        return(1);
    }
    return 0;
}
 


