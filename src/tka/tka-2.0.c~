#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <malloc.h>
#include "llsqwt.h"
#include "hdf5.h"
#include "minc2.h"
#include "fft.h"
#include "nrutil.h"
#include "mtga.h"
#include "minc_helper.h"
#include "pthread.h"
#include <unistd.h>

void useage();
int patlak_plot_multithread(data* pet, data* output, int* ref_mask_array, int*  brain_mask_array, double* integral_table, double* reference_tac,double *reference_tac_integral, float* time_widths, parameters*  user_options, int nthreads);
int logan_plot_multithread(data* pet,data* output, int* ref_mask_array, int*  brain_mask_array,double* integral_table, double* reference_tac,double *reference_tac_integral, float* time_widths,parameters*  user_options, int nthreads);
static double calcLinearRegression(int n, double* x, double* y);
int parse_input(int argc, char** argv,  parameters* user_options);
static int regression_choice(double *x, double *y, int nr, double *slope, double *ic, double *ssd, parameters* user_options);
int logan_plot(void* args);
int patlak_plot(void* args);


struct mtga_args {
    data* pet;
    data* output;
    int*  ref_mask_array;
    int* brain_mask_array;
    float* time_widths;
    parameters* user_options;
    double* integral_table;
    double* reference_tac;
    double* reference_tac_integral;
    int thread; 
    int nthreads;
};
/* ----------------------------- MNI Header -----------------------------------
@NAME       : calcTACIntegrals
@INPUT      : tacTable - the 2 dimensional TAC lookup table
              vol4d - the 4d volume data structure
s              ref_mask_array - the image mask used for this TAC
              brainMask - the brain mask used for this TAC
@OUTPUT     : integral_table - the completed integral table
@RETURNS    : (nothing)
@DESCRIPTION: Calculates the Time Activity Curve (TAC) integrals
@METHOD     :
@GLOBALS    :
@CALLS      :
@CREATED    : March 2008, Paul Lyon
@MODIFIED   : May 2011, Thomas Funck
---------------------------------------------------------------------------- */
void calcTACIntegrals( double *tacTable, int* brain_mask_array, int* ref_mask_array, double* integral_table, data* pet, int ref,int zstart, float* time_widths){
  int t, y, x,z, tmpi, tmpj, tmpk;
  int index3d, index4d,  lastIndex, maskIndex;
  int nframes=pet->tmax;
  int nslices = pet->zmax;
  int nrows   = pet->ymax;
  int ncols   = pet->xmax;
  int zend;
    if(ref==2) zend=zstart+1;
    else zend=nslices;
    /* Calculate the integral for the zero index in order to       */
    /* initialise integrals for the average ref. region            */
    /*                                                             */
    /* The integral for time t is equivalent to the area under the */
    /* TAC, calculated as area of the rectangle having height of   */
    /* the intensity mid-way between the difference from t to t+1, */
    /* and of width timeWidth[t]                                   */
    /* Remember to add the previous as integrals are cumulative    */

    /* Iterate through the voxel data, slice, row, then column at a time... */

  
	if( ref == 1  )	/*	Calculate Integrals for reference region:				   */
	{
		for (t=0; t<nframes; t++)		/* ref reg voxels are averaged together, only have t dim */
        {
			if (t>0) integral_table[t] = integral_table[t-1] + time_widths[t] *  ((tacTable[t-1]+tacTable[t])/2.0);	
			/*Integral=last value + frame duration * average value 	*/
      		else integral_table[t] = time_widths[t] * (tacTable[t]/2.0);															
			/*For first frame, just do frame duration * voxel/2		*/
	 	}
	
    } 
  else
  {

   for (z=zstart; z<zend; z++){   
    for (y=0; y<nrows; y++){
      for (x=0; x<ncols; x++){
        for (t=0; t<nframes; t++){	
			/*Set index for mask according to t, y, x dimensions 		*/
             index4d = ((nslices*nrows*ncols)*t)+z*ncols*nrows + (ncols*y)+x; 																		
            /*Set index for slab according to t, z, y, x dimensions 		*/
             index3d  = ((nrows*ncols)*z)+(ncols*y)+x;                        						     
            /*Set index for mask according to z, y, x dimensions 		*/

	        if ( brain_mask_array[index3d] == 1 && ref_mask_array[index3d] == 0 ) {
                if (t>0) integral_table[index4d] = integral_table[lastIndex] + time_widths[t] *  ((tacTable[lastIndex]+tacTable[index4d])/2.0);			   
		/*Integral=last value + frame duration * average value 	*/
                else integral_table[index4d] = time_widths[t] * (tacTable[index4d]/2.0);								
		/*For first frame, just do frame duration * voxel/2		*/
	      }
          else  integral_table[index4d]=0.0; 
/* Ignore if outside the brain mask or in the ref. reg   */
	lastIndex=index4d;																										/* current tableIndex gets stored for next integral calc	*/
         }
       }
     }
   }          
  }
  

}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : patlakPlot
@INPUT      : tacTable - the 2 dimensional TAC lookup table
              integral_table - the 2 dimensional integral table
              vol4d - the 4d volume data structure
              ref_mask_array - the image mask used for this TAC
              brainMask - the brain mask used for this TAC
@OUTPUT     : patlakTable - the completed 3 dimensional Logan table
@RETURNS    : (nothing)
@DESCRIPTION: Generates a Logan plot for this image 
              Performs linear regression to fit the plot
              (Logan used as PK11195 is reversible)
@METHOD     :
@GLOBALS    :
@CALLS      :
@CREATED    : March 2008, Paul Lyon
@MODIFIED   :
---------------------------------------------------------------------------- */

int patlak_plot_multithread(data* pet, data* output, int* ref_mask_array, int*  brain_mask_array, double* integral_table,  double* reference_tac,double *reference_tac_integral,float* time_widths, parameters*  user_options, int nthreads){
    int rc;
    pthread_t threads[nthreads];
    struct mtga_args thread_args[nthreads];
    for(int t=0; t< nthreads; t++){
        thread_args[t].pet=pet;
        thread_args[t].output=output;
        thread_args[t].ref_mask_array=ref_mask_array;
        thread_args[t].brain_mask_array=brain_mask_array;
        thread_args[t].time_widths=time_widths;
        thread_args[t].user_options=user_options;
        thread_args[t].thread=t;
        thread_args[t].nthreads=nthreads;
        thread_args[t].integral_table=integral_table;
        thread_args[t].reference_tac=reference_tac ;
        thread_args[t].reference_tac_integral=reference_tac_integral ;
        rc= pthread_create(&threads[t], NULL, patlak_plot, (void *) &thread_args[t] ); 
        if(rc!=0) pexit("Error: creating thread","" , 1);
    }
    for(int t=0; t<nthreads;t++){
        rc = pthread_join(threads[t], NULL);
        if(rc!=0) pexit("Error: joining thread", "", 1);
    }
    return(0);
}

int patlak_plot(void *args){
    data* pet = ((struct mtga_args*) args)->pet ; 
    data* output = ((struct mtga_args*) args)->output ; 
    int* ref_mask_array= ((struct mtga_args*) args)->ref_mask_array ;
    int* brain_mask_array= ((struct mtga_args*) args)->brain_mask_array;
    float* time_widths= ((struct mtga_args*) args)->time_widths;
    parameters* user_options= ((struct mtga_args*) args)->user_options;
    double* reference_tac = ((struct mtga_args*) args)->reference_tac ;
    double* reference_tac_integral = ((struct mtga_args*) args)->reference_tac_integral ;
    double* integral_table=((struct mtga_args*) args)->integral_table ;
    int thread = ((struct mtga_args*) args)->thread; 
    int nthreads = ((struct mtga_args*) args)->nthreads ;
    int thread_step=nthreads-1;
    int thrown      = 0;
    int npoints     = 0;
    double slope    = 0.0, ssd, ic;
    int t, x, y, z;
    int *tableIndex, index3d, index4d;
    int spatialIndex;
    int nframes = pet->tmax;
    int nslices = pet->zmax;
    int nrows   = pet->ymax;
    int ncols   = pet->xmax;
    double refIntegrals[nframes];      
    float total;
    int trigger=0;
    double min=0, max=0;
    double voxIntegrals[nframes];
    double Y[nframes];
    double X[nframes];
    double Y_all[nframes];
    double X_all[nframes];
            
    /* Calculate the number of points to include based on the cutoff */
    tableIndex=&index4d ;

    /* Calculate the number of points to include based on the cutoff */
      

    for (t=0; t<nframes; t++)
      { 
        if ( t-user_options->start_frame >= 0.0  )
        { 
          npoints=nframes-t;
          break; 
        }
      }
    thrown=t;
 /* Iterate through the voxel data, slice, row, then column at a time... */
 for(z=thread; z<nslices; z+=thread_step){   
     calcTACIntegrals(pet->data, brain_mask_array, ref_mask_array, integral_table, pet, 2,z, time_widths );   
    for (y=0; y<nrows; y++) {
      for (x=0; x<ncols; x++) {    
		  spatialIndex = ((nrows*ncols)*z)+(ncols*y)+x;
        /* If we are inside the brain mask but outside the ref. region... */ 
	    if (brain_mask_array[spatialIndex] == 1 && ref_mask_array[spatialIndex] == 0){
          /* Set up the x-axis points using the ref. region data => Intgrl(REF)/REF */
          /* Iterate through each frame, excluding the last...   */
            total=0;
                for (t=0; t<nframes; t++){
                    index4d = ((nrows*ncols*nslices)*t) +  ((nrows*ncols)*z) +(ncols*y)+x;
                    total +=  pet->data[*tableIndex];
                    //printf("%f\n", tacTable[index4d]);
                        if (pet->data[*tableIndex] == 0.0 || ref_mask_array[t] == 0.0)  { 
                           X_all[t] = 0.0; 
                           Y_all[t] = 0.0;
                        } 
            		    else {
                             Y_all[t]=pet->data[*tableIndex]/reference_tac[t];   
            			     X_all[t]=reference_tac_integral[t]/reference_tac[t];
                        }
            		    if (t >= thrown) { 
                            X[t-thrown]=X_all[t];  
                            Y[t-thrown]=Y_all[t];
                        } 
			        /* Store ratio if we are past the points thrown away*/
		            }
	/* Store ratio if we are past the points thrown away*/
	regression_choice(/* inputs */ X, Y, npoints,/* outpus */ &slope, &ic, &ssd, user_options); /* Calc slope of the linear regression and store it */

    output->data[spatialIndex]=(float) slope;
    //printf("%f\n", slope);
    if( max < slope ) max=slope;
	if( min > slope ) min=slope;
        float time=0; 
	    //if( isinf(slope) ||  isnan(slope) )   
        if(slope > 0 && trigger ==0){
			printf("time\tREF\tREF INT\tROI\tROI INT\tX-AXIS\tY-AXIS\n");

            trigger++;
            for(t=0; t<nframes; t++ ){
                index4d = ((nrows*ncols*nslices)*t) +  ((nrows*ncols)*z) +(ncols*y)+x;
                time += time_widths[t];
                printf("%f\t%f\t%f\t%f\t%f\t%f\t%f\n", time_widths[t], ref_mask_array[t], refIntegrals[t], pet->data[index4d], integral_table[index4d],  X_all[t],   Y_all[t] );
            }
            printf("Slope=%f\n", slope);


		} 
    } 
	else  	output->data[spatialIndex]=0.0;   /* Ignore everthing else */  
   }
  }
 }  
    if(strcmp(user_options->output_unit, "CMRg") == 0){  
        for (z=thread; z < nslices; z+=thread_step) for (y=0; y < nrows; y++) for (x=0; x < ncols; x++){ 
            output->data[z*nrows*ncols+y*ncols+x] *= (user_options->PG/user_options->LC);
        }
    }

  /* Clean up memory allocation */
    return(0);
}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : //loganPlot
@INPUT      : tacTable - the 2 dimensional TAC lookup table
              integral_table - the 2 dimensional integral table
              vol4d - the 4d volume data structure
              ref_mask_array - the image mask used for this TAC
              brainMask - the brain mask used for this TAC
@OUTPUT     : loganTable - the completed 2 dimensional Logan table
@RETURNS    : (nothing)
@DESCRIPTION: Generates a Logan plot for this image 
              Performs linear regression to fit the plot
              (Logan used as PK11195 is reversible)
@METHOD     :
@GLOBALS    :
@CALLS      :
@CREATED    : March 2008, Paul Lyon
@MODIFIED   :
---------------------------------------------------------------------------- */

int logan_plot_multithread(data* pet,data* output, int* ref_mask_array, int*  brain_mask_array,double* integral_table,  double* reference_tac,double *reference_tac_integral,float* time_widths,parameters*  user_options, int nthreads){
    int rc;
    pthread_t threads[nthreads];
    struct mtga_args thread_args[nthreads];
    for(int t=0; t< nthreads; t++){
        thread_args[t].pet=pet;
        thread_args[t].output=output;
        thread_args[t].ref_mask_array=ref_mask_array;
        thread_args[t].brain_mask_array=brain_mask_array;
        thread_args[t].time_widths=time_widths;
        thread_args[t].user_options=user_options;
        thread_args[t].thread=t;
        thread_args[t].nthreads=nthreads;
        thread_args[t].integral_table=integral_table;
        thread_args[t].reference_tac=reference_tac ;
        thread_args[t].reference_tac_integral=reference_tac_integral ;
        rc= pthread_create(&threads[t], NULL, logan_plot, (void *) &thread_args[t] ); 
        if(rc!=0) pexit("Error: creating thread","" , 1);
    }
    for(int t=0; t<nthreads;t++){
        rc = pthread_join(threads[t], NULL);
        if(rc!=0) pexit("Error: joining thread", "", 1);
    }
    return(0);
}


int logan_plot(void* args) {
    data* pet = ((struct mtga_args*) args)->pet ; 
    data* output = ((struct mtga_args*) args)->output ; 
    int* ref_mask_array= ((struct mtga_args*) args)->ref_mask_array ;
    int* brain_mask_array= ((struct mtga_args*) args)->brain_mask_array;
    float* time_widths= ((struct mtga_args*) args)->time_widths;
    parameters* user_options= ((struct mtga_args*) args)->user_options;
    double* integral_table =((struct mtga_args*) args)->integral_table ;
    double* reference_tac = ((struct mtga_args*) args)->reference_tac ;
    double* reference_tac_integral = ((struct mtga_args*) args)->reference_tac_integral ;
    int thread = ((struct mtga_args*) args)->thread; 
    int nthreads = ((struct mtga_args*) args)->nthreads ;
    int thread_step=nthreads;//-1;
    int thrown      = 0;
    int npoints     = 0;
    double slope    = 0.0, ic, ssd;
    int t, x, y, z;
    int index3d, index4d;
    int spatialIndex;
    int nframes = pet->tmax;
    int nslices = pet->zmax;
    int nrows   = pet->ymax;
    int ncols   = pet->xmax; 
    double refIntegrals[nframes];
    int trigger=0;
    double voxIntegrals[nframes]; // = malloc((nframes) * sizeof(double));
    double Y[nframes];//    = malloc((npoints) * sizeof(double));
    double X[nframes];//    = malloc((npoints) * sizeof(double));
    double Y_all[nframes];// = malloc((nframes) * sizeof(double));
    double X_all[nframes];// = malloc((nframes) * sizeof(double));
    double min=0, max=0;

    /* Calculate the number of points to include based on the cutoff */
 
    for (t=0; t<nframes; t++){ 
        if ( t - user_options->start_frame  >= 0.0 ){
          npoints=nframes-t;
          break; 
        }
    }
    thrown=t;
    /* If we are keeping only 3 or less points then continuing is futile */
      if (npoints <= 3)
      {
        fprintf(stderr,
                "Specified Logan time cutoff (%d) only leaves %d data points!\n",
                (int) user_options->start_frame, npoints);
        exit(EXIT_FAILURE);
      }

                    
      /* Iterate through the voxel data, slice, row, then column at a time... */
        for (z=thread; z < nslices; z+=thread_step){ 
            calcTACIntegrals(pet->data, brain_mask_array, ref_mask_array, integral_table, pet, 2,z, time_widths );
            for (y=0; y < nrows; y++){
                for (x=0; x < ncols; x++){
                    /* Ingore index 0 of loganTable so that we can map to tables directly */
                    spatialIndex = ((nrows*ncols)*z)+(ncols*y)+x;
                    /* If we are inside the brain mask but outside the ref. region... */ 
                    if (brain_mask_array[spatialIndex] == 1 && ref_mask_array[spatialIndex] == 0){
                    /* Set up the x-axis points using the ref. region data */
                    /* Iterate through each frame, excluding the last...   */
                        for (t=0; t<nframes; t++){
                        index4d = ((nrows*ncols*nslices)*t) +  ((nrows*ncols)*z) + (ncols*y) + x;
                        //printf("%d %f %f %f %f\n", t, reference_tac_integral[t],  reference_tac[t],pet->data[index4d], integral_table[index4d] ); 
                            if (pet->data[index4d] == 0.0) { 
                               X_all[t] = 0.0; 
                               Y_all[t] = 0.0;
                            } 
                            else{    
                                Y_all[t]=integral_table[index4d]/pet->data[index4d];  
                                if (user_options->arterial_fn != NULL) X_all[t]=(reference_tac_integral[t]/pet->data[index4d]); 
                                else  X_all[t]=(reference_tac_integral[t]+(reference_tac[t] / user_options->k2))/pet->data[index4d]; 	         
                            }
                            if (t >= thrown){ 
                                X[t-thrown]=X_all[t];  
                                Y[t-thrown]=Y_all[t];
                            } 
                        }				
                     
                        regression_choice(/* inputs */ X, Y, npoints,/* outputs */ &slope, &ic, &ssd, user_options);                     
                        /* Calc slope of the linear regression and store it */
                        output->data[spatialIndex]=  slope;

                   } 
                   else  output->data[spatialIndex]=0.0;   /* Ignore everthing else */    
                }
            }
        }

        if(strcmp(user_options->output_unit, "bp")==0){  
            for (z=thread; z < nslices; z+=thread_step){ 
                for (y=0; y < nrows; y++){
                    for (x=0; x < ncols; x++){
                        spatialIndex = ((nrows*ncols)*z)+(ncols*y)+x;
                        if (brain_mask_array[spatialIndex] == 1 && ref_mask_array[spatialIndex] == 0){
                            output->data[z*nrows*ncols+y*ncols+x] -= 1;
                        }
                    }
                }
            }
        }           

    return(0);
}




static int regression_choice(double *x, double *y, int nr, double *slope, double *ic, double *ssd, parameters* user_options){

    //if (user_options->regression_type=="P") {llsqperp3(x,y, nr, slope, ic, ssd);}
    /* add more regression options here...*/
    *slope=calcLinearRegression(nr,x,y); 

return 0;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : calcLinearRegression
@INPUT      : n - the total number of points
@OUTPUT     : (nothing)
@RETURNS    : The gradient of the linear regression
@DESCRIPTION: Calculates the slope of a simple linear regression performed on
              the set of coordinates {x,y} using least squares fitting.
@METHOD     :
@GLOBALS    :
@CALLS      :
@CREATED    : March 2008, Paul Lyon
@MODIFIED   :
---------------------------------------------------------------------------- */
static double calcLinearRegression(int n, double* x, double* y)
{
  int i;
  double meanX    = 0.0;
  double meanY    = 0.0;
  double top      = 0.0;
  double bottom   = 0.0;


  /* Calculate meanX and meanY */
  for (i=0; i<n; i++)
  {
    meanX+=x[i];
    meanY+=y[i];
  }
  meanX/=(double)n;
  meanY/=(double)n;

  /* Now iterate over i to calculate the simple linear regression using the  *
   * method of ordinary least squares which is illustrated as follows:       *
   *                                                                         *
   *   slope = sumi( (xi-meanX).(yi-meanY) ) / sumi(xi-meanX)^2              *
   *                                                                         *
   * where sumi is the sum over index i of variable x or y                   */

  for (i=0; i<n; i++)
  {
    top+=(x[i]-meanX)*(y[i]-meanY);
    bottom+=(x[i]-meanX)*(x[i]-meanX);
  }

	

  return (double)top/bottom;
}


float* get_time_widths(int nframes,  char* inputfile){
    char command[200];
    FILE *input_f=fopen(inputfile, "rt");
    char buffer[100], dlm[]=" ", *token;
    int t;
    int unit;
    float* time_widths=malloc( nframes *sizeof(*time_widths));

    //if(userOptions.min==TRUE) unit=60.0;
    //else unit=1.0;

    for(t=0; t<nframes; t++){
            fgets(buffer, sizeof(buffer), input_f);
            time_widths[t]=atof(buffer);
    }

    fclose(input_f);
    return(time_widths);
}


double** importArterialTAC (char* arterial_filename, int* points, parameters* userOptions){
    char buffer[100];
    char *dlm="	";
    int counter=0;
    double value, time;
    FILE* arterial_file=fopen(arterial_filename, "rt");
    double** arterialTAC=NULL; // 2D array for time [0] and radioactivity [1]
    //printf("Time\tValue\n");
    while(fgets(buffer, sizeof(buffer), arterial_file) != NULL){
        counter++;
        arterialTAC=realloc(arterialTAC, counter * sizeof(*arterialTAC));
        arterialTAC[counter-1]=malloc(sizeof(**arterialTAC)*2);

        arterialTAC[counter-1][0]=atof(strtok(buffer, dlm));
        //if( userOptions->arterial_convert_to_seconds ==TRUE)  arterialTAC[counter-1][0] *= 60;

        arterialTAC[counter-1][1]=atof(strtok(NULL, dlm));

        //printf("%s\t:\t%f\t%f\n", buffer, arterialTAC[counter-1][0], arterialTAC[counter-1][1]);
    }
    *points=counter;
    //printf("Imported arterial TAC file.\n");

    fclose(arterial_file);
    return(arterialTAC);
}

double* get_arterial_tac(int nframes,float* times, char* arterial_filename, parameters* userOptions  ){
    double* referenceTAC=malloc(nframes * sizeof(*referenceTAC) );
    int npoints=0;
    double** arterialTAC=importArterialTAC( arterial_filename, &npoints, userOptions); //read the arterial input from text file
        
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
                double d=arterialTAC[max][0] - arterialTAC[min][0];
                double d0 = 1-(times[i] - arterialTAC[min][0]) / d;
                double d1= 1-(arterialTAC[max][0] - times[i]) / d;
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

double* calculate_ref_tac(data* pet, int* ref_mask_array, float* time_widths, parameters* user_options){
    double* reference_tac;
    float sum=0;
    
    //OPTION 1:
    //Get reference TAC from arterial input
    if( user_options->arterial_fn != NULL){
        float times[pet->tmax];
        for(int t=0; t<pet->tmax; t++){
            times[t]=sum + time_widths[t]*0.5;
            sum+=time_widths[t];
        } 
        reference_tac=get_arterial_tac(pet->tmax, times, user_options->arterial_fn, user_options  );
        
    } else{
        //OPTION 2:
        //Get reference TAC from reference region in PET image
        int n=0;
        reference_tac=calloc(pet->tmax, sizeof(*reference_tac));
        for(int t=0; t< pet->tmax; t++){
            n=0;
            for(int z=0; z< pet->zmax; z++){
                for(int y=0; y< pet->ymax; y++){
                    for(int x=0; x< pet->xmax; x++){ 
                        int index= z*pet->ymax*pet->xmax+y*pet->xmax+x;
                        if( ref_mask_array[index] != 0){
                            int index4d=t * pet->zmax * pet->ymax * pet->xmax + z * pet->ymax * pet->xmax+ y * pet->xmax+x; 
                            reference_tac[t] += pet->data[ index4d ];
                            n++;
                        }
                    }
                }
            }
            reference_tac[t] /= n;
        }
    }

    return(reference_tac);
}

int main(int argc, char** argv){
    if(argc == 1 || strcmp(argv[1], "-help") == 0 ) useage();
    int* ref_mask_array, *brain_mask_array;
    parameters user_options;
    data *output, pet, brain_mask, ref_mask;
    double* integral_table;
    double* reference_tac;
    double* reference_tac_integral;
    int nthreads=sysconf(_SC_NPROCESSORS_ONLN);

    //Set default options
    user_options.time_unit="sec";    
    user_options.ref_mask_fn=NULL;
    user_options.brain_mask_fn=NULL;
    user_options.arterial_fn=NULL;
    user_options.output_unit="";

    if( parse_input(argc, argv, &user_options) != 0)  useage();
    
    //Read PET file
    pet.filename=user_options.pet_fn;
    pet.data = (double*) readVolume(&pet, 1, MI_TYPE_DOUBLE );
    //Read Reference Mask
    if( user_options.ref_mask_fn != NULL ){
        ref_mask.filename=user_options.ref_mask_fn;
        ref_mask_array =(int*) readVolume(&ref_mask, 1, MI_TYPE_INT );
    }
    //Read Brain Mask
    if( user_options.brain_mask_fn != NULL ){
        brain_mask.filename=user_options.brain_mask_fn;
        brain_mask_array = (int*) readVolume(&brain_mask, 1, MI_TYPE_INT );
    }

    //Get time-widths
    float* time_widths=get_time_widths(pet.tmax, user_options.time_widths_fn);

    if( strcmp(user_options.analysis, "pp") ||  strcmp(user_options.analysis, "lp") ){
        //Create 3D outputfile
        output=copyVolume(&pet);
        output->tmax=1;
        output->filename=user_options.output_fn;
        output->n=output->zmax*output->ymax*output->xmax;
        output->start[0]=pet.zstart; 
        output->start[1]=pet.ystart;
        output->start[2]=pet.xstart;
        output->step[0]=pet.zstep;
        output->step[1]=pet.ystep;
        output->step[2]=pet.xstep;
        output->wcount[0]=pet.zmax;
        output->wcount[1]=pet.ymax;
        output->wcount[2]=pet.xmax;
        output->data=calloc(output->n, sizeof(*output->data));

        integral_table=calloc(pet.n, sizeof(*integral_table));
        //Get reference TAC
        reference_tac=calculate_ref_tac(&pet, ref_mask_array, time_widths, &user_options);
        //Allocate memory for reference integral
        reference_tac_integral = calloc(pet.tmax, sizeof(*reference_tac_integral));
        //Find tac integral for reference region
        calcTACIntegrals(reference_tac, brain_mask_array, ref_mask_array,  reference_tac_integral,  &pet, 1, 0, time_widths);
    }
    else if( strcmp(user_options.analysis, "suv" ) ){
        //Create 4D output file
        output=copyVolume(&pet);
        output->data=malloc(output->n*sizeof(*output->data)); 
    }   
    else return 1;

    //////////////////
    // Run analysis //
    //////////////////
    if(strcmp(user_options.analysis, "lp") ==0 ){ //Logan Plot
        logan_plot_multithread(&pet, output, ref_mask_array, brain_mask_array, integral_table, reference_tac, reference_tac_integral, time_widths, &user_options, nthreads);
    }
    else if (strcmp(user_options.analysis, "pp")==0){ //Patlak Plot
        patlak_plot_multithread(&pet, output, ref_mask_array, brain_mask_array, integral_table, reference_tac, reference_tac_integral, time_widths,  &user_options, nthreads);
    }
    //else if (strcmp(user_options.analysis, "suv") == TRUE){ //SUV
    //  suv_multithread(&pet, output,  ref_mask_array, brain_mask_array, &user_options); 
    //}
    else return 1;
    writeVolume( output->filename, output->data, output->start, output->step, output->wcount, MI_TYPE_DOUBLE  );

    return(0);
}
 

int parse_input(int argc, char** argv,  parameters* user_options){
    int i, temp_i=0;
    int er=0;
    char* linear="-linear";
    char* nearest="-nearest";
    char* analysis_list[]={"lp", "pp", "suv"};
    char* output_units[]={"bp", "Vt", "K", "CMRg" };
    char* time_units[]={"sec", "min"};
    int set_pet=FALSE;
    int set_analysis=FALSE;
    int set_ref_mask=FALSE;
    int set_brain_mask=TRUE;
    int set_k2 =FALSE;
    int set_arterial=FALSE;
    int set_body_weight=FALSE;
    int set_injected_dose=FALSE;
    int set_output_fn=FALSE;
    int set_LC=FALSE;   
    int set_PG=FALSE;
    int set_time_widths=FALSE;
    int set_input_dose=FALSE;
    int set_output_unit=FALSE;
    int set_start_frame;

    
    for( i=1; i < argc; i++ ){
        //printf("%d\n", i);
        if(strcmp(argv[i], "-pet" ) == 0 ){
            i++;
            user_options->pet_fn=argv[i];
            set_pet=TRUE;
        }
        else if(strcmp(argv[i], "-ref_mask" ) ==0 ){
            i++;
            user_options->ref_mask_fn=argv[i];
            set_ref_mask=TRUE; 
        }
        else if(strcmp(argv[i], "-brain_mask" )==0 ){
            i++;
            user_options->brain_mask_fn=argv[i];
            set_brain_mask=TRUE; 
        }
        else if(strcmp(argv[i], "-o" )==0 ){
            i++;
            user_options->output_fn=argv[i];
            set_output_fn=TRUE; 
        }
        else if(strcmp(argv[i], "-arterial" )==0 ){
            i++; //
            user_options->arterial_fn=argv[i]; 
            set_arterial=TRUE;
        }
        else if(strcmp(argv[i], "-start_frame" )==0 ){
            i++; //
            user_options->start_frame=atoi(argv[i]); 
            set_start_frame=TRUE;
        }
        else if(strcmp(argv[i], "-time_widths" )==0 ){
            i++; //
            user_options->time_widths_fn=argv[i]; 
            set_time_widths=TRUE;
        }
        else if(strcmp(argv[i], "-k2" )==0 ){
            i++; //
            user_options->k2=atof(argv[i]); 
            set_k2=TRUE;
        }
        else if(strcmp(argv[i], "-PG" )==0 ){
            i++; //
            user_options->PG=atof(argv[i]); 
            set_PG=TRUE;
        }
        else if(strcmp(argv[i], "-LC" )==0 ){
            i++; //
            user_options->LC=atof(argv[i]); 
            set_LC=TRUE;
        }
        else if(strcmp(argv[i], "-body_weight" )==0 ){
            i++; //
            user_options->body_weight=atof(argv[i]); 
            set_body_weight=TRUE;
        }
        else if(strcmp(argv[i], "-injected_dose" )==0 ){
            i++; //
            user_options->injected_dose=atof(argv[i]); 
            set_input_dose=TRUE;
        }
        else if(strcmp(argv[i], "-time_unit" )==0 ){
            i++; //
            while(time_units[i] != NULL){
                if( strcmp(argv[i], time_units[i])  ){ 
                    user_options->time_unit=argv[i];
                    break;
                }
            }
        }
        else if(strcmp(argv[i], "-output_unit" )==0 ){
            i++; //
            for(int j=0; analysis_list[j] != NULL; j++){
                if( strcmp(argv[i], analysis_list[j])  ){ 
                    user_options->output_unit=argv[i];
                    set_output_unit=TRUE;
                    break;
                }
            }
        }
        else if(strcmp(argv[i], "-analysis" )==0 ) {
            i++;
            for(int j=0; analysis_list[j] != NULL; j++){
                if( strcmp(argv[i], analysis_list[j])  ){
                    user_options->analysis=argv[i];
                    set_analysis=TRUE;
                    break;
                }
            }
        }
    }


    if(set_analysis == FALSE){ 
        er=1;
        printf("Could not find analysis type, specify it with: -analysis <lp/pp>\n");
        return(er);
    }
    if(set_output_fn == FALSE){ 
        er=1;
        printf("Could not find output file type, specify it with: -o <output minc>\n");
        return(er);
    }
    if(set_start_frame==FALSE){
        printf("Could not find start frame, speficy with: -start_frame <integer>");
        er=1;
    }
    else if(strcmp(user_options->analysis, "lp") == TRUE){ //Logan Plot
        if(set_ref_mask == TRUE && set_k2 != TRUE) er=1;
        if(set_arterial==FALSE && set_ref_mask == FALSE) er=1; 
        if(set_start_frame == FALSE) er=1;
        if(set_time_widths == FALSE) er=1;
    }
    else if(strcmp(user_options->analysis, "pp") == TRUE){ //Patlak Plot
        if(set_arterial==FALSE && set_ref_mask == FALSE) er=1; 
        if(set_start_frame == FALSE) er=1;
        if( strcmp(user_options->output_unit, "CMRg")==TRUE){
            if(set_PG != TRUE || set_LC != TRUE) er=1;
        }
        if(set_time_widths == FALSE) er=1;
    }
    else if(strcmp(user_options->analysis, "sv") == TRUE){ //SUV 
        if(set_body_weight == FALSE) er=1;
        if(set_input_dose == FALSE) er=1;
    }
return(er);
}
 



void useage(){
    printf("\nInvalid parameters were specified!\nPlease put desired options before the analysis type.\n");
    printf("Usage:\n");
    printf("[-Options] -lp <k2> <input_image.mnc> [<ref_mask.mnc> or plasma_input_TAC.txt]\n\t <brain_mask.mnc> <start frame> <time widths file> <outputfile>\n");
    printf("[-Options] -pp <input_image.mnc> <ref_mask.mnc> <brain_mask.mnc> <start frame>  <time widths file> <outputfile>\n");
    printf("[-Options] -sr <input_image.mnc> <ref_mask.mnc> <brain_mask.mnc> <start frame> <time widths file> <outputfile>\n");
    printf("[-Options] -roi <input_image.mnc> <ROI_mask.mnc> <brainmask.mnc>\n");
    printf("[-Options] -sv  <injected dose> <patient's body weight> <start_frame> <end_frame> <input_image.mnc> <brainmask.mnc>\n");
    printf("General Options:\n");
    printf(" -arterial		Arterial input TAC, use .txt file instead of .mnc ref mask\n");
    printf(" -autok2		Automatically get tracer information from database\n");
    printf(" -bp 		Calculate binding potential\n");
    printf(" -RAM		Set maximum amount of RAM that will be used.\n");
    printf(" -min		Perform analysis in minutes.\n");
    printf(" -testStartTime     Test start time for Patlak/Logan Plot.\n");
    //printf(" -P  		Perpendicular regression model \n");   

    printf("Where:\n");
    printf(" -lp performs a Logan plot\n");
    printf(" -pp performs a Patlak plot\n");
    printf(" -sr performs the simplified reference tissue model plots\n");
    printf(" -sv calcultes standardized uptake value for each voxel.\n");   
    printf("Parameters are as follows:\n");
    printf("  <input_image.mnc>   - the name of the input file\n");
    printf("  <ref_mask.mnc>      - the reference region mask file\n");
    printf("  <brain_mask.mnc>    - the brain mask file\n");
    exit(1);
}
