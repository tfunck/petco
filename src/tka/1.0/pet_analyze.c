#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <malloc.h>
#include "llsqwt.h"
#include "hdf5.h"
#include "minc2.h"
#include "mtga.h"
#include "srtm.h"
#include "fft.h"

static void  patlakPlot(float *tacTable, char* petimage_file, float *outputTable, imageData imageInfo, misize_t *rstart, misize_t* rcount,   int startTime, parameters  userOptions,float*  max, float *min);

static void calcTACIntegrals(float *tacTable, float *integralTable, imageData *imageInfo, int ref, int RAMexceeded, int zstart)  ;  
void loganPlot(float* tacTable, char* petimage_file, float*  outputTable, imageData imageInfo,misize_t rstart[4], misize_t rcount[4], int startTime, parameters  userOptions, float *max ,float*  min, float *k2) ;
static float  calcLinearRegression(int n, double* x, double* y) ; 

static float testStartTime(float* tacTable, float*  outputTable, imageData imageInfo, misize_t rstart[4], misize_t rcount[4], parameters  userOptions, float *k2, float *refIntegrals, float* integralTable );

static int regression_choice(double *x, double *y, int nr, double *slope, double *ic, double *ssd, parameters user_options) ;

/* ----------------------------- MNI Header -----------------------------------
@NAME       : processSlices
@INPUT      : 
@OUTPUT     : 
@RETURNS    : (nothing)
@DESCRIPTION: 
@METHOD     :
@GLOBALS    :
@CALLS      :
@CREATED    : 
@MODIFIED   : 
---------------------------------------------------------------------------- */


int processSlices(imageData imageInfo, char* petimage_file,  parameters userOptions, int startFrame)
{
int x, y, z=0, t=0;
float *hyperslab; 
float *outputTable;       
/*The final results from the model are stored in the outputTable*/
float max=0.0, min=0.0;   
/*Maximum and minimum values for slice            */
misize_t nframes=imageInfo.sizes[0];      
/*Number of time frames, t                            */
misize_t nslices=imageInfo.sizes[1];      
/*Number of slices, z                             */
misize_t nrows=imageInfo.sizes[2];			
/*Number of rows, y                               */
misize_t ncols=imageInfo.sizes[3];			
/*Number of columns, x                            */
misize_t  rstart[4];  
/*Start values for reading data in a given dimension, 4D	          */
misize_t  rcount[4];  
/*Total length of data that will be read for x, y, z, t dimensions, 4D    */
misize_t  *wstart;  
/*Start values for writing data in a given dimension, 3D             */
misize_t  *wcount;  
/*Total length of data that will be written for x, y, z, t dimensions, 4D    */
double  RAMtotal; /* Total, approximated, amount of RAM used*/
double totalLogan=0.0, nLogan=0.0;
double *averageTable;
int *n;
int i, tableindex, maskindex;
tracer* k2parameters;
 
  if(userOptions.analysis_type != 3) 
  /*The analysis type is not SUV -> the output is a 3d image*/ 
  { 
    outputTable=malloc( nslices * nrows * ncols * sizeof(float));
    printf("Sizes: %d %d %d, Print %d\n",  nrows, ncols, nslices, outputTable);
    wcount=malloc(3*sizeof(misize_t));
    wstart=malloc(3*sizeof(misize_t));
    
   	wstart[0]=wstart[1]=wstart[2]=(misize_t) 0; 
	/* set start values for z, y and x to 0			*/
	wcount[0]=(misize_t) nslices; 				
	/* output write 1 voxel in z dim.           	*/
	wcount[1]=(misize_t) nrows;				
	/* output write "nrows" voxels in y dim.     	*/
	wcount[2]=(misize_t) ncols;   				
    /* output write "ncols" voxels in x dim.     	*/
  }    
  else  /*i.e., the image is a 4d output*/
  {
      outputTable=malloc( nframes * nrows * ncols * nslices * sizeof(float));
        printf("Sizes: %d %d %d %d\n", nframes, nrows, ncols, nslices); exit(0);
      wcount=malloc(4*sizeof(misize_t));
      wstart=malloc(4*sizeof(misize_t));    

       wstart[0]=wstart[1]=wstart[2]=wstart[3]=0; 
	/* set start values for z, y and x to 0			*/
	   wcount[0]=(misize_t) nframes; 				
	/* output write 1 voxel in z dim.           	*/
	   wcount[1]=(misize_t) nslices;				
	/* output write "nrows" voxels in y dim.     	*/
	   wcount[2]=(misize_t) nrows;   				
    /* output write "ncols" voxels in x dim.     	*/   
       wcount[3]=(misize_t) ncols;
  }

    k2parameters=imageInfo.tracerInfo;
	rstart[0]=rstart[1]=rstart[2]=rstart[3]= 0; 
	/*set start values for t, y and x to 0 			*/
	rcount[0]=(misize_t) nframes;				
	/* output write "nframes" voxels in t dim.   	*/
	rcount[1]=(misize_t) nslices;	   				
	/* output write 1 voxels in z dim.				*/
	rcount[2]=(misize_t) nrows;   				
	/* output read "nrows" voxels in y dim.      	*/
	rcount[3]=(misize_t) ncols;  				
	/* output read "ncols" voxels in x dim.      	*/	

	/*rough calculation of RAM usage based on size of arrays that have to be allocated*/
    RAMtotal=(double) 9.313 * pow(10, -10) * (sizeof(float)*nslices*nframes*nrows*ncols + 3*sizeof(float)*nrows*ncols*nslices);
       	/*additional memory is needed for logan or patlak plot, I dont remember for what though (maybe reference regio averages?*/
    if(userOptions.analysis_type == 0 || userOptions.analysis_type==1 ) RAMtotal +=(double) 9.313*pow(10, -10) * sizeof(float) * nframes*nrows*ncols*nslices ;    
    printf("total total %f\n", RAMtotal); 
    
    if(RAMtotal < userOptions.RAMthreshold ) {
        printf("\nRAM not exceeded. A little over %f GB of RAM will be used for this program.\n", RAMtotal);
        userOptions.RAMexceeded=FALSE;
        hyperslab=  malloc(nrows*nslices*nframes*ncols*sizeof(float));
                //miget_real_value_hyperslab(image, MI_TYPE_FLOAT, rstart, rcount, hyperslab) ;
        loadpetimage( z, imageInfo.sizes, rstart, rcount, petimage_file, hyperslab, imageInfo.brainmask, imageInfo.referencemask, userOptions);     
        //for(int i=0; i < nslices*ncols*nrows*nframes; i++) if(hyperslab[i] > 4) printf("%f\n", hyperslab[i]);               
    }
    else {    
        printf("RAM threshold exceeded. A little over %f GB of RAM will be used for this program.\n", RAMtotal);   
        userOptions.RAMexceeded=TRUE;   
        	if (hyperslab=(float *) malloc(nframes*nrows*ncols*sizeof(float)) == NULL) {
			fprintf(stderr,"Insufficient memory for hyperslab.\n");  
        	exit(0);
	    }
	    rcount[1]=1;  
    }   

    //Scale values by conversion factor
    if( userOptions.conversion_factor != 1) printf("Applying conversion factor: %f\n", userOptions.conversion_factor);
    for(int i=0; i< nframes * ncols * nslices *nrows; i++) hyperslab[i] *= userOptions.conversion_factor;

    printf("Applying model to PET image.\n");      		
    		
       switch (userOptions.analysis_type)
       {
           case 0: {
                    float tempk2=k2parameters->k2mean;           
		            loganPlot(hyperslab, petimage_file, outputTable, imageInfo, rstart, rcount, startFrame, userOptions, &max , &min, &tempk2);
                   }
                    break;
           case 1: patlakPlot(hyperslab, petimage_file, outputTable, imageInfo, rstart, rcount, startFrame,  userOptions,  &max, &min );
                    break;
           default: break;        
       }

      if(userOptions.write == TRUE){
        printf("max=%f, min=%f\n", max, min);
        //miopen_volume("3D_image.mnc", MI2_OPEN_READ, &outputimage);
        // imageInfo.outputimage = &outputimage;
        mihandle_t outputimage;
        printf("%s\n",imageInfo.outputimage_file );
        if(miopen_volume(imageInfo.outputimage_file  , MI2_OPEN_RDWR, &outputimage) != MI_NOERROR )	fprintf(stderr, "\tError opening\n");
        if( miset_volume_range ( outputimage, max, min) != MI_NOERROR ) fprintf(stderr, "\tError setting range.\n");
	    if ( miset_real_value_hyperslab( outputimage, MI_TYPE_FLOAT, wstart, wcount, outputTable) != MI_NOERROR )	fprintf(stderr, "\tError setting hyperslab");
        miclose_volume( outputimage );
    }


return 0;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : calcTACIntegrals
@INPUT      : tacTable - the 2 dimensional TAC lookup table
              vol4d - the 4d volume data structure
s              imageInfo.referencemask - the image mask used for this TAC
              brainMask - the brain mask used for this TAC
@OUTPUT     : integralTable - the completed integral table
@RETURNS    : (nothing)
@DESCRIPTION: Calculates the Time Activity Curve (TAC) integrals
@METHOD     :
@GLOBALS    :
@CALLS      :
@CREATED    : March 2008, Paul Lyon
@MODIFIED   : May 2011, Thomas Funck
---------------------------------------------------------------------------- */
void calcTACIntegrals(float *tacTable, float *integralTable, imageData *imageInfo, int ref, int RAMexceeded, int zstart)
{
  int t, y, x,z, tmpi, tmpj, tmpk;
  int *tableIndex, index3d, index4d,  lastIndex, maskIndex, indexTYX;
  int nframes = imageInfo->sizes[0];
  int nslices = imageInfo->sizes[1];
  int nrows   = imageInfo->sizes[2];
  int ncols   = imageInfo->sizes[3];
  char filestr[100]="";

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
			if (t>0) integralTable[t] = integralTable[t-1] + imageInfo->timeWidths[t] *  ((tacTable[t-1]+tacTable[t])/2.0);	
			/*Integral=last value + frame duration * average value 	*/
      		else integralTable[t] = imageInfo->timeWidths[t] * (tacTable[t]/2.0);															
			/*For first frame, just do frame duration * voxel/2		*/
	 	}
	
    } 
  else
  {
  tableIndex=&index4d;   
  if(RAMexceeded==TRUE) { tableIndex=&indexTYX; nslices = zstart+1; }  

   for (z=zstart; z<nslices; z++)
   {    
    for (y=0; y<nrows; y++) 
    {
      for (x=0; x<ncols; x++)
      {
        for (t=0; t<nframes; t++) 
        {	
			indexTYX= ((nrows*ncols)*t)+(ncols*y)+x; 
			/*Set index for mask according to t, y, x dimensions 		*/
             index4d = ((nslices*nrows*ncols)*t)+z*ncols*nrows + (ncols*y)+x; 																		
		/*Set index for slab according to t, z, y, x dimensions 		*/
	     index3d  = ((nrows*ncols)*z)+(ncols*y)+x;                        						     
		/*Set index for mask according to z, y, x dimensions 		*/
      	imageInfo->brainmask[0]; 
	if ( imageInfo->brainmask[index3d] == 1.0 && imageInfo->referencemask[index3d] == 0.0 )
          {
            if (t>0) integralTable[*tableIndex] = integralTable[lastIndex] + imageInfo->timeWidths[t] *  ((tacTable[lastIndex]+tacTable[*tableIndex])/2.0);			   
		/*Integral=last value + frame duration * average value 	*/
            else integralTable[*tableIndex] = imageInfo->timeWidths[t] * (tacTable[*tableIndex]/2.0);								
		/*For first frame, just do frame duration * voxel/2		*/
       //printf("%d %d %d: %f %f\n", y, x, t, tacTable[tableIndex], integralTable[tableIndex]); 
	      }
          else  integralTable[*tableIndex]=0.0; 
/* Ignore if outside the brain mask or in the ref. reg   */
	lastIndex=*tableIndex;																										/* current tableIndex gets stored for next integral calc	*/
         }
       }
     }
   }          
  }
  

}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : patlakPlot
@INPUT      : tacTable - the 2 dimensional TAC lookup table
              integralTable - the 2 dimensional integral table
              vol4d - the 4d volume data structure
              imageInfo.referencemask - the image mask used for this TAC
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


static void  patlakPlot(float *tacTable, char* petimage_file,  float *outputTable, imageData imageInfo, misize_t *rstart, misize_t* rcount,   int startFrame, parameters  userOptions,float*  max, float *min)
{
  int thrown      = 0;
  int npoints     = 0;
  double slope    = 0.0;
  double ssd, ic;
  double duration = 0.0;
  double* voxIntegrals;
  double* voxPoints;
  double* refPoints;
  double* voxPointsAll;
  double* refPointsAll;
  char filestr[100];
  int t, x, y, z;
  int *tableIndex, index3d, index4d;
  int spatialIndex;
  int nframes = imageInfo.sizes[0];
  int nslices = imageInfo.sizes[1];
  int nrows   = imageInfo.sizes[2];
  int ncols   = imageInfo.sizes[3];
  float *integralTable; 
  float refIntegrals[imageInfo.sizes[0]];      
    float total;
    int trigger=0;
            
  /* Calculate the number of points to include based on the cutoff */
     
 //if(userOptions.RAMexceeded == 0) 
 //{ 
        integralTable=malloc(nslices*nrows*ncols*nframes*sizeof(float)); 
        tableIndex=&index4d ;
 //}
 //else
// {
 //    integralTable=malloc(nframes*nrows*ncols*sizeof(float) );
 //    tableIndex=&index3d;
// }

 //calcTACIntegrals(tacTable, integralTable,  &imageInfo, 0, userOptions.RAMexceeded, 0); 
     calcTACIntegrals(tacTable, integralTable,  &imageInfo, 0, userOptions.RAMexceeded, 0); 
/* Calculate the number of points to include based on the cutoff */
 calcTACIntegrals(imageInfo.referenceArray, refIntegrals, &imageInfo, 1, userOptions.RAMexceeded, 0 );         

for (t=0; t<nframes; t++)
  { 
    if ( t-startFrame >= 0.0  )
    { 
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
            (int) startFrame,npoints);
    exit(EXIT_FAILURE);
  }

  /* Allocate the memory for the logan plot data points */
  voxIntegrals = malloc((nframes) * sizeof(double));
  voxPoints    = malloc((npoints) * sizeof(double));
  refPoints    = malloc((npoints) * sizeof(double));
  voxPointsAll = malloc((nframes) * sizeof(double));
  refPointsAll = malloc((nframes) * sizeof(double));

  if(voxIntegrals == NULL ||
     voxPoints    == NULL || refPoints    == NULL ||
     voxPointsAll == NULL || refPointsAll == NULL )
  {
    fprintf(stderr,"Unable to allocate enough memory for Logan plot!\n");
    exit(EXIT_FAILURE);
  }
  

//for(int i=0; i<nframes;i++) printf("imageInfo.arterialTAC[%d]=%f, integral=%f\n",i, imageInfo.referenceArray[i], refIntegrals[i]  ) ; 
//if(userOptions.testStartTime==TRUE)  startFrame= testStartTime( tacTable,  outputTable, imageInfo, rstart,  rcount, userOptions, 0, refIntegrals, integralTable );    
 /* Iterate through the voxel data, slice, row, then column at a time... */
 for(z=0; z<nslices; z++){    
    for (y=0; y<nrows; y++) {
      for (x=0; x<ncols; x++) {    
		  spatialIndex = ((nrows*ncols)*z)+(ncols*y)+x;
        /* If we are inside the brain mask but outside the ref. region... */ 
	    if (imageInfo.brainmask[spatialIndex] == 1 && imageInfo.referencemask[spatialIndex] == 0){
          /* Set up the x-axis points using the ref. region data => Intgrl(REF)/REF */
          /* Iterate through each frame, excluding the last...   */
            total=0;
                for (t=0; t<nframes; t++){
                    index4d = ((nrows*ncols*nslices)*t) +  ((nrows*ncols)*z) +(ncols*y)+x;
                    total +=  tacTable[*tableIndex];
                    //printf("%f\n", tacTable[index4d]);
                        if (tacTable[*tableIndex] == 0.0 || imageInfo.referenceArray[t] == 0.0)  { 
                           refPointsAll[t] = 0.0; 
                           voxPointsAll[t] = 0.0;
                        } 
            		    else {
                             voxPointsAll[t]=tacTable[*tableIndex]/imageInfo.referenceArray[t];   
            			     refPointsAll[t]=refIntegrals[t]/imageInfo.referenceArray[t];
                        }
            		    if (t >= thrown) { 
                            refPoints[t-thrown]=refPointsAll[t];  
                            voxPoints[t-thrown]=voxPointsAll[t];
                        } 
			        /* Store ratio if we are past the points thrown away*/
		            }
	/* Store ratio if we are past the points thrown away*/
	regression_choice(/* inputs */ refPoints, voxPoints, npoints,/* outpus */ &slope, &ic, &ssd, userOptions); /* Calc slope of the linear regression and store it */
    if(userOptions.CMRg == 1) slope = slope * userOptions.PG / userOptions.LC;
    outputTable[spatialIndex]=(float) slope;
    //printf("%f\n", slope);
    if( *max < slope ) *max=slope;
	if( *min > slope ) *min=slope;
        float time=0; 
	    //if( isinf(slope) ||  isnan(slope) )   
        if(slope > 0 && trigger ==0){
			printf("time\tREF\tREF INT\tROI\tROI INT\tX-AXIS\tY-AXIS\n");

            trigger++;
            for(t=0; t<nframes; t++ ){
                index4d = ((nrows*ncols*nslices)*t) +  ((nrows*ncols)*z) +(ncols*y)+x;
                time += imageInfo.timeWidths[t];
                printf("%f\t%f\t%f\t%f\t%f\t%f\t%f\n", imageInfo.timeWidths[t], imageInfo.referenceArray[t], refIntegrals[t], tacTable[index4d], integralTable[index4d],  refPointsAll[t],   voxPointsAll[t] );
            }
            printf("Slope=%f\n", slope);


		} 
    } 
	else  	outputTable[spatialIndex]=0.0;   /* Ignore everthing else */  
   }
  }
 }
  

  /* Clean up memory allocation */
  free(voxIntegrals);
  free(voxPoints);
  free(voxPointsAll);
  free(refPointsAll);
}


 /* ----------------------------- MNI Header -----------------------------------
@NAME       : //testStartTime   
-------------------------------------------------------------------------------*/   

static float testStartTime(float* tacTable, float*  outputTable, imageData imageInfo, misize_t rstart[4], misize_t rcount[4], parameters  userOptions, float *k2, float *refIntegrals, float* integralTable )   
{
  int thrown      = 0;
  int npoints     = 0;
  double slope    = 0.0;
  double ssd, ic;
  double duration = 0.0;
  double* voxIntegrals;
  double* voxPoints;
  double* refPoints;
  double* voxPointsAll;
  double* refPointsAll;
  char filestr[100];
  int t, x, y, z;
  int *tableIndex;
  int index3d, index4d;
  int spatialIndex;
  int nframes = imageInfo.sizes[0];
  int nslices = imageInfo.sizes[1];
  int nrows   = imageInfo.sizes[2];
  int ncols   = imageInfo.sizes[3];
  float startFrame;
  char userAnswer;

   if(userOptions.RAMexceeded == 0) tableIndex=&index4d;
   else tableIndex=&index3d;
   
  //calcTACIntegrals(tacTable, integralTable,  &imageInfo, 0, userOptions.RAMthreshold, 0); 
  /* Calculate the number of points to include based on the cutoff */
    //calcTACIntegrals(imageInfo.referenceArray, refIntegrals, &imageInfo, 1, userOptions.RAMexceeded, 0 );         

  /* Allocate the memory for the logan plot data points */
  voxIntegrals = malloc((nframes) * sizeof(double));
  voxPoints    = malloc((nframes) * sizeof(double));
  refPoints    = malloc((nframes) * sizeof(double));
  voxPointsAll = malloc((nframes) * sizeof(double));
  refPointsAll = malloc((nframes) * sizeof(double));
  
  if(voxIntegrals == NULL ||
     voxPoints    == NULL || refPoints    == NULL ||
     voxPointsAll == NULL || refPointsAll == NULL )
  {
    fprintf(stderr,"Unable to allocate enough memory for Logan plot!\n");
    exit(EXIT_FAILURE);
  }
  /* Iterate through the voxel data, slice, row, then column at a time... */
   for (z=nslices-1; z != 0; z--)  
    {    
	    if(userOptions.RAMexceeded == 1) {                      
                rstart[1]=z; 
                loadpetimage( z, imageInfo.sizes, rstart, rcount, imageInfo.smoothedImageFile, tacTable, imageInfo.brainmask, imageInfo.referencemask, userOptions);
                calcTACIntegrals(tacTable, integralTable,  &imageInfo, 0, userOptions.RAMexceeded, 0);  
        }      
        

        for (y=nrows-1; y != 0; y--)   
        {
            for (x=ncols-1; x !=0; x--)
            {  
            /* Ingore index 0 of loganTable so that we can map to tables directly */
		    spatialIndex = ((nrows*ncols)*z)+(ncols*y)+x;
            /* If we are inside the brain mask but outside the ref. region... */ 
	            if (spatialIndex == nrows*ncols*((int)nslices/2) + ((int)ncols*nrows/2) + ((int) ncols/2)   )  
                {
                    while(userAnswer != 'n' || userAnswer != 'N')
                    {    
                        printf("Please Enter a Start Time: ");
                        scanf("%f", &startFrame);
                          
                        for (t=0; t<nframes; t++)
                                    { 
                                        if ( t -startFrame  >= 0.0 )
                                        {
                                        npoints=nframes-t;
                                        break; 
                                        }
                                    }        

                /* Set up the x-axis points using the ref. region data */
                /* Iterate through each frame, excluding the last...   */
         	            for (t=0; t<nframes; t++)
         	            {
		                    index3d = ((nrows*ncols)*t) + (ncols*y) + x;
                            index4d = ((nrows*ncols*nslices)*t) +  ((nrows*ncols)*z) +(ncols*y)+x; 
                            if (tacTable[*tableIndex] == 0.0) { 
                                refPointsAll[t] = 0.0; 
                                voxPointsAll[t] = 0.0;
                            } 
            		        else{    
                                if(userOptions.analysis_type=0) /*Do Logan*/
                                {
                                    voxPointsAll[t]=integralTable[*tableIndex]/tacTable[*tableIndex];  
            			            if (userOptions.arterial==TRUE) refPointsAll[t]=(refIntegrals[t]/tacTable[*tableIndex]); 
                                    else  refPointsAll[t]=(refIntegrals[t]+(imageInfo.referenceArray[t]/ *k2))/tacTable[*tableIndex]; 	         
                                }
                                else if (userOptions.analysis_type=1) /*Do Patlak*/
                                {
                                voxPointsAll[t]=tacTable[*tableIndex]/imageInfo.referenceArray[t];   
            			        refPointsAll[t]=refIntegrals[t]/imageInfo.referenceArray[t];       
                                }
                            }

            		       
                                refPoints[t-thrown]=refPointsAll[t];  
                                voxPoints[t-thrown]=voxPointsAll[t];
                             
			            /* Store ratio if we are past the points thrown away*/
		               }							
	         

                    regression_choice(/* inputs */ refPoints, voxPoints, npoints,/* outpus */ &slope, &ic, &ssd, userOptions); /* Calc slope of the linear regression and store it */
                   
		            printf("Time\t\tReference TAC\tVoxel TAC\t\tROI Int.\tRef. Int.\tY-axis\tX-axis\n");

                    float time=0;
                    for(t=0; t<nframes; t++ ){
                        index4d = ((nrows*ncols)*t)+ z*nrows*ncols +(ncols*y)+x;
                        time += imageInfo.timeWidths[t];
                        printf("%f\t%f\t%f\t%f\t%f\t%f%f\t%f\n", time, imageInfo.referenceArray[t], integralTable[t], refIntegrals[t], tacTable[index4d], refPointsAll[t],   voxPointsAll[t] );
                    }
                    printf("Number of points included: %d\n ", npoints); 
			        printf("Slope: %f\tBP: %f\n\n\n", slope, slope -1 ); 

                    printf("Do you wish to use %f as your start time? If not we will try a different start time. (y/Y/n/N)\n", startFrame);                      
                    scanf("%s", &userAnswer);
                    
                    if(userAnswer == 'y' || userAnswer == 'Y') break;

			        }  
                }
             } 
         }
    }

 //printf("Slice %d:\tErrors=%d\tMax=%f\tMin=%f\n", z, errors, *max, *min );
  /* Clean up memory allocation */
  free(voxIntegrals);
  free(voxPoints);
  free(voxPointsAll);
  free(refPointsAll);
  free(refPoints);   
  
     return startFrame;
}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : //loganPlot
@INPUT      : tacTable - the 2 dimensional TAC lookup table
              integralTable - the 2 dimensional integral table
              vol4d - the 4d volume data structure
              imageInfo.referencemask - the image mask used for this TAC
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

void loganPlot(float* tacTable, char* petimage_file, float*  outputTable, imageData imageInfo, misize_t rstart[4], misize_t rcount[4], int startFrame, parameters  userOptions, float *max ,float*  min, float *k2) 
{
  int thrown      = 0;
  int npoints     = 0;
  double slope    = 0.0;
  double ssd, ic;
  double duration = 0.0;
  double* voxIntegrals;
  double* voxPoints;
  double* refPoints;
  double* voxPointsAll;
  double* refPointsAll;
  char filestr[100];
  int t, x, y, z;
  int *tableIndex;
  int index3d, index4d;
  int spatialIndex;
  int nframes = imageInfo.sizes[0];
  int nslices = imageInfo.sizes[1];
  int nrows   = imageInfo.sizes[2];
  int ncols   = imageInfo.sizes[3];
  float* integralTable; 
  float refIntegrals[imageInfo.sizes[0]];
  int trigger=0;

        integralTable=malloc(nslices*nrows*ncols*nframes*sizeof(float)); 
        tableIndex=&index4d ;

    calcTACIntegrals(tacTable, integralTable,  &imageInfo, 0, userOptions.RAMexceeded, 0); 
  /* Calculate the number of points to include based on the cutoff */
   calcTACIntegrals(imageInfo.referenceArray, refIntegrals, &imageInfo, 1, userOptions.RAMexceeded, 0 );         
 

for (t=0; t<nframes; t++)
  { 
    if ( t -startFrame  >= 0.0 )
    {
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
            (int) startFrame, npoints);
    exit(EXIT_FAILURE);
  }

  /* Allocate the memory for the logan plot data points */
  voxIntegrals = malloc((nframes) * sizeof(double));
  voxPoints    = malloc((npoints) * sizeof(double));
  refPoints    = malloc((npoints) * sizeof(double));
  voxPointsAll = malloc((nframes) * sizeof(double));
  refPointsAll = malloc((nframes) * sizeof(double));
  if(voxIntegrals == NULL ||
     voxPoints    == NULL || refPoints    == NULL ||
     voxPointsAll == NULL || refPointsAll == NULL )
  {
    fprintf(stderr,"Unable to allocate enough memory for Logan plot!\n");
    exit(EXIT_FAILURE);
  }
                
  /* Iterate through the voxel data, slice, row, then column at a time... */


    for (z=0; z < nslices; z++){   
        for (y=0; y < nrows; y++){
            for (x=0; x < ncols; x++){  

                /* Ingore index 0 of loganTable so that we can map to tables directly */
                spatialIndex = ((nrows*ncols)*z)+(ncols*y)+x;
                /* If we are inside the brain mask but outside the ref. region... */ 

                if (imageInfo.brainmask[spatialIndex] == 1 && imageInfo.referencemask[spatialIndex] == 0){
                /* Set up the x-axis points using the ref. region data */
                /* Iterate through each frame, excluding the last...   */
         	        for (t=0; t<nframes; t++){
                    index4d = ((nrows*ncols*nslices)*t) +  ((nrows*ncols)*z) +(ncols*y)+x; 
                    
 
                        if (tacTable[*tableIndex] == 0.0) { 
                           refPointsAll[t] = 0.0; 
                           voxPointsAll[t] = 0.0;
                        } 
            		    else{    
                            voxPointsAll[t]=integralTable[*tableIndex]/tacTable[*tableIndex];  
            			    if (userOptions.arterial==TRUE) refPointsAll[t]=(refIntegrals[t]/tacTable[*tableIndex]); 
				            else  refPointsAll[t]=(refIntegrals[t]+(imageInfo.referenceArray[t]/ *k2))/tacTable[*tableIndex]; 	         
                        }
            		    if (t >= thrown){ 
                            refPoints[t-thrown]=refPointsAll[t];  
                            voxPoints[t-thrown]=voxPointsAll[t];
                        } 
		            }							
	             
                regression_choice(/* inputs */ refPoints, voxPoints, npoints,/* outpus */ &slope, &ic, &ssd, userOptions); /* Calc slope of the linear regression and store it */
                
               if( slope > 6 && trigger==0){
               //if( z == 154 && y == 135 &&  x == 53 ){
                    printf("Index=%d\n", spatialIndex); 
	                printf("Time\t\tReference TAC\tRef. Int.\tVoxel TAC\tVoxel Int.\tX-axis\t\tY-axis\n");
                    float time=0;
                    trigger=1;
                    for(t=0; t<nframes; t++ ){
                        index4d = ((nrows*ncols*nslices)*t) +  ((nrows*ncols)*z) +(ncols*y)+x;
                        time += imageInfo.timeWidths[t];
                        printf("%f\t%f\t%f\t%f\t%f\t%f\t%f\n", time, imageInfo.referenceArray[t], refIntegrals[t], tacTable[index4d], integralTable[index4d],  refPointsAll[t],   voxPointsAll[t] );
                    }
                    printf("Slope=%f\n", slope);
                }

                if(userOptions.bp==1) slope -= 1; 
                outputTable[spatialIndex]= (float) slope;
                //printf("%f\t", slope);
               
    	        if( *max < slope && isinf(slope) == FALSE ) *max=slope;
		        if( *min > slope && isinf(slope) == FALSE) *min=slope;

		           /* if( isinf(slope) ||  isnan(slope) )   */
                    /* if(   ((nrows*ncols)*z) +(ncols*y)+x == 400841  )   
		            {printf("%d %d %d\n", z, y, x);
		            printf("time\trefIntegrals\timageInfo.referenceArray\t\trefPointsAll\ttacTable\tvoxPointsAll\n");
			            for(t=0; t<nframes; t++ )
			            {
			            index4d = ((nrows*ncols)*t)+ z*nrows*ncols +(ncols*y)+x;
				        printf("%d\t\t%f\t\t%f\t\t%f\t\t%f\t\t%f\n", t,refIntegrals[t], imageInfo.referenceArray[t], refPointsAll[t], tacTable[index4d],  voxPointsAll[t] );
			            
                        }
                    printf("npoints: %d\n ", npoints); 
			        printf("slope: %f\n\n\n", slope );
			        }*/  
               } 
               else  outputTable[spatialIndex]=0.0;   /* Ignore everthing else */    
   }
  }
 }  
 //printf("Slice %d:\tErrors=%d\tMax=%f\tMin=%f\n", z, errors, *max, *min );
  /* Clean up memory allocation */
  free(voxIntegrals);
  free(voxPoints);
  free(voxPointsAll);
  free(refPointsAll);
  free(integralTable);
}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : regression_choice
@INPUT      :
@OUTPUT     : 
@RETURNS    : (nothing)
@DESCRIPTION: 
@METHOD     :
@GLOBALS    :
@CALLS      :
@CREATED    : 
@MODIFIED   : 
---------------------------------------------------------------------------- */


static int regression_choice(double *x, double *y, int nr, double *slope, double *ic, double *ssd, parameters user_options)
{

if (user_options.regression_type=="P") {llsqperp3(x,y, nr, slope, ic, ssd);}
/* add more regression options here...*/
else *slope=calcLinearRegression(nr,x,y); 

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
static float calcLinearRegression(int n, double* x, double* y)
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

	

  return (float)top/bottom;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : roi
@INPUT      : tacTable - the 2 dimensional TAC lookup table
              vol4d - the 4d volume data structure complete with time info
              imageInfo.referencemask - the image mask used for this TAC
              brainMask - the brain mask used for this TAC
@OUTPUT     : (nothing)
@RETURNS    : Produces the SRTM plots
@DESCRIPTION: Calculates the average and standard deviation for an ROI

@METHOD     :
@GLOBALS    :
@CALLS      :
@CREATED    : Thomas Funck
@MODIFIED   :
---------------------------------------------------------------------------- */


void roi(misize_t *sizes, float* tacTable, float* referencemask, parameters userOptions)
{
int x, y,z, t=0;
int index; 
int nframes=sizes[0];
int nslices=sizes[1];
int nrows=sizes[2];
int ncols=sizes[3];
float average;
int n=0;
	      
  for(z=0; z<nslices; z++)
  {    
	for(y=0; y<nrows; y++)
	{
		for (x=0; x<ncols; x++)
		{
			index = ((nrows*ncols)*z) +(ncols*y)+x;	
            //printf("%f\n", referencemask[index] );
				if (referencemask[index] == 1.0) 
				{
				average += tacTable[index]; 
				n += 1;
			   
                printf("%d\t%f\t%d\n", index, tacTable[index], n );
				}
			}
		}	
	}
  
   
	printf("\n\nFrame\tAverage\t\tTotal\t\tNpoints\n"); 
	/*Print out the averages...*/
		for(t=0; t<nframes; t++)
		{
		printf("%d\t%f\t%f\t%d\n", n, average/  n, average, n );
		}
	printf("\n\n");
   

}
