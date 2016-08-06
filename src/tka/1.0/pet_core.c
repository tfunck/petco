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


static int processImage(char* petimage_file, char* referencemask_file, char* brainmask_file, char* time_widths_file, char* outputfile,  int startFrame, parameters userOptions)
{
  int t, i, tacMaxIndex;
  int nframes, nslices, nrows, ncols;
  float ***filterdata;           /* Holds the raw filter data          */
  misize_t pad_nframes,      /* Variable for the padded dimensions */
               pad_nrows,
               pad_ncols;
  mihandle_t petimage, outputimage, mask, smoothedImage, *image_ptr; 
  midimhandle_t  dimensions[4];
  imageData imageInfo;
  int  ***brainmask_array,  ***referencemask_array, z, zmax, tmax, ymax, xmax, x, y,  counter; 	
  misize_t  sizes[4], tempSizes[3];
  misize_t  location[3];
  float  *referenceArray; 
  double step[4];
  double separations[4];
  double starts[4];
  //float* timeWidths;
  //float *times;
  int *brainmask;
  float *referencemask;
  float *reference_ptr;

/* Load the 4d image file into a Volume4d data structure */
printf("Starting processing of the input file '%s'...\n",petimage_file);

	if(userOptions.input4d==1) {
	    getdimensions(petimage_file, &petimage, 4, sizes, dimensions, step, separations, starts);
        miclose_volume(petimage); 
	 }	
	else{
	    getdimensions(petimage_file, &petimage, 3, tempSizes, dimensions, step, separations, starts);
	    //miclose_volume(petimage); 
	    sizes[0]=1;
	    sizes[1]=tempSizes[0];
	    sizes[2]=tempSizes[1];
	    sizes[3]=tempSizes[2];	
	}

imageInfo.sizes=sizes;
nframes=sizes[0]; 
nslices=sizes[1]; 
nrows=sizes[2]; 
ncols=sizes[3]; /*makes it easier to read the code if the sizes are stored as individually named variables*/


loadmask_hyperslab(brainmask_file, &(imageInfo.brainmask));
    if(userOptions.input4d==TRUE){
        float* timeWidths=malloc(nframes*sizeof(float));
	    getTimeWidths(nframes, startFrame, timeWidths, time_widths_file, userOptions); 
 	    imageInfo.timeWidths=timeWidths;
    	imageInfo.times=calloc(nframes,sizeof(*imageInfo.times));
        for(int i=0; i<nframes;i++){ 
            if(i==0) imageInfo.times[i]=imageInfo.timeWidths[i]*0.5+imageInfo.times[i];
	        else imageInfo.times[i]=imageInfo.times[i-1]+0.5*(imageInfo.timeWidths[i-1]+imageInfo.timeWidths[i]);
        }
        //for(int i=0; i<nframes; i++) printf("Times: %f\n", imageInfo.times[i]);
    }
    if (userOptions.arterial == TRUE) {
	    imageInfo.referenceArray = get_arterial_TAC( nframes, imageInfo.times, referencemask_file, &userOptions);

        referencemask=malloc(nrows*ncols*nslices*sizeof(float));
	    for(i=0; i<nrows*ncols*nslices; i++) referencemask[i]=0;
	}
    else if( userOptions.analysis_type == 3   ) { 
        referencemask=malloc(nrows*ncols*nslices*sizeof(float));    
        for(i=0; i<nrows*ncols*nslices; i++) referencemask[i]=0;    
    }
	else loadmask_hyperslab(referencemask_file, &referencemask); 
	imageInfo.referencemask=referencemask; 



    if( (userOptions.analysis_type == 0 || userOptions.analysis_type==2)){ 
	    if (userOptions.autok2 == TRUE  ) imageInfo.tracerInfo=getTracerParameterValues(&petimage, &userOptions);
	    else {
	        imageInfo.tracerInfo=malloc(sizeof(tracer)); 
	        imageInfo.tracerInfo->k2mean=userOptions.k2;
	    }
    }
   
    if  (userOptions.analysis_type != 3) {   
        	createImage( outputfile, 3, sizes, separations, starts ); //Create 3D output image         
            //miclose_volume(outputimage);
    }
	else if  (userOptions.analysis_type == 3) {    
        createImage( outputfile, 4, sizes, separations, starts ); //Create 4D output image
	}   
	
    imageInfo.outputimage_file=outputfile;
	if(userOptions.arterial == FALSE && userOptions.write == TRUE && userOptions.analysis_type != 3) find_reference_values( petimage_file ,  sizes, &(imageInfo.referenceArray), imageInfo.referencemask); 
    for(int i=0; i<nframes; i++) printf("%f\t%f\n", imageInfo.times[i], imageInfo.referenceArray[i]);
    //processSlices is the key data processing step, where the data is read and analyzed according to the desired model

    processSlices(imageInfo, petimage_file, userOptions, startFrame);

 printf("\nCompleted processing of the input file '%s'.\n", petimage_file);
 return(EXIT_SUCCESS);
}

int main(int argc, char** argv)
{
  char* petimage_file=NULL;
  char* referencemask_file=NULL;
  char* brainMaskFile=NULL;
  char* outputfile=NULL;
  char* time_widths_file;
  int startFrame;
  int badparams;
  int splitframes;
  int genmask;
  int x0,y0,z0,x1,y1,z1;
  int frame_no;
  int result;
  int tempIcvId;
  int tempMincId;
  int expected_parameters=1;/*starts at 1 --> program command is always first input in terminal*/
  int counter=0;
  int arg=0;
  parameters user_options;
  user_options.RAMthreshold=0.5; /*Set default max RAM usage to 0.5GB, or 500MB*/
  user_options.CMRg=0;
  user_options.arterial=FALSE;
  user_options.autok2=FALSE;
  user_options.testStartTime=FALSE;
  user_options.arterial_convert_to_seconds=FALSE;
  user_options.conversion_factor=1;
    /* Ensure we have some parameters specified */
  badparams = (argc == 1);

  if (badparams == FALSE)
  {

    user_options.analysis_type=-1;//The -1 was changed from NULL
   
	//Check parameters for options and set expected number of parameters. 
	//This is needed so that more options can be added later and keep track of how
	//many parameters the program can expect to get.
	while (badparams == 0 && argv[counter] != NULL)
	{
	

            if (strcmp(argv[counter], "-autok2") == 0)  //check for autok2 option
			{
			expected_parameters += 1;
			user_options.autok2=1;
			} 
			
			 if (strcmp(argv[counter], "-RAM") == 0)  //change RAM threshold usage (Default == 500MB)
			{
			expected_parameters += 2;
			user_options.RAMthreshold=atof(argv[++counter]);
			} 
			
            if (strcmp(argv[counter], "-arterial") == 0) {
                if( strcmp( argv[counter+1], "-min") == 0 ) {
                    printf("Converting minutes to seconds.\n");
                    user_options.arterial_convert_to_seconds=TRUE;
                    //expected_parameters++;
                }else {printf("Not Converting\n");}
			expected_parameters += 1;
			user_options.arterial=TRUE;
			}
            if ( strcmp( argv[counter], "-CMRg") == 0 ){
                expected_parameters+=3;
                user_options.CMRg=1;
                counter += 1;
                user_options.PG=atof(argv[counter]);
                counter += 1;
                user_options.LC=atof(argv[counter]);
            }

			if (strcmp(argv[counter], "-bp") == 0)
			{
			expected_parameters += 1;
			user_options.bp=TRUE;
			}

			if (strcmp(argv[counter], "-min") == 0)
			{
			expected_parameters += 1;
			user_options.min=TRUE;
			}

  			if (strcmp(argv[counter], "-conversion_factor") == 0)
			{
			    expected_parameters += 2;
                counter += 1;
			    user_options.conversion_factor=atof(argv[counter]);
			}
            if (strcmp(argv[counter], "-testStartTime") == 0)
			{
			expected_parameters += 1;
			user_options.testStartTime=TRUE;
			}   

			if (strcmp(argv[counter], "-test") == 0)
			{
			expected_parameters += 1;
			user_options.test=TRUE;
			}
			if (strcmp(argv[counter], "-pp") == 0)
			{
			user_options.integrate=TRUE;
			user_options.input3d=FALSE;
			user_options.input4d=TRUE;
			user_options.write=TRUE;
			user_options.analysis_type=1;
			petimage_file=argv[counter+1];
 			referencemask_file=argv[counter+2];
 			brainMaskFile=argv[counter+3];
			startFrame=atof(argv[counter+4]);
            time_widths_file=argv[counter+5];
            outputfile=argv[counter+6];
			expected_parameters += 7; 
			}
			if (strcmp(argv[counter], "-lp") == 0){
			    user_options.integrate=TRUE;
			    user_options.input3d=FALSE;
			    user_options.input4d=TRUE;
			    user_options.write=TRUE;
			    user_options.analysis_type=0;
			    arg=1;
		    	//if using ref region and not detecting k2 from header, must be given after -lp
			    if (user_options.arterial == FALSE &&  user_options.autok2==FALSE){ 
			        user_options.k2=atof(argv[counter+arg++]); 		
			    }
			    petimage_file=argv[counter+arg++];
 			    referencemask_file=argv[counter+arg++];
 			    brainMaskFile=argv[counter+arg++];
			    startFrame=atof(argv[counter+arg++]);
                time_widths_file=argv[counter+arg++];
                outputfile=argv[counter+arg++];
			    expected_parameters += arg; 
			}
			if (strcmp(argv[counter], "-sr") == 0)
			{ 
			user_options.integrate=FALSE;
			user_options.input3d=FALSE;
			user_options.input4d=TRUE;
			user_options.write=TRUE;
			petimage_file=argv[counter+1];
 			referencemask_file=argv[counter+2];
 			brainMaskFile=argv[counter+3];
			startFrame=atof(argv[counter+4]);
            outputfile=argv[counter+5];
			expected_parameters += 6; 
			user_options.analysis_type=2;
			}
			if (strcmp(argv[counter], "-sv") == 0)  
			{
		   	user_options.integrate=FALSE;
			user_options.input3d=FALSE;
			user_options.input4d=TRUE;
			user_options.write=TRUE;     
 			user_options.analysis_type=3; 
			user_options.injected_dose=atof(argv[counter+1]);
            user_options.patient_weight=atof(argv[counter+2]);
			user_options.start_frame=atof(argv[counter+3]);
			user_options.end_frame=atof(argv[counter+4]);
			petimage_file=argv[counter+5];
 			brainMaskFile=argv[counter+6];
            outputfile=argv[counter+7];
            expected_parameters += 8; 

            }
			
		counter++;
	}
 }
  


printf("%d parameters are expected and %d parameters were received.\n", expected_parameters, argc);
if (expected_parameters != argc) badparams=TRUE;
	
  /* Display utsage message if we have bad parameters */
  if (badparams)
  {
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
    exit(EXIT_FAILURE);
  }



    result = processImage(petimage_file, referencemask_file, brainMaskFile, time_widths_file, outputfile,  startFrame, user_options);
  
    return(0);
}
