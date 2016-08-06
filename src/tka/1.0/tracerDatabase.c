#include <stdio.h>
#include <stdlib.h>
#include <minc2.h>
#include <math.h>
#include "mtga.h"
#include <string.h>
#ifndef TRUE
#   define TRUE 1
#   define FALSE 0
#endif

 /************
 *  Definitions
 * *********/  




int searchAutomatic(int *i, int numTracers, tracer* tracerDatabase, mihandle_t *petimage, char *tracerString, parameters* userOptions ); 

int openDatabase(tracer **tracerDatabase, mihandle_t* petimage, FILE *inputfile, char *filename  ); 

tracer* getTracerParameterValues(mihandle_t *petimage, parameters* userOptions);    

int readNewTracer( tracer* newTracer    , FILE* inputfile) ; 


int searchManual(int *i, int  numTracers,tracer* tracerDatabase,char  *tracerString); 


/******************
 *   FUNCTIONS
 * ***************/
/*int main()
{
mihandle_t petimage;
char* image_file="anno-08_2.mnc";


if ( miopen_volume( image_file, MI2_OPEN_READ, &petimage  ) != MI_NOERROR) 
    {
    fprintf(stderr,"Error Opening %s.\n", image_file );    
    exit(1);    
    }


getTracerParameterValues( &petimage);


return 0;
} */



tracer* getTracerParameterValues(mihandle_t *petimage, parameters* userOptions)
{
tracer *tracerDatabase;
tracer *outputTracer, tempoutputTracer;
char halflife_string[20], tracerString[30] ;
double halflife;  
int i=0, j;
int numTracers=0;
char *filename="databaseFile.txt";
FILE inputfile;

  numTracers=openDatabase(&tracerDatabase,petimage, &inputfile, filename);
                        //remember to close the file once we are done using it!!!
  if ( searchAutomatic(&i, numTracers, tracerDatabase, petimage, tracerString, userOptions) == 0) {outputTracer=&(tracerDatabase[i]); } /*check if we can find a tracer name that matches the tracer in the header*/
  else if ( searchManual(&i, numTracers, tracerDatabase, tracerString) == 0 ) { outputTracer=&(tracerDatabase[i]);  }
  else if ( saveNewTracer(&inputfile, filename, &outputTracer) != 0 ) 
    { 
        printf("Error: No tracer parameters were entered.\n"); 
        exit(1);
    }    
   

  miget_attr_values(*petimage, MI_TYPE_STRING, "/ecat-main/", "Isotope_Halflife", 20,  halflife_string);  
  halflife=atof(halflife_string);
  outputTracer->lambda=(log(2)/halflife);
  
  
  printf("\nWARNING!\nUsing the following parameters:\nName\t\t%s\nk2 Mean\t\t%f\nk2 Min\t\t%f\nk2 Max\t\t%f\nDecay Constant\t\t%f\n",outputTracer->name, outputTracer->k2mean, outputTracer->k2min, outputTracer->k2max, outputTracer->lambda );
  
//  for(j=0; j< numTracers; j++) if(j != i ) free(tracerDatabase[j]); /*free the tracers we are not using*/
  //free(tracerDatabase);
  
 

  return outputTracer;
}


int openDatabase(tracer **tracerDatabase, mihandle_t* petimage, FILE *inputfile, char *filename  )
{
int i=0, j;
float newLambda;
float newK2max;
float newK2min;
char newTracerName;


inputfile=fopen(filename, "r+");  

*tracerDatabase=malloc(sizeof(tracer));
   
    while( readNewTracer(  &((*tracerDatabase)[i]), inputfile ) != 1) 
    {      
    i++;
    *tracerDatabase=realloc(*tracerDatabase, (i+1)* sizeof(tracer));//Check this, no idea if this is how its done
    }
fclose(inputfile);
return i;   
}


int readNewTracer( tracer* newTracer    , FILE* inputfile)
{
char buffer[30];
char* dlm="";
char* token;
int   counter;
int i;
char halflife[100];


    if  ( fgets(buffer, sizeof(buffer), inputfile) != NULL )
    {   
        for(i=0;i < sizeof(buffer); i++) if(buffer[i]== '\n') buffer[i]='\0';
        strcpy(newTracer->name, buffer);
        printf("%s\t", newTracer->name );
         
 
        fgets(buffer, sizeof(buffer), inputfile);
        newTracer->k2mean=atof(buffer);   
        printf("%f\t", newTracer->k2mean );
    
        fgets(buffer, sizeof(buffer), inputfile); 
        newTracer->k2min=atof(buffer);           
        printf("%f\t", newTracer->k2min  );
        
        fgets(buffer, sizeof(buffer), inputfile); 
        newTracer->k2max= (double) atof(buffer); 
        printf("%f\n",newTracer->k2max);  
        fgets(buffer, sizeof(buffer), inputfile);   
    }
    else return 1;
 
  return 0;  
}

int checkTracer(char *tracerString)
{
char userAnswer;

    
  printf("Do you wish to use %s as your tracer? (y/Y/n/N)\n", tracerString); 
  scanf("%c", &userAnswer);
      if(userAnswer == 'y' || userAnswer == 'Y')   return 0;
      else if (userAnswer == 'n' || userAnswer == 'N') return 1;  
      else return checkTracer(tracerString);
}


int searchAutomatic(int *i, int numTracers, tracer* tracerDatabase, mihandle_t *petimage, char *tracerString, parameters* userOptions )
{
 
 miget_attr_values(*petimage, MI_TYPE_STRING,"/acquisition/","tracer", 30, tracerString);        
 for(*i=0; *i<numTracers; (*i)++ )   if(strcmp(tracerString, tracerDatabase[*i].name) == 0) 
                                    { 
                                           if(userOptions->autok2 == 0) return checkTracer(tracerString);
                                           else return 0;
                                    }   

 return 1;
          
}

int searchManual(int *i, int  numTracers,tracer* tracerDatabase,char  *tracerString) 
{
int k;

printf("Tracer was not found automatically, we will try to find it manually from a list of tracers. \n");
for(k=0; k<numTracers; k++) printf("%d\t%s\n", k,  tracerDatabase[k].name);
printf("%d\tTracer Not Listed\n", k);
*i=k+1;    
    while(*i>k || *i < 0 )
    {
        printf("Please enter the number of the option you wish to select: ");
        scanf("%d", i);
        printf("\n");
    }

    if(*i != k) return checkTracer(tracerDatabase[*i].name);
        else return 1;
}

int saveNewTracer(FILE *inputfile, char* filename, tracer**  outputTracer)
{
char userAnswer;
char newname[30];
double newk2mean;
double newk2min;
double newk2max;
char outVar[250];
int i;

*outputTracer=malloc(sizeof(tracer));

inputfile=fopen(filename, "r+");      
userAnswer='a';  

printf("No tracer was selected from the tracer database.\n"); 
printf("Would you like to add a new tracer to this database? (y/Y/n/N)\n"); fflush(stdout);  
     
    while( (userAnswer == 'y' || userAnswer == 'Y' || userAnswer == 'n' || userAnswer == 'N') == FALSE) scanf("%c", &userAnswer);
 
    if (userAnswer == 'n' || userAnswer == 'N') return 1;
   
userAnswer='a';      

while( (userAnswer == 'y' || userAnswer == 'Y' ) == FALSE ) 
    {   fflush(stdin);
    printf("Name of New Tracer: ");//fflush(stdout); fflush(stdin);
    getchar();
    fgets( newname, 30, stdin );  //scanf("%s", newname);
    for(i=0; i<30; i++) if(newname[i] == '\n') newname[i] = '\0'; 
    printf("\nMean k2 of %s:\n", newname);
    scanf("%lf", &newk2mean);   
    
    printf("Max k2 of %s:\n", newname);
    scanf("%lf", &newk2max );   

    printf("Min k2 of %s:\n", newname);
    scanf("%lf", &newk2min);

    printf("The following tracer and parameters will be added to the tracer database:\nName:\t%s\nMean k2\t%f\nMax k2\t%f\nMin k2\t%f\n", newname, newk2mean, newk2max, newk2min);

       printf("Is this correct? (y/Y/n/N)\n");    userAnswer='a';   
       while( (userAnswer == 'y' || userAnswer == 'Y' || userAnswer == 'n' || userAnswer == 'N') == FALSE) scanf("%c", &userAnswer);

    }  ; 

    sprintf(outVar,"\n\n%s\n%f\n%f\n%f\n", newname, newk2mean, newk2min, newk2max);
    fseek(inputfile, -1, SEEK_END);
    fputs(outVar, inputfile);
    printf("%s added to tracer database. If you wish to edit or rename a tracer, simply edit %s\n", newname, filename);

(*outputTracer)->k2mean=newk2mean;
(*outputTracer)->k2min=newk2min;
(*outputTracer)->k2max=newk2max;
strcpy((*outputTracer)->name, newname);

fclose(inputfile);
return 0;

}   
