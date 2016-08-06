#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <minc2.h>
#include <hdf5.h>
#include <math.h>


int main(int argc, char** argv){
 
 mihandle_t petimage;
 mitype_t datatype;
 char string[100];
 int length;
 char* attribute="time";
 char* path="/data";


printf("Opening %s\n", argv[1]);
 if ( miopen_volume(argv[1], MI2_OPEN_RDWR, &petimage) != MI_NOERROR) printf("Problem Opening File\n");
 //if ( micreate_group(petimage, "/helloworld", "test") != MI_NOERROR) printf("Problem creating group\n");
 
    //SET VALUES
 //if ( miset_attr_values ( petimage, MI_TYPE_STRING, "Hello/World/", "Hiya", sizeof(char)*3, values) != MI_NOERROR) printf("Problem creating group2\n");


 miget_attr_length(petimage, path, attribute, &length );
 printf("Length: %d\n", length);
 
 miget_attr_type(petimage, path, attribute, &datatype );
 if(datatype == MI_TYPE_STRING) printf("String\n");
 if(datatype == MI_TYPE_DOUBLE) printf("Double\n");
 if(datatype == MI_TYPE_FLOAT) printf("Float\n");
 if(datatype == MI_TYPE_INT) printf("Int\n");
 if(datatype == MI_TYPE_SHORT) printf("Short\n");
 if(datatype == MI_TYPE_UBYTE) printf("Ubyte\n");
 printf("Size: %d\n", (int) sizeof(datatype));

 if(datatype == MI_TYPE_STRING) {
    char data[length];
    miget_attr_values(petimage, datatype, path, attribute, length,  data); 
    printf("%s\n", data );
 }

 if(datatype == MI_TYPE_DOUBLE) { 
  double data; 
     miget_attr_values(petimage, datatype, path, attribute, length,  &data); 
     printf("%f\n",(double) data);
 }
 
 if(datatype == MI_TYPE_FLOAT) {
     float data; 
     miget_attr_values(petimage, datatype, path, attribute, length,  &data); 
     printf("%f\n", (float) data);
 }

 if(datatype == MI_TYPE_INT) { 
     int data; 
     miget_attr_values(petimage, datatype, path, attribute, length,  &data); 
     printf("%d\n", (int)  data);
 }

 return 0;
}
