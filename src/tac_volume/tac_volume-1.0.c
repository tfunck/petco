#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <omp.h>
#include "minc2.h"
#include "hdf5.h"
#include "minc_helper.h"
#include "pthread.h"

void useage();

struct roi_tac{
    unsigned char* roi;
    float** tac;
    int nroi;
    int nframes;
};

struct args{
    unsigned char* image_vol;
    float* new;
    unsigned char* roi;
    int nroi;
    float** tac;
    float* min;
    float* max;
    unsigned long xmax; 
    unsigned long ymax; 
    unsigned long zmax; 
    int thread;
    int nthreads;
    int time;
};


char *inputString(FILE* fp, size_t size){
//The size is extended by the input with the value of the provisional
    char *str;
    int ch;
    size_t len = 0;
    str = realloc(NULL, sizeof(char)*size);//size is start size
    if(!str)return str;
    while(EOF!=(ch=fgetc(fp)) && ch != '\n'){
        str[len++]=ch;
        if(len==size){
            str = realloc(str, sizeof(char)*(size+=16));
            if(!str)return str;
        }
    }
    str[len++]='\0';
    printf("Length: %d\n", len);
    return realloc(str, sizeof(char)*len);
}

int is_in(void* value, void* list, int n, int dtype ){
    unsigned int i=0;
    if( dtype == MI_TYPE_INT){ 
        for(i=0; i< n; i++){ 
            if( *((int*) value) == ((int*) list)[i]) return(i);
        }
    }
    else if( dtype == MI_TYPE_FLOAT){ 
        for(i=0; i< n; i++){ 
            if( *((float*) value) == ((float*) list)[i]) return(i);
        }
    }
    else if( dtype == MI_TYPE_UBYTE){
        //printf("UBYTE\n");
        for(i=0; i< n; i++){ 
            if( (int) *((unsigned char*) value)  == ((unsigned char*) list)[i] - '0' ){
                //printf("%d,%d\t%d\n",i,(int) *((unsigned char*) value), ((unsigned char*) list)[i] - '0'  );
                return(i);
            }
        }
    }
    else if( dtype ==  MI_TYPE_DOUBLE){ 
        for(i=0; i< n; i++){ 
            if( *((double*) value) == ((double*) list)[i]) return(i);
        }
    }
    //printf("%d\t%d\t%d\t%d\n", *((unsigned char*) value), i, dtype, MI_TYPE_UBYTE);
    return(-1);
}

char *getStrFromFile(FILE *fpin){
    int i=0, ch; 
    int chunk=sizeof(char);
    int size=2*chunk;
    if(EOF==(ch=fgetc(fpin)))return NULL;
    ungetc(ch, fpin);
    char *line = malloc(size);

    while(EOF!=(ch=fgetc(fpin)) && ch != '\n'){//Do not include '\n'
        line[i] = ch;
        i++;
        if(i==size-1){
            size += chunk;
            line = realloc(line, size);
        }
    }
    line[i] = '\0';
    //for(int j=0; j<i; j++) free(line[j]);
    return line;
}


struct roi_tac* read_csv(char* filename ){
    int i, j, t; 
    FILE* file=fopen(filename, "rt");
    char* buffer; //=inputString(file, 100);
    char *dlm=",";
    int nroi=0;
    unsigned char* roi=NULL;
    float** tac=NULL; // 2D array for time [0] and radioactivity [1]
    float value, time;
    char* token;
    //char** buff_list=NULL;

    struct roi_tac* output=malloc(sizeof(*output));

//    while(fgets(buffer, sizeof(buffer), file) != NULL){
    while( NULL != (buffer=getStrFromFile(file))  ){
        nroi++;
        tac=realloc(tac, nroi * sizeof(*tac));
        roi=realloc(roi, nroi * sizeof(*roi));
        tac[nroi-1]=NULL;
        token=strtok(buffer, dlm);
        
        //buff_list=realloc(buff_list, nroi * sizeof(*buff_list) );
        //buff_list[nroi-1]=buffer; 
        roi[nroi-1]=*token;

        printf("%c:\t", roi[nroi-1]); 
        token=strtok(NULL, dlm);
        t=0;
        while( token != '\0'){
            t++;
            tac[nroi-1]=realloc(tac[nroi-1], t * sizeof(**tac) );
            tac[nroi-1][t-1]=atof(token);
            printf("%3.1f ", tac[nroi-1][t-1]);
            token=strtok(NULL, dlm);
        }
        printf("\n");
        free(buffer);
    }
    
    //for(int i=0; i<nroi; i++) free(buff_list[i]);
    fclose(file);
    output->nframes=t;
    output->nroi=nroi;
    output->roi=roi;
    output->tac=tac;
    return(output);
}

int map(unsigned char val, unsigned char* list ){
    int i;
    for( i=0; list[i] - '0' !=  val; i++);
    //printf("%d\t%d\n", i, (int) val);
    return(i);
}

int main(int argc, char** argv){
    if(argc != 4 || strcmp(argv[1], "-help") == 0 ){
        useage();
    }

    data image, mask;

    //readVolume(data* volume, int read_hyperslab /*1 for read, o.w. don't read*/, mitype_t dtype)
    unsigned char* image_vol;
    int nthreads=sysconf(_SC_NPROCESSORS_ONLN);
    float* new; 
    image.filename=argv[1];
    char* tac_csv=argv[2];
    char* output_filename=argv[3];
    struct roi_tac *in = read_csv(tac_csv);
    unsigned char* roi=in->roi;
    float** tac=in->tac;
    int nroi=in->nroi;
    unsigned long tmax, zmax, ymax, xmax, ndim, nvox;
    float min=0, max=0;
    misize_t *wstarts=NULL, *wcount=NULL, *wcount_ptr=NULL;
    mihandle_t img;
    image_vol=(unsigned char*)readVolume(&image, 1, MI_TYPE_UBYTE);
    if(image_vol == NULL) pexit("Could not allocate input array. Image too large, not enough memory.\n","", 1);
    zmax=image.zmax;
    ymax=image.ymax;
    xmax=image.xmax;
    nvox=zmax*ymax*xmax;
    tmax=in->nframes;
    if(tmax > 1){ 
        ndim=4;
    }
    else{
        ndim=3;
    }

    wstarts=malloc(ndim * sizeof(*wstarts));
    wcount=malloc(ndim * sizeof(*wcount));

    if( ndim == 4){
        wcount[0]=tmax; 
        wcount[1]=zmax; 
        wcount[2]=ymax; 
        wcount[3]=xmax;
        wstarts[0]=wstarts[1]=wstarts[2]=wstarts[3]=0;
    }
    else if( ndim == 3){
        wcount[0]=zmax; 
        wcount[1]=ymax; 
        wcount[2]=xmax;
        wstarts[0]=wstarts[1]=wstarts[2]=0;
    }

    for(int i=0; i<nroi; i++) for(int t=0; t<tmax; t++) if(max < tac[i][t]) max=tac[i][t];

    new=calloc( nvox, sizeof(*new));

    if(new == NULL) pexit("Could not allocate memory for array of new image", "", 1);

    createVolume(output_filename, ndim, wcount, image.step, image.start, MI_TYPE_FLOAT);
    if(ndim==4) wcount[0]=1;
    miopen_volume(output_filename, MI2_OPEN_RDWR, &img);

    for( int t=0; t<tmax; t++) {
        if(ndim==4) {wstarts[0]=t;}
        for(unsigned long z=0; z< zmax; z++){
            for(unsigned long y=0; y< ymax; y++){
                for(unsigned long x=0; x < xmax; x++){
                    unsigned long index=z*ymax*xmax+y*xmax + x;
                    if( image_vol[index] != 0 ){
                        int i=is_in( &(image_vol[index]), roi, nroi, MI_TYPE_UBYTE);
                        if( i >= 0){
                            new[index]=tac[i][t];
                        } 
                    }

                    //new[index]=(int) image_vol[index];
                } 
            }
        }
        int status=miset_real_value_hyperslab(img, MI_TYPE_FLOAT, wstarts, wcount, new);
    }


    /**/
    miset_volume_range ( img, max,  min);
    miclose_volume(img);
    free(in);
    free(new);
    return(0);
}


void useage(){
    printf("Useage: image.mnc tac.csv tac.mnc\n");
    printf("Purpose: Put time activity curves (TAC) into a minc file based on image with integer labels\n");
    printf("Format of csv files:\n");
    printf("ROI, frame value 1, frame value2, frame value 3, ..., frame value n\n");
    printf("1, 50.0, 100.6, 400.8, 300.1, 250.0...\n");
    printf("2, 10.0, 200.0, 700.9, 400.2, 350.7...\n");
    exit(1);
}
