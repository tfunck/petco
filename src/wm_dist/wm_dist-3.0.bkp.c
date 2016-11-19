#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <malloc.h>
#include <math.h>
#include <pthread.h>
#include <unistd.h>
#include "minc2.h"
#include "minc_helper.h"
#define TRUE 1
#define FALSE 0
#define MAX 99999 //FIXME: Only makes sense for brains at the mm scale, may have to be adapted for other applications
int VERBOSE=FALSE;
const float pseudo_inf=9999999;
float dx, dy, dz, dx2, dy2, dz2;
int xmax, ymax, zmax, xymax;

struct node{
    float* dist;
    int index, z,y,x;
    struct node* left, *right;
};

void* wm_dist_threaded(void*);
void useage();

int length(char *array[]) {
    int i;
    for(i=0; array[i] != NULL; i++);
    return(i);
}

int check_input_files(char *file_inputs[]){
    int n_file_inputs=length(file_inputs);
    int i;
    if(VERBOSE) printf("Number of file inputs: %d\n", n_file_inputs);
    for(i=0; i < n_file_inputs; i++){
        if( access( file_inputs[i], R_OK ) != -1 ) {
            if(VERBOSE) printf("Can read file: %s\n", file_inputs[i]);
        } else {
            printf("Error: could not access %s\n", file_inputs[i]);
            return 1;
        }
    }
    return 0 ;
}

struct wm_vol_args{
    data* img;
    int** gm_border;
    int* img_vol;
    float* mat; 
    int label;
    int thread;
    int nthreads;
    int n;
};



int check_for_wm(int* img_vol, int z, int y, int x, int zmax, int ymax, int xmax, int WM){
    for( int i=-1; i <= 1; i++ )
        for( int j=-1; j <= 1; j++ ) 
            for( int k=-1; k <= 1; k++ ){
                int zi=z+i;
                int yi=y+j;
                int xi=x+k; 
                if( img_vol[zi*ymax*xmax+yi*xmax+xi ] == WM) return(1);
            }
    return(0);
}

int** wm_gm_border(data* img, int GM, int WM, int* img_vol, int* n){
    int zmax=img->zmax;
    int ymax=img->ymax;
    int xmax=img->xmax;
    int** border=NULL;
    *n=0;
    /*****************
    *   Find border  * 
    ******************/
    for(int z=0; z<zmax; z++ ){
        for(int y=0; y<ymax; y++ ){
            for(int x=0; x<xmax; x++ ){
                int index=z*ymax*xmax+y*xmax+x;
                int val=img_vol[index];
                if(val==GM){
                    if( check_for_wm(img_vol, z, y, x, zmax, ymax, xmax, WM) == 1){
                        *n += 1;
                        border=realloc(border, sizeof(*border) * *n);
                        border[*n-1]=malloc( sizeof(**border) * 4);
                        border[*n-1][0]=index;
                        border[*n-1][1]=z;
                        border[*n-1][2]=y;
                        border[*n-1][3]=x;
                        //printf("%d: %d %d %d\n", *n-1, z, y, x );
                    }
                } 
            }
        }
    }



    return(border);
}


void update(struct node* considered, _Bool* fixed, float* distances){
    int xpp, xmm, ypp, ymm, zpp, zmm;
    int ixp, ixpp, ixm, ixmm, iyp, iypp, iym, iymm, izp, izpp, izm, izmm;
    int xp, xm, yp, ym, zp, zm;
    float xdp, xdpp, xdm, xdmm,  ydp, ydpp, ydm, ydmm, zdm, zdmm, zdp, zdpp;
    float alist[6], blist[6], clist[6];
    float alta[3], altb[3], altc[3];
    float altdisc[3], altsols[3][2];
    float altsol[3];
    float sol, sol1, sol2;
    float alpha, disc, t, a,b, c; 
    float d=*(considered->dist); 
    int counter=0;
    int index = considered->index;
    int z=considered->z;
    int y=considered->y;
    int x=considered->x;
    xp=x+1;
    xpp=xp+1;
    xm=x-1;
    xmm=xm-1;
    yp=y+1;
    ypp=yp+1;
    ym=y-1;
    ymm=ym-1;
    zp=z+1;
    zpp=zp+1;
    zm=z-1;
    zmm=zm-1;
    ixp=z*xymax+y*xmax+xp; 
    ixpp=z*xymax+y*xmax+xpp;
    ixm=z*xymax+y*xmax+xm;
    ixmm=z*xymax+y*xmax+xmm;
    iyp=z*xymax+yp*xmax+x;
    iypp=z*xymax+ypp*xmax+x;
    iym=z*xymax+ym*xmax+x;
    iymm=z*xymax+ymm*xmax+x;
    izp=zp*xymax+y*xmax+x;
    izpp=zpp*xymax+y*xmax+x;
    izm=zm*xymax+y*xmax+x;
    izmm=zmm*xymax+y*xmax+x;
    xdp=distances[ixp];
    xdpp=distances[ixpp];
    xdm=distances[ixm];
    xdmm=distances[ixmm];
    ydp=distances[iyp];
    ydpp=distances[iypp];
    ydm=distances[iym];
    ydmm=distances[iymm];
    zdm=distances[izm];
    zdmm=distances[izmm];
    zdp=distances[izp];
    zdpp=distances[izpp];

    alist[0]=alist[1]=alist[2]=alist[3]=alist[4]=alist[5]=0;
    blist[0]=blist[1]=blist[2]=blist[3]=blist[4]=blist[5]=0;
    clist[0]=clist[1]=clist[2]=clist[3]=clist[4]=clist[5]=0;
    /*
    * Calculate new distance based on Eikonald equation
    */

    if( fixed[ixm]==TRUE) {//X-- second order
        if( fixed[ixmm]==TRUE){ 
            alpha = 9 / (4 * dx2);
            t=(4*xdm-xdmm)/3;
            alist[0] = alpha;
            blist[0] = -2*alpha*t;
            clist[0] = alpha*t*t;
        }
        else{ //X- first order
            alpha = 1/dx2;
            alist[0] = alpha;
            blist[0] = -2*xdm*alpha;
            clist[0] = xdm*xdm*alpha;
        }
    } 
    if( fixed[ixp]==TRUE) {//X++ second order
        if( fixed[ixpp]==TRUE){ 
            alpha = 9 / (4 * dx2);
            t=(xdpp-4*xdp)/3;
            alist[1] = alpha;
            blist[1] = 2*alpha*t;
            clist[1] = alpha*t*t;
        }
        else{ //X+ first order
            alpha = 1/dx2;
            alist[1] = alpha;
            blist[1] = -2*xdp*alpha;
            clist[1] = xdp*xdp*alpha;
        }
    }
    if( fixed[iym]==TRUE) {//Y-- second order
        if( fixed[iymm]==TRUE){ 
            alpha = 9 / (4 * dy2);
            t=(4*ydm-ydmm)/3;
            alist[2] = alpha;
            blist[2] = -2*alpha*t;
            clist[2] = alpha*t*t;
        }
        else{ //Y- first order
            alpha = 1/dy2;
            alist[2] = alpha;
            blist[2] = -2*ydm*alpha;
            clist[2] = ydm*ydm*alpha;
        }
    } 
    if(fixed[iyp]==TRUE) {//Y++ second order
        if( fixed[iypp]==TRUE){ 
            alpha = 9 / (4 * dy2);
            t=(ydpp-4*ydp)/3;
            alist[3] = alpha;
            blist[3] = 2*alpha*t;
            clist[3] = alpha*t*t;
        }
        else{//Y+ first order
            alpha = 1/dy2;
            alist[3] = alpha;
            blist[3] = -2*ydp*alpha;
            clist[3] = ydp*ydp*alpha;
        }
    }        
    
    if(fixed[izm]==TRUE) {//Z-- second order
        if( fixed[izmm]==TRUE){ 
            alpha = 9 / (4 * dz2);
            t=(4*zdm-zdmm)/3;
            alist[4] = alpha;
            blist[4] = -2*alpha*t;
            clist[4] = alpha*t*t;
        }
        else{ //Z- first order
            alpha = 1/dz2;
            alist[4] = alpha;
            blist[4] = -2*zdm*alpha;
            clist[4] = zdm*zdm*alpha;
        }
    } 
    if(fixed[izp]==TRUE) {//Z++ second order

        if( fixed[izpp]==TRUE){
            alpha = 9 / (4 * dz2);
            t=(zdpp-4*zdp)/3;
            alist[5] = alpha;
            blist[5] = 2*alpha*t;
            clist[5] = alpha*t*t;
        }
        else{ //Z+ first order
            alpha = 1/dz2;
            alist[5] = alpha;
            blist[5] = -2*zdp*alpha;
            clist[5] = zdp*zdp*alpha;
        }
    }    
    /***************************
    *Try to find a 3D solution*
    ***************************/ 
    a=alist[0]+alist[1]+alist[2]+alist[3]+alist[4]+alist[5]; 
    b=blist[0]+blist[1]+blist[2]+blist[3]+blist[4]+blist[5];
    c=-1+clist[0]+clist[1]+clist[2]+clist[3]+clist[4]+clist[5];
    disc=b*b - 4*a*c;

    float disc2=sqrt(disc);
    float check;
    float tolerance=0.01;
    sol1=(-b + disc2)/(2*a);
    sol2=(-b - disc2)/(2*a);
    sol=(sol1 > sol2) ?  sol1 : sol2;

    if(disc <0 || sol <= 0 || disc < 0|| isnan(sol) || isinf(sol) ){
        //printf("Trying 2D solution\n");
        /***************************
         *Try to find a 2D solution*
         ***************************/
        altsol[0]=altsol[1]=altsol[2]=0;
        alta[0]=alist[2]+alist[3]+alist[4]+alist[5];
        alta[1]=alist[0]+alist[1]+alist[4]+alist[5];
        alta[2]=alist[0]+alist[1]+alist[2]+alist[3];
        altb[0]=blist[2]+blist[3]+blist[4]+blist[5];
        altb[1]=blist[0]+blist[1]+blist[4]+blist[5];
        altb[2]=blist[0]+blist[1]+blist[2]+blist[3];
        altc[0]=-1+clist[2]+clist[3]+clist[4]+clist[5];
        altc[1]=-1+clist[0]+clist[1]+clist[4]+clist[5];
        altc[2]=-1+clist[0]+clist[1]+clist[2]+clist[3];
        sol=pseudo_inf;
        for(int k=0; k < 3; k++){
            altdisc[k]=altb[k]*altb[k] - 4*alta[k]*altc[k];
            if(altdisc[k] >= 0){
                altsols[k][0]=( - altb[k] + sqrt(altdisc[k]))/(2*alta[k]);
                altsols[k][1]=( - altb[k] - sqrt(altdisc[k]))/(2*alta[k]);
                altsol[k] = (altsols[k][0] > altsols[k][1]) ? altsols[k][0] : altsols[k][1];
                //printf("2D Solution[%d] %f\n",k, altsol[k]);
                if( isnan(altsol[k]) == 0){
                    if(sol==pseudo_inf || /*(altsol[k]>d &&*/ sol > altsol[k]){
                        sol=altsol[k];
                        a=alta[k];
                        b=altb[k];
                        c=altc[k];
                    }
                }
            } 
        }
        //printf("2D Solution %f\n", sol);
        if ( sol ==pseudo_inf || sol <=0 ){
            altsol[0]=altsol[1]=altsol[2]=0;
            /***************************
             *Try to find a 1D solution*
             ***************************/ 
            //printf("Trying 1D solution\n");
            alta[0]=alist[2]+alist[3];
            alta[1]=alist[0]+alist[1];
            alta[2]=alist[4]+alist[5];
            altb[0]=blist[2]+blist[3];
            altb[1]=blist[0]+blist[1];
            altb[2]=blist[4]+blist[5];
            altc[0]=-1+clist[2]+clist[3];
            altc[1]=-1+clist[0]+clist[1];
            altc[2]=-1+clist[4]+clist[5];
            sol=pseudo_inf;
            for(int k=0; k < 3; k++){
                altdisc[k]=altb[k]*altb[k] - 4*alta[k]*altc[k];
                //print("alt disc: %f\n", altdisc[k]);
                if(altdisc[k] >= 0){
                    altsols[k][0]=( - altb[k] + sqrt(altdisc[k]))/(2*alta[k]);
                    altsols[k][1]=( - altb[k] - sqrt(altdisc[k]))/(2*alta[k]);
                    altsol[k] = (altsols[k][0] > altsols[k][1]) ? altsols[k][0] : altsols[k][1];
                    //printf("1D Solution[%d] %f\n",k, altsol[k]);
                    float test_sol;
                    
                    if( isnan(altsol[k]) == 0){
                       // printf("1D: %f %f %f\n", sol, altsol[k], d);
                       if(sol==pseudo_inf || /*(altsol[k]>d &&*/ sol > altsol[k]/*)*/){
                            sol=altsol[k];
                            a=alta[k];
                            b=altb[k];
                            c=altc[k];
                        }
                    }
                } 
            }//printf("1D Solution %f\n", sol);
        }
    } //else printf("3D Solution %f\n", sol);

    check=fabs(a*sol*sol + b*sol + c);
     if( sol <= 0 ||  isnan(sol) || isinf(sol) || sol >= 13333311.0 ){
        printf("ixm: %d\n",ixm);
        printf("Fixed points\n");
        printf("X: %d\t%d\t%d\t%d\t%d\n", fixed[ixmm], fixed[ixm], fixed[index], fixed[ixp], fixed[ixpp] );
        printf("Y: %d\t%d\t%d\t%d\t%d\n", fixed[iymm], fixed[iym], fixed[index], fixed[iyp], fixed[iypp] );
        printf("Z: %d\t%d\t%d\t%d\t%d\n", fixed[izmm], fixed[izm], fixed[index], fixed[izp], fixed[izpp] );

        printf("Distances\n");
        printf("X: %7.2f\t%7.2f\t%7.2f\t%7.2f\t%7.2f\n", xdmm, xdm, d, xdp, xdpp );
        printf("Y: %7.2f\t%7.2f\t%7.2f\t%7.2f\t%7.2f\n", ydmm, ydm, d, ydp, ydpp );
        printf("Z: %7.2f\t%7.2f\t%7.2f\t%7.2f\t%7.2f\n", zdmm, zdm, d, zdp, zdpp );
        printf("a: ");
        for(int k=0; k<6; k++) printf("%f ", alist[k]); printf("\n");
        printf("b: ");
        for(int k=0; k<6; k++) printf("%f ", blist[k]); printf("\n");
        printf("c: ");
        for(int k=0; k<6; k++) printf("%f ", clist[k]); printf("\n");
        printf("a=%f, b=%f, c=%f\n", a, b, c);
        printf("Discriminant: %f\n", disc);
        printf("Check: %f\n", check );
        printf("%f --> %f %f --> %f\n", distances[index], sol1, sol2, sol);
        if(isnan(sol)) pexit("Solution is NaN!", "", 1);
        else if(isinf(sol)) pexit("Solution is inf!", "", 1);
        else if(sol<=0) pexit("Solution is 0.", "", 1);
        //else if(check > tolerance) pexit("Bad solution to Eikonal update step\n", "", 1);
        else if(distances[index] > sol) pexit("New distance is less than previous", "", 1);
        else pexit("Problem with solution of Eikonal equation", "", 1);
    }
    //printf("%f ---> %f\n", distances[index], sol);
    distances[index]=sol;
    if( considered->left != NULL) update(considered->left, fixed, distances);
    if( considered->right != NULL) update(considered->right, fixed, distances);
}


void swap_node(struct node* a, struct node* b){
    float* tdist=a->dist;
    int tindex=a->index;
    int z=a->z;
    int y=a->y;
    int x=a->x;
    a->dist=b->dist;
    a->index=b->index;
    a->z=b->z;
    a->y=b->y;
    a->x=b->x;

    b->dist=tdist;
    b->index=tindex;
    b->z=z;
    b->y=y;
    b->x=x;

}

struct node* min(struct node** in, int depth){
    struct node* smallest=NULL;
    struct node* child=NULL;

    if( (*in)->left == NULL){//First node is the smallest
        smallest=*in;
        child=(*in)->right;
        if( child != NULL){//Check if right node exists
           (*in)=child; //If the right node exists, link input node to right node
        }
        else{
            //Both left and right nodes are NULL, set in to NULL
            (*in)=NULL;
        }
    } 
    else if ((*in)->left != NULL && (*in)->left->left == NULL ){
        //Found second to last node
        //Left node is the smallest node
        smallest=(*in)->left;
        //the smallest node may have a right node (child), check if it exists
        child=smallest->right;
        if(child != NULL){
            //      15              15
            //     /  \            /  \
            //   10   17   -->   11   17
            //  /  \             /  
            //NULL  11        NULL  10 
            //
            //Move the right node up a level
            (*in)->left=child;
        }
        else{
            //If the smallest node does not have a child, 
            //sever the link between current node (in) and smallest
            (*in)->left=NULL;
        }
    }
    else if((*in)->left != NULL){
        smallest=min(&((*in)->left), depth+1);
    }
    if( *(smallest->dist) == pseudo_inf ){ 
        printf("Smallest equals pseudo_inf %d %d %d\n", smallest, *in, smallest->right);
        printf("current dist: %f, smallest: %f\n", *((*in)->dist), *smallest->dist ); exit(0);
    }
    return(smallest);
}


struct node* max(struct node* in){
    if(in->right != NULL) max(in->right);
    return(in);
}

int delete(struct node* in, float val, int index ){
    if( *(in->dist) == val && in->index == index){
        free(in);
        return(1);
    }
    else if( *(in->dist) <= val){
        if(delete(in->left, val, index)==1) in->left=NULL;

    }
    else if ( *(in->dist) > val){
        if(delete(in->right, val, index)==1) in->right=NULL;
    }
    return(0);
}
void search(struct node* in, int depth){
    if(in != NULL) {
        printf("%d: %d %d %d %d %f\n",depth, in->index, in->z, in->y, in->x, *(in->dist) );
        if( in->left != NULL) { printf("Left: "); search(in->left, depth+1);}
        if( in->right != NULL) { printf("Right: "); search(in->right, depth+1);}
    }
}

void free_nodes(struct node* in){
    if(in->right != NULL){ 
        free_nodes(in->right);
        free(in->right);
    }
    if(in->left != NULL){ 
        free_nodes(in->left);
        free(in->left);
    }
}

struct node* newNode(float* dist, int index, int z, int y, int x){
    struct node* new=malloc(sizeof(*new));
    new->dist=dist;
    new->index=index;
    new->z=z;
    new->y=y;
    new->x=x;
    new->left=NULL;
    new->right=NULL;
    return(new);

}

struct node* insert(struct node* in, float* dist, int index, int z, int y, int x){
    if(in==NULL){
        struct node* new=newNode(dist, index, z, y, x);
        in=new;
    }
    else if( *dist <= *(in->dist) ){
        in->left=insert(in->left, dist, index, z, y, x) ;
    } else{
        in->right=insert(in->right, dist, index, z,y,x);
    }
    return(in);
    
}

void wm_dist(data* img,int* img_vol, int** gm_border, int n, int label,  int start, int step){
    //Calculate inner and outer border voxels
    //Calculate distances to outer border voxels
    //Add border voxels to inner region
    int zmax=img->zmax;
    int ymax=img->ymax;
    int xmax=img->xmax;
    int xymax=xmax*ymax;
    float dx=img->xstep, dy=img->ystep, dz=img->zstep;
    int zi, yi, xi, z1, y1, x1, xp, xm, yp, ym, zp, zm;
    int i0, i1, i2, i3, i4, i5;
    int max =img->n3d;
    _Bool* fixed=calloc(max, sizeof(*fixed));
    _Bool* considered_array=calloc(max, sizeof(*considered_array));
    struct node* considered=NULL;
    float* distances=malloc(max*sizeof(*distances));
    int nfixed, nconsidered, nouter;
    int index;


    char testfn[100];
    float a=0,b=0,c=-1, alpha, t;
    int x, y, z;
    float disc, sol1, sol2, sol;
    float altsol[3];
    int wm_total=0;

    for(int i=0; i<max; i++){ 
        distances[i]=pseudo_inf; //Initialize starting values to infinity (pseudo_inf)
        fixed[i]=considered_array[i]=FALSE;
        if(img_vol[i]==label) wm_total++;
    }

    for(int i=start; i < n; i += step){ //Iterate over nodes on WM-GM border
        printf("i = %d\n", i); 
        index=gm_border[i][0];
        fixed[index]=TRUE; //Set the starting fixed node
        nfixed=1; 
        nconsidered=nouter=0; 
        z=gm_border[i][1];
        y=gm_border[i][2];
        x=gm_border[i][3];
        distances[index]=0;

        /************************************
         * Initialize the considered points *
         ************************************/
        xp=x+1;
        xm=x-1;
        yp=y+1;
        ym=y-1;
        zp=z+1;
        zm=z-1;
        i0=z*xymax+y*xmax+xp;
        i1=z*xymax+y*xmax+xm;
        i2=z*xymax+yp*xmax+x;
        i3=z*xymax+ym*xmax+x;
        i4=zp*xymax+y*xmax+x;
        i5=zm*xymax+y*xmax+x;

        if (img_vol[i0]==label && fixed[i0]==FALSE && considered_array[i0] == FALSE ){ 
                considered=insert(considered, &(distances[i0]), i0, z, y, xp);
                considered_array[i0] = TRUE;
                nconsidered++;
        }
        if (img_vol[i1]==label && fixed[i1]==FALSE  && considered_array[i1] == FALSE ){ 
                considered=insert(considered, &(distances[i1]), i1, z, y, xm);
                considered_array[i1] = TRUE;
                nconsidered++;
        }
        if (img_vol[i2]==label && fixed[i2]==FALSE && considered_array[i2] == FALSE  ){ 
                considered=insert(considered, &(distances[i2]), i2, z, yp, x);
                considered_array[i2] = TRUE;
                nconsidered++;
        }
        if (img_vol[i3]==label && fixed[i3]==FALSE && considered_array[i3] == FALSE  ){ 
                considered=insert(considered, &(distances[i3]), i3, z, ym, x);
                considered_array[i3] = TRUE;
                nconsidered++;
        }
        if (img_vol[i4]==label && fixed[i4]==FALSE && considered_array[i4] == FALSE  ){ 
                considered=insert(considered, &(distances[i4]), i4, zp, y, x);     
                considered_array[i4] = TRUE;
                nconsidered++;
        }
        if (img_vol[i5]==label && fixed[i5]==FALSE && considered_array[i5] == FALSE  ){
                considered=insert(considered, &(distances[i5]), i5, zm, y, x); 
                considered_array[i5] = TRUE;
                nconsidered++;
        }

        while(nconsidered > 0){
             //if(nfixed % 100 == 0) printf("%f %d %d\r", (float) 100.0 * nfixed/wm_total, nfixed, nconsidered);
            /*******************************
             * 1. Update considered points *    
             *******************************/
            update(considered, fixed, distances);
            /***************************************************
             * 2. Find smallest distance within border points  *
             ****************************************************/
            //search(considered, 1);
            struct node* minNode=min(&considered, 1);
            if( *(minNode->dist) == pseudo_inf) {
                printf("Dist is pseudo inf!\n"); 
                search(considered, 1);
                exit(0);
            }
            if( considered == NULL && nconsidered > 1){ 
                printf("nconsidered ==%d, considered=%d\n", nconsidered, considered);
                exit(0);   
            }
            /********************************************************
             * 3. Add mininum distance point to list of fixed points
             ********************************************************/
            index=minNode->index;
            z=minNode->z;
            y=minNode->y;
            x=minNode->x;
            //printf("Minimum: %d %d %d %d %f\n",index,z,y,x,*(minNode->dist) );
            //search(considered, 1);
            //printf("End second search\n");
            

            fixed[index]=TRUE;
            nfixed++;
            free(minNode);
            //Decrease the number of considered points
            nconsidered--;
    
            /***********************************************************************
             * Add points surrounding new fixed point to list of considered points
             ***********************************************************************/
            xp=x+1;
            xm=x-1;
            yp=y+1;
            ym=y-1;
            zp=z+1;
            zm=z-1;
            i0=z*xymax+y*xmax+xp;
            i1=z*xymax+y*xmax+xm;
            i2=z*xymax+yp*xmax+x;
            i3=z*xymax+ym*xmax+x;
            i4=zp*xymax+y*xmax+x;
            i5=zm*xymax+y*xmax+x;

            if (img_vol[i0]==label && fixed[i0]==FALSE && considered_array[i0] != TRUE ){ 
                    considered=insert(considered, &(distances[i0]), i0, z, y, xp);
                    considered_array[i0] = TRUE;
                    nconsidered++;
            }
            if (img_vol[i1]==label && fixed[i1]==FALSE && considered_array[i1] != TRUE ){ 
                    considered=insert(considered, &(distances[i1]), i1, z, y, xm);
                    considered_array[i1] = TRUE;
                    nconsidered++;
            }
            if (img_vol[i2]==label && fixed[i2]==FALSE  && considered_array[i2] != TRUE){ 
                    considered=insert(considered, &(distances[i2]), i2, z, yp, x);
                    considered_array[i2] = TRUE;
                    nconsidered++;
            }
            if (img_vol[i3]==label && fixed[i3]==FALSE  && considered_array[i3] != TRUE){ 

                    considered=insert(considered, &(distances[i3]), i3, z, ym, x);
                    considered_array[i3] = TRUE;
                    nconsidered++;
            }
            if (img_vol[i4]==label && fixed[i4]==FALSE && considered_array[i4] != TRUE ){ 
                    considered=insert(considered, &(distances[i4]), i4, zp, y, x);
                    considered_array[i4] = TRUE;
                    nconsidered++;
            }
            if (img_vol[i5]==label && fixed[i5]==FALSE && considered_array[i5] != TRUE ){
                    considered=insert(considered, &(distances[i5]), i5, zm, y, x); 
                    considered_array[i5] = TRUE;
                    nconsidered++;
            }

            //printf("%d\t%d\t%d\n",i,nfixed, nconsidered);
            //if( nfixed > 10) break;
        }
        if(i==8){ 
            for(int i=0; i<max; i++) distances[i] *= fixed[i];
            sprintf(testfn,"test_wm.mnc", i );
            writeVolume(testfn, distances, img->start, img->step, img->wcount, MI_TYPE_FLOAT );
            exit(0);
         }

        if(considered != NULL){
            free_nodes(considered);
            free(considered);
        }
        considered=NULL;

        for(int j=0; j<max; j++){ 
            distances[j]=pseudo_inf;
            fixed[j]=FALSE;
        }

    }


}

void* wm_dist_threaded(void* args){
    data* img= ((struct wm_vol_args*) args)->img;
    int label= ((struct wm_vol_args*) args)->label;
    int* img_vol=((struct wm_vol_args*) args)->img_vol;
    float* mat=((struct wm_vol_args*) args)->mat;
    int** gm_border=((struct wm_vol_args*) args)->gm_border;
    int n=((struct wm_vol_args*) args)->n;
    int thread= ((struct wm_vol_args*) args)->thread;
    int nthreads= ((struct wm_vol_args*) args)->nthreads;
    wm_dist(img, img_vol, gm_border, n, label, thread, nthreads);
}

int wm_dist_multithreaded(data* img, int* img_vol, int** gm_border, int label, float* mat, int n, int nthreads){
    int rc;
    pthread_t threads[nthreads];
    struct wm_vol_args thread_args[nthreads];
    for(int t=0; t< nthreads; t++){
        thread_args[t].img_vol=img_vol;
        thread_args[t].img=img;
        thread_args[t].label=label;
        thread_args[t].mat=mat;
        thread_args[t].n=n;
        thread_args[t].gm_border=gm_border;
        thread_args[t].thread=t;
        thread_args[t].nthreads=nthreads;
        rc= pthread_create(&threads[t], NULL, wm_dist_threaded, (void *) &thread_args[t] ); 
        if(rc!=0) pexit("Error: creating thread","" , 1);
    }
    for(int t=0; t<nthreads;t++){
        rc = pthread_join(threads[t], NULL);
        if(rc!=0) pexit("Error: joining thread", "", 1);
    }
    return(0);
}






/***************************************************************************************************
 *Author: Thomas Funck 
 *Contact: thomas.funck@mail.mcgill.ca
 *Date: January 15, 2016
 *Location: Montreal Neurological Institute, Montreal QC
 *
 *Name: surf_dist  
 *Synopsis: - program to calculate minimum distance from a region to any point within anatomically constrained space
 *
 *Inputs: 
 *          1) Mesh composed of vertices in .obj format
 *          2) Binary mask file (.txt) containing "1" for vertices inside region
 *
 *Outputs: 
 *          1) File (.txt) containing distance values at each node
 * 
 *
 *Description:
 *  Calculates distances from 3D region within anatomically constrained region. 
 *
 *
 *Version History:
 *
 *January, 2016
 *   Version-1.0: 
 *
 *
 * ****************************************************************/
int main(int argc, char** argv){
    if (argc < 4 || strcmp(argv[0],"-help") ==0  ) useage();
    
    int i=1;
    float dt=0.1;
    if(strcmp(argv[i++],"-dt") == 0){
        dt=atof(argv[i++]);
    }
    else i=1;
    data img;
    double* dist_vol;
    img.filename=argv[i++];
    int GM=atoi(argv[i++]);
    int WM=atoi(argv[i++]);
    int nvertices;
    int* img_vol;
    char *file_inputs[]={img.filename,  NULL}; //={mesh_filename, node_values_filename};
    int nthreads=1; //sysconf(_SC_NPROCESSORS_ONLN);
    float* mat=NULL;
    int n;
    VERBOSE=FALSE;

    //Check input files to make sure they exist
    if (check_input_files(file_inputs) != 0) exit(1);
    
    img_vol=(int*) readVolume(&img, 1, MI_TYPE_INT);
    dx=img.xstep;
    dx2=dx*dx;
    dy=img.ystep;
    dy2=dy*dy;
    dz=img.zstep;
    dz2=dz*dz;
    xmax=img.xmax;
    ymax=img.ymax;
    zmax=img.zmax;
    xymax=xmax*ymax;
    for(int z=0; z<img.zmax; z++ )
        for(int y=0; y<img.ymax; y++ )
            for(int x=0; x<img.xmax; x++ ) 
                if(z<=1 || y <= 1 || x<=1 || z >= img.zmax-2 || y >= img.ymax-2 || x >= img.xmax-2){
                    img_vol[z*img.ymax*img.xmax+y*img.xmax+x]=0;
                }

    /*Find GM-WM Border*/
    int** gm_border=wm_gm_border(&img, GM, WM, img_vol,  &n);
    //wm_dist_multithreaded(&img, img_vol, gm_border, WM, mat, n, nthreads);
    wm_dist(&img, img_vol, gm_border, n, WM, 0,1);
    return 0;
}

void useage(){
    printf("Name:\nwm_dist ~ calculate distances in wm.\n");
    printf("Description:\n");
    printf("Useage:\n<-dt increment> input.mnc GM WM \n");
    exit(1);
}
