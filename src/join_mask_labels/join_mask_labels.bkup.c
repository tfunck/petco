#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <omp.h>
#include "minc2.h"
#include "hdf5.h"
#include "minc_helper.h"

void useage();

float inverse3x3(double* m){
    //[ 0, 1, 2 ]
    //[ 3, 4, 5 ]
    //[ 6, 7, 8 ]
    //Get minors
    double d0 = m[4] * m[8] - m[7] * m[5];
    double d1 = m[3] * m[8] - m[6] * m[5];
    double d2 = m[3] * m[7] - m[6] * m[4];
    double d3 = m[1] * m[8] - m[7] * m[2];
    double d4 = m[0] * m[8] - m[6] * m[2];
    double d5 = m[0] * m[7] - m[6] * m[1];
    double d6 = m[1] * m[5] - m[4] * m[2];
    double d7 = m[0] * m[5] - m[3] * m[2];
    double d8 = m[0] * m[4] - m[3] * m[1];
    double det = m[0] * d0 - m[3] * d3 + m[6] * d6;
    //printf("%f %f\n", d0, m[0]);
    //printf("%f = %f - %f + %f\n", det, m[0] * d0, m[3] * d3, m[6] * d6);
    if(det != 0){
        //[ 0, -3, 6 ]
        //[ -1, 4, -7 ]
        //[ 2, -5, 8 ]
        m[0]= d0 / det;
        m[1]=-d3 / det;
        m[2]= d6 / det;
        m[3]=-d1 / det;
        m[4]= d4 / det;
        m[5]=-d7 / det;
        m[6]= d2 / det;
        m[7]=-d5 / det;
        m[8]= d8 / det;

    }
    return(det);
}
int min(int a, int b){
    if( a < b){ 
        return a;
    }
    else return b;
}

int max(int a, int b){
    if( a > b) return a;
    else return b;
}


int get_start(int x, int y){
    int diff=x-y;

    if( diff < 0 ) return(0);
    
    return(diff);

}

int fz(int i /*y*/, int j /*x*/, int k /*z*/, int* dim){ return(k*dim[1]*dim[0] + i*dim[0] + j);}
int fy(int i /*z*/, int j /*x*/, int k /*y*/, int* dim){ return(i*dim[2]*dim[1] + k*dim[1] + j);}
int fx(int i /*z*/, int j /*y*/, int k /*x*/, int* dim){ return(i*dim[1]*dim[2] + j*dim[2] + k);}

void write_z(float i /*y*/, float j /*x*/, float k /*z*/, FILE* fn){ fprintf(fn,"%f,%f,%f", k, i, j  ); }
void write_y(float i /*z*/, float j /*x*/, float k /*y*/, FILE* fn){ fprintf(fn,"%f,%f,%f", i, k, j  ); }
void write_x(float i /*z*/, float j /*y*/, float k /*x*/, FILE* fn){ fprintf(fn,"%f,%f,%f", i, j, k  ); }

int find_min(int* slab, int kernel_offset, int voxel_id,int* dim, int* cvox, int* k, float* step, float* start, int (*f)(), float*** knots, int vox, int* nknots){
    int *cur_i=&(cvox[0]), *cur_j=&(cvox[1]),*cur_k=&(cvox[2]);
    int i_start=get_start(*cur_i, kernel_offset);
    int j_start=get_start(*cur_j, kernel_offset);
    int i_end=i_start+kernel_offset+1;
    int j_end=j_start+kernel_offset+1;
    int dim1=dim[0], dim2=dim[1], dim3=dim[2];
    int index, index_min=0;
    float min=dim[0]*step[0], d;
    for(int i=i_start; i < i_end; i++){
        for(int j=j_start; j<j_end; j++){
            //printf("\t\t%d\n", j);
            index=(*f)(i,j, cvox[2],dim);
            int value=slab[index];


            //if(*nknots > 10) printf("\t%d %d = %d\n", i, j, value);
            if( value == voxel_id  ){
                //Calculates distance of test voxel from projection of current voxel onto current plane
                d=step[0] * step[0] * (*cur_i-i) * (*cur_i-i) + step[1] * step[1] * (*cur_j-j) * (*cur_j-j) ; 
                if( d < min){ 
                    min=d;
                }
                //if(*nknots > 10) printf("%d %d %d --> %f\n", i, j, *cur_k, d);
                
                if( d <= min ){
                    d = min;
                    if(i_start == 0) *cur_i=i;
                    else *cur_i= i ;
                    if(j_start == 0) *cur_j=j;
                    else *cur_j= j ;
                    *k=*cur_k;
                    index_min=index;
                }
            }
        }
    }
    if(index_min != 0){

        slab[index_min] *= -1;
        *nknots += 1;
        knots[vox]=realloc(knots[vox], sizeof(**knots) * *nknots );
        knots[vox][*nknots-1]=NULL;
        knots[vox][*nknots-1]=realloc(knots[vox][*nknots-1], sizeof(***knots) * 3 );
        knots[vox][*nknots-1][0]=*cur_i * step[0]+start[0];
        knots[vox][*nknots-1][1]=*cur_j * step[1]+start[1];
        knots[vox][*nknots-1][2]=*cur_k * step[2]+start[2];
        //printf("%d. New knot (world): %3.3f %3.3f %3.3f\n", *nknots, knots[vox][*nknots-1][0], knots[vox][*nknots-1][1], knots[vox][*nknots-1][2]);
        //printf("%d. New knot (voxel): %d %d %d\n", *nknots, *cur_i, *cur_j, *cur_k);
    }
    return(0);
}

int search_slab(int* slab,int kernel_offset, int* cvox, int* k, int* dim, int voxel_id, float* step, float* start, int (*f)(int, int, int, int* ), float*** knots, int vox, int* nknots, int search_depth ){
    int max_slice=*k+search_depth;
    //printf("%d %d %d\n", cvox[0], cvox[1], cvox[2]);
    if(cvox[2]==dim[2] || cvox[2] >  max_slice  ) return(0);

    find_min(slab, kernel_offset, voxel_id, dim, cvox, k, step, start, f, knots, vox, nknots);
    cvox[2] += 1;
    search_slab(slab,  kernel_offset,  cvox, k, dim, voxel_id, step, start, f, knots,vox, nknots, search_depth);
    
    return(0);
}

int set_parameters(int* dim, int amax, int bmax, int cmax, float*start, float astart, float bstart, float cstart, float* step, float astep, float bstep, float cstep){
    dim[0]=amax;
    dim[1]=bmax;
    dim[2]=cmax;

    start[0]=astart;
    start[1]=bstart;
    start[2]=cstart;

    step[0]=astep;
    step[1]=bstep;
    step[2]=cstep;
    return(0);
}



int get_lin_domain(double i0, double j0, double k0, double i1, double j1, double k1, double t, double* cur_i, double* cur_j, double* cur_k){
    //printf("y0: %f, y1: %f, t: %f\n", y0, y1, t);
    *cur_i = i0 + (i1 - i0) * t;
    *cur_j = j0 + (j1 - j0) * t;
    *cur_k = k0 + (k1 - k0) * t;
   return(0); 
}

int fill_voxels(float*** knots, int* nknots, int k, int n, int* dim, float* step, float* start, int* voxel_ids,float* distances,int* mask_vol, int (*findex)(int, int, int, int* ) ){
    double delta = sqrt( step[0]*step[0] + step[1]*step[1] + step[2]*step[2]  )/9;
    double x, y, d, z, z2;
    double sign=(double) fabs(start[0])/ start[0];
    double a;
    double b;
    double c;
    int nset=0;

    
    for(int e=0; e < n; e++){
        //printf("e=%d\n", e);
        for(int f=0; f < nknots[e]-1; f++){
            //if( f < nknots[e] -1){
            //    a=eqs[e][f][0];
            //    b=eqs[e][f][1];
            //    c=eqs[e][f][2];
            //}
            //printf("[a, b, c] = %f %f %f\n", a, b, c);
            double iw=(double) knots[e][f][0];
            double jw=(double) knots[e][f][1];
            double kw=(double) knots[e][f][2];
            double iwp=(double)knots[e][f+1][0];
            double jwp=(double)knots[e][f+1][1];
            double kwp=(double)knots[e][f+1][2];
            double ciw=iw, cjw=jw, ckw=kw;

            int i=real2voxel( (float) iw, start[0], step[0]);
            int j=real2voxel( (float) jw, start[1], step[1]);
            int k=real2voxel( (float) kw, start[2], step[2]);
            int ip=real2voxel( (float) iwp, start[0], step[0] );
            int jp=real2voxel( (float) jwp, start[1], step[1] );
            int kp=real2voxel( (float) kwp, start[2], step[2] );
            int ci=i;//z
            int cj=j;//x
            int ck=k;//y
            double t = 0; //hyp( cjw, ckw); //double t2= hyp( jwp, kwp);
            double max_dist=sqrt((iw-iwp)*(iw-iwp)+(jw-jwp)*(jw-jwp)+(kw-kwp)*(kw-kwp));
            //printf("\tFrom: %3.3f, %3.3f, %3.3f\n", iw, jw, kw );
            //printf("\tTo: %3.3f, %3.3f, %3.3f\n", iwp, jwp, kwp );
            //printf("\tdelta: %f\n", delta);
            int counter=0;
            double nsteps=0.2;
            while( t < 1 ){
                t += nsteps;

                get_lin_domain(iw, jw, kw, iwp, jwp, kwp, t, &ciw, &cjw, &ckw); 

                ci=real2voxel((float) ciw, start[0], step[0]);
                cj=real2voxel((float) cjw, start[1], step[1]);
                ck=real2voxel((float) ckw, start[2], step[2]);
                double i_vox= ci*step[0]+start[0];
                double j_vox= cj*step[1]+start[1];
                double k_vox= ck*step[2]+start[2];

                d= ( i_vox - ciw  )*( i_vox - ciw  );
                d+=( j_vox-(double) cjw )*( j_vox-(double) cjw );
                d+=( k_vox -(double) ckw )*( k_vox -(double) ckw );

                int index=(*findex)(ci,cj,ck,dim);
                //printf("%f %f\n", d ,distances[index]);
                if( d < distances[index] ){ 
                    distances[index]=d;
                    mask_vol[index]=30; //voxel_ids[e]; //voxel_ids[e];
                    //printf("Set! %d\n", voxel_ids[e]);
                    nset++;
                }

                counter++;
                if( counter > 50) break;

            } 
        } 
    }
    printf("\t%d", nset);

    return(0);
}

int get_knots(int* dim, int* mask_vol, int k, float**** knots, int** nknots, int** voxel_ids, int* n, int search_depth, float* step, float* start, int search_radius, int (*f)(int, int, int, int* ) ){
    //int max_knots=100; 
    int kernel_offset=round((search_radius + 1) / 2);
    int m;
    int cvox[3];
    int* tmp_nknots=*nknots;
    int* tmp_voxel_ids=*voxel_ids;
    float*** tmp_knots=*knots;

    for(int i=0; i < dim[0]; i++){
        for(int j=0; j < dim[1]; j++){
            int index = (*f)(i,j,k,dim);
            int voxel_id=mask_vol[index];

            //if( i==326 && j==1164){
            if( voxel_id  > 0.1  ){
                //mask_vol[index] *= -1;

                *n += 1;
                
                if( *n > 1){
                    tmp_nknots=realloc( tmp_nknots, sizeof(*tmp_nknots) * *n );
                    tmp_voxel_ids=realloc( tmp_voxel_ids, sizeof(*tmp_voxel_ids) * *n );
                    tmp_knots=realloc(tmp_knots, sizeof(*tmp_knots) * *n  );
                }
                tmp_voxel_ids[*n-1]=voxel_id;
                
                m=1;

                tmp_knots[ *n -1]=malloc(sizeof(**tmp_knots) * *n );
                tmp_knots[ *n -1][0]=malloc(sizeof(***tmp_knots) * 3);
                tmp_knots[ *n -1][0][0]=i*step[0]+start[0];
                tmp_knots[ *n -1][0][1]=j*step[1]+start[1];
                tmp_knots[ *n -1][0][2]=k*step[2]+start[2];
                //printf("1. Voxel %d %d %d\n",  i, j, k); 
                //printf("1. World %f %f %f\n", tmp_knots[ *n -1][0][0], tmp_knots[ *n -1][0][1],tmp_knots[ *n -1][0][2] );
                cvox[0]=i;
                cvox[1]=j;
                cvox[2]=k+1;
                search_slab(mask_vol, kernel_offset, cvox, &k, dim, voxel_id, step, start, f, tmp_knots, *n-1, &m, search_depth );
                tmp_nknots[ *n - 1]=m;

            }
            //}
        }
       //if( *n == max_knots) break;
    }
    *nknots=tmp_nknots;
    *voxel_ids=tmp_voxel_ids;
    *knots=tmp_knots;
    return(0);
}


int set_dim_parameters(data* mask, char* search_dim, int* dim, float* start, float* step  ){
    //Figure out dimension along which we will use the splines
    if( strcmp(search_dim, "z") == 0){
        set_parameters(dim, mask->ymax, mask->xmax, mask->zmax, start, mask->ystart, mask->xstart, mask->zstart, step, mask->ystep, mask->xstep, mask->zstep);
    }
    else if ( strcmp(search_dim, "y") == 0){
        set_parameters(dim, mask->zmax, mask->xmax, mask->ymax, start, mask->zstart, mask->xstart, mask->ystart, step, mask->zstep, mask->xstep, mask->ystep);
    }
    else if ( strcmp(search_dim, "x") == 0){
        set_parameters(dim, mask->zmax, mask->ymax, mask->xmax, start, mask->zstart, mask->ystart, mask->xstart, step, mask->zstep, mask->ystep, mask->xstep);
    }
    else { 
        printf("Error could not recognize dimension: %s. Please run with <z/y/x>.\n", search_dim);
        exit(0);
    }

    return(0);
}

int main(int argc, char** argv){
    if(argc != 7 || strcmp(argv[1], "-help") == 0 ){
        useage();
    }

    data mask;
    int* mask_vol;
    mask.filename=argv[1];
    char* search_dim=argv[2];
    int k_start=atoi(argv[3]);
    int search_radius=atoi(argv[4]);
    int search_depth=atoi(argv[5]);
    char* output_fn=argv[6];
    int kernel_offset;
    int i,j;
    int n;
    int max;
    int dim[3];
    float start[3];
    float step[3];
    float* distances;
    int* voxel_ids;
    float*** knots;
    int* nknots;
    int (*f)(int,int,int,int*);
    void (*write)(float, float, float, FILE*);

    if( strcmp(search_dim, "z") == 0) f=&fz;
    else if ( strcmp(search_dim, "y") == 0) f=&fy;
    else if ( strcmp(search_dim, "x") == 0) f=&fx;
    else useage();

    mask_vol=(int*)readVolume(&mask, 1, MI_TYPE_INT);
    distances=malloc( mask.n3d * sizeof(*distances));
    for(int i=0; i < mask.n3d; i++){
        distances[i]=mask.n3d * (step[0]+step[1]+step[2])/3;//set a big number for default distance
    }
    set_dim_parameters(&mask, search_dim, dim, start, step);
    for(int k=k_start; k<dim[2]; k++){
        voxel_ids=malloc(sizeof(*voxel_ids));
        knots=malloc(sizeof(*knots));
        nknots=malloc(sizeof(*nknots));
        n=0;
        get_knots(dim, mask_vol, k, &knots,  &nknots, &voxel_ids, &n, search_depth, step, start, search_radius, f );

        printf("k=%d : %d", k, n);
        //interpolate(dim, mask_vol, k, knots, eqs, nknots, n, step, start);
        fill_voxels(knots, nknots, k, n, dim, step, start,voxel_ids, distances, mask_vol, f );
        printf("\n");
        for(int i=0; i < n; i++){
            for(int k=0; k < nknots[i]; k++){ 
                free(knots[i][k]); 
            }
            free(knots[i]);
        }
        free(knots);
        free(voxel_ids);
        free(nknots);
        if( k > k_start + 10 ) break;
    }

    for(int i=0; i < mask.n3d; i++) if(mask_vol[i] < 0) mask_vol[i] *= -1;
    writeVolume(output_fn, mask_vol, mask.start, mask.step, mask.wcount, MI_TYPE_INT  );
    printf("Done.\n"); 

    return(0);
}


void useage(){
    printf("Useage: mask.mnc <axis: z, y, x> <axis start value> <voxel search radius> output.mnc  \n");
    printf("Purpose: Join mask labels across slices. The axis start value refers to the \n");
    printf("Example: join_mask_labels mymask.mnc z 0 10 mymask_interpol.mnc\n");
    exit(1);
}

/*int interpolate(int* dim, int* mask_vol,int k,float*** knots, float*** eqs, int* nknots, int n, float* step, float* start   ){
    double mat[3][3];
    double* mat_ptr=&mat;
    double a, b, c;
    double pi;
    double pj;
    double pk;
                    
    double qi;
    double qj;
    double qk;

    double ri;
    double rj;
    double rk;

    double x0;
    double x1;
    double x2;
    for(int d=0; d<n; d++){
        for(int e=0; e<n; e++){
            for(int f=0; f<3; f++) printf("%f, ", knots[d][e][f]);
            printf("\n");
        } printf("\n");
    } 
    for(int d=0; d<n; d++){
        if( nknots[d] >= 3){
            for(int e=0; e < nknots[d]-1; e++){
                mat[0][0]=mat[1][0]=mat[2][0]=1;
                //y = a + b x + c x^2

                if( e < nknots[d]-2 ){
                    pi=(double) knots[d][e][0];
                    pj=(double) knots[d][e][1];
                    pk=(double) knots[d][e][2];
                    
                    qi=(double) knots[d][e+1][0];
                    qj=(double) knots[d][e+1][1];
                    qk=(double) knots[d][e+1][2];

                    ri=(double) knots[d][e+2][0];
                    rj=(double) knots[d][e+2][1];
                    rk=(double) knots[d][e+2][2];

                    x0=(double) hyp(pj,pk);
                    x1=(double) hyp(qj,qk);
                    x2=(double) hyp(rj,rk);
                   
                    mat[0][1]=x0;
                    mat[1][1]=x1;
                    mat[2][1]=x2;

                    mat[0][2]=x0*x0;
                    mat[1][2]=x1*x1;
                    mat[2][2]=x2*x2;
                    printf("%f %f %f\n", pi, pj, pk);
                    printf("%f %f %f\n", qi, qj, qk);
                    printf("%f %f %f\n", ri, rj, rk);
                    printf("[ %3.1f\t\t%3.1f\t\t%3.1f ]\t[a] = [ %3.3f ]\n",mat[0][0],mat[0][1],mat[0][2], pi);
                    printf("[ %3.1f\t\t%3.1f\t\t%3.1f ]\t[b] = [ %3.3f ]\n",mat[1][0],mat[1][1],mat[1][2], qi);
                    printf("[ %3.1f\t\t%3.1f\t\t%3.1f ]\t[c] = [ %3.3f ]\n",mat[2][0],mat[2][1],mat[2][2], ri);
                    
                    inverse3x3(mat_ptr);
                    a=mat[0][0] * pi + mat[0][1] * qi + mat[0][2] * ri;
                    b=mat[1][0] * pi + mat[1][1] * qi + mat[1][2] * ri;
                    c=mat[2][0] * pi + mat[2][1] * qi + mat[2][2] * ri;
                    
                    
                    //printf("%3f = %3.3f + %3.3f (%3.3f) + %3.3f (%3.3f)^2 = %f\n", pi, a, b, x0, c, x0, a + b*x0 + c*x0*x0);
                    //printf("%3f = %3.3f + %3.3f (%3.3f) + %3.3f (%3.3f)^2 = %f\n", qi, a, b, x1, c, x1, a + b*x1 + c*x1*x1);
                    //printf("%3f = %3.3f + %3.3f (%3.3f) + %3.3f (%3.3f)^2 = %f\n", ri, a, b, x2, c, x2, a + b*x2 + c*x2*x2);
            
                }
                eqs[d][e][0]=a;
                eqs[d][e][1]=b;
                eqs[d][e][2]=c;
                //printf("%d. [a, b, c] = %f %f %f\n", d, a, b, c); 
                /*double x0=hyp(pj,pk);
                double x1=hyp(qj,qk);
                double x2=hyp(rj,rk);
                //double y=a+b*x+c*x*x;
                printf("%d %d\t%4.3f = %4.3f + %f * %4.3f + %4.3f * %4.3f = %4.3f\n",d,e,pi, a,b,x0,c, x0*x0, (a+b*x0+c*x0*x0) );
                printf("%d %d\t%4.3f = %4.3f + %f * %4.3f + %4.3f * %4.3f = %4.3f\n",d,e,qi, a,b,x1,c, x2*x1, (a+b*x1+c*x1*x1) );
                printf("%d %d\t%4.3f = %4.3f + %f * %4.3f + %4.3f * %4.3f = %4.3f\n",d,e,ri, a,b,x2,c, x2*x2, (a+b*x2+c*x2*x2) );
                //printf("Equations (%d, %d) : [a, b, c] = %f %f %f\n",d, e, a, b, c);
             }
        }
    } 

    return(0);
}
*/

/*
 int get_quad_domain(double a, double b, double c, double* x_ptr, double* y_ptr, double t, double z, double threshold, int max_iter){
    double y = *y_ptr ;
    double x = *x_ptr;
    double z1;
    double z2;
    double sign=fabs(*x_ptr) / *x_ptr ;
    //printf("x %f\n", *x_ptr);
    printf("\tError\t\tZ2\t\tGuess\t\tt\n");
    for(int i=0; i < max_iter; i++){
        z1 = quad(a,b,c,y,t,z);
        z2 = quad_diff(a,b,c,y,t,z);
        print("%d.\t%f\t%f\t%f\t%f\n",i, fabs(z1), z2, y, t);
        if ( z2 == 0) break;
        y = y - z1/z2; 
        if (fabs(z1) < threshold) break;
    }   
    
    x = sign * sqrt(t*t - y*y);
    //printf("x=%f\n", x);
    *x_ptr=x;
    *y_ptr=y;

    return(0);
}

double hyp(double x,double y){
    return(sqrt(x*x + y*y));
}



double quad(double a, double b, double c, double y, double t, double z  ){
    return( a + b * sqrt( sqrt(t*t - y*y) + y*y ) + c * (sqrt(t*t-y*y) + y*y) - z); 
}

double quad_diff(double a, double b, double c, double y, double t, double z  ){
    return(0.5*b*y*(2 - 1/sqrt(t*t - y*y ))/ sqrt(sqrt(t*t-y*y)+y*y) + c*y*( 2 - 1/sqrt(t*t-y*y)  )); 
}
   */
