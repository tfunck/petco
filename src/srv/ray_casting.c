#include <stdio.h>
#include <stdlib.h>
#include <malloc.h> 
#include <math.h>
#include "minc2.h"
#include "meshPVC.h"

//int local_VERBOSE = 1;
int local_VERBOSE = 0;



plane* calcLines2(vector* point1, vector* point2, vector* point3, vector* point4 )
{
 
 plane* newPlane=malloc(sizeof(plane));
//printf("%f\t%f\t%f\n%f\t%f\t%f\n%f\t%f\t%f\n\n", point1.x, point2.x, point3.x, point1.y, point2.y, point3.y, point1.z, point2.z, point3.z  );  
    newPlane->line1.x=point2->x - point1->x ;
    newPlane->line1.y=point2->y - point1->y ;
    newPlane->line1.z=point2->z - point1->z ;

    newPlane->line2.x=point4->x - point1->x ;
    newPlane->line2.y=point4->y - point1->y ;
    newPlane->line2.z=point4->z - point1->z ;

    newPlane->line3.x=point2->x - point3->x ;
    newPlane->line3.y=point2->y - point3->y ;
    newPlane->line3.z=point2->z - point3->z ;
   
    newPlane->line4.x=point4->x - point3->x ;
    newPlane->line4.y=point4->y - point3->y ;
    newPlane->line4.z=point4->z - point3->z ;
     
    newPlane->mag1= pow(newPlane->line1.x, 2) + pow(newPlane->line1.y, 2)  + pow(newPlane->line1.z, 2);
    newPlane->mag2= pow(newPlane->line2.x, 2) + pow(newPlane->line2.y, 2)  + pow(newPlane->line2.z, 2);
    newPlane->mag3= pow(newPlane->line3.x, 2) + pow(newPlane->line3.y, 2)  + pow(newPlane->line3.z, 2);

    newPlane->offset.x=point1->x;  newPlane->offset.y=point1->y;   newPlane->offset.z=point1->z;
    newPlane->offset2.x=point2->x;  newPlane->offset2.y=point2->y; newPlane->offset2.z=point2->z;
    newPlane->offset3.x=point3->x;  newPlane->offset3.y=point3->y; newPlane->offset3.z=point3->z; 
    newPlane->offset4.x=point4->x;  newPlane->offset4.y=point4->y; newPlane->offset4.z=point4->z; 
if(local_VERBOSE == 1){ 
    printf("(%1.3f)\t(%1.3f)\t(%1.3f)\t(%1.3f)\n", point1->x, point2->x, point3->x, point4->x);
    printf("(%1.3f)\t(%1.3f)\t(%1.3f)\t(%1.3f)\n", point1->y, point2->y, point3->y, point4->y);
    printf("(%1.3f)\t(%1.3f)\t(%1.3f)\t(%1.3f)\n", point1->z, point2->z, point3->z, point4->z);
    printf("line1\t\tline2\t\tline3\t\tline4\n");
printf("(%1.3f) \t(%1.3f) \t(%1.3f) \t(%1.3f)\n", newPlane->line1.x,newPlane->line2.x,newPlane->line3.x, newPlane->line4.x);
printf("(%1.3f) \t(%1.3f) \t(%1.3f) \t(%1.3f)\n", newPlane->line1.y,newPlane->line2.y,newPlane->line3.y,newPlane->line4.y);
printf("(%1.3f) \t(%1.3f) \t(%1.3f) \t(%1.3f)\n\n",newPlane->line1.z,newPlane->line2.z,newPlane->line3.z,newPlane->line4.z);
   }
    return newPlane;
}


plane* calcLines(vector* point1, vector* point2, vector* point3 )
{

 plane* newPlane=malloc(sizeof(plane));
//printf("%f\t%f\t%f\n%f\t%f\t%f\n%f\t%f\t%f\n\n", point1.x, point2.x, point3.x, point1.y, point2.y, point3.y, point1.z, point2.z, point3.z  );  
    newPlane->line1.x=point2->x - point1->x;
    newPlane->line1.y=point2->y - point1->y;
    newPlane->line1.z=point2->z - point1->z;

    newPlane->line2.x=point3->x - point1->x;
    newPlane->line2.y=point3->y - point1->y;
    newPlane->line2.z=point3->z - point1->z;
   
    newPlane->line3.x=point3->x - point2->x;
    newPlane->line3.y=point3->y - point2->y;
    newPlane->line3.z=point3->z - point2->z;

    newPlane->mag1= pow(newPlane->line1.x, 2) + pow(newPlane->line1.y, 2)  + pow(newPlane->line1.z, 2);
    newPlane->mag2= pow(newPlane->line2.x, 2) + pow(newPlane->line2.y, 2)  + pow(newPlane->line2.z, 2);
    newPlane->mag3= pow(newPlane->line3.x, 2) + pow(newPlane->line3.y, 2)  + pow(newPlane->line3.z, 2);

    newPlane->offset.x=point1->x;
    newPlane->offset.y=point1->y;
    newPlane->offset.z=point1->z;
  if(local_VERBOSE == 1){ 
    printf("(%1.3f)  (%1.3f) (%1.3f)\n", point1->x, point2->x, point3->x );
    printf("(%1.3f)  (%1.3f) (%1.3f)\n", point1->y, point2->y, point3->y );
    printf("(%1.3f)  (%1.3f) (%1.3f)\n", point1->z, point2->z, point3->z );
    printf("(%1.3f)+  (%1.3f)     (%1.3f)\n", point1->x, newPlane->line1.x, newPlane->line2.x);
    printf("(%1.3f)+a*(%1.3f) + b*(%1.3f)\n", point1->y, newPlane->line1.y, newPlane->line2.y);
    printf("(%1.3f)+  (%1.3f)     (%1.3f)\n\n", point1->z, newPlane->line1.z, newPlane->line2.z);
   }
    return newPlane;
}


plane** findPlanes(vector **inner_mesh, vector** outer_mesh,  int n_meshPts )
{
 int i, j=0, k=0;
 plane **local=malloc(sizeof(plane*) * n_meshPts * 3 );
 //REMEMBER: the first vertex in the inner and outer meshes is in the center of the surrounding neighbours
 //array of pointers to 3d vetors defining normals for triangles around central vertex
//printf("n_meshPts = %d\n", n_meshPts);
    for(i =1; i < n_meshPts; i++) //skip first inner mesh vertex by starting at 1
    {
        //find the vectors that define the planes that make up the inner surface
        if( i + 1 < n_meshPts)  local[j++]=calcLines( inner_mesh[0], inner_mesh[i], inner_mesh[i+1] );  
        else local[j++]=calcLines( inner_mesh[0], inner_mesh[i], inner_mesh[1] ); 
    }   

    for(i =1; i < n_meshPts; i++){ //skip first outer vertex by starting at 1
        //find the vectors that define the planes that make up the outer surface
        if( i + 1 < n_meshPts) local[j++]=calcLines( outer_mesh[0], outer_mesh[i], outer_mesh[i+1] ); 
        else local[j++]=calcLines( outer_mesh[0], outer_mesh[i], outer_mesh[1] ); 
    }   

    for(i =1; i < n_meshPts; i++) //skip first inner vertex by starting at 1
    {
         //find the vectors that define the planes that make up the walls between inner and outer surfaces 
        if( i + 1 < n_meshPts) local[j++]=calcLines2( inner_mesh[i], inner_mesh[i+1] ,  outer_mesh[i+1],outer_mesh[i] ); 
        else local[j++]=calcLines2(  inner_mesh[i], inner_mesh[1], outer_mesh[1] , outer_mesh[i]);
 	// if( i + 1 < n_meshPts) local[j++]=calcLines( inner_mesh[i], inner_mesh[i+1] ,  outer_mesh[i] ); 
        //else local[j++]=calcLines(  inner_mesh[i], inner_mesh[1] , outer_mesh[i]); 
    }   

//printf("Total number of planes: %d.\nShould equal %d.\n", j , (n_meshPts-1) * 3 );
    
    return local;
}

double dot( double A[3], double B[3])
{
return A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
}


int intersection( vector *x1,  vector* u, vector* v, vector* x0, vector** inner_mesh, vector** outer_mesh,  int meshPoint, int n_meshPoints, double rotate, double** matrix, int rect2 ) 
{
 double a, b, t;
 int i, k;
 vector x3;
 double au[3];
 double bv[3];
 x3.x=x1->x - x0->x;
 x3.y=x1->y - x0->y;
 x3.z=x1->z - x0->z;
 
//Define a vector r that goes to some arbitrary vector outside the mesh.
//by defuault rotate = 1, but if this is inconclusive then shift the ray by giving
//a positive value for r other than 1

matrix[0][0]=u->x; matrix[0][1]=v->x; matrix[0][2]=-x1->x * rotate; 	        matrix[0][3]=x3.x;
matrix[1][0]=u->y; matrix[1][1]=v->y; matrix[1][2]=-x1->y * pow(rotate, 2); 	matrix[1][3]=x3.y;
matrix[2][0]=u->z; matrix[2][1]=v->z; matrix[2][2]=-x1->z * pow(rotate, 3); 	matrix[2][3]=x3.z;
 rref( matrix, 4, 3);
/*
for( k =0; k < 3; k++) {
        for( i = 0; i < 4; i++ ) {printf("%1.3f\t", matrix[k][i]);
       
        } printf("\n");
    }
     printf("\n"); 

/**/

for(k=0; k<3; k++){
    for(i=0; i< 3; i++){
        if( matrix[k][i]== 1.0) {
            switch (i){
                case 0: { a=matrix[k][3]; break;}
                case 1: { b=matrix[k][3]; break;}
                case 2: { t=matrix[k][3]; break;}
                default: {printf("Error identifying plane variables.\n"); /*exit(0);*/ };       
            }
            break;
        }
        /* if we get to i=3 without finding any of the variables, that means all the values so far are 0.
         * so the linear combination of a bunch of zeros must also be zero (0*a+0*b+0*t = 0). if not, then
         * we have an inconsistent system and therefore there cannot be an intersection*/
        if(i == 2 && matrix[k][3] != 0) { 
		if (local_VERBOSE) printf("Inconsistent matrix.\n"); 
		return 0;
	}
    }
}

double s, q, r;  

if(local_VERBOSE) printf("Mesh Point %d:\ta= %f\tb=%f\tt=%f\n", meshPoint + 1, a, b, t);

if(t >= 0 && a >= 0 && b >= 0 && b <= 1  && a <= 1 ) 
	if( (rect2 == 0 && a != 1 && b != 1) || (rect2 == 1) )	{
    //this if-statement is to check if a and b are within the correct bounds for triagular plane segments.
    //to prevent overlap of points on u or v, dont allow b to == 0
  	if( a == 0 && b == 0 && t > 0 ) { 
		if(local_VERBOSE) printf("Hit central vertex. Only counting it once.\n");
	return -3;
	}
	if(  b != 0 ){
             au[0]= a*u->x; bv[0]=b*v->x; 
             au[1]= a*u->y; bv[1]=b*v->y;
             au[2]= a*u->z; bv[2]=b*v->z;
		matrix[0][0]= au[0]+bv[0]; matrix[0][1]=-v->x + u->x; matrix[0][2]=u->x; 
		matrix[1][0]= au[1]+bv[1]; matrix[1][1]=-v->y + u->y; matrix[1][2]=u->y; 
		matrix[2][0]= au[2]+bv[2]; matrix[2][1]=-v->z + u->z; matrix[2][2]=u->z; 

		if(local_VERBOSE){ 
		printf("%1.3f %1.3f %1.3f\n", matrix[0][0], matrix[0][1], matrix[0][2] );
		printf("%1.3f %1.3f %1.3f\n", matrix[1][0], matrix[1][1], matrix[1][2] );
		printf("%1.3f %1.3f %1.3f\n", matrix[2][0], matrix[2][1], matrix[2][2] );
		}	
		//printf("Target point = %1.10f %1.10f %1.10f\n", x0->x + matrix[0][0] -x1->x*t , x0->y+ matrix[1][0]-x1->y*t , x0->z + matrix[2][0]-x1->z*t   );
		//printf("Target point = %f %f %f\n", x1->x *(t +1)  , x1->y * (t+1), x1->z *( t + 1));
		rref( matrix, 3, 3);
		if(local_VERBOSE){ 
		printf("%1.2f %1.2f %1.12f\n", matrix[0][0], matrix[0][1], matrix[0][2] );
		printf("%1.2f %1.2f %1.12f\n", matrix[1][0], matrix[1][1], matrix[1][2] );
		printf("%1.2f %1.2f %1.12f\n", matrix[2][0], matrix[2][1], matrix[2][2] );
		}
		
		q=r=0;
		for(k=0; k<3; k++){
    			for(i=0; i< 3; i++){
        			if( matrix[k][i]== 1.0){ 
        				switch (i){
                			  case 0: { q=matrix[k][2]; break; }
                			  case 1: { r=matrix[k][2]; break; }
                			  case 2: break;
					  default: {printf("Error identifying plane variables.\n"); /*exit(0);*/ };       
           				}	
					//Below are two sets of termination conditions. The first make sure that that
					//the lines have the correct size. The second checks for an inconsistent
					//matrix which would mean that there is no intersection point and thus that 
					//the point falls inside the triangular segment.
					//if the magnitude of s*n is greater than the magnitude of n, then the point where
					//the intersection takes place is further than the point defined by p + n. and
					if( i==2 && matrix[k][2] != 0){
					if(local_VERBOSE) printf("Inconsistent Matrix for triangular segment\n");
					return 0;
					}
					if( q > 1)  { 
						if(t == 0 ) {
							if(local_VERBOSE) printf("Border point\n");
							return -1;
						 }
						if(local_VERBOSE) printf("Triangular Segment 1.\n"); 
					return 1;
					}
					else if(i == 2 && matrix[k][2] != 0) {
						if(local_VERBOSE) printf("Triangular Segment 2.\n"); 
					return 1;
					}
    					break;	
				}
			}
		}
		return 0; //the point is not in triangular segment
	}
}

return 0;
}

int pointinMesh(vector* point, vector** inner_mesh, vector** outer_mesh, int n_meshPts){
  plane **local;
  double rotate =1;
  double result=-2;
  //FILE* intersectionsFile=fopen("intersections.txt", "wt+");

	local=findPlanes(inner_mesh, outer_mesh,  n_meshPts);
	while(result == -2){
		result=traceRay( point, inner_mesh, outer_mesh, n_meshPts, local, rotate);
		rotate *= -1.1;
	if( rotate > 2) exit(0);
	};
return result;
}


int traceRay(vector* point, vector** inner_mesh, vector** outer_mesh, int n_meshPts, plane** local, double rotate )
{
 int i;
 int intersections = 0;
 vector r;
 int result;
 double extra=0;
//FILE* uvectorfile=fopen("uvector.txt", "wt+");
//FILE* vvectorfile=fopen("vvector.txt", "wt+");
//FILE* wvectorfile=fopen("wvector.txt", "wt+");
double** matrix=malloc(3* sizeof(double*)) ;
for( i=0; i<3; i++ )  matrix[i]=malloc(sizeof(double) *4);

    //Test if the ray can be constructed from the vectors that make up the planes of the surface  
    for(i=0; i < (n_meshPts-1) * 3; i++){
/*	if( i <( n_meshPts-1)*2){	
		fprintf( uvectorfile, "%f %f %f %f %f %f\n", local[i]->offset.x, local[i]->offset.y, local[i]->offset.z, local[i]->line1.x, local[i]->line1.y, local[i]->line1.z  );
		fprintf( vvectorfile, "%f %f %f %f %f %f\n", local[i]->offset.x, local[i]->offset.y, local[i]->offset.z, local[i]->line2.x, local[i]->line2.y, local[i]->line2.z  );
		fprintf( wvectorfile, "%f %f %f %f %f %f\n", local[i]->offset.x, local[i]->offset.y, local[i]->offset.z, local[i]->line3.x, local[i]->line3.y, local[i]->line3.z  );
        } 
	else {
		fprintf( uvectorfile, "%f %f %f %f %f %f\n", local[i]->offset.x, local[i]->offset.y, local[i]->offset.z, local[i]->line1.x, local[i]->line1.y, local[i]->line1.z  );
		fprintf( vvectorfile, "%f %f %f %f %f %f\n", local[i]->offset2.x, local[i]->offset2.y, local[i]->offset2.z, local[i]->line2.x, local[i]->line2.y, local[i]->line2.z  );
		fprintf( wvectorfile, "%f %f %f %f %f %f\n", local[i]->offset.x, local[i]->offset.y, local[i]->offset.z, local[i]->line3.x, local[i]->line3.y, local[i]->line3.z  );
	}
/**/

	/*minus 1 because exclude first point*/
	if( i <( n_meshPts-1)*2 )
	result=intersection( point, &local[i]->line1, &local[i]->line2, &local[i]->offset, inner_mesh, outer_mesh, i, n_meshPts-1 , rotate, matrix, 0);
	else{
	result=intersection( point, &local[i]->line1, &local[i]->line2, &local[i]->offset, inner_mesh, outer_mesh, i, n_meshPts-1 , rotate, matrix, 0);
		if(  result  == 1) intersections++;
        	else if (result == -1) { 
		//printf("found border\n"); 
		return 1;
		} //a -1 means that the point is on the border with the mesh and so we will count it as being inside. 
		else if (result == -2){
		return -2;
		}    
		else if (result == -3) extra=1; //one of the points is a vertex center
	result=intersection( point, &local[i]->line3, &local[i]->line4, &local[i]->offset3, inner_mesh, outer_mesh, i, n_meshPts-1 , rotate, matrix, 1);
	}

	if(  result  == 1) intersections++;
        else if (result == -1) { 
	//printf("found border\n"); 
	return 1;
	} //a -1 means that the point is on the border with the mesh and so we will count it as being inside. 
	else if (result == -2){
	return -2;
	}    
	else if (result == -3) extra=1; //one of the points is a vertex center
}
for( i=0; i<3; i++ )  free(matrix[i]);
free(matrix);
for( i=0; i < (n_meshPts-1)*3; i++) free(local[i]);
free(local);
//fclose(uvectorfile);
//fclose(vvectorfile);
//fclose(wvectorfile);
    intersections += extra;
   if(local_VERBOSE)  printf("Total of %d intersections\n", intersections); 
        
    if ( intersections%2 != 0) return 1; //the vertex is within the mesh because the number of intersections is odd
    else  return 0; //the vertex is outside of the mesh
   
  

}

