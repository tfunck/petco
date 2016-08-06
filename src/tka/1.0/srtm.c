#include <stdio.h>
#include <math.h>
#include "convolv.h" 
#include "minc2.h" 
#include "srtm.h"  
#include "mtga.h"




void srtmBasisFns(float* Cr, int nframes, float* weights, float **basis, tracer* k2parameters, int basisMax, vectorArg **Q_array, vectorArg **R_array, float **theta3Array, imageData imageInfo)
{
int i, j, k;// t;
int time;
float *times=imageInfo.times;
float *timeWidths=imageInfo.timeWidths;
float *referenceArray=imageInfo.referenceArray;
float *response;
float *ans;
//float theta1, theta2,
float theta3;
//float output =0;
//float min=0;
//float R1;
//float BP;
//float totaltime=0.0;
float rate;
//float BPmax=3;
//float time2=0;
//float k2; 
float k2Min= (float)k2parameters->k2min , k2Max=(float)k2parameters->k2max;
float theta3max, theta3min;
float lambda= (float)k2parameters->lambda ;
float e_theta3;
vectorArg argsQ, argsR, argsA;
int nsamples;
float *contREF;
float average;
float *convolution;
nsamples=(int)times[nframes-1]+ (int)timeWidths[nframes-1];    
contREF=malloc(nsamples*sizeof(float));  
response=malloc(nsamples*sizeof(float));   
ans= malloc(2*nsamples*sizeof(float));
 
 *theta3Array=malloc(basisMax*sizeof(float));
 *Q_array=malloc(basisMax*sizeof(vectorArg) );
 *R_array=malloc(basisMax*sizeof(vectorArg) );
    
    argsA.matrix=malloc(2 * sizeof(float *));
    for(j=0; j<2; j++) argsA.matrix[j]=malloc( nframes * sizeof(float) ); 
     
//printf("\tCr\n");
//for (i=0 ;i<nframes ;i++ ) printf("\t%f\n", Cr[i] ); 

   /* theta3Max = k2Max+lambda, theta3Min= k2Min/(1+BPMax)*/
   theta3max=k2Max+lambda;
   theta3min=lambda;
   rate=(theta3max - theta3min)/basisMax;
   
   //rate=exp(rate);

   printf("theta3min=%f\ttheta3max=%f\n", theta3min, theta3max);  
//printf("Rate=%f\tlambda=%f\ttheta3min=%f\ttheta3Max=%f\n", rate, lambda, theta3min, theta3max);
    
   printf("nsamples: %d", nsamples);
   for(k=0; k< nsamples; k++) 
       for(j=0; j<nframes; j++) 
        {
            if (times[j] > k) {contREF[k]=referenceArray[j];  break;} 
            else if ( times[nframes-1] < k) {contREF[k]=referenceArray[nframes-1]; break;} 
        }

    for(i=1; i <= basisMax; i++)
	{ 

    argsR.matrix=malloc(2 * sizeof(float*));
    for(j=0; j< 2; j++) argsR.matrix[j]=malloc(2 * sizeof(float) );        
    argsQ.matrix=malloc(2 * sizeof(float*));
    for(j=0; j< 2; j++) argsQ.matrix[j]=malloc(nframes * sizeof(float) );  
                         
   /*1: define theta3*/
   theta3= i*rate + theta3min;
   (*theta3Array)[i-1]=theta3;
  // printf("%d: Theta3=%f\n", i, theta3);
  
 
   /*Create response theta3 function*/
   for( j=0; j< nsamples; j++ )  { response[j]=exp(-1*theta3*j);   /*printf("response[%d]=%f\n", j, response[j]);*/ }   
   
   /*Create continuous version of discrete reference data*/

     //for(k=0; k< nsamples; k++) printf("test %f\n", contREF[k]);
   e_theta3=exp(-1*theta3);    

   /*calulcate convolution according to Gunn, et al. 1997 method (from Gunn matlab code)*/
   ans[0]=contREF[0]*(1-e_theta3)/theta3; 
   for(k=1; k < nsamples; k++) ans[k]=ans[k-1]*e_theta3+contREF[k]*(1-e_theta3)/theta3;    

   /*average output of pseudo continuous convolution into discrete data according to frames*/
   convolution=malloc(nframes*sizeof(float));
   for(j=0; j<nframes; j++) 
   {     
         average=0;
         for(k=(int) referenceArray[j]; k< (int) (referenceArray[j]+timeWidths[j]); k++) {average += ans[k];
            
         }              
         convolution[j]=average/timeWidths[j];
         printf("%d\t%f\t%f\t%f\n", j, Cr[j], convolution[j], average);
   }
   /*4: set matrix A with reference tissue and output function from convolution */   
   for(j=0; j<nframes; j++) 
   {
       argsA.matrix[0][j]=Cr[j]  * sqrt(weights[j]);
       argsA.matrix[1][j]=convolution[j] * sqrt(weights[j]);
       //printf("%f %f\n", Cr[j]*sqrt(weights[j]) , convolution[j]*sqrt(weights[j]) );
   }

   argsQ.start=0; argsQ.row=nframes; argsQ.col=2;
   argsR.start=0; argsR.row=2; argsR.col=2;                
   /*5: perform qr decomposition on matrix A*/
   qrdecomp(&argsA, &argsQ, &argsR, 2, nframes);
         
   (*Q_array)[i-1]=argsQ;
   (*R_array)[i-1]=argsR;
   basis[i-1]=convolution;
   }
          
for(i=0; i<2; i++) free(argsA.matrix[i]);
free(argsA.matrix);
}


float srtmBasisResult(float *Ct, float *Cr,  int nframes, float* weights,  int basisMax,double  lambda, vectorArg *Q_array, vectorArg *R_array, float **basis, float *theta3Array, vectorArg *argsB, vectorArg *argsQt, int *numMax )  
{
int i, j, k, t, iMin=0;
int lowestIndex;                             
float min;
float R1, k2, BP;
float SSR=0;
float CtWeighted[nframes];
float theta1, theta2, theta3, Ftheta1, Ftheta2, Ftheta3;
   for(j=0; j< nframes; j++ ) CtWeighted[j]=Ct[j]* sqrt(weights[j]);
 solve_qr( Q_array[0].matrix, R_array[0].matrix, CtWeighted, &theta1, &theta2, 2, nframes, argsB, argsQt); 
 for(j=0; j< nframes; j++ )   SSR +=  weights[j] *    pow(Ct[j]- (  theta1*Cr[j] + theta2* basis[0][j] )  , 2); 
// printf("%f = %f -  %f x  %f + %f x %f\n", SSR,  Ct[j], theta1, Cr[j], theta2, basis[0][j]);      
                             
 min=SSR;
 Ftheta1=theta1;
 Ftheta2=theta2;
 Ftheta3=theta3Array[0];

 //printf("%d\tSSR: %f\tTheta3: %f\n", 0, SSR, theta3Array[0]);       
 SSR=0;
 for(i=1; i<basisMax; i++)
     {  
        solve_qr(Q_array[i].matrix, R_array[i].matrix, CtWeighted, &theta1, &theta2, 2, nframes, argsB, argsQt); 
        /*7: calculate sum of squared differences */
        for(j=0; j< nframes; j++ )   SSR +=  weights[j] *  pow(Ct[j]- (theta1*Cr[j] +theta2*basis[i][j] ) , 2);
           //printf("%d\t%f = %f -  %f x  %f + %f x %f\n", j, SSR,  Ct[j], theta1, Cr[j], theta2, basis[i][j]);
        //printf("%d\tSSR: %f\tTheta3: %f\n", i, SSR, theta3Array[i]);
        if (SSR < min) 
        {   // printf("New Min: %d %f %f %f\n", i, SSR, min, theta3Array[i] );
            Ftheta1=theta1;
            Ftheta2=theta2;
            Ftheta3=theta3Array[i];
            iMin=i;
            min=SSR; 
        }
        SSR=0;     
     }
     
     R1=Ftheta1;
     k2=Ftheta2+R1*(Ftheta3- ( (float)lambda / 1) );
     //printf("%f / ( %f - %f) - 1\n", k2, Ftheta3, ( (float) lambda / 60 )  );
     BP=( (k2/(Ftheta3-( (float)lambda/ 1) ) ) -1);        


      if(iMin==basisMax-1) *numMax += 1;

   //printf("results after i=%d\n", iMin );  
   //printf("i:%d\ttheta3=%f\tBP=%f\ttheta1=%f\ttheta2=%f\tk2=%f\n\n\n\n", iMin, Ftheta3, BP, Ftheta1, Ftheta2, k2);     
   
  //    for(j=0; j< nframes; j++ )
  //     {
  //     SSR +=  weights[j] *   pow(Ct[j]- (Ftheta1*Cr[j] +Ftheta2*basis[iMin][j] ) , 2);        
  //    printf("%d\t%f = %f -  %f x  %f + %f x %f\n", j, SSR,  Ct[j], Ftheta1, Cr[j], Ftheta2, basis[iMin][j]);
  //   }
        
   // printf("i:%d\ttheta3=%f\tBP=%f\ttheta1=%f\ttheta2=%f\tk2=%f\n\n\n\n", iMin, Ftheta3, BP, Ftheta1, Ftheta2, k2); 
       
    return  BP;
}

void qrdecomp(vectorArg *input, vectorArg *Q, vectorArg *R, int col, int row)
{
  int i, j, k, x;
  vectorArg tempVector;
  vectorArg tempQ;    

tempVector.matrix=malloc(sizeof(float*) ); 
tempVector.matrix[0]=malloc(row * sizeof(float));

tempQ.matrix=malloc(sizeof(float *));
tempQ.matrix[0]=malloc(row*sizeof(float));   
    
for(i=0; i < col; i++)
    {
        for(x=0; x<row; x++) 
        {
            tempVector.matrix[0][x]=0.0;
            tempQ.matrix[0][x]=0.0;
        }
    /*Initialize temporary vector to 0*/
    //printf("\n\ni=%d\n",  i);
   
    /*Set the first non-zero value of R for column j, to the magnitude of the input column for col j*/
         
        if(i==0)
        {
            R->matrix[i][i]=magnitude(input->matrix[i], row);
            //printf("R.matrix[%d][%d]=%f\n", i, i, R->matrix[i][i] );
            input->row=row;
            input->col=i;
            input->start=i;
            Q->row=row;
            Q->col=i;
            Q->start=i;
            scalar_mult( (float) 1/R->matrix[i][i], input, Q); /*First column of Q*/
            //for(k=0; k<row; k++) printf("\tQ[0][%d]= %f = %f x %f\n ", i, Q->matrix[i][k], (float) 1/R->matrix[i][i],  input.matrix[0][k]);    
        }
        else
        {    
        tempVector.row=row;
        tempVector.col=0;
        tempVector.start=0;   

        tempQ.row=row;
        tempQ.col=0;
        tempQ.start=0;     
       
            for(j=0; j < i; j++  )
            {
              //  printf("j=%d\n", j);
            
                /*Allocate memory for temporary Q column*/
                Q->row=row;
                Q->col=j;
                Q->start=j;


                /* Multiply r[current col][current row] * Q[current col]  */
                scalar_mult( R->matrix[i][j], Q, &tempQ);

                /*Add the temporary Q, which has been multiplied by r scalar, to the sum of Q vectors*/
                matrix_add( &tempVector, &tempQ,  &tempVector, 1);
            //        for(k=0; k<row; k++) printf("3. tempVector[0][%d]  %f\n", k, tempVector.matrix[0][k] );  
            }  
        
        input->row=row;
        input->col=i;
        input->start=i;

        Q->row=row;
        Q->col=i;
        Q->start=i;

        matrix_add(input, &tempVector,  Q, 0);
//        for(k=0; k<col; k++) printf("After add, before mult: q[%d][%d]=%f\n", i, k, Q->matrix[i][k]);         
                   
        R->matrix[i][i]=magnitude(Q->matrix[i], row);    
        //printf("R.matrix[%d][%d]=%f\n", i, i, R->matrix[i][i] );
        scalar_mult( (float) 1.0/R->matrix[i][i], Q, Q); 
        /*for(k=0; k<row; k++) printf("\tq[%d][%d]=%f\n", i, k, Q->matrix[i][k]);*/
        }
       
        for(k=0; k<col; k++)
        {   
            if(k<i) R->matrix[k][i]=0.0;
            else if (k>i)
            {   
               // printf("k=%d\t", k );
                input->start=k;
                input->col=k;
                /*for(x=0; x<row; x++  ) printf("test input %f\n", input.matrix[k][x]  );*/
                R->matrix[k][i]=dot(Q, input);
          //  printf("R.matrix[%d][%d]=%f\n", i, i, R->matrix[k][i] );
            } 
    //    printf("\tR[%d][%d]=%f\t", i, k, R->matrix[i][k]); 
        }
    //printf("\n");  
    /*printf("\nQ\n, %d column:",i);
       for(k=0; k<=i; k++) 
        {
             printf("\n"); 
             for(j=0; j<row; j++) printf("Q[k][j]%f\n", Q->matrix[k][j] );  */ 
       // }       
    }
 
//   printf("\nR");
///        for(i=0; i<col; i++) 
//        {
//            printf("\n"); 
//            for(k=0; k<col; k++) printf("%f\t", R->matrix[k][i] );   
//        }  
//   printf("\n");

 
free(tempVector.matrix[0]);
free( tempVector.matrix );

free(tempQ.matrix[0]);
free( tempQ.matrix );
} 
         
void solve_qr(float** Q, float** R, float *Ct,  float* theta1, float* theta2, int col, int row, vectorArg *argsB, vectorArg *argsQt)
{
int i, j;
vectorArg argsCt;

argsB->row=row;
argsB->col=col;
argsB->start=0;

//argsB->matrix=malloc( sizeof(float *) );
//argsB->matrix[0]=malloc( 2 * sizeof(float) );
//printf("col %d\trow %d\n", col, row);

//argsQt->matrix=malloc(row *  sizeof(float *));
//for(i=0; i< row; i++) argsQt->matrix[i]=malloc( col * sizeof(float) );
argsQt->row=col;
argsQt->col=row;
argsQt->start=0;

//argsCt.matrix=malloc( sizeof(float*));
argsCt.matrix=&Ct;
argsCt.row=row;
argsCt.col=1;
argsCt.start=0;


transpose(Q, col, row, argsQt->matrix);
matrix_mult(argsQt, &argsCt, argsB);

/*
printf("\n\nQ"); 
for(i=0; i<row; i++) 
	{
		
		for(j=0; j< col; j++) printf("%f ", Q[j][i]);
 printf("\n");  
	} 
 	
printf("\n\nR"); 
for(i=0; i<col; i++) 
	{
	printf("\n"); 
	for(j=0; j< col; j++) printf(" %f ", R[j][i]);
	} 
printf("\n");
   
printf("\n\nargsQt->"); 
for(i=0; i<col; i++) 
	{
	printf("\n"); 
	for(j=0; j< row; j++) printf("%f ", argsQt->matrix[j][i]);
	}

    printf("\nB\n%f\n%f\n\n", argsB->matrix[0][0], argsB->matrix[0][1] );
  */

*theta2= argsB->matrix[0][1] / R[col-1][col-1];
//printf("theta2 %f / %f = %f\n",  argsB->matrix[0][1],  R[col-1][col-1],  argsB->matrix[0][1] / R[col-1][col-1] );

*theta1= (argsB->matrix[0][0]-R[1][0] * (*theta2) ) / R[0][0];
//printf("theta1  ( %f -  %f * %f ) / %f = %f\n", argsB->matrix[0][0], R[1][0], (*theta2),  R[0][0],  (argsB->matrix[0][0]-R[1][0] * (*theta2) ) / R[0][0]  );  

//printf("col %d\trow %d\n", col, row);




}

void transpose(float **A, int col, int row, float** B)
{
int i, j;

for(i=0; i<col; i++) for(j=0; j<row; j++)  B[j][i]=A[i][j];  
}


float magnitude(float *A, int row)
{
int i;
float sum=0.0;

    for(i=0; i< row; i++) {sum += A[i]*A[i]; /*printf("magnitude[%d]=%f. Matrix[%d]=%f\n", i, sum, i, A[i] ); */}
    //printf("Magnitude^2=%f\n", sum);
    sum=sqrt(sum); 
   // printf("Magnitude=%f\n", sum);
return  sum;
}



void normalize(float **A,int row, int col, float *sum )
{
  int i;
*sum=0;

for(i=0; i< row; i++) *sum += A[col][i]*A[col][i];
*sum=sqrt(*sum);

for(i=0; i<row; i++) A[col][i] /= *sum;

} 



float dot(vectorArg *argsA, vectorArg *argsB)
{
  int i;
  float product=0;
  float **A=argsA->matrix; 
  float **B=argsB->matrix;    
  int startA=argsA->start;
  int startB=argsB->start;
  int n=argsA->row;

    for (i=0; i<n; i++) 
    {
    product += A[startA][i]*B[startB][i];
  //  printf("\nA[%d][%d] %f * B[%d][%d] %f = %f\n", startA, i, A[startA][i], startB, i, B[startB][i], A[startA][i]*B[startB][i]  );
    }
return product;
}   

void matrix_mult(vectorArg *argsA, vectorArg *argsB, vectorArg *argsC)
{
  int i, j, k=0;
  float **C=argsC->matrix;
  float **A=argsA->matrix;
  float **B=argsB->matrix;  
  float product;
  int colA=argsA->col;
  int rowA=argsA->row;
  int colB=argsB->col;
  int rowB=argsB->row;
  

if(colA != rowB) {fprintf(stderr,"Error: rowA does not equal colB!\n"); exit(1);}
//printf("rowA=%d\tcolB=%d\trowB=%d\n", argsA.row, argsB.col, argsB.row );  
for(j=0; j<rowA; j++)
  {
    for(k=0; k<colB; k++)
    {
        for(i=0; i < colA; i++) 
        {
            product += A[i][j] * B[k][i]; 
            //printf("A: %d x %d\tB: %d x %d = %f\n", i, j, j,k, product  );
        }
        argsC->matrix[k][j]=product;
        product=0;
    }
  }

}

void scalar_mult(float scalar, vectorArg *argsA, vectorArg *argsB)
{
  int i, j, k;    
  float **A=argsA->matrix;
  float **B=argsB->matrix;  
  int colA=argsA->col;
  int rowA=argsA->row;
  int colB=argsB->col;
  int rowB=argsB->row;
  int startA=argsA->start;
  int startB=argsB->start;

 // printf("startA: %d, colA: %d\n", startA, colA );                
  for(i=startA, k=startB; i<colA+1; k++, i++)
  { 
       for(j=0; j<rowA; j++) 
       {
           B[k][j]=scalar * A[i][j]; 
     //      printf("scalar_mult: %d %d %f %f\n", i, j, argsA.matrix[i][j], scalar*A[i][j] );
       }
  }
 
}


void matrix_add(vectorArg *argsA, vectorArg *argsB, vectorArg *argsC, int operation)
{
  int i, j, l, k=0;
  float product;
  float **C=argsC->matrix;
  float **A=argsA->matrix;
  float **B=argsB->matrix;  
  int colA=argsA->col;
  int rowA=argsA->row;
  int colB=argsB->col;
  int rowB=argsB->row;    
  int startA=argsA->start;
  int startB=argsB->start;    
  int startC=argsC->start;

  if(colA-startA != colB-startB || rowA != rowB) {fprintf(stderr,"Error: two matrices are not of same dimensions and cannot be added.\n"); exit(1);}

  for(i=startA, l=startB,  k=startC; i<=colA; i++, k++, l++)
  {  // printf("i=%d, k=%d\n",i,k );
    for(j=0; j<rowB; j++)
    {
        if(operation == 1) C[k][j] = A[i][j] + B[l][j];
        else C[k][j]=A[i][j]-B[l][j];
    //printf("Matrix Add - C[%d][%d]=%f, %f and %f\n", i, j, C[k][j], A[i][j], B[l][j] );
    }
  }

} 

