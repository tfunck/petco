#include <stdio.h> 
#include <stdlib.h>
#include <math.h>



/*int main()
{
int colSize=4, rowSize=3;
int i;
double** matrix=malloc(rowSize * sizeof(double*))  ;
for( i=0;i < rowSize; i++) matrix[i]=malloc(colSize * sizeof(double));
matrix[0][0]=0; matrix[0][1]=-1; matrix[0][2]=-9; matrix[0][3]=0; 
matrix[1][0]=0; matrix[1][1]=-1; matrix[1][2]=-9; matrix[1][3]=0; 
matrix[2][0]=0; matrix[2][1]=-1; matrix[2][2]=-4.5; matrix[2][3]=-.5; 

int kRow, nCol;
for( kRow =0; kRow < rowSize; kRow++)
{
     for( nCol = 0; nCol < colSize; nCol ++ ) 
     {
          printf("%1.3f\t", matrix[kRow][nCol]);
     }
printf("\n");
}

    
rref( matrix , colSize, rowSize);
return 0;
}*/ 
double** rref(double** matrix, int colSize, int rowSize)

{
int  nRow, kRow, nCol, iRow, divisorCol;
double divisor, multiplicand;

 /*   for( kRow =0; kRow < rowSize; kRow++){
        for( nCol = 0; nCol < colSize; nCol ++ ) printf("%1.3f\t", matrix[kRow][nCol]);
        printf("\n");
    }
     printf("\n"); 
*/	    
	for( nRow =0; nRow < rowSize; nRow++){
        for( nCol = 0; nCol < colSize; nCol ++ ){
            divisor = matrix[nRow][nCol];
            divisorCol=nCol;
            if( divisor != 0) break;
        }
        //printf("Divisor %f %f\n", divisor, fabs(divisor) );
       
        if(  fabsf(divisor) > .00001  ){   
            //Divide the elements of row i by the first element in the row
            for( nCol = 0; nCol < colSize; nCol ++ ) matrix[nRow][nCol]= matrix[nRow][nCol] / divisor ;
            for( kRow=0; kRow < rowSize; kRow++ ){
                multiplicand = matrix[kRow][divisorCol];
                if( kRow != nRow  ) for( nCol = 0; nCol < colSize; nCol ++ ) matrix[kRow][nCol] = matrix[kRow][nCol] - multiplicand*matrix[nRow][nCol]; 
            }
        }
    }

 /*   for( kRow =0; kRow < rowSize; kRow++)
        {
            for( nCol = 0; nCol < colSize; nCol ++ ) 
            {
                printf("%1.3f\t", matrix[kRow][nCol]);
            }
        printf("\n");
        }

/*printf("\n\n\n"); */
return matrix;
}

