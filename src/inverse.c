#include <stdio.h>
#include <stdlib.h>

void printMat(float *mat){
	for(int y=0; y<3; y++ ){
		for(int x=0; x< 3; x++){
			printf("%f\t", mat[y*3+x]);
		}
	printf("\n");
	}	
printf("\n");
}

int isin(int number, int* array, int dim){

    for(int i=0; i<dim; i++) if( number==array[i]) return 1;

return 0;

}

int inverse2(double* matrix, double* out, int rows, int* leading_rows,  int depth){
    int i,j,index, leading_index;
    int cols=rows;
    float alpha;
    int leading_row=-1;
    if(depth==0){
        for(i=0; i<9; i++) out[i]=0;
        out[0]=out[4]=out[8]=1;
    }
    //printf("\nDEPTH: %d\n", depth);
    if(depth < rows){
        for(i=0; i<rows; i++){//offset by depth so we start searching plus one down diagonal
            index=i*cols+depth;
            //printf("%d\n", isin(i, leading_rows, rows) );
            if ( matrix[index] != 0 && isin(i, leading_rows, rows)==0 ){
                //find leading term, call it alpha
                alpha=matrix[index];
                leading_row=i;
                leading_rows[depth]=leading_row;
                //printf("Alpha %f\n", alpha);
                //find first row with non-zero element, this is the leading row
                break;
            }
        }
        
        //if entire row is equal to 0, then we have a lin. dep. matrix
        //printf("Leading row: %d\n", leading_row);
        if(leading_row == -1){ 
            fprintf(stderr, "Error: Singular matrix.\n");
            exit(1);
        }
        //Divide leading row by leading term alpha
        for(j=0; j<cols; j++)   {
            index=leading_row*cols+j;
            out[index] /= alpha;
            matrix[index] /= alpha;
        }
        //printf("Inverse:\n");
        //printMat(out);
        //printf("matrix:\n");
        //printMat(matrix);
        for(i=0; i<rows; i++){
            if(i != leading_row){
                //first element of current row is beta
                index=i*cols+depth;
                float beta=matrix[index];    
                //printf("Beta: %f\n", beta);
                //subtract elements of row j by corresponding element in row i times beta
                
                for(j=0; j<cols; j++){
                    index=i*cols+j;
                    leading_index=leading_row*cols+j;
                    //printf("\t%d %d: %f %f\n",i,j, beta, matrix[leading_index]  );
                    out[index] -= beta * out[leading_index];
                    matrix[index] -= beta * matrix[leading_index];
                    
                }
            }
        }
    
    /*printf("inverse:\n");
    printMat(out);
    printf("matrix:\n");
    printMat(matrix);
    */

    return inverse2(matrix, out, rows, leading_rows, depth+1); 
    } 
    else{
     /*   alpha=matrix[(depth)*(depth+1)+depth]; //assign bottom right element of matrix to alpha
        printf("Alpha: %f\n", alpha); 
        for(j=0; j < rows; j++){
            matrix[(depth)*(depth+1)+j] /= alpha;
            out[(depth)*(depth+1)+j] /= alpha;
        }
                    printf("inverse:\n");
    printMat(out);
    printf("matrix:\n");
    printMat(matrix);

        

        for(i=0; i < depth; i++){
            index=i*rows+(cols-1) ;
            out[index] -= matrix[index]*out[depth*rows+depth] ;
            matrix[index] -= matrix[index]*matrix[depth*rows+depth] ;  // matrix[depth+1][j];
            
            printf("%f x %f\n",matrix[index], out[depth*rows+depth] );
            
        }
            printf("inverse:\n");
    printMat(out);
    printf("matrix:\n");
    printMat(matrix);*/

       return 0;
    }


}

/*int main(){

float data[]={1, 2, 3, 4, 5, 6, 7, 8, 12};
float output[]={1,0,0,0,1,0,0,0,1};

//for(int i=0; i < 400*400*500; i++) 
inverse2(data, 3, 0 );
printMat(data);


return(0);
}*/
