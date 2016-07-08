#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "math/pprz_simple_matrix.h"
#include "math/pprz_matrix_decomp_float.h"

#define numcola  3 // global niggu ?!?!?!

// void MAT_MUL(int i, int k, int j, float **C, float **A, float **B);

// Actually pass the pointer of the first element to the array ( float(*A)[numcols]), which points to multiple arrays
void print_matrix(float A[][numcola], int m, int n);


int main(){
	puts("Lets compute some matrices esse");
	float A[3][3] = {{1, 1, 1}, {2, 2, 2}, {3, 3, 3}};
	float B[3][3] = {{1, 1, 1}, {1, 1, 1}, {1, 1, 1}};
	int i = 3;
	int j = 3;
	int k = 3;
	float a = 2;
	float b = 3;
	float c;
	// float **C = new float * [3];
	//for(int ii = 0; ii < 3; ii++)
		//C[ii] = new float[3];
//	float ** C;
	c = pythag(a,b);
//	C = calloc(j*i, sizeof(float*));
// #anton
//	float **own = malloc(j*i * sizeof(float * ) );
//	for (i = 0; i < ; i++ ) C[i] = calloc(3, sizeof(float));
float ( *C)[i] = calloc( i*j, i*sizeof(float));
	// C = calloc(j*i, sizeof(float));
	MAT_MUL(i, k ,j, C, A, B);
	printf("%8.3g \n",C[1][1]);
	//for(int ii = 0; ii < 3; ii++)
	//	delete[] C[ii];
	//delete[] C;
//	print_matrix(C, i, j);
return 0;
}

void print_matrix(float A[][3], int m, int n)
/* Display contents of an m x n matrix. */
{
	  int i, j;
	    for (i=0 ; i<m ; i++) {
		 if (i==0) printf("[ ");
		   else printf("  ");
		for (j=0 ; j<n ; j++)
			printf("%8.3g ",A[i][j]);
			if (i < m-1) printf("\n");
			else printf("]\n");
	}
}
