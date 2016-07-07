#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "math/pprz_simple_matrix.h"

// void MAT_MUL(int i, int k, int j, float **C, float **A, float **B);

int main(){
	puts("Lets compute some matrices esse");
	float A[3][3] = {{1, 1, 1}, {2, 2, 2}, {3, 3, 3}};
	float B[3][3] = {{1, 1, 1}, {1, 1, 1}, {1, 1, 1}};
	int i = 3;
	int j = 3;
	int k = 3;
	// float **C = new float * [3];
	//for(int ii = 0; ii < 3; ii++)
		//C[ii] = new float[3];

	
	float ** C;
	C = calloc(j*i, sizeof(float*));
// #anton
//	int **own = malloc( mem_size * sizeof( int * ) );
//	for ( i = 0; i < mem_size; i++ ) own[i] = calloc( 3, sizeof( int ) );
	// C = calloc(j*i, sizeof(float));
	MAT_MUL(i, k ,j, C, A, B);

	//for(int ii = 0; ii < 3; ii++)
	//	delete[] C[ii];
	//delete[] C;
return 0;
}

