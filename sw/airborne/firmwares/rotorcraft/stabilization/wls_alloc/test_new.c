#include "wls_alloc.h"
#include "stdio.h"

#define M  3
#define N  4
void test(void); // test the previous main function in paparazzi


int main(){
	test();
	return 0;
}

void test(){
	// m = size(v), n = size(u)
	// int m = 3 ; int n = 4;

    float Wv[M] = {3, 3, 1};
    float B_tmp[M][N] = {{-21.5189e-3, 21.5189e-3, 21.5189e-3, -21.5189e-3},{14.3894e-3, 14.3894e-3, -14.3894e-3, -14.3894e-3},{ 1.2538e-3,  -1.2538e-3, 1.2538e-3, -1.2538e-3}};
			  
    float** B = (float**)calloc(M, sizeof(float*));
    for (int i = 0; i < M; i++) {
        B[i] = (float*)calloc(N, sizeof(float*));
        for (int j = 0; j < N; j++) B[i][j] = B_tmp[i][j];
    }
    float umax[N] = {500, 500, 500, 500};
    float umin[N] = {-500, -500, -500, -500};

    float v[M] = {0.00055796,-3.5578,2.535};
    float u[N] = {0, 0, 0, 0};

    wls_alloc(u,v,umin,umax,B,N,M,0,0,Wv,0,0,1000,100);

    for(int i = 0; i < N; i++)
        printf("%.2f\n", u[i]);
}
