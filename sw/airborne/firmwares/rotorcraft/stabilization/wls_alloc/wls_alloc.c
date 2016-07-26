// Script written by Anton Naruta && Daniel Hoppener 2016


#include "wls_alloc.h"
#include <stdio.h>

// the wrapper can use any solver function
float* qr_solve_wrapper(int m, int n, float** A, float* b) {
    float in[m * n];
    // convert A to 1d array 
    int k = 0;
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < m; i++) {
            in[k++] = A[i][j];
        }
    }
    // use solver
    return qr_solve(m, n, in, b);
}

int wls_alloc(float* u, float* v, float* umin, float* umax, float** B,
              int n_u, int n_w, float* u_guess, float* W_init, float* Wv,
              float* Wu, float* ud, float gamma_sq, int imax) {
    // allocate variables, use defaults where parameters are set to 0 
    if(!gamma_sq) gamma_sq = 100000;
    if(!imax) imax = 100;
    int n_c = n_u + n_w;

    float** A = (float**)calloc(n_c, sizeof(float*));
    float** A_free = (float**)calloc(n_c, sizeof(float*));

    float b[n_c];
    float d[n_c];

    int free_index[n_u];
    int free_index_lookup[n_u];
    int n_free = 0;
    int free_chk = -1;

    // provide loop feedback 
    int verbose = 1;

    int iter = 0;
    float* p_free;
    float p[n_u];
    float u_opt[n_u];
    int infeasible_index[n_u];
    int n_infeasible = 0;
    float lambda[n_u];
    float W[n_u];

    if (!u_guess) {
        for (int i = 0; i < n_u; i++) u[i] = (umax[i] + umin[i]) * 0.5;
    } else {
        for (int i = 0; i < n_u; i++) u[i] = u_guess[i];
    }
    W_init ? memcpy(W, W_init, n_u * sizeof(float))
           : memset(W, 0, n_u * sizeof(float));

    memset(free_index_lookup, -1, n_u * sizeof(float));


    // find free indices
    for (int i = 0; i < n_u; i++) {
        if (W[i] == 0) {
            free_index_lookup[i] = n_free;
            free_index[n_free++] = i;
        }
    }

    // allocate and fill up A, A_free, b and d
    for (int i = 0; i < n_w; i++) {
        A[i] = (float*)calloc(n_u, sizeof(float));
        A_free[i] = (float*)calloc(n_u, sizeof(float));
        b[i] = Wv ? gamma_sq * Wv[i] * v[i] : gamma_sq * v[i];
        d[i] = b[i];
        for (int j = 0; j < n_u; j++) {
            A[i][j] = Wv ? gamma_sq * Wv[i] * B[i][j] : gamma_sq * B[i][j];
            d[i] -= A[i][j] * u[j];
        }
    }
    for (int i = n_w; i < n_c; i++) {
        A[i] = (float*)calloc(n_u, sizeof(float));
        A_free[i] = (float*)calloc(n_u, sizeof(float));
        memset(A[i], 0, n_u * sizeof(float));
        A[i][i - n_w] = Wu ? Wu[i - n_w] : 1.0;
        b[i] = ud ? (Wu ? Wu[i] * ud[i] : ud[i]) : 0;
        d[i] = b[i] - A[i][i - n_w] * u[i - n_w];
    }

    while (iter++ < imax) {
        // clear u, copy u to u_opt
        memset(p, 0, n_u * sizeof(float));
        memcpy(u_opt, u, n_u * sizeof(float));

        if (free_chk != n_free) {
            for (int i = 0; i < n_c; i++) {
                for (int j = 0; j < n_free; j++) {
                    A_free[i][j] = A[i][free_index[j]];
                }
            }
            free_chk = n_free;
        }
	// print iteration
	if (verbose) {
		printf("Iteration %d \n",iter);
		printf("u = \n");
		for(int i = 0; i < n_u; i++)
        		printf("%.2f\n", u[i]);
    		}

	if (!n_free) {
                // No free variables left, all actuators saturated
                for (int i = 0; i < n_c; i++) {
                    free(A[i]);
                    free(A_free[i]);
                }
                free(A);
                free(A_free);
                return iter;
            }

        // use a solver to find solution to A_free*p_free = d
        p_free = qr_solve_wrapper(n_c, n_free, A_free, d);
        for (int i = 0; i < n_free; i++) {
            p[free_index[i]] = p_free[i];
            u_opt[free_index[i]] += p_free[i];
        }
        // check limits
        n_infeasible = 0;
        for (int i = 0; i < n_u; i++) {
            if (u_opt[i] > umax[i] || u_opt[i] < umin[i]) {
                infeasible_index[n_infeasible++] = i;
            }
        }
        if (n_infeasible == 0) {
            // all variables are within limits
            memcpy(u, u_opt, n_u * sizeof(float));
            memset(lambda, 0, n_u * sizeof(float));

            // d = d + A_free*p_free; lambda = A*d;
            for (int i = 0; i < n_c; i++) {
                for (int k = 0; k < n_free; k++) {
                    d[i] += A_free[i][k] * p_free[k];
                }
                for (int k = 0; k < n_u; k++) {
                    lambda[k] += A[i][k] * d[i];
                }
            }
            bool break_flag = true;
            int i_neg;

            // lambda = lambda x W;
            for (int i = 0; i < n_u; i++) {
                lambda[i] *= W[i];
                // if any lambdas are negative, keep looking for solution
                if (lambda[i] < -DBL_EPSILON) {
                    break_flag = false;
                    W[i] = 0;
                    // add a free index
                    if (free_index_lookup[i] < 0) {
                        free_index_lookup[i] = n_free;
                        free_index[n_free++] = i;
                    }
                }
            }
            if (break_flag) {
                // if solution is found de-allocate A, A_free and return number of iterations
                for (int i = 0; i < n_c; i++) {
                    free(A[i]);
                    free(A_free[i]);
                }
                free(A);
                free(A_free);
                return iter;
            }
        } else {
            float alpha = INFINITY;
            float alpha_tmp;
            int id_alpha;

            // find the lowest distance from the limit among the free variables
            for (int i = 0; i < n_free; i++) {
                int id = free_index[i];
                alpha_tmp = (p[id] < 0) ? (umin[id] - u[id]) / p[id]
                                        : (umax[id] - u[id]) / p[id];
                if (alpha_tmp < alpha) {
                    alpha = alpha_tmp;
                    id_alpha = id;
                }
            }

            // update input u = u + alpha*p
            for (int i = 0; i < n_u; i++) {
                u[i] += alpha * p[i];
            }
            // update d = d-alpha*A*p_free
            for (int i = 0; i < n_c; i++) {
                for (int k = 0; k < n_free; k++) {
                    d[i] -= A_free[i][k] * alpha * p_free[k];
                }
            }
            // get rid of a free index
            W[id_alpha] = (p[id_alpha] > 0) ? 1.0 : -1.0;

            free_index[free_index_lookup[id_alpha]] = free_index[--n_free];
            free_index_lookup[free_index[free_index_lookup[id_alpha]]] =
                free_index_lookup[id_alpha];
            free_index_lookup[id_alpha] = -1;
        }
    }
    // solution failed, de-allocate A and A_free
    for (int i = 0; i < n_c; i++) {
        free(A[i]);
        free(A_free[i]);
    }
    free(A);
    free(A_free);
    // return negative one to indicate failure
    return -1;
}
