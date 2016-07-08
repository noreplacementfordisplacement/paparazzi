/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
WLSc_r control allocator

By: Daan Hoppener 2016

Script to be executed from command terminal:

Compile: gcc wlsc_paparazzi2.c -llapack -lblas -lm -o wlsc_paparazzi
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
#include <stdio.h>
#include <stdlib.h>
#include <math.h> /* sqrt */
#define PREC 1e-8

/*compile on linux so with UNDERSCORE*/
#define SGEMM sgemm_
#define SGEMV sgemv_
#define SCOPY scopy_
#define SGELS sgels_
#define SSCAL sscal_
/*#else
#define SGEMM dgemm
#define SGEMV dgemv
#define SCOPY dcopy
#define SGELS dgels
#define SSCAL dscal
#endif*/

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * ~~~          Subroutine prototypes          ~~~
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
void insert_Wu(float *Wu, float *A, int k, int m);
int copy_columns(float *A, float *B, float *W, int m, int n);
void add_vector_cond(float *a, float *b, float *W, int m);
int is_feasible(float *u, float *umin, float *umax, float *W, int m);
int lagrange_mult(float *lambda_ns, float *W, int m);
float max_step(float *u, float *p_free, float *umin, float *umax,
	       float *W, int m, int *i_block, int *if_block);
void print_vector(float *b, int n);
void print_matrix(float *A, int m, int n);
void wls_solve(float *B, float *v, float *umin, float *umax, float *Wv, float *Wu, float *ud, float gam, float *u, float *W, int imax);


int main(){
 	float v[3] = {3,3,3};
	
	/* B = G1+G2 */
	float B[12] =  {-21.5198e-3, 14.3894e-3, 1.2538e-3, 21.5198e-3, 14.3894e-3, -1.2538e-3, 21.5198e-3, -14.3894e-3, 1.2538e-3,-21.5198e-3, -14.3894e-3, -1.2538e-3};

	float Wv[9] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 5.0};
	
/*	float Wu[4] = {1.0, 1.0, 1.0, 1.0};*/
	float Wu[16] = {1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0};
	float umin[4] = {-1000, -1000, -1000, -1000};
	float umax[4] = {500, 500, 500, 500};
	float ud[4] = {0.0, 0.0, 0.0, 0.0};

	float gam = 100000;

	float u[4] = {-1000, 0.0, 0.0, 0.0};
	float W[4] = {-1.0, 0.0, 0.0, 0.0};
	int imax = 100;

// Call WLS solve:	
	wls_solve(B,v,umin,umax,Wv,Wu,ud,gam,u,W,imax);
return 0;
}
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * ~~~              Core function              ~~~
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
void wls_solve(float *B, float *v, float *umin, float *umax, float *Wv, float *Wu, float *ud, float gam, float *u, float *W, int imax){
    /* Input variables */ 
	/*float *B,*v,*umin,*umax,*Wv,*Wu,*ud,gam,*u0,*W0;*/

  /* Output variables */
	/*	float *u,*W;*/

  /* Other variables */
	int i;
	/*int i,k,m;*/
	int m = 4;
	int k = 3;
	/* Variables */
	 float mu, *A, *d, *A_free, *p_free, *u_opt, *lambda_ns, alpha, *v_r;
 	 int i_block, if_block, i_neg, iter, m_free, verbose = 0;

  /* Variables used in BLAS calls */
  float done = 1.0, dzero = 0.0, dmone = -1.0, *work, *A_tmp;
  int km = k+m, ione = 1, info, lwork, tmp;
  char *cht = "T", *chn = "N", *chl = "L";
  if (imax<0) {
    /* Turn on verbosity (output info to command window). */
    imax = -imax;
    verbose = 1; /* verbose = 2 for even more details */
  }

  verbose = 2;
  /* Allocate variables. */
  A = (float *)calloc((k+m)*m, sizeof(float));
  d = (float *)calloc(k+m, sizeof(float));
  v_r = (float *)calloc(k, sizeof(float));
  A_free = (float *)calloc((k+m)*m, sizeof(float));
  A_tmp = (float *)calloc((k+m)*m, sizeof(float));
  p_free = (float *)calloc(k+m, sizeof(float));
  u_opt = (float *)calloc(m, sizeof(float));
  lambda_ns = (float *)calloc(m, sizeof(float));

  /* Work area for least squares solves SGELS */
  tmp = 64; /* Generous block size, see
	       http://www.netlib.org/lapack/lug/node120.html#secblocksize*/
  lwork = m+km*tmp;
  work = (float *)calloc(lwork,sizeof(float));

  /* Reformulate problem:
   *    min ||Wu(u-ud)||^2 + gam ||Wv(Bu-v)||^2 <=>
   *    min ||mu Wv(Bu-v)|| = min ||Au-b||
   *        ||   Wu(u-ud)||
   * where mu^2 = gam, A = [mu Wv*B ; Wu], b = [mu Wv*v ; Wu*ud] 
   *
   * Change of variables: u = u0+p -> ||Au-b|| = ||Ap-d||
   * where d = b - A*u0
   */
  
  mu = sqrt(gam);
  
  /* Create A = [mu*Wv*B ; Wu] */
  SGEMM(chn,chn,&k,&m,&k,&mu,Wv,&k,B,&k,&dzero,A,&km);
  insert_Wu(Wu,A,k,m)
	
	print_matrix(A, km, m); /* verified! CHECK!!! */

  /* Create b = [mu Wv*v ; Wu*ud], store in d  */
  SGEMV(chn,&k,&k,&mu,Wv,&k,v,&ione,&dzero,d,&ione);
  SGEMV(chn,&m,&m,&done,Wu,&m,ud,&ione,&dzero,&d[k],&ione);

	print_matrix(d, km, ione);  /* Verified!! */ 

  
  /* Create d = b - A*u0 */
  SGEMV(chn,&km,&m,&dmone,A,&km,u,&ione,&done,d,&ione);
/*	puts("Holla at ya boi d = b - A*u0  is now:" );
	print_matrix(d, km, ione); */ /* Verified!! */
  /* ----------------------------------------------------------
   *  Iterate until optimum is found, or the maximum number of
   *  iterations is reached.
   * ---------------------------------------------------------- */
  
  for (iter=1 ; iter<=imax ; iter++) {

    if (verbose) {
      /* Output iteration info */
      printf("------------------------------------------\n");
      printf("Iteration: %d\n",iter);
      printf(" u = ");
      print_vector(u,m);
      printf(" W = ");
      print_vector(W,m);
    }
    
    /* ----------------------------------------
     *  Compute optimal perturbation vector p.
     * ---------------------------------------- */

    /* Eliminate saturated variables (create A_free). */
    m_free = copy_columns(A,A_free,W,km,m);
    if (verbose > 1) {
      printf(" A_free =\n");
      print_matrix(A_free,km,m_free);
      printf(" d = ");
      print_vector(d,km);
    }

    /* Solve the reduced optimization problem for the free variables, 
       i.e., compute least squares solution p_free = A_free\d. */
    SCOPY(&km,d,&ione,p_free,&ione);      /* p_free = d */
    tmp = km*m_free;
    SCOPY(&tmp,A_free,&ione,A_tmp,&ione); /* A_tmp = A_free */
                                           /* p_free = A_tmp\p_free */
    SGELS(chn,&km,&m_free,&ione,A_tmp,&km,p_free,&km,work,&lwork,&info);

    if (verbose) {  
      printf(" p_free = ");
      print_vector(p_free,m_free);
    }
    
    /* ----------------------------
     *  Is the new point feasible?
     * ---------------------------- */

    /* Create u_opt = u+p */
    SCOPY(&m,u,&ione,u_opt,&ione);    /* u_opt = u */
    add_vector_cond(p_free,u_opt,W,m); /* u_opt += p */
    if (verbose) {
      printf(" u + p = ");
      print_vector(u_opt,m);
    }

    if (is_feasible(u_opt,umin,umax,W,m)) {
      
       /* ----------------------------
       *  Yes, check for optimality.
       * ---------------------------- */

      if (verbose) {
	printf("Feasible w.r.t. constraints.\n");
      }

      /* Update point and residuals. */
      SCOPY(&m,u_opt,&ione,u,&ione); /* u = u_opt */
				      /* d = d - A_free*p_free */
      SGEMV(chn,&km,&m_free,&dmone,A_free,&km,p_free,&ione,&done,d,&ione);

      /* Compute Lagrange multipliers, not considering the sign due to
	 upper/lower bound. lambda_ns = A'*d */
      SGEMV(cht,&km,&m,&done,A,&km,d,&ione,&dzero,lambda_ns,&ione);
      /* i_neg is the index of the most negative lambda = W.*lambda_ns */
      i_neg = lagrange_mult(lambda_ns,W,m);
      
      if (i_neg<0) { /* All Lagrange multipliers positive */
	/*------------------------ \
	| Optimum found, bail out. |
	\ ------------------------*/
	if (verbose) {
	  printf("Optimum reached. Optimal solution:\n u = ");
	  print_vector(u,m);
	  printf("Control Objective:\n v = ");
	  print_vector(v,k);
	  printf("Realized Control Objective:\n v_r = ");
	  SGEMV(chn,&k,&m,&done,B,&k,u,&ione,&dzero,v_r,&ione);
	  print_vector(v_r,k);
	  printf("------------------------------------------\n");
	}

	/* Clean up */
	free(A);
	free(d);
	free(A_free);
	free(A_tmp);
	free(p_free);
	free(u_opt);
	free(lambda_ns);
	free(work);
	
	/* Return nr of iterations */
	return iter;
    }
    
      /* --------------------------------------------------
       *  Optimum not found, remove one active constraint.
       * -------------------------------------------------- */

      /* Remove constraint with most negative Lagrange multiplier
	 lambda from working set. */
      W[i_neg] = 0;
      if (verbose) {
	printf("Removing variable %d from working set.\n",i_neg+1);
      }
      
    } else /* is_feasible */ {
      
      /* ---------------------------------------
       *  No, find primary blocking constraint.
       * --------------------------------------- */

      if (verbose) {
	printf("Infeasible w.r.t. constraints.\n");
      }
      
      /* Compute primary blocking constraint (i_block, if_block) and
	 the maximum feasible step length along p (alpha). */
      alpha = max_step(u,p_free,umin,umax,W,m,&i_block,&if_block);
      
      /* Update point and residual. */
      SSCAL(&m_free,&alpha,p_free,&ione); /* p_free = alpha*p_free */
      add_vector_cond(p_free,u,W,m);       /* u += p */
      SGEMV(chn,&km,&m_free,&dmone,A_free,&km,
	     p_free,&ione,&done,d,&ione);  /* d -= A_free*p_free */
      
      /* Include the blocking constraint in the working set. 
	 Lower bound: -1, upper bound: +1 */
      W[i_block] = (p_free[if_block]>0) ? 1 : -1;
      
      if (verbose) {
	if (verbose>1) {
	  printf("Scaled perturbation: p_free = ");
	  print_vector(p_free,m_free);
	}
	printf("Adding variable %d (%s) to working set.\n",i_block+1,
	       (W[i_block]==1 ? "max" : "min"));
	printf("Reduced step length: alpha = %.3g\n",alpha);
      }

    } /* is_feasible */

  } /* for-loop */

  if (verbose) {
    printf("Max nr of iterations reached. Suboptimal solution:\n u = ");
    print_vector(u,m);
    printf("------------------------------------------\n");
  }

  /* Free allocated variables */
  free(A);
  free(d);
  free(A_free);
  free(A_tmp);
  free(p_free);
  free(u_opt);
  free(lambda_ns);
  free(work);
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * ~~~             Gateway routine             ~~~
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

/*void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]) */

  /* Check nr of inputs */
 /* if (nrhs != 11) {
    puts("Wrong number of input arguments.");
  }*/

  /* Create a pointers to input variables */
 /* B	  = mxGetPr(prhs[0]);
  v	  = mxGetPr(prhs[1]);
  umin	  = mxGetPr(prhs[2]);
  umax	  = mxGetPr(prhs[3]);
  Wv	  = mxGetPr(prhs[4]);
  Wu	  = mxGetPr(prhs[5]);
  ud	  = mxGetPr(prhs[6]);
  gam	  = mxGetScalar(prhs[7]);
  u0	  = mxGetPr(prhs[8]);
  W0	  = mxGetPr(prhs[9]);
  imax	  = mxGetScalar(prhs[10]);*/
  
  /* Number of variables */
/*  k = mxGetM(prhs[0]);*/
/*  m = mxGetN(prhs[0]);*/

  /* Allocate space for output variables */
/*  plhs[0] = mxCreateDoubleMatrix(m, 1, mxREAL);
  plhs[1] = mxCreateDoubleMatrix(m, 1, mxREAL);
  plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);*/

  /* Create pointers to output variables */
  /*u    = mxGetPr(plhs[0]);
  W    = mxGetPr(plhs[1]);
  iter = mxGetPr(plhs[2]);*/

  /* Copy initial values of u and W */
 /* SCOPY(&m,u0,&ione,u,&ione);*/ /* u = u0 */
 /* SCOPY(&m,W0,&ione,W,&ione);*/ /* W = W0 */

  /* Call wls_solve */
  /**iter = wls_solve(B,v,umin,umax,Wv,Wu,ud,gam,u,W,imax,k,m);*/

/*}*/

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * ~~~               Subroutines               ~~~
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void insert_Wu(float *Wu, float *A, int k, int m)
/* Insert Wu into A */ 
{
  int i, j, idx_A, idx_Wu;
  
  idx_Wu = 0;
  idx_A = k;

  for (j = 0 ; j < m ; j++) {
    for (i = 0 ; i < m ; i++) {
      A[idx_A] = Wu[idx_Wu];
      idx_A++;
      idx_Wu++;
    }
    idx_A += k;
  }
}  

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

int copy_columns(float *A, float *B, float *W, int m, int n)
/* Copy columns with zero entry in W from A (m x n) to B.
   Return number of columns copied = # zero elements in W */
{
  int i, idx_A = 0, idx_B = 0, ione = 1, cols = 0;
  
  for (i=0 ; i<n ; i++) {
    /* For each column of A */
    if (!W[i]) {
      /* Copy i:th column of A to B */
      SCOPY(&m,&A[idx_A],&ione,&B[idx_B],&ione);
      idx_B += m;
      cols++;
    }
    idx_A += m;
  }
  return cols;
}
 
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void add_vector_cond(float *a, float *b, float *W, int m)
/* Set b(not(W)) += a */
{
  int i, j=-1;

  for (i=0 ; i<m ; i++) {
    if (!W[i]) {
      j++;
      b[i] += a[j];
    }
  }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

int is_feasible(float *u, float *umin, float *umax, float *W, int m)
/* Check if umin <= u <= umax for free components of u */
{
  int i;
  
  for (i=0 ; i<m ; i++) {
    if (W[i]==0) /* Free variable */
      if ((u[i] < umin[i]) || (u[i] > umax[i]))
	return 0;
  }
  return 1;
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

int lagrange_mult(float *lambda_ns, float *W, int m)
/* Return the index of the most negative multiplier. Return -1 if all
   multipliers are positive. */
{
  float lambda, lambda_min;
  int i, i_min;
  
  /* Check if all lambda corresponding to active constraints are all
     positive. If not, find the most negative one. */
  lambda_min = -PREC;
  i_min = -1;
  for (i=0 ; i<m ; i++) {
    if (W[i]) { /* active constraint */
      lambda = W[i]*lambda_ns[i];
      if (lambda < lambda_min) {
	lambda_min = lambda;
	i_min = i;
      }
    }
  }

  /* Return index of constraint to be removed from working set. */
  return i_min;
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

float max_step(float *u, float *p_free, float *umin, float *umax,
	       float *W, int m, int *i_block, int *if_block)
/* Compute the maximum step length, a, such that u+a*p_free is
   feasible w.r.t. umin and umax. Store index of blocking constraint
   in i_block and if_block. */
{
  int i, j = -1;
  float a, a_min = 1;

  /* Determine step length a_min and the corresponding variable index. */
  for (i=0 ; i<m ; i++) {
    if (W[i]==0) { /* free variable */
      j++; /* index in p_free */

      /* Check the direction of the i:th perturbation component. */
      if (p_free[j]>0) {
	a = (umax[i]-u[i])/p_free[j]; /* step length to upper bound */
      } else {
	if(p_free[j]<0) {
	  a = (umin[i]-u[i])/p_free[j]; /* step length to lower bound */
	} else { /* p_free[j]=0, degenerate case */
	  a = 2; /* arbitrary number > 1 */
	}
      }

      /* Most restrictive step length so far? */
      if (a<a_min) {
	/* Yes, store step length and indeces. */
	a_min = a;
	*i_block = i;
	*if_block = j;
      }
    }
  }
  
  /* Return maximum feasible step length */
  return a_min;
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void print_vector(float *b, int n)
/* Display contents of an n vector. */
{
  int i;
  
  printf("[ ");
  for (i=0 ; i<n ; i++)
    printf("%.3g ",b[i]);
  printf("]\n");
}

void print_matrix(float *A, int m, int n)
/* Display contents of an m x n matrix. */
{
  int i, j;

  for (i=0 ; i<m ; i++) {
    if (i==0) printf("[ ");
    else printf("  ");

    for (j=0 ; j<n ; j++)
      printf("%8.3g ",A[i+j*m]);
    if (i<m-1) printf("\n");
    else printf("]\n");
  }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

