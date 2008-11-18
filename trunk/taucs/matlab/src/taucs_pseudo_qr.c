/***************************************************************
 * Matlab Interface for BoomerAMG 
 * Written by Haim Avron, June 2006.
 **************************************************************/

#include "mex.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <sys/timeb.h>
#include <sys/time.h>
#include <sys/resource.h>


#include "taucs.h"

extern taucs_multiqr_factor* taucs_ccs_factor_pseudo_qr(taucs_ccs_matrix* A, int *column_order,
							double max_kappa_R,
							int keep_q, void *B, int nrhs);

extern taucs_ccs_matrix *taucs_multiqr_get_R(taucs_multiqr_factor *F, int **column_order);

extern taucs_ccs_matrix *taucs_get_column_order(taucs_multiqr_factor *F, int **column_order);

double taucs_multiqr_get_perb_value(taucs_multiqr_factor *F);

int taucs_multiqr_get_perb_number(taucs_multiqr_factor *F);

void taucs_multiqr_get_perb_indices(taucs_multiqr_factor *F, int *indices);

void taucs_multiqr_factor_free(taucs_multiqr_factor* F);

extern void taucs_ccs_order(taucs_ccs_matrix* m, 
			    int** perm, int** invperm,
			    char* which);

extern void taucs_profile_report();

static void free_ccs_matrix(taucs_ccs_matrix *a)
{
  if (a)
    {
      free(a->rowind);
      free(a->colptr);
      free(a->values.d/*taucs_values*/);
      free(a);
    }
}

double get_cpu_time()
{
  struct rusage a;
  
  getrusage(RUSAGE_SELF,&a);
  
  return (double) 
    (double) ((a.ru_utime).tv_sec +(a.ru_stime).tv_sec ) +
    (double) ((a.ru_utime).tv_usec+(a.ru_stime).tv_usec) * 1.0e-6;

}

typedef struct 
{
  double part_goal;
  int seed;
  bool stretch;
  bool profile;
  bool amwb_force;
} params_struct;

/* ordering for qsort */
int *compare_vals;
int compar(const void *a, const void *b)
{
  int val_a = compare_vals[*(const int *)a];
  int val_b = compare_vals[*(const int *)b];

  return val_a - val_b;
}


/* Fill taucs_ccs_matrix from a matlab's mxArray */
void fill_taucs_ccs_matrix(const mxArray *matlab_A, taucs_ccs_matrix *taucs_A)
{
  taucs_A->n = mxGetN(matlab_A);
  taucs_A->m = mxGetM(matlab_A);
  taucs_A->colptr = mxGetJc(matlab_A);
  taucs_A->rowind = mxGetIr(matlab_A);
  taucs_A->values.d = mxGetPr(matlab_A);
  taucs_A->flags = TAUCS_DOUBLE; 
}


/* Fill taucs_ccs_matrix from a matlab's mxArray */
mxArray *convert_taucs_ccs_2_matlab(taucs_ccs_matrix *M)
{
  int nnz = M->colptr[M->n];
  mxArray *matlab_M;
  int *tmp = (int *)mxMalloc(M->n * sizeof(int));
  int i, j;

  matlab_M = mxCreateSparse(M->m, M->n, nnz, mxREAL);
  
  memcpy(mxGetJc(matlab_M), M->colptr, (M->n + 1) * sizeof(int));
  
  /* The following is complicated because the row indications in each col must be sorted */
  for(i = 0; i < M->n; i++)
  {
    /* Find order */
    compare_vals = M->rowind + M->colptr[i];
    for(j = 0; j < M->colptr[i+1] - M->colptr[i]; j++)
      tmp[j] = j;
    qsort(tmp, M->colptr[i+1] - M->colptr[i], sizeof(int), compar);

    /* copy vals */
    for(j = 0; j < M->colptr[i+1] - M->colptr[i]; j++)
    {
      *(mxGetIr(matlab_M) + M->colptr[i] + j) = M->rowind[M->colptr[i] + tmp[j]];
      *(mxGetPr(matlab_M) + M->colptr[i] + j) = M->values.d[M->colptr[i] + tmp[j]];
    }
  }

  /* 
  memcpy(mxGetIr(matlab_M), M->rowind, nnz * sizeof(int));
  memcpy(mxGetPr(matlab_M), M->values.d, nnz * sizeof(double)); */

  mxFree(tmp);

  return matlab_M;
}

/* Actual mex function */
void mexFunction(int nargout, mxArray *argout[], int nargin, const mxArray *argin[]) 
{
  const mxArray *matlab_A;
  taucs_ccs_matrix A, *R;
  taucs_multiqr_factor *F;
  int *column_order, *icol_order;
  double *pr, *B, max_kappa_R;
  int i, nrhs;

  if (nargin > 2)
    taucs_logfile("/tmp/taucs_qr.log");

  /* Make sure A is sparse */
  matlab_A = argin[0];
  if (!mxIsSparse(matlab_A))
  {
    mexPrintf("input matrix A must be sparse\n");
    return;
  }

  max_kappa_R = mxGetScalar(argin[1]);

  /* Get B */
  if (nargin > 2 && !mxIsEmpty(argin[2]) && mxIsDouble(argin[2])) 
  {
    argout[2] = mxDuplicateArray(argin[2]);
    B = mxGetPr(argout[2]);
    nrhs = mxGetN(argout[2]);
  } else
  {
    if (nargout > 2)
      argout[2] = mxCreateDoubleMatrix(0, 0, mxREAL);
    B = NULL;
    nrhs = 0;
  }

  /* now we can convert to matrix */
  fill_taucs_ccs_matrix(matlab_A, &A);

  /* Run factorizatrization */
  double tstart = get_cpu_time();
  taucs_ccs_order(&A, &column_order, &icol_order, "colamd");
  F = taucs_ccs_factor_pseudo_qr(&A, column_order, max_kappa_R, 0, B, nrhs);
  free(column_order);
  free(icol_order);
  if (F == NULL)
  {
    mexPrintf("Factorization failed\n");fflush(stdout);
    return;
  }
  if (nargin > 3)
    mexPrintf("Actual factor time is %.2es\n", get_cpu_time() - tstart);

  /* Copy to matlab allocated structure (just to make sure no problems with free */
  tstart = get_cpu_time();

  R = taucs_multiqr_get_R(F, &column_order);
  argout[0] = convert_taucs_ccs_2_matlab(R);

  argout[1] =  mxCreateDoubleMatrix(1, A.n, mxREAL);
  pr = mxGetPr(argout[1]);
  for(i = 0; i < A.n; i++)
    pr[i] = column_order[i] + 1;

  free(column_order);

  if (nargin > 3)
    mexPrintf("Getting result time is %.2es\n", get_cpu_time() - tstart);

  /* Free TAUCS format */
  free_ccs_matrix(R);

  if (nargout > 3) 
  {
    int perb_number =  taucs_multiqr_get_perb_number(F);
    argout[3] =  mxCreateDoubleMatrix(1, perb_number, mxREAL);
    pr = mxGetPr(argout[3]);
    int *indices = malloc(perb_number * sizeof(int));
    taucs_multiqr_get_perb_indices(F, indices);
    for(i = 0; i < perb_number; i++)
      pr[i] = (int)indices[i] + (indices[i] > 0 ? 1 : -1);
  }

  if (nargout > 4)
    argout[4] = mxCreateDoubleScalar(taucs_multiqr_get_perb_value(F));

  taucs_multiqr_factor_free(F);

  if (nargin > 3) 
  {
    taucs_profile_report();
    mexPrintf("Look at more profiling at /tmp/taucs_qr.log\n");
  }

}

