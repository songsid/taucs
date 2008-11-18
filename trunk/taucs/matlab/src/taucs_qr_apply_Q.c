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
  taucs_multiqr_factor *F;

  F = *(taucs_multiqr_factor **)mxGetData(argin[0]);

  argout[0] = mxDuplicateArray(argin[1]);
  taucs_multiqr_apply_many_Q(F, mxGetN(argin[1]), mxGetPr(argout[0]), mxGetM(argin[1]), mxGetPr(argin[1]), mxGetM(argin[1]));
}

