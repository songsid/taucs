/***************************************************************
 * MEX interfaces for TAUCS
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


/* Fill taucs_ccs_matrix from a matlab's mxArray */
void mex_fill_taucs_ccs_matrix(const mxArray *matlab_A, taucs_ccs_matrix *taucs_A)
{
  taucs_A->n = mxGetN(matlab_A);
  taucs_A->m = mxGetM(matlab_A);
  taucs_A->colptr = mxGetJc(matlab_A);
  taucs_A->rowind = mxGetIr(matlab_A);
  taucs_A->values.d = mxGetPr(matlab_A);
  taucs_A->flags = TAUCS_DOUBLE; 
}

/* ordering for qsort */
static int *compare_vals;
static int compar(const void *a, const void *b)
{
  int val_a = compare_vals[*(const int *)a];
  int val_b = compare_vals[*(const int *)b];

  return val_a - val_b;
}

/* Fill taucs_ccs_matrix from a matlab's mxArray */
mxArray *mex_convert_taucs_ccs_2_matlab(taucs_ccs_matrix *M)
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

  mxFree(tmp);

  return matlab_M;
}

/* free a css matrix */
void mex_free_ccs_matrix(taucs_ccs_matrix *a)
{
  if (a)
    {
      free(a->rowind);
      free(a->colptr);
      free(a->values.d/*taucs_values*/);
      free(a);
    }
}

/* get cpu time */
double mex_get_cpu_time()
{
  struct rusage a;
  
  getrusage(RUSAGE_SELF,&a);
  
  return (double) 
    (double) ((a.ru_utime).tv_sec +(a.ru_stime).tv_sec ) +
    (double) ((a.ru_utime).tv_usec+(a.ru_stime).tv_usec) * 1.0e-6;

}



