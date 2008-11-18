/***************************************************************
 * MEX for solving overdetermined LS equations using TAUCS.
 **************************************************************/

#include "taucs_mex.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <sys/timeb.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/resource.h>

/* mex (main) function */
void mexFunction(int nargout, mxArray *argout[], int nargin, const mxArray *argin[]) 
{
  const mxArray *matlab_A, *matlab_B;
  taucs_ccs_matrix A, *R;
  taucs_multiqr_factor *F;
  int *column_order, *icol_order;
  double *pr, *B, *X, *QTX;
  int i, j, nrhs, ld;
  double max_kappa_R;


  /* There is some randomness in the perturbation algorithm (not really important).
     The following will make this function determenistic */
  srand(2233);

  /* Make sure A is sparse */
  matlab_A = argin[0];
  if (!mxIsSparse(matlab_A))
  {
    mexErrMsgTxt("input matrix is not sparse\n");
    return;
  }
  mex_fill_taucs_ccs_matrix(matlab_A, &A);

  /* Copy B */
  matlab_B = argin[1];
  B = mxMalloc(mxGetM(matlab_B) * mxGetN(matlab_B) * sizeof(double));
  memcpy(B, mxGetPr(matlab_B), mxGetM(matlab_B) * mxGetN(matlab_B) * sizeof(double));
  nrhs = mxGetN(matlab_B);

  if (nargin > 2)
  {
    max_kappa_R = mxGetScalar(argin[2]);
    if (max_kappa_R < 0)
      mexErrMsgTxt("max kappa R is smaller than 0");
  }
    else
    max_kappa_R = -1;

  /* Run factorization */
  taucs_ccs_order(&A, &column_order, &icol_order, "colamd");
  if (max_kappa_R < 0)
    F = taucs_ccs_factor_qr(&A, column_order, 0, B, nrhs);
  else
    F = taucs_ccs_factor_pseudo_qr(&A, column_order, max_kappa_R, 0, B, nrhs);
  free(column_order);
  free(icol_order);
  if (F == NULL)
    mexErrMsgTxt("factorization failed");

  /* Solve */
  ld = mxGetN(matlab_A);
  QTX = mxMalloc(ld * nrhs * sizeof(double));
  taucs_multiqr_solve_many_R(F, nrhs, QTX, mxGetN(matlab_A), B, mxGetM(matlab_A));

  
  argout[0] = mxCreateDoubleMatrix(ld, nrhs, mxREAL);
  X = mxGetPr(argout[0]);

  column_order = taucs_multiqr_get_column_order(F);
  for(j = 0 ; j < mxGetN(matlab_A); j++)
    for (i = 0; i < nrhs; i++)
      X[column_order[j] + ld * i] = QTX[j + ld * i]; 
   

  /* Free */
  mxFree(B);
  free(column_order);
  taucs_multiqr_factor_free(F);
  mxFree(QTX);
}

