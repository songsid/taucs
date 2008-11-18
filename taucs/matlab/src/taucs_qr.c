/***************************************************************
 * MATLAB interface to TAUCS QR
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

/* simple mode */
void simple_mode(int nargout, mxArray *argout[], int nargin, const mxArray *argin[]) 
{
  const mxArray *matlab_A;
  taucs_ccs_matrix A, *R;
  taucs_multiqr_factor *F;
  int *column_order, *icol_order;
  double *pr, *B, max_kappa_R;
  int i, nrhs;

  if (nargin > 3)
    taucs_logfile("/tmp/taucs_qr.log"); 


  /* Getting parameters */
  matlab_A = argin[0];
  if (nargin > 1 && !mxIsEmpty(argin[1]) &&  mxIsDouble(argin[1])) 
  {
    argout[2] = mxDuplicateArray(argin[1]);
    B = mxGetPr(argout[2]);
    nrhs = mxGetN(argout[2]);
  } else
  {
    if (nargout > 2)
      argout[2] = mxCreateDoubleMatrix(0, 0, mxREAL);
    B = NULL;
    nrhs = 0;
  }
  if (nargin > 2)
  {
    max_kappa_R = mxGetScalar(argin[2]);
    if (max_kappa_R < 0)
      mexErrMsgTxt("max kappa R is smaller than 0");
  }
    else
    max_kappa_R = -1;

  /* Convert matrix */
  mex_fill_taucs_ccs_matrix(matlab_A, &A);

  /* Run factorizatrization */
  double tstart = mex_get_cpu_time();
  taucs_ccs_order(&A, &column_order, &icol_order, "colamd");
  if (max_kappa_R < 0)
    F = taucs_ccs_factor_qr(&A, column_order, 0, B, nrhs);
  else
    F = taucs_ccs_factor_pseudo_qr(&A, column_order, max_kappa_R, 0, B, nrhs);
  free(column_order);
  free(icol_order);
  if (F == NULL)
    mexErrMsgTxt("factorization failed");

  /* Copy to matlab allocated structure (just to make sure no problems with free */
  tstart = mex_get_cpu_time();

  R = taucs_multiqr_get_R(F, &column_order);
  argout[0] = mex_convert_taucs_ccs_2_matlab(R);

  argout[1] =  mxCreateDoubleMatrix(1, A.n, mxREAL);
  pr = mxGetPr(argout[1]);
  for(i = 0; i < A.n; i++)
    pr[i] = column_order[i] + 1;

  free(column_order);


  /* Free TAUCS format */
  mex_free_ccs_matrix(R);

  if (nargout > 3) 
  {
    if (max_kappa_R < 0)
      argout[3] =  mxCreateDoubleMatrix(0, 0, mxREAL);
    else 
    {
      int perb_number =  taucs_multiqr_get_perb_number(F);
      argout[3] =  mxCreateDoubleMatrix(1, perb_number, mxREAL);
      pr = mxGetPr(argout[3]);
      int *indices = malloc(perb_number * sizeof(int));
      taucs_multiqr_get_perb_indices(F, indices);
      for(i = 0; i < perb_number; i++)
	pr[i] = (int)indices[i] + (indices[i] > 0 ? 1 : -1);
    }
  }

  if (nargout > 4) 
  {
    if (max_kappa_R < 0)
      argout[4] = mxCreateDoubleScalar(0);
    else
      argout[4] = mxCreateDoubleScalar(taucs_multiqr_get_perb_value(F));
  }

  taucs_multiqr_factor_free(F);
}


void robust_factor(int nargout, mxArray *argout[], int nargin, const mxArray *argin[]) 
{
  const mxArray *matlab_A;
  taucs_ccs_matrix A;
  taucs_multiqr_factor *F;
  int *column_order, *icol_order;
  double max_kappa_R;
  
  
  matlab_A = argin[1];
  mex_fill_taucs_ccs_matrix(matlab_A, &A);
  
  if (nargin > 2)
  {
    max_kappa_R = mxGetScalar(argin[2]);
    if (max_kappa_R < 0)
      mexErrMsgTxt("max kappa R is smaller than 0");
  }
  else
    max_kappa_R = -1;

  taucs_ccs_order(&A, &column_order, &icol_order, "colamd");
  if (max_kappa_R < 0)
    F = taucs_ccs_factor_qr(&A, column_order, 1, NULL, 0);
  else
    F = taucs_ccs_factor_pseudo_qr(&A, column_order, max_kappa_R, 1, NULL, 0);
  if (F == NULL)
    mexErrMsgTxt("factorization failed");
  
  /* Keep factor pointer as long number */
  argout[0] = mxCreateNumericMatrix(1, 1, (sizeof(unsigned long) == 4 ? mxUINT32_CLASS : mxUINT64_CLASS), 0);
  *(taucs_multiqr_factor **)mxGetData(argout[0]) = F;


  /* Free TAUCS format */
  free(column_order);
  free(icol_order);   
}

void robust_get_R(int nargout, mxArray *argout[], int nargin, const mxArray *argin[]) 
{
  taucs_multiqr_factor *F;
  taucs_ccs_matrix *R;
  int *column_order;
  double *pr;
  int i;
  
  F = *(taucs_multiqr_factor **)mxGetData(argin[1]);
  
  R = taucs_multiqr_get_R(F, &column_order);
  argout[0] = mex_convert_taucs_ccs_2_matlab(R);
  argout[1] =  mxCreateDoubleMatrix(1, R->n, mxREAL);
  pr = mxGetPr(argout[1]);
  for(i = 0; i < R->n; i++)
    pr[i] = column_order[i] + 1;

  free(column_order);
  free_ccs_matrix(R);

  if (nargout > 2) 
  {
    int perb_number =  taucs_multiqr_get_perb_number(F);
     
    if (perb_number == 0)
      argout[2] =  mxCreateDoubleMatrix(0, 0, mxREAL);
    else 
    {
      argout[2] =  mxCreateDoubleMatrix(1, perb_number, mxREAL);
      pr = mxGetPr(argout[2]);
      int *indices = malloc(perb_number * sizeof(int));
      taucs_multiqr_get_perb_indices(F, indices);
      for(i = 0; i < perb_number; i++)
	pr[i] = (int)indices[i] + (indices[i] > 0 ? 1 : -1);
    }
  }

  if (nargout > 3) 
  {
    int perb_number =  taucs_multiqr_get_perb_number(F);
    if (perb_number == 0)
      argout[3] = mxCreateDoubleScalar(0);
    else
      argout[3] = mxCreateDoubleScalar(taucs_multiqr_get_perb_value(F));
  }
}


void robust_apply_Qt(int nargout, mxArray *argout[], int nargin, const mxArray *argin[]) 
{
  taucs_multiqr_factor *F;

  F = *(taucs_multiqr_factor **)mxGetData(argin[1]);

  const mxArray *b = argin[2];
  argout[0] = mxDuplicateArray(b);
  mxArray *x = argout[0];
  taucs_multiqr_apply_many_Qt(F, mxGetN(b), mxGetPr(x), mxGetM(b), mxGetPr(b), mxGetM(b));
}

void robust_apply_Q(int nargout, mxArray *argout[], int nargin, const mxArray *argin[]) 
{
  taucs_multiqr_factor *F;

  F = *(taucs_multiqr_factor **)mxGetData(argin[1]);

  const mxArray *b = argin[2];
  argout[0] = mxDuplicateArray(b);
  mxArray *x = argout[0];
  taucs_multiqr_apply_many_Q(F, mxGetN(b), mxGetPr(x), mxGetM(b), mxGetPr(b), mxGetM(b));
}

void robust_free(int nargout, mxArray *argout[], int nargin, const mxArray *argin[]) 
{
  taucs_multiqr_factor *F;

  F = *(taucs_multiqr_factor **)mxGetData(argin[1]);
  taucs_multiqr_factor_free(F);
}

/* robust mode */
void robust_mode(int nargout, mxArray *argout[], int nargin, const mxArray *argin[]) 
{
  if (nargin < 2)
    mexErrMsgTxt("only operation parameter");

  int operation = (int)mxGetScalar(argin[0]);
  switch(operation) 
  {
  case 0:
    robust_factor(nargout, argout, nargin, argin);
    break;
  case 1:
    robust_get_R(nargout, argout, nargin, argin);
    break;
  case 2:
    robust_apply_Qt(nargout, argout, nargin, argin);
    break;
  case 3:
    robust_apply_Q(nargout, argout, nargin, argin);
    break;
  case 4:
    robust_free(nargout, argout, nargin, argin);
    break;
  default:
    mexErrMsgTxt("unknown operation");
    break;
  }
}

/* mex (main) function */
void mexFunction(int nargout, mxArray *argout[], int nargin, const mxArray *argin[])
{
  /* There is some randomness in the perturbation algorithm (not really important).
     The following will make this function determenistic */
  srand(2233);

  if (nargin < 1)
    mexErrMsgTxt("zero parameters");


  if (mxIsSparse(argin[0]))
    simple_mode(nargout, argout, nargin, argin);
  else
    robust_mode(nargout, argout, nargin, argin);
}
