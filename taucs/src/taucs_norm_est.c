/*********************************************************/
/* TAUCS                                                 */
/* Author: Sivan Toledo                                  */
/*                                                       */
/* Functions for estimating norm of sparse matrices      */
/* Contributed by Haim Avron                             */
/*********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "taucs.h"

#ifndef TAUCS_CORE
#error "This is a TAUCS core file: you must define a primitive data type"
#endif


// TODO: Complex
#define RNDM ((double)rand()/(double)RAND_MAX);


#define NUM_POWER_ITS 4

#ifdef TAUCS_CORE_GENERAL

taucs_double taucs_gross_largest_sv_explicit(int m, int n,
					     taucs_ccs_matrix *A,
					     taucs_ccs_matrix *At, /* IF NULL WILL BE CALCULATED */
					     void **left_sv, void **right_sv)
{
#ifdef TAUCS_CONFIG_DREAL
  if (A->flags & TAUCS_DOUBLE)
    return taucs_dgross_largest_sv_explicit(m, n, A, At, left_sv, right_sv);
#endif

#ifdef TAUCS_CONFIG_SREAL
  if (A->flags & TAUCS_SINGLE)
    return taucs_sgross_largest_sv_explicit(m, n, A, At, left_sv, right_sv);
#endif

#ifdef TAUCS_CONFIG_DCOMPLEX
  if (A->flags & TAUCS_DCOMPLEX)
    return taucs_zgross_largest_sv_explicit(m, n, A, At, left_sv, right_sv);
#endif

#ifdef TAUCS_CONFIG_SCOMPLEX
  if (A->flags & TAUCS_SCOMPLEX)
    return taucs_cgross_largest_sv_explicit(m, n, A, At, left_sv, right_sv);
#endif

  return taucs_get_nan(); 
}

#endif

#ifndef TAUCS_CORE_GENERAL

taucs_double taucs_dtl(gross_largest_sv_explicit)(int m, int n,
						  taucs_ccs_matrix *A,
						  taucs_ccs_matrix *At, /* IF NULL WILL BE CALCULATED */
						  void **_p_left_sv, void **_p_right_sv)
{
  int i, iter;

  taucs_datatype **p_left_sv = (taucs_datatype **)_p_left_sv;
  taucs_datatype **p_right_sv = (taucs_datatype **)_p_right_sv;
  taucs_datatype *mv = taucs_malloc(m * sizeof(taucs_datatype));
  taucs_datatype *sv;

  sv = taucs_malloc(n * sizeof(taucs_datatype));;

  if (At == NULL)
    At = taucs_ccs_transpose(A);

  /* Init with random normalized vector */
  for(i = 0; i < n; i++)
    sv[i] = RNDM;
  taucs_double norm = taucs_vec_norm2(n, TAUCS_CORE_DATATYPE, sv);
  for(i = 0; i < n; i++)
    sv[i] /= norm;

  for(iter = 0; iter < NUM_POWER_ITS; iter++)
  { 
    taucs_ccs_times_vec(A, sv, mv);
    taucs_ccs_times_vec(At, mv, sv);
    norm = taucs_vec_norm2(n, TAUCS_CORE_DATATYPE, sv);
    for(i = 0; i < n; i++)
      sv[i] /= norm;
  }

  // Calc left from right svs
  if (p_left_sv != NULL) 
  {
    taucs_ccs_times_vec(A, sv, mv);
    taucs_double norm2 = taucs_vec_norm2(m, TAUCS_CORE_DATATYPE, mv); /* In infinite accuracy we can use norm */
    for(i = 0; i < m; i++)
      mv[i] /= norm2;
    *p_left_sv = mv;
  }
  else
    taucs_free(mv);

  if (p_right_sv != NULL)
    *p_right_sv = sv;
  else
    taucs_free(sv);

  return sqrt(norm);  
}

#endif


#ifdef TAUCS_CORE_GENERAL
taucs_double taucs_gross_largest_sv_implicit(int m, int n, int flags,
					     void matrix_func(void *, void *, void *), 
					     void trans_matrix_func(void *, void *, void *),
					     void *params, 
					     void **left_sv, void **right_sv)
{
#ifdef TAUCS_CONFIG_DREAL
  if (flags & TAUCS_DOUBLE)
    return taucs_dgross_largest_sv_implicit(m, n, matrix_func, trans_matrix_func, params, left_sv, right_sv);
#endif

#ifdef TAUCS_CONFIG_SREAL
  if (flags & TAUCS_SINGLE)
    return taucs_sgross_largest_sv_implicit(m, n, matrix_func, trans_matrix_func, params, left_sv, right_sv);
#endif

#ifdef TAUCS_CONFIG_DCOMPLEX
  if (flags & TAUCS_DCOMPLEX)
    return taucs_zgross_largest_sv_implicit(m, n, matrix_func, trans_matrix_func, params, left_sv, right_sv);
#endif

#ifdef TAUCS_CONFIG_SCOMPLEX
  if (flags & TAUCS_SCOMPLEX)
    return taucs_cgross_largest_sv_implicit(m, n, matrix_func, trans_matrix_func, params, left_sv, right_sv);
#endif

  return taucs_get_nan(); 
}

#endif

#ifndef TAUCS_CORE_GENERAL

taucs_double taucs_dtl(gross_largest_sv_implicit)(int m, int n,
						  void matrix_func(void *, void *, void *), 
						  void trans_matrix_func(void *, void *, void *),
						  void *params, 
						  void **_p_left_sv, void **_p_right_sv)
{
  int i, iter;

  taucs_datatype **p_left_sv = (taucs_datatype **)_p_left_sv;
  taucs_datatype **p_right_sv = (taucs_datatype **)_p_right_sv;
  taucs_datatype *mv = taucs_malloc(m * sizeof(taucs_datatype));
  taucs_datatype *sv;

  sv = taucs_malloc(n * sizeof(taucs_datatype));

  /* Init with random normalized vector */
  for(i = 0; i < n; i++)
    sv[i] = RNDM;
  taucs_double norm = taucs_vec_norm2(n, TAUCS_CORE_DATATYPE, sv);
  for(i = 0; i < n; i++)
    sv[i] /= norm;

  for(iter = 0; iter < NUM_POWER_ITS; iter++)
  { 
    matrix_func(sv, mv, params);
    trans_matrix_func(mv, sv, params);
    norm = taucs_vec_norm2(n, TAUCS_CORE_DATATYPE, sv);
    for(i = 0; i < n; i++)
      sv[i] /= norm;
  }

  // Calc left from right svs
  if (p_left_sv != NULL) 
  {
    matrix_func(sv, mv, params);
    taucs_double norm2 = taucs_vec_norm2(m, TAUCS_CORE_DATATYPE, mv); /* In infinite accuracy we can use norm */
    for(i = 0; i < m; i++)
      mv[i] /= norm2;
    *p_left_sv = mv;
  }
  else
    taucs_free(mv);

  if (p_right_sv != NULL)
    *p_right_sv = sv;
  else
    taucs_free(sv);

  return sqrt(norm);  
}

#endif

#ifdef TAUCS_CORE_GENERAL

taucs_double taucs_norm_1(taucs_ccs_matrix *A)
{
#ifdef TAUCS_CONFIG_DREAL
  if (A->flags & TAUCS_DOUBLE)
    return taucs_dnorm_1(A);
#endif

#ifdef TAUCS_CONFIG_SREAL
  if (A->flags & TAUCS_SINGLE)
    return taucs_snorm_1(A);
#endif

#ifdef TAUCS_CONFIG_DCOMPLEX
  if (A->flags & TAUCS_DCOMPLEX)
    return taucs_znorm_1(A);
#endif

#ifdef TAUCS_CONFIG_SCOMPLEX
  if (A->flags & TAUCS_SCOMPLEX)
    return taucs_cnorm_1(A);
#endif

  return taucs_get_nan(); 
}

#endif


#ifndef TAUCS_CORE_GENERAL
taucs_double taucs_dtl(norm_1)(taucs_ccs_matrix *A)
{
  taucs_double norm1 = -1;
  for(int i = 0; i < A->n; i++)
  {
    taucs_double this = 0.0;
    for(int j = A->colptr[i]; j < A->colptr[i+1]; j++)
      this += taucs_abs(A->taucs_values[j]);
    norm1 = max(this, norm1);
  }

  return norm1;
	 
}

#endif
