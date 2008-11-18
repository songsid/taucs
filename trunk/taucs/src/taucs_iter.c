/*********************************************************/
/* TAUCS                                                 */
/* Author: Sivan Toledo                                  */
/*********************************************************/

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include "taucs.h"

#ifdef TAUCS_CONFIG_PFUNC
#include "pfunc.h"
#endif

#ifdef TAUCS_CORE_DOUBLE

/*********************************************************/
/* utilities                                             */
/*********************************************************/
/*extern int _isnan(double);*/

static int element_size(int flags)
{
  if (flags & TAUCS_SINGLE)   return sizeof(taucs_single);
  if (flags & TAUCS_DOUBLE)   return sizeof(taucs_double);
  if (flags & TAUCS_SCOMPLEX) return sizeof(taucs_scomplex);
  if (flags & TAUCS_DCOMPLEX) return sizeof(taucs_dcomplex);
  if (flags & TAUCS_INT)      return sizeof(int);
  assert(0);
  return -1;
}

static double dotprod(int n, double* v, double* u)
{
  double x;
  int i;

  for (i=0, x=0.0; i<n; i++) x += v[i]*u[i];

  return x;
}

#ifdef TAUCS_CONFIG_PFUNC
static double parallel_dotprod(double *v, double *w, int sv, int ev, double *S, int tid, pfunc_handle_t handle, int nproc)
{
  //  static int k[2] = {0, 0};

  int i = 0;
  *(S + tid) = 0.0;
  for (i = sv; i < ev; i++)
    *(S + tid) += v[i] * w[i];

  //double hit = taucs_wtime();
  pfunc_barrier();
  //double exit = taucs_wtime();

  //taucs_printf("(tid: %d) K = %d, Hit = %f, Exit = %f.\n", tid, k[tid], hit, exit);
  //k[tid]++;

  double r = 0.0;
  for(i = 0; i < nproc; i++)
    r += S[i];



  return r;
}

#endif

static double twonorm(int n, double* v)
{
  /*
  double norm;
  int i;

  for (i=0, norm=0.0; i<n; i++) norm += v[i]*v[i];

  norm = sqrt(norm);
  return norm;
  */

  double ssq, scale, absvi;/*norm omer*/
  int i;

  if (n==1) return fabs(v[0]);

  scale = 0.0;
  ssq   = 1.0;

  for (i=0; i<n; i++) {
    if ( v[i] != 0 ) {
      absvi = fabs(v[i]);
      if (scale < absvi) {
	ssq   = 1.0 + ssq * (scale/absvi)*(scale/absvi);
	scale = absvi;
      } else
	ssq   = ssq + (absvi/scale)*(absvi/scale);
    }
  }
  return scale * sqrt( ssq );
}

#ifdef TAUCS_CONFIG_PFUNC

/* Map and reduce implementation of twonorm */
static void twonorm_map(double* v, int sv, int ev, double *scale, double *ssq)
{
  double absvi;
  int i;

  if (ev == sv + 1) {
    *scale = fabs(v[sv]);
    *ssq = 1;
  }

  double scale_ = 0.0;
  double ssq_ = 1.0;

  for (i = sv; i < ev; i++) {
    if ( v[i] != 0 ) {
      absvi = fabs(v[i]);
      if (scale_ < absvi) {
	ssq_   = 1.0 + ssq_ * (scale_/absvi)*(scale_/absvi);
	scale_ = absvi;
      } else
	ssq_ = ssq_ + (absvi/scale_)*(absvi/scale_);
    }
  }
  
  *scale = scale_;
  *ssq = ssq_;
}

static double twonorm_reduce(double *S1, double *S2, int nproc)
{
  double scale = S1[0];
  double ssq = S2[0];

  for(int i = 1; i < nproc; i++) {
    if (S1[i] < scale)
      ssq += (S1[i] / scale) * (S1[i] / scale) * S2[i];
    else {
      ssq = S2[i] + ssq * (scale / S1[i]) * (scale / S1[i]);
      scale = S1[i];
    }
  }
  
  return scale * sqrt(ssq);
}

static double parallel_twonorm(double *v, int sv, int ev, double *S1, double *S2, int tid, pfunc_handle_t handle, int nproc)
{
  //static int k[2] = {0, 0};

  twonorm_map(v, sv, ev, S1 + tid, S2 + tid);
  //double hit = taucs_wtime();
  pfunc_barrier();
  //double exit = taucs_wtime();
  //taucs_printf("(tid: %d) K = %d, Hit = %f, Exit = %f.\n", tid, k[tid], hit, exit);
  //k[tid]++;

  return twonorm_reduce(S1, S2, nproc);
}

#endif 
 
/*********************************************************/
/* conjugate gradients                                   */
/*********************************************************/

#ifdef TAUCS_CONFIG_PFUNC
static void taucs_pfunc_conjugate_gradients(void *args);

/* from metis.h */
typedef int idxtype; 
void METIS_PartGraphRecursive(int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, int *, idxtype *);

int taucs_parallel_conjugate_gradients(taucs_ccs_matrix* A,
				       int               (*precond_fn)(void*,void* x,void* b),
				       void*             precond_args,
				       void*             vX,
				       void*             vB,
				       int               itermax,
				       double            convergetol,
				       int               nproc
				       )
{
  int i, j, ip;
  char **args = taucs_malloc(nproc * sizeof(void *));
  pfunc_handle_t *pfunc_handles = taucs_malloc(nproc * sizeof(pfunc_handle_t));

  /* Partition the matrix */
  if (!(A->flags & TAUCS_SYMMETRIC) && !(A->flags & TAUCS_HERMITIAN)) {
    taucs_printf("taucs_parallel_conjugate_gradients: Matrix must be symmetric \n");
    return -1;
  }
  /* this routine may actually work on UPPER as well */
  if (!(A->flags & TAUCS_LOWER)) {
    taucs_printf("taucs_parallel_conjugate_gradients: Lower part of the matrix must be represented \n");
    return -1;
  }

  taucs_printf("taucs_parallel_conjugate_gradients: Calling METIS to partition the matrix.\n");

  int n   = A->n;
  int nnz = (A->colptr)[n];
  
  int *partition    = (int*) taucs_malloc(n * sizeof(int));

  idxtype *xadj = (int*) taucs_malloc((n+1) * sizeof(int));
  idxtype *adj  = (int*) taucs_malloc(2* nnz * sizeof(int));

  if (!partition || !xadj || !adj) {
    taucs_free(partition);
    taucs_free(xadj);
    taucs_free(adj);
    return -1;
  }

  int *ptr = partition;
  int *len = partition;

  for (i=0; i<n; i++) len[i] = 0;

  for (j=0; j<n; j++) {
    for (ip = (A->colptr)[j]; ip < (A->colptr)[j+1]; ip++) {
      i = (A->rowind)[ip];
      if (i != j) {
	len[i] ++;
	len[j] ++;
      }
    }
  }

  xadj[0] = 0;
  for (i=1; i<=n; i++) xadj[i] = xadj[i-1] + len[i-1];
  
  for (i=0; i<n; i++) ptr[i] = xadj[i];

  for (j=0; j<n; j++) {
    for (ip = (A->colptr)[j]; ip < (A->colptr)[j+1]; ip++) {
      i = (A->rowind)[ip];
      if (i != j) {
	adj[ ptr[i] ] = j;
	adj[ ptr[j] ] = i;
	ptr[i] ++;
	ptr[j] ++;
      }
    }
  }

  
  int wgtflag = 0;
  int numflag = 0;
  int nparts = nproc;
  int options_flag = 0;
  int edgecut;
  METIS_PartGraphRecursive(&n,
			   xadj,adj,
			   NULL, NULL,
			   &wgtflag, &numflag,
			   &nparts, &options_flag,
			   &edgecut, partition);
  


  taucs_printf("taucs_parallel_conjugate_gradients: METIS finished partitioning the matrix.\n");

  taucs_printf("taucs_parallel_conjugate_gradients: Reordering the matrix for parallelism.\n");

  /* Reorder the matrix by groups */
  int *group_sizes = (int *)taucs_malloc(sizeof(int) * nproc);
  int *boundary = (int *)taucs_malloc(sizeof(int) * n);
  int *boundary_sizes = (int *)taucs_malloc(sizeof(int) * nproc);
  memset(boundary, 0, n * sizeof(int));
  memset(boundary_sizes, 0, nproc * sizeof(int));
  memset(group_sizes, 0, nproc * sizeof(int));
  for(j = 0; j < n; j++) {
    group_sizes[partition[j]]++;
    for (int ip = (A->colptr)[j]; ip < (A->colptr[j+1]); ip++) {
      int k  = (A->rowind)[ip];
      if (partition[k] != partition[j]) {
	boundary[j] = 1;
	boundary[k] = 1;
      }
    }
  }

  for(j = 0; j < n; j++)
    if (boundary[j])
      boundary_sizes[partition[j]]++;
  int *group_starts = (int *)taucs_malloc(sizeof(int) * (nproc + 1));
  int *boundary_starts = (int *)taucs_malloc(sizeof(int) * nproc);
  group_starts[0] = 0;
  for(i = 1; i <= nproc; i++) {
    group_starts[i] = group_starts[i - 1] + group_sizes[i - 1];
    boundary_starts[i - 1] = group_starts[i] - boundary_sizes[i - 1];
  }
  int *perm = (int *)taucs_malloc(sizeof(int) * n);
  int *invperm = (int *)taucs_malloc(sizeof(int) * n);
  for(i = 0; i < n; i++) {
    if (!boundary[i]) {
      perm[group_starts[partition[i]]] = i;
      invperm[i] = group_starts[partition[i]];
      group_starts[partition[i]]++;
    } else {
      perm[boundary_starts[partition[i]]] = i;
      invperm[i] = boundary_starts[partition[i]];
      boundary_starts[partition[i]]++;
    }
  }
  group_starts[0] = 0;
  for(i = 1; i <= nproc; i++) {
    group_starts[i] = group_starts[i - 1] + group_sizes[i - 1];
    boundary_starts[i - 1] = group_starts[i] - boundary_sizes[i - 1];
  }
  taucs_free(partition);

  taucs_ccs_matrix *PAPT = taucs_ccs_permute_symmetrically(A, perm,invperm);
  void *vPB = vX;
  taucs_vec_permute(n, A->flags, vB, vPB, perm);
  void *vPX = (void*)taucs_malloc(element_size(A->flags) * n);
  memset(vPX, 0, element_size(A->flags) * n);

  taucs_printf("taucs_parallel_conjugate_gradients: Finished reordering.\n");
  
  /* Allocate shared spaces */
  double *P = (double*)taucs_malloc(n * sizeof(double));
  double *R = (double*)taucs_malloc(n * sizeof(double));
  double *Q = (double*)taucs_malloc(n * sizeof(double));
  double *Z = (double*)taucs_malloc(n * sizeof(double));
  double *S1 = (double*)taucs_malloc(3 * nproc * sizeof(double));
  double *S2 = (double*)taucs_malloc(nproc * sizeof(double));
  double *S3 = (double*)taucs_malloc(n * nproc * sizeof(double));
  memset(S3, 0, n * nproc * sizeof(double));

  /* Create group */
  pfunc_group_t pfunc_group;
  pfunc_group_init(&pfunc_group); 
  pfunc_group_set(&pfunc_group, GROUP_ID, 1234); 
  pfunc_group_set(&pfunc_group, GROUP_SIZE, nproc); 
  pfunc_group_set(&pfunc_group, GROUP_BARRIER_TYPE, BARRIER_SPIN);


  for(i = 0; i < nproc; i++) {
    pfunc_handle_init(&pfunc_handles[i]);

    pfunc_pack(&args[i], 
	       "void*, void *, void*, void*, void*, int, double, int, int, int, double*, double*, double*, double*, double*, double*, double*, int, void *, int", 
	       PAPT, precond_fn, precond_args, 
	       vPX, vPB, itermax, convergetol, group_starts[i], boundary_starts[i], group_starts[i + 1],
	       P, R, Q, Z, S1, S2, S3, i, pfunc_handles[i], nproc);

    pfunc_run (&pfunc_handles[i], 
	       PFUNC_ATTR_DEFAULT,
	       pfunc_group,
	       taucs_pfunc_conjugate_gradients,
	       args[i]);
  }
  pfunc_wait_all(pfunc_handles, nproc);

  taucs_vec_permute(n, A->flags, vPX, vPX, invperm);  
  
  for(i = 0; i < nproc; i++)
    pfunc_handle_clear(pfunc_handles[i]);

  taucs_free(pfunc_handles);
  taucs_free(args);
  taucs_free(group_starts);
  taucs_free(perm);
  taucs_free(invperm);
  taucs_free(vPX);

  taucs_free(P);
  taucs_free(R);
  taucs_free(Q);
  taucs_free(Z);
  taucs_free(S1);
  taucs_free(S2);
  
  return 0;
}

static void taucs_pfunc_conjugate_gradients(void *args)
{
  /* Parameters */
  taucs_ccs_matrix* A;
  int (*precond_fn)(void*,void* x,void* b);
  void* precond_args;
  void* vX;
  void* vB;
  int   itermax;
  double convergetol;
  int sv, mv, ev;
  int tid, nproc;
  int i;
  double *P, *R, *Q, *Z;
  double *S1, *S2, *S3;
  pfunc_handle_t handle; 

  pfunc_unpack(args, 
	       "void*, void *, void*, void*, void*, int, double, int, int, int, double*, double*, double*, double*, double*, double*, double*, int, void *, int", 
	       (void *)&A, (void *)&precond_fn, &precond_args, 
	       &vX, &vB, &itermax, &convergetol, &sv, &mv, &ev, 
	       &P, &R, &Q, &Z, &S1, &S2, &S3, &tid, &handle, &nproc);


  taucs_printf("taucs_pfunc_conjugate_gradients: (tid %d) Starting to work. sv = %d, mv = %d, ev = %d. \n", tid, sv, mv, ev);

  double* X = (double*) vX;
  double* B = (double*) vB;

  double Alpha, Beta, Rho, Init_norm, ratio, Res_norm, Rtmp ;
  double Rho0 = 0.0; /* warning */
  double Tiny = 0.1e-28;
  int    Iter;

  int n = A->n;

  double ct;
  if (tid == 0) 
    ct = -taucs_wtime();

  taucs_ccs_parallel_times_vec(A, X, R, sv, mv, ev, tid, S3, nproc, handle);

  for (i = sv; i < ev; i++)
    R[i] = B[i] - R[i]; 

  Res_norm = Init_norm = parallel_twonorm(R, sv, ev, S1, S2, tid, handle, nproc);

  if (tid == 0)
    printf("(tid: %d) two norm of initial residual %.2e\n", tid, Init_norm);

  if ( Init_norm == 0.0 ) 
    Init_norm = 1.0;

  ratio = 1.0;
 
  Iter = 0;

  while ( ratio > convergetol && Iter <= itermax ) {
    Iter++;
    
    if (precond_fn)
      (*precond_fn)(precond_args,Z,R);
    else
      for (i = sv; i < ev; i++)
	Z[i] = R[i];

    Rho = parallel_dotprod(R, Z, sv, ev, S1 + 2 * nproc, tid, handle, nproc);

    if ( Iter == 1 ) {
      for (i = sv; i < ev; i++) 
	P[i] = Z[i];
    } else {
      Beta = Rho /(Rho0 + Tiny);
      for (i = sv; i < ev; i++) 
	P[i] = Z[i] + Beta * P[i];
    };

    taucs_ccs_parallel_times_vec(A ,P, Q, sv, mv, ev, tid, S3, nproc, handle); /* Q = A*P */

    Rtmp = parallel_dotprod(P, Q, sv, ev, S1, tid, handle, nproc);

    Alpha = Rho/(Rtmp+Tiny);

    for (i = sv; i < ev; i++) {
      X[i] = X[i] + Alpha * P[i];
      R[i] = R[i] - Alpha * Q[i];
    }

    Rho0  = Rho;

    Res_norm = parallel_twonorm(R, sv, ev, S1 + nproc, S2, tid, handle, nproc);
    
    ratio = Res_norm/Init_norm;
    if (tid == 0 && Iter % 25 == 0) 
      taucs_printf("cg: n=%d at iteration %d the convergence ratio is %.2e, Rnorm %.2e\n", 
		   A->n,Iter, ratio,Res_norm) ;
  }
  pfunc_barrier();

  if (tid == 0) {
    ct += taucs_wtime();
    taucs_printf("Iteration time = %.2es\n", ct);
  }

  if (tid == 0 && Iter > 0) {
    taucs_printf("cg: n=%d iterations = %d Reduction in residual norm %.2e, Rnorm %.2e\n", 
		 A->n,Iter,ratio,Res_norm) ;
    taucs_ccs_times_vec(A,X,R);
    for (i=0; i<n; i++) R[i] = B[i] - R[i];
    taucs_printf("cg: true relative residual norm %.2e\n",twonorm(n,R) / Init_norm);
  }

  // THE FOLLOWING LINES ARE DEBUG
  pfunc_barrier();
  fflush(stdout);
  taucs_printf("taucs_pfunc_conjugate_gradients: (tid %d) Finished working.\n", tid);
}

#endif

int 
taucs_conjugate_gradients(taucs_ccs_matrix* A,
			  int               (*precond_fn)(void*,void* x,void* b),
			  void*             precond_args,
			  void*             vX,
			  void*             vB,
			  int               itermax,
			  double            convergetol
			  )
{
  double* X = (double*) vX;
  double* B = (double*) vB;
  double *P, *R, *Q, *Z ;
  double Alpha, Beta, Rho, Init_norm, ratio, Res_norm, Rtmp ;
  double Rho0 = 0.0; /* warning */
  /*double t1, t2,  cpus[9] ; omer*/
  /*
  double one[2] = {1.0, 0.0};
  double zero[2] = {0.0, 0.0} ;
  */
  /*  double Tiny = 0.0;*/
  double Tiny = 0.1e-28;
  int    Iter;
  /*int    stats[6] ; omer*/
  int    i,n;

#define RESVEC_NO
#ifdef RESVEC
  FILE* f;
  double* resvec = (double*) taucs_malloc((itermax+2) * sizeof(double));
  assert(resvec);
  for (i=0; i<=itermax; i++) {
    /*double inf = 1.0/0.0; omer*/
    double nan = taucs_get_nan()/*inf - inf; omer*/
    assert(taucs_isnan(nan));
    resvec[i] = nan;
  }
#endif

  n = A->n;
 
  P = (double*) taucs_malloc(n * sizeof(double));
  R = (double*) taucs_malloc(n * sizeof(double));
  Q = (double*) taucs_malloc(n * sizeof(double));
  Z = (double*) taucs_malloc(n * sizeof(double));

#define TAUCS_REMOVE_CONST_NO
#ifdef TAUCS_REMOVE_CONST
    {
      double s;
      for (i=0, s=0.0; i<n; i++) s += B[i];
      for (i=0, s=0.0; i<n; i++) B[i] -= s;
    }
#endif

  /*
  for (i=0; i<n; i++) X[i] = 0;
  for (i=0; i<n; i++) R[i] = B[i];
  */

  double ct = -taucs_wtime();

  taucs_ccs_times_vec(A,X,R);

  for (i=0; i<n; i++) R[i] = B[i] - R[i];

  Res_norm = Init_norm = twonorm(n,R);
  printf("two norm of initial residual %.2e\n",Init_norm);
  if ( Init_norm == 0.0 ) Init_norm = 1.0;
  ratio = 1.0;
 
  Iter = 0;
 
#ifdef RESVEC
  resvec[Iter] = Res_norm;
#endif

  while ( ratio > convergetol && Iter <= itermax ) {
    Iter++;
    
    if (precond_fn)
      (*precond_fn)(precond_args,Z,R);
    else
      for (i=0; i<n; i++) Z[i] = R[i];

    for (i=0,Rho=0.0; i<n; i++) Rho += R[i] * Z[i];


    if ( Iter == 1 ) {
      for (i=0; i<n; i++) P[i] = Z[i];
    } else {
      Beta = Rho /(Rho0 + Tiny);
      for (i=0; i<n; i++) P[i] = Z[i] + Beta * P[i];
    };

    taucs_ccs_times_vec(A,P,Q); /* Q = A*P */

    for (i=0,Rtmp=0.0; i<n; i++) Rtmp += P[i] * Q[i];
  
    Alpha = Rho/(Rtmp+Tiny);

    for (i=0; i<n; i++) X[i] = X[i] + Alpha * P[i];

    for (i=0; i<n; i++) R[i] = R[i] - Alpha * Q[i];

#ifdef TAUCS_REMOVE_CONST
    {
      double s;
      for (i=0, s=0.0; i<n; i++) s += R[i];
      for (i=0, s=0.0; i<n; i++) R[i] -= s;
    }
#endif


    Rho0  = Rho;

    Res_norm = twonorm(n,R);

#if 0
    taucs_ccs_times_vec(A,X,R);
    for (i=0; i<n; i++) R[i] -= B[i];
    Res_norm = twonorm(n,R);
#endif

#ifdef RESVEC
  resvec[Iter] = Res_norm;
#endif

    ratio = Res_norm/Init_norm;
    if (Iter % 25 == 0) 
      taucs_printf("cg: n=%d at iteration %d the convergence ratio is %.2e, Rnorm %.2e\n", 
		   A->n,Iter, ratio,Res_norm) ;
  }

  ct += taucs_wtime();
  taucs_printf("Iteration time = %.2es\n", ct);

  if (Iter > 0) {
    taucs_printf("cg: n=%d iterations = %d Reduction in residual norm %.2e, Rnorm %.2e\n", 
		 A->n,Iter,ratio,Res_norm) ;
    taucs_ccs_times_vec(A,X,R);
    for (i=0; i<n; i++) R[i] = B[i] - R[i];
    taucs_printf("cg: true residual norm %.2e\n",twonorm(n,R));
  }

  taucs_free(P) ;
  taucs_free(R) ;
  taucs_free(Q) ;
  taucs_free(Z) ;
 
#ifdef RESVEC
  f=fopen("resvec","a");
  assert(f);
  for (i=0; i<=itermax && !taucs_isnan(resvec[i]); i++) {
    fprintf(f,"%.3e\n",resvec[i]);
  }
  fclose(f);
  taucs_free(resvec);
#endif

  return 0; 
}                                                                             

/*********************************************************/
/* minres                                                */
/*********************************************************/

int 
taucs_minres(taucs_ccs_matrix*  A,
	     int                (*precond_fn)(void*,void* x,void* b),
	     void*              precond_args,
	     void*              vX,
	     void*              vB,
	     int                itermax,
	     double             convergetol)
{
  double* X = (double*) vX;
  double* B = (double*) vB;

  double *Xcg, *R, *V, *VV, *Vold, *Volder, *M, *Mold, *Molder;
  double tolb, normr, alpha, beta, beta1, betaold;
  double gamma, gammabar, delta, deltabar, epsilon;
  double cs,sn,snprod, numer, denom;
  int    Iter;
  int    i,n;

  n = A->n;
 
  R      = (double*) taucs_malloc(n * sizeof(double));
  Xcg    = (double*) taucs_malloc(n * sizeof(double));
  VV     = (double*) taucs_malloc(n * sizeof(double));
  V      = (double*) taucs_malloc(n * sizeof(double));
  Vold   = (double*) taucs_malloc(n * sizeof(double));
  Volder = (double*) taucs_malloc(n * sizeof(double));
  M      = (double*) taucs_malloc(n * sizeof(double));
  Mold   = (double*) taucs_malloc(n * sizeof(double));
  Molder = (double*) taucs_malloc(n * sizeof(double));

  tolb = convergetol * twonorm(n,B);
  taucs_printf("minres: residual convergence tolerance %.1e\n",tolb);
 
  for (i=0; i<n; i++) X[i] = 0;    /* x = 0 */
  for (i=0; i<n; i++) R[i] = B[i]; /* r = b-A*x */

  normr = twonorm(n,R);
  if ( normr == 0.0 ) {
    taucs_printf("minres: initial residual == 0\n");
    return -1;
  }

  for (i=0; i<n; i++) V[i]    = R[i];    /* v = r */
  for (i=0; i<n; i++) Vold[i] = R[i];    /* vold = r */
  
  if (precond_fn)
    (*precond_fn)(precond_args,V,Vold);
  else
    for (i=0; i<n; i++) V[i] = Vold[i];
  
  beta1 = dotprod(n,Vold,V);
  if (beta1 < 0.0) {
    taucs_printf("minres: error (1)\n");
    return -1;
  }
  beta1 = sqrt(beta1);

  { int flag = 0;
    for (i=0; i<n; i++) {
      if (taucs_isnan(V[i]) && flag < 10) 
	taucs_printf("minres: V has nan's in position %d\n",i);
      flag++;
    }
  }


  snprod = beta1;
  taucs_printf(">>> %e %e %e\n",beta1,snprod,normr);


  for (i=0; i<n; i++) VV[i] = V[i] / beta1;
  
  taucs_ccs_times_vec(A,VV,V); /* V = A*VV */
  
  alpha = dotprod(n,VV,V);
  
  for (i=0; i<n; i++) V[i] -= (alpha/beta1) * Vold[i];
  
  /* local reorthogonalization */

  numer = dotprod(n,VV,V);
  denom = dotprod(n,VV,VV);

  for (i=0; i<n; i++) V[i] -= (numer/denom) * VV[i];

  for (i=0; i<n; i++) Volder[i] = Vold[i];
  for (i=0; i<n; i++) Vold[i]   = V[i];
  
  if (precond_fn)
    (*precond_fn)(precond_args,V,Vold);
  else
    for (i=0; i<n; i++) V[i] = Vold[i];
  
  betaold = beta1;
  beta = dotprod(n,Vold,V);
  if (beta < 0.0) {
    taucs_printf("minres: error (2)\n");
    return -1;
  }
  beta = sqrt(beta);
  
  gammabar = alpha;
  epsilon = 0.0;
  deltabar = beta;
  gamma = sqrt(gammabar*gammabar + beta*beta);


  for (i=0; i<n; i++) Mold[i] = 0.0;
  for (i=0; i<n; i++) M[i]    = VV[i] / gamma;

  cs = gammabar / gamma;
  sn = beta / gamma;


  for (i=0; i<n; i++) X[i] += snprod*cs*M[i];
  snprod = snprod * sn;

  /* generate CG iterates */
  for (i=0; i<n; i++) Xcg[i] = X[i] + snprod*(sn/cs)*M[i];

  /* compute residual again */
  
  taucs_ccs_times_vec(A,X,R); 
  for (i=0; i<n; i++) R[i] = B[i] - R[i];  /* r = b - A*x */
  normr = twonorm(n,R);

  taucs_printf("minres: starting iterations, residual norm is %.1e\n",normr);
  
  for ( Iter=1; Iter <= itermax; Iter++ ) {

    for (i=0; i<n; i++) VV[i] = V[i] / beta;
    taucs_ccs_times_vec(A,VV,V); 
    for (i=0; i<n; i++) V[i] -= (beta/betaold) * Volder[i];
    alpha = dotprod(n,VV,V);
    for (i=0; i<n; i++) V[i] -= (alpha/beta) * Vold[i];

    for (i=0; i<n; i++) Volder[i] = Vold[i];
    for (i=0; i<n; i++) Vold  [i] = V   [i];
    
    if (precond_fn)
      (*precond_fn)(precond_args,V,Vold);
    else
      for (i=0; i<n; i++) V[i] = Vold[i];

    betaold = beta;
    beta = dotprod(n,Vold,V);
    if (beta < 0.0) {
      taucs_printf("minres: error (3)\n");
      return -1;
    }
    beta = sqrt(beta);

    delta = cs*deltabar + sn*alpha;
    for (i=0; i<n; i++) Molder[i] = Mold[i];
    for (i=0; i<n; i++) Mold  [i] = M   [i];
    for (i=0; i<n; i++) M[i] = VV[i] - delta*Mold[i] - epsilon*Molder[i];
    gammabar = sn*deltabar - cs*alpha;
    epsilon = sn*beta;
    deltabar = -cs*beta;
    gamma = sqrt(gammabar*gammabar + beta*beta);
    for (i=0; i<n; i++) M[i] = M[i]/ gamma;
    cs = gammabar / gamma;
    sn = beta / gamma;

    /* stagnation test; skipped */
    
    for (i=0; i<n; i++) X[i] += snprod*cs*M[i];
    snprod = snprod*sn;
    for (i=0; i<n; i++) Xcg[i] = X[i] + snprod*(sn/cs)*M[i];
    
    if (precond_fn) {
      taucs_ccs_times_vec(A,X,R); 
      for (i=0; i<n; i++) R[i] = B[i] - R[i];  /* r = b - A*x */
      normr = twonorm(n,R);
    } else {
      normr = fabs(snprod); 
      if (normr <= tolb) {
	/* double check */
	taucs_ccs_times_vec(A,X,R); 
	for (i=0; i<n; i++) R[i] = B[i] - R[i];  /* r = b - A*x */
	normr = twonorm(n,R);
      }
    }

    if (Iter > -1)
      taucs_printf("minres: n=%d iterations = %d residual norm %12.4e\n", A->n,Iter,normr);

    if (normr <= tolb) break;
  }

  taucs_printf("minres: done. n=%d iterations = %d residual norm %12.4e\n", A->n,Iter,normr);
 
  taucs_free(Molder) ;
  taucs_free(Mold) ;
  taucs_free(M) ;
  taucs_free(Volder) ;
  taucs_free(Vold) ;
  taucs_free(V) ;
  taucs_free(VV) ;
  taucs_free(Xcg) ;
  taucs_free(R) ;
 
  return 0; 
}                                                                             

#endif /* TAUCS_CORE_DOUBLE */

