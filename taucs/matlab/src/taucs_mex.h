#include "mex.h"
#include "taucs.h"


void mex_fill_taucs_ccs_matrix(const mxArray *matlab_A, taucs_ccs_matrix *taucs_A);
mxArray *mex_convert_taucs_ccs_2_matlab(taucs_ccs_matrix *M);
void mex_free_ccs_matrix(taucs_ccs_matrix *a);
double mex_get_cpu_time();
