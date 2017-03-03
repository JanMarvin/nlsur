#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP nlsur_arma_reshape(SEXP, SEXP);
extern SEXP nlsur_calc_reg(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP nlsur_calc_robust(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP nlsur_calc_ssr(SEXP, SEXP, SEXP);
extern SEXP nlsur_wt_mean(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"nlsur_arma_reshape", (DL_FUNC) &nlsur_arma_reshape, 2},
  {"nlsur_calc_reg",     (DL_FUNC) &nlsur_calc_reg,     7},
  {"nlsur_calc_robust",  (DL_FUNC) &nlsur_calc_robust,  5},
  {"nlsur_calc_ssr",     (DL_FUNC) &nlsur_calc_ssr,     3},
  {"nlsur_wt_mean",      (DL_FUNC) &nlsur_wt_mean,      2},
  {NULL, NULL, 0}
};

void R_init_nlsur(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
