#ifndef LME4_ST_H
#define LME4_ST_H

#include <R.h>
/* Rdefines.h includes Rinternals.h (for SEXP, REAL, etc.) and defines
 * GET_SLOT, MAKE_CLASS, NEW_OBJECT, SET_SLOT, etc. */
#include <Rdefines.h>
#ifdef	__cplusplus
extern "C" {
#endif

SEXP mer_ST_chol(SEXP x);
SEXP mer_ST_getPars(SEXP x);
SEXP mer_ST_initialize(SEXP ST, SEXP Gp, SEXP Zt);
SEXP mer_ST_setPars(SEXP x, SEXP pars);

SEXP spR_optimize(SEXP x, SEXP verbP);
SEXP spR_update_mu(SEXP x);

#ifdef	__cplusplus
}
#endif

#endif /* LME4_ST_H */
