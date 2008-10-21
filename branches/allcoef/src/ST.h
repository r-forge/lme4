#ifndef LME4_ST_H
#define LME4_ST_H

#include <R.h>
/* Rdefines.h includes Rinternals.h (for SEXP, REAL, etc.) and defines
 * GET_SLOT, MAKE_CLASS, NEW_OBJECT, SET_SLOT, etc. */
#include <Rdefines.h>
#ifdef	__cplusplus
extern "C" {
#endif

    SEXP ST_generate_A(SEXP Zt, SEXP Gp, SEXP ST);
    SEXP ST_update_A(SEXP x);
    SEXP mer_ST_chol(SEXP x);
    SEXP mer_ST_getPars(SEXP x);
    SEXP mer_ST_initialize(SEXP ST, SEXP Gp, SEXP Zt);
    SEXP mer_ST_setPars(SEXP x, SEXP pars);

#ifdef	__cplusplus
}
#endif

#endif /* LME4_ST_H */
