#ifndef LME4_ST_H
#define LME4_ST_H

#include <R.h>
/* Rdefines.h includes Rinternals.h (for SEXP, REAL, etc.) and defines
 * GET_SLOT, MAKE_CLASS, NEW_OBJECT, SET_SLOT, etc. */
#include <Rdefines.h>
#include <R_ext/Lapack.h>        // for Lapack (dpotrf, etc.) and BLAS
#include "Matrix.h"

#ifdef	__cplusplus
extern "C" {
#endif

    SEXP ST_Lambda(SEXP x);
    SEXP ST_Tmatrix(SEXP x);
    SEXP ST_bounds(SEXP x);
    SEXP ST_chol(SEXP x);
    SEXP ST_create_A(SEXP x, SEXP rho);
    SEXP ST_create_ranef(SEXP ST, SEXP u, SEXP perm);
    SEXP ST_getPars(SEXP x);
    SEXP ST_initialize(SEXP x, SEXP rho);
    SEXP ST_postVar(SEXP x, SEXP L, SEXP perm, SEXP flist, SEXP which);
    SEXP ST_setPars(SEXP x, SEXP pars, SEXP rho);
    SEXP ST_update_A(SEXP ST, SEXP rho);
    SEXP ST_validate(SEXP x);

#ifdef	__cplusplus
}
#endif /* __cplusplus */

#endif /* LME4_ST_H */
