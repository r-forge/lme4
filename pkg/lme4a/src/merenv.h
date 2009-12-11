#ifdef	__cplusplus
#define NO_C_HEADERS
#include <cstdio>
#include <cmath>
#include <cstring>
extern "C" {
#endif 

#include <R.h>
// Rdefines.h includes Rinternals.h (for SEXP, REAL, etc.) and defines
// GET_SLOT, MAKE_CLASS, NEW_OBJECT, SET_SLOT, etc.
#include <Rdefines.h>
#include "Matrix.h"

    SEXP glmer_linkinv(SEXP rho);

    SEXP lme4_dup_env_contents(SEXP dest, SEXP src, SEXP nms);

    SEXP lmerenv_deviance(SEXP rho, SEXP newth);
    SEXP lmerenv_validate(SEXP rho);

    SEXP merenvtrms_condVar(SEXP rho, SEXP scale);
    SEXP merenvtrms_show(SEXP rho);
    SEXP merenvtrms_validate(SEXP rho);

#ifdef	__cplusplus
}
#endif 
