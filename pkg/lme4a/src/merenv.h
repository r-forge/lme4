#ifndef LME4_MERENV_H
#define LME4_MERENV_H

#ifdef	__cplusplus
#define NO_C_HEADERS
#include <cstdio>
#include <cmath>
#include <cstring>
#endif 

#include <R.h>
// Rdefines.h includes Rinternals.h (for SEXP, REAL, etc.) and defines
// GET_SLOT, MAKE_CLASS, NEW_OBJECT, SET_SLOT, etc.
#include <Rdefines.h>
#ifdef	__cplusplus
extern "C" {
#endif

    SEXP lmerenv_deviance(SEXP rho, SEXP thnew);
    SEXP lmerenv_validate(SEXP rho);

#ifdef	__cplusplus
}
#endif /* __cplusplus */

#endif /* LME4_MERENV_H */
