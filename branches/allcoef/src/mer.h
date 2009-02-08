#ifndef LME4_LMER_H
#define LME4_LMER_H

#ifdef	__cplusplus
#define NO_C_HEADERS
#include <cstdio>
#include <cmath>
#include <cstring>
#endif 

#include <R.h>
/* Rdefines.h includes Rinternals.h (for SEXP, REAL, etc.) and defines
 * GET_SLOT, MAKE_CLASS, NEW_OBJECT, SET_SLOT, etc. */
#include <Rdefines.h>
#ifdef	__cplusplus
extern "C" {
#endif

// SEXP lme4_rWishart(SEXP ns, SEXP dfp, SEXP scal); 
SEXP lme4_ghq(SEXP np);

//SEXP merMCMC_VarCorr(SEXP x, SEXP typ);
//SEXP merMCMC_validate(SEXP x);

SEXP mer_A_to_U(SEXP rho);
SEXP mer_MCMCsamp(SEXP x, SEXP fm);
SEXP mer_PIRLS(SEXP rho);
//SEXP mer_postVar(SEXP x, SEXP which);
SEXP mer_update_dev(SEXP rho);
SEXP mer_update_mu(SEXP rho);
SEXP mer_validate(SEXP rho);

//SEXP spR_optimize(SEXP x, SEXP verbP);
//SEXP spR_update_mu(SEXP x);

#ifdef	__cplusplus
}

#endif /* __cplusplus */

#endif /* LME4_LMER_H */
