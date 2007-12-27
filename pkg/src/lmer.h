#ifndef LME4_LMER_H
#define LME4_LMER_H

#include <R.h>
/* Rdefines.h includes Rinternals.h (for SEXP, REAL, etc.) and defines
 * GET_SLOT, MAKE_CLASS, NEW_OBJECT, SET_SLOT, etc. */
#include <Rdefines.h>

SEXP lme4_rWishart(SEXP ns, SEXP dfp, SEXP scal);

SEXP mer_MCMCsamp(SEXP x, SEXP savebp, SEXP nsampp,
		  SEXP transp, SEXP verbose, SEXP deviance); 
SEXP mer_ST_getPars(SEXP x);
SEXP mer_ST_initialize(SEXP ST, SEXP Gp, SEXP Zt);
SEXP mer_ST_setPars(SEXP x, SEXP pars);
SEXP mer_condMode(SEXP x, SEXP verbP);
SEXP mer_create_L(SEXP Vt);
SEXP mer_create_Vt(SEXP Zt, SEXP ST, SEXP Gp);
SEXP mer_optimize(SEXP x, SEXP verb);
SEXP mer_postVar(SEXP x, SEXP useScale);
SEXP mer_profiled_deviance(SEXP x);
SEXP mer_sigma(SEXP x, SEXP which);
SEXP mer_update_L(SEXP x);
SEXP mer_update_dev(SEXP x);
SEXP mer_update_effects(SEXP x);
SEXP mer_update_mu(SEXP x);
SEXP mer_validate(SEXP x);

SEXP nlmer_create_A(SEXP Vt, SEXP sP);

SEXP pedigree_chol(SEXP x, SEXP ans);

#endif /* LME4_LMER_H */
