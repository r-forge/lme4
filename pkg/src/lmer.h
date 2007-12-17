#ifndef LME4_LMER_H
#define LME4_LMER_H

#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <Rversion.h>
#include <R_ext/Lapack.h>
#include <R_ext/stats_package.h>
#include "Matrix.h"

#ifdef HAVE_VISIBILITY_ATTRIBUTE
# define attr_hidden __attribute__ ((visibility ("hidden")))
#else
# define attr_hidden
#endif

#ifdef __GNUC__
# undef alloca
# define alloca(x) __builtin_alloca((x))
#else
/* this is necessary (and sufficient) for Solaris 10: */
# ifdef __sun
#  include <alloca.h>
# endif
#endif


extern
#include "Syms.h"

#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("lme4", String)
#else
#define _(String) (String)
#endif

SEXP lme4_rWishart(SEXP ns, SEXP dfp, SEXP scal);

SEXP lmer_profiled_deviance(SEXP x);

SEXP mer_MCMCsamp(SEXP x, SEXP savebp, SEXP nsampp,
		  SEXP transp, SEXP verbose, SEXP deviance); 
SEXP mer_ST_getPars(SEXP x);
SEXP mer_ST_initialize(SEXP ST, SEXP Gp, SEXP Zt);
SEXP mer_ST_setPars(SEXP x, SEXP pars);
SEXP mer_condMode(SEXP x, SEXP verbP);
SEXP mer_create_L(SEXP Vt);
SEXP mer_create_Vt(SEXP Zt, SEXP ST, SEXP Gp);
SEXP mer_eta(SEXP x);
SEXP mer_optimize(SEXP x, SEXP verb);
SEXP mer_pwrss(SEXP x);
SEXP mer_postVar(SEXP x, SEXP useScale);
SEXP mer_sigma(SEXP x, SEXP which);
SEXP mer_update_L(SEXP x);
SEXP mer_update_dev(SEXP x);
SEXP mer_update_effects(SEXP x);
SEXP mer_update_eta(SEXP x);
SEXP mer_update_mu(SEXP x);
SEXP mer_update_swts(SEXP x);
SEXP mer_validate(SEXP x);

SEXP nlmer_create_A(SEXP Vt, SEXP sP);

SEXP pedigree_chol(SEXP x, SEXP ans);

#endif /* LME4_LMER_H */
