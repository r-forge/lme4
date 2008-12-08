#include "mer.h"
#include "ST.h"
#include <R_ext/Rdynload.h>
#include "Matrix.h"
#include "Syms.h" 

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

static R_CallMethodDef CallEntries[] = {
  
    CALLDEF(ST_Lambda, 1),
    CALLDEF(ST_Tmatrix, 1),
    CALLDEF(ST_bounds, 1),
    CALLDEF(ST_chol, 1),
    CALLDEF(ST_create_A, 2),
    CALLDEF(ST_create_ranef, 3),
    CALLDEF(ST_getPars, 1),
    CALLDEF(ST_initialize, 2),
    CALLDEF(ST_setPars, 3),
    CALLDEF(ST_update_A, 2),
    CALLDEF(ST_validate, 1),

    CALLDEF(lme4_ghq, 1),

    CALLDEF(merMCMC_VarCorr, 2),
    CALLDEF(merMCMC_validate, 1),
    
    CALLDEF(mer_A_to_U, 1),
    CALLDEF(mer_MCMCsamp, 2),
    CALLDEF(mer_PIRLS, 1),
    CALLDEF(mer_postVar, 2),
    CALLDEF(mer_update_dev, 1),
    CALLDEF(mer_update_mu, 1),
    CALLDEF(mer_validate, 1),
    
/*    CALLDEF(spR_optimize, 2), */
/*    CALLDEF(spR_update_mu, 1), */

  {NULL, NULL, 0}
};

/** cholmod_common struct local to the package */
cholmod_common c;

/** Initializer for lme4, called upon loading the package.
 *
 *  Register routines that can be called directly from R.
 *  Initialize CHOLMOD and require the LL' form of the factorization.
 *  Install the symbols to be used by functions in the package.
 */
#ifdef HAVE_VISIBILITY_ATTRIBUTE
__attribute__ ((visibility ("default")))
#endif
void R_init_lme4(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);


    M_R_cholmod_start(&c);
    c.final_ll = 1;	    /* LL' form of simplicial factorization */

    lme4_ASym = install("A");
    lme4_DimSym = install("Dim");
    lme4_GpSym = install("Gp");
    lme4_LSym = install("L");
    lme4_RXSym = install("RX");
    lme4_RZXSym = install("RZX");
    lme4_STSym = install("ST");
    lme4_XSym = install("X");
    lme4_ZtSym = install("Zt");
    lme4_devianceSym = install("deviance");
    lme4_dimsSym = install("dims");
    lme4_etaGammaSym = install("etaGamma");
    lme4_etaSym = install("eta");
    lme4_fixefSym = install("fixef");
    lme4_fl1Sym = install("fl1");
    lme4_flistSym = install("flist");
    lme4_ghwSym = install("ghw");
    lme4_ghxSym = install("ghx");
    lme4_gradientSym = install("gradient");
    lme4_iSym = install("i");
    lme4_merSym = install("mer");
    lme4_muEtaSym = install("muEta");
    lme4_muSym = install("mu");
    lme4_nlenvSym = install("nlenv");
    lme4_nlmodelSym = install("nlmodel");
    lme4_offsetSym = install("offset");
    lme4_pSym = install("p");
    lme4_pWtSym = install("pWt");
    lme4_permSym = install("perm");
    lme4_ranefSym = install("ranef");
    lme4_residSym = install("resid");
    lme4_rhoSym = install("rho");
    lme4_sigmaSym = install("sigma");
    lme4_sqrtrWtSym = install("sqrtrWt");
    lme4_uSym = install("u");
    lme4_varSym = install("var");
    lme4_xSym = install("x");
    lme4_ySym = install("y");
}

/** Finalizer for lme4 called upon unloading the package.
 *
 */
void R_unload_lme4(DllInfo *dll){
    M_cholmod_finish(&c);
}
