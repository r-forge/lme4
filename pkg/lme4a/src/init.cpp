#include "merenv.h"
#include "Matrix.h"
#include "Syms.h"
#include <R_ext/Rdynload.h>

extern "C" SEXP LMMdeviance(SEXP);
extern "C" SEXP PIRLS(SEXP,SEXP,SEXP,SEXP);
extern "C" SEXP feSetBeta(SEXP,SEXP);
extern "C" SEXP reUpdateLambda(SEXP,SEXP);
extern "C" SEXP updateRzxRx(SEXP);
extern "C" SEXP updateDc(SEXP);

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

static R_CallMethodDef CallEntries[] = {
    CALLDEF(LMMdeviance, 1),
    CALLDEF(PIRLS, 4),

    CALLDEF(feSetBeta, 2),

    CALLDEF(glmer_IRLS, 1),
    CALLDEF(glmer_PIRLS, 1),
    CALLDEF(glmer_PIRLSbeta, 1),
    CALLDEF(glmer_update_RX, 1),


    CALLDEF(lme4_ghq, 1),
    CALLDEF(lme4_dup_env_contents, 3),

    CALLDEF(lmer_deviance, 2),
    CALLDEF(lmer_validate, 1),

    CALLDEF(merenv_update_Lambda_Ut, 2),

    CALLDEF(merenvtrms_condVar, 2),
    CALLDEF(merenvtrms_show, 1),
    CALLDEF(merenvtrms_validate, 1),

    CALLDEF(reUpdateLambda, 2),
    CALLDEF(updateRzxRx, 1),
    CALLDEF(updateDc, 1),
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
extern "C"
void R_init_lme4a(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, (Rboolean)FALSE);


    M_R_cholmod_start(&c);
    c.final_ll = 1;	    /* LL' form of simplicial factorization */

    lme4_DimSym = install("Dim");
    lme4_LSym = install("L");
    lme4_LambdaSym = install("Lambda");
    lme4_LindSym = install("Lind");
    lme4_RXSym = install("RX");
    lme4_RZXSym = install("RZX");
    lme4_UtSym = install("Ut");
    lme4_XSym = install("X");
    lme4_ZtSym = install("Zt");
    lme4_betaSym = install("beta");
    lme4_familySym = install("family");
    lme4_flistSym = install("flist");
    lme4_gammaSym = install("gamma");
    lme4_ghwSym = install("ghw");
    lme4_ghxSym = install("ghx");
    lme4_iSym = install("i");
    lme4_ldL2Sym = install("ldL2");
    lme4_muEtaSym = install("muEta");
    lme4_muSym = install("mu");
    lme4_nlmodelSym = install("nlmodel");
    lme4_offsetSym = install("offset");
    lme4_pSym = install("p");
    lme4_pwrssSym = install("pwrss");
    lme4_sqrtrwtSym = install("sqrtrWt");/* capital 'W' - as 'always' in "mer" class */
    lme4_sqrtXwtSym = install("sqrtXWt");
    lme4_thetaSym = install("theta");
    lme4_uSym = install("u");
    lme4_varSym = install("var");
    lme4_weightsSym = install("weights");
    lme4_xSym = install("x");
    lme4_ySym = install("y");
}

/** Finalizer for lme4 called upon unloading the package.
 *
 */
extern "C"
void R_unload_lme4(DllInfo *dll){
    M_cholmod_finish(&c);
}
