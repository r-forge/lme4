#include "Matrix.h"
#include <R_ext/Rdynload.h>

extern "C" SEXP glmIRLS(SEXP,SEXP);
extern "C" SEXP lme4_ghq(SEXP);
//extern "C" SEXP lme4_PermChk(SEXP,SEXP);
extern "C" SEXP merDeviance(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
extern "C" SEXP reTrmsCondVar(SEXP,SEXP);
extern "C" SEXP testFam(SEXP,SEXP,SEXP);
extern "C" SEXP updateDc(SEXP,SEXP,SEXP,SEXP);

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

static R_CallMethodDef CallEntries[] = {
    CALLDEF(glmIRLS, 2),
    CALLDEF(lme4_ghq, 1),
//    CALLDEF(lme4_PermChk, 2),
    CALLDEF(merDeviance, 6),
    CALLDEF(reTrmsCondVar, 2),
    CALLDEF(testFam, 3),
    CALLDEF(updateDc, 4),
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
}

/** Finalizer for lme4 called upon unloading the package.
 *
 */
extern "C"
void R_unload_lme4a(DllInfo *dll){
    M_cholmod_finish(&c);
}
