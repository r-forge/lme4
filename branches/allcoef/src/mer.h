#ifndef LME4_LMER_H
#define LME4_LMER_H

#include <R.h>
/* Rdefines.h includes Rinternals.h (for SEXP, REAL, etc.) and defines
 * GET_SLOT, MAKE_CLASS, NEW_OBJECT, SET_SLOT, etc. */
#include <Rdefines.h>
#ifdef	__cplusplus
extern "C" {
#endif

// SEXP lme4_rWishart(SEXP ns, SEXP dfp, SEXP scal); 
SEXP lme4_ghq(SEXP np);

SEXP merMCMC_VarCorr(SEXP x, SEXP typ);
SEXP merMCMC_validate(SEXP x);

SEXP mer_MCMCsamp(SEXP x, SEXP fm);
SEXP mer_PIRLS(SEXP x);
SEXP mer_create_L(SEXP CmP);
SEXP mer_optimize(SEXP x);
SEXP mer_postVar(SEXP x, SEXP which);
SEXP mer_update_dev(SEXP x);
SEXP mer_update_mu(SEXP x);
SEXP mer_validate(SEXP x);

//SEXP spR_optimize(SEXP x, SEXP verbP);
//SEXP spR_update_mu(SEXP x);

#ifdef	__cplusplus
}

#include "Matrix.h"
#include "lme4utils.hpp"

class mer {
public:
    mer(SEXP x);		//< external, R-level object
    ~mer(){
	delete L; delete A; 
	d[wrss_POS] = wrss;
	d[usqr_POS] = usqr;
	d[pwrss_POS] = pwrss;
    }
/**
 * Update the u and fixef slots in an mer object
 * @return updated deviance
 */
    double PIRLS();
/**
 * Update the conditional mean, mu
 * @return penalized, weighted residual sum of squares
 *
 * \note This function uses the existing weights without updating
 * them. Reweighting must be done explicitly in PIRLS.
 */
    double update_mu();
/**
 * Evaluate the deviance, possibly using AGQ.
 */
    double update_dev();

private:
    static int i1;
    static const int BUF_SIZE = 127, CM_MAXITER = 300;
    static double mone, one, zero;
    static const double CM_TOL, CM_SMIN, GHQ_EPS,
	LTHRESH, MLTHRESH, MPTHRESH, PTHRESH, INVEPS;

    int *dims, *perm, N, n, p, q, s;
    double *RX, *RZX, *V, *X, *d, *eta, *fixef, *etaGamma, *mu,
	*muEta, *offset, *pWt, *srwt, *res, *u, *var, *y,
	*ghx, *ghw, wrss, pwrss, usqr;
    SEXP flistP, nlmodel, pnames, rho;
    CHM_FR L;
    CHM_SP A;

    void extractL(SEXP x);
    void extractA(SEXP x);
    double *apply_perm(double *dest, const double *src)
    {
	for (int i = 0; i < q; i++) dest[i] = src[perm ? perm[i] : i];
	return dest;
    }

    double* apply_iperm(double *dest, const double *src)
    {
	for (int i = 0; i < q; i++) dest[perm ? perm[i] : i] = src[i];
	return dest;
    }
/**
 * Create a sparse matrix U from A, etaGamma, muEta, srwt and Perm
 * @return a freshly allocated q by n sparse matrix
 */
    CHM_SP A_to_U();
/**
 * Fill in the V matrix using X, etaGamma, muEta and srwt.
 * @param V pointer to storage for N * p doubles
 */
    double* X_to_V(double *V);
    void eval_nonlin(const double *tmp);
    void eval_muEta();
    void eval_varFunc();
    double* eval_devResid(double *ans, const int *Grps);
};

#endif /* __cplusplus */

#endif /* LME4_LMER_H */
