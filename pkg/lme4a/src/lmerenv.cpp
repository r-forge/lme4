#include "lmerenv.h"
#include "lme4utils.hpp"

lmerenv::lmerenv(SEXP rho) : merenv(rho)
{
    REML = asLogical(findVarBound(rho, install("REML")));
    if (N != n)
	error(_("nrow(X) = %d must match length(y) = %d for lmer"),
	      N, n);
    Zty = VAR_dMatrix_x(rho, install("Zty"), q, 1);
    ldRX2 = VAR_REAL_NULL(rho, install("ldRX2"), 1, TRUE);
    if (X) {			// dense X
	RX = VAR_dMatrix_x(rho, lme4_RXSym, p, p);
	RZX = VAR_dMatrix_x(rho, lme4_RZXSym, q, p);
	XtX = VAR_dMatrix_x(rho, install("XtX"), p, p);
	ZtX = VAR_dMatrix_x(rho, install("ZtX"), q, p);
	Xty = VAR_REAL_NULL(rho, install("Xty"), p);
    } else {
	error(_("code not yet written"));
    }
}

/**
 * Update to new value of theta and evaluate the deviance or REML
 * criterion.
 *
 * @param thnew pointer to an numeric vector theta
 *
 * @return deviance or REML criterion according to the value of REML
 */
double lmerenv::update_dev(SEXP thnew) {
    double mone[2] = {-1,0}, one[2] = {1,0};
    CHM_DN cZty = N_AS_CHM_DN(Zty, q, 1);
    int info;

    update_Lambda_Ut(thnew);
    M_cholmod_factorize_p(Ut, one, (int*)NULL, (size_t)0, L, &c);
    *ldL2 = M_chm_factor_ldetL2(L);
    if (X) {
	CHM_DN cZtX = N_AS_CHM_DN(ZtX, q, p), cu, tmp1, tmp2;
				// create cu and RZX
	tmp1 = M_cholmod_copy_dense(cZty, &c);
	crossprod_Lambda(cZty, tmp1);
	tmp2 = M_cholmod_solve(CHOLMOD_P, L, tmp1, &c);
	M_cholmod_free_dense(&tmp1, &c);
	cu = M_cholmod_solve(CHOLMOD_L, L, tmp2, &c);
	M_cholmod_free_dense(&tmp2, &c);
	tmp1 = M_cholmod_copy_dense(cZtX, &c);
	crossprod_Lambda(cZtX, tmp1);
	tmp2 = M_cholmod_solve(CHOLMOD_P, L, tmp1, &c);
	M_cholmod_free_dense(&tmp1, &c);
	tmp1 = M_cholmod_solve(CHOLMOD_L, L, tmp2, &c);
	M_cholmod_free_dense(&tmp2, &c);
	dble_cpy(RZX, (double*)(tmp1->x), q * p);
	M_cholmod_free_dense(&tmp1, &c);
				// downdate and factor XtX, solve for fixef
	dble_cpy(RX, XtX, p * p);
	F77_CALL(dsyrk)("U", "T", &p, &q, mone, RZX, &q, one, RX, &p);
	dble_cpy(fixef, Xty, p);
	F77_CALL(dgemv)("T", &q, &p, mone, RZX, &q, (double*)(cu->x),
			&i1, one, fixef, &i1);
	F77_CALL(dposv)("U", &p, &i1, RX, &p, fixef, &p, &info);
	if (info)
	    error(_("Downdated X'X is not positive definite, %d."), info);
				// evaluate ldRX2
	*ldRX2 = 0;
	for (int i = 0; i < p; i++) *ldRX2 += 2 * log(RX[i * (p + 1)]);
				// solve for u
	F77_CALL(dgemv)("N", &q, &p, mone, RZX, &q, fixef, &i1, one,
			(double*)(cu->x), &i1);
	tmp1 = M_cholmod_solve(CHOLMOD_Lt, L, cu, &c);
	M_cholmod_free_dense(&cu, &c);
	tmp2 = M_cholmod_solve(CHOLMOD_Pt, L, tmp1, &c);
	M_cholmod_free_dense(&tmp1, &c);
	dble_cpy(u, (double*)(tmp2->x), q);
	M_cholmod_free_dense(&tmp2, &c);
    } else {
	error(_("code not yet written"));
// Structures that are of different class with a sparseX are
// X (dgCMatrix), XtX (dpoMatrix), ZtX (dgCMatrix), RX (dtCMatrix),
// RZX (dgCMatrix) 
    }
    *prss = sqr_length(u, q);
    update_eta();
    for (int i = 0; i < n; i++) {
	double resi = y[i] - eta[i];
	*prss += resi * resi;
    }
    if (REML) {
	double nmp = (double)(n - p);
	return *ldL2 + *ldRX2 + nmp * (1 + log(2 * PI * (*prss)/nmp));
    }				       
    return *ldL2 + n * (1 + log(2 * PI * (*prss)/((double)n)));
}

/* Externally callable functions */

/**
 * Evaluate the deviance or REML criterion
 *
 * @param rho pointer to an merenv environment
 * @param thnew pointer to an numeric vector theta
 *
 * @return deviance value
 */
SEXP lmerenv_deviance(SEXP rho, SEXP thnew) {
    return ScalarReal(lmerenv(rho).update_dev(thnew));
}

/**
 * Check validity of an merenv environment
 *
 * @param x pointer to an merenv environment
 */
SEXP lmerenv_validate(SEXP rho) {
    return ScalarLogical(lmerenv(rho).validate());
}
