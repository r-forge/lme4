#include "lmerenv.h"
#include "lme4utils.hpp"

lmerenv::lmerenv(SEXP rho) : merenv(rho)
{
    REML = asLogical(findVarBound(rho, install("REML")));
    if (N != n)
	error(_("nrow(X) = %d must match length(y) = %d for lmer"),
	      N, n);
    Xty = VAR_REAL_NULL(rho, install("Xty"), p);
    Zty = VAR_dMatrix_x(rho, install("Zty"), q, 1);
    ldRX2 = VAR_REAL_NULL(rho, install("ldRX2"), 1, TRUE);
    if (X) {			// dense X
	RX = VAR_dMatrix_x(rho, lme4_RXSym, p, p);
	RZX = VAR_dMatrix_x(rho, lme4_RZXSym, q, p);
	XtX = VAR_dMatrix_x(rho, install("XtX"), p, p);
	ZtX = VAR_dMatrix_x(rho, install("ZtX"), q, p);
    } else {
	sRX = VAR_CHM_SP(rho, lme4_RXSym, p, p);
	sRZX = VAR_CHM_SP(rho, lme4_RZXSym, p, p);
	sXtX = VAR_CHM_SP(rho, install("XtX"), p, p);
	sZtX = VAR_CHM_SP(rho, install("ZtX"), p, p);
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
    CHM_DN cZty = N_AS_CHM_DN(Zty, q, 1), cu, Dtmp1, Dtmp2;
    int info;

    update_Lambda_Ut(thnew);
    M_cholmod_factorize_p(Ut, one, (int*)NULL, (size_t)0, L, &c);
    *ldL2 = M_chm_factor_ldetL2(L);
				// create cu 
    Dtmp1 = M_cholmod_copy_dense(cZty, &c);
    crossprod_Lambda(cZty, Dtmp1);
    Dtmp2 = M_cholmod_solve(CHOLMOD_P, L, Dtmp1, &c);
    M_cholmod_free_dense(&Dtmp1, &c);
    cu = M_cholmod_solve(CHOLMOD_L, L, Dtmp2, &c);
    M_cholmod_free_dense(&Dtmp2, &c);
    if (X) {			// update dense RZX
	CHM_DN cZtX = N_AS_CHM_DN(ZtX, q, p);
	Dtmp1 = M_cholmod_copy_dense(cZtX, &c);
	crossprod_Lambda(cZtX, Dtmp1);
	Dtmp2 = M_cholmod_solve(CHOLMOD_P, L, Dtmp1, &c);
	M_cholmod_free_dense(&Dtmp1, &c);
	Dtmp1 = M_cholmod_solve(CHOLMOD_L, L, Dtmp2, &c);
	M_cholmod_free_dense(&Dtmp2, &c);
	dble_cpy(RZX, (double*)(Dtmp1->x), q * p);
	M_cholmod_free_dense(&Dtmp1, &c);
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
	Dtmp1 = M_cholmod_solve(CHOLMOD_Lt, L, cu, &c);
	Dtmp2 = M_cholmod_solve(CHOLMOD_Pt, L, Dtmp1, &c);
	M_cholmod_free_dense(&Dtmp1, &c);
	dble_cpy(u, (double*)(Dtmp2->x), q);
	M_cholmod_free_dense(&Dtmp2, &c);
    } else {
	CHM_SP Stmp1, Stmp2;
	error(_("code not yet written"));
// Structures that are of different class with a sparseX are
// X (dgCMatrix), XtX (dpoMatrix), ZtX (dgCMatrix), RX (dtCMatrix),
// RZX (dgCMatrix) 

// Need a crossprod_Lambda for sparse matrices.  Should have a safe
// copy of the contents of a sparse matrix to another (safe in the
// sense that it first checks for consistency).
    }
    M_cholmod_free_dense(&cu, &c);
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
