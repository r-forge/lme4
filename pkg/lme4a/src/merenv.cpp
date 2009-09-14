#include "merenv.h"

#include "lme4utils.hpp"

int merenv::i1 = 1;

// Definition of methods for the mer class

merenv::merenv(SEXP rho)
{
    // Extract slots that must have positive length.
    // Get dimensions of the problem
    SEXP sl = findVarBound(rho, lme4_ySym);
    if (!(n = LENGTH(sl)) || !isReal(sl)) // n = length of response
	error(_("Response vector y must be numeric (double)"));
    y = REAL(sl);

    sl = findVarBound(rho, lme4_fixefSym);
    if (!(p = LENGTH(sl)) || !isReal(sl)) // p = length of fixef
	error(_("fixef vector must be numeric (double)"));
    fixef = REAL(sl);

    sl = findVarBound(rho, lme4_uSym);
    if (!(q = LENGTH(sl)) || !isReal(sl)) // q = length of u
	error(_("u vector must be numeric (double)"));
    u = REAL(sl);

    sl = findVarBound(rho, install("theta"));
    if (!(nth = LENGTH(sl)) || !isReal(sl)) // nth = length of theta
	error(_("theta vector must be numeric (double)"));
    theta = REAL(sl);

    sl = findVarBound(rho, install("Lind"));
    if (!(nLind = LENGTH(sl)) || !isInteger(sl))
	error(_("Lind vector must be integer"));
    Lind = INTEGER(sl);
    
    // allow for Z and X to have a multiple of n rows.
    Zt = VAR_CHM_SP(rho, lme4_ZtSym, q, 0);
    N = (int)Zt->ncol;
    Ut = VAR_CHM_SP(rho, lme4_UtSym, q, N);
    L = new cholmod_factor;
    M_as_cholmod_factor(L, findVarInFrame(rho, lme4_LSym));
				// versions of Lambda
//FIXME: Use S4 class to determine sparseX and diagonalLambda
    if (asLogical(findVarBound(rho, install("diagonalLambda")))) {
	Lambdax = VAR_dMatrix_x(rho, lme4_LambdaSym, q, q);
	if (nLind != q)
	    error(_("Lind should be of length q = %d for diagonal Lambda"),
		  q);
	Lambda = (CHM_SP) NULL;
    } else {
	Lambda = VAR_CHM_SP(rho, lme4_LambdaSym, q, q);
	int nnz = M_cholmod_nnz(Lambda, &c);
	if (nnz != nLind)
	    error(_("length(Lind) = %d should be  %d"), nLind, nnz);
	Lambdax = (double*)(Lambda->x);
    }
    eta = VAR_REAL_NULL(rho, install("fitted"), N);
    ldL2 = VAR_REAL_NULL(rho, install("ldL2"), 1);
    prss = VAR_REAL_NULL(rho, install("prss"), 1);
    weights = VAR_REAL_NULL(rho, install("weights"), n, TRUE);
    offset = VAR_REAL_NULL(rho, lme4_offsetSym, N, TRUE);
}

void merenv::update_Lambda_Ut(SEXP thnew) {
				// check and copy thnew
    if (!isReal(thnew) || LENGTH(thnew) != nth)
	error(_("theta must be numeric and of length %d"), nth);
    dble_cpy(theta, REAL(thnew), nth);
				// update Lambda
    for (int i = 0; i < nLind; i++) Lambdax[i] = theta[Lind[i] - 1];
				// update Ut from Lambda and Zt
    if (Lambda) {
	CHM_SP Lamtr = M_cholmod_transpose(Lambda, TRUE/*values*/, &c),
	    tmp = M_cholmod_ssmult(Lamtr, Zt, 0/*stype*/,
				   TRUE/*values*/, TRUE/*sorted*/, &c);
	// Should we store t(Lambda), instead of Lambda?
	int *it = (int*)(tmp->i), *iu = (int*)(Ut->i),
	    *pt = (int*)(tmp->p), *pu = (int*)(Ut->p);
	double *xt = (double*)(tmp->x), *xu = (double*)(Ut->x); 

	M_cholmod_free_sparse(&Lamtr, &c);
	for (int j = 0; j <= n; j++)
	    if (pt[j] != pu[j])
		error(_("Ut is not consistent with Lambda %*% Zt"));
	for (int j = 0; j < pu[n]; j++) {
	    if (it[j] != iu[j])
		error(_("Ut is not consistent with Lambda %*% Zt"));
	    xu[j] = xt[j];
	}
	M_cholmod_free_sparse(&tmp, &c);
    } else {			// special case for diagonal Lambda
	int *iz = (int*)(Zt->i), nnz = ((int*)(Zt->p))[n];
	double *xu = (double*)(Ut->x), *xz = (double*)(Zt->x);
	for (int i = 0; i < nnz; i++) xu[i] = xz[i] * Lambdax[iz[i]];
    }
}

void merenv::update_eta() {    
    if (offset)			// initialize to offset if used
	dble_cpy(eta, offset, N);
    else dble_zero(eta, N);	// otherwise to zero
    CHM_DN ceta = N_AS_CHM_DN(eta, N, 1), cu = N_AS_CHM_DN(u, q, 1);
    double one[2] = {1,0};
    M_cholmod_sdmult(Ut, 1/*transpose*/, one, one, cu, ceta, &c);
}

void mersparse::update_eta() {
// call method for base class
    double one[2] = {1,0};
    CHM_DN cfixef = N_AS_CHM_DN(fixef, p, 1);
    M_cholmod_sdmult(sX, 0/*no transpose*/, one, one, cfixef, ceta, &c);
}

void merdense::update_eta() {
// call method for base class
    double one[2] = {1,0};
    F77_CALL(dgemv)("N", &n, &p, one, X, &n, fixef, &i1, one, eta, &i1);
}

CHM_DN merenv::crossprod_Lambda(CHM_DN rhs, CHM_DN ans) {
    if (((int)(rhs->nrow)) != q)
	error(_("in crossprod_Lambda, rhs->nrow = %d, should be %d"),
	      rhs->nrow, q);
    if (Lambda) {
	double one[2] = {1,0}, zero[2] = {0,0};
	M_cholmod_sdmult(Lambda, 1/*transpose*/, one, zero, rhs, ans, &c);
    } else {
	double *ax = (double*)(ans->x), *rx = (double*)(rhs->x);
	int nc = rhs->ncol;
    	for (int j = 0; j < nc * q; j++) ax[j] = rx[j] * Lambdax[j % q];
    }
    return ans;
}

merdense::merdense(SEXP rho) : merenv(rho) {
    X = VAR_dMatrix_x(rho, lme4_XSym, N, p);
    RX = VAR_dMatrix_x(rho, lme4_RXSym, p, p);
    RZX = VAR_dMatrix_x(rho, lme4_RZXSym, q, p);
}

mersparse::mersparse(SEXP rho) : merenv(rho) {
    X = VAR_CHM_SP(rho, lme4_XSym, N, p);
    RX = VAR_CHM_SP(rho, lme4_RXSym, p, p);
    RZX = VAR_CHM_SP(rho, lme4_RZXSym, q, p);
}

lmer::lmer(SEXP rho)
{
    REML = asLogical(findVarBound(rho, install("REML")));
    if (N != n)
	error(_("nrow(X) = %d must match length(y) = %d for lmer"),
	      N, n);
    Xty = VAR_REAL_NULL(rho, install("Xty"), p);
    Zty = VAR_dMatrix_x(rho, install("Zty"), q, 1);
    ldRX2 = VAR_REAL_NULL(rho, install("ldRX2"), 1, TRUE);
}

lmerdense::lmerdense(SEXP rho) : merdense(rho) : lmer(rho) {
    ZtX = VAR_dMatrix_x(rho, install("ZtX"), q, p);
    XtX = VAR_dMatrix_x(rho, install("XtX"), p, p);
}

lmersparse::lmersparse(SEXP rho) : mersparse(rho) : lmer(rho) {
    ZtX = VAR_CHM_SP(rho, install("ZtX"), q, p);
    XtX = VAR_CHM_SP(rho, install("XtX"), p, p);
}

/**
 * Update to new value of theta and evaluate the deviance or REML
 * criterion.
 *
 * @param thnew pointer to an numeric vector theta
 *
 * @return deviance or REML criterion according to the value of REML
 */
double lmerdense::update_dev(SEXP thnew) {
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
