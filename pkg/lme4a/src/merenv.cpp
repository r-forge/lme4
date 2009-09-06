#include "merenv.h"

#include "Matrix.h"
#include "lme4utils.hpp"

class merenv {
public:
    merenv(SEXP rho);		//< instantiate from an environment
    ~merenv(){
	delete L;
	delete Zt;
	delete Ut;
	if (!dLam) delete Lambda;
    }
    double update_dev(SEXP thnew);
    int validate();

private:
    static int i1;
    int *Lind, *perm, REML, dLam, n, nLind, nmp, nth,
	p, q, sparseX;
    double *Lambdax, *RX, *RZX, *X, *XtX, *Xty, *ZtX, *Zty,
	*beta, *eta, *ldL2, *ldRX2, *offset, *pWt,
	*prss, *theta, *u, *weights, *y;
    CHM_FR L;
    CHM_SP Lambda, Ut, Zt, ZtsX, sX, sXtX, sRX, sRZX;
    void update_Lambda();
    void update_Ut();
    void update_eta();
    CHM_DN crossprod_Lambda(CHM_DN rhs, CHM_DN ans);
};

int merenv::i1 = 1;

// Inlined utilties

/**
 * Copy the first nn elements of src to dest
 *
 * @param src source vector
 * @param dest destination vector
 * @param nn number of elements in src and dest
 *
 * @return dest
 */
inline double *dble_cpy(double *dest, const double *src, int nn)
{
    for (int i = 0; i < nn; i++)
	dest[i] = src[i];
    return dest;
}

/**
 * Zero the first nn elements of double pointer dest
 *
 * @param dest vector
 * @param nn number of elements in dest
 *
 * @return dest
 */
inline double *dble_zero(double *dest, int nn)
{
    for (int i = 0; i < nn; i++)
	dest[i] = 0.;
    return dest;
}

/**
 * Zero the first nn elements of int pointer dest
 *
 * @param dest vector
 * @param nn number of elements in dest
 *
 * @return dest
 */
inline int *int_zero(int *dest, int nn)
{
    for (int i = 0; i < nn; i++)
	dest[i] = 0;
    return dest;
}

/**
 * Evaluate the squared length of the first nn elements of x
 *
 * @param x vector
 * @param nn number of elements in x
 *
 * @return the squared length of x
 */
inline double sqr_length(const double *x, int nn) {
    double ans = 0;
    for (int i = 0; i < nn; i++)
	ans += x[i] * x[i];
    return ans;
}

inline SEXP findVarBound(SEXP rho, SEXP nm) {
    SEXP var = findVarInFrame(rho, nm);
    if (var == R_UnboundValue)
	error(_("object named '%s' not found in environment"),
	      CHAR(PRINTNAME(nm)));
    return var;
}
    
/** Non-inlined utilities
 */

CHM_SP VAR_CHM_SP(SEXP rho, SEXP nm, int nrow, int ncol)
{
    CHM_SP ans = new cholmod_sparse;
    SEXP var = findVarBound(rho, nm);
    M_as_cholmod_sparse(ans, var, (Rboolean)TRUE, (Rboolean)FALSE);
    if (((int)(ans->nrow)) != nrow)
	error(_("Number of rows of %s is %d, should be %d"),
	      CHAR(PRINTNAME(nm)), ans->nrow, nrow);
    if (((int)(ans->ncol)) != ncol)
	error(_("Number of columns of %s is %d, should be %d"),
	      CHAR(PRINTNAME(nm)), ans->ncol, ncol);
    return ans;
}

double *VAR_dMatrix_x(SEXP rho, SEXP nm, int nrow, int ncol)
{    
    SEXP var = findVarBound(rho, nm);

    // FIXME: Should check here to ensure that the object "is" a dMatrix
    int *dims = INTEGER(GET_SLOT(var, lme4_DimSym));
    if (dims[0] != nrow || dims[1] != ncol)
	error(_("object named '%s' should be %d by %d"),
	      CHAR(PRINTNAME(nm)), nrow, ncol);
    return REAL(GET_SLOT(var, lme4_xSym));
}


// Definition of methods for the mer class

merenv::merenv(SEXP rho)
{
    // Extract slots that must have positive length.
    // Get dimensions of the problem
    SEXP sl = findVarBound(rho, lme4_ySym);
    if (!(n = LENGTH(sl)) || !isReal(sl)) // n = length of response
	error(_("Response vector y must be numeric (double)"));
    y = REAL(sl);

    sl = findVarBound(rho, lme4_betaSym);
    if (!(p = LENGTH(sl)) || !isReal(sl)) // p = length of beta
	error(_("beta vector must be numeric (double)"));
    beta = REAL(sl);

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
    
    Zt = VAR_CHM_SP(rho, lme4_ZtSym, q, n);
    Ut = VAR_CHM_SP(rho, lme4_UtSym, q, n);
    L = new cholmod_factor;
    M_as_cholmod_factor(L, findVarInFrame(rho, lme4_LSym));
				// scalar logicals
    REML = asLogical(findVarBound(rho, install("REML")));
    dLam = asLogical(findVarBound(rho, install("diagonalLambda")));
    sparseX = asLogical(findVarBound(rho, install("sparseX")));
				// versions of Lambda
    if (dLam) {
	Lambdax = VAR_dMatrix_x(rho, lme4_LambdaSym, q, q);
	if (nLind != q)
	    error(_("Lind should be of length q = %d for diagonal Lambda"),
		  q);
    } else {
	Lambda = VAR_CHM_SP(rho, lme4_LambdaSym, q, q);
	int nnz = M_cholmod_nnz(Lambda, &c);
	if (nnz != nLind)
	    error(_("length(Lind) = %d should be  %d"), nLind, nnz);
	Lambdax = (double*)(Lambda->x);
    }
    if (sparseX) {
	ZtsX = VAR_CHM_SP(rho, install("ZtX"), q, p);
	sX = VAR_CHM_SP(rho, lme4_XSym, n, p);
	sXtX = VAR_CHM_SP(rho, install("XtX"), p, p);
	sRX = VAR_CHM_SP(rho, lme4_RXSym, p, p);
	sRZX = VAR_CHM_SP(rho, lme4_RZXSym, q, p);
    } else {
	RX = VAR_dMatrix_x(rho, lme4_RXSym, p, p);
	RZX = VAR_dMatrix_x(rho, lme4_RZXSym, q, p);
	X = VAR_dMatrix_x(rho, lme4_XSym, n, p);
	XtX = VAR_dMatrix_x(rho, install("XtX"), p, p);
	ZtX = VAR_dMatrix_x(rho, install("ZtX"), q, p);
    }
    Xty = VAR_REAL_NULL(rho, install("Xty"), p);
    eta = VAR_REAL_NULL(rho, install("fitted"), n);
    ldL2 = VAR_REAL_NULL(rho, install("ldL2"), 1);
    ldRX2 = VAR_REAL_NULL(rho, install("ldRX2"), 1);
    prss = VAR_REAL_NULL(rho, install("prss"), 1);
    Zty = VAR_dMatrix_x(rho, install("Zty"), q, 1);
    weights = VAR_REAL_NULL(rho, install("weights"), n, TRUE);
    offset = VAR_REAL_NULL(rho, lme4_offsetSym, q, TRUE);
}

void merenv::update_Lambda() {
    for (int i = 0; i < nLind; i++) Lambdax[i] = theta[Lind[i] - 1];
}	

void merenv::update_Ut() {
    if (dLam) {
	int *iz = (int*)(Zt->i), nnz = ((int*)(Zt->p))[n];
	double *xu = (double*)(Ut->x), *xz = (double*)(Zt->x);
	for (int i = 0; i < nnz; i++) xu[i] = xz[i] * Lambdax[iz[i]];
    } else {
				// store Lambda transpose, not Lambda?
	CHM_SP Lamtr = M_cholmod_transpose(Lambda, TRUE/*values*/, &c),
	    tmp = M_cholmod_ssmult(Lamtr, Zt, 0/*stype*/,
				   TRUE/*values*/, TRUE/*sorted*/, &c);
	int *it = (int*)(tmp->i), *iu = (int*)(Ut->i),
	    *pt = (int*)(tmp->p), *pu = (int*)(Ut->p);
	int nnz = pu[n];
	double *xt = (double*)(tmp->x), *xu = (double*)(Ut->x); 

	M_cholmod_free_sparse(&Lamtr, &c);
	for (int j = 0; j <= n; j++)
	    if (pt[j] != pu[j])
		error(_("Ut is not consistent with Lambda %*% Zt"));
	for (int j = 0; j < nnz; j++) {
	    if (it[j] != iu[j])
		error(_("Ut is not consistent with Lambda %*% Zt"));
	    xu[j] = xt[j];
	}
	M_cholmod_free_sparse(&tmp, &c);
    }
}

void merenv::update_eta() {
    CHM_DN ceta = N_AS_CHM_DN(eta, n, 1), cu = N_AS_CHM_DN(u, q, 1);
    double one[2] = {1,0};
    if (offset)			// initialize to offset if used
	dble_cpy(eta, offset, n);
    else dble_zero(eta, n);	// otherwise to zero
    if (sparseX) {
	CHM_DN cbeta = N_AS_CHM_DN(beta, p, 1);
	M_cholmod_sdmult(sX, 0/*transpose*/, one, one, cbeta, ceta, &c);
    } else
	F77_CALL(dgemv)("N", &n, &p, one, X, &n, beta, &i1, one, eta, &i1);
    M_cholmod_sdmult(Ut, 1/*transpose*/, one, one, cu, ceta, &c);
}

CHM_DN merenv::crossprod_Lambda(CHM_DN rhs, CHM_DN ans) {
    int nc = rhs->ncol;
    double one[2] = {1,0}, zero[2] = {0,0};
    
    if (((int)(rhs->nrow)) != q)
	error(_("in crossprod_Lambda, rhs->ncol = %d should be %d"),
	      rhs->ncol, q);
    if (dLam) {
	double *ax = (double*)(ans->x), *rx = (double*)(rhs->x);
	for (int j = 0; j < nc * q; j++) ax[j] = rx[j] * Lambdax[j % q];
    } else
	M_cholmod_sdmult(Lambda, 1/*transpose*/, one, zero, rhs, ans, &c);
    return ans;
}

/**
 * Update to new value of theta and evaluate the deviance or REML
 * criterion.
 *
 * @param thnew pointer to an numeric vector theta
 *
 * @return deviance or REML criterion according to the value of REML
 */
double merenv::update_dev(SEXP thnew) {
    double mone[2] = {-1,0}, one[2] = {1,0};
    CHM_DN cZty = N_AS_CHM_DN(Zty, q, 1);
    int info;

    if (!isReal(thnew) || LENGTH(thnew) != nth)
	error(_("theta must be numeric and of length %d"), nth);
    dble_cpy(theta, REAL(thnew), nth);
    update_Lambda();
    update_Ut();
    M_cholmod_factorize_p(Ut, one, (int*)NULL, (size_t)0, L, &c);
    *ldL2 = M_chm_factor_ldetL2(L);
    if (sparseX) {
	error(_("code not yet written"));
    } else {
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
				// downdate and factor XtX, solve for beta
	dble_cpy(RX, XtX, p * p);
	F77_CALL(dsyrk)("U", "T", &p, &q, mone, RZX, &q, one, RX, &p);
	dble_cpy(beta, Xty, p);
	F77_CALL(dgemv)("T", &q, &p, mone, RZX, &q, (double*)(cu->x),
			&i1, one, beta, &i1);
	F77_CALL(dposv)("U", &p, &i1, RX, &p, beta, &p, &info);
	if (info)
	    error(_("Downdated X'X is not positive definite, %d."), info);
				// evaluate ldRX2
	*ldRX2 = 0;
	for (int i = 0; i < p; i++) *ldRX2 += 2 * log(RX[i * (p + 1)]);
				// solve for u
	F77_CALL(dgemv)("N", &q, &p, mone, RZX, &q, beta, &i1, one,
			(double*)(cu->x), &i1);
	tmp1 = M_cholmod_solve(CHOLMOD_Lt, L, cu, &c);
	M_cholmod_free_dense(&cu, &c);
	tmp2 = M_cholmod_solve(CHOLMOD_Pt, L, tmp1, &c);
	M_cholmod_free_dense(&tmp1, &c);
	dble_cpy(u, (double*)(tmp2->x), q);
	M_cholmod_free_dense(&tmp2, &c);
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

int merenv::validate() {
    return 1;			// checking is done in constructor
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
    return ScalarReal(merenv(rho).update_dev(thnew));
}

/**
 * Check validity of an merenv environment
 *
 * @param x pointer to an merenv environment
 */
SEXP lmerenv_validate(SEXP rho) {
    return ScalarLogical(merenv(rho).validate());
}
    
