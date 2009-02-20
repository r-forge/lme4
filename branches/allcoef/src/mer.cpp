#include "mer.h"

#include "Matrix.h"
#include "lme4utils.hpp"

class mer {
public:
    mer(SEXP rho);		//< instantiate from an environment
    ~mer(){
	delete L; delete A; 
	d[wrss_POS] = wrss;
	d[usqr_POS] = usqr;
	d[pwrss_POS] = pwrss;
	d[sigmaML_POS] = sigmaML;
	d[sigmaREML_POS] = sigmaREML;
	d[ldL2_POS] = ldL2;
	d[ldRX2_POS] = ldRX2;
    }
/**
 * Create a sparse matrix U from A, etaGamma, muEta, srwt and Perm
 * @return a freshly allocated q by n sparse matrix
 */
    CHM_SP A_to_U();
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
    double *RX, *RZX, *V, *X, *beta0, *d, *eta, *fixef, *etaGamma,
	*mu, *muEta, *nvec, *offset, *pWt, *srwt, *res, *u, *u0, *var, *y,
	*ghx, *ghw, ldL2, ldRX2, pwrss, sigmaML, sigmaREML,
	usqr, wrss;
    SEXP flistP, nlmodel, pnames, nlenv;
    CHM_FR L;
    CHM_SP A;

    void extractL(SEXP rho);
/**
 * Fill in the V matrix using X, etaGamma, muEta and srwt.
 */
    void X_to_V();
    void eval_nonlin(const double *tmp);
    void eval_muEta();
    void eval_varFunc();
    double* eval_m2lcond(double *ans, const int *Grps);
    double* eval_devResid(double *ans, const int *Grps);
};

double mer::mone = -1;		// These don't really change but they
double mer::one = 1;		// are passed to Fortran code as
double mer::zero = 0;		// pointers. FIXME: Add const to those
				// declarations.
int mer::i1 = 1;

const double mer::CM_TOL = 1e-12;
const double mer::CM_SMIN = 1e-5;
const double mer::LTHRESH = 30;
const double mer::MLTHRESH = -30;
const double mer::MPTHRESH = qnorm5(DOUBLE_EPS, 0, 1, 1, 0);
const double mer::PTHRESH = -MPTHRESH;
const double mer::INVEPS = 1/DOUBLE_EPS;
    
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
inline double sqr_length(const double *x, int nn)
{
    double ans = 0;
    for (int i = 0; i < nn; i++)
	ans += x[i] * x[i];
    return ans;
}

/**
 * Evaluate y * log(y/mu) with the correct limiting value at y = 0.
 *
 * @param y 
 * @param mu
 *
 * @return y * log(y/mu) for y > 0, 0 for y == 0.
 */
inline double y_log_y(double y, double mu)
{
    return (y) ? (y * log(y/mu)) : 0;
}

double *VAR_REAL_NULL(SEXP rho, SEXP nm, int len, int nullOK, int absentOK)
{
    const char *pn = CHAR(PRINTNAME(nm));
    SEXP var = findVarInFrame(rho, nm);
    if (var == R_UnboundValue) {
	if (absentOK) {
	    defineVar(nm, allocVector(REALSXP, len), rho);
	    return(REAL(findVarInFrame(rho, nm)));
	} else error(_("object named '%s' not found in environment"), pn);
    }
    int ll = LENGTH(var);
    
    if (!ll) {
	if (nullOK) return (double*) NULL;
	error(_("numeric object '%s' may not have length 0"), pn);
    }
    if (len && ll != len)
	error(_("Expected numeric object '%s' to be length %d, got %d"),
	      pn, len, ll);
    if (!isReal(var))
	error(_("numeric object '%s' not found in env"), pn);
    return REAL(var);
}

double *VAR_REAL_NULL(SEXP rho, SEXP nm, int len, int nullOK)
{
    return VAR_REAL_NULL(rho, nm, len, nullOK, FALSE);
}

double *VAR_REAL_NULL(SEXP rho, SEXP nm, int len)
{
    return VAR_REAL_NULL(rho, nm, len, FALSE, FALSE);
}

// Definition of methods for the mer class

void mer::extractL(SEXP rho)
{
    SEXP LL = findVarInFrame(rho, lme4_LSym);
    if (LL == R_UnboundValue) {	// need to create L
	double d1[2] = {1,0};

	CHM_SP U = A_to_U();
	L = M_cholmod_analyze(U, &c);
	Memcpy(perm, (int*)(L->Perm), q);
	M_cholmod_free_sparse(&U, &c);
	M_cholmod_free_factor(&L, &c);
				//set the natural ordering
	int nm = c.nmethods, ord0 = c.method[0].ordering, posto = c.postorder;
	c.nmethods = 1; c.method[0].ordering = CHOLMOD_NATURAL;
	c.postorder = FALSE;

	U = A_to_U();
	L = M_cholmod_analyze(U, &c);
				// force a numeric factorization
	if (!M_cholmod_factorize_p(U, d1, (int*)NULL, 0 /*fsize*/, L, &c))
	    error(_("cholmod_factorize_p failed: status %d, minor %d from ncol %d"),
		  c.status, L->minor, L->n);
	M_cholmod_free_sparse(&U, &c);
				//restore the former settings
	c.nmethods = nm; c.method[0].ordering = ord0; c.postorder = posto;
				//define object L in rho
	LL = M_chm_factor_to_SEXP(L, 0);
	defineVar(lme4_LSym, LL, rho);
	M_cholmod_free_factor(&L, &c);
	LL = findVarInFrame(rho, lme4_LSym);
    } 
    L = new cholmod_factor;
    M_as_cholmod_factor(L, LL);
    if (!(L->is_ll))	   // should never happen, but check anyway
	error(_("L must be LL', not LDL'"));
}

mer::mer(SEXP rho)
{
    // Extract slots that must have positive length.
    // Get dimensions of the problem
    SEXP sl = findVarInFrame(rho, lme4_ySym);
    if (sl == R_UnboundValue)
	error(_("Environment must contain response y"));
    if (!(n = LENGTH(sl)) || !isReal(sl)) // n = length of response
	error(_("Response vector y must be numeric (double)"));
    y = REAL(sl);
    sl = findVarInFrame(rho, lme4_fixefSym);
    if (sl == R_UnboundValue)
	error(_("Environment must contain numeric vector fixef"));
    if (!(p = LENGTH(sl)) || !isReal(sl)) // p = length of fixef
	error(_("Fixef vector must be numeric (double)"));
    fixef = REAL(sl);
    sl = findVarInFrame(rho, lme4_uSym);
    if (sl == R_UnboundValue)
	error(_("Environment must contain numeric vector u"));
    if (!(q = LENGTH(sl)) || !isReal(sl)) // q = length of u
	error(_("u vector must be numeric (double)"));
    u = REAL(sl);
    sl = findVarInFrame(rho, lme4_etaGammaSym);
    if (sl == R_UnboundValue)
	error(_("Environment must contain numeric matrix etaGamma"));
    s = LENGTH(sl);		// this may be zero
    if (s % n)
	error(_("length(etaGamma) = %d is not a multiple of length(y) = %d"),
	      s, n);
    s = s / n;
    if (s) {
	if (!isReal(sl) || !isMatrix(sl))
	    error(_("Non-null etaGamma must be a numeric matrix"));
	etaGamma = REAL(sl);
	pnames = VECTOR_ELT(getAttrib(sl, R_DimNamesSymbol), 1);
	if (!isString(pnames) || LENGTH(pnames) != s)
	    error(_("Slot V must be a matrix with %d named columns"), s);
    } else etaGamma = (double*) NULL;
    
    N = (s ? s * n : n);	// number of rows in X (cols in A)

    A = new cholmod_sparse;
    M_as_cholmod_sparse(A, findVarInFrame(rho, lme4_ASym),
			(Rboolean)TRUE, (Rboolean)FALSE);
    if (((int)(A->ncol)) != N)
	error(_("Number of columns of A = %d is not (n = %d)*(s = %d)"),
	      A->ncol, n, s);
    sl = findVarInFrame(rho, lme4_permSym);
    if (!isInteger(sl) || LENGTH(sl) != q)
	error(_("\"perm\" must be an integer vector of length %d"), q);
    perm = INTEGER(sl);
    extractL(rho);

    RX = VAR_REAL_NULL(rho, lme4_RXSym, p * p);
    RZX = VAR_REAL_NULL(rho, lme4_RZXSym, q * p);
    X = VAR_REAL_NULL(rho, lme4_XSym, N * p);
    beta0 = VAR_REAL_NULL(rho, lme4_startSym, p, TRUE);
    d = VAR_REAL_NULL(rho, lme4_devianceSym, NULLdev_POS + 1);
    dims = INTEGER(findVarInFrame(rho, lme4_dimsSym));
    eta = VAR_REAL_NULL(rho, lme4_etaSym, n, FALSE, TRUE);
    mu = VAR_REAL_NULL(rho, lme4_muSym, n, FALSE, TRUE);
    offset = VAR_REAL_NULL(rho, lme4_offsetSym, n, TRUE);
    pWt = VAR_REAL_NULL(rho, lme4_weightsSym, n, TRUE);
    res = VAR_REAL_NULL(rho, lme4_residSym, n, FALSE, TRUE);
    int ncv = dims[vTyp_POS] != 1; // non-constant variance function
    var = VAR_REAL_NULL(rho, lme4_varSym, n, !ncv, ncv);
    int nidl = dims[lTyp_POS] != 5; // non-identity link
    muEta = VAR_REAL_NULL(rho, lme4_muEtaSym, n, !nidl, nidl);
    u0 = VAR_REAL_NULL(rho, install("u0"), q, TRUE);
    nvec = (double*)NULL;
    sl = findVarInFrame(rho, install("n"));
    if (sl != R_UnboundValue) {
	if (!isReal(sl)) error(_("Vector n must be numeric (double)"));
	if (LENGTH(sl) != n) error(_("length(n) = %d != length(y) = %d"),
				   LENGTH(sl), n);
	nvec = REAL(sl);
    }
    int reswt = pWt || var;	// force non-null srwt if TRUE
    srwt = VAR_REAL_NULL(rho, lme4_sqrtrWtSym, n, !reswt, reswt);
    nlmodel = findVarInFrame(rho, lme4_nlmodelSym);
    nlenv = findVarInFrame(rho, lme4_nlenvSym);
}

CHM_SP mer::A_to_U()
{
    CHM_TR At = M_cholmod_sparse_to_triplet(A, &c);
    int *Ati = (int*)(At->i), *Atj = (int*)(At->j), nz = At->nnz,
	*iperm = new int[q];
    double *Atx = (double*)(At->x);
	
    for (int j = 0; j < q; j++) iperm[perm[j]] = j;
    for (int p = 0; p < nz; p++) {
	int j = Atj[p], jj;
	    
	Atj[p] = jj = j % n;
	Atx[p] *= (etaGamma ? etaGamma[j] : 1) *
	    (muEta ? muEta[jj] : 1) * (srwt ? srwt[jj] : 1);
	Ati[p] = iperm[Ati[p]];
    }
    At->ncol = n;
    CHM_SP U = M_cholmod_triplet_to_sparse(At, nz, &c);
    M_cholmod_free_triplet(&At, &c);
    delete[] iperm;
    return U;
}

void mer::eval_nonlin(const double *tmp2)
{
    for (int i = 0; i < s; i++) { // par. vals. into environment
	SEXP vv = findVarInFrame(nlenv, install(CHAR(STRING_ELT(pnames, i))));
	if (!isReal(vv) || LENGTH(vv) != n)
	    error(_("Parameter %s in the environment must be a length %d numeric vector"),
		  CHAR(STRING_ELT(pnames, i)), n);
	dble_cpy(REAL(vv), tmp2 + i * n, n);
    }
    SEXP vv = PROTECT(eval(nlmodel, nlenv));
    if (!isReal(vv) || LENGTH(vv) != n)
	error(_("evaluated model is not a numeric vector of length %d"), n);
    SEXP gg = getAttrib(vv, lme4_gradientSym);
    if (!isReal(gg) || !isMatrix(gg))
	error(_("gradient attribute of evaluated model must be a numeric matrix"));
    int *gdims = INTEGER(getAttrib(gg, R_DimSymbol));
    if (gdims[0] != n ||gdims[1] != s)
	error(_("gradient matrix must be of size %d by %d"), n, s);
    // colnames of the gradient corresponding to the order of the
    // pnames has been checked
    dble_cpy(eta, REAL(vv), n);
    dble_cpy(etaGamma, REAL(gg), n * s);
    UNPROTECT(1);
}

void mer::X_to_V()
{
    dble_zero(V, n * p);	// zero the array
    for (int j = 0; j < p; j++) {
	for (int i = 0; i < N; i++) {
	    int ii = i % n;
	    V[ii + j * n] +=
		X[i + j * N] * (etaGamma ? etaGamma[i] : 1) *
		(muEta ? muEta[ii] : 1) * (srwt ? srwt[ii] : 1);
	}
    }
}

/**
 * Update the eta, etaGamma, mu, muEta, res and var members from the
 * current values of the fixef and u.  Also evaluate d[wrss_POS] using
 * the current sqrtrWt slot.  The sqrtrWt slot is changed in update_L.
 *
 * @return penalized, weighted residual sum of squares
 */
double mer::update_mu()
{
    double *tmp = new double[N], *pu = new double[q], d1[2] = {1,0};
    CHM_DN ctmp = N_AS_CHM_DN(tmp, N, 1), cpu = N_AS_CHM_DN(pu, q, 1);

				// tmp := offset or tmp := 0
    for (int i = 0; i < N; i++) tmp[i] = offset ? offset[i] : 0;
				// tmp := tmp + X beta
    F77_CALL(dgemv)("N", &N, &p, d1, X, &N, fixef, &i1, d1, tmp, &i1);
				// tmp := tmp + A'P'u
    for (int j = 0; j < q; j++) pu[perm[j]] = u[j];
    if (!M_cholmod_sdmult(A, 1 /* trans */, d1, d1, cpu, ctmp, &c))
	error(_("cholmod_sdmult error returned"));
    delete[] pu;
				// fill in eta
    if (etaGamma) eval_nonlin(tmp); else dble_cpy(eta, tmp, n);
    delete[] tmp;
				// inverse link
    if (muEta) {eval_muEta(); eval_varFunc();} else dble_cpy(mu, eta, n);
    
    wrss = 0;		
    for (int i = 0; i < n; i++) { // update res and wrss
	double wtres = (res[i] = y[i] - mu[i]) * (srwt ? srwt[i] : 1);
	wrss += wtres * wtres;
    }
    usqr = sqr_length(u, q);
    return pwrss = wrss + usqr;
}

double mer::PIRLS()
{
    int cvg, info, verb = dims[verb_POS];
    double *betaold = new double[p],
	*cbeta = new double[p], *tmp = new double[q], *uold = new double[q],
	*wtres = new double[n], cfac = ((double)n)/((double)(q+p)),
	crit, pwrss_old, step, d1[2] = {1,0}, d0[2] = {0,0};
    CHM_DN SOL, cRZX = N_AS_CHM_DN(RZX, q, p), cV,
	cwtres = N_AS_CHM_DN(wtres, n, 1), ctmp = N_AS_CHM_DN(tmp, q, 1);
    CHM_SP U;
    R_CheckStack();

    if (verb < 0) Rprintf("cfac = %g\n", cfac);

    // reset u and fixef to initial values.  This can result in more
    // iterations but it gives a repeatable function evaluation for
    // the optimizer.
    dble_zero(u, q); 
    if (u0) Memcpy(u, u0, q);
    dble_zero(fixef, p);
    if (beta0) Memcpy(fixef, beta0, p);

    V = new double[n * p];
    cV = N_AS_CHM_DN(V, n, p);

    cvg = FALSE;
    update_mu();
    for (int i = 0; i < CM_MAXITER; i++) {
	dble_cpy(uold, u, q);	// record current coefficients
	dble_cpy(betaold, fixef, p);

	if (srwt) {	  // Update the weights and weighted residuals
	    for (int j = 0; j < n; j++)
		wtres[j] = res[j] *
		    (srwt[j] = sqrt((pWt? pWt[j] : 1) * (var ? var[j] : 1)));
	} else dble_cpy(wtres, res, n);
	pwrss_old = sqr_length(wtres, n) + sqr_length(u, q);			 
	U = A_to_U();		// create U
	if (!M_cholmod_factorize_p(U, d1, (int*)NULL, 0 /*fsize*/, L, &c))
	    error(_("cholmod_factorize_p failed: status %d, minor %d from ncol %d"),
		  c.status, L->minor, L->n);
	X_to_V();		   // update V
	if (!M_cholmod_sdmult(U, 0/*no transpose*/, d1, d0, cV, cRZX, &c))
	    error(_("cholmod_sdmult failed: status %d"), c.status);
	if (!(SOL = M_cholmod_solve(CHOLMOD_L, L, cRZX, &c)))
	    error(_("cholmod_solve (CHOLMOD_L) failed: status %d"), c.status);
	dble_cpy(RZX, (double*)(SOL->x), q * p);
	M_cholmod_free_dense(&SOL, &c);
				// solve for RX in downdated V'V
	F77_CALL(dsyrk)("U", "T", &p, &n, d1, V, &n, d0, RX, &p); //V'V
	F77_CALL(dsyrk)("U", "T", &p, &q, &mone, RZX, &q, d1, RX, &p);
	F77_CALL(dpotrf)("U", &p, RX, &p, &info);
	if (info)
	    error(_("Downdated V'V is not positive definite, %d."), info);
				// tmp := U %*% wtdResid 
	M_cholmod_sdmult(U, 0 /* notrans */, d1, d0, cwtres, ctmp, &c);
	for (int j = 0; j < q; j++) tmp[j] -= u[j]; // tmp := tmp - u
	M_cholmod_free_sparse(&U, &c);
				// solve L %*% tmp = tmp
	if (!(SOL = M_cholmod_solve(CHOLMOD_L, L, ctmp, &c)))
	    error(_("cholmod_solve (CHOLMOD_L) failed"));
	dble_cpy(tmp, (double*)(SOL->x), q);
	M_cholmod_free_dense(&SOL, &c);
				// solve RX'cbeta = V'wtres - RZX'cu
	F77_CALL(dgemv)("T", &n, &p, d1, V, &n, wtres, &i1,
			d0, cbeta, &i1);
	F77_CALL(dgemv)("T", &q, &p, &mone, RZX, &q, tmp, &i1,
			d1, cbeta, &i1);
	F77_CALL(dtrsv)("U", "T", "N", &p, RX, &p, cbeta, &i1);
				// evaluate convergence criterion
	double cul2 = sqr_length(tmp, q), cbetal2 = sqr_length(cbeta, p);
	crit = cfac * (cul2 + cbetal2)/ pwrss_old;
	if (verb < 0) Rprintf("cul2 = %g, cbetal2 = %g\n", cul2, cbetal2);
	if (crit < CM_TOL) {	// don't do needless evaluations 
	    cvg = TRUE;
	    break;
	}
				// solve for delta-beta
	F77_CALL(dtrsv)("U", "N", "N", &p, RX, &p, cbeta, &i1);
				// solve t(L) %*% SOL = tmp - RZX cbeta
	F77_CALL(dgemv)("N", &q, &p, &mone, RZX, &q, cbeta, &i1,
			d1, tmp, &i1);
	if (!(SOL = M_cholmod_solve(CHOLMOD_Lt, L, ctmp, &c)))
	    error(_("cholmod_solve (CHOLMOD_Lt) failed"));
	dble_cpy(tmp, (double*)(SOL->x), q);
	M_cholmod_free_dense(&SOL, &c);
	
	for (step = 1; step > CM_SMIN; step /= 2) { // step halving 
	    for (int j = 0; j < q; j++)
		u[j] = uold[j] + step * tmp[j];
	    for (int j = 0; j < p; j++)
		fixef[j] = betaold[j] + step * cbeta[j];
	    pwrss = update_mu();
	    if (verb < 0)
		Rprintf("%2d,%8.6f,%12.4g: %15.6g %15.6g %15.6g %15.6g %15.6g\n",
			i, step, crit, pwrss, pwrss_old, fixef[0], u[0], u[1]);
	    if (pwrss < pwrss_old) {
		pwrss_old = pwrss;
		break;
	    }
	}
	if (step <= CM_SMIN) break;
	if (!(muEta || etaGamma)) { // Gaussian linear mixed models require
	    cvg = TRUE;		    // only 1 iteration
	    break;
	}
    }
    delete[] V;
    delete[] betaold;
    delete[] cbeta;
    delete[] tmp;
    delete[] uold;
    delete[] wtres;
    ldRX2 = 0;
    for (int j = 0; j < p; j++) ldRX2 += 2 * log(RX[j * (p + 1)]);
    ldL2 = M_chm_factor_ldetL2(L);
    sigmaML = sqrt(pwrss/(srwt ? sqr_length(srwt, n) : (double) n));
    sigmaREML = (etaGamma || muEta) ? NA_REAL :
	sigmaML * sqrt((((double) n)/((double)(n - p))));

    return update_dev();
}

/**
 * Evaluate the sum of the deviance residuals for a GLM family
 *
 * @param ans pointer to vector of partial sums
 * @param fac indices associating observations with partial sums 
 *            (may be (int*)NULL)
 * @return ans
 */
double* mer::eval_devResid(double *ans, const int *fac)
{
    const int vTyp = dims[vTyp_POS];
    for (int i = 0; i < n; i++) {
	double mui = mu[i], wi = pWt ? pWt[i] : 1, yi = y[i];
	double ri = yi - mui;
	int ai = fac ? (fac[i] - 1) : 0;
	switch(vTyp) {
	case 1:			/* constant variance */
	    ans[ai] += wi * ri * ri;
	    break;
	case 2:			/* mu(1-mu) variance */
	    ans[ai] += 2 * wi *
		(y_log_y(yi, mui) + y_log_y(1 - yi, 1 - mui));
	    break;
	case 3:			/* mu variance */
	    ans[ai] += 2 * wi * (y_log_y(yi, mui) - (yi - mui));
	    break;
	case 4:			/* mu^2 variance */
	    ans[ai] += 2 * wi * (y_log_y(yi, mui) - (yi - mui)/mui);
	    break;
	case 5:			/* mu^3 variance */
	    ans[ai] += wi * (ri * ri)/(yi * mui * mui);
	    break;
	default:
	    error(_("Unknown vTyp value %d"), vTyp);
	}
    }
    return ans;
}

/**
 * Evaluate the inverse link and the derivative d mu/d eta for the GLM
 * link function of type dims[lTyp_POS]
 */
void mer::eval_muEta()
{
    int lTyp = dims[lTyp_POS];
    for (int i = 0; i < n; i++) {
	double etai = eta[i], tmp, t2;
	switch(lTyp) {
	case 1:		// logit
	    tmp = (etai < MLTHRESH) ? DOUBLE_EPS : ((etai > LTHRESH) ?
						    INVEPS : exp(etai));
	    mu[i] = tmp/(1 + tmp);
	    muEta[i] = mu[i] * (1 - mu[i]);
	    break;
	case 2:		// probit 
	    mu[i] = (etai < MPTHRESH) ? DOUBLE_EPS :
	    ((etai > PTHRESH) ? 1 - DOUBLE_EPS : pnorm5(etai, 0, 1, 1, 0));
	    tmp = dnorm4(eta[i], 0, 1, 0);
	    muEta[i] = (tmp < DOUBLE_EPS) ? DOUBLE_EPS : tmp;
	    break;
	case 3:		// cauchit 
	    error(_("cauchit link not yet coded"));
	    break;
	case 4:		// cloglog
	    tmp = (etai < MLTHRESH) ? DOUBLE_EPS : ((etai > LTHRESH) ?
						    INVEPS : exp(etai));
	    t2 = -expm1(-tmp);
	    mu[i] = (t2 < DOUBLE_EPS) ? DOUBLE_EPS : t2;
	    muEta[i] = tmp * exp(-tmp);
	    break;
	case 5:		// identity
	    mu[i] = etai;
	    muEta[i] = 1.;
	    break;
	case 6:		// log 
	    tmp = exp(etai);
	    muEta[i] = mu[i] = (tmp < DOUBLE_EPS) ? DOUBLE_EPS : tmp;
	    break;
	case 7:		// sqrt
	    mu[i] = etai * etai;
	    muEta[i] = 2 * etai;
	    break;
	default:
	    error(_("General form of glmer_linkinv not yet written"));
	}
    }
}

/**
 * Evaluate the GLM variance function of type vTyp given mu
 *
 * vTyp type of variance function: the 1-based index into
 *        c("constant", "mu(1-mu)", "mu", "mu^2", "mu^3")
 *
 */
void mer::eval_varFunc()
{
    int vTyp = dims[vTyp_POS];
    for (int i = 0; i < n; i++) {
	double mui = mu[i];
	switch(vTyp) {
	case 1:			/* constant variance */
	    var[i] = 1.;
	    break;
	case 2:			/* mu(1-mu) variance */
	    if (mui <= 0 || mui >= 1)
		error(_("mu[i] must be in the range (0,1): mu = %g, i = %d"),
		      mu, i);
	    var[i] = mui * (1 - mui);
	    break;
	case 3:			/* mu variance */
	    if (mui <= 0)
		error(_("mu[i] must be positive: mu = %g, i = %d"), mu, i);
	    var[i] = mui;
	    break;
	case 4:			/* mu^2 variance */
	    if (mui <= 0)
		error(_("mu[i] must be positive: mu = %g, i = %d"), mu, i);
	    var[i] = mui * mui;
	    break;
	case 5:			/* mu^3 variance */
	    if (mui <= 0)
		error(_("mu[i] must be positive: mu = %g, i = %d"), mu, i);
	    var[i] = mui * mui * mui;
	    break;
	default:
	    error(_("Unknown vTyp value %d"), vTyp);
	}
    }
}

/**
 * Evaluate -2*log(cond dens)
 *
 * @param ans double precision vector to hold the partial sums
 * @param Grps integer vector of groups.  If (int*)NULL, no groups are used
 *
 * @return ans
 *
 * \note: The ans vector *must* be zeroed before passing to this method.
 *
 */
double* mer::eval_m2lcond(double *ans, const int *Grps)
{
    const int fTyp = dims[fTyp_POS], vTyp = dims[vTyp_POS];
    if ((fTyp == 4 && vTyp != 3) || (fTyp == 1 && vTyp != 2))
	error(_("Variance function is not consistent with poisson or binomial family"));
    for (int i = 0; i < n; i++) {
	double mui = mu[i], ni = nvec[i], wi = pWt ? pWt[i] : 1, yi = y[i];
	int ai = Grps ? (Grps[i] - 1) : 0;
	switch(fTyp) {
	case 1:			// binomial family
	    ans[ai] -= 2 * dbinom(yi, ni, mui, 1);
	case 5:			// Poisson family
	    ans[ai] -= 2 * dpois(yi, mui, 1) * wi;
	    break;
	default:
	    error(_("Unknown fTyp value %d"), fTyp);
	}
    }
    return ans;
    
}

double mer::update_dev()
{
    int nAGQ = dims[nAGQ_POS];
    double dn = (double)n;

    d[ML_POS] = ldL2;

    if (dims[fTyp_POS] == 2 && dims[lTyp_POS] == 5 &&
	dims[vTyp_POS] == 1) {	// Gaussian linear mixed models
	double dnmp = (double)(n - p);

	d[REML_POS] = ldL2 + ldRX2 +
	    dnmp * (1. + log(pwrss) + log(2. * PI / dnmp));
	d[ML_POS] += dn*(1. + log(pwrss) + log(2. * PI / dn));
	return d[dims[isREML_POS] ? REML_POS : ML_POS];
    }

    if (nAGQ == 1) {		// Laplace evaluation
	double ans = 0;
	eval_m2lcond(&ans, (int*) NULL);
	d[disc_POS] = ans;
	d[ML_POS] += d[disc_POS] + d[usqr_POS];
	return d[ML_POS];
    } else {			// Adaptive Gauss-Hermite quadrature
				// Single grouping factor has been checked.
	const int nl = nlevels(VECTOR_ELT(flistP, 0));
	const int nre = q / nl;	// number of random effects per level
	int *fl0 = INTEGER(VECTOR_ELT(flistP, 0)),
	    *pointer = Alloca(nre, int);
	double *uold = new double[q], *tmp = new double[nl],
	    w_pro = 1, z_sum = 0; // values needed in AGQ evaluation 
				// constants used in AGQ evaluation
	const double sigma = (dims[useSc_POS] == 1)?d[sigmaML_POS]:1;
	const double factor = - 0.5 / (sigma * sigma);
	R_CheckStack();
	
	dble_cpy(uold, u, q); // store conditional mode
	int_zero(pointer, nre);
	dble_zero(tmp, nl);
	
	while(pointer[nre - 1] < nAGQ){
	    double *z = new double[q];	  // current abscissas
	    // current penalized residuals in different levels 
	    double *ans = new double[nl];
	    
	    for(int i = 0; i < nre; ++i){ // update abscissas and weights
		for(int j = 0; j < nl; ++j)
		    z[i + j * nre] = ghx[pointer[i]];
		w_pro *= ghw[pointer[i]];
		if (!muEta) z_sum += z[pointer[i]] * z[pointer[i]];
	    }
	    
	    CHM_DN cz = N_AS_CHM_DN(z, q, 1), sol;
	    if(!(sol = M_cholmod_solve(CHOLMOD_L, L, cz, &c)))
		error(_("cholmod_solve(CHOLMOD_L) failed"));
	    dble_cpy(z, (double *)sol->x, q);
	    M_cholmod_free_dense(&sol, &c);
	    
	    for(int i = 0; i < q; ++i) u[i] = uold[i] + sigma * z[i];
	    update_mu();
	    
	    dble_zero(ans, nl);
	    eval_devResid(ans, fl0);
	    
	    for(int i = 0; i < nre; ++i)
		for(int j = 0; j < nl; ++j)
		    ans[j] += u[i + j * nre] * u[i + j * nre];
	    
	    for(int i = 0; i < nl; ++i)
		tmp[i] += exp( factor * ans[i] + z_sum) * w_pro / sqrt(PI);
	    // move pointer to next combination of weights and abbsicas
	    int count = 0;
	    pointer[count]++;
	    while(pointer[count] == nAGQ && count < nre - 1){
		pointer[count] = 0;
		pointer[++count]++;
	    }
	    
	    w_pro = 1;
	    z_sum = 0;
	    
	    delete[] z; delete[] ans;
	}
	
	for(int j = 0; j < nl; ++j) d[ML_POS] -= 2 * log(tmp[j]);
	dble_cpy(u, uold, q);
	update_mu();
	d[ML_POS] += muEta ? 0 : (dn * log(2*PI*d[pwrss_POS]/dn));
	delete[] tmp; delete[] uold;
    }
    return d[ML_POS];
}

/* Externally callable functions */

/**
 * Generate zeros and weights of Hermite polynomial of order N, for
 * the AGQ method.
 *
 * Derived from Fortran code in package 'glmmML'
 *
 * @param N order of the Hermite polynomial
 * @param x zeros of the polynomial, abscissas for AGQ
 * @param w weights used in AGQ
 */

static void internal_ghq(int N, double *x, double *w)
{
    const double GHQ_EPS = 1e-15;
    const int GHQ_MAXIT = 40;
    int NR, IT, I, K, J;
    double Z = 0, HF = 0, HD = 0;
    double Z0, F0, F1, P, FD, Q, WP, GD, R, R1, R2;
    double HN = 1/(double)N;
    double *X = Calloc(N + 1, double), *W = Calloc(N + 1, double);

    for(NR = 1; NR <= N / 2; NR++){
	if(NR == 1)
	    Z = -1.1611 + 1.46 * sqrt((double)N);
	else
	    Z -= HN * (N/2 + 1 - NR);
	for (IT = 0; IT <= GHQ_MAXIT; IT++) {
	    Z0 = Z;
	    F0 = 1.0;
	    F1 = 2.0 * Z;
	    for(K = 2; K <= N; ++K){
		HF = 2.0 * Z * F1 - 2.0 * (double)(K - 1.0) * F0;
		HD = 2.0 * K * F1;
		F0 = F1;
		F1 = HF;
	    }
	    P = 1.0;
	    for(I = 1; I <= NR-1; ++I){
		P *= (Z - X[I]);
	    }
	    FD = HF / P;
	    Q = 0.0;
	    for(I = 1; I <= NR - 1; ++I){
		WP = 1.0;
		for(J = 1; J <= NR - 1; ++J){
		    if(J != I) WP *= ( Z - X[J] );
		}
		Q += WP;
	    }
	    GD = (HD-Q*FD)/P;
	    Z -= (FD/GD);
	    if (fabs((Z - Z0) / Z) < GHQ_EPS) break;
	}

	X[NR] = Z;
	X[N+1-NR] = -Z;
	R=1.0;
	for(K = 1; K <= N; ++K){
	    R *= (2.0 * (double)K );
	}
	W[N+1-NR] = W[NR] = 3.544907701811 * R / (HD*HD);
    }

    if( N % 2 ){
	R1=1.0;
	R2=1.0;
	for(J = 1; J <= N; ++J){
	    R1=2.0*R1*J;
	    if(J>=(N+1)/2) R2 *= J;
	}
	W[N/2+1]=0.88622692545276*R1/(R2*R2);
	X[N/2+1]=0.0;
    }

    dble_cpy(x, X + 1, N);
    dble_cpy(w, W + 1, N);

    if(X) Free(X);
    if(W) Free(W);
}

/**
 * Return zeros and weights of Hermite polynomial of order n as a list
 *
 * @param np pointer to a scalar integer SEXP
 * @return a list with two components, the abscissas and the weights.
 *
 */
SEXP lme4_ghq(SEXP np)
{
    int n = asInteger(np);
    SEXP ans = PROTECT(allocVector(VECSXP, 2));

    if (n < 1) n = 1;
    SET_VECTOR_ELT(ans, 0, allocVector(REALSXP, n ));
    SET_VECTOR_ELT(ans, 1, allocVector(REALSXP, n ));
    
    internal_ghq(n, REAL(VECTOR_ELT(ans, 0)), REAL(VECTOR_ELT(ans, 1)));
    UNPROTECT(1);
    return ans;
}

/**
 * Generate a Markov-chain Monte Carlo sample from an mer object
 *
 * @param x pointer to an merMCMC object
 * @param fm pointer to an mer object
 *
 * @return x with samples filled in
 */
SEXP mer_MCMCsamp(SEXP x, SEXP fm)
{
    return x;			// This is, obviously, a stub.
}

/**
 * Update the deviance vector in GLMMs, NLMMs and GNLMMs
 * If nAGQ > 1, adaptive Gauss-Hermite quadrature is applied.
 *
 * @param x pointer to an mer environment
 *
 * @return deviance value
 */
SEXP mer_update_dev(SEXP x) {return ScalarReal(mer(x).update_dev());}

/**
 * Update the u and fixef slots in an mer object to their conditional
 * modes using PIRLS (penalized, iteratively reweighted least squares).
 *
 * @param x an mer object
 * 
 * @return deviance value
 */
SEXP mer_PIRLS(SEXP x) {return ScalarReal(mer(x).PIRLS());}

/**
 * Create the U matrix from the A matrix
 *
 * @param x an mer environment
 * 
 * @return U
 */
SEXP mer_A_to_U(SEXP x) {
    CHM_SP U = mer(x).A_to_U();
    SEXP ans = CHM_SP2SEXP(U, "dgCMatrix");
    M_cholmod_free_sparse(&U, &c);
    return ans;
}

/**
 * Externally callable update_mu.
 *
 * Update the eta, v, mu, resid and var slots according to the current
 * values of the parameters and u.  Also evaluate d[wrss_POS] using
 * the current contents of sqrtrWt.  The sqrtrWt slot is updated in
 * update_L.
 *
 * @param x pointer to an mer object
 *
 * @return penalized, weighted residual sum of squares
 */
SEXP mer_update_mu(SEXP x){return ScalarReal(mer(x).update_mu());}

/**
 * Check validity of an mer object
 *
 * @param x Pointer to an mer object
 * @return TRUE, the real checking is in mer::mer
 */
SEXP mer_validate(SEXP x)
{
    return ScalarLogical(1);
}

