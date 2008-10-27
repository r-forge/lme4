#include "mer.h"

double mer::mone = -1;		// These don't really change but they
double mer::one = 1;		// are passed to Fortran code as
double mer::zero = 0;		// pointers 
int mer::i1 = 1;

const double mer::CM_TOL = 1e-10;
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
    for (int i = 0; i < nn; i++) dest[i] = src[i];
    return dest;
}

/**
 * Zero the first nn elements of dest
 *
 * @param dest vector
 * @param nn number of elements in dest
 */
inline void dble_zero(double *dest, int nn)
{
    for (int i = 0; i < nn; i++) dest[i] = 0.;
}

inline void int_zero(int *dest, int nn)
{
    for (int i = 0; i < nn; i++) dest[i] = 0.;
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
    for (int i = 0; i < nn; i++) ans += x[i] * x[i];
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

/* Constants */
    
#ifndef BUF_SIZE
/** size of buffer for an error message */
#define BUF_SIZE 127
#endif	

// Definition of methods for the mer class

void mer::extractA(SEXP x)
{
    // The struct in the result of AS_CHM_SP is allocated with
    // alloca and will be freed when it goes out of scope.  The
    // contents of the pointers do not.  We copy the pointers to a
    // struct we will explicitly free in the destructor.

    CHM_SP Acp = AS_CHM_SP(x);
    A = new cholmod_sparse;
    memcpy(A, Acp, sizeof(cholmod_sparse));
}

void mer::extractL(SEXP x)
{
    // The struct in the result of AS_CHM_FR is allocated with
    // alloca and will be freed when it goes out of scope.  The
    // contents of the pointers do not.  We copy the pointers to a
    // struct we will explicitly free in the destructor.

    CHM_FR Lcp = AS_CHM_FR(x);
    if (!(Lcp->is_ll))	   /* should never happen, but check anyway */
	error(_("L must be LL', not LDL'"));
    L = new cholmod_factor;
    memcpy(L, Lcp, sizeof(cholmod_factor));
}

mer::mer(SEXP x)
{
    // Extract slots that must be positive length.
    // Get dimensions of the problem
    SEXP sl = GET_SLOT(x, lme4_ySym);
    n = LENGTH(sl);  	// number of observations
    y = REAL(sl);
    sl = GET_SLOT(x, lme4_fixefSym);
    p = LENGTH(sl);		// number of fixed effects
    fixef = REAL(sl);
    sl = GET_SLOT(x, lme4_uSym);
    q = LENGTH(sl);		// number of random effects
    u = REAL(sl);
    sl = GET_SLOT(x, lme4_etaGammaSym);
    s = LENGTH(sl);
    if (s % n)
	error(_("length(etaGamma) = %d is not a multiple of length(y) = %d"),
	      s, n);
    s = s / n;
    if (s) {
	pnames = VECTOR_ELT(sl, 1);
	if (!isString(pnames) || LENGTH(pnames) != s)
	    error(_("Slot V must be a matrix with %d named columns"), s);
	etaGamma = REAL(sl);
    } else etaGamma = (double*) NULL;
    
    N = (s ? s * n : n);	// number of rows in X (cols in A)

    extractA(GET_SLOT(x, lme4_ASym));
    if (((int)(A->ncol)) != N)
	error(_("Number of columns of A = %d is not (n = %d)*(s = %d)"),
	      A->ncol, N);
    extractL(GET_SLOT(x, lme4_LSym));

    RX = SLOT_REAL_NULL(x, lme4_RXSym);
    RZX = SLOT_REAL_NULL(x, lme4_RZXSym);
    X = SLOT_REAL_NULL(x, lme4_XSym);
    d = REAL(GET_SLOT(x, lme4_devianceSym));
    dims = INTEGER(GET_SLOT(x, lme4_dimsSym));
    eta = REAL(GET_SLOT(x, lme4_etaSym));
    etaGamma = SLOT_REAL_NULL(x, lme4_etaGammaSym);
    mu = REAL(GET_SLOT(x, lme4_muSym));
    muEta = SLOT_REAL_NULL(x, lme4_muEtaSym);
    offset = SLOT_REAL_NULL(x, lme4_offsetSym);
    pWt = SLOT_REAL_NULL(x, lme4_pWtSym);
    perm = INTEGER(GET_SLOT(x, lme4_permSym));
    res = SLOT_REAL_NULL(x, lme4_residSym);
    srwt = SLOT_REAL_NULL(x, lme4_sqrtrWtSym);
    var = SLOT_REAL_NULL(x, lme4_varSym);
    if ((pWt || var) & !srwt)
	error(_("0 = length(srwt) != max(length(pwt) = %d, length(var) = %d)"),
	      LENGTH(GET_SLOT(x, lme4_pWtSym)), LENGTH(GET_SLOT(x, lme4_varSym)));

    nlmodel = GET_SLOT(x, lme4_nlmodelSym);
    rho = GET_SLOT(x, lme4_envSym);
}


CHM_SP mer::A_to_U()
{
    CHM_TR At = M_cholmod_sparse_to_triplet(A, &c);
    int *Ati = (int*)(At->i), *Atj = (int*)(At->j), nz = At->nnz;
    double *Atx = (double*)(At->x);
    CHM_SP U;
	
    for (int p = 0; p < nz; p++) {
	int j = Atj[p], jj;
	    
	Atj[p] = jj = j % n;
	Atx[p] *= (etaGamma ? etaGamma[j] : 1) *
	    (muEta ? muEta[jj] : 1) * (srwt ? srwt[jj] : 1);
	Ati[p] = perm[Ati[p]];
    }
    At->ncol = N;
    U = M_cholmod_triplet_to_sparse(At, nz, &c);
    M_cholmod_free_triplet(&At, &c);
    return U;
}

void mer::eval_nonlin(const double *tmp2)
{
    for (int i = 0; i < s; i++) { /* par. vals. into env. */
	SEXP vv = findVarInFrame(rho, install(CHAR(STRING_ELT(pnames, i))));
	if (!isReal(vv) || LENGTH(vv) != n)
	    error(_("Parameter %s in the environment must be a length %d numeric vector"),
		  CHAR(STRING_ELT(pnames, i)), n);
	dble_cpy(REAL(vv), tmp2 + i * n, n);
    }
    SEXP vv = PROTECT(eval(nlmodel, rho));
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
    dble_cpy(etaGamma, REAL(gg), n * p);
    UNPROTECT(1);
}

double* mer::X_to_V(double *V)
{
    for (int i = 0; i < N*p; i++) V[i] = 0.; // zero the array
    for (int j = 0; j < p; j++) {
	for (int i = 0; i < N; i++) {
	    int ii = i % n;
	    V[ii + j * n] +=
		X[i + j * N] * (etaGamma ? etaGamma[i] : 1) *
		(muEta ? muEta[ii] : 1) * (srwt ? srwt[ii] : 1);
	}
    }
    return V;
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
    double *tmp = new double[N], d1[2] = {1,0};
    CHM_DN ctmp = N_AS_CHM_DN(tmp, N, 1), cu = N_AS_CHM_DN(u, q, 1);

				// tmp := offset or tmp := 0
    for (int i = 0; i < N; i++) tmp[i] = offset ? offset[i] : 0;
				// tmp := tmp + X beta
    F77_CALL(dgemv)("N", &N, &p, d1, X, &N, fixef, &i1, d1, tmp, &i1);
				// tmp := tmp + A' u
    if (!M_cholmod_sdmult(A, 1 /* trans */, d1, d1, cu, ctmp, &c))
	error(_("cholmod_sdmult error returned"));
				// fill in eta
//FIXME: Do the check for etaGamma in eval_nonlin?
    if (etaGamma) eval_nonlin(tmp); else dble_cpy(eta, tmp, n);
    delete[] tmp;
				// inverse link
    if (muEta) {eval_muEta(); eval_varFunc();} else dble_cpy(mu, eta, n);
    
    wrss = 0;		
    for (int i = 0; i < n; i++) { // update res and wrss
	double wtres = (res[i] = y[i] - mu[i]) * (srwt ? srwt[i] : 1);
	wrss += wtres * wtres;
    }
    usqr = sqr_length(u, q);	// store u'u 
    return pwrss = wrss + usqr;
}

double mer::PIRLS()
{
    int cvg, info, verb = dims[verb_POS];
    double *V = new double[n * p], *betaold = new double[p],
	*cbeta = new double[p], *tmp = new double[q], *uold = new double[q],
	*wtres = new double[q], cfac = ((double)n)/((double)(q+p)),
	crit, pwrss_old, step, d1[2] = {1,0}, d0[2] = {0,0};
    CHM_SP U;
    CHM_DN SOL, cRZX = N_AS_CHM_DN(RZX, q, p), cV = N_AS_CHM_DN(V, n, p),
	cwtres = N_AS_CHM_DN(wtres, n, 1), ctmp = N_AS_CHM_DN(tmp, q, 1);
    R_CheckStack();

    dble_zero(u, q); // resetting u to zero at the beginning of each
		     // PIRLS evaluation requires more iterations but
		     // is necessary to obtain a repeatable
		     // evaluation.  If this is not done the
		     // optimization algorithm can take wild steps.
    cvg = FALSE;
    update_mu();
    for (int i = 0; i < CM_MAXITER; i++) {
	dble_cpy(uold, u, q);	// record current coefficients
	dble_cpy(betaold, fixef, p);

	if (srwt) {    /* Update the weights and weighted residuals */
	    for (int j = 0; j < n; j++)
		wtres[j] = res[j] *
		    (srwt[j] = sqrt((pWt ? pWt[j] : 1.0) * (var ? var[j] : 1.0)));
	} else dble_cpy(wtres, res, n);
	pwrss_old = sqr_length(wtres, n) + sqr_length(u, q);			 
	U = A_to_U();		// create U
	X_to_V(V);		// update V
	if (!M_cholmod_factorize_p(U, d1, (int*)NULL, 0 /*fsize*/, L, &c)) // L
	    error(_("cholmod_factorize_p failed: status %d, minor %d from ncol %d"),
		  c.status, L->minor, L->n);
	if (!M_cholmod_sdmult(U, 0/*no transpose*/, d1, d0, cV, cRZX, &c)) // RZX
	    error(_("cholmod_sdmult failed: status %d"), c.status);
	if (!(SOL = M_cholmod_solve(CHOLMOD_L, L, cRZX, &c)))
	    error(_("cholmod_solve (CHOLMOD_P) failed: status %d"), c.status);
	dble_cpy(RZX, (double*)(SOL->x), q * p);
	M_cholmod_free_dense(&SOL, &c);
				// solve for RX in downdated V'V*/
	F77_CALL(dsyrk)("U", "T", &p, &n, d1, V, &n, d0, RX, &p); /* V'V */
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
	/* evaluate convergence criterion */
	crit = cfac * (sqr_length(tmp, q) + sqr_length(cbeta, p))/ pwrss_old;
	if (crit < CM_TOL) { /* don't do needless evaluations */
	    cvg = TRUE;
	    break;
	}
	/* solve for delta-beta */
	F77_CALL(dtrsv)("U", "N", "N", &p, RX, &p, cbeta, &i1);
	/* solve t(L) %*% SOL = tmp - RZX cbeta */
	F77_CALL(dgemv)("N", &q, &p, &mone, RZX, &q, cbeta, &i1,
			d1, tmp, &i1);
	if (!(SOL = M_cholmod_solve(CHOLMOD_Lt, L, ctmp, &c)))
	    error(_("cholmod_solve (CHOLMOD_Lt) failed"));
	dble_cpy(tmp, (double*)(SOL->x), q);
	M_cholmod_free_dense(&SOL, &c);
	
	for (step = 1; step > CM_SMIN; step /= 2) { /* step halving */
	    for (int j = 0; j < q; j++) u[j] = uold[j] + step * tmp[j];
	    for (int j = 0; j < p; j++) fixef[j] = betaold[j] + step * cbeta[j];
	    pwrss = update_mu();
	    if (verb < 0)
		Rprintf("%2d,%8.6f,%12.4g: %15.6g %15.6g %15.6g %15.6g\n",
			i, step, crit, pwrss, pwrss_old, u[1], u[2]);
	    if (pwrss < pwrss_old) {
		pwrss_old = pwrss;
		break;
	    }
	}
	if (step <= CM_SMIN) break;
	if (!(muEta || etaGamma)) { /* linear mixed models require */
	    cvg = TRUE;		    /* only 1 iteration */
	    break;
	}
    }
    delete[] V; delete[] betaold; delete[] cbeta; delete[] tmp;
    delete[] uold; delete[] wtres;

    return update_dev();
}



/* Stand-alone utility functions (sorted by function name) */

/**
 * Check that slot sym of object x is a numeric matrix of dimension nr
 * by nc.
 *
 * @param buf character buffer of length nb + 1
 * @param nb number of writable positions in the buffer
 * @param x pointer to an mer object
 * @param sym name (symbol, actually) of the slot to check
 * @param nr expected number of rows
 * @param nc expected number of columns
 *
 * @return 0 for success, number of characters written to buf on failure
 */
static int chkDims(char *buf, int nb, SEXP x, SEXP sym, int nr, int nc)
{
    SEXP MP = GET_SLOT(x, sym);
    int *dm = isMatrix(MP) ? INTEGER(getAttrib(MP, R_DimSymbol)) : (int*) NULL;
    if (!dm || !isReal(MP) || dm[0] != nr || dm[1] != nc)
	return snprintf(buf, BUF_SIZE,
			_("Slot %s must be a numeric matrix of size %d by %d"),
			CHAR(PRINTNAME(sym)), nr, nc);
    return 0;
}

/** Check that the length of the sym slot in x is len or, possibly, zero.
 *
 * @param buf character buffer of length nb + 1
 * @param nb number of writable positions in the buffer
 * @param x pointer to an mer object
 * @param sym name (symbol, actually) of the slot to check
 * @param len expected length
 * @param zerok is a length of zero allowed?
 *
 * @return 0 for success, number of characters written to buf on failure

 */
static int chkLen(char *buf, int nb, SEXP x, SEXP sym, int len, int zerok)
{
    int ll;

    if (!(ll = LENGTH(GET_SLOT(x, sym))) == len && zerok && !ll)
	return snprintf(buf, BUF_SIZE, _("Slot %s must have length %d."),
			CHAR(PRINTNAME(sym)), len);
    return 0;
}

/**
 * Evaluate the sum of the deviance residuals for a GLM family
 *
 * @param ans pointer to vector of partial sums
 * @param fac indices associating observations with partial sums 
 *            (may be (int*)NULL)
  * @return the sum of the deviance residuals
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
 *
 * @param mu pointer to the mu vector
 * @param muEta pointer to the muEta vector
 * @param eta pointer to the eta vector
 * @param n length of mu, muEta and eta
 * @param lTyp type of link: the 1-based index into
 *        c("logit", "probit", "cauchit", "cloglog", "identity",
 *          "log", "sqrt", "1/mu^2", "inverse")
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

double mer::update_dev()
{
    int nAGQ = dims[nAGQ_POS];
    double dn = (double) n;

    d[ML_POS] = d[ldL2_POS];

    if (nAGQ == 1) {		// Laplace evaluation
	double ans = 0;
	eval_devResid(&ans, (int*) NULL);
	d[disc_POS] = ans;
//FIXME: Need to take the definition of the AIC function from the
//R glm family and use it to evaluate the deviance.
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
 * Create and initialize L
 *
 * @param CmP pointer to the model matrix for the orthogonal random
 * effects (transposed)
 *
 * @return L
 */
SEXP mer_create_L(SEXP CmP)
{
    double one[] = {1, 0};
    CHM_SP Cm = AS_CHM_SP(CmP);
    CHM_FR L;
    R_CheckStack();

    L = M_cholmod_analyze(Cm, &c);
    if (!M_cholmod_factorize_p(Cm, one, (int*)NULL, 0 /*fsize*/, L, &c))
	error(_("cholmod_factorize_p failed: status %d, minor %d from ncol %d"),
	      c.status, L->minor, L->n);

    return M_chm_factor_to_SEXP(L, 1);
}

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
    return x; /* This is (obviously) a stub.  The real function is below. */
}

/**
 * Optimize the profiled deviance of an lmer object or the Laplace
 * approximation to the deviance of a nlmer or glmer object.
 *
 * @param x pointer to an mer object
 *
 * @return R_NilValue
 */
SEXP mer_optimize(SEXP x)
{
    return R_NilValue;
}

/**
 * Update the deviance vector in GLMMs, NLMMs and GNLMMs
 * If nAGQ > 1, adaptive Gauss-Hermite quadrature is applied.
 *
 * @param x pointer to an mer object
 *
 * @return R_NilValue
 */
SEXP mer_update_dev(SEXP x)
{
    return ScalarReal(mer(x).update_dev());
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

SEXP merMCMC_VarCorr(SEXP x, SEXP typP)
{
    SEXP ST = GET_SLOT(x, lme4_STSym),
	ncs = GET_SLOT(x, install("nc"));
    int *Gp = Gp_SLOT(x), *Sd = INTEGER(GET_DIM(ST)), *nc = INTEGER(ncs);
    int maxnc = 0, nt = LENGTH(ncs), np = Sd[0], nsamp = Sd[1], pos;
    double *sig = SLOT_REAL_NULL(x, lme4_sigmaSym);
    SEXP ans = PROTECT(allocMatrix(REALSXP, nsamp, np + (sig ? 1 : 0)));
    double *av = REAL(ans), *STx = REAL(ST);
    double *as = av + nsamp * np, *t1, *t2, var;
    int *nlev = Alloca(nt, int);
    R_CheckStack();
    
    for (int j = 0; j < nt; j++) {
	nlev[j] = (Gp[j + 1] - Gp[j])/nc[j];
	if (maxnc < nc[j]) maxnc = nc[j];
    }
    if (maxnc > 1) {
	t1 = Alloca(maxnc * maxnc, double);
	t2 = Alloca(maxnc * maxnc, double);
	R_CheckStack();
    }
    
    for (int i = 0; i < nsamp; i++) {
	var = 1; pos = 0;
	if (sig) var = as[i] = sig[i] * sig[i];
	for (int k = 0; k < nt; k++) {
	    if (nc[k] < 2) {
		double sd = STx[pos + i * np] * sig[i];
		av[i + nsamp * pos++] = sd * sd;
	    }
	    else error(_("Code not yet written"));
	}
    }
    
    UNPROTECT(1);
    return ans;
}

/**
 * Check validity of an mer object
 *
 * @param x Pointer to an mer object
 *
 * @return TRUE if the object is a valid mer object, otherwise a string
 *         that describes the violation.
 */
SEXP mer_validate(SEXP x)
{
    SEXP devianceP = GET_SLOT(x, lme4_devianceSym),
	dimsP = GET_SLOT(x, lme4_dimsSym),
	flistP = GET_SLOT(x, lme4_flistSym), asgnP;
    int *dd = INTEGER(dimsP), *asgn;
    const int n = dd[n_POS], nAGQ = dd[nAGQ_POS],
	nt = dd[nt_POS], nfl = LENGTH(flistP),
	p = dd[p_POS], q = dd[q_POS], s = dd[s_POS];
    int nq, nv = n * s;
    CHM_SP Zt = Zt_SLOT(x), A =  A_SLOT(x);
    CHM_FR L = L_SLOT(x);
    char *buf = Alloca(BUF_SIZE + 1, char);
    R_CheckStack();
				/* check lengths */
    asgnP = getAttrib(flistP, install("assign"));
    if (!isInteger(asgnP) || LENGTH(asgnP) != nt)
	return mkString(_("Slot flist must have integer attribute 'assign' of length dims['nt']"));
    asgn = INTEGER(asgnP);
    if (nAGQ < 1)
	return mkString(_("nAGQ must be positive"));
    if ((nAGQ > 1) & (nfl != 1))
	return mkString(_("AGQ method requires a single grouping factor"));

    for (int i = 0; i < nt; i++)
	if (asgn[i] <= 0 || asgn[i] > nfl)
	    return mkString(_("All elements of the assign attribute must be in [1,length(flist)]"));

    if (LENGTH(devianceP) != (NULLdev_POS + 1) ||
	LENGTH(getAttrib(devianceP, R_NamesSymbol)) != (NULLdev_POS + 1))
	return mkString(_("deviance slot not named or incorrect length"));
    if (LENGTH(dimsP) != (cvg_POS + 1) ||
	LENGTH(getAttrib(dimsP, R_NamesSymbol)) != (cvg_POS + 1))
	return mkString(_("dims slot not named or incorrect length"));
    if (((int)(L->n)) != q || !L->is_ll || !L->is_monotonic)
	return mkString(_("Slot L must be a monotonic LL' factorization of size dims['q']"));
    if (((int)(Zt->nrow)) != q || ((int)(Zt->ncol)) != nv)
	return mkString(_("Slot Zt must by dims['q']  by dims['n']*dims['s']"));
    if (((int)(A->nrow)) != q || ((int)(A->ncol)) != nv)
	return mkString(_("Slot A must be dims['q']  by dims['n']*dims['s']"));
    if (chkLen(buf, BUF_SIZE, x, lme4_etaSym, n, 0)) return(mkString(buf));
    if (chkLen(buf, BUF_SIZE, x, lme4_etaGammaSym, nv, 1)) return(mkString(buf));
    if (chkLen(buf, BUF_SIZE, x, lme4_fixefSym, p, 0)) return(mkString(buf));
    if (chkLen(buf, BUF_SIZE, x, lme4_muEtaSym, n, 1)) return(mkString(buf));
    if (chkLen(buf, BUF_SIZE, x, lme4_muSym, n, 0)) return(mkString(buf));
    if (chkLen(buf, BUF_SIZE, x, lme4_offsetSym, n, 1)) return(mkString(buf));
    if (chkLen(buf, BUF_SIZE, x, lme4_pWtSym, n, 1)) return(mkString(buf));
    if (chkLen(buf, BUF_SIZE, x, lme4_ranefSym, q, 0)) return(mkString(buf));
    if (chkLen(buf, BUF_SIZE, x, lme4_residSym, n, 0)) return(mkString(buf));
    if (chkLen(buf, BUF_SIZE, x, lme4_sqrtrWtSym, n, 1)) return(mkString(buf));
    if (chkLen(buf, BUF_SIZE, x, lme4_uSym, q, 0)) return(mkString(buf));
    if (chkLen(buf, BUF_SIZE, x, lme4_varSym, n, 1)) return(mkString(buf));
    if (chkLen(buf, BUF_SIZE, x, lme4_ySym, n, 0)) return(mkString(buf));
    if (chkDims(buf, BUF_SIZE, x, lme4_XSym, nv, p)) return(mkString(buf));
    if (chkDims(buf, BUF_SIZE, x, lme4_RZXSym, q, p)) return(mkString(buf));
    if (chkDims(buf, BUF_SIZE, x, lme4_RXSym, p, p)) return(mkString(buf));

    nq = 0;
    for (int i = 0; i < LENGTH(flistP); i++) {
	SEXP fli = VECTOR_ELT(flistP, i);
	if (!isFactor(fli))
	    return mkString(_("flist must be a list of factors"));
    }
    return ScalarLogical(1);
}

/**
 * Check validity of an merMCMC object
 *
 * @param x Pointer to an merMCMC object
 *
 * @return TRUE if the object is a valid merMCMC object, otherwise a string
 *         that describes the violation.
 */
SEXP merMCMC_validate(SEXP x)
{
    SEXP GpP = GET_SLOT(x, lme4_GpSym),
	devianceP = GET_SLOT(x, lme4_devianceSym),
	dimsP = GET_SLOT(x, lme4_dimsSym);
    int *Gp = INTEGER(GpP), *dd = INTEGER(dimsP);
    int nt = dd[nt_POS], np = dd[np_POS], nsamp = LENGTH(devianceP),
	p = dd[p_POS], q = dd[q_POS];
    char *buf = Alloca(BUF_SIZE + 1, char);
    R_CheckStack();
				/* check lengths */
    if (nsamp <= 0)
	return mkString(_("number of samples must be positive"));
    if (LENGTH(dimsP) != (cvg_POS + 1) ||
	LENGTH(getAttrib(dimsP, R_NamesSymbol)) != (cvg_POS + 1))
	return mkString(_("dims slot not named or incorrect length"));
    if (LENGTH(GpP) != nt + 1)
	return mkString(_("Slot Gp must have length dims['nt'] + 1"));
    if (Gp[0] != 0 || Gp[nt] != q)
	return mkString(_("Gp[1] != 0 or Gp[dims['nt'] + 1] != dims['q']"));

//     if (chkLen(buf, BUF_SIZE, x, lme4_ncSym, nt, 0))
// 	return(mkString(buf));
    if (chkLen(buf, BUF_SIZE, x, lme4_sigmaSym, nsamp, !dd[useSc_POS]))
	return(mkString(buf));
    if (chkDims(buf, BUF_SIZE, x, lme4_STSym, np, nsamp))
	return(mkString(buf));
    if (chkDims(buf, BUF_SIZE, x, lme4_fixefSym, p, nsamp))
	return(mkString(buf));
    if (LENGTH(GET_SLOT(x, lme4_ranefSym)))
	if (chkDims(buf, BUF_SIZE, x, lme4_ranefSym, q, nsamp))
	    return(mkString(buf));
    return ScalarLogical(1);
}

/**
 * Update the u and fixef slots in an mer object
 *
 * Update the u and fixef slots in an mer object to their conditional
 * modes using PIRLS (penalized, iteratively reweighted least squares).
 *
 * @param x an mer object
 * 
 * @return not sure yet
 */
SEXP mer_PIRLS(SEXP x) {return ScalarReal(mer(x).PIRLS());}

#if 0


/**
 * Update the eta, mu, resid and var slots in a sparseRasch object
 * from the current values of the model parameters in the beta slot.
 *
 * @param x pointer to an sparseRasch object
 *
 * @return R_NilValue
 */
SEXP spR_update_mu(SEXP x)
{
    int *dims = DIMS_SLOT(x);
    int n = dims[n_POS];
    double *d = DEV_SLOT(x), *eta = Calloc(n, double), *mu = MU_SLOT(x),
	*offset = OFFSET_SLOT(x), *srwt = SRWT_SLOT(x),
	*res = RESID_SLOT(x), *y = Y_SLOT(x), one[] = {1,0};
    CHM_SP Zt = Zt_SLOT(x);
    CHM_DN cbeta = AS_CHM_DN(GET_SLOT(x, lme4_fixefSym)),
	ceta = N_AS_CHM_DN(eta, n, 1);
    R_CheckStack();
    
    for (int i = 0; i < n; i++) eta[i] = offset ? offset[i] : 0;
    if (!M_cholmod_sdmult(Zt, 1 /* trans */, one, one, cbeta, ceta, &c))
	error(_("cholmod_sdmult error returned"));
    lme4_muEta(mu, MUETA_SLOT(x), eta, n, dims[lTyp_POS]);
    lme4_varFunc(VAR_SLOT(x), mu, n, dims[vTyp_POS]);

    d[wrss_POS] = 0;		/* update resid slot and d[wrss_POS] */
    for (int i = 0; i < n; i++) {
	res[i] = (y[i] - mu[i]) * srwt[i];
	d[wrss_POS] += res[i] * res[i];
    }

    Free(eta);
    return R_NilValue;
}

/**
 * Optimize the log-likelihood for the sparse representation of a
 * Rasch model.
 *
 * @param x pointer to a sparseRasch object
 * @param verbP pointer to indicator of verbose output
 *
 * @return R_NilValue
 */
SEXP spR_optimize(SEXP x, SEXP verbP)
{
    int *zp, m, n, nnz, verb = asInteger(verbP);
    double *Zcp, *Ztx, *d = DEV_SLOT(x), *fixef = FIXEF_SLOT(x),
	*mu = MU_SLOT(x), *muEta = MUETA_SLOT(x),
	*pWt = PWT_SLOT(x), *res = RESID_SLOT(x),
	*srwt = SRWT_SLOT(x), *tmp, *tmp2,
	*var = VAR_SLOT(x), *y = Y_SLOT(x), cfac, crit,
	one[] = {1, 0}, step, wrss_old, zero[] = {0, 0};
    CHM_SP Zt = Zt_SLOT(x);
    CHM_FR L = L_SLOT(x);
    CHM_DN cres = N_AS_CHM_DN(res, Zt->ncol, 1), ctmp, sol;
    R_CheckStack();
    
    zp = (int*)Zt->p;
    m = Zt->nrow;
    n = Zt->ncol;
    nnz = zp[n];
    Zcp = Calloc(nnz, double);
    Ztx = (double*)Zt->x;
    tmp = Calloc(m, double);
    tmp2 = Calloc(m, double);
    ctmp = N_AS_CHM_DN(tmp, m, 1);
    cfac = ((double)n) / ((double)m);
    R_CheckStack();

    spR_update_mu(x);
    for (int i = 0; ; i++) {
	d[wrss_POS] = 0;
	for (int j = 0; j < n; j++) {
	    srwt[j] = sqrt((pWt ? pWt[j] : 1) / var[j]);
	    res[j] = srwt[j] * (y[j] - mu[j]);
	    d[wrss_POS] += res[j] * res[j];
	    for (int p = zp[j]; p < zp[j + 1]; p++)
		Zcp[p] = srwt[j] * muEta[j] * Ztx[p];
	}
	Zt->x = (void*)Zcp;
	if (!M_cholmod_factorize_p(Zt, zero, (int*)NULL, 0 /*fsize*/, L, &c))
	    error(_("cholmod_factorize_p failed: status %d, minor %d from ncol %d"),
		  c.status, L->minor, L->n);
	wrss_old = d[wrss_POS];
				/* tmp := Zt %*% muEta %*% var %*% wtdResid */
	M_cholmod_sdmult(Zt, 0 /* notrans */, one, zero, cres, ctmp, &c);
	dble_cpy(tmp2, tmp, m);
	apply_perm(tmp, tmp2, (int*)L->Perm, m);
				/* solve L %*% sol = tmp */
	if (!(sol = M_cholmod_solve(CHOLMOD_L, L, ctmp, &c)))
	    error(_("cholmod_solve (CHOLMOD_L) failed"));
	dble_cpy(tmp, (double*)(sol->x), m);
	M_cholmod_free_dense(&sol, &c);
				/* evaluate convergence criterion */
	crit = cfac * sqr_length(tmp, m) / wrss_old;
	if (crit < CM_TOL) break; /* don't do needless evaluations */
				/* solve t(L) %*% sol = tmp */
	if (!(sol = M_cholmod_solve(CHOLMOD_Lt, L, ctmp, &c)))
	    error(_("cholmod_solve (CHOLMOD_Lt) failed"));
	dble_cpy(tmp, (double*)(sol->x), m);
	M_cholmod_free_dense(&sol, &c);

	dble_cpy(tmp2, fixef, m);
	for (step = 1; step > CM_SMIN; step /= 2) { /* step halving */
	    for (int j = 0; j < m; j++) fixef[j] = tmp2[j] + step * tmp[j];
	    spR_update_mu(x);
	    if (verb < 0)
		Rprintf("%2d,%8.6f,%12.4g: %15.6g %15.6g %15.6g %15.6g\n",
			i, step, crit, d[wrss_POS], wrss_old, fixef[1], fixef[2]);
	    if (d[wrss_POS] < wrss_old) {
		wrss_old = d[wrss_POS];
		break;
	    }
	}
	if (step <= CM_SMIN || i > CM_MAXITER) return 0;
    }
    Free(tmp2);
    Free(Zcp);
    Free(tmp);
    return R_NilValue;
}

/**
 * Update the ST list of arrays and the sparse model matrix A
 *
 * @param x an mer object
 * @param sigma current standard deviation of the per-observation
 *        noise terms.
 * @param vals pointer to memory in which to store the updated values
 *        of the ST parameters
 */
static void MCMC_ST(SEXP x, double sigma, double *vals)
{
    int *Gp = Gp_SLOT(x), *dims = DIMS_SLOT(x), *perm = PERM_VEC(x);
    int nt = dims[nt_POS], q = dims[q_POS], pos = 0;
    double *Ptu = Calloc(q, double), *u = U_SLOT(x);
    double **st = Alloca(nt, double*);
    int *nc = Alloca(nt, int), *nlev = Alloca(nt, int);
    R_CheckStack();

    ST_nc_nlev(GET_SLOT(x, lme4_STSym), Gp, st, nc, nlev);
				/* inverse permutation of u */
    for (int j = 0; j < q; j++) Ptu[perm[j]] = u[j];

    for (int i = 0; i < nt; i++) {
	int qi = nc[i], nl = nlev[i];
	double *sti = st[i], *ui = Ptu + Gp[i];

	if (qi == 1) {
	    *sti *= sqrt(sqr_length(ui, nl)/rchisq((double)nl))/sigma;
	    vals[pos++] = *sti;
	} else {
	    int info, qisq = qi * qi, qip1 = qi + 1;
	    double invsigsq = 1/(sigma * sigma), one[] = {1,0}, zer[] = {0,0};

	    double *scal = Alloca(qisq, double), *tmp = Alloca(qisq, double);
	    R_CheckStack();
				/* scal := crossprod(U_i)/sigma^2 */
	    F77_CALL(dsyrk)("L", "T", &qi, &nl, &invsigsq, ui, &nl,
			    zer, scal, &qi); 
	    F77_CALL(dpotrf)("L", &qi, scal, &qi, &info);
	    if (info)
		error(_("scal matrix %d is not positive definite"), i + 1);
	    /* solve W_i' B_i = R_i for W_i a random std Wishart factor */
	    F77_CALL(dtrsm)("L", "L", "T", "N", &qi, &qi, one,
			    std_rWishart_factor((double)(nl - qi + 1), qi, 0, tmp),
			    &qi, scal, &qi);
				/* B_i := T(theta) %*% S(theta) %*% B_i */
	    for (int k = 0; k < qi; k++)
		for (int j = 0; j < qi; j++)
		    scal[j * qi + k] *= st[i][k * qip1];
	    F77_CALL(dtrmm)("L", "L", "N", "U", &qi, &qi, one, st[i], &qi, scal, &qi);
				/* Create the lower Cholesky factor */
	    F77_CALL(dsyrk)("L", "T", &qi, &qi, one, scal, &qi, zer, tmp, &qi); 
	    F77_CALL(dpotrf)("L", &qi, tmp, &qi, &info);
	    if (info)
		error(_("crossproduct matrix %d is not positive definite"), i + 1);
				/* Update S_i and T_i */
	    for (int j = 0; j < qi; j++) vals[pos++] = sti[j * qip1] = tmp[j * qip1];
	    for (int j = 0; j < qi; j++)
		for (int k = j + 1; k < qi; k++)
		    vals[pos++] = sti[j * qi + k] = tmp[j * qi + k]/sti[j * qip1];
	}
    }
    update_A(x);
    Free(Ptu);
}



/* Functions for sampling from a Wishart distribution */
/* FIXME: Move these to the R sources */

/**
 * Simulate the Cholesky factor of a standardized Wishart variate with
 * dimension p and nu degrees of freedom.
 *
 * @param nu degrees of freedom
 * @param p dimension of the Wishart distribution
 * @param upper if 0 the result is lower triangular, otherwise upper
                triangular
 * @param ans array of size p * p to hold the result
 *
 * @return ans
 */
static double
*std_rWishart_factor(double nu, int p, int upper, double ans[])
{
    int pp1 = p + 1;

    if (nu < (double) p || p <= 0)
	error("inconsistent degrees of freedom and dimension");

    AZERO(ans, p * p);
    for (int j = 0; j < p; j++) {	/* jth column */
	ans[j * pp1] = sqrt(rchisq(nu - (double) j));
	for (int i = 0; i < j; i++) {
	    int uind = i + j * p, /* upper triangle index */
		lind = j + i * p; /* lower triangle index */
	    ans[(upper ? uind : lind)] = norm_rand();
	    ans[(upper ? lind : uind)] = 0;
	}
    }
    return ans;
}

/**
 * Simulate a sample of random matrices from a Wishart distribution
 *
 * @param ns Number of samples to generate
 * @param nuP Degrees of freedom
 * @param scal Positive-definite scale matrix
 *
 * @return
 */
SEXP
lme4_rWishart(SEXP ns, SEXP nuP, SEXP scal)
{
    SEXP ans;
    int *dims = INTEGER(getAttrib(scal, R_DimSymbol)), info,
	n = asInteger(ns), psqr;
    double *scCp, *ansp, *tmp, nu = asReal(nuP), one = 1, zero = 0;

    if (!isMatrix(scal) || !isReal(scal) || dims[0] != dims[1])
	error("scal must be a square, real matrix");
    if (n <= 0) n = 1;
    psqr = dims[0] * dims[0];
    tmp = Alloca(psqr, double);
    scCp = Alloca(psqr, double);
    R_CheckStack();

    dble_cpy(scCp, REAL(scal), psqr);
    AZERO(tmp, psqr);
    F77_CALL(dpotrf)("U", &(dims[0]), scCp, &(dims[0]), &info);
    if (info)
	error("scal matrix is not positive-definite");
    PROTECT(ans = alloc3DArray(REALSXP, dims[0], dims[0], n));
    ansp = REAL(ans);
    GetRNGstate();
    for (int j = 0; j < n; j++) {
	double *ansj = ansp + j * psqr;
	std_rWishart_factor(nu, dims[0], 1, tmp);
	F77_CALL(dtrmm)("R", "U", "N", "N", dims, dims,
			&one, scCp, dims, tmp, dims);
	F77_CALL(dsyrk)("U", "T", &(dims[1]), &(dims[1]),
			&one, tmp, &(dims[1]),
			&zero, ansj, &(dims[1]));

	for (int i = 1; i < dims[0]; i++)
	    for (int k = 0; k < i; k++)
		ansj[i + k * dims[0]] = ansj[k + i * dims[0]];
    }

    PutRNGstate();
    UNPROTECT(1);
    return ans;
}


/**
 * Permute the vector src according to the inverse of perm into dest
 *
 * @param dest destination
 * @param src source
 * @param perm NULL or 0-based permutation of length n
 * @param n length of src, dest and perm
 *
 * @return dest
 *
 * \note If perm is NULL the first n elements of src are copied to dest.
 */
static R_INLINE double*
apply_iperm(double *dest, const double *src, const int *perm, int n)
{
    for (int i = 0; i < n; i++) dest[perm ? perm[i] : i] = src[i];
    return dest;
}

/**
 * Update the L, RZX, RX, sqrtrWt and resid slots.  It is assumed that
 * update_mu has already been called at the current values of u and
 * the model parameters and that A has been updated.
 *
 * @param x pointer to an mer object
 *
 * @return penalized weighted residual sum of squares
 *
 */
static double update_L(SEXP x)
{
    int *dims = DIMS_SLOT(x);
    int n = dims[n_POS], p = dims[p_POS], s = dims[s_POS];
    double *RZX = RZX_SLOT(x), *RX = RX_SLOT(x), *V = etaGamma_SLOT(x),
	*X = X_SLOT(x), *cx = Cx_SLOT(x), *d = DEV_SLOT(x),
	*res = RESID_SLOT(x), *mu = MU_SLOT(x), *muEta = MUETA_SLOT(x),
	*pwt = PWT_SLOT(x), *sXwt = SXWT_SLOT(x), *srwt = SRWT_SLOT(x),
	*var =  VAR_SLOT(x), *y = Y_SLOT(x), one[] = {1,0};
    CHM_SP A = A_SLOT(x);
    CHM_FR L = L_SLOT(x);
    R_CheckStack();
    
    if (var || pwt) {	   /* Update srwt and res. Reevaluate wrss. */
	d[wrss_POS] = 0;
	for (int j = 0; j < n; j++) {
	    srwt[j] = sqrt((pwt ? pwt[j] : 1.0) / (var ? var[j] : 1.0));
	    res[j] = srwt[j] * (y[j] - mu[j]);
	    d[wrss_POS] += res[j] * res[j];
	}
    }
/* Reconsider having separate Cx and C slots.  Change names to U and V? */
    if (sXwt) {			/* Update sXwt and C */
	int *ai = (int*)A->i, *ap = (int*)A->p;
	double *ax = (double*)(A->x), *WX, *wy;
	CHM_SP C = A;

	for (int j = 0; j < s; j++) {
	    for (int i = 0; i < n; i++) {
		int ja = i + j * n;
		sXwt[ja] = (srwt ? srwt[i] : 1) *
		    (muEta ? muEta[i] : 1) * (V ? V[ja] : 1);
	    }
	}
	if (cx) {		/* C is a scaled version of A */
	    for (int j = 0; j < n; j++)
		for (int p = ap[j]; p < ap[j + 1]; p++)
		    cx[p] = ax[p] * sXwt[j];
	    A->x = (void*)cx;
	} else {
	    int *ci, *cp;
	    C = Cm_SLOT(x);
	    R_CheckStack();

	    ci = (int*)C->i; cp = (int*)C->p; cx = (double*)C->x;
	    AZERO(cx, cp[n]);
	    for (int j = 0; j < s; j++)
		for (int i = 0; i < n; i++) {
		    int ja = i + j * n, pc;
		    for (int pa = ap[ja]; pa < ap[ja + 1]; pa++) {
			for (pc = cp[i]; pc < cp[i + 1]; pc++)
			    if (ci[pc] == ai[pa]) break;
			if (pc >= cp[i + 1])
			    error(_("Structure of Cm and A are not consistent"));
			cx[pc] += ax[pa] * sXwt[ja];
		    }
		}
	    A = C;
	}
    		     /* Replace X and y by weighted X and y */
	WX = Calloc(n * p, double);
	wy = Calloc(n, double);
	
	AZERO(WX, n * p);
	for (int i = 0; i < n; i++) {
	    wy[i] = sXwt[i] * y[i];
	    for (int j = 0; j < p; j++)
		WX[i + j * n] += sXwt[i] * X[i + j * n];
	}
	X = WX;
	y = wy;
    }
    if (!M_cholmod_factorize_p(A, one, (int*)NULL, 0 /*fsize*/, L, &c))
	error(_("cholmod_factorize_p failed: status %d, minor %d from ncol %d"),
	      c.status, L->minor, L->n);

    d[ldL2_POS] = M_chm_factor_ldetL2(L);
    d[pwrss_POS] = d[usqr_POS] + d[wrss_POS];
    d[sigmaML_POS] = sqrt(d[pwrss_POS]/
			  (srwt ? sqr_length(srwt, n) : (double) n));
    d[sigmaREML_POS] = (V || muEta) ? NA_REAL :
	d[sigmaML_POS] * sqrt((((double) n)/((double)(n - p))));
    return d[pwrss_POS];
}

/**
 * Update the RZX and RX slots in an mer object. update_L should be
 * called before update_RX
 *
 * @param x pointer to an mer object
 *
 * @return profiled deviance or REML deviance
 */
static double update_RX(SEXP x)
{
    int *dims = DIMS_SLOT(x), info;
    int n = dims[n_POS], p = dims[p_POS], q = dims[q_POS], s = dims[s_POS];
    double *cx = Cx_SLOT(x), *d = DEV_SLOT(x),
	*RZX = RZX_SLOT(x), *RX = RX_SLOT(x), *sXwt = SXWT_SLOT(x),
	*WX = (double*) NULL, *X = X_SLOT(x);
    CHM_SP A = A_SLOT(x);
    CHM_FR L = L_SLOT(x);
    CHM_DN cRZX = N_AS_CHM_DN(RZX, q, p), ans;
    R_CheckStack();

    if (sXwt) {			/* Create W^{1/2}GHX in WX */
	WX = Calloc(n * p, double);

	AZERO(WX, n * p);
	for (int j = 0; j < p; j++)
	    for (int k = 0; k < s; k++)
		for (int i = 0; i < n; i++)
		    WX[i + j * n] +=
			sXwt[i + k * n] * X[i + n * (k + j * s)];
	X = WX;
	/* Replace A by C, either just the x component or the entire matrix */
	if (cx) A->x = (void*)cx;
	else {
	    A = Cm_SLOT(x);
	    R_CheckStack();
	}
    }
				/* solve L %*% RZX = PAW^{1/2}GHX */
    P_sdmult(RZX, (int*)L->Perm, A, X, p); /* right-hand side */
    ans = M_cholmod_solve(CHOLMOD_L, L, cRZX, &c); /* solution */
    dble_cpy(RZX, (double*)(ans->x), q * p);
    M_cholmod_free_dense(&ans, &c);
    				/* downdate X'X and factor  */
    F77_CALL(dsyrk)("U", "T", &p, &n, one, X, &n, zero, RX, &p); /* X'X */
    F77_CALL(dsyrk)("U", "T", &p, &q, mone, RZX, &q, one, RX, &p);
    F77_CALL(dpotrf)("U", &p, RX, &p, &info);
    if (info)
	error(_("Downdated X'X is not positive definite, %d."), info);
				/* accumulate log(det(RX)^2)  */
    d[ldRX2_POS] = 0;
    for (int j = 0; j < p; j++) d[ldRX2_POS] += 2 * log(RX[j * (p + 1)]);

    if (WX) Free(WX);
    return d[ML_POS];
}


/**
 * Update the projections of the response vector onto the column
 * spaces of the random effects and the fixed effects.  This function
 * is needed separately for the one-argument form of the anova function.
 *
 * @param x an mer object
 * @param pb position to store the random-effects projection
 * @param pbeta position to store the fixed-effects projection
 */
static void lmm_update_projection(SEXP x, double *pb, double *pbeta)
{
    int *dims = DIMS_SLOT(x), i1 = 1;
    int n = dims[n_POS], p = dims[p_POS], q = dims[q_POS];
    double *WX = (double*) NULL, *X = X_SLOT(x),
	*d = DEV_SLOT(x), *RZX = RZX_SLOT(x), *RX = RX_SLOT(x),
	*sXwt = SXWT_SLOT(x),
	*wy = (double*)NULL, *y = Y_SLOT(x),
	mone[] = {-1,0}, one[] = {1,0}, zero[] = {0,0};
    CHM_SP A = A_SLOT(x);
    CHM_FR L = L_SLOT(x);
    CHM_DN cpb = N_AS_CHM_DN(pb, q, 1), sol;
    R_CheckStack();
	
    if (sXwt) {		     /* Replace X and y by weighted X and y */
	WX = Calloc(n * p, double);
	wy = Calloc(n, double);
	
	AZERO(WX, n * p);
	for (int i = 0; i < n; i++) {
	    wy[i] = sXwt[i] * y[i];
	    for (int j = 0; j < p; j++)
		WX[i + j * n] += sXwt[i] * X[i + j * n];
	}
	X = WX;
	y = wy;
    }
				/* solve L del1 = PAy */
    P_sdmult(pb, (int*)L->Perm, A, y, 1);
    sol = M_cholmod_solve(CHOLMOD_L, L, cpb, &c);
    dble_cpy(pb, (double*)sol->x, q);
    M_cholmod_free_dense(&sol, &c);
				/* solve RX' del2 = X'y - RZX'del1 */
    F77_CALL(dgemv)("T", &n, &p, one, X, &n,
		    y, &i1, zero, pbeta, &i1);
    F77_CALL(dgemv)("T", &q, &p, mone, RZX, &q,
		    pb, &i1, one, pbeta, &i1);
    F77_CALL(dtrsv)("U", "T", "N", &p, RX, &p, pbeta, &i1);
    d[pwrss_POS] = sqr_length(y, n)
	- (sqr_length(pbeta, p) + sqr_length(pb, q));
    if (d[pwrss_POS] < 0)
	error(_("Calculated PWRSS for a LMM is negative"));
    if (wy) Free(wy);
    if (WX) Free(WX);
}


/**
 * Determine the conditional estimates of the fixed effects
 * and the conditional mode of u for a linear mixed model
 *
 * @param x an mer object
 */
static double lmm_update_fixef_u(SEXP x)
{
    double ans = NA_REAL;
    if (!etaGamma_SLOT(x) && !MUETA_SLOT(x)) { /* linear mixed model */
	int *dims = DIMS_SLOT(x), i1 = 1;
	int n = dims[n_POS], p = dims[p_POS], q = dims[q_POS];
	double *d = DEV_SLOT(x), *fixef = FIXEF_SLOT(x),
	    *srwt = SRWT_SLOT(x), *u = U_SLOT(x),
	    dn = (double)n, dnmp = (double)(n-p),
	    mone[] = {-1,0}, one[] = {1,0};
	CHM_FR L = L_SLOT(x);
	CHM_DN cu = N_AS_CHM_DN(u, q, 1), sol;
	R_CheckStack();

	lmm_update_projection(x, u, fixef);
				/* solve RX beta-hat = del2 */
	F77_CALL(dtrsv)("U", "N", "N", &p, RX_SLOT(x), &p, fixef, &i1);
	F77_CALL(dgemv)("N", &q, &p, mone, RZX_SLOT(x), &q, fixef,
			&i1, one, u, &i1);
	if (!(sol = M_cholmod_solve(CHOLMOD_Lt, L, cu, &c)))
	    error(_("cholmod_solve (CHOLMOD_Lt) failed:"));
	dble_cpy(u, (double*)(sol->x), q);
	M_cholmod_free_dense(&sol, &c);

	d[usqr_POS] = sqr_length(u, q);
	d[wrss_POS] = d[disc_POS] = d[pwrss_POS] - d[usqr_POS];
	d[ML_POS] = d[ldL2_POS] +
	    dn * (1 + log(d[pwrss_POS]) + log(2 * PI / dn));
	d[REML_POS] = d[ldL2_POS] + d[ldRX2_POS] +
	    dnmp * (1. + log(d[pwrss_POS]) + log(2. * PI / dnmp));
	d[sigmaML_POS] = sqrt(d[pwrss_POS]/
			      (srwt ? sqr_length(srwt, n) : dn));
	d[sigmaREML_POS] = d[sigmaML_POS] * sqrt(dn/dnmp);

	ans = d[dims[isREML_POS] ? REML_POS : ML_POS];
    }
    return ans;
}


/**
 * Iterate to determine the conditional modes of the random effects.
 *
 * @param x pointer to an mer object
 *
 * @return number of iterations to convergence (0 for non-convergence) 
 */
static int update_u(SEXP x)
{
    int *dims = DIMS_SLOT(x);
    int i, n = dims[n_POS], q = dims[q_POS], verb = dims[verb_POS];
    double *Cx = Cx_SLOT(x), *sXwt = SXWT_SLOT(x),
	*res = RESID_SLOT(x), *u = U_SLOT(x),
	cfac = ((double)n) / ((double)q), 
	crit, pwrss, pwrss_old, step;
    double *tmp = Alloca(q, double), *tmp1 = Alloca(q, double),
	*uold = Alloca(q, double), one[] = {1,0}, zero[] = {0,0};
    CHM_FR L = L_SLOT(x);
    CHM_DN cres = N_AS_CHM_DN(res, n, 1),
	ctmp = N_AS_CHM_DN(tmp, q, 1), sol;
    CHM_SP C = Cm_SLOT(x);
    R_CheckStack();
    
    if (!sXwt) return(0);	/* nothing to do for LMMs */
    if (!(L->is_ll)) error(_("L must be LL', not LDL'"));
    if (q > n) error(_("q = %d > n = %d"), q, n);
    if (Cx) {		    /* A and C have the same structure */
	C = A_SLOT(x);
	R_CheckStack();
	C->x = (void*)Cx;
    }
    
    /* resetting u to zero at the beginning of each evaluation
     * requires more iterations but is necessary to obtain a
     * repeatable evaluation.  If this is not done the optimization
     * algorithm can take wild steps. */
    AZERO(u, q);
    update_mu(x);
    for (i = 0; ; i++) {
	
	dble_cpy(uold, u, q);
	pwrss_old = update_L(x);
				/* tmp := PC %*% wtdResid */
	M_cholmod_sdmult(C, 0 /* notrans */, one, zero, cres, ctmp, &c);
	dble_cpy(tmp1, tmp, q);
	apply_perm(tmp, tmp1, (int*)L->Perm, q);
				/* tmp := tmp - u */
	for (int j = 0; j < q; j++) tmp[j] -= u[j];
				/* solve L %*% sol = tmp */
	if (!(sol = M_cholmod_solve(CHOLMOD_L, L, ctmp, &c)))
	    error(_("cholmod_solve (CHOLMOD_L) failed"));
	dble_cpy(tmp, (double*)(sol->x), q);
	M_cholmod_free_dense(&sol, &c);
				/* evaluate convergence criterion */
	crit = cfac * sqr_length(tmp, q) / pwrss_old;
	if (crit < CM_TOL) break; /* don't do needless evaluations */
				/* solve t(L) %*% sol = tmp */
	if (!(sol = M_cholmod_solve(CHOLMOD_Lt, L, ctmp, &c)))
	    error(_("cholmod_solve (CHOLMOD_Lt) failed"));
	dble_cpy(tmp, (double*)(sol->x), q);
	M_cholmod_free_dense(&sol, &c);

	for (step = 1; step > CM_SMIN; step /= 2) { /* step halving */
	    for (int j = 0; j < q; j++) u[j] = uold[j] + step * tmp[j];
	    pwrss = update_mu(x);
	    if (verb < 0)
		Rprintf("%2d,%8.6f,%12.4g: %15.6g %15.6g %15.6g %15.6g\n",
			i, step, crit, pwrss, pwrss_old, u[1], u[2]);
	    if (pwrss < pwrss_old) {
		pwrss_old = pwrss;
		break;
	    }
	}
	if (step <= CM_SMIN || i > CM_MAXITER) return 0;
    }
    return i;
}

/**
 * Update the fixed effects and the orthogonal random effects in an
 * MCMC sample from an mer object.
 *
 * @param x an mer object
 * @param sigma current standard deviation of the per-observation
 *        noise terms.
 * @param fvals pointer to memory in which to store the updated beta
 * @param rvals pointer to memory in which to store the updated b (may
 *              be (double*)NULL)
 */
static void MCMC_beta_u(SEXP x, double sigma, double *fvals, double *rvals)
{
    int *dims = DIMS_SLOT(x);
    int i1 = 1, p = dims[p_POS], q = dims[q_POS];
    double *V = etaGamma_SLOT(x), *fixef = FIXEF_SLOT(x), *muEta = MUETA_SLOT(x),
	*u = U_SLOT(x), mone[] = {-1,0}, one[] = {1,0};
    CHM_FR L = L_SLOT(x);
    double *del1 = Alloca(q, double), *del2 = Alloca(p, double);
    CHM_DN sol, rhs = N_AS_CHM_DN(del1, q, 1);
    R_CheckStack();

    if (V || muEta) {
	error(_("Update not yet written"));
    } else {			/* Linear mixed model */
	update_L(x);
	update_RX(x);
	lmm_update_fixef_u(x);
				/* Update beta */
	for (int j = 0; j < p; j++) del2[j] = sigma * norm_rand();
	F77_CALL(dtrsv)("U", "N", "N", &p, RX_SLOT(x), &p, del2, &i1);
	for (int j = 0; j < p; j++) fixef[j] += del2[j];
				/* Update u */
	for (int j = 0; j < q; j++) del1[j] = sigma * norm_rand();
	F77_CALL(dgemv)("N", &q, &p, mone, RZX_SLOT(x), &q,
			del2, &i1, one, del1, &i1);
	sol = M_cholmod_solve(CHOLMOD_Lt, L, rhs, &c);
	for (int j = 0; j < q; j++) u[j] += ((double*)(sol->x))[j];
	M_cholmod_free_dense(&sol, &c);
	update_mu(x);	     /* and parts of the deviance slot */
    }
    dble_cpy(fvals, fixef, p);
    if (rvals) {
	update_ranef(x);
	dble_cpy(rvals, RANEF_SLOT(x), q);
    }
}

/**
 * Update the theta_T parameters from the ST arrays in place.
 *
 * @param x an mer object
 * @param sigma current standard deviation of the per-observation
 *        noise terms.
 */
/* FIXME: Probably should fold this function into MCMC_S */
static void MCMC_T(SEXP x, double sigma)
{
    int *Gp = Gp_SLOT(x), nt = (DIMS_SLOT(x))[nt_POS];
    double **st = Alloca(nt, double*);
    int *nc = Alloca(nt, int), *nlev = Alloca(nt, int);
    R_CheckStack();

    if (ST_nc_nlev(GET_SLOT(x, lme4_STSym), Gp, st, nc, nlev) < 2) return;
    error("Code for non-trivial theta_T not yet written");
}

/**
 * Update the theta_S parameters from the ST arrays in place.
 *
 * @param x an mer object
 * @param sigma current standard deviation of the per-observation
 *        noise terms.
 */
static void MCMC_S(SEXP x, double sigma)
{
    CHM_SP A = A_SLOT(x), Zt = Zt_SLOT(x);
    int *Gp = Gp_SLOT(x), *ai = (int*)(A->i),
	*ap = (int*)(A->p), *dims = DIMS_SLOT(x), *perm = PERM_VEC(x);
    int annz = ap[A->ncol], info, i1 = 1, n = dims[n_POS],
	nt = dims[nt_POS], ns, p = dims[p_POS], pos,
	q = dims[q_POS], znnz = ((int*)(Zt->p))[Zt->ncol];
    double *R, *ax = (double*)(A->x), *b = RANEF_SLOT(x),
	*eta = ETA_SLOT(x), *offset = OFFSET_SLOT(x),  
	*rr, *ss, one = 1, *u = U_SLOT(x), *y = Y_SLOT(x);
    int *nc = Alloca(nt, int), *nlev = Alloca(nt, int),
	*spt = Alloca(nt + 1, int);
    double **st = Alloca(nt, double*);
    R_CheckStack();

    ST_nc_nlev(GET_SLOT(x, lme4_STSym), Gp, st, nc, nlev);
    ns = 0;			/* ns is length(theta_S) */
    spt[0] = 0;			/* pointers into ss for terms */
    for (int i = 0; i < nt; i++) {
	ns += nc[i];
	spt[i + 1] = spt[i] + nc[i];
    }

    if (annz == znnz) { /* Copy Z' to A unless A has new nonzeros */
	dble_cpy(ax, (double*)(Zt->x), znnz);
    } else error("Code not yet written for MCMC_S with NLMMs");
				/* Create T'Zt in A */
    Tt_Zt(A, Gp, nc, nlev, st, nt); 
				/* Create P'u in ranef slot */
    for (int i = 0; i < q; i++) b[perm[i]] = u[i];
				/* Create X\beta + offset in eta slot */
    for (int i = 0; i < n; i++) eta[i] = offset ? offset[i] : 0;
    F77_CALL(dgemv)("N", &n, &p, &one, X_SLOT(x), &n,
		    FIXEF_SLOT(x), &i1, &one, eta, &i1);
				/* Allocate R, rr and ss */
    R = Alloca(ns * ns, double); /* crossproduct matrix then factor */
    rr = Alloca(ns, double);	 /* row of model matrix for theta_S */
    ss = Alloca(ns, double);	 /* right hand side, then theta_S */
    R_CheckStack();
    AZERO(R, ns * ns);
    AZERO(ss, ns);
    /* Accumulate crossproduct from pseudo-data part of model matrix */
    for (int i = 0; i < q; i++) {
	int sj = theta_S_ind(i, nt, Gp, nlev, spt);
	AZERO(rr, ns);
	rr[sj] = b[i];
	F77_CALL(dsyr)("U", &ns, &one, rr, &i1, R, &ns);
    }	
    /* Accumulate crossproduct and residual product of the model matrix. */
    /* This is done one row at a time.  Rows of the model matrix
     * correspond to columns of T'Zt */
    for (int j = 0; j < n; j++) { /* jth column of T'Zt */
	AZERO(rr, ns);
	for (int p = ap[j]; p < ap[j + 1]; p++) {
	    int i = ai[p];	/* row in T'Zt */
	    int sj = theta_S_ind(i, nt, Gp, nlev, spt);

	    rr[sj] += ax[p] * b[i];
	    ss[sj] += rr[sj] * (y[j] - eta[j]);
	}
	F77_CALL(dsyr)("U", &ns, &one, rr, &i1, R, &ns);
    }
    F77_CALL(dposv)("U", &ns, &i1, R, &ns, ss, &ns, &info);
    if (info)
	error(_("Model matrix for theta_S is not positive definite, %d."), info);
    for (int j = 0; j < ns; j++) rr[j] = sigma * norm_rand();
    /* Sample from the conditional Gaussian distribution */
    F77_CALL(dtrsv)("U", "N", "N", &ns, R, &ns, rr, &i1);
    for (int j = 0; j < ns; j++) ss[j] += rr[j];
    /* Copy positive part of solution onto diagonals of ST */
    pos = 0;
    for (int i = 0; i < nt; i++) {
	for (int j = 0; j < nc[i]; j++) {
	    st[i][j * (nc[i] + 1)] = (ss[pos] > 0) ? ss[pos] : 0;
	    pos++;
	}
    }
    update_A(x);
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
    SEXP devsamp = GET_SLOT(x, lme4_devianceSym);
    int *dims = DIMS_SLOT(x), nsamp = LENGTH(devsamp);
    int n = dims[n_POS], np = dims[np_POS],
	p = dims[p_POS], q = dims[q_POS];
    double
	*STsamp = REAL(GET_SLOT(x, lme4_STSym)),
	*d = DEV_SLOT(fm), *dev = REAL(devsamp),
	*sig = SLOT_REAL_NULL(x, lme4_sigmaSym),
	*fixsamp = FIXEF_SLOT(x), *resamp = RANEF_SLOT(x);

    GetRNGstate();
    /* The first column of storage slots contains the fitted values */
    for (int i = 1; i < nsamp; i++) {
/* FIXME: This is probably wrong for a model with weights. */
	if (sig) 		/* update and store sigma */
	    sig[i] = sqrt(d[pwrss_POS]/rchisq((double)(n + q)));
			/* update L, RX, beta, u, eta, mu, res and d */
	MCMC_beta_u(fm, sig ? sig[i] : 1, fixsamp + i * p,
		    resamp + (resamp ? i : 0) * q); 
	dev[i] = d[ML_POS];
				/* update theta_T, theta_S and A */
 	MCMC_T(fm, sig ? sig[i] : 1);
 	MCMC_S(fm, sig ? sig[i] : 1);
	ST_getPars(fm, STsamp + i * np); /* record theta */
    }
    PutRNGstate();
		/* Restore pars from the first columns of the samples */
    dble_cpy(FIXEF_SLOT(fm), fixsamp, p);
    ST_setPars(fm, STsamp);
    update_ranef(fm);

    return x;
}

/**
 * Externally callable update_L.
 * Update the A, L, sqrtrWt and resid slots.  It is assumed that
 * update_mu has already been called at the current values of u and
 * the model parameters.
 *
 * @param x an mer object
 *
 * @return penalized weighted residual sum of squares
 *
 */
SEXP mer_update_L(SEXP x){return ScalarReal(update_L(x));}

/**
 * Externally callable update_u.
 *
 * Iterate to determine the conditional modes of the random effects.
 *
 * @param x pointer to an mer object
 *
 * @return number of iterations to convergence (0 for non-convergence) 
 */
SEXP mer_update_u(SEXP x){return ScalarInteger(update_u(x));}

/**
 * Externally callable lmm_update_projection.
 * Create the projections onto the column spaces of the random effects
 * and the fixed effects.
 *
 * @param x an mer object
 *
 * @return a list with two elements, both REAL vectors
 */
SEXP mer_update_projection(SEXP x)
{
    SEXP ans = PROTECT(allocVector(VECSXP, 2));
    int *dims = DIMS_SLOT(x);

    SET_VECTOR_ELT(ans, 0, allocVector(REALSXP, dims[q_POS]));
    SET_VECTOR_ELT(ans, 1, allocVector(REALSXP, dims[p_POS]));
    lmm_update_projection(x, REAL(VECTOR_ELT(ans, 0)),
			  REAL(VECTOR_ELT(ans, 1)));

    UNPROTECT(1);
    return ans;
}

/**
 * Externally callable update_RX
 *
 * @param x pointer to an mer object
 *
 * @return profiled deviance or REML deviance
 */
SEXP
mer_update_RX(SEXP x)
{
    return ScalarReal(update_RX(x));
}

/**
 * Create PAX in dest.
 *
 * @param dest values to be calculated
 * @param perm NULL or a 0-based permutation vector defining P
 * @param A sparse matrix
 * @param X dense matrix
 * @param nc number of columns in X
 *
 */
static void
P_sdmult(double *dest, const int *perm, const CHM_SP A,
	 const double *X, int nc)
{
    int *ai = (int*)(A->i), *ap = (int*)(A->p), m = A->nrow, n = A->ncol;
    double *ax = (double*)(A->x), *tmp = Alloca(m, double);
    R_CheckStack();

    for (int k = 0; k < nc; k++) {
	AZERO(tmp, m);
	for (int j = 0; j < n; j++) {
	    for (int p = ap[j]; p < ap[j + 1]; p++)
		tmp[ai[p]] += X[j + k * n] * ax[p];
	}
	apply_perm(dest + k * m, tmp, perm, m);
    }
}

SEXP mer_optimize(SEXP x)
{
    SEXP ST = GET_SLOT(x, lme4_STSym);
    int *dims = DIMS_SLOT(x);
    int lmm = !(MUETA_SLOT(x) || ETAGAMMA_SLOT(x)), nt = dims[nt_POS];
    int nv = dims[np_POS] + (lmm ? 0 : dims[p_POS]), verb = dims[verb_POS];
    int liv = S_iv_length(OPT, nv), lv = S_v_length(OPT, nv);
    double *g = (double*)NULL, *h = (double*)NULL, fx = R_PosInf;
    double *fixef = FIXEF_SLOT(x);
    int *iv = Alloca(liv, int);
    double *b = Alloca(2 * nv, double), *d = Alloca(nv, double),
	*v = Alloca(lv, double), *xv = Alloca(nv, double);
    R_CheckStack();

    ST_getPars(x, xv);
    if (!lmm) {
	double eta = 1.e-5; /* estimated rel. error on computed lpdisc */

	dble_cpy(xv + dims[np_POS], fixef, dims[p_POS]);
	v[31] = eta;		/* RFCTOL */
	v[36] = eta;		/* SCTOL */
	v[41] = eta;		/* ETA0 */
    }
				/* initialize the state vectors v and iv */
    S_Rf_divset(OPT, iv, liv, lv, v);
    iv[OUTLEV] = (verb < 0) ? -verb : verb;
    iv[MXFCAL] = dims[mxfn_POS];
    iv[MXITER] = dims[mxit_POS];
				/* set the bounds to plus/minus Infty  */
    for (int i = 0; i < nv; i++) {
	b[2*i] = R_NegInf; b[2*i+1] = R_PosInf; d[i] = 1;
    }
				/* reset lower bounds on elements of theta_S */
    for (int i = 0, pos = 0; i < nt; i++) {
	int nc = *INTEGER(getAttrib(VECTOR_ELT(ST, i), R_DimSymbol));
	for (int j = 0; j < nc; j++) b[pos + 2*j] = 0;
	pos += nc * (nc + 1);
    }

    do {
	if (!lmm)	  /* Some time soon we will get rid of this */
	    dble_cpy(fixef, xv + dims[np_POS], dims[p_POS]);
	ST_setPars(x, xv);	/* update ST and A etc. */
	fx = update_dev(x);
	S_nlminb_iterate(b, d, fx, g, h, iv, liv, lv, nv, v, xv);
    } while (iv[0] == 1 || iv[0] == 2);
    ST_setPars(x, xv);
    update_ranef(x);
    dims[cvg_POS] = iv[0];
    return R_NilValue;
}

/* this is the nonlinear part of update_dev.  It should all be folded
 * into one function and use the definition of the likelihood methods
 * here instead of duplicating all this code.
 */

    else {  /* NLMM */
	double dn = (double) dims[n_POS];
	
	d[disc_POS] = d[wrss_POS];
	if (nAGQ > 1) {
	    /* Adaptive Gauss-Hermite quadrature */
	    /* Single grouping factor has been checked. */
	    const int nl = nlevels(VECTOR_ELT(flistP, 0));
	    const int nre = q / nl;
	    int *fl0 = INTEGER(VECTOR_ELT(flistP, 0)), *pointer = Alloca(nre, int);
	    double *ghw = GHW_SLOT(x), *ghx = GHX_SLOT(x),
		*res = RESID_SLOT(x), *tmp = Calloc(nl, double),
		/* store conditional mode */
		*uold = dble_cpy(Calloc(q, double), u, q),
		w_pro = 1, z_sum = 0;           /* values needed in AGQ evaluation */
	    const double sigma = d[sigmaML_POS];   /* MLE of sigma */
	    const double factor = - 1 / (2 * sigma * sigma);

	    R_CheckStack();

	    AZERO(pointer, nre);
	    AZERO(tmp, nl);
	    
	    d[ML_POS] = dn * log(2*PI*d[pwrss_POS]/dn) + d[ldL2_POS];
	    
	    /* implementation of AGQ method (Laplacian will be a trivial case) */
	    AZERO(pointer, nre);                    /* assign initial pointers, all 0 */
	    AZERO(tmp, nl);

	    /* add accuracy to integration approximation */
	    while(pointer[nre - 1] < nAGQ){
		double *z = Calloc(q, double);       /* current abscissas */
		double *presid = Calloc(nl, double); /* current penalized residuals in different levels */
		
		/* update abscissas and weights */
		for(int i = 0; i < nre; ++i){
		    for(int j = 0; j < nl; ++j){
			z[i + j * nre] = ghx[pointer[i]];
		    }
		    z_sum += ghx[pointer[i]] * ghx[pointer[i]];
		    w_pro *= ghw[pointer[i]];
		}
		CHM_DN cz = N_AS_CHM_DN(z, q, 1), sol;
		if(!(sol = M_cholmod_solve(CHOLMOD_L, L, cz, &c)))
		    error(_("cholmod_solve(CHOLMOD_L) failed"));
		dble_cpy(z, (double *)sol->x, q);
		M_cholmod_free_dense(&sol, &c);
		for(int i = 0; i < q; ++i){
		    u[i] = uold[i] + sigma * z[i];
		}
		update_mu(x);
		
		AZERO(presid, nl);
		for(int i = 0; i < dims[n_POS]; ++i){
		    presid[fl0[i]-1] += ( res[i] * res[i] );
		}
		
		for(int i = 0; i < nre; ++i){
		    for(int j = 0; j < nl; ++j)
			presid[j] += u[i + j * nre] * u[i + j * nre];
		}
		
		for(int j = 0; j < nl; ++j){
		    tmp[j] += exp(factor * presid[j] + z_sum) * w_pro / sqrt(PI);
		}
		
		/* move pointer to next combination of weights and abbsicas */
		int count = 0;
		pointer[count]++;
		while(pointer[count] == nAGQ && count < nre - 1){
		    pointer[count] = 0;
		    pointer[++count]++;
		}
		if(z) Free(z);
		w_pro = 1;
		z_sum = 0;
		if(presid) Free(presid);
	    }
	    
	    for(int j = 0; j < nl; ++j){
		d[ML_POS] -= ( 2 * log(tmp[j]) );
	    }
	    
	    dble_cpy(u, uold, q);
	    update_mu(x);
	    if(tmp)   Free(tmp);
	    if(uold)  Free(uold);
	    
	}
	else{
	    d[ML_POS] = dn*(1 + log(d[pwrss_POS]) + log(2*PI/dn)) + d[ldL2_POS];
        }
    }



#endif
