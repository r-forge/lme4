#include "mer.h"

using namespace Rcpp;
using namespace MatrixNs;

static void showdbl(const double* x, const char* nm, int n) {
    int n5 = (n < 5) ? n : 5;
    Rprintf("%s: %g", nm, x[0]);
    for (int i = 1; i < n5; i++) Rprintf(", %g", x[i]);
    Rprintf("\n");
}

#if 0
static void showint(const int* x, const char* nm, int n) {
    int n20 = (n < 20) ? n : 20;
    Rprintf("%s: %d", nm, x[0]);
    for (int i = 1; i < n20; i++) Rprintf(", %d", x[i]);
    Rprintf("\n");
}

static void showCHM_DN(const CHM_DN x, const std::string &nm) {
    Rprintf("%s: nrow = %d, ncol = %d, nzmax = %d, d = %d, xtype = %d, dtype = %d\n",
	    nm.c_str(), x->nrow, x->ncol, x->nzmax, x->d, x->xtype, x->dtype);
    showdbl((double*)x->x, "x", x->nzmax);
}

static void showCHM_FR(const CHM_FR x, const std::string &nm) {
    Rprintf("%s: n = %d, minor = %d, xtype = %d, itype = %d, dtype = %d\n",
	    nm.c_str(), x->n, x->minor, x->xtype, x->itype, x->dtype);
    Rprintf("ordering = %d, is_ll = %d, is_super = %d, is_monotonic = %d, nzmax = %d\n",
	    x->ordering, x->is_ll, x->is_super, x->is_monotonic, x->nzmax);
    int *pp = (int*)x->p;
    showint(pp, "p", (x->n) + 1);
    showint((int*)x->Perm, "Perm", x->n);
    if (pp[x->n] > 0) {
	showint((int*)x->i, "i", x->n);
	showdbl((double*)x->x, "x", x->n);
    } else Rprintf("nnz = 0\n");
}

static void showCHM_SP(const CHM_SP x, const std::string &nm) {
    Rprintf("%s: nrow = %d, ncol = %d, xtype = %d, stype = %d, itype = %d, dtype = %d\n",
	    nm.c_str(), x->nrow, x->ncol, x->xtype, x->stype,
	    x->itype, x->dtype);
    int nc = x->ncol, nnz = M_cholmod_nnz(x, &c);
    showint((int*)x->p, "p", nc + 1);
    if (nnz > 0) {
	showint((int*)x->i, "i", nnz);
	showdbl((double*)x->x, "x", nnz);
    } else Rprintf("nnz = 0\n");
}
#endif

namespace mer{

    reModule::reModule(S4 xp) :
	d_L(S4(xp.slot("L"))),
	d_Lambda(S4(xp.slot("Lambda"))),
	d_Ut(S4(xp.slot("Ut"))),
	d_Zt(S4(xp.slot("Zt"))),
	Lind(xp.slot("Lind")),
	lower(xp.slot("lower")),
	theta(xp.slot("theta")),
	d_ldL2(NumericVector(xp.slot("ldL2")).begin()),
	d_u(xp.slot("u"))
    {}
    
    void reModule::setU(const double *nu) {
	std::copy(nu, nu + d_u.size(), d_u.begin());
    }
    
/** 
 * Lambda@x[] <- theta[Lind]; Ut <- crossprod(Lambda, Zt);
 * update(L,Ut,1); ldL2
 * 
 * @param nt New value of theta
 */
    void reModule::updateTheta(NumericVector const &nt) {
	if (nt.size() != theta.size())
	    ::Rf_error("length(theta) = %d != length(newtheta) = %d",
		       theta.size(), nt.size());
	std::copy(nt.begin(), nt.end(), theta.begin());
	double *Lamx = (double*)d_Lambda.x, *nnt = nt.begin();
	int *Li = Lind.begin();
	for (int i = 0; i < Lind.size(); i++) Lamx[i] = nnt[Li[i] - 1];

	CHM_SP LamTrZt = d_Lambda.crossprod(d_Zt);
	d_Ut.update(*LamTrZt);
	M_cholmod_free_sparse(&LamTrZt, &c);
	d_L.update(d_Ut, 1.);
	*d_ldL2 = M_chm_factor_ldetL2(&d_L);
    }

    void reModule::rwUpdateL(NumericMatrix const &Xwt,
			     NumericVector const &wtres,
			     NumericVector const &u,
			     double *ans) {
	double one = 1., mone = -1;
	CHM_SP wtUt = M_cholmod_copy_sparse(&d_Ut, &c);
	const chmDn cXwt(Xwt);
	if (d_Ut.ncol == d_Zt.ncol) {
	    M_cholmod_scale(&cXwt, CHOLMOD_COL, wtUt, &c);
	} else Rf_error("Multiple columns in Xwt");
	M_cholmod_factorize_p(wtUt, &one, (int*)NULL, (size_t)0,
			      &d_L, &c);
	*d_ldL2 = M_chm_factor_ldetL2(&d_L);
	std::copy(u.begin(), u.end(), ans);
	chmDn cans(ans, d_Ut.nrow, 1);
	const chmDn cwtres(wtres);
	M_cholmod_sdmult(wtUt, 0/*trans*/, &one, &mone, &cwtres,
			 &cans, &c);
	M_cholmod_free_sparse(&wtUt, &c);
	CHM_DN incr = M_cholmod_solve(CHOLMOD_A, &d_L, &cans, &c);
	double *incx = (double*)incr->x;
	std::copy(incx, incx + d_Ut.nrow, ans);
	M_cholmod_free_dense(&incr, &c);
    }

/** 
 * @return squared length of u
 */
    double reModule::sqLenU() const {
	return std::inner_product(d_u.begin(), d_u.end(),
				  d_u.begin(), double());
    }

/** 
 * Solve for u given cu from the response module
 * 
 * @param resp Response module
 */
    void reModule::updateU(const merResp &resp) {
	chmDn cu(resp.cu);
	CHM_DN t1 = d_L.solve(CHOLMOD_Lt, &cu);
	CHM_DN t2 = d_L.solve(CHOLMOD_Pt, t1);
	::M_cholmod_free_dense(&t1, &c);
	double *t2b = (double*)t2->x;
	std::copy(t2b,  t2b + d_u.size(), d_u.begin());
	::M_cholmod_free_dense(&t2, &c);
    }

/** 
 * Add the contribution of the random effects to the linear predictor
 * vector.
 * 
 * @param gam Current gamma vector
 */
    void reModule::incGamma(double *gam) const {
	chmDn cgam(gam, d_Ut.ncol, 1);
	d_Ut.dmult('T', 1., 1., chmDn(d_u), cgam);
    }
	
    void reModule::incGamma(NumericVector &gam) const {
	chmDn cgam(gam);
	d_Ut.dmult('T', 1., 1., chmDn(d_u), cgam);
    }
	
/** 
 * Update the destination dense matrix as
 * solve(L, solve(L, Lambda %*% src, system = "P"), system = "L")
 * 
 * @param src source dense matrix
 * @param dest destination dense matrix
 */
    void reModule::DupdateL(chmDn const &src, chmDn &dest) const {
	d_Lambda.dmult('T', 1., 0., src, dest);
	CHM_DN t1 = d_L.solve(CHOLMOD_P, &dest);
	CHM_DN t2 = d_L.solve(CHOLMOD_L, t1);
	::M_cholmod_free_dense(&t1, &c);

	double *t2b = (double*)t2->x, *db = (double*)(dest.x);
	int sz = src.nr() * src.nc();
	std::copy(t2b, t2b + sz, db);
	::M_cholmod_free_dense(&t2, &c);
    }

/** 
 * Return
 * solve(L, solve(L, Lambda %*% src, system = "P"), system = "L")
 * 
 * @param src source sparse matrix
 * @return as above
 */
    CHM_SP reModule::SupdateL(chmSp const &src) const {
	CHM_SP t1 = d_Lambda.crossprod(src);
	CHM_SP t2 = d_L.spsolve(CHOLMOD_P, t1);
	::M_cholmod_free_sparse(&t1, &c);
	t1 = d_L.spsolve(CHOLMOD_L, t2);
	::M_cholmod_free_sparse(&t2, &c);
	return t1;
    }
    
    double reModule::ldL2() const { return *d_ldL2; }

    merResp::merResp(S4 xp)
	: d_wrss(NumericVector(xp.slot("wrss")).begin()),
	  d_offset(xp.slot("offset")),
	  d_sqrtrwt(xp.slot("sqrtrwt")),
	  d_wtres(xp.slot("wtres")),
	  mu(xp.slot("mu")),
	  weights(xp.slot("weights")),
	  y(xp.slot("y")),
	  Utr(xp.slot("Utr")),
	  Vtr(xp.slot("Vtr")),
	  cbeta(xp.slot("cbeta")),
	  cu(xp.slot("cu"))
    {
	int n = y.size(), os = d_offset.size();

	if (mu.size() != n || d_wtres.size() != n ||
	    weights.size() != n || d_sqrtrwt.size() != n)
	    ::Rf_error("y, mu, sqrtrwt, wtres and weights slots must have equal lengths");
	if (os < 1 || os % n)
	    ::Rf_error("length(offset) must be a positive multiple of length(y)");
	if (cbeta.size() != Vtr.size())
	    ::Rf_error("cbeta and Vtr slots must have equal lengths");
	if (cu.size() != Utr.size())
	    ::Rf_error("cu and Utr slots must have equal lengths");
    }

    double merResp::wrss() const { return *d_wrss; }
    const NumericVector &merResp::sqrtrwt() const { return d_sqrtrwt; }
    const NumericVector &merResp::offset() const { return d_offset; }
    const NumericVector &merResp::wtres() const { return d_wtres; }
/** 
 * Update cu using Lambda and L
 *  cu <- solve(L, solve(L, crossprod(Lambda, Utr),
 *                       sys = "P"), 
 *              sys = "L")
 * 
 * @param re random-effects module
 */
    void merResp::updateL(reModule const &re) {
	chmDn ccu(cu);
	re.DupdateL(chmDn(Utr), ccu);
    }

    void merResp::updateMu(NumericVector const &gamma) {
	std::copy(gamma.begin(), gamma.end(), mu.begin());
    }

    void glmerResp::updateMu(NumericVector const &gamma) {
	std::copy(gamma.begin(), gamma.end(), d_gamma.begin());
	linkInv();
	MuEta();
    }

    void nlmerResp::updateMu(NumericVector const &gamma) {
	int n = y.size();
	double *gg = gamma.begin();

	showdbl(gamma.begin(), "gamma in resp::updateMu", gamma.size());
	for (int p = 0; p < pnames.size(); p++) {
	    std::string pn(pnames[p]);
	    Rprintf("p = %d, parameter name = %s\n", p, &pn[0]);
	    NumericVector pp = nlenv.get(pn);
	    showdbl(pp.begin(), "current", pp.size());
	    std::copy(gg + n * p, gg + n * (p + 1), pp.begin());
	    showdbl(pp.begin(), "updated", pp.size());
	}
	Rprintf("\nBefore evaluation\n");
	NumericVector rr = nlmod.eval(nlenv);
	showdbl(rr.begin(), "Model evaluation", rr.size());
	if (rr.size() != n)
	    Rf_error("Length mu = %d, expected %d", rr.size(), n);
	std::copy(rr.begin(), rr.end(), mu.begin());
	NumericMatrix rrg = rrg.attr("gradient");
	std::copy(rrg.begin(), rrg.end(), d_sqrtXwt.begin());
	showdbl(rrg.begin(), "Gradient evaluation", rrg.size());
    }

/** 
 * Update the wtres vector and return the sum of squares
 *   wtres <- sqrtrwts * (y - mu)
 *   return(wrss <- sum(wtres^2))
 *
 * @return Updated weighted residual sum of squares
 */
    double merResp::updateWrss() {
				// wtres <- y - mu
	std::transform(y.begin(), y.end(), mu.begin(),
		       d_wtres.begin(), std::minus<double>());
				// wtres <- wtres * sqrtrwt
	std::transform(d_wtres.begin(), d_wtres.end(), d_sqrtrwt.begin(),
		       d_wtres.begin(), std::multiplies<double>());
	*d_wrss = std::inner_product(d_wtres.begin(), d_wtres.end(),
				     d_wtres.begin(), double());
	return *d_wrss;
    }

    rwResp::rwResp(S4 xp)
	: merResp(xp),
	  d_gamma(xp.slot("gamma")),
	  d_sqrtXwt(SEXP(xp.slot("sqrtXwt"))) {
    }

    const NumericMatrix &rwResp::sqrtXwt() const {
	return d_sqrtXwt;
    }

    double *rwResp::gamma() {
	return d_gamma.begin();
    }

    glmerResp::glmerResp(S4 xp)
	: rwResp(xp),
	  family(SEXP(xp.slot("family"))),
	  muEta(xp.slot("muEta")),
	  n(xp.slot("n")),
	  d_var(xp.slot("var")) {
    }

    nlmerResp::nlmerResp(S4 xp)
	: rwResp(xp),
	  nlenv(SEXP(xp.slot("nlenv"))),
	  nlmod(SEXP(xp.slot("nlmod"))),
	  pnames(xp.slot("pnames")) {
    }

    feModule::feModule(S4 xp)
	: d_beta(xp.slot("beta")),
	  d_ldRX2(NumericVector(xp.slot("ldRX2")).begin()) {
    }

    double feModule::ldRX2() const {
	return *d_ldRX2;
    }

    void feModule::setBeta(NumericVector const &nbeta) {
	if (nbeta.size() != d_beta.size())
	    Rf_error("length(beta) = %d should be %d",
		     nbeta.size(), d_beta.size());
	std::copy(nbeta.begin(), nbeta.end(), d_beta.begin());
    }

    const NumericVector &feModule::beta() const {
	return d_beta;
    }

    deFeMod::deFeMod(S4 xp) :
	feModule(xp),
	d_X(S4(xp.slot("X"))),
	d_RZX(S4(xp.slot("RZX"))),
	d_RX(S4(xp.slot("RX")))
    { }

    lmerDeFeMod::lmerDeFeMod(S4 xp) :
	deFeMod(xp),
	ZtX(S4(xp.slot("ZtX"))),
	XtX(S4(xp.slot("XtX")))
    { }

/** 
 * Update RZX and RX
 *   RZX <<- solve(L, solve(L, crossprod(Lambda, ZtX), 
 *                          sys = "P"),
 *                 sys = "L")
 *   RX <<- chol(XtX - crossprod(RZX))
 * 
 * @param re a random-effects module
 */
    void lmerDeFeMod::updateRzxRx(reModule const &re) {
	chmDn cRZX(d_RZX);
	re.DupdateL(chmDn(ZtX), cRZX);
	d_RX.update('T', -1., d_RZX, 1., XtX);
	*d_ldRX2 = d_RX.logDet2();
    }

 /** 
  * Update beta
  * 	resp@cbeta <- Vtr - crossprod(RZX, cu)
  *	beta <- solve(RX, solve(t(RX), resp@cbeta))
  *	resp@cu <- resp@cu - RZX %*% beta
  * 
  * @param resp response module
  */
    void lmerDeFeMod::updateBeta(merResp &resp) {
	std::copy(resp.Vtr.begin(), resp.Vtr.end(), resp.cbeta.begin());
	d_RZX.dgemv('T', -1., resp.cu, 1., resp.cbeta);
	std::copy(resp.cbeta.begin(), resp.cbeta.end(), d_beta.begin());
	d_RX.dpotrs(d_beta);
	d_RZX.dgemv('N', -1., d_beta, 1., resp.cu);
    }

/** 
 *  gamma += crossprod(Ut, u)
 * 
 * @param gam Value of gamma to be updated
 */
    void deFeMod::incGamma(double *gam) const {
	d_X.dgemv('N', 1., d_beta, 1., gam);
    }
    
    void deFeMod::incGamma(NumericVector &gam) const {
	d_X.dgemv('N', 1., d_beta, 1., gam.begin());
    }
    
    rwDeFeMod::rwDeFeMod(S4 xp)
	: deFeMod(xp),
	  d_V(S4(xp.slot("V"))) {
    }

    const dgeMatrix &rwDeFeMod::V() const {
	return d_V;
    }

    spFeMod::spFeMod(Rcpp::S4 xp) :
	feModule(xp),
	X(S4(xp.slot("X"))),
	RZX(S4(xp.slot("RZX"))),
	RX(S4(xp.slot("RX")))
    { }

    lmerSpFeMod::lmerSpFeMod(Rcpp::S4 xp) :
	spFeMod(xp),
	ZtX(Rcpp::S4(SEXP(xp.slot("ZtX")))),
	XtX(Rcpp::S4(SEXP(xp.slot("XtX"))))
    { }

/** 
 * Update RZX and RX
 *   RZX <<- solve(L, solve(L, crossprod(Lambda, ZtX), 
 *                          sys = "P"),
 *                 sys = "L")
 *   RX <- chol(XtX - crossprod(RZX))
 * 
 * @param re a random-effects module
 */
    void lmerSpFeMod::updateRzxRx(reModule const &re) {
	double mone[] = {-1.,0}, one[] = {1.,0};
	CHM_SP t2 = re.SupdateL(ZtX);

	RZX.update(*t2);
	M_cholmod_free_sparse(&t2, &c);	
	
	CHM_SP t1 = RZX.crossprod();
	t2 = M_cholmod_add(&XtX, t1, one, mone, 1/*values*/,
			     1/*sorted*/, &c);
	M_cholmod_free_sparse(&t1, &c);
	RX.update(*t2);

	M_cholmod_free_sparse(&t2, &c);

	*d_ldRX2 = M_chm_factor_ldetL2(&RX);
    }

 /** 
  * Update beta
  * 	resp@cbeta <- Vtr - crossprod(RZX, cu)
  *	beta <- solve(RX, solve(t(RX), resp@cbeta))
  *	resp@cu <- resp@cu - RZX %*% beta
  * 
  * @param resp response module
  */
    void lmerSpFeMod::updateBeta(merResp &resp) {
	std::copy(resp.Vtr.begin(), resp.Vtr.end(), resp.cbeta.begin());
	chmDn ccbeta(resp.cbeta), ccu(resp.cu);
	RZX.dmult('T', -1., 1., ccu, ccbeta);
	std::copy(resp.cbeta.begin(), resp.cbeta.end(), d_beta.begin());
	chmDn cbeta(d_beta);
	CHM_DN t1 = RX.solve(CHOLMOD_A, &cbeta);
	double *t1b = (double*)t1->x;
	std::copy(t1b,  t1b + d_beta.size(), d_beta.begin());
	M_cholmod_free_dense(&t1, &c);
	RZX.dmult('N', -1., 1., cbeta, ccu);
    }

    void spFeMod::incGamma(NumericVector &gam) const {
	chmDn bb(d_beta), gg(gam);
	X.dmult('N', 1., 1., bb, gg);
    }
    
    static inline double sqrtquotient(double x, double y) {
	return sqrt(x/y);
    }

    void glmerResp::updateSqrtRWt() {
	std::transform(weights.begin(), weights.end(), d_var.begin(),
		       d_sqrtrwt.begin(), sqrtquotient);
    }

    void glmerResp::updateSqrtXWt() {
	MuEta();
	std::transform(d_sqrtrwt.begin(), d_sqrtrwt.end(),
		       muEta.begin(), d_sqrtXwt.begin(),
		       std::multiplies<double>());
    }
    
    glmer::glmer(S4 xp) :
	re(S4(xp.slot("re"))),
	resp(S4(xp.slot("resp"))) {
    }
    
    glmerDe::glmerDe(S4 xp) :
	glmer(xp),
	fe(Rcpp::S4(xp.slot("fe"))) {
    }

    glmerSp::glmerSp(S4 xp) :
	glmer(xp),
	fe(S4(xp.slot("fe"))) {
    }

    void glmerDe::updateGamma() {
	const NumericVector u = re.u(), offset = resp.offset();
	NumericVector b(u.size()), gamma(offset.size());
	std::copy(offset.begin(), offset.end(), gamma.begin());
	chmDn bb(b), gg(gamma);

	re.Lambda().dmult('N', 1., 0., chmDn(u), bb);
	re.Zt().dmult('T', 1., 1., bb, gg);
	fe.X().dmult('N', 1., 1., chmDn(fe.beta()), gg);
	resp.updateMu(gamma);
    }

    void rwDeFeMod::updateV(NumericMatrix const &sqrtXwt) {
	int nc = sqrtXwt.ncol();
	if (nc != 1) Rf_error("code for nlmer not yet written");
// FIXME: Rewrite this using std::transform
	std::copy(d_X.x.begin(), d_X.x.end(), d_V.x.begin());
	int m = d_X.nrow(), n = d_X.ncol();
	double *vv = d_V.x.begin(), *ww = sqrtXwt.begin();
	for (int j = 0; j < n; j++)
	    for (int i = 0; i < m; i++) vv[i + j * m] *= ww[i];
    }

    void rwDeFeMod::updateRX(bool useRZX) {
	if (useRZX) throw std::range_error("should not happen");
	d_RX.update(d_V);
    }

#define CM_TOL 1.e-5
#define CM_MAXITER 30
#define CM_SMIN 1.e-4

/**
 * Determine the Euclidean distance between two vectors relative to the
 * square root of the product of their lengths
 */
    static double compare_vec(const double *v1, const double *v2, int n) {
	double num = 0., d1 = 0., d2 = 0.;
	for (int i = 0; i < n; i++) {
	    double diff = v1[i] - v2[i];
	    num += diff * diff;
	    d1 += v1[i] * v1[i];
	    d2 += v2[i] * v2[i];
	}
	num = sqrt(num); d1 = sqrt(d1); d2 = sqrt(d2);
	if (d1 == 0) {
	    if (d2 == 0) return num;
	    else return num/d2;
	}
	if (d2 == 0) return num/d1;
	return num/sqrt(d1 * d2);
    }

    double glmerDe::IRLS(int verbose) {
	const NumericVector var = resp.var();
	std::vector<double> varold(var.size());
	NumericVector bb = fe.beta();
	int p = bb.size();
	NumericVector betabase(p), newbeta(p), incr(p); 
	double crit, step, pwrss0, pwrss1;

	crit = 10. * CM_TOL;
	    // weights, var and sqrtrwt are assumed established.
	for (int i = 0; crit >= CM_TOL && i < CM_MAXITER; i++) {
				// store copies of var and beta
	    std::copy(bb.begin(), bb.end(), betabase.begin());
	    std::copy(var.begin(), var.end(), varold.begin());
	    updateGamma();     // using current beta
	    resp.linkInv();    // mu
	    pwrss0 = resp.updateWrss() + re.sqLenU();
	    resp.updateSqrtXWt();      // muEta and sqrtXwt
	    fe.updateV(resp.sqrtXwt());	// derive V from X
				// Vtr <- crossprod(V, wtres)
	    fe.V().dgemv('T', 1., resp.wtres(), 0., resp.Vtr);
	    fe.updateRX(false);	// factor V'V
	    std::copy(resp.Vtr.begin(), resp.Vtr.end(), incr.begin());
	    fe.RX().dpotrs(incr); // evaluate the increment
	    if (verbose > 1) showdbl(&incr[0], "deltab", p);
	    pwrss1 = pwrss0;	// force one evaluation of the loop
	    for (step = 1.; pwrss0 <= pwrss1 && step > CM_SMIN; step /= 2.) {
				// newbeta <- betabase + incr * mult
		std::transform(incr.begin(), incr.end(), newbeta.begin(),
		       std::bind2nd(std::multiplies<double>(), step));
		std::transform(betabase.begin(), betabase.end(),
			       newbeta.begin(), newbeta.begin(),
			       std::plus<double>());
		fe.setBeta(newbeta);
		updateGamma();
		resp.linkInv();
		pwrss1 = resp.updateWrss();
		NumericVector bb = fe.beta();
		if (verbose > 1) {
		    Rprintf("step = %g, pwrss0 = %g; pwrss1 = %g\n",
			    step, pwrss0, pwrss1, *bb.begin());	
		    showdbl(bb.begin(), "beta", p);
		}
	    }
	    resp.variance();
	    resp.updateSqrtRWt();
	    crit = compare_vec(&varold[0], var.begin(), varold.size());
	    if (verbose > 1) Rprintf("convergence criterion: %g\n", crit);
	}
	return Laplace();
    } // IRLS

    double glmer::Laplace() {
	NumericVector devRes = resp.devResid();
	return std::accumulate(devRes.begin(), devRes.end(),
			       double()) +
	    re.ldL2() + re.sqLenU();
    }
	
    double glmerDe::PIRLS(int verbose) {
	const NumericVector var = resp.var();
	std::vector<double> varold(var.size());
	NumericVector uu = re.u();
	int q = uu.size();
	std::vector<double> incr(q), ubase(q), newU(q);
	double crit, step, pwrss0, pwrss1;

	crit = 10. * CM_TOL;
	// weights, var and sqrtrwt are assumed established.
	for (int i = 0; crit >= CM_TOL && i < CM_MAXITER; i++) {
				// store copies of var and u
	    std::copy(uu.begin(), uu.end(), ubase.begin());
	    std::copy(var.begin(), var.end(), varold.begin());
	    updateGamma();	// linear predictor
	    resp.linkInv();	// mu
	    pwrss0 = resp.updateWrss() + re.sqLenU();
	    resp.updateSqrtXWt();      // muEta and sqrtXwt
				// solve for incr
	    re.rwUpdateL(resp.sqrtXwt(), resp.wtres(), uu, &incr[0]);
	    pwrss1 = pwrss0;	// force one evaluation of the loop
	    for (step = 1.; pwrss0 <= pwrss1 && step > CM_SMIN; step /= 2.) {
		std::transform(incr.begin(), incr.end(), newU.begin(),
			       std::bind2nd(std::multiplies<double>(), step));
		std::transform(ubase.begin(), ubase.end(), newU.begin(), newU.begin(),
			       std::plus<double>());
		re.setU(&newU[0]);
		updateGamma();
		resp.linkInv();
		pwrss1 = resp.updateWrss() + re.sqLenU();
		if (verbose > 1) {
		    Rprintf("step = %g, pwrss0 = %g; pwrss1 = %g\n",
			    step, pwrss0, pwrss1, uu[0]);
		    showdbl(uu.begin(), "u", uu.size());
		}
	    }
	    resp.variance();
	    resp.updateSqrtRWt();
	    crit = compare_vec(&varold[0], var.begin(), varold.size());
	    if (verbose > 1) Rprintf("convergence criterion: %g\n", crit);
	}
	return Laplace();
    } // PIRLS

    nlmer::nlmer(S4 xp) :
	re(S4(xp.slot("re"))),
	resp(S4(xp.slot("resp"))) {
    }

    nlmerDe::nlmerDe(S4 xp) :
	nlmer(xp),
	fe(S4(xp.slot("fe"))) {
    }

    void nlmerDe::updateMu() {
	const NumericVector u = re.u(), offset = resp.offset();
	NumericVector b(u.size()), gamma(offset.size());
	std::copy(offset.begin(), offset.end(), gamma.begin());
	chmDn bb(b), gg(gamma);

	showdbl(gamma.begin(), "gamma after offset", gamma.size());
	re.Lambda().dmult('N', 1., 0., chmDn(u), bb);
	showdbl(u.begin(), "u", u.size());
	showdbl(b.begin(), "b", b.size());
	re.Zt().dmult('T', 1., 1., bb, gg);
	showdbl(gamma.begin(), "gamma after reUpdate", gamma.size());
	fe.X().dmult('N', 1., 1., chmDn(fe.beta()), gg);
	showdbl(gamma.begin(), "gamma after feUpdate", gamma.size());
	resp.updateMu(gamma);
    }

} // namespace mer

RCPP_FUNCTION_2(double,lmerDeUpdate,S4 xp,NumericVector nt) {
    mer::lmer<mer::lmerDeFeMod> lm(xp);
    return lm.updateTheta(nt);
}

RCPP_FUNCTION_2(double,lmerSpUpdate,S4 xp,NumericVector nt) {
    mer::lmer<mer::lmerSpFeMod> lm(xp);
    return lm.updateTheta(nt);
}

RCPP_FUNCTION_VOID_2(reModUpdate,S4 xp,NumericVector nt) {
    mer::reModule re(xp);
    re.updateTheta(nt);
}

RCPP_FUNCTION_1(double,lmerDeDeviance,S4 xp) {
    mer::lmer<mer::lmerDeFeMod> lm(xp);
    return lm.deviance();
}

RCPP_FUNCTION_1(double,lmerSpDeviance,S4 xp) {
    mer::lmer<mer::lmerSpFeMod> lm(xp);
    return lm.deviance();
}

RCPP_FUNCTION_2(double,glmerDeIRLS,S4 xp,IntegerVector verb) {
    mer::glmerDe glmr(xp);
    return glmr.IRLS(verb[0]);
}

RCPP_FUNCTION_2(double,glmerDePIRLS,S4 xp,IntegerVector verb) {
    mer::glmerDe glmr(xp);
    return glmr.PIRLS(verb[0]);
}

RCPP_FUNCTION_1(double,glmerLaplace,S4 xp) {
    mer::glmer glmr(xp);
    return glmr.Laplace();
}

RCPP_FUNCTION_VOID_2(feModSetBeta,S4 xp,NumericVector nbeta) {
    mer::feModule fe(xp);
    fe.setBeta(nbeta);
}

RCPP_FUNCTION_VOID_1(nlmerDeEval,S4 xp) {
    mer::nlmerDe nlmr(xp);
    Rprintf("past initialization of nlmerDe object\n");
    nlmr.updateMu();
}
