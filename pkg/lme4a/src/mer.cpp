#include "mer.h"

using namespace Rcpp;
using namespace MatrixNs;

static void showdbl(const double* x, const char* nm, int n) {
    int n5 = (n < 5) ? n : 5;
    Rprintf("%s: %g", nm, x[0]);
    for (int i = 1; i < n5; i++) Rprintf(", %g", x[i]);
    Rprintf("\n");
}

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

#if 0
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
#endif

static void showCHM_SP(const CHM_SP x, const std::string &nm) {
    Rprintf("%s: nrow = %d, ncol = %d, xtype = %d, stype = %d, itype = %d, dtype = %d\n",
	    nm.c_str(), x->nrow, x->ncol, x->xtype, x->stype,
	    x->itype, x->dtype);
    int nc = x->ncol, nnz = ::M_cholmod_nnz(x, &c);
    showint((int*)x->p, "p", nc + 1);
    if (nnz > 0) {
	showint((int*)x->i, "i", nnz);
	showdbl((double*)x->x, "x", nnz);
    } else Rprintf("nnz = 0\n");
}

namespace mer{

    reModule::reModule(S4 xp) :
	L(S4(xp.slot("L"))),
	Lambda(S4(xp.slot("Lambda"))),
	Ut(S4(xp.slot("Ut"))),
	Zt(S4(xp.slot("Zt"))),
	Lind(xp.slot("Lind")),
	lower(xp.slot("lower")),
	theta(xp.slot("theta")),
	d_ldL2(NumericVector(xp.slot("ldL2")).begin()),
	d_u(xp.slot("u"))
    {}
    
    void reModule::setU(const double *nu) {
	std::copy(nu, nu + d_u.size(), d_u.begin());
    }
    
    const NumericVector &reModule::u() const {return d_u;}

/** 
 * Lambda@x[] <- theta[Lind]; Ut <- crossprod(Lambda, Zt);
 * update(L,Ut,1); ldL2
 * 
 * @param nt New value of theta
 */
    void reModule::updateTheta(const NumericVector &nt) {
	if (nt.size() != theta.size())
	    ::Rf_error("length(theta) = %d != length(newtheta) = %d",
		       theta.size(), nt.size());
	std::copy(nt.begin(), nt.end(), theta.begin());
	double *Lamx = (double*)Lambda.x, *nnt = nt.begin();
	int *Li = Lind.begin();
	for (int i = 0; i < Lind.size(); i++) Lamx[i] = nnt[Li[i] - 1];

	CHM_SP LamTrZt = Lambda.crossprod(Zt);
	Ut.update(*LamTrZt);
	M_cholmod_free_sparse(&LamTrZt, &c);
	L.update(Ut, 1.);
	*d_ldL2 = M_chm_factor_ldetL2(&L);
    }

    void reModule::rwUpdateL(dgeMatrix const &Xwt,
			     NumericVector const &wtres,
			     double *ans) {
	double one = 1., zero = 0.;
	CHM_SP wtUt = M_cholmod_copy_sparse(&Ut, &c);
	showCHM_SP(wtUt, "Ut");
	chmDn cXwt(Xwt);
	showCHM_DN(&cXwt, "Xwt");
	if (Xwt.ncol() == 1) {
	    M_cholmod_scale(&cXwt, CHOLMOD_COL, wtUt, &c);
	}
	showCHM_SP(wtUt, "wtUt");
	M_cholmod_factorize_p(wtUt, &one, (int*)NULL, (size_t)0,
			      &L, &c);
	*d_ldL2 = M_chm_factor_ldetL2(&L);
	chmDn cans(ans, Ut.nrow, 1), cwtres(wtres);
	M_cholmod_sdmult(wtUt, 0/*trans*/, &one, &zero, &cwtres,
			 &cans, &c);
	M_cholmod_free_sparse(&wtUt, &c);
	CHM_DN incr = M_cholmod_solve(CHOLMOD_A, &L, &cans, &c);
	double *incx = (double*)incr->x;
	std::copy(incx, incx + Ut.nrow, ans);
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
	CHM_DN t1 = L.solve(CHOLMOD_Lt, &cu);
	CHM_DN t2 = L.solve(CHOLMOD_Pt, t1);
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
    void reModule::incGamma(NumericVector &gam) const {
	chmDn cgam(gam);	// gamma <- gamma + crossprod(Ut, uu)
	Ut.dmult('T', 1., 1., chmDn(d_u), cgam);
    }
	
/** 
 * Update the destination dense matrix as
 * solve(L, solve(L, Lambda %*% src, system = "P"), system = "L")
 * 
 * @param src source dense matrix
 * @param dest destination dense matrix
 */
    void reModule::DupdateL(chmDn const &src, chmDn &dest) const {
	Lambda.dmult('T', 1., 0., src, dest);
	CHM_DN t1 = L.solve(CHOLMOD_P, &dest);
	CHM_DN t2 = L.solve(CHOLMOD_L, t1);
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
	CHM_SP t1 = Lambda.crossprod(src);
	CHM_SP t2 = L.spsolve(CHOLMOD_P, t1);
	::M_cholmod_free_sparse(&t1, &c);
	t1 = L.spsolve(CHOLMOD_L, t2);
	::M_cholmod_free_sparse(&t2, &c);
	return t1;
    }
    
    double reModule::ldL2() const { return *d_ldL2; }

    merResp::merResp(S4 xp) :
	Utr(SEXP(xp.slot("Utr"))),
	Vtr(SEXP(xp.slot("Vtr"))),
	cbeta(SEXP(xp.slot("cbeta"))),
	cu(SEXP(xp.slot("cu"))),
	mu(SEXP(xp.slot("mu"))),
	offset(SEXP(xp.slot("offset"))),
	wtres(SEXP(xp.slot("wtres"))),
	weights(SEXP(xp.slot("weights"))),
	y(SEXP(xp.slot("y"))),
	cUtr(Utr),
	ccu(cu)
    {
	NumericVector wrssVec(SEXP(xp.slot("wrss")));
	wrss = wrssVec.begin();

	int n = y.size(), os = offset.size();

	if (mu.size() != n || wtres.size() != n ||
	    weights.size() != n)
	    ::Rf_error("y, mu, wtres and weights slots must have equal lengths");
	if (os < 1 || os % n)
	    ::Rf_error("length(offset) must be a positive multiple of length(y)");
	if (cbeta.size() != Vtr.size())
	    ::Rf_error("cbeta and Vtr slots must have equal lengths");
	if (cu.size() != Utr.size())
	    ::Rf_error("cu and Utr slots must have equal lengths");
    }

/** 
 * Update cu using Lambda and L
 *  cu <- solve(L, solve(L, crossprod(Lambda, Utr),
 *                       sys = "P"), 
 *              sys = "L")
 * 
 * @param re random-effects module
 */
    void merResp::updateL(reModule const &re) {
	re.DupdateL(cUtr, ccu);
    }

    void merResp::updateMu(NumericVector const &gamma) {
	if (gamma.size() != mu.size())
	    throw std::range_error("dimensions of gamma and mu");
	std::copy(gamma.begin(), gamma.end(), mu.begin());
    }

/** 
 * Update the wtres vector and return the sum of squares
 *   wtres <- sqrtrwts * (y - mu)
 *   return(wrss <- sum(wtres^2))
 *
 * @return Updated weighted residual sum of squares
 */
// FIXME: Move the sqrtrwt vector to the merResp class and use the
// general form here too
    double merResp::updateWrss() {
	int n = y.size();
	double *mm = mu.begin(),
	    *rr = wtres.begin(),
	    *ww = weights.begin(),
	    *yy = y.begin(),
	    ans = 0;
	for (int i = 0; i < n; i++) {
	    rr[i] = yy[i] - mm[i];
	    ans += rr[i] * rr[i] * ww[i];
	}
	*wrss = ans;
	return ans;
    }

    rwResp::rwResp(S4 xp) :
	merResp(xp),
	gamma(SEXP(xp.slot("gamma"))),
	sqrtrwt(SEXP(xp.slot("sqrtrwt"))),
	sqrtXwt(Rcpp::S4(SEXP(xp.slot("sqrtXwt")))) {
    }

// FIXME: This should be the method for merResp.  Move sqrtrwt to the merResp class    
    double rwResp::updateWrss() {
	std::transform(y.begin(), y.end(), mu.begin(), // wtres <- y - mu
		       wtres.begin(), std::minus<double>());
				// wtres <- wtres * sqrtrwt
	std::transform(wtres.begin(), wtres.end(), sqrtrwt.begin(),
		       wtres.begin(), std::multiplies<double>());
	*wrss = std::inner_product(wtres.begin(), wtres.end(),
				   wtres.begin(), double());
	return *wrss;
    }

    feModule::feModule(S4 xp) :
	d_beta(xp.slot("beta")),
	d_ldRX2(xp.slot("ldRX2"))
    { }

    double feModule::ldRX2() const {
	return *d_ldRX2.begin();
    }

    void feModule::setBeta(const double *nb) {
	std::copy(nb, nb + d_beta.size(), d_beta.begin());
    }
#if 0    
    void feModule::setBeta(NumericVector const &BB, double mm) {
	if (bb.size() != BB.size())
	    throw std::range_error("setBeta dimension mismatch");
	std::transform(BB.begin(), BB.end(), bb.begin(),
		       std::bind2nd(std::multiplies<double>(), mm));
    }
    
    void feModule::setBeta(double vv) {
	std::fill(bb.begin(), bb.end(), vv);
    }
#endif    
    const NumericVector& feModule::beta() const {
	return d_beta;
    }

    deFeMod::deFeMod(S4 xp) :
	feModule(xp),
	X(S4(xp.slot("X"))),
	RZX(S4(xp.slot("RZX"))),
	d_RX(S4(xp.slot("RX")))
    { }

    const MatrixNs::Cholesky &deFeMod::RX() const {
	return d_RX;
    }

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
	chmDn cRZX(RZX);
	re.DupdateL(chmDn(ZtX), cRZX);
	d_RX.update('T', -1., RZX, 1., XtX);
	*d_ldRX2.begin() = d_RX.logDet2();
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
	RZX.dgemv('T', -1., resp.cu, 1., resp.cbeta);
	std::copy(resp.cbeta.begin(), resp.cbeta.end(), d_beta.begin());
	d_RX.dpotrs(d_beta);
	RZX.dgemv('N', -1., d_beta, 1., resp.cu);
    }

/** 
 *  gamma += crossprod(Ut, u)
 * 
 * @param gam Value of gamma to be updated
 */
    void deFeMod::incGamma(NumericVector &gam) const {
	X.dgemv('N', 1., d_beta, 1., gam);
    }
    
    rwDeFeMod::rwDeFeMod(S4 xp) :
	deFeMod(xp),
	V(Rcpp::S4(SEXP(xp.slot("V"))))
    { }

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

	*d_ldRX2.begin() = M_chm_factor_ldetL2(&RX);
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
	chmDn ccbeta(resp.cbeta);
	RZX.dmult('T', -1., 1., resp.ccu, ccbeta);
	std::copy(resp.cbeta.begin(), resp.cbeta.end(), d_beta.begin());
	chmDn cbeta(d_beta);
	CHM_DN t1 = RX.solve(CHOLMOD_A, &cbeta);
	double *t1b = (double*)t1->x;
	std::copy(t1b,  t1b + d_beta.size(), d_beta.begin());
	M_cholmod_free_dense(&t1, &c);
	RZX.dmult('N', -1., 1., cbeta, resp.ccu);
    }

    void spFeMod::incGamma(NumericVector &gam) const {
	chmDn bb(d_beta), gg(gam);
	X.dmult('N', 1., 1., bb, gg);
    }

    void glmerResp::updateSqrtRWt() {
	std::transform(weights.begin(), weights.end(), var.begin(),
		       sqrtrwt.begin(), std::divides<double>());
	std::transform(sqrtrwt.begin(), sqrtrwt.end(),
		       sqrtrwt.begin(), sqrt);
    }

    void glmerResp::updateSqrtXWt() {
	MuEta();
	std::transform(sqrtrwt.begin(), sqrtrwt.end(), muEta.begin(),
		       sqrtXwt.x.begin(), std::multiplies<double>());
    }

    glmerDe::glmerDe(Rcpp::S4 xp) :
	glmer(xp),
	fe(Rcpp::S4(SEXP(xp.slot("fe")))) {
    }

    void glmerDe::updateGamma() {
	std::copy(resp.offset.begin(), resp.offset.end(),
		  resp.gamma.begin());
	fe.incGamma(resp.gamma);
	re.incGamma(resp.gamma);
    }

    void rwDeFeMod::updateV(rwResp const &resp) {
	int nc = resp.sqrtXwt.ncol();
	if (nc != 1) Rf_error("code for nlmer not yet written");
	std::copy(X.x.begin(), X.x.end(), V.x.begin());
	int m = X.nrow(), n = X.ncol();
	double *vv = V.x.begin(), *ww = resp.sqrtXwt.x.begin();
	for (int j = 0; j < n; j++)
	    for (int i = 0; i < m; i++) vv[i + j * m] *= ww[i];
    }

    void rwDeFeMod::updateRX(bool useRZX) {
	if (useRZX) throw std::range_error("should not happen");
	d_RX.update(V);
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

    double glmerDe::IRLS() {
	std::vector<double> varold(resp.var.size());
	NumericVector bb = fe.beta();
	int p = bb.size();
	std::vector<double> betabase(p), newbeta(p), incr(p); 
	double crit, step, pwrss0, pwrss1;

	crit = 10. * CM_TOL;
	    // weights, var and sqrtrwt are assumed established.
	for (int i = 0; crit >= CM_TOL && i < CM_MAXITER; i++) {
				// store copies of var and beta
	    std::copy(bb.begin(), bb.end(), betabase.begin());
	    std::copy(resp.var.begin(), resp.var.end(), varold.begin());
	    updateGamma();     // using current beta
	    resp.linkInv();    // mu
	    pwrss0 = resp.updateWrss() + re.sqLenU();
	    resp.updateSqrtXWt();      // muEta and sqrtXwt
	    fe.updateV(resp);	       // derive V from X
				// Vtr <- crossprod(V, wtres)
	    fe.V.dgemv('T', 1., resp.wtres, 0., resp.Vtr);
	    fe.updateRX(false);	// factor V'V
	    std::copy(resp.Vtr.begin(), resp.Vtr.end(), incr.begin());
	    fe.RX().dpotrs(incr); // evaluate the increment
	    pwrss1 = pwrss0;	// force one evaluation of the loop
	    for (step = 1.; pwrss0 <= pwrss1 && step > CM_SMIN; step /= 2.) {
				// newbeta <- betabase + incr * mult
		std::transform(incr.begin(), incr.end(), newbeta.begin(),
		       std::bind2nd(std::multiplies<double>(), step));
		std::transform(betabase.begin(), betabase.end(),
			       newbeta.begin(), newbeta.begin(),
			       std::plus<double>());
		fe.setBeta(&newbeta[0]);
		updateGamma();
		resp.linkInv();
		pwrss1 = resp.updateWrss();
		NumericVector bb = fe.beta();
		Rprintf("step = %g, pwrss0 = %g; pwrss1 = %g\n",
			step, pwrss0, pwrss1, *bb.begin());
		showdbl(bb.begin(), "u", bb.size());
	    }
	    resp.variance();
	    resp.updateSqrtRWt();
	    crit = compare_vec(&varold[0], &resp.var[0], varold.size());
	    Rprintf("convergence criterion: %g\n", crit);
	}
	NumericVector devRes =
	    resp.family.devResid(resp.mu, resp.weights, resp.y);
	return std::accumulate(devRes.begin(), devRes.end(), double());
    } // IRLS

    double glmerDe::PIRLS() {
	std::vector<double> varold(resp.var.size());
	NumericVector uu = re.u();
	int q = uu.size();
	std::vector<double> incr(q), ubase(q), newU(q);
	double crit, step, pwrss0, pwrss1;

	crit = 10. * CM_TOL;
	// weights, var and sqrtrwt are assumed established.
	for (int i = 0; crit >= CM_TOL && i < CM_MAXITER; i++) {
				// store copies of var and u
	    std::copy(uu.begin(), uu.end(), ubase.begin());
	    std::copy(resp.var.begin(), resp.var.end(), varold.begin());
	    updateGamma();	// linear predictor
	    resp.linkInv();	// mu
	    pwrss0 = resp.updateWrss() + re.sqLenU();
	    resp.updateSqrtXWt();      // muEta and sqrtXwt
	    re.rwUpdateL(resp.sqrtXwt, resp.wtres, &incr[0]);
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
		Rprintf("step = %g, pwrss0 = %g; pwrss1 = %g\n",
			step, pwrss0, pwrss1, uu[0]);
		showdbl(uu.begin(), "u", uu.size());
	    }
	    resp.variance();
	    resp.updateSqrtRWt();
	    crit = compare_vec(&varold[0], &resp.var[0], varold.size());
	    Rprintf("convergence criterion: %g\n", crit);
	}
	NumericVector devRes =
	    resp.family.devResid(resp.mu, resp.weights, resp.y);
	return std::accumulate(devRes.begin(), devRes.end(), double());
    } // PIRLS

} // namespace mer

RCPP_FUNCTION_2(double,lmerDeUpdate, S4 xp,NumericVector nt) {
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

RCPP_FUNCTION_1(double,glmerDeIRLS,S4 xp) {
    mer::glmerDe glmr(xp);
    return glmr.IRLS();
}

RCPP_FUNCTION_1(double,glmerDePIRLS,S4 xp) {
    mer::glmerDe glmr(xp);
    return glmr.PIRLS();
}
