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


/** 
 * Update the destination dense matrix as
 * solve(L, solve(L, Lambda %*% src, system = "P"), system = "L")
 * 
 * @param Lambda a sparse q by q matrix
 * @param L a q by q sparse Cholesky factor
 * @param src source dense matrix
 * @param dest destination dense matrix
 */
    void DupdateL(chmSp const &Lambda, chmFr const &L,
		  chmDn const &src, chmDn &dest) {
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
 * @param Lambda a sparse q by q matrix
 * @param L a q by q sparse Cholesky factor
 * @param src source sparse matrix
 * @return as above
 */
    CHM_SP SupdateL(chmSp const &Lambda, chmFr const &L,
		    chmSp const &src) {
	CHM_SP t1 = Lambda.crossprod(src);
	CHM_SP t2 = L.spsolve(CHOLMOD_P, t1);
	::M_cholmod_free_sparse(&t1, &c);
	t1 = L.spsolve(CHOLMOD_L, t2);
	::M_cholmod_free_sparse(&t2, &c);
	return t1;
    }
    
    reModule::reModule(S4 xp)
	: d_L(S4(xp.slot("L"))),
	  d_Lambda(S4(xp.slot("Lambda"))),
	  d_Ut(S4(xp.slot("Ut"))),
	  d_Zt(S4(xp.slot("Zt"))),
	  Lind(xp.slot("Lind")),
	  lower(xp.slot("lower")),
	  theta(xp.slot("theta")),
	  d_Utr(xp.slot("Utr")),
	  d_cu(d_Utr.size()),
	  d_u(xp.slot("u")),
	  d_ldL2(NumericVector(xp.slot("ldL2")).begin()),
	  d_sqlLenU(std::inner_product(d_u.begin(), d_u.end(),
				       d_u.begin(), double())) {
    }
    
    void reModule::setU(const double *nu) {
	std::copy(nu, nu + d_u.size(), d_u.begin());
	d_sqlLenU = std::inner_product(d_u.begin(), d_u.end(),
				       d_u.begin(), double());
    }
    
/** 
 * Check and install new value of theta.  
 *
 * Update Lambda using theta and Lind.  Update Ut and L.
 * Evaluate ldL2.  Evaluate cu from Utr.
 * 
 * @param nt New value of theta
 */
    void reModule::updateTheta(NumericVector const &nt) {
	if (nt.size() != theta.size())
	    ::Rf_error("length(theta) = %d != length(newtheta) = %d",
		       theta.size(), nt.size());
//FIXME: Check values of nt against the lower bounds in lower.
//Perhaps use the Rcpp::any_if algorithm

	double *Lamx = (double*)d_Lambda.x, *th = theta.begin();
	std::copy(nt.begin(), nt.end(), th);
	int *Li = Lind.begin(), Lis = Lind.size();
	for (int i = 0; i < Lis; i++) Lamx[i] = th[Li[i] - 1];

	CHM_SP LamTrZt = d_Lambda.crossprod(d_Zt);
	d_Ut.update(*LamTrZt);
	M_cholmod_free_sparse(&LamTrZt, &c);
	d_L.update(d_Ut, 1.);
	*d_ldL2 = M_chm_factor_ldetL2(&d_L);

	chmDn ccu(d_cu);
	DupdateL(d_Lambda, d_L, chmDn(d_Utr), ccu);
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
 * Solve for u given the updated cu
 * 
 * @param cu 
 */
    void reModule::updateU(chmDn const &cu) {
	CHM_DN t1 = d_L.solve(CHOLMOD_Lt, &cu);
	CHM_DN t2 = d_L.solve(CHOLMOD_Pt, t1);
	::M_cholmod_free_dense(&t1, &c);
	setU((double*)t2->x);
	::M_cholmod_free_dense(&t2, &c);
    }

    merResp::merResp(S4 xp)
	: d_wrss(NumericVector(xp.slot("wrss")).begin()),
	  d_offset(xp.slot("offset")),
	  d_sqrtrwt(xp.slot("sqrtrwt")),
	  d_wtres(xp.slot("wtres")),
	  mu(xp.slot("mu")),
	  weights(xp.slot("weights")),
	  y(xp.slot("y"))
    {
	int n = y.size(), os = d_offset.size();

	if (mu.size() != n || d_wtres.size() != n ||
	    weights.size() != n || d_sqrtrwt.size() != n)
	    ::Rf_error("y, mu, sqrtrwt, wtres and weights slots must have equal lengths");
	if (os < 1 || os % n)
	    ::Rf_error("length(offset) must be a positive multiple of length(y)");
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

    double merResp::updateMu(NumericVector const &gamma) {
	std::copy(gamma.begin(), gamma.end(), mu.begin());
	return updateWrss();
    }

    double glmerResp::updateMu(NumericVector const &gamma) {
	std::copy(gamma.begin(), gamma.end(), eta.begin());
	linkInv();
	return updateWrss();
    }

    double nlmerResp::updateMu(NumericVector const &gamma) {
	int n = y.size();
	double *gg = gamma.begin();

	for (int p = 0; p < pnames.size(); p++) {
	    std::string pn(pnames[p]);
	    NumericVector pp = nlenv.get(pn);
	    std::copy(gg + n * p, gg + n * (p + 1), pp.begin());
	}
	NumericVector rr = nlmod.eval(SEXP(nlenv));
	if (rr.size() != n)
	    Rf_error("Length mu = %d, expected %d", rr.size(), n);
	std::copy(rr.begin(), rr.end(), mu.begin());
	NumericMatrix rrg = rr.attr("gradient");
	std::copy(rrg.begin(), rrg.end(), d_sqrtXwt.begin());
	return updateWrss();
    }

    lmerResp::lmerResp(S4 xp)
	: merResp(xp),
	  d_reml(*LogicalVector(xp.slot("REML")).begin()) {
    }

    rwResp::rwResp(S4 xp)
	: merResp(xp),
	  d_sqrtXwt(SEXP(xp.slot("sqrtXwt"))) {
    }

    glmerResp::glmerResp(S4 xp)
	: rwResp(xp),
	  family(SEXP(xp.slot("family"))),
	  eta(xp.slot("eta")),
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
	  d_Vtr(xp.slot("Vtr")),
	  d_ldRX2(NumericVector(xp.slot("ldRX2")).begin()) {
    }

    void feModule::setBeta(NumericVector const &nbeta) {
	if (nbeta.size() != d_beta.size())
	    Rf_error("length(beta) = %d should be %d",
		     nbeta.size(), d_beta.size());
	std::copy(nbeta.begin(), nbeta.end(), d_beta.begin());
    }

    deFeMod::deFeMod(S4 xp)
	: feModule(xp),
	  d_X(S4(xp.slot("X"))),
	  d_RZX(S4(xp.slot("RZX"))),
	  d_UtV(S4(xp.slot("UtV"))),
	  d_V(S4(xp.slot("V"))),
	  d_VtV(S4(xp.slot("VtV"))),
	  d_RX(S4(xp.slot("RX"))) {
    }

/** 
 * Update RZX and RX
 *   RZX <<- solve(L, solve(L, crossprod(Lambda, ZtX), 
 *                          sys = "P"),
 *                 sys = "L")
 *   RX <<- chol(XtX - crossprod(RZX))
 * 
 * @param Lambda 
 * @param L
 */
    void deFeMod::updateRzxRx(chmSp const &Lambda, chmFr const &L) {
	chmDn cRZX(d_RZX);
	DupdateL(Lambda, L, chmDn(d_UtV), cRZX);
	d_RX.update('T', -1., d_RZX, 1., d_VtV);
	*d_ldRX2 = d_RX.logDet2();
    }

/** 
 * Update Vtr
 */
    void deFeMod::updateVtr(NumericVector const &wtres) {
	d_V.dgemv('T', 1., wtres, 0., d_Vtr);
    }

    void deFeMod::updateV(NumericMatrix const &sqrtXwt) {
	int nc = sqrtXwt.ncol();
	if (nc != 1) Rf_error("code for nlmer not yet written");
// FIXME: Rewrite this using std::transform
	std::copy(d_X.x.begin(), d_X.x.end(), d_V.x.begin());
	int m = d_X.nrow(), n = d_X.ncol();
	double *vv = d_V.x.begin(), *ww = sqrtXwt.begin();
	for (int j = 0; j < n; j++)
	    for (int i = 0; i < m; i++) vv[i + j * m] *= ww[i];
    }
    
    void deFeMod::updateUtV(chmSp const &Ut) {
	chmDn cUtV(d_UtV);
	Ut.dmult('N', 1., 0., chmDn(d_V), cUtV);
    }

 /** 
  * Update beta
  *	beta <- solve(RX, solve(t(RX), Vtr - crossprod(RZX, cu)))
  * 
  * @param cu intermediate solution of random-effects.  Updated to
  *     cu - RZX %*% beta on return
  */
    NumericVector deFeMod::updateBeta(NumericVector const &cu) {
	NumericVector ans(cu.size());
	std::copy(cu.begin(), cu.end(), ans.begin());
	std::copy(d_Vtr.begin(), d_Vtr.end(), d_beta.begin());
	
	d_RZX.dgemv('T', -1., ans, 1., d_beta);
	d_RX.dpotrs(d_beta);
	d_RZX.dgemv('N', -1., d_beta, 1., ans);
	return ans;
    }

    spFeMod::spFeMod(Rcpp::S4 xp)
	: feModule(xp),
	  d_RZX(S4(xp.slot("RZX"))),
	  d_UtV(S4(xp.slot("UtV"))),
	  d_V(S4(xp.slot("V"))),
	  d_VtV(S4(xp.slot("VtV"))),
	  d_X(S4(xp.slot("X"))),
	  d_RX(S4(xp.slot("RX"))) {
    }

/** 
 * Update RZX and RX
 *   RZX <<- solve(L, solve(L, crossprod(Lambda, ZtX), 
 *                          sys = "P"),
 *                 sys = "L")
 *   RX <- chol(XtX - crossprod(RZX))
 * 
 * @param re a random-effects module
 */
    void spFeMod::updateRzxRx(chmSp const &Lambda, chmFr const &L) {
	double mone[] = {-1.,0}, one[] = {1.,0};
	CHM_SP t2 = SupdateL(Lambda, L, d_UtV);

	d_RZX.update(*t2);
	M_cholmod_free_sparse(&t2, &c);	
	
	CHM_SP t1 = d_RZX.crossprod();
	t2 = M_cholmod_add(&d_VtV, t1, one, mone, 1/*values*/,
			     1/*sorted*/, &c);
	M_cholmod_free_sparse(&t1, &c);
	d_RX.update(*t2);

	M_cholmod_free_sparse(&t2, &c);

	*d_ldRX2 = M_chm_factor_ldetL2(&d_RX);
    }

 /** 
  * Update beta
  * 	resp@cbeta <- Vtr - crossprod(RZX, cu)
  *	beta <- solve(RX, solve(t(RX), resp@cbeta))
  *	resp@cu <- resp@cu - RZX %*% beta
  * 
  * @param resp response module
  */
    NumericVector spFeMod::updateBeta(NumericVector const &cu) {
	NumericVector ans(cu.size());
	std::copy(cu.begin(), cu.end(), ans.begin());
	std::copy(d_Vtr.begin(), d_Vtr.end(), d_beta.begin());
	chmDn cbeta(d_beta), cans(ans);
	d_RZX.dmult('T', -1., 1., chmDn(cu), cbeta);
	CHM_DN t1 = d_RX.solve(CHOLMOD_A, &cbeta);
	double *t1b = (double*)t1->x;
	std::copy(t1b,  t1b + d_beta.size(), d_beta.begin());
	M_cholmod_free_dense(&t1, &c);
	d_RZX.dmult('N', -1., 1., cbeta, cans);
	return ans;
    }

    static inline double sqrtquotient(double x, double y) {
	return sqrt(x/y);
    }

    void glmerResp::updateSqrtRWt() {
	variance();
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

    double glmerDe::updateMu() {
	const NumericVector u = re.u(), offset = resp.offset();
	NumericVector b(u.size()), gamma(offset.size());
	std::copy(offset.begin(), offset.end(), gamma.begin());
	chmDn bb(b), gg(gamma);

	re.Lambda().dmult('N', 1., 0., chmDn(u), bb);
	re.Zt().dmult('T', 1., 1., bb, gg);
	fe.X().dmult('N', 1., 1., chmDn(fe.beta()), gg);
	return resp.updateMu(gamma);
    }

    void deFeMod::updateRX(bool useRZX) {
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
	const NumericVector var = resp.var(), u = re.u();
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
				// using current beta and u
	    pwrss0 = updateMu() + re.sqrLenU();
	    resp.updateSqrtXWt();	// muEta and sqrtXwt
	    fe.updateV(resp.sqrtXwt());	// derive V from X
	    fe.updateVtr(resp.wtres());	// crossprod(V, wtres)
	    fe.updateRX(false);		// factor V'V
	    std::copy(fe.Vtr().begin(), fe.Vtr().end(), incr.begin());
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
		pwrss1 = updateMu() + re.sqrLenU();
		if (verbose > 1) {
		    NumericVector bb = fe.beta();
		    Rprintf("step = %g, pwrss0 = %g; pwrss1 = %g\n",
			    step, pwrss0, pwrss1, *bb.begin());	
		    showdbl(bb.begin(), "beta", p);
		}
	    }
	    resp.updateSqrtRWt();
	    crit = compare_vec(&varold[0], var.begin(), varold.size());
	    if (verbose > 1) Rprintf("convergence criterion: %g\n", crit);
	}
	return Laplace();
    } // IRLS

    double glmer::Laplace() {
	const NumericVector devRes = resp.devResid();
	return re.ldL2() + re.sqrLenU() +
	    std::accumulate(devRes.begin(), devRes.end(), double());
    }
	
    double glmerDe::PIRLS(int verbose) {
	const NumericVector var = resp.var(), u = re.u();
	std::vector<double> varold(var.size());
	int q = u.size();
	NumericVector incr(q), ubase(q), newU(q);
	double crit, step, pwrss0, pwrss1;

	crit = 10. * CM_TOL;
	// weights, var and sqrtrwt are assumed established.
	for (int i = 0; crit >= CM_TOL && i < CM_MAXITER; i++) {
				// store copies of var and u
	    std::copy(u.begin(), u.end(), ubase.begin());
	    std::copy(var.begin(), var.end(), varold.begin());
	    updateMu();	// linear predictor
	    pwrss0 = resp.updateWrss() + re.sqrLenU();
	    resp.updateSqrtXWt();      // muEta and sqrtXwt
				// solve for incr
	    re.rwUpdateL(resp.sqrtXwt(), resp.wtres(), u, &incr[0]);
	    pwrss1 = pwrss0;	// force one evaluation of the loop
	    for (step = 1.; pwrss0 <= pwrss1 && step > CM_SMIN; step /= 2.) {
		std::transform(incr.begin(), incr.end(), newU.begin(),
			       std::bind2nd(std::multiplies<double>(),
					    step));
		std::transform(ubase.begin(), ubase.end(),
			       newU.begin(), newU.begin(),
			       std::plus<double>());
		re.setU(newU.begin());
		updateMu();
		pwrss1 = resp.updateWrss() + re.sqrLenU();
		if (verbose > 1) {
		    Rprintf("step = %g, pwrss0 = %g; pwrss1 = %g\n",
			    step, pwrss0, pwrss1);
		    showdbl(u.begin(), "u", u.size());
		}
	    }
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

    double nlmerDe::updateMu() {
	const NumericVector u = re.u(), offset = resp.offset();
	NumericVector b(u.size()), gamma(offset.size());
	std::copy(offset.begin(), offset.end(), gamma.begin());
	chmDn bb(b), gg(gamma);

	re.Lambda().dmult('N', 1., 0., chmDn(u), bb);
	re.Zt().dmult('T', 1., 1., bb, gg);
	fe.X().dmult('N', 1., 1., chmDn(fe.beta()), gg);
	return resp.updateMu(gamma);
    }

} // namespace mer

RCPP_FUNCTION_2(double,lmerDeUpdate,S4 xp,NumericVector nt) {
    mer::mer<mer::deFeMod,mer::lmerResp> lm(xp);
    lm.updateTheta(nt);
    lm.solveBetaU();
    lm.updateMu();
    return lm.reml() ? lm.profREML() : lm.profDev();
}

RCPP_FUNCTION_VOID_1(lmerDeUpdateMu,S4 xp) {
    mer::mer<mer::deFeMod,mer::lmerResp> lm(xp);
    lm.updateMu();
}

RCPP_FUNCTION_2(double,lmerSpUpdate,S4 xp,NumericVector nt) {
    mer::mer<mer::spFeMod,mer::lmerResp> lm(xp);
    lm.updateTheta(nt);
    lm.solveBetaU();
    lm.updateMu();
    return lm.reml() ? lm.profREML() : lm.profDev();
}

RCPP_FUNCTION_VOID_1(lmerSpUpdateMu,S4 xp) {
    mer::mer<mer::spFeMod,mer::lmerResp> lm(xp);
    lm.updateMu();
}

RCPP_FUNCTION_VOID_2(reModUpdate,S4 xp,NumericVector nt) {
    mer::reModule re(xp);
    re.updateTheta(nt);
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

RCPP_FUNCTION_1(double,nlmerDeEval,S4 xp) {
    mer::nlmerDe nlmr(xp);
    return nlmr.updateMu();
}
