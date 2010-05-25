#include "mer.h"

using namespace Rcpp;
using namespace MatrixNs;

namespace mer{

    void showdbl(const double* x, const char* nm, int n) {
	if (n < 1) {
	    Rprintf("%s[%d]:\n", nm, n);
	    return;
	}
	int n5 = (n < 5) ? n : 5;
	Rprintf("%s[1:%d]: %g", nm, n, x[0]);
	for (int i = 1; i < n5; i++) Rprintf(", %g", x[i]);
	if (n > 5) Rprintf(", ...");
	Rprintf("\n");
    }

    void showincr(double step, double c0, double c1,
		  NumericVector const& incr, const char* nm) {
	Rprintf("step = %8.5f, pwrss0 = %12g, pwrss1 = %12g\n",
		step, c0, c1);
	showdbl(incr.begin(), nm, incr.size());
    }

    void showint(const int* x, const char* nm, int n) {
	if (n < 1) {
	    Rprintf("%s[%d]:\n", nm, n);
	    return;
	}
	int n20 = (n < 20) ? n : 20;
	Rprintf("%s[1:%d]: %d", nm, n, x[0]);
	for (int i = 1; i < n20; i++) Rprintf(", %d", x[i]);
	if (n > 20) Rprintf(", ...");
	Rprintf("\n");
    }

    void showCHM_DN(const CHM_DN x, const std::string &nm) {
	Rprintf("%s: nrow = %d, ncol = %d, nzmax = %d, d = %d, xtype = %d, dtype = %d\n",
		nm.c_str(), x->nrow, x->ncol, x->nzmax, x->d, x->xtype, x->dtype);
	showdbl((double*)x->x, "x", x->nzmax);
    }

    void showCHM_FR(const CHM_FR x, const std::string &nm) {
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

    void showCHM_SP(const CHM_SP x, const std::string &nm) {
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

/** 
 * Determine the weighted Euclidean distance between two vectors,
 * relative to the square root of the product of their lengths
 * 
 * @param v1 First vector to compare
 * @param v2 Second vector to compare
 * @param wt square root of the weights
 * 
 * @return relative difference between the matrices
 */
    double compareVecWt(NumericVector const& v1,
			NumericVector const& v2,
			NumericVector const& wt) {
	int n = v1.size();
	double num, d1, d2;
	if (v2.size() != n || wt.size() != n)
	    Rf_error("%s: size mismatch, %d != %d or != %d\n",
		     "compareVecWt", n, v2.size(), wt.size());
	std::vector<double> a(n);

	std::transform(v1.begin(), v1.end(), v2.begin(),
		       a.begin(), std::minus<double>());
	std::transform(a.begin(), a.end(), wt.begin(), a.begin(),
		       std::multiplies<double>());
	num = std::inner_product(a.begin(), a.end(), a.begin(), double());

	std::transform(v1.begin(), v1.end(), wt.begin(), a.begin(),
		       std::multiplies<double>());
	d1 = std::inner_product(a.begin(), a.end(), a.begin(), double());
	std::transform(v2.begin(), v2.end(), wt.begin(), a.begin(),
		       std::multiplies<double>());
	d2 = std::inner_product(a.begin(), a.end(), a.begin(), double());

	if (d1 == 0) {
	    return (d2 == 0) ? 0. : sqrt(num/d2);
	}
	return (d2 == 0) ? sqrt(num/d1) : sqrt(num/sqrt(d1 * d2));
    }

/** 
 * Update the destination dense matrix as
 * solve(L, solve(L, Lambda %*% src, system = "P"), system = "L")
 * 
 * @param Lambda a sparse q by q matrix
 * @param L a q by q sparse Cholesky factor
 * @param src source dense matrix
 * @param dest destination dense matrix
 */
    void DupdateL(chmSp const& Lambda,
		  chmFr const&      L,
		  chmDn const&    src,
		  chmDn          dest) {
	Lambda.dmult('T', 1., 0., src, dest);
	NumericMatrix
	    ans = L.solve(CHOLMOD_L, L.solve(CHOLMOD_P, &dest));
	std::copy(ans.begin(), ans.end(), (double*)dest.x);
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
    CHM_SP SupdateL(chmSp const& Lambda,
		    chmFr const&      L,
		    chmSp const&    src) {
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
	  d_Lind(xp.slot("Lind")),
	  d_Utr(xp.slot("Utr")),
	  d_cu(d_Utr.size()),
	  d_lower(xp.slot("lower")),
	  d_theta(xp.slot("theta")),
	  d_u(xp.slot("u")),
	  d_ldL2(NumericVector(xp.slot("ldL2")).begin()),
	  d_sqlLenU(std::inner_product(d_u.begin(), d_u.end(),
				       d_u.begin(), double())) {
    }
    
    void reModule::setU(NumericVector const &ubase,
			NumericVector const &incr, double step) {
	int q = d_u.size();
	if (ubase.size() != d_u.size())
	    Rf_error("%s: expected %s.size() = %d, got %d",
		     "reModule::setU", "ubase", q, ubase.size());
	if (step == 0.) {
	    std::copy(ubase.begin(), ubase.end(), d_u.begin());
	} else {
	    if (incr.size() != q)
		Rf_error("%s: expected %s.size() = %d, got %d",
			 "reModule::setU", "incr", q, incr.size());
	    std::transform(incr.begin(), incr.end(), d_u.begin(),
			   std::bind2nd(std::multiplies<double>(), step));
	    std::transform(ubase.begin(), ubase.end(), d_u.begin(),
			   d_u.begin(), std::plus<double>());
	}
	d_sqlLenU = std::inner_product(d_u.begin(), d_u.end(),
				       d_u.begin(), double());
    }
    
/** 
 * Check and install new value of theta.  Update Lambda.
 *
 * @param nt New value of theta.
 */
    void reModule::updateTheta(NumericVector const& nt) {
	if (nt.size() != d_theta.size())
	    Rf_error("%s: %s[1:%d], expected [1:%d]",
		     "updateTheta", "newtheta",
		     nt.size(), d_theta.size());
//FIXME: Check values of nt against the lower bounds in lower.
//Perhaps use the Rcpp::any_if algorithm
//	NumericVector diffs(p);
//	std::transform(nt.begin(), nt.end(), d_lower.begin(),
				// store new theta
	std::copy(nt.begin(), nt.end(), d_theta.begin());
				// update Lambda from theta and Lind
	double *Lamx = (double*)d_Lambda.x, *th = d_theta.begin();
	int *Li = d_Lind.begin(), Lis = d_Lind.size();
	for (int i = 0; i < Lis; i++) Lamx[i] = th[Li[i] - 1];

	updateLcu();
    }
/*
 * Update Ut from Zt and sqrtXwt.
 * Update Utr from wtres and Ut.
 */ 
    void reModule::reweight(NumericMatrix const&   Xwt,
			    NumericVector const& wtres) {
	int Wnc = Xwt.ncol()
//	    ,Wnr = Xwt.nrow()
//	    ,Znc = d_Zt.ncol
//	    ,Znr = d_Zt.nrow
	    ;
	if (Wnc == 1) {
	    d_Ut.update(d_Zt);
	    d_Ut.scale(CHOLMOD_COL, chmDn(Xwt));
	} else Rf_error("Multiple columns in Xwt");
	chmDn cUtr(d_Utr);
	d_Ut.dmult('N', 1., 0., chmDn(wtres), cUtr);
    }

/*
 * Update L and ldL2.
 * Update cu from Utr and L.
 */ 
    void reModule::updateLcu() {
	double one = 1.;
	CHM_SP LtUt = d_Lambda.crossprod(d_Ut);
	M_cholmod_factorize_p(LtUt, &one, (int*)NULL, (size_t)0,
			      &d_L, &c);
	*d_ldL2 = M_chm_factor_ldetL2(&d_L);
	M_cholmod_free_sparse(&LtUt, &c);

	chmDn ccu(d_cu);
	DupdateL(d_Lambda, d_L, chmDn(d_Utr), d_cu);
    }

    NumericVector reModule::UIncr() {
	NumericMatrix ans = d_L.solve(CHOLMOD_A, d_Utr);
	return NumericVector(SEXP(ans));
    }

/** 
 * Solve for u given the updated cu
 * 
 * @param cu 
 */
    void reModule::updateU(chmDn const &cu) {
	NumericMatrix t1 = d_L.solve(CHOLMOD_Lt, &cu);
	const chmDn ct1(t1);
	setU(d_L.solve(CHOLMOD_Pt, &ct1));
    }

    merResp::merResp(S4 xp)
	: d_wrss(NumericVector(xp.slot("wrss")).begin()),
	  d_offset(xp.slot("offset")),
	  d_sqrtrwt(xp.slot("sqrtrwt")),
	  d_wtres(xp.slot("wtres")),
	  d_mu(xp.slot("mu")),
	  d_weights(xp.slot("weights")),
	  d_y(xp.slot("y")),
	  d_sqrtXwt(SEXP(xp.slot("sqrtXwt"))) {
	int n = d_y.size(), os = d_offset.size();

	if (d_mu.size() != n || d_wtres.size() != n ||
	    d_weights.size() != n || d_sqrtrwt.size() != n)
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
	std::transform(d_y.begin(), d_y.end(), d_mu.begin(),
		       d_wtres.begin(), std::minus<double>());
				// wtres <- wtres * sqrtrwt
	std::transform(d_wtres.begin(), d_wtres.end(), d_sqrtrwt.begin(),
		       d_wtres.begin(), std::multiplies<double>());
	*d_wrss = std::inner_product(d_wtres.begin(), d_wtres.end(),
				     d_wtres.begin(), double());
	return *d_wrss;
    }

    lmerResp::lmerResp(S4 xp)
	: merResp(xp),
	  d_reml(*IntegerVector(xp.slot("REML")).begin()) {
	std::fill(d_mu.begin(), d_mu.end(), double());
	updateWrss();
    }

    double lmerResp::Laplace(double ldL2, double ldRX2, double sqrL)const{
	double lnum = 2.* PI * (*d_wrss + sqrL), n = (double)d_y.size();
	if (d_reml == 0) return ldL2 + n * (1 + log(lnum / n));
	double nmp = n - d_reml;
	return ldL2 + ldRX2 + nmp * (1 + log(lnum / nmp));
    }

    double lmerResp::updateMu(NumericVector const &gamma) {
	std::copy(gamma.begin(), gamma.end(), d_mu.begin());
	return updateWrss();
    }

    glmerResp::glmerResp(S4 xp)
	: merResp(xp),
	  d_devres(NumericVector(xp.slot("devres")).begin()),
	  family(SEXP(xp.slot("family"))),
	  d_eta(xp.slot("eta")),
	  d_muEta(xp.slot("muEta")),
	  d_n(xp.slot("n")),
	  d_var(xp.slot("var")) {
	updateWts();
    }

    double glmerResp::Laplace(double  ldL2,
			      double ldRX2,
			      double  sqrL) const{
	NumericVector
	    dr = const_cast<glmFamily&>(family).devResid(     d_mu,
							 d_weights,
							       d_y);
	*d_devres = std::accumulate(dr.begin(), dr.end(), double());
	return ldL2 + sqrL + *d_devres;

    }
	
    double glmerResp::updateMu(NumericVector const &gamma) {
	std::copy(gamma.begin(), gamma.end(), d_eta.begin());
	family.linkInv(d_mu, d_eta);
	return updateWrss();
    }

    static inline double sqrtquot(double x, double y) {
	return sqrt(x / y);
    }

    double glmerResp::updateWts() {
	family.variance(d_var, d_mu);
	std::transform(d_weights.begin(), d_weights.end(),
		       d_var.begin(), d_sqrtrwt.begin(), sqrtquot);
	family.muEta(d_muEta, d_eta);
	std::transform(d_sqrtrwt.begin(), d_sqrtrwt.end(),
		       d_muEta.begin(), d_sqrtXwt.begin(),
		       std::multiplies<double>());
	return updateWrss();
    }

    nlmerResp::nlmerResp(S4 xp)
	: merResp(xp),
	  nlenv(SEXP(xp.slot("nlenv"))),
	  nlmod(SEXP(xp.slot("nlmod"))),
	  pnames(xp.slot("pnames")) {
    }

    double nlmerResp::Laplace(double  ldL2,
			      double ldRX2,
			      double  sqrL) const{
	double lnum = 2.* PI * (*d_wrss + sqrL),
	    n = (double)d_y.size();
	return ldL2 + n * (1 + log(lnum / n));
    }

    double nlmerResp::updateMu(NumericVector const &gamma) {
	int n = d_y.size();
	double *gg = gamma.begin();

	for (int p = 0; p < pnames.size(); p++) {
	    std::string pn(pnames[p]);
	    NumericVector pp = nlenv.get(pn);
	    std::copy(gg + n * p, gg + n * (p + 1), pp.begin());
	}
	NumericVector rr = nlmod.eval(SEXP(nlenv));
	if (rr.size() != n)
	    Rf_error("Length mu = %d, expected %d", rr.size(), n);
	std::copy(rr.begin(), rr.end(), d_mu.begin());
	NumericMatrix rrg = rr.attr("gradient");
	std::copy(rrg.begin(), rrg.end(), d_sqrtXwt.begin());
	return updateWrss();
    }

    feModule::feModule(S4 xp)
	: d_beta(xp.slot("beta")),
	  d_Vtr(xp.slot("Vtr")),
	  d_ldRX2(NumericVector(xp.slot("ldRX2")).begin()) {
    }

    void feModule::setBeta(NumericVector const &bbase,
			   NumericVector const &incr, double step) {
	int p = d_beta.size();
	if (bbase.size() != d_beta.size())
	    Rf_error("%s: expected %s.size() = %d, got %d",
		     "feModule::setBeta", "bbase", p, bbase.size());
	if (step == 0.) {
	    std::copy(bbase.begin(), bbase.end(), d_beta.begin());
	} else {
	    if (incr.size() != p)
		Rf_error("%s: expected %s.size() = %d, got %d",
			 "feModule::setBeta", "incr", p, incr.size());
	    std::transform(incr.begin(), incr.end(), d_beta.begin(),
			   std::bind2nd(std::multiplies<double>(), step));
	    std::transform(bbase.begin(), bbase.end(), d_beta.begin(),
			   d_beta.begin(), std::plus<double>());
	}
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

    void deFeMod::reweight(chmSp         const&    Ut,
			   NumericMatrix const&   Xwt,
			   NumericVector const& wtres) {
	if (Xwt.size() != d_X.nrow())
	    Rf_error("%s: dimension mismatch %s(%d,%d), %s(%d,%d)",
		     "deFeMod::reweight", "X", d_X.nrow(), d_X.ncol(),
		     "Xwt", Xwt.nrow(), Xwt.ncol());
	int Wnc = Xwt.ncol(), Wnr = Xwt.nrow(),
	    Xnc = d_X.ncol(), Xnr = d_X.nrow();
	double *V = d_V.x.begin(), *W = Xwt.begin(), *X = d_X.x.begin();
	NumericVector pr(Wnr);
				// Use the first column of W to initialize V
	for (int j = 0; j < Xnc; j++)
	    std::transform(X + j*Xnr, X + (j+1)*Xnr, W,
			   V + j * Wnr, std::multiplies<double>());
				// increment by other columns, if any
	for (int Wj = 1; Wj < Wnc; Wj++) {
	    double *Xpt = X + Wj * Wnr;
	    for (int j = 0; j < Xnc; j++) {
		std::transform(Xpt + j * Xnr, Xpt + (j+1)*Xnr, W,
			       pr.begin(), std::multiplies<double>());
		std::transform(pr.begin(), pr.end(),
			       V + j*Wnr, V + j*Wnr, std::plus<double>());
	    }
	}
	d_V.dgemv('T', 1., wtres, 0., d_Vtr);
	chmDn cUtV(d_UtV);
	Ut.dmult('N', 1., 0., chmDn(d_V), cUtV);
    }

    NumericVector deFeMod::BetaIncr() {
	d_RX.update(d_V);
	NumericVector ans(d_beta.size());
	std::copy(d_Vtr.begin(), d_Vtr.end(), ans.begin());
	d_RX.dpotrs(ans);
	return ans;
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
	NumericMatrix t1 = d_RX.solve(CHOLMOD_A, &cbeta);
	std::copy(t1.begin(),  t1.end(), d_beta.begin());
	d_RZX.dmult('N', -1., 1., cbeta, cans);
	return ans;
    }

    void spFeMod::reweight(chmSp         const&      Ut,
			   NumericMatrix const& sqrtXwt,
			   NumericVector const&   wtres) {
//FIXME: Code here not written yet.
	// updateV(sqrtXwt);
	// updateVtr(wtres);
	// d_RX.update(d_V);
    }

    NumericVector spFeMod::rwBetaIncr(NumericMatrix const& sqrtXwt,
			              NumericVector const&   wtres) {
	chmDn cbeta(d_beta);
	NumericMatrix ans = d_RX.solve(CHOLMOD_A, &cbeta);
	return NumericVector(SEXP(ans));
    }

} // namespace mer

RCPP_FUNCTION_2(double, lmerDeUpdate, S4 xp, NumericVector nt) {
    mer::mer<mer::deFeMod,mer::lmerResp> lm(xp);
    lm.updateTheta(nt);
    lm.solveBetaU();
    lm.updateMu();
    return lm.Laplace();
}

RCPP_FUNCTION_2(double, lmerSpUpdate, S4 xp, NumericVector nt) {
    mer::mer<mer::spFeMod,mer::lmerResp> lm(xp);
    lm.updateTheta(nt);
    lm.solveBetaU();
    lm.updateMu();
    return lm.Laplace();
}

RCPP_FUNCTION_2(double, glmerDeIRLS, S4 xp, int verb) {
    mer::mer<mer::deFeMod,mer::glmerResp> glmr(xp);
    return glmr.IRLS(verb);
}

RCPP_FUNCTION_3(double, glmerDePIRLS, S4 xp, NumericVector nt, int verb) {
    mer::mer<mer::deFeMod,mer::glmerResp> glmr(xp);
    glmr.updateTheta(nt);
    return glmr.PIRLS(verb);
}

RCPP_FUNCTION_VOID_1(glmerDeUpdateRzxRx,S4 xp) {
    mer::mer<mer::deFeMod,mer::glmerResp> glmr(xp);
    return glmr.updateRzxRx();
}

RCPP_FUNCTION_VOID_2(feSetBeta, S4 xp, NumericVector nbeta) {
    mer::feModule fe(xp);
    fe.setBeta(nbeta);
}

RCPP_FUNCTION_2(double, nlmerDeEval, S4 xp, int verb) {
    mer::mer<mer::deFeMod,mer::nlmerResp> nlmr(xp);
    return nlmr.IRLS(verb);
}

RCPP_FUNCTION_1(double,nlmerDeIRLS,S4 xp) {
    mer::mer<mer::deFeMod,mer::nlmerResp> nlmr(xp);
    return nlmr.updateMu();
}
