#include "mer.h"
#include <R_ext/BLAS.h>

using namespace Rcpp;
using namespace MatrixNs;
using namespace std;

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

    void showCHM_DN(const CHM_DN x, const string &nm) {
	Rprintf("%s: nrow = %d, ncol = %d, nzmax = %d, d = %d, xtype = %d, dtype = %d\n",
		nm.c_str(), x->nrow, x->ncol, x->nzmax, x->d, x->xtype, x->dtype);
	showdbl((double*)x->x, "x", x->nzmax);
    }

    void showCHM_FR(const CHM_FR x, const string &nm) {
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

    void showCHM_SP(const CHM_SP x, const string &nm) {
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
	vector<double> a(n);

	transform(v1.begin(), v1.end(), v2.begin(),
		       a.begin(), minus<double>());
	transform(a.begin(), a.end(), wt.begin(), a.begin(),
		       multiplies<double>());
	num = inner_product(a.begin(), a.end(), a.begin(), double());

	transform(v1.begin(), v1.end(), wt.begin(), a.begin(),
		       multiplies<double>());
	d1 = inner_product(a.begin(), a.end(), a.begin(), double());
	transform(v2.begin(), v2.end(), wt.begin(), a.begin(),
		       multiplies<double>());
	d2 = inner_product(a.begin(), a.end(), a.begin(), double());

	if (d1 == 0) {
	    return (d2 == 0) ? 0. : sqrt(num/d2);
	}
	return (d2 == 0) ? sqrt(num/d1) : sqrt(num/sqrt(d1 * d2));
    }

//FIXME: change d_Ut to a CHM_SP and write a separate utility to take
//sqrtXwt and a sparse matrix.
    reModule::reModule(S4 xp)
	: d_xp(       xp),
	  d_L(     S4(xp.slot("L"))),
	  d_Lambda(S4(xp.slot("Lambda"))),
	  d_Ut(    S4(xp.slot("Ut"))),
	  d_Zt(    S4(xp.slot("Zt"))),
	  d_Lind(     xp.slot("Lind")),
	  d_lower(    xp.slot("lower")),
	  d_u(        xp.slot("u")),
	  d_cu(d_u.size()),
	  d_sqrLenU(inner_product(d_u.begin(), d_u.end(), d_u.begin(), double())) {
    }

    /** 
     * Update L, Ut and cu for new weights.
     *
     * Update Ut from Zt and sqrtXwt, then L from Lambda and Ut
     * Update cu from wtres, Lambda and Ut.
     * 
     * @param Xwt Matrix of weights for the model matrix
     * @param wtres weighted residuals
     */
    void reModule::reweight(Rcpp::NumericMatrix const&   Xwt,
			    Rcpp::NumericVector const& wtres) {
	double mone = -1., one = 1.; 
	int Wnc = Xwt.ncol();
	if (Wnc == 1) {
	    d_Ut.update(d_Zt);	// copy Zt to Ut
	    d_Ut.scale(CHOLMOD_COL, chmDn(Xwt));
	} else {
	    int n = Xwt.nrow();
	    CHM_TR tr = M_cholmod_sparse_to_triplet(&d_Zt, &c);
	    int *j = (int*)tr->j, nnz = tr->nnz;
	    double *x = (double*)tr->x, *W = Xwt.begin();
	    for (int k = 0; k < nnz; k++) {
		x[k] *= W[j[k]];
		j[k] = j[k] % n;
	    }
	    tr->ncol = (size_t)n;

	    CHM_SP sp = M_cholmod_triplet_to_sparse(tr, nnz, &c);
	    M_cholmod_free_triplet(&tr, &c);
	    d_Ut.update(*sp);
	    M_cholmod_free_sparse(&sp, &c);
	}
				// update the factor L
	CHM_SP LambdatUt = d_Lambda.crossprod(d_Ut);
	d_L.update(*LambdatUt, 1.);
	d_ldL2 = d_L.logDet2();
				// update cu
	chmDn ccu(d_cu), cwtres(wtres);
	copy(d_u.begin(), d_u.end(), d_cu.begin());
	M_cholmod_sdmult(LambdatUt, 0/*trans*/, &one, &mone, &cwtres, &ccu, &c);
	M_cholmod_free_sparse(&LambdatUt, &c);
	NumericMatrix
	    ans = d_L.solve(CHOLMOD_L, d_L.solve(CHOLMOD_P, d_cu));
	copy(ans.begin(), ans.end(), d_cu.begin());
    }

    /** 
     * Install a new value of U, either a single vector or as a
     * combination of a base, an increment and a step factor.
     * 
     * @param ubase base value of u
     * @param incr increment relative to the base
     * @param step step fraction
     */
    void reModule::setU(Rcpp::NumericVector const &ubase,
			Rcpp::NumericVector const &incr, double step) {
	int q = d_u.size();
	if (ubase.size() != q)
	    Rf_error("%s: expected %s.size() = %d, got %d",
		     "reModule::setU", "ubase", q, ubase.size());
	if (step == 0.) {
	    copy(ubase.begin(), ubase.end(), d_u.begin());
	} else {
	    if (incr.size() != q)
		Rf_error("%s: expected %s.size() = %d, got %d",
			 "reModule::setU", "incr", q, incr.size());
	    transform(incr.begin(), incr.end(), d_u.begin(),
			   bind2nd(multiplies<double>(), step));
	    transform(ubase.begin(), ubase.end(), d_u.begin(),
			   d_u.begin(), plus<double>());
	}
	d_sqrLenU = inner_product(d_u.begin(), d_u.end(),
				       d_u.begin(), double());
    }

    /** 
     * Solve for u (or the increment for u) only.
     * 
     */  
    void reModule::solveU() {
	NumericMatrix ans = d_L.solve(CHOLMOD_A, d_cu);
	setU(NumericVector(SEXP(ans)));
    }

//    void reModule::updateDcmp(Rcpp::NumericVector &cmp) const {  // Need Matrix_0.999375-42 or later
    void reModule::updateDcmp(Rcpp::NumericVector &cmp) {
	cmp["ldL2"] = d_L.logDet2();
	cmp["ussq"] = inner_product(d_u.begin(), d_u.end(),
				    d_u.begin(), double());
    }

    /** 
     * Check and install new value of theta.  Update Lambda.
     * 
     * @param nt New value of theta
     */
    void reModule::updateLambda(NumericVector const& nt) {
	if (nt.size() != d_lower.size())
	    Rf_error("%s: %s[1:%d], expected [1:%d]",
		     "updateLambda", "newtheta",
		     nt.size(), d_lower.size());
	double *th = nt.begin(), *ll = d_lower.begin();
				// check for a feasible point
	for (int i = 0; i < nt.size(); i++)
	    if (th[i] < ll[i] || !R_finite(th[i]))
		Rf_error("updateLambda: theta not in feasible region");
				// store (a copy of) theta
	d_xp.slot("theta") = clone(nt);
				// update Lambda from theta and Lind
	double *Lamx = (double*)d_Lambda.x;
	int *Li = d_Lind.begin(), Lis = d_Lind.size();
	for (int i = 0; i < Lis; i++) Lamx[i] = th[Li[i] - 1];
    }

    /** 
     * Solve for u given the updated cu
     * 
     * @param cu 
     */
    void reModule::updateU(Rcpp::NumericVector const &cu) {
	NumericMatrix nu = d_L.solve(CHOLMOD_Pt, d_L.solve(CHOLMOD_Lt, cu));
	setU(NumericVector(SEXP(nu)));
    }

    /** 
     * Zero the contents of d_u and d_sqrLenU
     */
    void reModule::zeroU() {
	fill(d_u.begin(), d_u.end(), double());
	d_sqrLenU = 0.;
    }
    
    merResp::merResp(Rcpp::S4 xp)
	: d_xp(xp),
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
	    Rf_error("y, mu, sqrtrwt, wtres and weights slots must have equal lengths");
	if (os < 1 || os % n)
	    Rf_error("length(offset) must be a positive multiple of length(y)");
    }

    /** 
     * Update the wtres vector and return its sum of squares
     *   wtres <- sqrtrwts * (y - mu)
     *   return(wrss <- sum(wtres^2))
     *
     * @return Updated weighted residual sum of squares
     */
    double merResp::updateWrss() {
				// wtres <- y - mu
	transform(d_y.begin(), d_y.end(), d_mu.begin(),
		  d_wtres.begin(), minus<double>());
				// wtres <- wtres * sqrtrwt
	transform(d_wtres.begin(), d_wtres.end(), d_sqrtrwt.begin(),
		  d_wtres.begin(), multiplies<double>());
	d_wrss = inner_product(d_wtres.begin(), d_wtres.end(),
			       d_wtres.begin(), double());
	return d_wrss;
    }

    void merResp::updateDcmp(Rcpp::NumericVector& cmp) const {
	double wrss = inner_product(d_wtres.begin(), d_wtres.end(),
				    d_wtres.begin(), double());
	double n = (double)d_y.size(), ussq = cmp["ussq"];
	double pwrss = wrss + ussq;
	cmp["wrss"] = wrss;
	cmp["pwrss"] = pwrss;
	cmp["sigmaML"] = sqrt(pwrss/n);
    }

    lmerResp::lmerResp(Rcpp::S4 xp)
	: merResp(xp),
	  d_reml(*IntegerVector(xp.slot("REML")).begin()) {
	copy(d_offset.begin(), d_offset.end(), d_mu.begin());
	updateWrss();
    }

    double lmerResp::Laplace(double  ldL2,
			     double ldRX2,
			     double  sqrL) const {
	double lnum = 2.* PI * (d_wrss + sqrL), n = (double)d_y.size();
	if (d_reml == 0) return ldL2 + n * (1. + log(lnum / n));
	double nmp = n - d_reml;
	return ldL2 + ldRX2 + nmp * (1. + log(lnum / nmp));
    }

    double lmerResp::updateMu(Rcpp::NumericVector const &gamma) {
	copy(gamma.begin(), gamma.end(), d_mu.begin());
	return updateWrss();
    }
    
    glmerResp::glmerResp(Rcpp::S4 xp)
	: merResp(    xp),
	  family(SEXP(xp.slot("family"))),
	  d_eta(      xp.slot("eta")),
	  d_n(        xp.slot("n")) {
	updateWts();
    }
    
    Rcpp::NumericVector glmerResp::devResid() const {
	return family.devResid(d_mu, d_weights, d_y);
    }
	
    double glmerResp::Laplace(double  ldL2,
			      double ldRX2,
			      double  sqrL) const{
	NumericVector dr = devResid();
	return ldL2 + sqrL + accumulate(dr.begin(), dr.end(), double());
    }
	
    double glmerResp::updateMu(NumericVector const &gamma) {
	copy(gamma.begin(), gamma.end(), d_eta.begin());
	family.linkInv(d_mu, d_eta);
	return updateWrss();
    }

    static inline double sqrtquot(double x, double y) {
	return sqrt(x / y);
    }

    double glmerResp::updateWts() {
	NumericVector mueta = family.  muEta(d_eta),
	                 vv = family.variance(d_mu);
	transform(d_weights.begin(), d_weights.end(), vv.begin(),
		  d_sqrtrwt.begin(), sqrtquot);
	transform(d_sqrtrwt.begin(), d_sqrtrwt.end(), mueta.begin(),
		  d_sqrtXwt.begin(), multiplies<double>());
	return updateWrss();
    }

    void glmerResp::updateDcmp(Rcpp::NumericVector& cmp) const {
	merResp::updateDcmp(cmp);
	NumericVector dr = devResid();
	cmp["drsum"] = accumulate(dr.begin(), dr.end(), double());
    }
	
    nlmerResp::nlmerResp(S4 xp)
	: merResp(xp),
	  nlenv(SEXP(xp.slot("nlenv"))),
	  nlmod(SEXP(xp.slot("nlmod"))),
	  pnames(SEXP(xp.slot("pnames"))) {
    }

    double nlmerResp::Laplace(double  ldL2,
			      double ldRX2,
			      double  sqrL) const{
	double lnum = 2.* PI * (d_wrss + sqrL),
	    n = (double)d_y.size();
	return ldL2 + n * (1 + log(lnum / n));
    }

    double nlmerResp::updateMu(NumericVector const &gamma) {
	int n = d_y.size();
	double *gg = gamma.begin();

	for (int p = 0; p < pnames.size(); p++) {
	    string pn(pnames[p]);
	    NumericVector pp = nlenv.get(pn);
	    copy(gg + n * p, gg + n * (p + 1), pp.begin());
	}
	NumericVector rr = nlmod.eval(SEXP(nlenv));
	if (rr.size() != n)
	    Rf_error("Length mu = %d, expected %d", rr.size(), n);
	copy(rr.begin(), rr.end(), d_mu.begin());
	NumericMatrix rrg = rr.attr("gradient");
	copy(rrg.begin(), rrg.end(), d_sqrtXwt.begin());
	return updateWrss();
    }

    feModule::feModule(S4 xp)
	: d_beta(xp.slot("beta")),
	  d_Vtr(d_beta.size()) {
    }

    void feModule::setBeta(NumericVector const &bbase,
			   NumericVector const &incr, double step) {
	if (d_beta.size() == 0) return;
	int p = d_beta.size();
	if (bbase.size() != p)
	    Rf_error("%s: expected %s.size() = %d, got %d",
		     "feModule::setBeta", "bbase", p, bbase.size());
	if (step == 0.) {
	    copy(bbase.begin(), bbase.end(), d_beta.begin());
	} else {
	    if (incr.size() != p)
		Rf_error("%s: expected %s.size() = %d, got %d",
			 "feModule::setBeta", "incr", p, incr.size());
	    transform(incr.begin(), incr.end(), d_beta.begin(),
			   bind2nd(multiplies<double>(), step));
	    transform(bbase.begin(), bbase.end(), d_beta.begin(),
			   d_beta.begin(), plus<double>());
	}
    }

    deFeMod::deFeMod(S4 xp, int n)
	: feModule(        xp),
	  d_RZX(        S4(xp.slot("RZX"))),
	  d_X(          S4(xp.slot("X"))),
	  d_RX(         S4(xp.slot("RX"))),
	  d_UtV(d_RZX.nrow(), d_RZX.ncol()),
	  d_V(             n,   d_X.ncol()),
	  d_VtV(             d_beta.size()) {
    }

    /** 
     * Update V, UtV, VtV and Vtr
     * 
     * @param Ut from the reModule
     * @param Xwt square root of the weights for the model matrices
     * @param wtres weighted residuals
     */
    void deFeMod::reweight(MatrixNs::chmSp       const&    Ut,
			   Rcpp::NumericMatrix   const&   Xwt,
			   Rcpp::NumericVector   const& wtres) {
	if (d_beta.size() == 0) return;
	if (Xwt.size() != d_X.nrow())
	    Rf_error("%s: dimension mismatch %s(%d,%d), %s(%d,%d)",
		     "deFeMod::reweight", "X", d_X.nrow(), d_X.ncol(),
		     "Xwt", Xwt.nrow(), Xwt.ncol());
	int Wnc = Xwt.ncol(), Wnr = Xwt.nrow(),
	    Xnc = d_X.ncol(), Xnr = d_X.nrow();
	double *V = d_V.x().begin(), *X = d_X.x().begin();

	if (Wnc == 1) {
	    for (int j = 0; j < Xnc; j++) 
		transform(Xwt.begin(), Xwt.end(), X + j*Xnr,
			       V + j*Xnr, multiplies<double>());
	} else {
	    int i1 = 1;
	    double one = 1., zero = 0.;
	    NumericVector tmp(Xnr), mm(Wnc); 
	    fill(mm.begin(), mm.end(), 1.);
	    for (int j = 0; j < Xnc; j++) {
		transform(Xwt.begin(), Xwt.end(), X + j*Xnr,
			       tmp.begin(), multiplies<double>());
		F77_CALL(dgemv)("N", &Wnr, &Wnc, &one, tmp.begin(),
				&Wnr, mm.begin(), &i1, &zero,
				V + j * Wnr, &i1);
	    }
	}
	d_V.dgemv('T', 1., wtres, 0., d_Vtr);
	chmDn cUtV(d_UtV);
	Ut.dmult('N', 1., 0., chmDn(d_V), cUtV);
	d_VtV.dsyrk(d_V, 1., 0.);
    }

    /** 
     * Solve (V'V)beta = Vtr for beta.
     * 
     */
    void deFeMod::solveBeta() {
	int p = d_beta.size();
	if (p == 0) return;
	MatrixNs::Cholesky chol(d_V);
	copy(d_Vtr.begin(), d_Vtr.end(), d_beta.begin());
	chol.dpotrs(d_beta);
    }

    void deFeMod::updateDcmp(Rcpp::NumericVector& cmp) const {
	cmp["ldRX2"] = d_RX.logDet2();
    }

    /** 
     * Update RZX and RX
     *   RZX <<- solve(L, solve(L, crossprod(Lambda, ZtX), 
     *                          sys = "P"),
     *                 sys = "L")
     *   RX <<- chol(XtX - crossprod(RZX))
     * 
     * @param Lambda relative covariance factor for the random effects
     * @param L sparse Cholesky factor for the random effects
     */
    void deFeMod::updateRzxRx(MatrixNs::chmSp const& Lambda,
			      MatrixNs::chmFr const&      L) {
	if (d_beta.size() == 0) return;
	chmDn cRZX(d_RZX);
	Lambda.dmult('T', 1., 0., chmDn(d_UtV), cRZX);
	NumericMatrix
	    ans = L.solve(CHOLMOD_L, L.solve(CHOLMOD_P, &cRZX));
	d_RZX.setX(ans);
	d_RX.update('T', -1., d_RZX, 1., d_VtV);
	d_ldRX2 = d_RX.logDet2();
    }

    /** 
     * Update beta
     *	beta <- solve(RX, solve(t(RX), Vtr - crossprod(RZX, cu)))
     * 
     * @param cu intermediate solution of random-effects.
     * 
     * @return cu - RZX %*% beta
     */
    Rcpp::NumericVector deFeMod::updateBeta(Rcpp::NumericVector const& cu) {
	NumericVector ans(cu.size());
	copy(cu.begin(), cu.end(), ans.begin());
	if (d_beta.size() == 0) return ans;
	copy(d_Vtr.begin(), d_Vtr.end(), d_beta.begin());
	
	d_RZX.dgemv('T', -1., ans, 1., d_beta);
	d_RX.dpotrs(d_beta);
	d_RZX.dgemv('N', -1., d_beta, 1., ans);
	return ans;
    }

    spFeMod::spFeMod(Rcpp::S4 xp, int n)
	: feModule(    xp),
	  d_RZX(    S4(xp.slot("RZX"))),
	  d_X(      S4(xp.slot("X"))),
	  d_RX(     S4(xp.slot("RX"))) {
    }

    spFeMod::~spFeMod() {
	if (d_VtV) M_cholmod_free_sparse(&d_VtV, &c);
	if (d_V) M_cholmod_free_sparse(&d_V, &c);
	if (d_UtV) M_cholmod_free_sparse(&d_UtV, &c);
    }

    /** 
     * Update beta
     *	beta <- solve(RX, solve(t(RX), Vtr - crossprod(RZX, cu)))
     * 
     * @param cu intermediate solution of random-effects.
     * 
     * @return cu - RZX %*% beta
     */
    NumericVector spFeMod::updateBeta(NumericVector const &cu) {
	NumericVector ans(cu.size());
	copy(cu.begin(), cu.end(), ans.begin());
	copy(d_Vtr.begin(), d_Vtr.end(), d_beta.begin());
	chmDn cbeta(d_beta), cans(ans);
	d_RZX.dmult('T', -1., 1., chmDn(cu), cbeta);
	NumericMatrix t1 = d_RX.solve(CHOLMOD_A, &cbeta);
	copy(t1.begin(),  t1.end(), d_beta.begin());
	d_RZX.dmult('N', -1., 1., cbeta, cans);
	return ans;
    }


    /** 
     * Update RZX and RX
     *   RZX <<- solve(L, solve(L, crossprod(Lambda, ZtX), 
     *                          sys = "P"),
     *                 sys = "L")
     *   RX <<- chol(XtX - crossprod(RZX))
     * 
     * @param Lambda relative covariance factor for the random effects
     * @param L sparse Cholesky factor for the random effects
     */
    void spFeMod::updateRzxRx(chmSp const &Lambda, chmFr const &L) {
	double mone[] = {-1.,0}, one[] = {1.,0};
	CHM_SP t1 = Lambda.crossprod(d_UtV);
	CHM_SP t2 = L.spsolve(CHOLMOD_P, t1);
	M_cholmod_free_sparse(&t1, &c);
	t1 = L.spsolve(CHOLMOD_L, t2);
	M_cholmod_free_sparse(&t2, &c);
	d_RZX.update(*t1);
	M_cholmod_free_sparse(&t1, &c);	
	
	t1 = d_RZX.crossprod();
	t2 = M_cholmod_add(d_VtV, t1, one, mone, 1/*values*/,
			     1/*sorted*/, &c);
	M_cholmod_free_sparse(&t1, &c);
	d_RX.update(*t2);

	M_cholmod_free_sparse(&t2, &c);

	d_ldRX2 = d_RX.logDet2();
    }

    /** 
     * Reweight the feModule.
     *
     * Update V, UtV and Vtr
     * 
     * @param Ut from the reModule
     * @param Xwt square root of the weights for the model matrices
     * @param wtres weighted residuals
     */
    void spFeMod::reweight(chmSp         const&      Ut,
			   NumericMatrix const& sqrtXwt,
			   NumericVector const&   wtres) {
	double one = 1., zero = 0.;
	if (d_beta.size() == 0) return;
	int Wnc = sqrtXwt.ncol(), Wnr = sqrtXwt.nrow(),
	    Xnc = d_X.ncol, Xnr = d_X.nrow;
	if (sqrtXwt.size() != (int)d_X.nrow)
	    Rf_error("%s: dimension mismatch %s(%d,%d), %s(%d,%d)",
		     "deFeMod::reweight", "X", Xnr, Xnc,
		     "Xwt", Wnr, Wnc);
	if (Wnc == 1) {
	    if (d_V) M_cholmod_free_sparse(&d_V, &c);
	    d_V = M_cholmod_copy_sparse(&d_X, &c);
	    chmDn csqrtX(sqrtXwt);
	    M_cholmod_scale(&csqrtX, CHOLMOD_ROW, d_V, &c);
	} else throw runtime_error("spFeMod::reweight: multiple columns in sqrtXwt");
// FIXME rewrite this using the triplet representation
	
	if (d_UtV) M_cholmod_free_sparse(&d_UtV, &c);
	d_UtV = M_cholmod_ssmult(&Ut, d_V, 0/*styp*/,1/*vals*/,1/*srtd*/, &c);

	if (d_VtV) M_cholmod_free_sparse(&d_VtV, &c);
	CHM_SP t1 = M_cholmod_transpose(d_V, 1/*vals*/, &c);
	d_VtV = M_cholmod_ssmult(t1, d_V, 1/*styp*/,1/*vals*/,1/*srtd*/, &c);
	M_cholmod_free_sparse(&t1, &c);

	chmDn cVtr(d_Vtr);
	const chmDn cwtres(wtres);
	M_cholmod_sdmult(d_V, 'T', &one, &zero, &cwtres, &cVtr, &c);
    }

    /** 
     * Solve (V'V)beta = Vtr for beta.
     * 
     */
    void spFeMod::solveBeta() {
	d_RX.update(*d_VtV);
	NumericMatrix ans = d_RX.solve(CHOLMOD_A, d_Vtr);
	copy(ans.begin(), ans.end(), d_beta.begin());
    }

//    void spFeMod::updateDcmp(Rcpp::NumericVector& cmp) const {  // needs Matrix_0.999375-42 or later
    void spFeMod::updateDcmp(Rcpp::NumericVector& cmp) {
	cmp["ldRX2"] = d_RX.logDet2();
    }


} // namespace mer

RCPP_FUNCTION_1(double, LMMdeviance, S4 xp) {
    S4 fe(xp.slot("fe")), resp(xp.slot("resp"));
    if (!resp.is("lmerResp")) 
	throw runtime_error("LMMupdate on non-lmer object");
    if (fe.is("deFeMod")) {
	mer::mer<mer::deFeMod,mer::lmerResp> lm(xp);
	return lm.LMMdeviance();
    } else if (fe.is("spFeMod")) {
	mer::mer<mer::spFeMod,mer::lmerResp> lm(xp);
	return lm.LMMdeviance();
    } else throw runtime_error("fe slot is neither deFeMod nor spFeMod");
}

RCPP_FUNCTION_4(double, PIRLS, S4 xp, NumericVector u0, int verb, int alg) {
    if (alg < 1 || alg > 3) throw range_error("alg must be 1, 2 or 3");
    mer::Alg aa = (alg == 1) ? mer::Beta : ((alg == 2) ? mer::U : mer::BetaU);
    S4 fe(xp.slot("fe")), resp(xp.slot("resp"));
    bool de = fe.is("deFeMod");
    if (!de && !fe.is("spFeMod"))
	throw runtime_error("fe slot is neither deFeMod nor spFeMod");
    if (resp.is("glmerResp")) {
	if (de) {
	    mer::mer<mer::deFeMod,mer::glmerResp> glmr(xp);
	    return glmr.PIRLS(u0, verb, aa);
	} else {
	    mer::mer<mer::spFeMod,mer::glmerResp> glmr(xp);
	    return glmr.PIRLS(u0, verb, aa);
	}
    } else if (resp.is("nlmerResp")) {
	if (de) {
	    mer::mer<mer::deFeMod,mer::nlmerResp> nlmr(xp);
	    return nlmr.PIRLS(u0, verb, aa);
	} else {
	    mer::mer<mer::spFeMod,mer::nlmerResp> nlmr(xp);
	    return nlmr.PIRLS(u0, verb, aa);
	}
    }
    throw runtime_error("resp slot is not glmerResp or nlmerResp in PIRLS");
}

RCPP_FUNCTION_VOID_2(feSetBeta, S4 xp, NumericVector nbeta) {
    mer::feModule fe(xp);
    fe.setBeta(nbeta);
}

RCPP_FUNCTION_VOID_2(reUpdateLambda, S4 xp, NumericVector nth) {
    mer::reModule re(xp);
    re.updateLambda(nth);
}

RCPP_FUNCTION_VOID_1(updateRzxRx, S4 xp) {
    S4 fe(xp.slot("fe")), resp(xp.slot("resp"));
    bool de = fe.is("deFeMod");
    if (!de && !fe.is("spFeMod"))
	throw runtime_error("fe slot is neither deFeMod nor spFeMod");
    if (resp.is("glmerResp")) {
	if (de) {
	    mer::mer<mer::deFeMod,mer::glmerResp> glmr(xp);
	    glmr.updateRzxRx();
	} else {
	    mer::mer<mer::spFeMod,mer::glmerResp> glmr(xp);
	    glmr.updateRzxRx();
	}
    } else if (resp.is("nlmerResp")) {
	if (de) {
	    mer::mer<mer::deFeMod,mer::nlmerResp> nlmr(xp);
	    nlmr.updateRzxRx();
	} else {
	    mer::mer<mer::spFeMod,mer::nlmerResp> nlmr(xp);
	    nlmr.updateRzxRx();
	}
    } else 
	throw runtime_error("resp slot is not glmerResp or nlmerResp in updateRzxRx");
}

RCPP_FUNCTION_VOID_1(updateDc, S4 xp) {
    List ll(xp.slot("devcomp"));
    NumericVector cmp = ll["cmp"];
    IntegerVector dims = ll["dims"];
    S4 fe(xp.slot("fe")), re(xp.slot("re"));
    mer::merResp resp(S4(xp.slot("resp")));
    mer::reModule reM(re);
    reM.updateDcmp(cmp); // with Matrix_0.999375-42 or later can do this in one step
//    mer::reModule(re).updateDcmp(cmp); // need Matrix_0.999375-42 or later
    resp.updateDcmp(cmp);
    int n = resp.mu().size();
    if (fe.is("deFeMod")) mer::deFeMod(fe, n).updateDcmp(cmp);
    if (fe.is("spFeMod")) // mer::spFeMod(fe, n).updateDcmp(cmp);   // needs Matrix_0.999375-42 or later
    {
	mer::spFeMod spFe(fe, n);
	spFe.updateDcmp(cmp);
    }

    cmp["sigmaREML"] = cmp["sigmaML"] * sqrt((double)dims["n"]/(double)dims["nmp"]);

    ll["cmp"] = cmp;		// should this be necessary?
}
    
