#include "mer.h"

using namespace Rcpp;
using namespace MatrixNs;

namespace mer{

    void showdbl(const double* x, const char* nm, int n) {
	int n5 = (n < 5) ? n : 5;
	Rprintf("%s: %g", nm, x[0]);
	for (int i = 1; i < n5; i++) Rprintf(", %g", x[i]);
	Rprintf("\n");
    }

    void showint(const int* x, const char* nm, int n) {
	int n20 = (n < 20) ? n : 20;
	Rprintf("%s: %d", nm, x[0]);
	for (int i = 1; i < n20; i++) Rprintf(", %d", x[i]);
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
 * Determine the Euclidean distance between two vectors relative to the
 * square root of the product of their lengths
 * 
 * @param v1 First vector to compare
 * @param v2 Second vector to compare
 * 
 * @return relative difference between the vectors
 */
    double compare_vec(const NumericVector v1, const NumericVector v2) {
	int n = v1.size();
	if (v2.size() != n)
	    Rf_error("%s: size mismatch, %d != %d\n",
		     "compare_vec", n, v2.size());
	std::vector<double> diff(v1.size());
	std::transform(v1.begin(), v1.end(), v2.begin(),
		       diff.begin(), std::minus<double>());
	double d1 = sqrt(std::inner_product(v1.begin(), v1.end(),
					    v1.begin(), double())),
	    d2 = sqrt(std::inner_product(v1.begin(), v1.end(),
					 v1.begin(), double())),
	    num = sqrt(std::inner_product(diff.begin(), diff.end(),
					  diff.begin(), double()));
	if (d1 == 0) {
	    if (d2 == 0) return num;
	    else return num/d2;
	}
	if (d2 == 0) return num/d1;
	return num/sqrt(d1 * d2);
    }

    template<typename T> 
    static T myClone(T object) {
	T ans(object.size());
	std::copy(object.begin(), object.end(), ans.begin());
	return ans;
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
    
    void reModule::updateSqrLenU() {
	d_sqlLenU = std::inner_product(d_u.begin(), d_u.end(),
				       d_u.begin(), double());
    }

    void reModule::setU(const cholmod_dense* nu) {
	int q = d_u.size();
	if (nu->nrow != (size_t)q || nu->ncol != (size_t)1)
	    Rf_error("%s: expected dimensions (%d,1), got (%d,%d)",
		     "reModule::setU(chmDn)", q, nu->nrow, nu->ncol);
	double *uu = (double*)nu->x;
	std::copy(uu, uu + q, d_u.begin());
	updateSqrLenU();
    }
    
    void reModule::setU(NumericVector const &nu) {
	if (nu.size() != d_u.size())
	    Rf_error("%s: expected %d, got %d",
		     "reModule::setU(NumericVector)", d_u.size(), nu.size());
	std::copy(nu.begin(), nu.end(), d_u.begin());
	updateSqrLenU();
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

    NumericVector reModule::rwUpdateL(NumericMatrix const &Xwt,
				      NumericVector const &wtres) {
	double one = 1., mone = -1;
	CHM_SP wtUt = M_cholmod_copy_sparse(&d_Ut, &c);
	const chmDn cXwt(Xwt);
	if (d_Ut.ncol == d_Zt.ncol) {
	    M_cholmod_scale(&cXwt, CHOLMOD_COL, wtUt, &c);
	} else Rf_error("Multiple columns in Xwt");
	M_cholmod_factorize_p(wtUt, &one, (int*)NULL, (size_t)0,
			      &d_L, &c);
	*d_ldL2 = M_chm_factor_ldetL2(&d_L);
	NumericVector ans = myClone(d_u);
	chmDn cans(ans);
	const chmDn cwtres(wtres);
	M_cholmod_sdmult(wtUt, 0/*trans*/, &one, &mone, &cwtres,
			 &cans, &c);
	M_cholmod_free_sparse(&wtUt, &c);
	CHM_DN incr = M_cholmod_solve(CHOLMOD_A, &d_L, &cans, &c);
	double *incx = (double*)incr->x;
	std::copy(incx, incx + d_Ut.nrow, ans.begin());
	M_cholmod_free_dense(&incr, &c);
	return ans;
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
	setU(t2);
	::M_cholmod_free_dense(&t2, &c);
    }

    merResp::merResp(S4 xp)
	: d_wrss(NumericVector(xp.slot("wrss")).begin()),
	  d_offset(xp.slot("offset")),
	  d_sqrtrwt(xp.slot("sqrtrwt")),
	  d_wtres(xp.slot("wtres")),
	  mu(xp.slot("mu")),
	  weights(xp.slot("weights")),
	  y(xp.slot("y")),
	  d_sqrtXwt(SEXP(xp.slot("sqrtXwt"))) {
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

    NumericVector glmerResp::devResid() {
	return family.devResid(mu, weights, y);
    }

    double glmerResp::unscaledDeviance() {
	const NumericVector dr = devResid();
	return std::accumulate(dr.begin(), dr.end(), double());
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

    glmerResp::glmerResp(S4 xp)
	: merResp(xp),
	  family(SEXP(xp.slot("family"))),
	  eta(xp.slot("eta")),
	  muEta(xp.slot("muEta")),
	  n(xp.slot("n")),
	  d_var(xp.slot("var")) {
    }

    nlmerResp::nlmerResp(S4 xp)
	: merResp(xp),
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

    void deFeMod::updateRX(bool useRZX) {
	if (useRZX) throw std::range_error("should not happen");
	d_RX.update(d_V);
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

} // namespace mer

RCPP_FUNCTION_2(double,lmerDeUpdate,S4 xp,NumericVector nt) {
    mer::mer<mer::deFeMod,mer::lmerResp> lm(xp);
    lm.updateTheta(nt);
    lm.solveBetaU();
    lm.updateMu();
    return lm.reml() ? lm.profREML() : lm.profDev();
}

RCPP_FUNCTION_2(double,lmerSpUpdate,S4 xp,NumericVector nt) {
    mer::mer<mer::spFeMod,mer::lmerResp> lm(xp);
    lm.updateTheta(nt);
    lm.solveBetaU();
    lm.updateMu();
    return lm.reml() ? lm.profREML() : lm.profDev();
}

RCPP_FUNCTION_VOID_2(reUpdateTheta,S4 xp,NumericVector nt) {
    mer::reModule re(xp);
    re.updateTheta(nt);
}

RCPP_FUNCTION_2(double,glmerDeIRLS,S4 xp,IntegerVector verb) {
    mer::mer<mer::deFeMod,mer::glmerResp> glmr(xp);
    return glmr.IRLS(verb[0]);
}

RCPP_FUNCTION_2(double,glmerDePIRLS,S4 xp,IntegerVector verb) {
    mer::mer<mer::deFeMod,mer::glmerResp> glmr(xp);
    return glmr.PIRLS(verb[0]);
}

RCPP_FUNCTION_1(double,glmerDeLaplace,S4 xp) {
    mer::mer<mer::deFeMod,mer::glmerResp> glmr(xp);
    return glmr.Laplace();
}

RCPP_FUNCTION_VOID_1(glmerDeUpdateRzxRx,S4 xp) {
    mer::mer<mer::deFeMod,mer::glmerResp> glmr(xp);
    return glmr.updateRzxRx();
}

RCPP_FUNCTION_VOID_2(feSetBeta,S4 xp,NumericVector nbeta) {
    mer::feModule fe(xp);
    fe.setBeta(nbeta);
}

RCPP_FUNCTION_1(double,nlmerDeEval,S4 xp) {
    mer::mer<mer::deFeMod,mer::nlmerResp> nlmr(xp);
    return nlmr.updateMu();
}
