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
 * Determine the Euclidean distance between two vectors relative to the
 * square root of the product of their lengths
 * 
 * @param v1 First vector to compare
 * @param v2 Second vector to compare
 * 
 * @return relative difference between the vectors
 */
    double compare_vec(NumericVector const& v1, NumericVector const& v2) {
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
	NumericMatrix t1 = L.solve(CHOLMOD_P, &dest);
	chmDn ct1(t1);
	NumericMatrix t2 = L.solve(CHOLMOD_L, &ct1);
	std::copy(t2.begin(), t2.end(), (double*)dest.x);
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

    NumericVector reModule::rwUIncr(NumericMatrix const&   Xwt,
				    NumericVector const& wtres) {
	double one = 1., mone = -1;
	CHM_SP wtUt = M_cholmod_copy_sparse(&d_Ut, &c);
	const chmDn cXwt(Xwt);
	if (d_Ut.ncol == d_Zt.ncol) {
	    M_cholmod_scale(&cXwt, CHOLMOD_COL, wtUt, &c);
	} else Rf_error("Multiple columns in Xwt");
	M_cholmod_factorize_p(wtUt, &one, (int*)NULL, (size_t)0,
			      &d_L, &c);
	*d_ldL2 = M_chm_factor_ldetL2(&d_L);
	NumericVector ans = clone(d_u);
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
	  family(SEXP(xp.slot("family"))),
	  d_eta(xp.slot("eta")),
	  d_muEta(xp.slot("muEta")),
	  d_n(xp.slot("n")),
	  d_var(xp.slot("var")) {
	updateWts();
    }

    double glmerResp::Laplace(double ldL2, double ldRX2, double sqrL) const{
	const NumericVector dr =
	    const_cast<glmFamily&>(family).devResid(d_mu, d_weights, d_y);
	return ldL2 + sqrL +
	    std::accumulate(dr.begin(), dr.end(), double());
    }
	
    double glmerResp::updateMu(NumericVector const &gamma) {
	std::copy(gamma.begin(), gamma.end(), d_eta.begin());
	family.linkInv(d_mu, d_eta);
	return updateWrss();
    }

    static inline double sqrtquot(double x, double y){return sqrt(x / y);}

    void glmerResp::updateWts() {
	family.variance(d_var, d_mu);
	std::transform(d_weights.begin(), d_weights.end(),
		       d_var.begin(), d_sqrtrwt.begin(), sqrtquot);
	family.muEta(d_muEta, d_eta);
	std::transform(d_sqrtrwt.begin(), d_sqrtrwt.end(),
		       d_muEta.begin(), d_sqrtXwt.begin(),
		       std::multiplies<double>());
    }

    nlmerResp::nlmerResp(S4 xp)
	: merResp(xp),
	  nlenv(SEXP(xp.slot("nlenv"))),
	  nlmod(SEXP(xp.slot("nlmod"))),
	  pnames(xp.slot("pnames")) {
    }

    double nlmerResp::Laplace(double ldL2, double ldRX2, double sqrL)const{
	double lnum = 2.* PI * (*d_wrss + sqrL), n = (double)d_y.size();
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

    NumericVector deFeMod::rwBetaIncr(NumericMatrix const&sqrtXwt,
			              NumericVector const&  wtres) {
	if (sqrtXwt.size() != d_X.nrow())
	    Rf_error("%s: dimension mismatch %s(%d,%d), %s(%d,%d)",
		     "deFeMod::rwBetaIncr", "X", d_X.nrow(), d_X.ncol(),
		     "sqrtXwt", sqrtXwt.nrow(), sqrtXwt.ncol());
	int nc = sqrtXwt.ncol(), nr = sqrtXwt.nrow(), Xnc = d_X.ncol();
	double *V = d_V.x.begin(), *W = sqrtXwt.begin(), *X = d_X.x.begin();
	if (nc == 1) {
	    for (int j = 0; j < Xnc; j++)
		std::transform(X + j * nr, X + (j+1) * nr, W,
			       V + j * nr, std::multiplies<double>());
	} else {
	    Rf_error("code for nlmer not yet written");
	}
	d_V.dgemv('T', 1., wtres, 0., d_Vtr);
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

    NumericVector spFeMod::rwBetaIncr(NumericMatrix const& sqrtXwt,
			              NumericVector const&   wtres) {
//FIXME: Code here not written yet.
	// updateV(sqrtXwt);
	// updateVtr(wtres);
	// d_RX.update(d_V);
	chmDn cbeta(d_beta);
	NumericMatrix ans = d_RX.solve(CHOLMOD_A, &cbeta);
	return NumericVector(SEXP(ans));
    }

} // namespace mer

RCPP_FUNCTION_2(double,lmerDeUpdate,S4 xp,NumericVector nt) {
    mer::mer<mer::deFeMod,mer::lmerResp> lm(xp);
    lm.updateTheta(nt);
    lm.solveBetaU();
    lm.updateMu();
    return lm.Laplace();
}

RCPP_FUNCTION_2(double,lmerSpUpdate,S4 xp,NumericVector nt) {
    mer::mer<mer::spFeMod,mer::lmerResp> lm(xp);
    lm.updateTheta(nt);
    lm.solveBetaU();
    lm.updateMu();
    return lm.Laplace();
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
