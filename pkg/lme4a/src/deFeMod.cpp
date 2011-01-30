#include "feModule.h"
#include <R_ext/BLAS.h>
#include "utilities.h"

using namespace std;
using namespace MatrixNs;  // for chmDn, chmFr and chmSp

namespace mer {
    deFeMod::deFeMod(Rcpp::S4 xp, int n)
	: matMod::dPredModule(                  xp, n),
	  d_RZX(Rcpp::clone(Rcpp::S4(xp.slot("RZX")))),
	  d_UtV(d_RZX.nrow(), d_RZX.ncol()),
	  d_VtV(             d_coef.size()),
	  d_coef0(       d_coef.size(), 0.),
	  d_incr(            d_coef.size()) {
    }

    deFeMod::deFeMod(Rcpp::S4 mm, int n, int p, int q)
	: matMod::dPredModule(mm,n,p), d_RZX(q, p),
	  d_UtV(q, p), d_VtV(p, 'U'),
	  d_coef0(p), d_incr(p) {
    }

    /** 
     * Update UtV and VtV
     *
     * Note: In a dPredModule object there is no d_VtV member as the
     * d_fac member can be updated directly from d_V.  However, in the
     * deFeMod object it is necessary to form d_VtV separately.
     *
     * @param Ut from the reModule
     */
    void deFeMod::updateUtV(cholmod_sparse   const* Ut) {
	if (d_coef.size() == 0) return;
     	double one = 1., zero = 0.;
     	chmDn cUtV(d_UtV), cV(d_V);
     	M_cholmod_sdmult(Ut, 0/*trans*/, &one, &zero, &cV, &cUtV, &c);
	d_VtV.dsyrk(d_V, 1., 0.);
    }

    void deFeMod::updateUtVp(Rcpp::XPtr<cholmod_sparse> p) {
	updateUtV(p);
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
	if (d_coef.size() == 0) return;
	chmDn cRZX(d_RZX);
	Lambda.dmult('T', 1., 0., chmDn(d_UtV), cRZX);
	Rcpp::NumericMatrix
	    ans = L.solve(CHOLMOD_L, L.solve(CHOLMOD_P, &cRZX));
	d_RZX.setX(ans);
	d_fac.update('T', -1., d_RZX, 1., d_VtV);
	d_ldRX2 = d_fac.logDet2();
    }
    
    void deFeMod::updateRzxpRxpp(Rcpp::XPtr<MatrixNs::chmSp> Lambdap,
				 Rcpp::XPtr<MatrixNs::chmFr> Lp) {
	updateRzxRx(*Lambdap, *Lp);
    }

    // Need to copy this here, at least for the present.  A module
    // will not be able to find the definition in the dPredModule
    // class

    Rcpp::NumericVector deFeMod::linPred() const {
	Rcpp::NumericVector ans(d_X.nrow());
	d_X.dgemv('N', 1., d_coef, 0., ans);
	return ans;
    }

    Rcpp::NumericVector deFeMod::linPred1(double fac) const {
#ifdef USE_RCPP_SUGAR
	Rcpp::NumericVector cc = d_coef0 + fac * d_incr;
	copy(cc.begin(), cc.end(), d_coef.begin());
#else
	fill(d_coef.begin(), d_coef.end(), fac);
	transform(d_coef.begin(), d_coef.end(), d_incr.begin(),
		  d_coef.begin(), multiplies<double>());
	transform(d_coef.begin(), d_coef.end(), d_coef0.begin(),
		  d_coef.begin(), plus<double>());
#endif
	return linPred();
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
	Rcpp::NumericVector ans = clone(cu);
	if (d_coef.size() == 0) return ans;

	copy(d_Vtr.begin(), d_Vtr.end(), d_coef.begin());
	d_RZX.dgemv('T', -1., ans, 1., d_coef);
	d_fac.dtrtrs('T', d_coef.begin());
	d_CcNumer = sum(d_coef * d_coef);
	d_fac.dtrtrs('N', d_coef.begin());
	d_RZX.dgemv('N', -1., d_coef, 1., ans);
	return ans;
    }

    void deFeMod::solveIncr() {
	d_fac.update(d_VtV);
	copy(d_Vtr.begin(), d_Vtr.end(), d_incr.begin());
	d_fac.dtrtrs('T', d_incr.begin());
	d_fac.dtrtrs('N', d_incr.begin());
    }
    
    Rcpp::NumericVector deFeMod::updateIncr(Rcpp::NumericVector const& cu) {
	Rcpp::NumericVector ans = clone(cu);
	if (d_coef.size() == 0) return ans;
#ifdef LME4A_DEBUG
showdbl(ans, "cu on entry");
showdbl(d_Vtr, "Vtr on entry");
#endif
	copy(d_Vtr.begin(), d_Vtr.end(), d_incr.begin());
	d_RZX.dgemv('T', -1., ans, 1., d_incr);
#ifdef LME4A_DEBUG
showdbl(d_incr, "Vtr - crossprod(RZX, cu)");
#endif
	d_fac.dtrtrs('T', d_incr.begin());
#ifdef LME4A_DEBUG
showdbl(d_incr, "RX^{-T}(Vtr - crossprod(RZX, cu))");
#endif
#ifdef USE_RCPP_SUGAR
        d_CcNumer = sum(d_incr * d_incr);
#else
        Rcpp::NumericVector tmp(d_incr.size());
	transform(d_incr.begin(), d_incr.end(), d_incr.begin(),
		  tmp.begin(), multiplies<double>());
	d_CcNumer = accumulate(tmp.begin(), tmp.end(), double());
#endif
        d_fac.dtrtrs('N', d_incr.begin());
#ifdef LME4A_DEBUG
showdbl(d_incr, "RX^{-1}RX^{-T}(Vtr - crossprod(RZX, cu))");
#endif
	d_RZX.dgemv('N', -1., d_incr, 1., ans);
#ifdef LME4A_DEBUG
showdbl(ans, "updated cu");
#endif
	return ans;
    }
    
    void deFeMod::installCoef0() {
	copy(d_coef.begin(), d_coef.end(), d_coef0.begin());
    }

    void deFeMod::setCoef0 (const Rcpp::NumericVector& cc)
	throw (runtime_error) {
	if (cc.size() != d_coef0.size())
	    throw runtime_error("setCoef0: size mismatch");
	copy(cc.begin(), cc.end(), d_coef0.begin());
    }

    void deFeMod::setIncr (const Rcpp::NumericVector& ii)
	throw (runtime_error) {
	if (ii.size() != d_incr.size())
	    throw runtime_error("setIncr: size mismatch");
	copy(ii.begin(), ii.end(), d_incr.begin());
    }

    /** 
     * Update V, VtV and Vtr
     * 
     * @param Xwt square root of the weights for the model matrices
     * @param wtres weighted residuals
     */
    void deFeMod::reweight(Rcpp::NumericMatrix   const&   Xwt,
			   Rcpp::NumericVector   const& wtres)
	throw(runtime_error) {
	if (d_coef.size() == 0) return;
	chmDn cXwt(Xwt);
	if ((Xwt.rows() * Xwt.cols()) != d_X.nrow())
	    throw runtime_error("dimension mismatch");
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
	    Rcpp::NumericVector tmp(Xnr), mm(Wnc); 
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
	d_VtV.dsyrk(d_V, 1., 0.);
    }

}
