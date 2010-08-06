#include "feModule.h"

using namespace std;
using namespace MatrixNs;  // for chmDn, chmFr and chmSp

namespace mer {
    deFeMod::deFeMod(Rcpp::S4 xp, int n)
	: matMod::dPredModule(                  xp, n),
	  d_RZX(Rcpp::clone(Rcpp::S4(xp.slot("RZX")))),
	  d_UtV(d_RZX.nrow(), d_RZX.ncol()),
	  d_VtV(             d_coef.size()) {
    }

    /** 
     * Update UtV and VtV
     *
     * Note: In a dPredModule object there is no d_VtV member as the
     * d_fac member can be updated directly from d_V.  However, in the
     * deFeMod object it is necessary to form d_VtV separately.
     * Because a call to this member function follows a call to
     * reweight we do that here.  It may be better to create a
     * separate reweight method that first calls the reweight method
     * from the parent class.  Alternatively, we could include a d_VtV
     * member in the dPredModule class.
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
	d_fac.dpotrs(d_coef);
	d_RZX.dgemv('N', -1., d_coef, 1., ans);
	return ans;
    }
}
