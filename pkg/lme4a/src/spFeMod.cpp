#include "mer.h"

using namespace std;
using namespace MatrixNs;

namespace mer {
    spFeMod::spFeMod(Rcpp::S4 xp, int n)
	: feModule(      xp),
	  d_RZX(Rcpp::S4(xp.slot("RZX"))),
	  d_X(  Rcpp::S4(xp.slot("X"))),
	  d_RX( Rcpp::S4(xp.slot("RX"))) {
				// Not sure if this is necessary
	d_UtV = d_V = d_VtV = (CHM_SP)NULL;
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
    Rcpp::NumericVector spFeMod::updateBeta(Rcpp::NumericVector const &cu) {
	Rcpp::NumericVector ans = clone(cu);
	copy(cu.begin(), cu.end(), ans.begin());
	copy(d_Vtr.begin(), d_Vtr.end(), d_beta.begin());
	chmDn cbeta(d_beta), cans(ans);
	d_RZX.dmult('T', -1., 1., chmDn(cu), cbeta);
	Rcpp::NumericMatrix t1 = d_RX.solve(CHOLMOD_A, &cbeta);
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
     * @param sqrtXwt square root of the weights for the model matrices
     * @param wtres weighted residuals
     */
    void spFeMod::reweight(cholmod_sparse      const*      Ut,
			   Rcpp::NumericMatrix const& sqrtXwt,
			   Rcpp::NumericVector const&   wtres) {
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
	d_UtV = M_cholmod_ssmult(Ut, d_V, 0/*styp*/,1/*vals*/,1/*srtd*/, &c);

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
	Rcpp::NumericMatrix ans = d_RX.solve(CHOLMOD_A, d_Vtr);
	copy(ans.begin(), ans.end(), d_beta.begin());
    }

//    void spFeMod::updateDcmp(Rcpp::NumericVector& cmp) const {  // needs Matrix_0.999375-42 or later
    // void spFeMod::updateDcmp(Rcpp::List& ll) {
    // 	ll["ldRX2"] = d_ldRX2;
    // 	Rcpp::S4 RX = d_RX.S4();
    // 	ll["RX"] = RX;
    // }
}
