#include "feModule.h"

using namespace std;
using namespace MatrixNs;

namespace mer {
    spFeMod::spFeMod(Rcpp::S4 xp, int n)
	: matMod::sPredModule(    xp, n),
	  d_RZX(Rcpp::clone(Rcpp::S4(xp.slot("RZX")))) {
				// Not sure if this is necessary
	d_UtV = d_VtV = (CHM_SP)NULL;
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
    Rcpp::NumericVector
    spFeMod::updateBeta(Rcpp::NumericVector const &cu) {
	Rcpp::NumericVector ans(cu.size());
	copy(cu.begin(), cu.end(), ans.begin());
	copy(d_Vtr.begin(), d_Vtr.end(), d_coef.begin());
	chmDn ccoef(d_coef), cans(ans);
	d_RZX.dmult('T', -1., 1., chmDn(cu), ccoef);
	Rcpp::NumericMatrix t1 = d_fac.solve(CHOLMOD_A, &ccoef);
	copy(t1.begin(),  t1.end(), d_coef.begin());
	d_RZX.dmult('N', -1., 1., ccoef, cans);
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
	d_fac.update(*t2);

	M_cholmod_free_sparse(&t2, &c);

	d_ldRX2 = d_fac.logDet2();
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
    void spFeMod::updateUtV(cholmod_sparse   const* Ut) {
 	if (d_UtV) M_cholmod_free_sparse(&d_UtV, &c);
 	d_UtV = M_cholmod_ssmult(Ut, d_V, 0/*styp*/,
				 1/*vals*/,1/*srtd*/, &c);
 	if (d_VtV) M_cholmod_free_sparse(&d_VtV, &c);
 	CHM_SP t1 = M_cholmod_transpose(d_V, 1/*vals*/, &c);
 	d_VtV = M_cholmod_ssmult(t1, d_V, 1/*styp*/,1/*vals*/,
				 1/*srtd*/, &c);
 	M_cholmod_free_sparse(&t1, &c);
    }

    // /** 
    //  * Reweight the feModule.
    //  *
    //  * Update V, UtV and Vtr
    //  * 
    //  * @param Ut from the reModule
    //  * @param Xwt square root of the weights for the model matrices
    //  * @param wtres weighted residuals
    //  */
//     void spFeMod::reweight(cholmod_sparse      const*      Ut,
// 			   Rcpp::NumericMatrix const&     Xwt,
// 			   Rcpp::NumericVector const&   wtres) {
// 	double one = 1., zero = 0.;
// 	if (d_beta.size() == 0) return;
// 	int Wnc = Xwt.ncol(), Wnr = Xwt.nrow(),
// 	    Xnc = d_X.ncol, Xnr = d_X.nrow;
// 	if ((Xwt.rows() * Xwt.cols()) != (int)d_X.nrow)
// 	    Rf_error("%s: dimension mismatch %s(%d,%d), %s(%d,%d)",
// 		     "spFeMod::reweight", "X", Xnr, Xnc,
// 		     "Xwt", Wnr, Wnc);
// 	if (Wnc == 1) {
// 	    if (d_V) M_cholmod_free_sparse(&d_V, &c);
// 	    d_V = M_cholmod_copy_sparse(&d_X, &c);
// 	    chmDn csqrtX(Xwt);
// 	    M_cholmod_scale(&csqrtX, CHOLMOD_ROW, d_V, &c);
// 	} else throw runtime_error("spFeMod::reweight: multiple columns in Xwt");
// // FIXME rewrite this using the triplet representation
	
// 	if (d_UtV) M_cholmod_free_sparse(&d_UtV, &c);
// 	d_UtV = M_cholmod_ssmult(Ut, d_V, 0/*styp*/,1/*vals*/,1/*srtd*/, &c);

// 	if (d_VtV) M_cholmod_free_sparse(&d_VtV, &c);
// 	CHM_SP t1 = M_cholmod_transpose(d_V, 1/*vals*/, &c);
// 	d_VtV = M_cholmod_ssmult(t1, d_V, 1/*styp*/,1/*vals*/,1/*srtd*/, &c);
// 	M_cholmod_free_sparse(&t1, &c);

// 	chmDn cVtr(d_Vtr);
// 	const chmDn cwtres(wtres);
// 	M_cholmod_sdmult(d_V, 'T', &one, &zero, &cwtres, &cVtr, &c);
//     }

    // /** 
    //  * Solve (V'V)beta = Vtr for beta.
    //  * 
    //  */
//    double spFeMod::solveBeta() {
// FIXME: This may be a bad idea.  d_RX may have a more complicated
// structure than the factor of d_VtV.
    // 	d_RX.update(*d_VtV);
    // 	Rcpp::NumericMatrix c1 = d_RX.solve(CHOLMOD_L, d_RX.solve(CHOLMOD_P, d_Vtr));
    // 	double ans = inner_product(c1.begin(), c1.end(), c1.begin(), double());
    // 	Rcpp::NumericMatrix mm = d_RX.solve(CHOLMOD_Pt, d_RX.solve(CHOLMOD_Lt, c1));
    // 	copy(mm.begin(), mm.end(), d_beta.begin());
    // 	return ans;
    // }
}
