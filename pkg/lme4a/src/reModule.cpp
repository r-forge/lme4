#include "mer.h"

using namespace Rcpp;
using namespace MatrixNs;
using namespace std;

namespace mer{
    reModule::reModule(S4 xp)
	: d_xp(             xp),
	  d_L(     S4(clone(SEXP(xp.slot("L"))))),
	  d_Lambda(S4(clone(SEXP(xp.slot("Lambda"))))),
	  d_Zt(          S4(xp.slot("Zt"))),
	  d_Lind(           xp.slot("Lind")),
	  d_lower(          xp.slot("lower")),
	  d_u(              d_L.n),       
	  d_cu(             d_L.n) {
	d_Ut = (CHM_SP)NULL;
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
	if (d_Ut) M_cholmod_free_sparse(&d_Ut, &c);
	if (Wnc == 1) {
	    d_Ut = M_cholmod_copy_sparse(&d_Zt, &c);
	    chmDn csqrtX(Xwt);
	    M_cholmod_scale(&csqrtX, CHOLMOD_COL, d_Ut, &c);
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

	    d_Ut = M_cholmod_triplet_to_sparse(tr, nnz, &c);
	    M_cholmod_free_triplet(&tr, &c);
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
	NumericVector res = (step == 0.) ? ubase : ubase + incr * step;
	copy(res.begin(), res.end(), d_u.begin());
	d_sqrLenU = inner_product(d_u.begin(), d_u.end(), d_u.begin(), double());
    }

    /** 
     * Solve for u (or the increment for u) only.
     * 
     */  
    void reModule::solveU() {
	NumericMatrix ans = d_L.solve(CHOLMOD_A, d_cu);
	setU(NumericVector(SEXP(ans)));
    }

    // void reModule::updateDcmp(Rcpp::List& ll) const {
    // 	ll["ldL2"] = d_ldL2;
    // 	ll["ussq"] = d_sqrLenU;
    // 	Rcpp::S4 L = d_L.S4();
    // 	ll["L"]      = L;
    // }

    /** 
     * Check and install new value of theta.  Update Lambda.
     * 
     * @param nt New value of theta
     */
    void reModule::updateLambda(NumericVector const& nt) {
				// check that nt is feasible
	if (any(nt < d_lower).is_true())
	    throw runtime_error("updateLambda: theta not in feasible region");

	// R_len_t nth = d_lower.size(), Lis = d_Lind.size();
	// if (nt.size() != nth)
	//     throw runtime_error("size mismatch of nt and d_lower in updateLambda");
	// double *Lamx = (double*)d_Lambda.x, *ll = d_lower.begin(), *th = nt.begin();
	// for (R_len_t i = 0; i < nth; ++i)
	//     if (th[i] < ll[i])
	// 	throw runtime_error("updateLambda: theta not in feasible region");
				// update Lambda from theta and Lind
	double *Lamx = (double*)d_Lambda.x, *th = nt.begin();
	int *Li = d_Lind.begin();
	for (R_len_t i = 0; i < d_Lind.size(); i++) Lamx[i] = th[Li[i] - 1];
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
}