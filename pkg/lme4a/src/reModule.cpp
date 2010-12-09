#include "reModule.h"

using namespace Rcpp;
using namespace MatrixNs;
using namespace std;

namespace mer{
    reModule::reModule(Rcpp::S4 xp)
	: d_xp(             xp),
	  d_L(     S4(clone(SEXP(xp.slot("L"))))),
	  d_Lambda(S4(clone(SEXP(xp.slot("Lambda"))))),
	  d_Zt(          S4(xp.slot("Zt"))),
	  d_Lind(           xp.slot("Lind")),
	  d_lower(          xp.slot("lower")),
	  d_theta(          xp.slot("theta")),
	  d_u(              d_L.n),       
	  d_cu(             d_L.n) {
	d_Ut = (CHM_SP)NULL;
    }

    Rcpp::NumericVector reModule::linPred() const {
	NumericVector bb(d_u.size()), ans(d_Zt.nc());
	chmDn cans(ans), cbb(bb);
	d_Lambda.dmult('N',1.,0.,chmDn(d_u),cbb);
	d_Zt.dmult('T',1.,0.,cbb,cans);
	return ans;
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
    void reModule::setU(const Rcpp::NumericVector& ubase,
			const Rcpp::NumericVector& incr, double step) throw(std::runtime_error) {
	int q = d_u.size();
	if (ubase.size() != q)
	    throw std::runtime_error("size mismatch");
	    // Rf_error("%s: expected %s.size() = %d, got %d",
	    // 	     "reModule::setU", "ubase", q, ubase.size());
	NumericVector res = (step == 0.) ? ubase : ubase + incr * step;
	copy(res.begin(), res.end(), d_u.begin());
    }

    /** 
     * Solve for u (or the increment for u) only.
     * 
     */  
    double reModule::solveU() {
	NumericMatrix c1 = d_L.solve(CHOLMOD_L, d_L.solve(CHOLMOD_P, d_cu));
	double ans = inner_product(c1.begin(), c1.end(), c1.begin(), double());
	NumericMatrix mm = d_L.solve(CHOLMOD_Pt, d_L.solve(CHOLMOD_Lt, c1));
	copy(mm.begin(), mm.end(), d_u.begin());
	return ans;
    }

    /** 
     * Check and install new value of theta.  Update Lambda.
     * 
     * @param nt New value of theta
     */
    void reModule::setTheta(const Rcpp::NumericVector& nt) {
	R_len_t nth = d_lower.size();
	if (nt.size() != nth)
	    throw runtime_error("setTheta: size mismatch of nt and d_lower");
	if (any(nt < d_lower).is_true()) // check that nt is feasible
	    throw runtime_error("setTheta: theta not in feasible region");
	copy(nt.begin(), nt.end(), d_theta.begin());
				// update Lambda from theta and Lind
	double *Lamx = (double*)d_Lambda.x, *th = d_theta.begin();
	int *Li = d_Lind.begin(), nLind = d_Lind.size();
	for (R_len_t i = 0; i < nLind; i++) Lamx[i] = th[Li[i] - 1];
    }

    /** 
     * Solve for u given the updated cu
     * 
     * @param cu 
     */
    void reModule::updateU(const Rcpp::NumericVector& cu) {
	NumericMatrix nu = d_L.solve(CHOLMOD_Pt, d_L.solve(CHOLMOD_Lt, cu));
	setU(NumericVector(SEXP(nu)));
    }

}
