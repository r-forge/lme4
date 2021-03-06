#include "reModule.h"

using namespace Rcpp;
using namespace MatrixNs;
using namespace std;

namespace mer{
    reModule::reModule(Rcpp::S4 xp)
	: d_L(     S4(clone(SEXP(xp.slot("L"))))),
	  d_Lambda(S4(clone(SEXP(xp.slot("Lambda"))))),
	  d_Zt(          S4(xp.slot("Zt"))),
	  d_Lind(           xp.slot("Lind")),
	  d_lower(          xp.slot("lower")),
	  d_theta(          xp.slot("theta")),
	  d_u0(             d_L.n),
	  d_incr(           d_L.n),
	  d_u(              d_L.n),
	  d_cu(             d_L.n) {
	d_Ut = (CHM_SP)NULL;
    }

    reModule::reModule(Rcpp::S4 Zt, Rcpp::S4 Lambda, Rcpp::S4 L,
	Rcpp::IntegerVector Lind, Rcpp::NumericVector lower)
	throw (MatrixNs::wrongS4)
	: d_L(L), d_Lambda(Lambda), d_Zt(Zt), d_Lind(Lind),
	  d_lower(lower), d_theta(lower.size()),
	  d_u0(d_Lambda.nr(), 0.), d_incr(d_Lambda.nr()),
	  d_u(d_Lambda.nr()),d_cu(d_Lambda.nr()) {
	d_Ut = (CHM_SP)NULL;
    }

    Rcpp::NumericVector reModule::b() const {
	NumericVector ans(d_u.size());
	chmDn cans(ans);
	d_Lambda.dmult('N',1.,0.,chmDn(d_u),cans);
	return ans;
    }

    Rcpp::NumericVector reModule::linPred() const {
	NumericVector bb = b(), ans(d_Zt.nc());
	chmDn cans(ans), cbb(bb);
	d_Zt.dmult('T',1.,0.,cbb,cans);
	return ans;
    }

    Rcpp::XPtr<cholmod_sparse> reModule::Utp() const {
	Rcpp::XPtr<cholmod_sparse> p(d_Ut, false);
	return p;
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
	copy(d_u0.begin(), d_u0.end(), d_cu.begin());
	M_cholmod_sdmult(LambdatUt, 0/*trans*/, &one, &mone, &cwtres, &ccu, &c);
	M_cholmod_free_sparse(&LambdatUt, &c);
	NumericMatrix
	    ans = d_L.solve(CHOLMOD_L, d_L.solve(CHOLMOD_P, d_cu));
	copy(ans.begin(), ans.end(), d_cu.begin());
	d_CcNumer = inner_product(d_cu.begin(), d_cu.end(), d_cu.begin(), double());
    }

    /** 
     * Solve for the increment only.
     * 
     */  
    double reModule::solveIncr() {
	NumericMatrix mm = d_L.solve(CHOLMOD_Pt, d_L.solve(CHOLMOD_Lt, d_cu));
	copy(mm.begin(), mm.end(), d_incr.begin());
	return d_CcNumer;
    }

    double reModule::sqrLenU() const {
	return inner_product(d_u.begin(), d_u.end(), d_u.begin(), double());
    }
	
    void reModule::installU0 () {
	copy(d_u.begin(), d_u.end(), d_u0.begin());
    }

    void reModule::setU0 (const Rcpp::NumericVector& uu)
	throw (std::runtime_error) {
	if (uu.size() != d_u0.size())
	    throw runtime_error("size mismatch");
	copy(uu.begin(), uu.end(), d_u0.begin());
    }

    void reModule::setIncr (const Rcpp::NumericVector& ii)
	throw (std::runtime_error) {
	if (ii.size() != d_incr.size())
	    throw runtime_error("size mismatch");
	copy(ii.begin(), ii.end(), d_incr.begin());
    }

    /** 
     * Check and install new value of theta.  Update Lambda.
     * 
     * @param nt New value of theta
     */
    void reModule::setTheta(const Rcpp::NumericVector& nt)
	throw (std::runtime_error) {
	R_len_t nth = d_lower.size(), nLind = d_Lind.size();
	if (nt.size() != nth)
	    throw runtime_error("setTheta: size mismatch of nt and d_lower");
#ifdef USE_RCPP_SUGAR
	if (any(nt < d_lower).is_true()) // check that nt is feasible
	    throw runtime_error("setTheta: theta not in feasible region");
#else
	double *lp = d_lower.begin(), *ntp = nt.begin();
	for (R_len_t i = 0; i < nth; i++)
	    if (ntp[i] < lp[i])
		throw runtime_error("setTheta: theta not in feasible region");
#endif
	copy(nt.begin(), nt.end(), d_theta.begin());
				// update Lambda from theta and Lind
	double *Lamx = (double*)d_Lambda.x, *th = d_theta.begin();
	int *Li = d_Lind.begin();
	for (R_len_t i = 0; i < nLind; i++) Lamx[i] = th[Li[i] - 1];
    }

    /** 
     * Solve for the increment given the updated cu
     * 
     * @param cu 
     */
    void reModule::updateIncr(const Rcpp::NumericVector& cu) {
	NumericMatrix nu = d_L.solve(CHOLMOD_Pt, d_L.solve(CHOLMOD_Lt, cu));
	copy(nu.begin(), nu.end(), d_incr.begin());
    }
    
    Rcpp::NumericVector reModule::linPred1(double fac) {
#if USE_RCPP_SUGAR	
	d_u = d_u0 + fac * d_incr;
#else
	fill(d_u.begin(), d_u.end(), fac);
	transform(d_u.begin(), d_u.end(), d_incr.begin(),
		  d_u.begin(), multiplies<double>());
	transform(d_u.begin(), d_u.end(), d_u0.begin(),
		  d_u.begin(), plus<double>());
#endif
	return linPred();
    }
}
