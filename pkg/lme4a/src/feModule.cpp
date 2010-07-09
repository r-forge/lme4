#include "mer.h"

using namespace std;

namespace mer {
    feModule::feModule(Rcpp::S4 xp)
	: d_beta(Rcpp::clone(SEXP(xp.slot("beta")))),
	  d_Vtr(              d_beta.size()) {
    }

    void feModule::setBeta(Rcpp::NumericVector const& bbase,
			   Rcpp::NumericVector const&  incr,
			   double                      step) {
	R_len_t p = d_beta.size();
	if (bbase.size() != p)
	    throw runtime_error("feModule::setBeta size mismatch of beta and bbase");
	if (p == 0) return;
	if (step == 0.) {
	    std::copy(bbase.begin(), bbase.end(), d_beta.begin());
	} else {
//	    Rcpp::NumericVector res = cbase + incr * step;  // needs Rcpp_0.8.3
//	    copy(res.begin(), res.end(), d_coef.begin());
	    double *cb = bbase.begin(), *inc = incr.begin(), *cc = d_beta.begin();
	    for (R_len_t i = 0; i < p; i++) cc[i] = cb[i] + inc[i] * step;
	}
    }
}
