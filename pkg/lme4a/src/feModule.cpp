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
	if (p == 0) return;
	if (bbase.size() != p)
	    throw runtime_error("feModule::setBeta size mismatch of beta and bbase");

	Rcpp::NumericVector
	    res = (step == 0.) ? bbase : bbase + incr * step;
	copy(res.begin(), res.end(), d_beta.begin());
    }
}
