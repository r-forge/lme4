#include "mer.h"

using namespace std;

namespace mer{
    nlmerResp::nlmerResp(Rcpp::S4 xp)
	: merResp(xp),
	  nlenv(SEXP(xp.slot("nlenv"))),
	  nlmod(SEXP(xp.slot("nlmod"))),
	  pnames(SEXP(xp.slot("pnames"))) {
    }

    double nlmerResp::Laplace(double  ldL2,
			      double ldRX2,
			      double  sqrL) const{
	double lnum = 2.* PI * (d_wrss + sqrL),
	    n = (double)d_y.size();
	return ldL2 + n * (1 + log(lnum / n));
    }

    double nlmerResp::updateMu(Rcpp::NumericVector const &gamma) {
	int n = d_y.size();
	double *gg = gamma.begin();

	for (int p = 0; p < pnames.size(); p++) {
	    string pn(pnames[p]);
	    Rcpp::NumericVector pp = nlenv.get(pn);
	    copy(gg + n * p, gg + n * (p + 1), pp.begin());
	}
	Rcpp::NumericVector rr = nlmod.eval(SEXP(nlenv));
	if (rr.size() != n)
	    Rf_error("Length mu = %d, expected %d", rr.size(), n);
	copy(rr.begin(), rr.end(), d_mu.begin());
	Rcpp::NumericMatrix rrg = rr.attr("gradient");
	copy(rrg.begin(), rrg.end(), d_sqrtXwt.begin());
	return updateWrss();
    }

    void nlmerResp::updateDcmp(Rcpp::List& ll) const {
	merResp::updateDcmp(ll);
    }
}
