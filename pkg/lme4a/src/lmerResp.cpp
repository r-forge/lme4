#include "mer.h"

using namespace std;

namespace mer{
    lmerResp::lmerResp(Rcpp::S4 xp)
	: merResp(xp),
	  d_reml(*Rcpp::IntegerVector(xp.slot("REML")).begin()) {
	copy(d_offset.begin(), d_offset.end(), d_mu.begin());
	updateWrss();
	transform(d_weights.begin(), d_weights.end(), d_sqrtrwt.begin(), sqrtFun());
	copy(d_sqrtrwt.begin(), d_sqrtrwt.end(), d_sqrtXwt.begin());
    }

    double lmerResp::Laplace(double  ldL2,
			     double ldRX2,
			     double  sqrL) const {
	double lnum = 2.* PI * (d_wrss + sqrL), n = (double)d_y.size();
	if (d_reml == 0) return ldL2 + n * (1. + log(lnum / n));
	double nmp = n - d_reml;
	return ldL2 + ldRX2 + nmp * (1. + log(lnum / nmp));
    }

    double lmerResp::updateMu(Rcpp::NumericVector const &gamma) {
	copy(gamma.begin(), gamma.end(), d_mu.begin());
	return updateWrss();
    }
}
