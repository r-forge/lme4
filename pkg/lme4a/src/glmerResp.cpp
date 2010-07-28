#include "respModule.h"

using namespace std;

namespace mer{

    glmerResp::glmerResp(Rcpp::S4 xp)
	: merResp(    xp),
	  family(SEXP(xp.slot("family"))),
	  d_eta(      xp.slot("eta")),
	  d_n(        xp.slot("n")) {
	updateWts();
    }
    
    Rcpp::NumericVector glmerResp::devResid() const {
	return family.devResid(d_mu, d_weights, d_y);
    }
	
    double glmerResp::Laplace(double  ldL2,
			      double ldRX2,
			      double  sqrL) const{
	Rcpp::NumericVector dr = devResid();
	return ldL2 + sqrL + accumulate(dr.begin(), dr.end(), double());
    }
	
    double glmerResp::updateMu(Rcpp::NumericVector const &gamma) {
	transform(gamma.begin(), gamma.end(), d_offset.begin(),
		  d_eta.begin(), plus<double>());
	Rcpp::NumericVector li = family.linkInv(d_eta);
	copy(li.begin(), li.end(), d_mu.begin());
	return updateWrss();
    }

    struct sqrtquotFun : std::binary_function<double,double,double> {
	inline double operator() (double x, double y) {
	    return sqrt(x / y);
	}
    };

    double glmerResp::updateWts() {
	Rcpp::NumericVector mueta = family.  muEta(d_eta),
	                       vv = family.variance(d_mu);
	transform(d_weights.begin(), d_weights.end(), vv.begin(),
		  d_sqrtrwt.begin(), sqrtquotFun());
	transform(d_sqrtrwt.begin(), d_sqrtrwt.end(), mueta.begin(),
		  d_sqrtXwt.begin(), multiplies<double>());
	return updateWrss();
    }
}
