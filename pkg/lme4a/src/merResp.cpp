#include "mer.h"

using namespace std;

namespace mer{
    merResp::merResp(Rcpp::S4 xp)
	: d_xp(                      xp),
	  d_offset(                  xp.slot("offset")),
	  d_weights(                 xp.slot("weights")),
	  d_y(                       xp.slot("y")),
	  d_sqrtrwt(                d_y.size()),
	  d_wtres(                  d_y.size()),
	  d_mu(                     d_y.size()),
	  d_sqrtXwt(d_y.size(), d_offset.size()/d_y.size()) {
	int n = d_y.size(), os = d_offset.size();
	if (d_mu.size() != n || d_weights.size() != n || d_sqrtrwt.size() != n)
	    Rf_error("y, mu, sqrtrwt and weights slots must have equal lengths");
	if (os < 1 || os % n)
	    Rf_error("length(offset) must be a positive multiple of length(y)");
    }

    /** 
     * Update the wtres vector and return its sum of squares
     *   wtres <- sqrtrwts * (y - mu)
     *   return(wrss <- sum(wtres^2))
     *
     * @return Updated weighted residual sum of squares
     */
    double merResp::updateWrss() {
				// wtres <- y - mu
	transform(d_y.begin(), d_y.end(), d_mu.begin(),
		  d_wtres.begin(), minus<double>());
				// wtres <- wtres * sqrtrwt
	transform(d_wtres.begin(), d_wtres.end(), d_sqrtrwt.begin(),
		  d_wtres.begin(), multiplies<double>());
	d_wrss = inner_product(d_wtres.begin(), d_wtres.end(),
			       d_wtres.begin(), double());
	return d_wrss;
    }

    /** 
     * Update the devcomp$cmp vector in place by calculating and
     * installing the wrss (weighted residual sum of squares), pwrss
     * (penalized weighted ...) and sigmaML elements.
     * 
     * @param cmp the devcomp list to update
     */
    void merResp::updateDcmp(Rcpp::List& ll) const {
	Rcpp::List devcomp = ll["devcomp"];
	Rcpp::NumericVector cmp = devcomp["cmp"];

	double wrss = inner_product(d_wtres.begin(), d_wtres.end(),
				    d_wtres.begin(), double());
	double n = (double)d_y.size(), ussq = cmp["ussq"];
	double pwrss = wrss + ussq;
	cmp["wrss"] = wrss;
	cmp["pwrss"] = pwrss;
	cmp["sigmaML"] = sqrt(pwrss/n);
	ll["mu"] = d_mu;
    }

}
