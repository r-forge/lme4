#include "respModule.h"
#include <cmath>

using namespace Rcpp;

namespace mer {

    merResp::merResp(Rcpp::S4 xp)  throw (std::runtime_error)
	: d_y(                    xp.slot("y")),
	  d_weights(        xp.slot("weights")),
	  d_offset(          xp.slot("offset")),
	  d_mu(                 d_y.size(), 0.),
	  d_sqrtrwt(                d_y.size()),
	  d_wtres(                  d_y.size()),
	  d_sqrtXwt(d_y.size(), d_offset.size()/d_y.size()) {
	int n = d_y.size(), os = d_offset.size();
	if (d_mu.size() != n || d_weights.size() != n || d_sqrtrwt.size() != n)
	    throw std::runtime_error("y, mu, sqrtrwt and weights slots must have equal lengths");
	if (os < 1 || os % n)
	    throw std::runtime_error("length(offset) must be a positive multiple of length(y)");
	init();
    }

    merResp::merResp(Rcpp::NumericVector y)
	throw (std::runtime_error)
	: d_y(y), d_weights(y.size(), 1.0), d_offset(y.size()),
	  d_mu(y.size(), 0.), d_sqrtrwt(y.size(), 1.),
	  d_wtres(y.size()),
	  d_sqrtXwt(y.size(), 1) {
	init();
    }

    merResp::merResp(Rcpp::NumericVector y, Rcpp::NumericVector weights)
	throw (std::runtime_error)
	: d_y(y), d_weights(weights), d_offset(y.size()),
	  d_mu(y.size()), d_sqrtrwt(y.size()), d_wtres(y.size()),
	  d_sqrtXwt(y.size(), 1) {
	if (weights.size() != y.size())
	    throw std::runtime_error(
		"lengths of y and wts must agree");
	init();
    }

    merResp::merResp(Rcpp::NumericVector y, Rcpp::NumericVector weights,
	Rcpp::NumericVector offset) throw (std::runtime_error)
	: d_y(y), d_weights(weights), d_offset(offset),
	  d_mu(y.size()), d_sqrtrwt(y.size()), d_wtres(y.size()),
	  d_sqrtXwt(y.size(), 1) {
	int nn = y.size();
	if (weights.size() != nn || offset.size() != nn)
	    throw std::runtime_error(
		"lengths of y, weights and offset must agree");
	init();
    }

    void merResp::init() {
#ifdef USE_RCPP_SUGAR
	d_sqrtrwt = sqrt(d_weights);
#else
	std::transform(d_weights.begin(), d_weights.end(), d_sqrtrwt.begin(), &::sqrt);
#endif
	std::copy(d_sqrtrwt.begin(), d_sqrtrwt.end(), d_sqrtXwt.begin());
	updateWrss();
    }

    /** 
     * Update the wtres vector and return its sum of squares
     *   wtres <- sqrtrwt * (y - mu)
     *   return(wrss <- sum(wtres^2))
     *
     * @return Updated weighted residual sum of squares
     */
    double merResp::updateWrss() {
#ifdef USE_RCPP_SUGAR
	d_wtres = (d_y - d_mu) * d_sqrtrwt;
	d_wrss = sum(d_wtres * d_wtres);
#else
	NumericVector tmp(d_y.size());
	std::transform(d_y.begin(), d_y.end(), d_mu.begin(), tmp.begin(),
		       std::minus<double>());
	std::transform(tmp.begin(), tmp.end(), d_sqrtrwt.begin(),
		       d_wtres.begin(), std::multiplies<double>());
	std::transform(d_wtres.begin(), d_wtres.end(), d_wtres.begin(),
		       tmp.begin(), std::multiplies<double>());
	d_wrss = std::accumulate(tmp.begin(), tmp.end(), double());
#endif
	return d_wrss;
    }

    Rcpp::NumericVector merResp::devResid() const {
	return NumericVector(0);
    }

    lmerResp::lmerResp(Rcpp::S4 xp) throw (std::runtime_error)
	: merResp(xp),
	  d_reml(*Rcpp::IntegerVector(xp.slot("REML")).begin()) {
	std::copy(d_offset.begin(), d_offset.end(), d_mu.begin());
	updateWrss();
    }

    lmerResp::lmerResp(int rr, Rcpp::NumericVector y)
	throw (std::runtime_error)
	: merResp(y), d_reml(rr) {
    }

    lmerResp::lmerResp(int rr, Rcpp::NumericVector y, Rcpp::NumericVector weights)
	throw (std::runtime_error)
	: merResp(y, weights), d_reml(rr) {
    }

    lmerResp::lmerResp(int rr, Rcpp::NumericVector y, Rcpp::NumericVector weights,
	Rcpp::NumericVector offset) throw (std::runtime_error)
	: merResp(y, weights, offset), d_reml(rr) {
    }

    double lmerResp::Laplace(double  ldL2,
			     double ldRX2,
			     double  sqrL) const {
	double lnum = 2.* PI * (d_wrss + sqrL), n = (double)d_y.size();
	if (d_reml == 0) return ldL2 + n * (1. + log(lnum / n));
	double nmp = n - d_reml;
	return ldL2 + ldRX2 + nmp * (1. + log(lnum / nmp));
    }

    void lmerResp::setWeights(const Rcpp::NumericVector& ww)
	throw (std::runtime_error) {
	if (ww.size() != d_weights.size())
	    throw std::runtime_error("setWeights: Size mismatch");
	std::copy(ww.begin(), ww.end(), d_weights.begin());
    }
					       
    void lmerResp::setOffset(const Rcpp::NumericVector& oo)
	throw (std::runtime_error) {
	if (oo.size() != d_offset.size())
	    throw std::runtime_error("setOffset: Size mismatch");
	std::copy(oo.begin(), oo.end(), d_offset.begin());
    }
					       
    void lmerResp::setReml(int rr) throw (std::runtime_error) {
	if (rr < 0)
	    throw std::runtime_error("setReml: negative rr");
	d_reml = rr;
    }
					       
    double lmerResp::updateMu(const Rcpp::NumericVector& gamma) {
#ifdef USE_RCPP_SUGAR
	d_mu = d_offset + gamma;
#else
	std::transform(gamma.begin(), gamma.end(), d_offset.begin(),
		       d_mu.begin(), std::plus<double>());
#endif
	return updateWrss();
    }
    
    glmerResp::glmerResp(Rcpp::S4 xp) throw (std::runtime_error)
	: merResp(xp),
	  d_fam(SEXP(xp.slot("family"))),
	  d_n(       xp.slot("n")),
	  d_eta(     xp.slot("eta")) {
	updateWts();
    }

    glmerResp::glmerResp(Rcpp::List fam, Rcpp::NumericVector y)
	throw (std::runtime_error)
	: merResp(y), d_fam(fam), d_eta(y.size()) {
    }

    glmerResp::glmerResp(Rcpp::List fam, Rcpp::NumericVector y,
			 Rcpp::NumericVector weights)
	throw (std::runtime_error)
	: merResp(y, weights), d_fam(fam), d_eta(y.size()) {
    }

    glmerResp::glmerResp(Rcpp::List fam, Rcpp::NumericVector y,
			 Rcpp::NumericVector weights,
			 Rcpp::NumericVector offset)
	throw (std::runtime_error) 
	: merResp(y, weights, offset), d_fam(fam),
	  d_eta(y.size()) {
    }

    glmerResp::glmerResp(Rcpp::List fam, Rcpp::NumericVector y,
			 Rcpp::NumericVector weights,
			 Rcpp::NumericVector offset,
			 Rcpp::NumericVector n)
	throw (std::runtime_error)
 	: merResp(y, weights, offset), d_fam(fam), d_n(n),
	  d_eta(y.size()) {
	if (n.size() != y.size())
	    throw std::runtime_error("lengths of y and n must agree");
    }

    glmerResp::glmerResp(Rcpp::List fam, Rcpp::NumericVector y,
		     Rcpp::NumericVector weights,
		     Rcpp::NumericVector offset,
		     Rcpp::NumericVector n, Rcpp::NumericVector eta)
	throw (std::runtime_error) 
	: merResp(y, weights, offset), d_fam(fam), d_n(n),
	  d_eta(eta) {
	int nn = y.size();
	if (n.size() != nn || eta.size() != nn )
	    throw std::runtime_error(
		"lengths of y, n and eta must agree");
    }

    double glmerResp::residDeviance() const {
#ifdef USE_RCPP_SUGAR
	return sum(devResid());
#else
	NumericVector dd = devResid();
	return std::accumulate(dd.begin(), dd.end(), double());
#endif
    }

    double glmerResp::updateWts() {
#ifdef USE_RCPP_SUGAR
	d_sqrtrwt = sqrt(d_weights/variance());
	NumericVector tmp = muEta() * d_sqrtrwt;
#else
	NumericVector vv = variance();
	std::transform(d_weights.begin(), d_weights.end(), vv.begin(),
		       d_sqrtrwt.begin(), std::divides<double>());
	std::transform(d_sqrtrwt.begin(), d_sqrtrwt.end(),
		       d_sqrtrwt.begin(), &::sqrt);
	NumericVector tmp = muEta();
	std::transform(tmp.begin(), tmp.end(), d_sqrtrwt.begin(),
		       tmp.begin(), std::multiplies<double>());
#endif
	std::copy(tmp.begin(), tmp.end(), d_sqrtXwt.begin());
	
	return updateWrss();
    }

    Rcpp::NumericVector glmerResp::wrkResids() const {
#ifdef USE_RCPP_SUGAR
	return (d_y - d_mu)/muEta();
#else
	NumericVector rr(d_y.size()), me = muEta();
	std::transform(d_y.begin(), d_y.end(), d_mu.begin(),
		       rr.begin(), std::minus<double>());
	std::transform(rr.begin(), rr.end(), me.begin(),
		       rr.begin(), std::divides<double>());
	return rr;
#endif
    }

    Rcpp::NumericVector glmerResp::wrkResp() const {
#ifdef USE_RCPP_SUGAR
	return (d_eta - d_offset) + wrkResids();
#else
	NumericVector rr(d_eta.size()), ww = wrkResids();
	std::transform(d_eta.begin(), d_eta.end(), d_offset.begin(),
		       rr.begin(), std::minus<double>());
	std::transform(rr.begin(), rr.end(), ww.begin(), rr.begin(),
		       std::plus<double>());
	return rr;
#endif	
    }

    Rcpp::NumericVector glmerResp::sqrtWrkWt() const {
	NumericVector me = muEta();
#ifdef USE_RCPP_SUGAR
	return sqrt(d_weights * me * me / variance());
#else
	NumericVector vv = variance();
	std::transform(me.begin(), me.end(), me.begin(), me.begin(),
		       std::multiplies<double>());
	std::transform(me.begin(), me.end(), d_weights.begin(),
		       me.begin(), std::multiplies<double>());
	std::transform(me.begin(), me.end(), vv.begin(), me.begin(),
		       std::divides<double>());
	std::transform(me.begin(), me.end(), me.begin(), &::sqrt);
	return me;
#endif	
    }

    double glmerResp::Laplace(double  ldL2,
			      double ldRX2,
			      double  sqrL) const{
	return ldL2 + sqrL + residDeviance();
    }

    double glmerResp::updateMu(const Rcpp::NumericVector& gamma) {
#ifdef USE_RCPP_SUGAR
	d_eta = d_offset + gamma;
#else
	std::transform(d_offset.begin(), d_offset.end(), gamma.begin(),
		       d_eta.begin(), std::plus<double>());
#endif
	NumericVector mmu = d_fam.linkInv(d_eta);
	std::copy(mmu.begin(), mmu.end(), d_mu.begin());
	return updateWrss();
    }

    void glmerResp::setPwrss(double val) {
	d_pwrss = val;
    }

    nlmerResp::nlmerResp(Rcpp::S4 xp)
	: merResp(xp),
	  d_nlenv(SEXP(xp.slot("nlenv"))),
	  d_nlmod(SEXP(xp.slot("nlmod"))),
	  d_pnames(SEXP(xp.slot("pnames"))) {
    }

    double nlmerResp::Laplace(double  ldL2,
			      double ldRX2,
			      double  sqrL) const{
	double lnum = 2.* PI * (d_wrss + sqrL),
	    n = (double)d_y.size();
	return ldL2 + n * (1 + log(lnum / n));
    }

    double nlmerResp::updateMu(Rcpp::NumericVector const &gamma) throw(std::runtime_error) {
	int n = d_y.size();
#ifdef USE_RCPP_SUGAR
	Rcpp::NumericVector gam = gamma + d_offset;
#else
	NumericVector gam(d_offset.size());
	std::transform(gamma.begin(), gamma.end(), d_offset.begin(),
		       gam.begin(), std::plus<double>());
#endif
	double *gg = gam.begin();

	for (int p = 0; p < d_pnames.size(); p++) {
	    std::string pn(d_pnames[p]);
	    Rcpp::NumericVector pp = d_nlenv.get(pn);
	    std::copy(gg + n * p, gg + n * (p + 1), pp.begin());
	}
	NumericVector rr = d_nlmod.eval(SEXP(d_nlenv));
	if (rr.size() != n)
	    throw std::runtime_error("dimension mismatch");
	std::copy(rr.begin(), rr.end(), d_mu.begin());
	NumericMatrix rrg = rr.attr("gradient");
	std::copy(rrg.begin(), rrg.end(), d_sqrtXwt.begin());
	return updateWrss();
    }
}
