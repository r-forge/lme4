#include "respModule.h"

using namespace Rcpp;

namespace mer {

    merResp::merResp(Rcpp::S4 xp)  throw (std::runtime_error)
	: d_y(                       xp.slot("y")),
	  d_weights(                 xp.slot("weights")),
	  d_offset(                  xp.slot("offset")),
	  d_mu(                     d_y.size()),
	  d_sqrtrwt(                d_y.size()),
	  d_wtres(                  d_y.size()),
	  d_sqrtXwt(d_y.size(), d_offset.size()/d_y.size()) {
	int n = d_y.size(), os = d_offset.size();
	if (d_mu.size() != n || d_weights.size() != n || d_sqrtrwt.size() != n)
	    throw std::runtime_error("y, mu, sqrtrwt and weights slots must have equal lengths");
	if (os < 1 || os % n)
	    throw std::runtime_error("length(offset) must be a positive multiple of length(y)");
	d_sqrtrwt = sqrt(d_weights);
	std::copy(d_sqrtrwt.begin(), d_sqrtrwt.end(), d_sqrtXwt.begin());
    }

    merResp::merResp(Rcpp::NumericVector y)
	throw (std::runtime_error)
	: d_y(y), d_weights(y.size(), 1.0), d_offset(y.size()),
	  d_mu(y.size()), d_sqrtrwt(y.size()),
	  d_wtres(y.size()),
	  d_sqrtXwt(y.size(), 1.0) {
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
	d_sqrtrwt = sqrt(d_weights);
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
	d_wtres = (d_y - d_mu) * d_sqrtrwt;
	d_wrss = sum(d_wtres * d_wtres);
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

    double lmerResp::updateMu(const Rcpp::NumericVector& gamma) {
	d_mu = d_offset + gamma;
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

    Rcpp::NumericVector glmerResp::devResid() const {
	return d_fam.devResid(d_mu, d_weights, d_y);
    }

    Rcpp::NumericVector glmerResp::muEta() const {
	return d_fam.muEta(d_eta);
    }

    Rcpp::NumericVector glmerResp::variance() const {
	return d_fam.variance(d_mu);
    }

    double glmerResp::residDeviance() const {
	return sum(devResid());
    }

    double glmerResp::updateWts() {
	d_sqrtrwt = sqrt(d_weights/d_fam.variance(d_mu));
	NumericVector tmp = muEta() * d_sqrtrwt;
	std::copy(tmp.begin(), tmp.end(), d_sqrtXwt.begin());
	return updateWrss();
    }

    Rcpp::NumericVector glmerResp::wrkResids() const {
	return (d_y - d_mu)/muEta();
    }

    Rcpp::NumericVector glmerResp::wrkResp() const {
	return d_eta + wrkResids();
    }

    Rcpp::NumericVector glmerResp::sqrtWrkWt() const {
	Rcpp::NumericVector me = muEta();
	return sqrt(d_weights * me * me / variance());
    }

    double glmerResp::Laplace(double  ldL2,
			      double ldRX2,
			      double  sqrL) const{
	return ldL2 + sqrL + residDeviance();
    }
	
    double glmerResp::updateMu(const Rcpp::NumericVector& gamma) {
	d_eta = d_offset + gamma;
	NumericVector mmu = d_fam.linkInv(d_eta);
	std::copy(mmu.begin(), mmu.end(), d_mu.begin());
	return updateWrss();
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
	Rcpp::NumericVector gam = gamma + d_offset;
	double *gg = gam.begin();

	for (int p = 0; p < d_pnames.size(); p++) {
	    std::string pn(d_pnames[p]);
	    Rcpp::NumericVector pp = d_nlenv.get(pn);
	    std::copy(gg + n * p, gg + n * (p + 1), pp.begin());
	}
	Rcpp::NumericVector rr = d_nlmod.eval(SEXP(d_nlenv));
	if (rr.size() != n)
	    throw std::runtime_error("dimension mismatch");
	std::copy(rr.begin(), rr.end(), d_mu.begin());
	Rcpp::NumericMatrix rrg = rr.attr("gradient");
	std::copy(rrg.begin(), rrg.end(), d_sqrtXwt.begin());
	return updateWrss();
    }
}
