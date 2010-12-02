#include "RcppModels.h"

using namespace Rcpp;

// Irwls::Irwls(Rcpp::NumericMatrix xR, Rcpp::NumericVector yR) :
//     x(xR), y(yR) {
//     int n = x.nrow(), p = x.ncol();
//     if (y.size() != n)
// 	throw std::runtime_error("incompatible dimensions");
//     X = arma::mat(x.begin(), n, p, false, true);
//     Y = arma::vec(y.begin(), n, false, true);
// }

// arma::vec Irwls::fit(Rcpp::NumericVector wR) {
//     int n = x.nrow();
//     if (wR.size() != n)
// 	throw std::runtime_error("length(weights) != nrow(X)");
//     wrt = sqrt(wR);
//     arma::vec W(wrt.begin(), n, false, true);
//     return arma::solve(diagmat(W) * X, W % Y);
// }

namespace glm {
    lmResp::lmResp(Rcpp::NumericVector y)
	throw (std::runtime_error)
	: d_y(y), d_weights(y.size(), 1.0), d_offset(y.size()),
	  d_mu(y.size()), d_sqrtrwt(y.size()),
	  d_sqrtXwt(y.size(), 1.0) {
	init();
    }

    lmResp::lmResp(Rcpp::NumericVector y, Rcpp::NumericVector weights)
	throw (std::runtime_error)
	: d_y(y), d_weights(weights), d_offset(y.size()),
	  d_mu(y.size()), d_sqrtrwt(y.size()),
	  d_sqrtXwt(y.size(), 1) {
	if (weights.size() != y.size())
	    throw std::runtime_error(
		"lengths of y and wts must agree");
	init();
    }

    lmResp::lmResp(Rcpp::NumericVector y, Rcpp::NumericVector weights,
	Rcpp::NumericVector offset) throw (std::runtime_error)
	: d_y(y), d_weights(weights), d_offset(offset),
	  d_mu(y.size()), d_sqrtrwt(y.size()),
	  d_sqrtXwt(y.size(), 1) {
	int nn = y.size();
	if (weights.size() != nn || offset.size() != nn)
	    throw std::runtime_error(
		"lengths of y, weights and offset must agree");
	init();
    }

    void lmResp::init() {
	d_sqrtrwt = sqrt(d_weights);
	std::copy(d_sqrtrwt.begin(), d_sqrtrwt.end(), d_sqrtXwt.begin());
    }

    double lmResp::updateMu(const Rcpp::NumericVector& gamma) {
	d_mu = d_offset + gamma;
	return updateWrss();
    }
    
    double lmResp::updateWrss() {
	d_wtres = d_sqrtrwt * (d_y - d_mu);
	d_wrss = std::inner_product(d_wtres.begin(), d_wtres.end(),
				    d_wtres.begin(), double());
	return d_wrss;
    }	

    glmResp::glmResp(Rcpp::List fam, Rcpp::NumericVector y)
	throw (std::runtime_error)
	: lmResp(y), d_fam(fam), d_eta(y.size()) {
    }

    glmResp::glmResp(Rcpp::List fam, Rcpp::NumericVector y,
		     Rcpp::NumericVector weights)
	throw (std::runtime_error)
	: lmResp(y, weights), d_fam(fam), d_eta(y.size()) {
    }

    glmResp::glmResp(Rcpp::List fam, Rcpp::NumericVector y,
		     Rcpp::NumericVector weights,
		     Rcpp::NumericVector offset)
	throw (std::runtime_error) 
	: lmResp(y, weights, offset), d_fam(fam),
	  d_eta(y.size()) {
    }

    glmResp::glmResp(Rcpp::List fam, Rcpp::NumericVector y,
		     Rcpp::NumericVector weights,
		     Rcpp::NumericVector offset,
		     Rcpp::NumericVector n)
	throw (std::runtime_error)
 	: lmResp(y, weights, offset), d_fam(fam), d_n(n),
	  d_eta(y.size()) {
	if (n.size() != y.size())
	    throw std::runtime_error("lengths of y and n must agree");
    }

    glmResp::glmResp(Rcpp::List fam, Rcpp::NumericVector y,
		     Rcpp::NumericVector weights,
		     Rcpp::NumericVector offset,
		     Rcpp::NumericVector n, Rcpp::NumericVector eta)
	throw (std::runtime_error) 
	: lmResp(y, weights, offset), d_fam(fam), d_n(n),
	  d_eta(eta) {
	int nn = y.size();
	if (n.size() != nn || eta.size() != nn )
	    throw std::runtime_error(
		"lengths of y, n and eta must agree");
    }

    Rcpp::NumericVector glmResp::devResid() const {
	return d_fam.devResid(d_mu, d_weights, d_y);
    }

    Rcpp::NumericVector glmResp::muEta() const {
	return d_fam.muEta(d_eta);
    }

    Rcpp::NumericVector glmResp::variance() const {
	return d_fam.variance(d_mu);
    }
    

    double glmResp::residDeviance() const {
	return sum(devResid());
    }

    double glmResp::updateWts() {
	d_sqrtrwt = sqrt(d_weights/d_fam.variance(d_mu));
	NumericVector muEta = d_fam.muEta(d_eta);
	std::transform(muEta.begin(), muEta.end(), d_sqrtrwt.begin(),
		       d_sqrtXwt.begin(), std::multiplies<double>());
	return updateWrss();
    }

    Rcpp::NumericVector glmResp::wrkResids() const {
	return (d_y - d_mu)/muEta();
    }

    Rcpp::NumericVector glmResp::wrkResp() const {
	return d_eta + wrkResids();
    }

    Rcpp::NumericVector glmResp::sqrtWrkWt() const {
	Rcpp::NumericVector me = muEta();
	return sqrt(d_weights * me * me / variance());
    }

    double glmResp::updateMu(const Rcpp::NumericVector& gamma) {
	d_eta = d_offset + gamma;
	NumericVector mmu = d_fam.linkInv(d_eta);
	std::copy(mmu.begin(), mmu.end(), d_mu.begin());
	return updateWrss();
    }
    
    predMod::predMod(int p)
	: d_coef(p), d_coef0(p), d_delta(p),
	  a_coef (d_coef.begin(),  p, false, true),
	  a_coef0(d_coef0.begin(), p, false, true),
	  a_delta(d_delta.begin(), p, false, true) {
    }

    predMod::predMod(Rcpp::NumericVector coef0)
	: d_coef(coef0.size()), d_coef0(coef0), d_delta(coef0.size()),
	  a_coef(d_coef.begin(), d_coef.size(), false, true),
	  a_coef0(d_coef0.begin(), d_coef0.size(), false, true),
	  a_delta(d_delta.begin(), d_delta.size(), false, true) {
    }
    
   void predMod::setCoef0(const Rcpp::NumericVector& cc)
	throw (std::runtime_error) {
	if (cc.size() != d_coef0.size())
	    throw std::runtime_error("size mismatch in setCoef0");
	std::copy(cc.begin(), cc.end(), d_coef0.begin());
    }

    void predMod::installCoef0() {
	std::copy(d_coef.begin(), d_coef.end(), d_coef0.begin());
    }

    densePred::densePred(Rcpp::NumericMatrix x)
	throw (std::runtime_error)
	: predMod(x.ncol()), d_X(x),
	  a_X(x.begin(), x.nrow(), x.ncol(), false, true) {
    }

    densePred::densePred(Rcpp::NumericMatrix x,
			 Rcpp::NumericVector coef0)
	throw (std::runtime_error)
	: predMod(coef0), d_X(x),
	  a_X(x.begin(), x.nrow(), x.ncol(), false, true) {
	if (d_coef0.size() != d_X.ncol())
	    throw std::runtime_error("length(coef0) != ncol(X)");
    }

    void densePred::installCoef0() {
	std::copy(d_coef.begin(), d_coef.end(), d_coef0.begin());
    }

   void densePred::setCoef0(const Rcpp::NumericVector& cc)
	throw (std::runtime_error) {
	if (cc.size() != d_coef0.size())
	    throw std::runtime_error("size mismatch in setCoef0");
	std::copy(cc.begin(), cc.end(), d_coef0.begin());
    }

    
    // Returns the linear predictor from the full step
    Rcpp::NumericVector
    densePred::gamma(const Rcpp::NumericVector& xwts,
		     const Rcpp::NumericVector& wtres)
	throw (std::runtime_error) {
	int n = d_X.nrow();
	if (xwts.size() != n || wtres.size() != n)
	    throw
		std::runtime_error("length(xwts) or length(wtres) != nrow(X)");
	arma::vec a_xwts = arma::vec(xwts.begin(), n, false, true),
	    a_wtres = arma::vec(wtres.begin(), n, false, true);
	arma::vec tmp = solve(diagmat(a_xwts) * a_X, a_wtres);
	std::copy(tmp.begin(), tmp.end(), a_delta.begin());
	d_coef = d_coef0 + d_delta;
	return NumericVector(wrap(a_X * a_coef));
    }
    
    // Returns the linear predictor from a step of fac
    Rcpp::NumericVector densePred::dgamma(double fac) {
	return NumericVector(wrap(a_X * (a_coef0 + fac * a_delta)));
    }    
}    
