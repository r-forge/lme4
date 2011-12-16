// respModule.cpp: response modules using Eigen
//
// Copyright (C)       2011 Douglas Bates, Martin Maechler and Ben Bolker
//
// This file is part of lme4.

#include "respModule.h"
#include <cmath>

namespace lme4Eigen {
    using Eigen::VectorXd;

    using Rcpp::List;
    using Rcpp::NumericMatrix;
    using Rcpp::as;

    using std::copy;
    using std::string;
    using std::invalid_argument;

    typedef Eigen::Map<VectorXd>  MVec;

    lmResp::lmResp(SEXP y, SEXP weights, SEXP offset, SEXP mu, SEXP sqrtXwt,
		   SEXP sqrtrwt, SEXP wtres)
	: d_y(      as<MVec>(y)),
	  d_weights(as<MVec>(weights)),
	  d_offset( as<MVec>(offset)),
	  d_mu(     as<MVec>(mu)),
	  d_sqrtXwt(as<MVec>(sqrtXwt)),
	  d_sqrtrwt(as<MVec>(sqrtrwt)),
	  d_wtres(  as<MVec>(wtres)) {
	updateWrss();
    }

    double lmResp::updateMu(const VectorXd& gamma) {
	if (gamma.size() != d_offset.size())
	    throw invalid_argument("updateMu: Size mismatch");
	d_mu = d_offset + gamma;
	return updateWrss();
    }

    /** 
     * Update the wtres vector and return its sum of squares
     *   wtres <- sqrtrwt * (y - mu)
     *   return(wrss <- sum(wtres^2))
     *
     * @return Updated weighted residual sum of squares
     */
    double lmResp::updateWrss() {
	d_wtres = d_sqrtrwt.cwiseProduct(d_y - d_mu);
	d_wrss  = d_wtres.squaredNorm();
	return d_wrss;
    }

    void lmResp::setOffset(const VectorXd& oo) {
	if (oo.size() != d_offset.size())
	    throw invalid_argument("setOffset: Size mismatch");
	d_offset = oo;
    }

    void lmResp::setWeights(const VectorXd& ww) {
	if (ww.size() != d_weights.size())
	    throw invalid_argument("setWeights: Size mismatch");
	d_weights = ww;
    }

    lmerResp::lmerResp(SEXP y, SEXP weights, SEXP offset, SEXP mu,
		       SEXP sqrtXwt, SEXP sqrtrwt, SEXP wtres)
	: lmResp(y, weights, offset, mu, sqrtXwt, sqrtrwt, wtres),
	  d_reml(0) {
    }

    double lmerResp::Laplace(double ldL2, double ldRX2, double sqrL) const {
	double lnum = std::log(2.* M_PI * (d_wrss + sqrL));
	if (d_reml == 0) return ldL2 + d_y.size() * (1. + lnum - std::log(d_y.size()));
	double nmp = d_y.size() - d_reml;
	return ldL2 + ldRX2 + nmp * (1. + lnum - std::log(nmp));
    }

    void lmerResp::setReml(int rr) {
	if (rr < 0) throw invalid_argument("setReml: negative value for REML not meaningful");
	d_reml = rr;
    }
    
    glmResp::glmResp(List fam, SEXP y, SEXP weights, SEXP offset,
		     SEXP mu, SEXP sqrtXwt, SEXP sqrtrwt, SEXP wtres, SEXP eta, SEXP n)
	: lmResp(y, weights, offset, mu, sqrtXwt, sqrtrwt, wtres),
	  d_fam(fam),
	  d_eta(as<MVec>(eta)),
	  d_n(as<MVec>(n)) {
    }

    VectorXd glmResp::devResid() const {
	return d_fam.devResid(d_mu, d_weights, d_y);
    }

    VectorXd glmResp::muEta() const {
	return d_fam.muEta(d_eta);
    }

    VectorXd glmResp::variance() const {
	return d_fam.variance(d_mu);
    }

    VectorXd glmResp::wrkResids() const {
	return (d_y - d_mu).cwiseQuotient(muEta());
    }

    VectorXd glmResp::wrkResp() const {
	return (d_eta - d_offset) + wrkResids();
    }

    VectorXd glmResp::sqrtWrkWt() const {
	const VectorXd me(muEta());
	return d_weights.cwiseProduct(me).cwiseProduct(me).cwiseQuotient(variance()).cwiseSqrt();
    }

    double glmResp::Laplace(double ldL2, double ldRX2, double sqrL) const {
	return ldL2 + sqrL + resDev();
    }

    double glmResp::resDev() const {
	return devResid().sum();
    }

    double glmResp::updateMu(const VectorXd& gamma) {
	d_eta = d_offset + gamma; // lengths are checked here
	d_mu  = d_fam.linkInv(d_eta);
	return updateWrss();
    }

    double glmResp::updateWts() {
	d_sqrtrwt = d_weights.cwiseQuotient(variance()).cwiseSqrt();
	d_sqrtXwt = muEta().cwiseProduct(d_sqrtrwt);
	return updateWrss();
    }

    void glmResp::setN(const VectorXd& n) {
	if (n.size() != d_n.size())
	    throw invalid_argument("n size mismatch");
	d_n = n;
    }

    nlsResp::nlsResp(SEXP y, SEXP weights, SEXP offset, SEXP mu, SEXP sqrtXwt,
		     SEXP sqrtrwt, SEXP wtres, Language mm, Environment ee,
		     CharacterVector pp)
	: lmResp(y, weights, offset, mu, sqrtXwt, sqrtrwt, wtres),
	  d_nlenv(ee), d_nlmod(mm), d_pnames(pp) {
    }

    double nlsResp::Laplace(double ldL2, double ldRX2, double sqrL) const {
	double lnum = 2.* PI * (d_wrss + sqrL), n = d_y.size();
	return ldL2 + n * (1 + log(lnum / n));
    }

    double nlsResp::updateMu(VectorXd const &gamma) {
	int             n = d_y.size();
	VectorXd      gam = gamma + d_offset;
	const double  *gg = gam.data();

	for (int p = 0; p < d_pnames.size(); p++) {
	    string pn(d_pnames[p]);
	    NumericVector pp = d_nlenv.get(pn);
	    copy(gg + n * p, gg + n * (p + 1), pp.begin());
	}
	NumericVector  rr = d_nlmod.eval(SEXP(d_nlenv));
	if (rr.size() != n)
	    throw invalid_argument("dimension mismatch");
	copy(rr.begin(), rr.end(), d_mu.data());
	NumericMatrix  gr = rr.attr("gradient");
	copy(gr.begin(), gr.end(), d_sqrtXwt.data());
	return updateWrss();
    }

}
