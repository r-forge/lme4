#include "glmFamily.h"
#include <cmath>
#include <limits>

using namespace Rcpp;
using namespace std;

namespace glm {
    // Establish the values for the class constants
    double glmFamily::epsilon = numeric_limits<double>::epsilon();
    
    // initialize the function maps (i.e. associative arrays of functions)
    drmap glmFamily::devRes = drmap();

    fmap glmFamily::linvs = fmap();
    fmap glmFamily::lnks = fmap();
    fmap glmFamily::muEtas = fmap();
    fmap glmFamily::varFuncs = fmap();
    
    void glmFamily::initMaps() {
	// initialize the static maps.  The identity link is
	// guaranteed to be initialized if any maps are initialized
	    lnks["log"]                  = &log;
	    muEtas["log"] = linvs["log"] = &exp;
	    
	    lnks["sqrt"]                 = &sqrt;
	    linvs["sqrt"]                = &sqrf;
	    muEtas["sqrt"]               = &twoxf;
	    
	    lnks["identity"]             = &identf;
	    linvs["identity"]            = &identf;
	    muEtas["identity"]           = &onef;
	    
	    lnks["inverse"]              = &inversef;
	    linvs["inverse"]             = &inversef;
	    muEtas["inverse"]            = &invderivf;
	    
	    lnks["logit"]                = &logitLink;
	    linvs["logit"]               = &logitLinkInv;
	    muEtas["logit"]              = &logitMuEta;
	    
	    lnks["probit"]               = &probitLink;
	    linvs["probit"]              = &probitLinkInv;
	    muEtas["probit"]             = &probitMuEta;
	    
//	    lnks["cloglog"]              = &cloglogLink;
	    linvs["cloglog"]             = &cloglogLinkInv;
	    muEtas["cloglog"]            = &cloglogMuEta;
	    
	    devRes["Gamma"]              = &GammaDevRes;
	    varFuncs["Gamma"]            = &sqrf;   // x^2

	    devRes["binomial"]           = &BinomialDevRes;
	    varFuncs["binomial"]         = &x1mxf;  // x * (1 - x)

	    devRes["gaussian"]           = &GaussianDevRes;
	    varFuncs["gaussian"]         = &onef;   // 1

	    varFuncs["inverse.gaussian"] = &cubef;  // x^3

	    devRes["poisson"]            = &PoissonDevRes;
	    varFuncs["poisson"]          = &identf; // x
    }
    
    glmFamily::glmFamily(List ll) throw (std::runtime_error)
	: lst(ll),
	  // d_family(as<std::string>(wrap(ll["family"]))),
	  // d_link(as<std::string>(wrap(ll["link"]))),
// I haven't been able to work out an expression to initialize the
// Functions from list components.  This is a placeholder until I can
// do so.
	  d_devRes("c"), d_linkfun("c"), d_linkinv("c"),
	  d_muEta("c"), d_variance("c") {
	  // d_devRes(wrap(ll["dev.resids"])),
	  // d_linkfun(wrap(ll["linkfun"])),
	  // d_linkinv(wrap(ll["linkinv"])),
	  // d_muEta(wrap(ll["mu.eta"])),
	  // d_variance(wrap(ll["variance"])) {
#if RCPP_VERSION > Rcpp_Version(0,8,9)
	if (!lst.inherits("family"))
	    throw std::runtime_error("glmFamily requires a list of (S3) class \"family\"");
#else
	std::string fam = as<std::string>(lst.attr("class"));
	if (fam != "family")
	    throw std::runtime_error("glmFamily requires a list of (S3) class \"family\"");
#endif
 	CharacterVector ff = lst["family"], lnk = lst["link"];
 	d_family = as<std::string>(ff);
 	d_link = as<std::string>(lnk);
 	d_linkinv = ll["linkinv"];
 	d_linkfun = ll["linkfun"];
 	d_muEta = ll["mu.eta"];
 	d_variance = ll["variance"];
 	d_devRes = ll["dev.resids"];

	if (!lnks.count("identity")) initMaps();
    }

    Rcpp::NumericVector
    glmFamily::linkFun(Rcpp::NumericVector const &mu) const {
	if (lnks.count(d_link))
	    return NumericVector::import_transform(mu.begin(), mu.end(), lnks[d_link]);
#if RCPP_VERSION > Rcpp_Version(0,8,9)
	return d_linkfun(mu);
#else
	Function linkfun = ((const_cast<glmFamily*>(this))->lst)["linkfun"];
	return linkfun(mu);
#endif
    }
    
    Rcpp::NumericVector
    glmFamily::linkInv(Rcpp::NumericVector const &eta) const {
	if (linvs.count(d_link))
	    return NumericVector::import_transform(eta.begin(), eta.end(), linvs[d_link]);
#if RCPP_VERSION > Rcpp_Version(0,8,9)
	return d_linkinv(eta);
#else
	Function linkinv = ((const_cast<glmFamily*>(this))->lst)["linkinv"];
	return linkinv(eta);
#endif
    }

    Rcpp::NumericVector
    glmFamily::muEta(Rcpp::NumericVector const &eta) const {
	if (muEtas.count(d_link))
	    return NumericVector::import_transform(eta.begin(), eta.end(), muEtas[d_link]);
#if RCPP_VERSION > Rcpp_Version(0,8,9)
	return d_muEta(eta);
#else
	Function muEta = ((const_cast<glmFamily*>(this))->lst)["mu.eta"];
	return muEta(eta);
#endif
    }
    
    Rcpp::NumericVector
    glmFamily::variance(Rcpp::NumericVector const &mu) const {
	if (varFuncs.count(d_link))
	    return NumericVector::import_transform(mu.begin(), mu.end(), varFuncs[d_link]);
#if RCPP_VERSION > Rcpp_Version(0,8,9)
	return d_variance(mu);
#else
	Function varFunc = ((const_cast<glmFamily*>(this))->lst)["variance"];
	return varFunc(mu);
#endif
    }
    
    Rcpp::NumericVector
    glmFamily::devResid(Rcpp::NumericVector const &mu,
			Rcpp::NumericVector const &weights,
			Rcpp::NumericVector const &y) const {
	if (devRes.count(d_family)) {
	    double (*f)(double, double, double) = devRes[d_family];
	    int nobs = mu.size();
	    NumericVector ans(nobs);
	    double *mm = mu.begin(), *ww = weights.begin(),
		*yy = y.begin(), *aa = ans.begin();
	    for (int i = 0; i < nobs; i++)
		aa[i] = f(yy[i], mm[i], ww[i]);
	    return ans;
	}
#if RCPP_VERSION > Rcpp_Version(0,8,9)
	return d_devRes(y, mu, weights);
#else
	Function devRes = ((const_cast<glmFamily*>(this))->lst)["dev.resids"];
	return devRes(y, mu, weights);
#endif
    }
}
