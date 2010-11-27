#include "glmFamily.h"
#include <cmath>
#include <limits>

using namespace Rcpp;
using namespace std;

namespace glm {
    // Establish the values for the class constants
    double glmFamily::epsilon = numeric_limits<double>::epsilon();
//    double glmFamily::INVEPS = 1. / glmFamily::epsilon;
//    double glmFamily::LTHRESH = 30.;
//    double glmFamily::MLTHRESH = -30.;

    // initialize the function maps to an empty map
    fmap glmFamily::linvs = fmap();
    fmap glmFamily::lnks = fmap();
    fmap glmFamily::muEtas = fmap();
    fmap glmFamily::varFuncs = fmap();

    glmFamily::glmFamily(SEXP ll) : lst(ll) {
	if (as<string>(lst.attr("class")) != "family")
	    Rf_error("glmFamily only from list of (S3) class \"family\"");

	CharacterVector fam = lst["family"], llink = lst["link"];
	char *pt = fam[0]; family = string(pt);
	pt = llink[0]; link = string(pt);

	// initialize the static maps.  The identity link is
	// guaranteed to be initialized if any maps are initialized
	if (!lnks.count("identity")) {
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

	    varFuncs["Gamma"]            = &sqrf;   // x^2
	    varFuncs["binomial"]         = &x1mxf;  // x * (1 - x)
	    varFuncs["inverse.gaussian"] = &cubef;  // x^3
	    varFuncs["gaussian"]         = &onef;   // 1
	    varFuncs["poisson"]          = &identf; // x
	}
    }

    Rcpp::NumericVector
    glmFamily::linkFun(Rcpp::NumericVector const &mu) const {
	if (lnks.count(link)) {	// sapply the known scalar function
	    return NumericVector::import_transform(mu.begin(), mu.end(), lnks[link]);
	} else {		// use the R function
	    Function linkfun = ((const_cast<glmFamily*>(this))->lst)["linkfun"];
	    // The const_cast is needed so that this member function
	    // can be const and also use the extraction of a list
	    // component.
	    return linkfun(mu);
	}
    }

    Rcpp::NumericVector
    glmFamily::linkInv(Rcpp::NumericVector const &eta) const {
	if (linvs.count(link)) {
	    return NumericVector::import_transform(eta.begin(), eta.end(), linvs[link]);
	} else {
	    Function linkinv = ((const_cast<glmFamily*>(this))->lst)["linkinv"];
	    return linkinv(eta);
	}
    }

    Rcpp::NumericVector
    glmFamily::muEta(Rcpp::NumericVector const &eta) const {
	if (muEtas.count(link)) {
	    return NumericVector::import_transform(eta.begin(), eta.end(), muEtas[link]);
	}
	Function mu_eta = ((const_cast<glmFamily*>(this))->lst)["mu.eta"];
	return mu_eta(eta);
    }

    Rcpp::NumericVector
    glmFamily::variance(Rcpp::NumericVector const &mu) const {
	if (varFuncs.count(link)) {
	    return NumericVector::import_transform(mu.begin(), mu.end(), varFuncs[link]);
	}
	Function vv = ((const_cast<glmFamily*>(this))->lst)["variance"];
	return vv(mu);
    }

    Rcpp::NumericVector
    glmFamily::devResid(Rcpp::NumericVector const &mu,
			Rcpp::NumericVector const &weights,
			Rcpp::NumericVector const &y) const {
	Function devres = ((const_cast<glmFamily*>(this))->lst)["dev.resids"];
	return devres(y, mu, weights);
    }
}


RCPP_MODULE(glm) {

    class_<glm::glmFamily>( "glmFamily" )

//? not yet valid: .constructor(init_1<List>())

    .method("linkFun",        &glm::glmFamily::linkFun)
    .method("linkInv",        &glm::glmFamily::linkInv)
    .method("muEta",          &glm::glmFamily::muEta)
    .method("devResid",       &glm::glmFamily::devResid)
    .method("variance",       &glm::glmFamily::variance)
    ;

}
