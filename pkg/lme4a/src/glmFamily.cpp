#include "glmFamily.h"
#include <cmath>
#include <limits>

using namespace Rcpp;

double glmFamily::epsilon = std::numeric_limits<double>::epsilon();
double glmFamily::INVEPS = 1. / glmFamily::epsilon;
double glmFamily::LTHRESH = 30.;
double glmFamily::MLTHRESH = -30.;

fmap glmFamily::linvs = fmap();
fmap glmFamily::lnks = fmap();
fmap glmFamily::muEtas = fmap();
fmap glmFamily::varFuncs = fmap();

glmFamily::glmFamily(SEXP ll) : lst(ll) {
    if (as<std::string>(lst.attr("class")) != "family")
	Rf_error("glmFamily only from list of (S3) class \"family\"");

    CharacterVector fam = lst["family"], llink = lst["link"];
    char *pt = fam[0]; family = std::string(pt);
    pt = llink[0]; link = std::string(pt);

    if (!lnks.count("identity")) { // initialize the static maps
	lnks["log"] = &log;
	muEtas["log"] = linvs["log"] = &exp;

	lnks["sqrt"] = &sqrt;
	linvs["sqrt"] = &sqrf;
	muEtas["sqrt"] = &twoxf;

	lnks["identity"] = linvs["identity"] = &identf;
	muEtas["identity"] = &onef;

	lnks["inverse"] = linvs["inverse"] = &inversef;
	muEtas["inverse"] = &invderivf;

	lnks["logit"] = &logitLink;
	linvs["logit"] = &logitLinkInv;
	muEtas["logit"] = &logitMuEta;

	lnks["probit"] = &probitLink;
	linvs["probit"] = &probitLinkInv;
	muEtas["probit"] = &probitMuEta;

	varFuncs["Gamma"] = &sqrf;             // x^2
	varFuncs["binomial"] = &x1mxf;         // x * (1 - x)
	varFuncs["inverse.gaussian"] = &cubef; // x^3
	varFuncs["gaussian"] = &onef;	       // 1
	varFuncs["poisson"] = &identf;	       // x
    }
}

void
glmFamily::linkFun(NumericVector &eta, NumericVector const &mu) {
    if (lnks.count(link)) {
	std::transform(mu.begin(), mu.end(), eta.begin(), lnks[link]);
    } else {
	Function linkfun = lst["linkfun"];
	NumericVector ans = linkfun(mu);
	std::copy(ans.begin(), ans.end(), eta.begin());
    }
}

void
glmFamily::linkInv(NumericVector &mu, NumericVector const &eta) {
    if (linvs.count(link)) {
	std::transform(eta.begin(), eta.end(), mu.begin(), linvs[link]);
    } else {
	Function linkinv = lst["linkinv"];
	NumericVector ans = linkinv(eta);
	std::copy(ans.begin(), ans.end(), mu.begin());
    }
}

void
glmFamily::muEta(NumericVector &mueta, NumericVector const &eta) {
    if (muEtas.count(link)) {
	std::transform(eta.begin(), eta.end(), mueta.begin(),
		       muEtas[link]);
    } else {
	Function mm = lst["mu.eta"];
	NumericVector ans = mm(eta);
	std::copy(ans.begin(), ans.end(), mueta.begin());
    }
}

void
glmFamily::variance(NumericVector &vv, NumericVector const &mu) {
    if (varFuncs.count(link)) {
	std::transform(mu.begin(), mu.end(), vv.begin(), varFuncs[link]);
    } else {
	Function mm = lst["variance"];
	NumericVector ans = mm(mu);
	std::copy(ans.begin(), ans.end(), vv.begin());
    }
}

NumericVector
glmFamily::devResid(NumericVector const &mu,
		    NumericVector const &weights,
		    NumericVector const &y) {
    Function devres = lst["dev.resids"];
    NumericVector dd = devres(y, mu, weights);
    return dd;
}
