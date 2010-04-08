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

	varFuncs["gamma"] = &sqrf;
	varFuncs["gaussian"] = &onef;
	varFuncs["poisson"] = &identf;
    }
}

void
glmFamily::linkFun(NumericVector& eta, const NumericVector& mu) {
    if (lnks.count(link)) {
	std::transform(mu.begin(), mu.end(), eta.begin(), lnks[link]);
    } else {
	Function linkfun = lst["linkfun"];
	NumericVector ans = linkfun(mu);
	std::copy(ans.begin(), ans.end(), eta.begin());
    }
}

void
glmFamily::linkInv(NumericVector& mu, const NumericVector& eta) {
    if (linvs.count(link)) {
	std::transform(eta.begin(), eta.end(), mu.begin(), linvs[link]);
    } else {
	Function linkinv = lst["linkinv"];
	NumericVector ans = linkinv(eta);
	std::copy(ans.begin(), ans.end(), mu.begin());
    }
}

void
glmFamily::muEta(NumericVector& mueta, const NumericVector& eta) {
    if (muEtas.count(link)) {
	std::transform(eta.begin(), eta.end(), mueta.begin(), muEtas[link]);
    } else {
	Function mm = lst["mu.eta"];
	NumericVector ans = mm(eta);
	std::copy(ans.begin(), ans.end(), mueta.begin());
    }
}

void
glmFamily::variance(NumericVector& vv, const NumericVector& eta) {
    if (varFuncs.count(link)) {
	std::transform(eta.begin(), eta.end(), vv.begin(), varFuncs[link]);
    } else {
	Function mm = lst["variance"];
	NumericVector ans = mm(eta);
	std::copy(ans.begin(), ans.end(), vv.begin());
    }
}

double
glmFamily::devResid(const NumericVector& mu, const NumericVector& weights,
		    const NumericVector& y) {
    Function devres = lst["dev.resids"];
    NumericVector dd = devres(y, mu, weights);
    return std::accumulate(dd.begin(), dd.end(), double());
}


SEXP
glmFamily::show() {
    return List::create(
	_["family"] = family,
	_["link"] = link,
	_["compiledLink"] = bool(lnks.count(link)),
	_["compiledLinv"] = bool(linvs.count(link)),
	_["compiledmuEta"] = bool(muEtas.count(link)),
	_["compiledVar"] = bool(varFuncs.count(family))
	);
}

extern "C" SEXP 
family_link(SEXP family, SEXP pmu) {
    NumericVector mu(pmu);
    NumericVector eta(mu.size());

    glmFamily(family).linkFun(eta, mu);
    return eta;
}

extern "C" SEXP 
family_linkinv(SEXP family, SEXP peta) {
    NumericVector eta(peta);
    NumericVector mu(eta.size());

    glmFamily(family).linkInv(mu, eta);
    return mu;
}

extern "C" SEXP 
family_muEta(SEXP family, SEXP peta) {
    NumericVector eta(peta);
    NumericVector muEta(eta.size());

    glmFamily(family).muEta(muEta, eta);
    return muEta;
}

extern "C" SEXP
family_show(SEXP family) {
    return glmFamily(family).show();
}
