#include "glmFamily.h"
#include <cmath>
#include <iostream>

static double sqr(double x) {return x * x;}
static double twox(double x) {return 2 * x;}
static double ident(double x) {return x;}
static double one(double x) {return 1;}

using namespace Rcpp;

fmap glmFamily::linvs = fmap();
fmap glmFamily::lnks = fmap();
fmap glmFamily::muEtas = fmap();
fmap glmFamily::varFuncs = fmap();

glmFamily::glmFamily(List ll) : lst(ll) {
    if (as<std::string>(ll.attr("class")) != "family")
	Rf_error("glmFamily only from list of (S3) class \"family\"");
    CharacterVector fam = ll["family"], llink = ll["link"];
    char *fpt = fam[0];
    Rprintf("Assigned fam of length %d %d\n", fam.size());
    if (fam.size() > 0) Rprintf("fam[0] = %s", fpt);
    std::string ff( fam[0] ), lll( llink[0] );
    std::cout << "fam length = " << fam.size() << 
	", fam[0] = " << ff << std::endl;
    std::cout << "llink length = " << llink.size() <<
	", llink[0] = " << lll << std::endl;
    family = ff;
    link = lll;
    std::cout << family << " " << link << std::endl;

    if (!lnks.count("identity")) { // initialize the static maps
	lnks["log"] = &log;
	muEtas["log"] = linvs["log"] = &exp;
	lnks["sqrt"] = &sqrt;
	linvs["sqrt"] = &sqr;
	muEtas["sqrt"] = &twox;
	lnks["identity"] = linvs["identity"] = &ident;
	muEtas["identity"] = &one;
	varFuncs["gaussian"] = &one;
    }
}

void glmFamily::linkFun(NumericVector eta, const NumericVector mu) {
    if (lnks.count(link)) {
	std::transform(mu.begin(), mu.end(), eta.begin(), lnks[link]);
    } else {
	Function linkfun = lst["linkfun"];
	NumericVector ans(linkfun(mu));
	std::copy(ans.begin(), ans.end(), eta.begin());
    }
}

void glmFamily::linkInv(NumericVector mu, const NumericVector eta) {
    if (linvs.count(link)) {
	std::transform(eta.begin(), eta.end(), mu.begin(), linvs[link]);
    } else {
	Function linkinv = lst["linkinv"];
	NumericVector ans(linkinv(eta));
	std::copy(ans.begin(), ans.end(), mu.begin());
    }
}

void glmFamily::muEta(NumericVector mueta, const NumericVector eta) {
    if (muEtas.count(link)) {
	std::transform(eta.begin(), eta.end(), mueta.begin(), muEtas[link]);
    } else {
	Function mm = lst["mu.eta"];
	NumericVector ans(mm(eta));
	std::copy(ans.begin(), ans.end(), mueta.begin());
    }
}

void glmFamily::variance(NumericVector vv, const NumericVector eta) {
    if (varFuncs.count(link)) {
	std::transform(eta.begin(), eta.end(), vv.begin(), varFuncs[link]);
    } else {
	Function mm = lst["variance"];
	NumericVector ans(mm(eta));
	std::copy(ans.begin(), ans.end(), vv.begin());
    }
}

double glmFamily::devResid(const NumericVector mu, NumericVector weights,
			   const NumericVector y) {
    return 0;
}

SEXP glmFamily::show() {
    return wrap(List::create(
		    _["family"] = family,
		    _["link"] = link
//		    _["compiledLink"] = bool(lnks.count(link)),
//		    _["compiledLinv"] = bool(linvs.count(link)),
//		    _["compiledmuEta"] = bool(muEtas.count(link)),
//		    _["compiledVar"] = bool(varFuncs.count(link))
		    ));
}

extern "C" SEXP
family_show(SEXP family) {
    List ll(family);
    glmFamily ff(ll);
    return ff.show();
}
