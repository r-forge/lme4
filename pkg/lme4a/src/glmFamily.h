// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
#ifndef LME4_GLMFAMILY_H
#define LME4_GLMFAMILY_H

#include <Rcpp.h>
#include <Rmath.h>

namespace glm {
    typedef std::map<std::string, double(*)(double)> fmap;

    class glmFamily {
	std::string family, link;
	Rcpp::List           lst;
    public:
	glmFamily(SEXP);
	
	Rcpp::NumericVector  linkFun(Rcpp::NumericVector const&) const;
	Rcpp::NumericVector  linkInv(Rcpp::NumericVector const&) const;
	Rcpp::NumericVector devResid(
	    Rcpp::NumericVector const&,
	    Rcpp::NumericVector const&,
	    Rcpp::NumericVector const&) const;
	Rcpp::NumericVector    muEta(Rcpp::NumericVector const&) const;
	Rcpp::NumericVector variance(Rcpp::NumericVector const&) const;
    private:
	static fmap
 	    devs,			//< scalar deviance functions
	    lnks,			//< scalar link functions
	    linvs,			//< scalar linkinv functions
	    muEtas,			//< scalar muEta functions
	    varFuncs;			//< scalar variance functions
    
	static double epsilon, INVEPS, LTHRESH, MLTHRESH;

	static inline double         cubef(double x) {return x * x * x;}
	static inline double        identf(double x) {return x;}
	static inline double     invderivf(double x) {return -1/(x * x);}
	static inline double      inversef(double x) {return 1/x;}
	static inline double          onef(double x) {return 1;}
	static inline double          sqrf(double x) {return x * x;}
	static inline double         twoxf(double x) {return 2 * x;}
	static inline double         x1mxf(double x) {return std::max(epsilon, x * (1 - x));}
	static inline double     finitePos(double x) { // truncate to [eps, 1/eps]
	    return std::max(epsilon, std::min(INVEPS, x));
	}
	static inline double      finite01(double x) { // truncate to [eps, 1 - eps]
	    return std::max(epsilon, std::min(1. - epsilon, x));
	}
	static inline double  logitLinkInv(double x) {
	    return Rf_plogis(x, 0., 1., 1, 0);
	}
	static inline double     logitLink(double x) {
	    return Rf_qlogis(x, 0., 1., 1, 0);
	}
	static inline double    logitMuEta(double x) {
	    return Rf_dlogis(x, 0., 1., 0);
	}
	static inline double probitLinkInv(double x) {
	    return Rf_pnorm5(x, 0., 1., 1, 0);
	}
	static inline double    probitLink(double x) {
	    return Rf_qnorm5(x, 0., 1., 1, 0);
	}
	static inline double   probitMuEta(double x) {
	    return Rf_dnorm4(x, 0., 1., 0);
	}
    };
}
    
#endif /* LME4_GLMFAMILY_H */

