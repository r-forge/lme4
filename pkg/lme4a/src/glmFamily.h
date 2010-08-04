// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
#ifndef LME4_GLMFAMILY_H
#define LME4_GLMFAMILY_H

#include <Rcpp.h>
#include <Rmath.h>

namespace glm {
    // Map (associative array) of functions returning double from a double
    typedef std::map<std::string, double(*)(double)> fmap;

    class glmFamily {
	std::string family, link; // as in the R glm family
	Rcpp::List           lst; // original list from R
    public:
	glmFamily(SEXP);
				// Application of functions from the family
				// The scalar transformations use
				// compiled code when available 
	Rcpp::NumericVector  linkFun(Rcpp::NumericVector const&) const;
	Rcpp::NumericVector  linkInv(Rcpp::NumericVector const&) const;
	Rcpp::NumericVector devResid(
	    Rcpp::NumericVector const&,
	    Rcpp::NumericVector const&,
	    Rcpp::NumericVector const&) const;
	Rcpp::NumericVector    muEta(Rcpp::NumericVector const&) const;
	Rcpp::NumericVector variance(Rcpp::NumericVector const&) const;
    private:
	// Class (as opposed to instance) members storing the scalar functions
	static fmap
 	    devs,			//< scalar deviance functions
	    lnks,			//< scalar link functions
	    linvs,			//< scalar linkinv functions
	    muEtas,			//< scalar muEta functions
	    varFuncs;			//< scalar variance functions
    
	// Thresholds common to the class (FIXME: we should eliminate
	// most of these)
	static double epsilon, INVEPS, LTHRESH, MLTHRESH;

	// Scalar functions that will used in transform applications
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
	static inline double   logitLinkInv(double x) {
	    return Rf_plogis(x, 0., 1., 1, 0);
	}
	static inline double      logitLink(double x) {
	    return Rf_qlogis(x, 0., 1., 1, 0);
	}
	static inline double     logitMuEta(double x) {
	    return Rf_dlogis(x, 0., 1., 0);
	}
	static inline double  probitLinkInv(double x) {
	    return Rf_pnorm5(x, 0., 1., 1, 0);
	}
	static inline double     probitLink(double x) {
	    return Rf_qnorm5(x, 0., 1., 1, 0);
	}
	static inline double    probitMuEta(double x) {
	    return Rf_dnorm4(x, 0., 1., 0);
	}
	static inline double
	pgumbel(double q, double loc, double scale, int lower_tail) {
	    q = (q - loc) / scale;
	    q = -exp(q);
	    return lower_tail ? -expm1(q) : exp(q);
	}
	static inline double cloglogLinkInv(double x) {
	    return pgumbel(x, 0., 1., 1);
	}
	static inline double
	dgumbel(double x, double loc, double scale, int give_log) {
	    x = (x - loc) / scale;
	    x = -exp(-x) - x - log(scale);
	    return give_log ? x : exp(x);
	}
	static inline double   cloglogMuEta(double x) {
	    return dgumbel(x, 0., 1., 0);
	}
    };
}
    
#endif /* LME4_GLMFAMILY_H */

