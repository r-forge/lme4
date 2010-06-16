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
	
	void  linkFun(Rcpp::NumericVector&, Rcpp::NumericVector const&);
	void  linkInv(Rcpp::NumericVector&, Rcpp::NumericVector const&);
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
	    varFuncs;		//< scalar variance functions
    
	static double epsilon, INVEPS, LTHRESH, MLTHRESH;

	static double cubef(double x) {return x * x * x;}
	static double identf(double x) {return x;}
	static double invderivf(double x) {return -1/(x * x);}
	static double inversef(double x) {return 1/x;}
	static double onef(double x) {return 1;}
	static double sqrf(double x) {return x * x;}
	static double twoxf(double x) {return 2 * x;}
	static double x1mxf(double x) {return std::max(epsilon, x * (1 - x));}
	
	static double finitePos(double x) { // truncate to [eps, 1/eps]
	    return std::max(epsilon, std::min(INVEPS, x));
	}
	static double finite01(double x) { // truncate to [eps, 1 - eps]
	    return std::max(epsilon, std::min(1. - epsilon, x));
	}
	
	static double logitLinkInv(double x) {
	    double tmp = finitePos(exp(x));
	    return tmp/(1 + tmp);
	}
	static double logitLink(double x) {
	    double xx = finite01(x);
	    return log(xx / (1 - xx));
	}
	static double logitMuEta(double x) {
	    return x1mxf(logitLinkInv(x));
	}
	
	static double probitLink(double x) {
	    return Rf_qnorm5(x, 0, 1, 1, 0);
	}
	static double probitLinkInv(double x) {
	    return finite01(Rf_pnorm5(x, 0, 1, 1, 0));
	}
	static double probitMuEta(double x) {
	    return finitePos(Rf_dnorm4(x, 0, 1, 0));
	}
    };
}
    
#endif /* LME4_GLMFAMILY_H */

