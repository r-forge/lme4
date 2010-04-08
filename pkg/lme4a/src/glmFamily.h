// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
#ifndef LME4_GLMFAMILY_H
#define LME4_GLMFAMILY_H

#include <Rcpp.h>
#include <Rmath.h>

typedef std::map<std::string, double(*)(double)> fmap;

/// generalized linear model family
class glmFamily {
public:
    glmFamily(SEXP);
    void linkFun(Rcpp::NumericVector& eta,
		 const Rcpp::NumericVector& mu);
    void linkInv(Rcpp::NumericVector& mu,
		 const Rcpp::NumericVector& eta);
    void muEta(Rcpp::NumericVector& mueta,
	       const Rcpp::NumericVector& eta);
    void variance(Rcpp::NumericVector& vv,
		  const Rcpp::NumericVector& eta);
    double devResid(const Rcpp::NumericVector& mu,
		    const Rcpp::NumericVector& weights,
		    const Rcpp::NumericVector& y);
    SEXP show();
    
    std::string family, link;
private:
    static fmap
	lnks,			//< scalar link functions
	linvs,			//< scalar linkinv functions
	muEtas,			//< scalar muEta functions
	varFuncs;		//< scalar variance functions
    Rcpp::List lst;
    
    static double epsilon, INVEPS, LTHRESH, MLTHRESH;

    static double sqrf(double x) {return x * x;}
    static double twoxf(double x) {return 2 * x;}
    static double identf(double x) {return x;}
    static double onef(double x) {return 1;}
    static double inversef(double x) {return 1/x;}
    static double invderivf(double x) {return -1/(x * x);}

    static double safeExp(double x) {
	return (x < MLTHRESH) ? epsilon : ((x > LTHRESH) ? INVEPS : exp(x));
    }
    static double logitLinkInv(double x) {
	double tmp = safeExp(x);
	return tmp/(1 + tmp);
    }
    static double logitLink(double x) {return log(x/(1 - x));}
    static double logitMuEta(double x) {
	double tmp = safeExp(x);
	double mu = tmp / (1 + tmp);
	return mu * (1 - mu);
    }

    static double probitLink(double x) {return Rf_qnorm5(x, 0, 1, 1, 0);}
    static double probitLinkInv(double x) {
	return std::max(std::min(Rf_pnorm5(x, 0, 1, 1, 0), 1. - epsilon), epsilon);
    }
    static double probitMuEta(double x) {
	return std::max(Rf_dnorm4(x, 0, 1, 0), epsilon);
    }
};

#endif /* LME4_GLMFAMILY_H */

