#ifndef LME4_GLMFAMILY_H
#define LME4_GLMFAMILY_H

#include <Rcpp.h>

typedef std::map<std::string, double(*)(double)> fmap;

/// generalized linear model family
class glmFamily {
    static fmap
	lnks,			//< scalar link functions
	linvs,			//< scalar linkinv functions
	muEtas,			//< scalar muEta functions
	varFuncs;		//< scalar variance functions
    Rcpp::List lst;
public:
    glmFamily(Rcpp::List);
    void linkFun(Rcpp::NumericVector eta, const Rcpp::NumericVector mu);
    void linkInv(Rcpp::NumericVector mu, const Rcpp::NumericVector eta);
    void muEta(Rcpp::NumericVector mueta, const Rcpp::NumericVector eta);
    void variance(Rcpp::NumericVector vv, const Rcpp::NumericVector eta);
    double devResid(const Rcpp::NumericVector mu,
		    const Rcpp::NumericVector weights,
		    const Rcpp::NumericVector y);
    SEXP show();
    
    std::string family, link;
};

#endif /* LME4_GLMFAMILY_H */

