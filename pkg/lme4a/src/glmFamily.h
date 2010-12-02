// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
#ifndef LME4_GLMFAMILY_H
#define LME4_GLMFAMILY_H

#include <Rcpp.h>
#include <Rmath.h>

namespace glm {
    // Map (associative array) of functions returning double from a double
    typedef std::map<std::string, double(*)(double)> fmap;
    typedef std::map<std::string,
		     double(*)(double,double,double)> drmap;

    class glmFamily {
    protected:
	Rcpp::List           lst; // original list from R
	std::string d_family, d_link; // as in the R glm family
	Rcpp::Function d_devRes, d_linkfun, d_linkinv, d_muEta, d_variance;
    public:
	glmFamily(Rcpp::List) throw (std::runtime_error);

	const std::string& fam() const {return d_family;}
	const std::string& lnk() const {return d_link;}

	// initialize the associative arrays of scalar functions
	void initMaps();
	
	// Application of functions from the family
	// The scalar transformations use compiled code when available 
	Rcpp::NumericVector  linkFun(Rcpp::NumericVector const&) const;
	Rcpp::NumericVector  linkInv(Rcpp::NumericVector const&) const;
	Rcpp::NumericVector devResid(
	    Rcpp::NumericVector const&,
	    Rcpp::NumericVector const&,
	    Rcpp::NumericVector const&) const;
	Rcpp::NumericVector    muEta(Rcpp::NumericVector const&) const;
	Rcpp::NumericVector variance(Rcpp::NumericVector const&) const;
    private:
	// Class members that are maps storing the scalar functions
	static fmap
	    lnks,		// scalar link functions
	    linvs,		// scalar linkinv functions
	    muEtas,		// scalar muEta functions
	    varFuncs;		// scalar variance functions
	static drmap devRes;
	
	// Thresholds common to the class
	static double epsilon;

	// Utility for deviance residuals
	static inline double y_log_y(double y, double mu) {
	    return y ? y*log(y/mu) : 0;
	}

	// Scalar functions that will used in transform applications
	static inline double         cubef(double x) {return x * x * x;}
	static inline double        identf(double x) {return x;}
	static inline double     invderivf(double x) {return -1/(x * x);}
	static inline double      inversef(double x) {return 1/x;}
	static inline double          onef(double x) {return 1;}
	static inline double          sqrf(double x) {return x * x;}
	static inline double         twoxf(double x) {return 2 * x;}
	static inline double         x1mxf(double x) {
	    return std::max(epsilon, x * (1 - x));
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
	// cumulative probability function of the complement of the
	// Gumbel distribution (i.e. pgumbel(x,0.,1.,0) == 1-pgumbel2(-x,0.,1.,0))
	static inline double	
	pgumbel2(double q, double loc, double scale, int lower_tail) {
	    q = (q - loc) / scale;
	    q = -exp(q);
	    return lower_tail ? -expm1(q) : exp(q);
	}
	static inline double cloglogLinkInv(double x) {
	    return pgumbel2(x, 0., 1., 1);
	}

	//density of the complement of the Gumbel distribution
	static inline double
	dgumbel2(double x, double loc, double scale, int give_log) {
	    x = (x - loc) / scale;
	    x = x - exp(x) - log(scale);
	    return give_log ? x : exp(x);
	}
	static inline double   cloglogMuEta(double x) {
	    return dgumbel2(x, 0., 1., 0);
	}
	
	// deviance residuals functions
	static inline double
	BinomialDevRes(double y, double mu, double wt) {
	    return 2 * wt * (y_log_y(y, mu) + y_log_y(1 - y, 1 - mu));
	}
	static inline double logYMu(double y, double mu) {
	    return y ? log(y/mu) : 0;
	}
	static inline double
	GammaDevRes   (double y, double mu, double wt) {
	    return -2 * wt * (logYMu(y, mu) - (y - mu)/mu);
	}
	static inline double
	GaussianDevRes(double y, double mu, double wt) {
	    double res = y - mu;
	    return wt * res * res;
	}
	static inline double
	PoissonDevRes (double y, double mu, double wt) {
	    return 2 * wt * (y_log_y(y, mu) - (y - mu));
	}
    };
}
    
#endif /* LME4_GLMFAMILY_H */

