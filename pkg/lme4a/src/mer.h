// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
#ifndef LME4_MER_H
#define LME4_MER_H

#include <Rcpp.h>
#include "MatrixNs.h"

namespace mer {
    class merResp;		// forward declaration

    class reModule {
    public:
	reModule(Rcpp::S4);

	void updateTheta(const Rcpp::NumericVector&);
	double sqLenU() const {	// squared length of u
	    return std::inner_product(u.begin(), u.end(), u.begin(), double());
	}
	void updateU(const merResp&);
	void incGamma(Rcpp::NumericVector&);

	MatrixNs::chmFa L;
	MatrixNs::chmSp Lambda, Ut, Zt;
	Rcpp::IntegerVector Lind;
	Rcpp::NumericVector lower, theta, u, ldL2;
    };
    
    class merResp {
    public:
	merResp(Rcpp::S4);
	void updateL(reModule&);
	void updateWrss();	// resid := y - mu; wrss := sum((sqrtrwts * resid)^2)
	
	Rcpp::NumericVector Utr, Vtr, cbeta,
	    cu, mu, offset, resid, weights, wrss, y; 
	MatrixNs::chmDn cUtr, ccu;
    };

    class feModule {
    public:
	Rcpp::NumericVector beta;
	feModule(Rcpp::S4 xp) : beta(SEXP(xp.slot("beta"))) { }
    };
    
    class deFeMod : public feModule {
    public:
	deFeMod(Rcpp::S4 xp) :
	    feModule(xp),
	    X(Rcpp::S4(SEXP(xp.slot("X")))),
	    RZX(Rcpp::S4(SEXP(xp.slot("RZX")))),
	    RX(Rcpp::S4(SEXP(xp.slot("RX")))),
	    cRZX(RZX) { }
	
	MatrixNs::dgeMatrix X, RZX;
	MatrixNs::Cholesky RX;
	MatrixNs::chmDn cRZX;
    };

    class lmerDeFeMod : public deFeMod {
    public:
	lmerDeFeMod(Rcpp::S4);
	void updateL(reModule&);
	void updateBeta(merResp&);
	void incGamma(Rcpp::NumericVector &gam) {
	    X.dgemv('N', 1., beta, 1., gam);
	}
	
	MatrixNs::dgeMatrix ZtX;
	MatrixNs::dpoMatrix XtX;
	Rcpp::NumericVector ldR2;
	MatrixNs::chmDn cZtX;
    };

    class lmer {
    public:
	lmer(Rcpp::S4 xp) :
	    re(Rcpp::S4(SEXP(xp.slot("re")))),
	    resp(Rcpp::S4(SEXP(xp.slot("resp")))),
	    REML(SEXP(xp.slot("REML"))) {
	    reml = (bool)*REML.begin();
	}
    
	double deviance();

	reModule re;
	merResp resp;
	Rcpp::LogicalVector REML;
	bool reml;
    };
	
    class lmerDe : public lmer {
    public:
	lmerDe(Rcpp::S4 xp) :
	    lmer(xp),
	    fe(Rcpp::S4(SEXP(xp.slot("fe")))) {
	}
	double REMLcrit();
	double updateTheta(const Rcpp::NumericVector&);
	lmerDeFeMod fe;
    };
}

#endif

