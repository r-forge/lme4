// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
#ifndef LME4_MER_H
#define LME4_MER_H

#include <Rcpp.h>
#include "MatrixNs.h"

namespace mer {
    class reModule {
    public:
	reModule(Rcpp::S4);
	void updateTheta(const Rcpp::NumericVector&);
	double sqLenU();	// squared length of u

	MatrixNs::dCHMfactor L;
	MatrixNs::chmSp Lambda, Ut, Zt;
	Rcpp::IntegerVector Lind;
	Rcpp::NumericVector lower, theta, u, ldL2;
    };
    
    class merResp {
    public:
	merResp(Rcpp::S4);
	void updateL(reModule&);  // FIXME: see comments in mer.cpp
	
	Rcpp::NumericVector Utr, Vtr, cbeta,
	    cu, mu, offset, resid, weights, wrss, y; 
	MatrixNs::chmDn cUtr, ccu;
    };

    class feModule {
    public:
	feModule(Rcpp::S4);
	
	Rcpp::NumericVector beta;
    };
    
    class deFeMod : public feModule {
    public:
	deFeMod(Rcpp::S4);
	
	MatrixNs::dgeMatrix X, RZX;
	MatrixNs::Cholesky RX;
	MatrixNs::chmDn cRZX;
    };

    class lmerDeFeMod : public deFeMod {
    public:
	lmerDeFeMod(Rcpp::S4);
	void updateL(reModule&);
	void updateBeta(merResp&);
	
	MatrixNs::dgeMatrix ZtX;
	MatrixNs::dpoMatrix XtX;
	Rcpp::NumericVector ldR2;
	MatrixNs::chmDn cZtX;
    };

    class lmer {
    public:
	lmer(Rcpp::S4);
	double deviance();

	reModule re;
	merResp resp;
	Rcpp::LogicalVector REML;
	bool reml;
    };
	
    class lmerDe : public lmer {
    public:
	lmerDe(Rcpp::S4);
	double REMLcrit();
	double updateTheta(const Rcpp::NumericVector&);

	lmerDeFeMod fe;
    };
}

#endif

