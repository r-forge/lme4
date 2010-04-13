// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
#ifndef LME4_MER_H
#define LME4_MER_H

#include <Rcpp.h>
#include "Matrix_ns.h"

namespace mer {
    
    class merResp {
    public:
	merResp(Rcpp::S4);
	
	Rcpp::NumericVector Utr, Vtr, cbeta,
	    cu, mu, offset, resid, weights, wrss, y; 
    };

    class reModule {
    public:
	reModule(Rcpp::S4);
	void updateTheta(Rcpp::NumericVector);

	Matrix::dCHMfactor L;
	Matrix::dgCMatrix Lambda, Ut, Zt;
	Rcpp::IntegerVector Lind;
	Rcpp::NumericVector lower, theta, u, ldL2;
    };

    class feModule {
    public:
	feModule(Rcpp::S4);
	
	Rcpp::NumericVector beta;
    };
    
    class deFeMod : public feModule {
    public:
	deFeMod(Rcpp::S4);
	
	Matrix::dgeMatrix X, RZX;
	Matrix::Cholesky RX;
    };

    class lmerDeFeMod : public deFeMod {
    public:
	lmerDeFeMod(Rcpp::S4);
	
	Matrix::dgeMatrix ZtX;
//	Matrix::dpoMatrix XtX;
	Rcpp::NumericVector ldR2;
    };

}

#endif

