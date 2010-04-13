// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
#ifndef LME4_MER_H
#define LME4_MER_H

#include <Rcpp.h>
#include <Matrix.h>

namespace mer {
    
    class merResp {
    public:
	merResp(Rcpp::S4);
	
	Rcpp::NumericVector Utr, Vtr, cbeta,
	    cu, mu, offset, resid, weights, wrss, y; 
    };

    class dgCMatrix {
    public:
	dgCMatrix(Rcpp::S4);
	~dgCMatrix(){delete sp;}
	SEXP dims() const;
	void update(CHM_SP);
	
	Rcpp::IntegerVector i, p, Dim;
	Rcpp::List Dimnames, factors;
	Rcpp::NumericVector x;
	CHM_SP sp;
    };
    
    class reModule {
    public:
	reModule(Rcpp::S4);
	void updateTheta(Rcpp::NumericVector);

	Rcpp::S4 L;
	dgCMatrix Lambda, Ut, Zt;
	Rcpp::IntegerVector Lind;
	Rcpp::NumericVector lower, theta, u, ldL2;
    };

}

#endif

