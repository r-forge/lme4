// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
#ifndef LME4_MER_H
#define LME4_MER_H

#include <Rcpp.h>

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
	SEXP dims();
	
	Rcpp::IntegerVector i, p, Dim;
	Rcpp::List Dimnames, factors;
	Rcpp::NumericVector x;
    };
    
    class reModule {
    public:
	reModule(Rcpp::S4);

	Rcpp::S4 L;
	dgCMatrix Lambda, Ut, Zt;
	Rcpp::IntegerVector Lind;
	Rcpp::NumericVector lower, theta, u, ldL2;
    };

}

#endif

