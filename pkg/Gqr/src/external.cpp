// external.cpp: externally .Call'able functions in GaussQuad package
//
// Copyright (C)       2011 Douglas Bates, Martin Maechler and Ben Bolker
//
// This file is part of the GaussQuad package

#include <Rcpp.h>
#include "GQ.h"

extern "C" {
    using GaussQuad::GQ;

    SEXP Gauss_Quad(SEXP n, SEXP rule, SEXP a, SEXP b, SEXP alpha, SEXP beta) {
	BEGIN_RCPP;
	GaussQuad::Rule rr=GaussQuad::Legendre; // -Wall
	switch(::Rf_asInteger(rule)) {
	case 1: rr = GaussQuad::Legendre;  break;
	case 2: rr = GaussQuad::Chebyshev; break;
	case 3: rr = GaussQuad::Gegenbauer; break;
	case 4: rr = GaussQuad::Jacobi; break;
	case 5: rr = GaussQuad::Laguerre; break;
	case 6: rr = GaussQuad::Hermite; break;
	case 7: rr = GaussQuad::Exponential; break;
	case 8: rr = GaussQuad::Rational; break;
	case 9: rr = GaussQuad::Type2; break;
	default: throw std::invalid_argument("Unknown rule");
	}
	GQ gq(::Rf_asInteger(n),  rr, ::Rf_asReal(a),
	      ::Rf_asReal(b), ::Rf_asReal(alpha), ::Rf_asReal(beta));
	return Rcpp::List::create(Rcpp::Named("knots")   = gq.t(),
				  Rcpp::Named("weights") = gq.wts()
//				  , Rcpp::Named("diag")    = gq.diag()
//				  , Rcpp::Named("sub")     = gq.sub()
//				  , Rcpp::Named("zemu")    = ::Rf_ScalarReal(gq.zemu())
	    );
	END_RCPP;
    }
}

#include <R_ext/Rdynload.h>

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

static R_CallMethodDef CallEntries[] = {
    CALLDEF(Gauss_Quad, 6),
    {NULL, NULL, 0}
};

/** Initializer for GaussQuad, called upon loading the package.
 *
 *  Register routines that can be called directly from R.
 */
extern "C"
void R_init_Gqr(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, (Rboolean)FALSE);
}

