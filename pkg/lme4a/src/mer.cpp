#include "mer.h"

using namespace Rcpp;

namespace mer{
    
    merResp::merResp(S4 xp) :
	Utr(SEXP(xp.slot("Utr"))),
	Vtr(SEXP(xp.slot("Vtr"))),
	cbeta(SEXP(xp.slot("cbeta"))),
	cu(SEXP(xp.slot("cu"))),
	mu(SEXP(xp.slot("mu"))),
	offset(SEXP(xp.slot("offset"))),
	resid(SEXP(xp.slot("resid"))),
	weights(SEXP(xp.slot("weights"))),
	wrss(SEXP(xp.slot("wrss"))),
	y(SEXP(xp.slot("y")))
    {
	int n = y.size(), ws = weights.size();
	
	if (mu.size() != n || resid.size() != n)
	    ::Rf_error("y, mu and resid slots must have equal lengths");
	if (cbeta.size() != Vtr.size())
	    ::Rf_error("cbeta and Vtr slots must have equal lengths");
	if (cu.size() != Utr.size())
	    ::Rf_error("cu and Utr slots must have equal lengths");
	if (ws && ws != n)
	    ::Rf_error("weights slot must have length 0 or n");
    }

    dgCMatrix::dgCMatrix(S4 xp) :
	i(SEXP(xp.slot("i"))),
	p(SEXP(xp.slot("p"))),
	Dim(SEXP(xp.slot("Dim"))),
	Dimnames(SEXP(xp.slot("Dim"))),
	factors(SEXP(xp.slot("Dim"))),
	x(SEXP(xp.slot("x")))
    {
	int nnz = p[Dim[1]];
	if (i.size() != nnz || x.size() != nnz)
	    Rf_error("size of i and x must match p[Dim[2] + 1]");
    }

    SEXP dgCMatrix::dims() {
	return wrap(clone(Dim));
    }
    
    reModule::reModule(S4 xp) :
	L(S4(SEXP(xp.slot("L")))),
	Lambda(S4(SEXP(xp.slot("Lambda")))),
	Ut(S4(SEXP(xp.slot("Ut")))),
	Zt(S4(SEXP(xp.slot("Zt")))),
	Lind(SEXP(xp.slot("Lind"))),
	lower(SEXP(xp.slot("lower"))),
	theta(SEXP(xp.slot("theta"))),
	u(SEXP(xp.slot("u"))),
	ldL2(SEXP(xp.slot("ldL2")))
    {
    }
}

extern "C"
SEXP check_dgCMatrix(SEXP xp) {
    mer::dgCMatrix rr(S4(xp));
    return Rf_ScalarLogical(1);
// Somehow rr is not what I think it should be here.
//    return rr.dims();
}

extern "C"
SEXP check_resp(SEXP xp) {
    mer::merResp rr(S4(xp));
    return Rf_ScalarLogical(1);
}

extern "C"
SEXP check_reModule(SEXP xp) {
    mer::reModule rr(S4(xp));
    return Rf_ScalarLogical(1);
}
