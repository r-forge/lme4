#include "mer.h"
#include <cstring>

using namespace Rcpp;

extern cholmod_common c;

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

    reModule::reModule(S4 xp) :
	L(S4(SEXP(xp.slot("L")))),
	Lambda(S4(SEXP(xp.slot("Lambda")))),
	Ut(S4(SEXP(xp.slot("Ut")))),
	Zt(S4(SEXP(xp.slot("Zt")))),
	Lind(SEXP(xp.slot("Lind"))),
	lower(SEXP(xp.slot("lower"))),
	theta(SEXP(xp.slot("theta"))),
	u(SEXP(xp.slot("u"))),
	ldL2(SEXP(xp.slot("ldL2"))) {
    }

    void reModule::updateTheta(NumericVector nt) {
	if (nt.size() != theta.size())
	    ::Rf_error("length(theta) = %d != length(newtheta) = %d",
		       theta.size(), nt.size());
	std::copy(nt.begin(), nt.end(), theta.begin());
	double *LamxP = (Lambda.x).begin();
	for (int *Li = Lind.begin(); Li < Lind.end(); Li++)
	    *LamxP++ = nt[(*Li) - 1];

	CHM_SP LamTr = ::M_cholmod_transpose(Lambda.sp, 1/*values*/, &c);
	CHM_SP LamTrZt = ::M_cholmod_ssmult(LamTr, Zt.sp, 0/*stype*/,
					    1/*values*/, 1/*sorted*/, &c);
	::M_cholmod_free_sparse(&LamTr, &c);
	Ut.update(LamTrZt);
	::M_cholmod_free_sparse(&LamTrZt, &c);
	L.update(Ut.sp, 1.);
	*(ldL2.begin()) = ::M_chm_factor_ldetL2(L.fa);
    }
}

extern "C"
SEXP check_resp(SEXP xp) {
    S4 x4(xp);
    mer::merResp rr(x4);
    return Rf_ScalarLogical(1);
}

extern "C"
SEXP update_reModule(SEXP xp, SEXP ntheta) {
    S4 x4(xp);
    mer::reModule rr(x4);
    NumericVector nt(ntheta);
    rr.updateTheta(nt);
    return R_NilValue;
}

