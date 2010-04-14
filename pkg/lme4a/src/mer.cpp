#include "mer.h"

using namespace Rcpp;
using namespace Matrix;

extern cholmod_common c;

static double M_2PI = 6.283185307179586476925286766559;

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

    static void updateL(reModule &re, chmDn &src, chmDn &dest) {
	double one[] = {1, 0}, zero[] = {0, 0};
	int sz = (src.dn)->nrow * (src.dn)->ncol;

	M_cholmod_sdmult(re.Lambda.sp, 1/*trans*/, one, zero, src.dn, dest.dn, &c);
	CHM_DN t1 = M_cholmod_solve(CHOLMOD_P, re.L.fa, desy.dn, &c);
	CHM_DN t2 = M_cholmod_solve(CHOLMOD_L, re.L.fa, t1, &c);
	M_cholmod_free_dense(&t1, &c);
	std::copy((double*)t2->x, (double*)((t2->x) + .size()), cu.begin());
	M_cholmod_free_dense(&t2, &c);
    }

    void merResp::updateL(reModule &re) {
	chmDn uu(Utr), cuu(cu);
	updateL(re, uu, cuu);
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

    void reModule::updateTheta(const NumericVector &nt) {
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

    double reModule::sqLenU() {
	return std::inner_product(u.begin(), u.end(), u.begin(), double());
    }

    feModule::feModule(S4 xp) : beta(SEXP(xp.slot("beta"))) {
    }

    deFeMod::deFeMod(S4 xp) :
	feModule(xp),
	X(S4(SEXP(xp.slot("X")))),
	RZX(S4(SEXP(xp.slot("RZX")))),
	RX(S4(SEXP(xp.slot("RX")))) {
    }

    lmerDeFeMod::lmerDeFeMod(S4 xp) :
	deFeMod(xp),
	ZtX(S4(SEXP(xp.slot("ZtX")))),
	XtX(S4(SEXP(xp.slot("XtX")))),
	ldR2(SEXP(xp.slot("ldR2"))) {
    }
    
    void lmerDeFeMod::updateL(reModule &re) {
	chmDn ztx(ZtX), rzx(RZX);
	updateL(re, ztx, rzx);
    }

    lmer::lmer(S4 xp) :
	re(S4(SEXP(xp.slot("re")))),
	resp(S4(SEXP(xp.slot("resp")))),
	REML(SEXP(xp.slot("REML"))) {
	reml = (bool)*REML.begin();
    }
    
    lmerDe::lmerDe(S4 xp) :
	lmer(xp),
	fe(S4(SEXP(xp.slot("fe")))) {
    }

    double lmer::deviance() {
	return re.ldL2 + log(M_2PI * (re.sqLenU() + resp.wrss)) +
	    1./((double)resp.y.size());
    }

    double lmerDe::updateTheta(const NumericVector &nt) {
	re.updateTheta(nt);
	fe.updateL(re);
	resp.updateL(re);
	return 0.;
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

