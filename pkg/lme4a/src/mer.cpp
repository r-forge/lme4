#include "mer.h"

using namespace Rcpp;
using namespace MatrixNs;

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
	y(SEXP(xp.slot("y"))),
	ccu(cu),
	cUtr(Utr)
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

    // update dense src and dest objects with Lambda and L
    static void DupdateL(reModule &re, chmDn &src, chmDn &dest) {
	double one[] = {1, 0}, zero[] = {0, 0};
	Rprintf("update1: L of size %d\n", (re.L.fa)->n);
	Rprintf("Lambda(%d,%d), src(%d,%d), dest(%d,%d)\n",
		re.Lambda.nrow, re.Lambda.ncol,
		src.nr(), src.nc(), dest.nr(), dest.nc());
	int sz = src.nr() * src.nc();

	::M_cholmod_sdmult(&(re.Lambda), 1/*trans*/, one, zero, &src, &dest, &c);
	CHM_DN t1 = ::M_cholmod_solve(CHOLMOD_P, re.L.fa, &dest, &c);
	CHM_DN t2 = ::M_cholmod_solve(CHOLMOD_L, re.L.fa, t1, &c);
	::M_cholmod_free_dense(&t1, &c);
	double *t2b = (double*)t2->x, *db = (double*)(dest.x);
	std::copy(t2b, t2b + sz, db);
	::M_cholmod_free_dense(&t2, &c);
    }

    void merResp::updateL(reModule &re) {
	DupdateL(re, cUtr, ccu);
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
	double *Lamx = (double*)Lambda.x, *nnt = nt.begin();
	int *Li = Lind.begin();
	for (int i = 0; i < Lind.size(); i++) Lamx[i] = nnt[Li[i] - 1];

	CHM_SP LamTr = ::M_cholmod_transpose(&Lambda, 1/*values*/, &c);
	CHM_SP LamTrZt = ::M_cholmod_ssmult(LamTr, &Zt, 0/*stype*/,
					    1/*values*/, 1/*sorted*/, &c);
	::M_cholmod_free_sparse(&LamTr, &c);
	Ut.update(LamTrZt);
	::M_cholmod_free_sparse(&LamTrZt, &c);
	L.update(&Ut, 1.);
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
	RX(S4(SEXP(xp.slot("RX")))),
	cRZX(RZX) {
    }

    lmerDeFeMod::lmerDeFeMod(S4 xp) :
	deFeMod(xp),
	ZtX(S4(SEXP(xp.slot("ZtX")))),
	XtX(S4(SEXP(xp.slot("XtX")))),
	ldR2(SEXP(xp.slot("ldR2"))),
	cZtX(ZtX) {
    }
    
    void lmerDeFeMod::updateL(reModule &re) {
	DupdateL(re, cZtX, cRZX);
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
	double prss = re.sqLenU() + *(resp.wrss.begin());
	return *(re.ldL2.begin()) + log(M_2PI * prss) + 1./((double)resp.y.size());
    }

    double lmerDe::updateTheta(const NumericVector &nt) {
	Rprintf("begin re.updateTheta\n");
	re.updateTheta(nt);
	Rprintf("begin resp.updateL\n");
	resp.updateL(re);
	Rprintf("begin fe.updateL\n");
	fe.updateL(re);
	return 0.;
    }
}

extern "C"
SEXP update_lmer2(SEXP xp, SEXP ntheta) {
    S4 x4(xp);
    NumericVector nt(ntheta);
    mer::lmerDe lm(x4);

    return Rf_ScalarReal(lm.updateTheta(nt));
}

