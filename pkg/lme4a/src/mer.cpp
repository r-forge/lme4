#include "mer.h"

using namespace Rcpp;
using namespace MatrixNs;

extern cholmod_common c;

static double M_2PI = 6.283185307179586476925286766559;

namespace mer{

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
	Ut.update(*LamTrZt);
	::M_cholmod_free_sparse(&LamTrZt, &c);
	L.update(Ut, 1.);
	*(ldL2.begin()) = ::M_chm_factor_ldetL2(&L);
    }

    void reModule::updateU(const merResp &resp) {
	std::copy(resp.cu.begin(), resp.cu.end(), u.begin());
	chmDn cu(u);
	::M_cholmod_solve(CHOLMOD_Lt, &L, &cu, &c);
	::M_cholmod_solve(CHOLMOD_Pt, &L, &cu, &c);
    }

    void reModule::incGamma(NumericVector &gam) {
	double one[] = {1, 0}, zero[] = {0, 0};
	chmDn cu(u), gg(gam);
	::M_cholmod_sdmult(&Ut, 1/*trans*/, one, one, &cu, &gg, &c);
    }

    // Generic update of RZX from ZtX or cu from Zty using an reModule
    static void DupdateL(reModule &re, chmDn &src, chmDn &dest) {
	double one[] = {1, 0}, zero[] = {0, 0};
	int sz = src.nr() * src.nc();

	::M_cholmod_sdmult(&(re.Lambda), 1/*trans*/, one, zero, &src, &dest, &c);
	CHM_DN t1 = ::M_cholmod_solve(CHOLMOD_P, &(re.L), &dest, &c);
	CHM_DN t2 = ::M_cholmod_solve(CHOLMOD_L, &(re.L), t1, &c);
	::M_cholmod_free_dense(&t1, &c);
	double *t2b = (double*)t2->x, *db = (double*)(dest.x);
	std::copy(t2b, t2b + sz, db);
	::M_cholmod_free_dense(&t2, &c);
    }

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

    void merResp::updateL(reModule &re) {
	DupdateL(re, cUtr, ccu);
    }

    void merResp::updateWrss() {
	int n = y.size();
	double *mm = mu.begin(), *rr = resid.begin(), *yy = y.begin(),
	    rss = 0;
	if (weights.size()) Rf_error("At present, weights must be numeric(0)");
	for (int i = 0; i < n; i++) {
	    rr[i] = yy[i] - mm[i];
	    rss += rr[i] * rr[i];
	}
	*(wrss.begin()) = rss;
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
	RX.update('T', -1., RZX, 1., XtX);
// calculate ldR2 here
    }

    void lmerDeFeMod::updateBeta(merResp &resp) {
	std::copy(resp.Vtr.begin(), resp.Vtr.end(), beta.begin());
	RZX.dgemv('T', -1., resp.cu, 1., resp.cbeta);
	std::copy(resp.cbeta.begin(), resp.cbeta.end(), beta.begin());
	RX.dpotrs(beta);
	RZX.dgemv('N', -1., beta, 1., resp.cu);
    }

    double lmer::deviance() {
	double prss = re.sqLenU() + *(resp.wrss.begin());
	return *(re.ldL2.begin()) + log(M_2PI * prss) + 1./((double)resp.y.size());
    }

    double lmerDe::updateTheta(const NumericVector &nt) {
	re.updateTheta(nt);
	resp.updateL(re);
	fe.updateL(re);
	fe.updateBeta(resp);
	re.updateU(resp);
				// update resp.mu
	if (resp.offset.size())	// FIXME: define offset to be rep(0, n) when numeric(0)
	    std::copy(resp.offset.begin(), resp.offset.end(), resp.mu.begin());
	else std::fill(resp.mu.begin(), resp.mu.end(), double());
	re.incGamma(resp.mu);
	fe.incGamma(resp.mu);
				// update resp.resid and resp.wrss
	resp.updateWrss();
	return deviance();
    }
}

extern "C"
SEXP update_lmer2(SEXP xp, SEXP ntheta) {
    Rcpp::S4 x4(xp);
    Rcpp::NumericVector nt(ntheta);
    mer::lmerDe lm(x4);

    return Rf_ScalarReal(lm.updateTheta(nt));
}

