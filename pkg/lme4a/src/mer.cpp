#include "mer.h"

using namespace Rcpp;
using namespace MatrixNs;

extern cholmod_common c;

static double l2PI = log(2. * PI);

namespace mer{

    reModule::reModule(S4 xp) :
	L(S4(SEXP(xp.slot("L")))),
	Lambda(S4(SEXP(xp.slot("Lambda")))),
	Ut(S4(SEXP(xp.slot("Ut")))),
	Zt(S4(SEXP(xp.slot("Zt")))),
	Lind(SEXP(xp.slot("Lind"))),
	lower(SEXP(xp.slot("lower"))),
	theta(SEXP(xp.slot("theta"))),
	u(SEXP(xp.slot("u")))
    {
	NumericVector ldL2Vec(SEXP(xp.slot("ldL2")));
	ldL2 = ldL2Vec.begin();
    }

    void reModule::updateTheta(const NumericVector &nt) {
	if (nt.size() != theta.size())
	    ::Rf_error("length(theta) = %d != length(newtheta) = %d",
		       theta.size(), nt.size());
	std::copy(nt.begin(), nt.end(), theta.begin());
	double *Lamx = (double*)Lambda.x, *nnt = nt.begin();
	int *Li = Lind.begin();
	for (int i = 0; i < Lind.size(); i++) Lamx[i] = nnt[Li[i] - 1];

	CHM_SP LamTr = Lambda.transpose();
	CHM_SP LamTrZt = ::M_cholmod_ssmult(LamTr, &Zt, 0/*stype*/,
					    1/*values*/, 1/*sorted*/, &c);
	::M_cholmod_free_sparse(&LamTr, &c);
	Ut.update(*LamTrZt);
	::M_cholmod_free_sparse(&LamTrZt, &c);
	L.update(Ut, 1.);
	*ldL2 = ::M_chm_factor_ldetL2(&L);
    }

    void reModule::updateU(const merResp &resp) {
	chmDn cu(resp.cu);
	CHM_DN t1 = L.solve(CHOLMOD_Lt, &cu);
	CHM_DN t2 = L.solve(CHOLMOD_Pt, t1);
	::M_cholmod_free_dense(&t1, &c);
	double *t2b = (double*)t2->x;
	std::copy(t2b,  t2b + u.size(), u.begin());
	::M_cholmod_free_dense(&t2, &c);
    }

    void reModule::incGamma(NumericVector &gam) {
	chmDn cu(u), cgam(gam);
	Ut.dmult('T', 1., 1., cu, cgam);
    }

    //< Create RZX from ZtX or cu from Zty
    static void DupdateL(reModule &re, chmDn &src, chmDn &dest) {
	re.Lambda.dmult(1/*trans*/, 1., 0., src, dest);
	CHM_DN t1 = re.L.solve(CHOLMOD_P, &dest);
	CHM_DN t2 = re.L.solve(CHOLMOD_L, t1);
	::M_cholmod_free_dense(&t1, &c);

	double *t2b = (double*)t2->x, *db = (double*)(dest.x);
	int sz = src.nr() * src.nc();
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
	y(SEXP(xp.slot("y"))),
	cUtr(Utr),
	ccu(cu)
    {
	NumericVector wrssVec(SEXP(xp.slot("wrss")));
	wrss = wrssVec.begin();

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
	double *mm = mu.begin(), *rr = resid.begin(), *yy = y.begin();
	if (weights.size()) Rf_error("At present, weights must be numeric(0)");
	*wrss = 0;
	for (int i = 0; i < n; i++) {
	    rr[i] = yy[i] - mm[i];
	}
	*wrss = std::inner_product(resid.begin(), resid.end(), resid.begin(), double());
    }

    lmerDeFeMod::lmerDeFeMod(S4 xp) :
	deFeMod(xp),
	ZtX(S4(SEXP(xp.slot("ZtX")))),
	XtX(S4(SEXP(xp.slot("XtX")))),
	cZtX(ZtX)
    {
	NumericVector ldR2Vec(SEXP(xp.slot("ldR2")));
	ldR2 = ldR2Vec.begin();
    }
    
    void lmerDeFeMod::updateRzxRx(reModule &re) {
	DupdateL(re, cZtX, cRZX);
	RX.update('T', -1., RZX, 1., XtX);
	*ldR2 = RX.logDet2();
    }

    void lmerDeFeMod::updateBeta(merResp &resp) {
	std::copy(resp.Vtr.begin(), resp.Vtr.end(), resp.cbeta.begin());
	RZX.dgemv('T', -1., resp.cu, 1., resp.cbeta);
	std::copy(resp.cbeta.begin(), resp.cbeta.end(), beta.begin());
	RX.dpotrs(beta);
	RZX.dgemv('N', -1., beta, 1., resp.cu);
    }

    double lmer::deviance() {
	double nn = (double)resp.y.size(),
	    prss = re.sqLenU() + *resp.wrss;
	return *re.ldL2 + nn * (1 + l2PI + log(prss/nn));
    }

    double lmerDe::updateTheta(const NumericVector &nt) {
	re.updateTheta(nt);
	resp.updateL(re);
	fe.updateRzxRx(re);
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

