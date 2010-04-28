#include "mer.h"

using namespace Rcpp;
using namespace MatrixNs;

static double l2PI = log(2. * PI);

static void showdbl(const double* x, const char* nm, int n) {
    int n5 = (n < 5) ? n : 5;
    Rprintf("%s: %g", nm, x[0]);
    for (int i = 1; i < n5; i++) Rprintf(", %g", x[i]);
    Rprintf("\n");
}

static void showint(const int* x, const char* nm, int n) {
    int n20 = (n < 20) ? n : 20;
    Rprintf("%s: %d", nm, x[0]);
    for (int i = 1; i < n20; i++) Rprintf(", %d", x[i]);
    Rprintf("\n");
}

static void showCHM_SP(const CHM_SP x, const std::string nm) {
    Rprintf("%s: nrow = %d, ncol = %d, xtype = %d, stype = %d, itype = %d, dtype = %d\n",
	    nm.c_str(), x->nrow, x->ncol, x->xtype, x->stype, x->itype, x->dtype);
    int nc = x->ncol, nnz = ::M_cholmod_nnz(x, &c);
    showint((int*)x->p, "p", nc + 1);
    showint((int*)x->i, "i", nnz);
    showdbl((double*)x->x, "x", nnz);
}

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

	CHM_SP LamTrZt = Lambda.crossprod(Zt);
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
	re.Lambda.dmult('T', 1., 0., src, dest);
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
    
    void lmerDeFeMod::updateRzxRx(reModule &re) {
	chmDn cZtX(ZtX), cRZX(RZX);
	DupdateL(re, cZtX, cRZX);
	RX.update('T', -1., RZX, 1., XtX);
	*ldRX2 = RX.logDet2();
    }

    void lmerDeFeMod::updateBeta(merResp &resp) {
	std::copy(resp.Vtr.begin(), resp.Vtr.end(), resp.cbeta.begin());
	RZX.dgemv('T', -1., resp.cu, 1., resp.cbeta);
	std::copy(resp.cbeta.begin(), resp.cbeta.end(), beta.begin());
	RX.dpotrs(beta);
	RZX.dgemv('N', -1., beta, 1., resp.cu);
    }

    void lmerSpFeMod::updateRzxRx(reModule &re) {
	double mone[] = {-1.,0}, one[] = {1.,0};
//	showCHM_SP(&re.Lambda, "Lambda");
	CHM_SP t1 = re.Lambda.transpose();
//	showCHM_SP(t1, "t(Lambda)");
	CHM_SP t2 = ::M_cholmod_ssmult(t1, &ZtX, 0/*stype*/,
				       1/*values*/, 1/*sorted*/, &c);
//	showCHM_SP(t2, "crossprod(Lambda, ZtX)");
//	Rprintf("L.n = %d\n", re.L.n);
//	showint((int*)re.L.Perm, "re.L.Perm", re.L.n);
	::M_cholmod_free_sparse(&t1, &c);
	t1 = re.L.spsolve(CHOLMOD_P, t2);
//	showCHM_SP(t1, "solve(L, crossprod(Lambda, ZtX), \"P\"");

//	showint((int*)t1->p, "t1->p", (t1->ncol) + 1);
//	Rprintf("solve(crossprod(Lambda, ZtX), L, \"P\") (%d, %d), nnz = %d\n",
//		t1->nrow, t1->ncol, ::M_cholmod_nnz(t1, &c));
	::M_cholmod_free_sparse(&t2, &c);
	t2 = re.L.spsolve(CHOLMOD_L, t1);
//	Rprintf("solve(crossprod(Lambda, ZtX), L, \"L\") (%d, %d), nnz = %d\n",
//		t2->nrow, t2->ncol, ::M_cholmod_nnz(t2, &c));
	RZX.update(*t2);

	::M_cholmod_free_sparse(&t1, &c);
	t1 = RZX.transpose();
	::M_cholmod_free_sparse(&t2, &c);
	t2 = ::M_cholmod_aat(t1, (int*)NULL, 0/*fsize*/, 1/*mode*/, &c);
	::M_cholmod_free_sparse(&t1, &c);
// Need to copy and convert to a symmetric matrix; grumble, grumble.
	t1 = ::M_cholmod_copy(t2, 1/*stype*/, 1/*mode*/, &c);
	::M_cholmod_free_sparse(&t2, &c);

	t2 = ::M_cholmod_add(&XtX, t1, one, mone, 1/*values*/, 1/*sorted*/, &c);
//	showCHM_SP(t1, "t1");
	::M_cholmod_free_sparse(&t1, &c);
	RX.update(*t2, 0.);
	::M_cholmod_free_sparse(&t2, &c);
	*ldRX2 = ::M_chm_factor_ldetL2(&RX);
    }

    void lmerSpFeMod::updateBeta(merResp &resp) {
	std::copy(resp.Vtr.begin(), resp.Vtr.end(), resp.cbeta.begin());
	chmDn ccbeta(resp.cbeta);
	RZX.dmult('T', -1., 1., resp.ccu, ccbeta);
	std::copy(resp.cbeta.begin(), resp.cbeta.end(), beta.begin());
	chmDn cbeta(beta);
	CHM_DN t1 = RX.solve(CHOLMOD_A, &cbeta);
	double *t1b = (double*)t1->x;
	std::copy(t1b,  t1b + beta.size(), beta.begin());
	RZX.dmult('N', -1., 1., cbeta, resp.ccu);
    }

    double lmer::deviance() {
	double nn = (double)resp.y.size(),
	    prss = re.sqLenU() + *resp.wrss;
	return *re.ldL2 + nn * (1 + l2PI + log(prss/nn));
    }

    double lmerSp::reCrit() {
	double nmp = (double)(resp.y.size() - fe.beta.size()),
	    prss = re.sqLenU() + *resp.wrss;
	return *re.ldL2 + *fe.ldRX2 + nmp * (1 + l2PI + log(prss/nmp));
    }

    double lmerDe::reCrit() {
	double nmp = (double)(resp.y.size() - fe.beta.size()),
	    prss = re.sqLenU() + *resp.wrss;
	return *re.ldL2 + *fe.ldRX2 + nmp * (1 + l2PI + log(prss/nmp));
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
	return reml ? reCrit() : deviance();
    }

//FIXME:: factor out the parts that do not depend on the fe class (or make those virtual functions?)
    double lmerSp::updateTheta(const NumericVector &nt) {
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
	return reml ? reCrit() : deviance();
    }
}

extern "C"
SEXP update_lmerDe(SEXP xp, SEXP ntheta) {
    S4 x4(xp);
    NumericVector nt(ntheta);
    mer::lmerDe lm(x4);

    return Rf_ScalarReal(lm.updateTheta(nt));
}

extern "C"
SEXP update_lmerSp(SEXP xp, SEXP ntheta) {
    S4 x4(xp);
    NumericVector nt(ntheta);
    mer::lmerSp lm(x4);

    return Rf_ScalarReal(lm.updateTheta(nt));
}

extern "C"
SEXP deviance_lmerDe(SEXP xp) {
    S4 x4(xp);
    mer::lmerDe lm(x4);
    return Rf_ScalarReal(lm.deviance());
}

extern "C"
SEXP deviance_lmerSp(SEXP xp) {
    S4 x4(xp);
    mer::lmerSp lm(x4);
    return Rf_ScalarReal(lm.deviance());
}

