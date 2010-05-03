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
    if (nnz > 0) {
	showint((int*)x->i, "i", nnz);
	showdbl((double*)x->x, "x", nnz);
    } else Rprintf("nnz = 0\n");
}

static void showCHM_FR(const CHM_FR x, const std::string nm) {
    Rprintf("%s: n = %d, minor = %d, xtype = %d, itype = %d, dtype = %d\n",
	    nm.c_str(), x->n, x->minor, x->xtype, x->itype, x->dtype);
    Rprintf("ordering = %d, is_ll = %d, is_super = %d, is_monotonic = %d, nzmax = %d\n",
	    x->ordering, x->is_ll, x->is_super, x->is_monotonic, x->nzmax);
    int *pp = (int*)x->p;
    showint(pp, "p", (x->n) + 1);
    showint((int*)x->Perm, "Perm", x->n);
    if (pp[x->n] > 0) {
	showint((int*)x->i, "i", x->n);
	showdbl((double*)x->x, "x", x->n);
    } else Rprintf("nnz = 0\n");
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

    void rwReMod::incGamma(NumericVector &gam) {
	std::vector<double> uu(u.size());
	std::transform(u.begin(), u.end(), ubase.begin(), uu.begin(), std::plus<double>());
	chmDn cuu(uu), cgam(gam);
	Ut.dmult('T', 1., 1., cuu, cgam);
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

	int n = y.size(), os = offset.size();

	if (mu.size() != n || resid.size() != n ||
	    weights.size() != n || offset.size() != n)
	    ::Rf_error("y, mu, resid and weights slots must have equal lengths");
	if (os < 1 || os % n)
	    ::Rf_error("length(offset) must be a positive multiple of length(y)");
	if (cbeta.size() != Vtr.size())
	    ::Rf_error("cbeta and Vtr slots must have equal lengths");
	if (cu.size() != Utr.size())
	    ::Rf_error("cu and Utr slots must have equal lengths");
    }

    void merResp::updateL(reModule &re) {
	DupdateL(re, cUtr, ccu);
    }

    double merResp::updateWrss() {
	int n = y.size();
	double *mm = mu.begin(),
	    *rr = resid.begin(),
	    *ww = weights.begin(),
	    *yy = y.begin(),
	    ans = 0;
	for (int i = 0; i < n; i++) {
	    rr[i] = yy[i] - mm[i];
	    ans += rr[i] * rr[i] * ww[i];
	}
	*wrss = ans;
	return ans;
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

    void rwDeFeMod::incGamma(NumericVector &gam) {
// Should use std::vector<double> but need to add appropriate dgeMatrix::dgemv method
// It is not clear to me what the best design for that is.  I could
// blow up the number of possible combinations or I could coerce
// everything to a base type.  I don't know if as< std::vector<double> >
// is sufficiently lightweight to use std::vector<double> as the base type.
	NumericVector bb(beta.size()); 
	std::transform(beta.begin(), beta.end(), betabase.begin(),
		       bb.begin(), std::plus<double>());
	X.dgemv('N', 1., bb, 1., gam);
    }
    
    void lmerSpFeMod::updateRzxRx(reModule &re) {
	double mone[] = {-1.,0}, one[] = {1.,0};

	CHM_SP t2 = re.Lambda.crossprod(ZtX);
//	showCHM_FR(&(re.L), "L");
	CHM_SP t1 = re.L.spsolve(CHOLMOD_P, t2);
//	showCHM_SP(t1, "solve(L, crossprod(Lambda, ZtX), \"P\")");
	::M_cholmod_free_sparse(&t2, &c);
	t2 = re.L.spsolve(CHOLMOD_L, t1);
//	showCHM_SP(t2, "solve(crossprod(Lambda, ZtX), L, \"L\")");
	RZX.update(*t2);
	
	::M_cholmod_free_sparse(&t1, &c);
	t1 = RZX.crossprod();
//	showCHM_SP(t1, "RZX.crossprod()");
	t2 = ::M_cholmod_add(&XtX, t1, one, mone, 1/*values*/, 1/*sorted*/, &c);
	::M_cholmod_free_sparse(&t1, &c);
//	showCHM_SP(t2, "XtX-RZX.crossprod");

//	showCHM_FR(&RX, "RX before update");
	RX.update(*t2);
//	showCHM_FR(&RX, "RX after update");
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

    void glmerResp::updateSqrtRWt() {
	double *rr = sqrtrwt.begin(), *ww = weights.begin(), *vv = var.begin();
	int n = var.size();

	for (int i = 0; i < n; i++) rr[i] = sqrt(ww[i] / vv[i]);
    }

    void glmerResp::updateSqrtXWt() {
	std::transform(sqrtrwt.begin(), sqrtrwt.end(), muEta.begin(),
		       sqrtXwt.x.begin(), std::multiplies<double>());
    }

    void glmerDe::updateGamma() {
	std::copy(resp.offset.begin(), resp.offset.end(), resp.gamma.begin());
	fe.incGamma(resp.gamma);
	re.incGamma(resp.gamma);
    }

//    void rwDeFeMod::updateV(rwResp const &resp) {
    void rwDeFeMod::updateV(rwResp &resp) {
	int nc = resp.sqrtXwt.ncol();
	if (nc != 1) Rf_error("code for nlmer not yet written");
	std::copy(X.x.begin(), X.x.end(), V.x.begin());
	int m = X.nrow(), n = X.ncol();
	double *vv = V.x.begin(), *ww = resp.sqrtXwt.x.begin();
	for (int j = 0; j < n; j++)
	    for (int i = 0; i < m; i++) vv[i + j * m] *= ww[i];
    }

#define CM_TOL 1.e-5
#define CM_MAXITER 30

    double glmerDe::IRLS() {
	std::vector<double> varold(resp.var.size()), betaold(fe.beta.size());
	double crit, step, wrss0, wrss1;

	resp.updateSqrtRWt();
	crit = 10. * CM_TOL;
	for (int i = 0; crit >= CM_TOL && i < CM_MAXITER; i++) {
	    std::copy(resp.var.begin(), resp.var.end(), varold.begin());
	    std::copy(fe.beta.begin(), fe.beta.end(), betaold.begin());
	    updateGamma();	// linear predictor
	    resp.linkInv();	// mu
	    resp.MuEta();
	    resp.updateSqrtXWt();
	    wrss0 = resp.updateWrss();
	    
	    fe.updateV(resp);	// derive V from X
	    fe.RX.update(fe.V);	// factor V'V
	    
#if 0
	    fe.V.dgemv(
	    V->drmult(1/*transpose*/, 1, 0, cwtres, Vtwr);
	    deltaf = RXp->solveA(Vtwr);
	    if (verbose > 1) showdbl((double*)deltaf->x, "deltaf", p);
	    wrss1 = wrss0;		// force one evaluation of the loop
	    for (step = 1.; wrss0 <= wrss1 && step > CM_SMIN; step /= 2.) {
		for (int k = 0; k < p; k++)
		    beta[k] = betaold[k] + step * ((double*)deltaf->x)[k];
		update_gamma();
		linkinv();
		wrss1 = update_wtres();
		if (verbose > 1) {
		    Rprintf("step = %g, wrss0 = %g; wrss1 = %g\n",
			    step, wrss0, wrss1, beta[0]);
		    showdbl(beta, "beta", p);
		}
	    }
	    update_sqrtrwt();
	    crit = compare_vec(varold, var, n);
	    if (verbose > 1) Rprintf("convergence criterion: %g\n", crit);
	    VtV->freeA(); delete VtV;
	    V->freeA(); delete V;
	    M_cholmod_free_dense(&deltaf, &c);
#endif 
	}
//	return devres;
	return 0.;
    } // IRLS
}

RCPP_FUNCTION_2(double,lmerDeUpdate,S4 xp, NumericVector nt) {
    mer::lmerDe lm(xp);
    return lm.updateTheta(nt);
}

RCPP_FUNCTION_2(double,lmerSpUpdate,S4 xp, NumericVector nt) {
    mer::lmerSp lm(xp);
    return lm.updateTheta(nt);
}

RCPP_FUNCTION_1(double,lmerDeDeviance,S4 xp) {
    mer::lmerDe lm(xp);
    return lm.deviance();
}

RCPP_FUNCTION_1(double,lmerSpDeviance,S4 xp) {
    mer::lmerSp lm(xp);
    return lm.deviance();
}

RCPP_FUNCTION_1(double,glmerDeIRLS,S4 xp) {
    mer::glmerDe glmr(xp);
    return 0.;
}
