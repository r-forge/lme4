#include "merenv.h"
#include "merenv.hpp"

#include "lme4utils.hpp"

#if 0
static void showdbl(const double* x, const char* nm, int n) {
    int n5 = (n < 5) ? n : 5;
    Rprintf("%s: %g", nm, x[0]);
    for (int i = 1; i < n5; i++) Rprintf(", %g", x[i]);
    Rprintf("\n");
}
#endif

// Definition of methods for the merenv class
merenv::merenv(SEXP rho)
{
    if (!isEnvironment(rho))
	error(_("argument rho must be an environment"));

    // Extract slots that must have positive length.
    // Get dimensions of the problem
    SEXP sl = findVarBound(rho, lme4_ySym);
    if (!(n = LENGTH(sl)) || !isReal(sl)) // n = length of response
	error(_("Response vector y must be numeric (double) and non-empty"));
    y = REAL(sl);

    sl = findVarBound(rho, lme4_fixefSym);
    if (!isReal(sl))		// p = length of fixef
	error(_("fixef vector must be numeric (double)"));
    p = LENGTH(sl);
    fixef = REAL(sl);

    sl = findVarBound(rho, lme4_uSym);
    if (!(q = LENGTH(sl)) || !isReal(sl)) // q = length of u
	error(_("u vector must be numeric (double)"));
    u = REAL(sl);

    sl = findVarBound(rho, install("theta"));
    if (!(nth = LENGTH(sl)) || !isReal(sl)) // nth = length of theta
	error(_("theta vector must be numeric (double)"));
    theta = REAL(sl);

    sl = findVarBound(rho, lme4_LindSym);
    if (!(nLind = LENGTH(sl)) || !isInteger(sl))
	error(_("Lind vector must be integer"));
    Lind = INTEGER(sl);
    
    // allow for Z and X to have a multiple of n rows.
    Zt = VAR_CHM_SP(rho, lme4_ZtSym, q, 0);
    N = (int)Zt->ncol;
    Ut = VAR_CHM_SP(rho, lme4_UtSym, q, N);
    L = VAR_CHM_FR(rho, lme4_LSym, q);
				// versions of Lambda
///FIXME: Use S4 class to determine diagonalLambda
    diagonalLambda = asLogical(findVarBound(rho, install("diagonalLambda")));
    if (diagonalLambda) {
	Lambdax = VAR_dMatrix_x(rho, lme4_LambdaSym, q, q);
	Lambda = M_cholmod_speye((size_t)q, (size_t) q, CHOLMOD_REAL, &c);
	dble_cpy((double*)Lambda->x, Lambdax, q);
    } else {
	Lambda = VAR_CHM_SP(rho, lme4_LambdaSym, q, q);
    }
    int nnz = M_cholmod_nnz(Lambda, &c);
    if (nnz != nLind)
	error(_("length(Lind) = %d should be  %d"), nLind, nnz);

    gam = VAR_REAL_NULL(rho, install("gamma"), N);
    mu = VAR_REAL_NULL(rho, lme4_muSym, N);
    ldL2 = VAR_REAL_NULL(rho, install("ldL2"), 1);
    pwrss = VAR_REAL_NULL(rho, install("pwrss"), 1);
    weights = VAR_REAL_NULL(rho, install("weights"), n, TRUE);
    sqrtrwt = VAR_REAL_NULL(rho, install("sqrtrwt"), n, TRUE);
    sqrtXwt = VAR_REAL_NULL(rho, install("sqrtXwt"), n, TRUE);
    offset = VAR_REAL_NULL(rho, lme4_offsetSym, N, TRUE);
    if ((sparseX = asLogical(findVarBound(rho, install("sparseX"))))) {
	Xp = new CHM_rs(VAR_CHM_SP(rho, lme4_XSym, N, p));
	RXp = new Cholesky_rs(VAR_CHM_FR(rho, lme4_RXSym, p));
	RZXp = new CHM_rs(VAR_CHM_SP(rho, lme4_RZXSym, q, p));
    } else {
	Xp = new CHM_rd(VAR_CHM_DN(rho, lme4_XSym, N, p));
	RXp = new Cholesky_rd(findVarBound(rho, lme4_RXSym), p);
	RZXp = new CHM_rd(VAR_CHM_DN(rho, lme4_RZXSym, q, p));
    }	
}

merenv::~merenv(){
    if (diagonalLambda) {
	dble_cpy(Lambdax, (double*)Lambda->x, q);
	M_cholmod_free_sparse(&Lambda, &c);
    }
    delete Xp;
    delete RXp;
    delete RZXp;
}

CHM_r *merenv::Vp() {
    CHM_r *ans;
    if (sparseX) {
	CHM_SP A = M_cholmod_copy_sparse((dynamic_cast<CHM_rs*>(Xp))->A, &c);
	if (sqrtXwt) {
	    M_cholmod_scale(N_AS_CHM_DN(sqrtXwt,A->nrow,1), CHOLMOD_ROW, A, &c);
	}
	ans = new CHM_rs(A);
    } else {
	CHM_DN A = M_cholmod_copy_dense((dynamic_cast<CHM_rd*>(Xp))->A, &c);
	if (sqrtXwt) {
	    int m = (int)A->nrow, n = (int)A->ncol;
	    double *ax = (double*)A->x;
	    for (int j = 0; j < n; j++)
		for (int i = 0; i < m; i++) ax[i + j*m] *= sqrtXwt[j];
	}
	ans = new CHM_rd(A);
    }
    return ans;
}
void merenv::update_gamma() {
    static double one[] = {1,0};
    CHM_DN cu = N_AS_CHM_DN(u, q, 1),
	cgamma = N_AS_CHM_DN(gam, N, 1),
	cfixef = N_AS_CHM_DN(fixef, p, 1);
    if (offset)			// initialize to offset if used
	dble_cpy(gam, offset, N);
    else			// otherwise initialize to zero
	dble_zero(gam, N);
    M_cholmod_sdmult(Ut, 1/*transpose*/, one, one, cu, cgamma, &c);
    Xp->drmult(0/*transpose*/, 1, 1, cfixef, cgamma);
}

double merenv::update_pwrss() {
    *pwrss = sqr_length(u, q);
    for (int i = 0; i < n; i++) {
	double resi = (y[i] - mu[i]) * (sqrtrwt ? sqrtrwt[i] : 1);
	*pwrss += resi * resi;
    }
    return *pwrss;
}

void merenv::update_Lambda_Ut(SEXP thnew) {
				// check and copy thnew
    if (!isReal(thnew) || LENGTH(thnew) != nth)
	error(_("theta must be numeric and of length %d"), nth);
    dble_cpy(theta, REAL(thnew), nth);
				// update Lambda
    double *Lambdax = (double*)Lambda->x;
    for (int i = 0; i < nLind; i++) Lambdax[i] = theta[Lind[i] - 1];
				// update Ut from Lambda and Zt
    CHM_SP tmp;
    int m = (int)Ut->nrow, n = (int)Ut->ncol;
    if (diagonalLambda) {
	tmp = M_cholmod_ssmult(Lambda, Zt, 0/*stype*/, TRUE/*values*/,
			       TRUE/*sorted*/, &c);
    } else {
	CHM_SP Lamtr = M_cholmod_transpose(Lambda, 1/*values*/, &c);
	tmp = M_cholmod_ssmult(Lamtr, Zt, 0/*stype*/, TRUE/*values*/,
			       TRUE/*sorted*/, &c);
	M_cholmod_free_sparse(&Lamtr, &c);
    }
    if (!(m == (int)tmp->nrow && n == (int)tmp->ncol))
	error(_("Lambda'Zt (%d by %d) doesn't match Ut (%d by %d)"),
	      tmp->nrow, tmp->ncol, m, n);
    int *di = (int*)(Ut->i), *si = (int*)(tmp->i),
	*dp = (int*)(Ut->p), *sp = (int*)(tmp->p);
    double *dx = (double*)(Ut->x), *sx = (double*)(tmp->x); 
    
    for (int j = 0; j <= n; j++)
	if (dp[j] != sp[j])
	    error(_("Ut->p[%d] not consistent with Lambda'Zt"), j);
    for (int j = 0; j < dp[n]; j++) {
	if (di[j] != si[j])
	    error(_("Ut->i[%d] not consistent with Lambda'Zt"), j);
	dx[j] = sx[j];
    }
    M_cholmod_free_sparse(&tmp, &c);
// Note: we do not update L here because there may be weights to consider
// Well, that depends on how U is defined.  For that we may need sqrtXwts.
}

CHM_DN merenv::crossprod_Lambda(CHM_DN rhs, CHM_DN ans) {
    static double one[] = {1,0}, zero[] = {0,0};
    M_cholmod_sdmult(Lambda, 1/*transpose*/, one, zero, rhs, ans, &c);
    return ans;
}

CHM_SP merenv::spcrossprod_Lambda(CHM_SP src) {
    CHM_SP Lamtr = M_cholmod_transpose(Lambda, TRUE/*values*/, &c),
	tmp = M_cholmod_ssmult(Lamtr, Zt, 0/*stype*/,
			       TRUE/*values*/, TRUE/*sorted*/, &c);
    M_cholmod_free_sparse(&Lamtr, &c);
    return tmp;
}

CHM_DN merenv::solvePL(CHM_DN src) {
    CHM_DN tmp1, tmp2;

    tmp1 = M_cholmod_copy_dense(src, &c);
    crossprod_Lambda(src, tmp1);
    tmp2 = M_cholmod_solve(CHOLMOD_P, L, tmp1, &c);
    M_cholmod_free_dense(&tmp1, &c);
    tmp1 = M_cholmod_solve(CHOLMOD_L, L, tmp2, &c);
    M_cholmod_free_dense(&tmp2, &c);
    return tmp1;
}

CHM_SP merenv::solvePL(CHM_SP src) {
    CHM_SP tmp1, tmp2;

    tmp1 = spcrossprod_Lambda(src);
    tmp2 = M_cholmod_spsolve(CHOLMOD_P, L, tmp1, &c);
    M_cholmod_free_sparse(&tmp1, &c);
    tmp1 = M_cholmod_spsolve(CHOLMOD_L, L, tmp2, &c);
    M_cholmod_free_sparse(&tmp2, &c);
    return tmp1;
}

// definitions of methods for the merenvtrms classes

merenvtrms::merenvtrms(SEXP rho) : merenv(rho) {
    flist = findVarBound(rho, lme4_flistSym);
    if (!isNewList(flist))
	error(_("Object \"%s\" must be a list"), "flist");
    nfac = LENGTH(flist);

    SEXP cnms = findVarBound(rho, install("cnms"));
    if (!isNewList(cnms))
	error(_("Object \"cnms\" must be a list"));
    ntrm = LENGTH(cnms);

    nc = (int*)R_alloc(sizeof(int), ntrm);
    nl = (int*)R_alloc(sizeof(int), nfac);
    apt = (int*)R_alloc(sizeof(int), nfac + 1);
    apt[0] = 0;

    SEXP asgn = getAttrib(flist, install("assign"));
    if (LENGTH(asgn) != ntrm)
	error(_("length(attr(flist, \"assign\")) != length(cnms)"),
	      LENGTH(asgn), ntrm);
    int *assign = INTEGER(asgn);

    for (int i = 0; i < nfac; i++) {
	SEXP ff = VECTOR_ELT(flist, i);
	if (!isFactor(ff))
	    error(_("Element %d of flist is not a factor"), i + 1);
	if (LENGTH(ff) != n)
	    error(_("Element %d has length %d; should be %d"),
		  i + 1, LENGTH(ff), n);
	nl[i] = LENGTH(getAttrib(ff, R_LevelsSymbol));
    }
				// check range in assign, store nc and apt
    for (int i = 0; i < ntrm; i++) {
	int ii = assign[i];
	if (ii < 1 || nfac < ii)
	    error(_("assign attribute els must be in [1,%d]"), nfac);
	if (i > 0 && assign[i - 1] > assign[i])
	    error(_("assign attribute must be non-decreasing"));
	nc[i] = LENGTH(VECTOR_ELT(cnms, i));
	apt[ii] = i + 1;
    }
    if (apt[nfac] != ntrm)
	error(_("assign attribute does not agree with cnms"));
    for (int i = 0; i < nfac; i++)
	if (apt[i] >= apt[i + 1])
	    error(_("assign attribute missing index %d"), i + 1);
}

void merenvtrms::show() {
    Rprintf("merenvtrms object: ntrm = %d, nfac = %d\n nl:", ntrm, nfac);
    for (int i = 0; i < nfac; i++) Rprintf(" %d", nl[i]);
    Rprintf("\n apt:");
    for (int i = 0; i <= nfac; i++) Rprintf(" %d", apt[i]);
    Rprintf("\n nc:");
    for (int i = 0; i < ntrm; i++) Rprintf(" %d", nc[i]);
    Rprintf("\n");
}
    
SEXP merenvtrms::condVar(double scale) {
    SEXP ans = PROTECT(allocVector(VECSXP, nfac));
    int offset = 0;
    
    if (scale < 0 || !R_FINITE(scale))
	error(_("scale must be non-negative and finite"));

    double scsqr = scale * scale;
    for (int i = 0; i < nfac; i++) {
	int api = apt[i], nli = nl[i],
	    ntrm = apt[i + 1] - api; // number of terms for this factor
	int *ncol = new int[ntrm], *cumcol = new int[ntrm + 1];
	cumcol[0] = 0;
	for (int j = 0; j < ntrm; j++)
	    cumcol[j + 1] = cumcol[j] + (ncol[j] = nc[api + j]);
	int nct = cumcol[ntrm];
	int *cset = new int[nct], nct2 = nct * nct;
	
	SET_VECTOR_ELT(ans, i, alloc3DArray(REALSXP, nct, nct, nli));
	double *ai = REAL(VECTOR_ELT(ans, i));

	for (int j = 0; j < nli; j++) {
	    for (int jj = 0; jj < ntrm; jj++)
		for (int k = 0; k < ncol[jj]; k++)
		    cset[k] = offset + j * ncol[jj] + cumcol[jj] * nli + k;

	    CHM_SP cols =
		M_cholmod_submatrix(Lambda, (int*)NULL, -1, cset, nct,
				    1/*values*/, 1/*sorted*/, &c);
	    CHM_SP sol = M_cholmod_spsolve(CHOLMOD_A, L, cols, &c);
	    CHM_SP tcols = M_cholmod_transpose(cols, 1/*values*/, &c);
	    M_cholmod_free_sparse(&cols, &c);
	    CHM_SP var = M_cholmod_ssmult(tcols, sol, 0/*stype*/,
					  1/*values*/, 1/*sorted*/, &c);
	    M_cholmod_free_sparse(&sol, &c);
	    M_cholmod_free_sparse(&tcols, &c);
	    CHM_DN dvar = M_cholmod_sparse_to_dense(var, &c);
	    M_cholmod_free_sparse(&var, &c);
	    Memcpy(ai + j * nct2, (double*)dvar->x, nct2);
	    M_cholmod_free_dense(&dvar, &c);
	}
	for (int k = 0; k < nct2 * nli; k++) ai[k] *= scsqr;
	offset += nli * nct;
	delete[] cset;
	delete[] cumcol;
	delete[] ncol;
    }

    UNPROTECT(1);
    return ans;
}

lmer::lmer(SEXP rho) : merenv(rho) {
    if (N != n)
	error(_("nrow(X) = %d must match length(y) = %d for lmer"),
	      N, n);
    REML = asLogical(findVarBound(rho, install("REML")));
    ldRX2 = VAR_REAL_NULL(rho, install("ldRX2"), 1, TRUE);
    Xty = VAR_REAL_NULL(rho, install("Xty"), p);
    Zty = VAR_dMatrix_x(rho, install("Zty"), q, 1);
    if (sparseX) {
	XtXp = new CHM_rs(VAR_CHM_SP(rho, install("XtX"), p, p));
	ZtXp = new CHM_rs(VAR_CHM_SP(rho, install("ZtX"), q, p));
    } else {
	XtXp = new dpoMatrix(findVarBound(rho, install("XtX")));
	if (!(XtXp->nrow() == p && XtXp->ncol() == p))
	    error(_("Dimensions of %s are %d by %d, should be %d by %d"),
		  "XtX", XtXp->nrow(), XtXp->ncol(), p, p);
	ZtXp = new CHM_rd(VAR_CHM_DN(rho, install("ZtX"), q, p));
    }
}

lmer::~lmer() {
    delete XtXp;
    delete ZtXp;
}

double lmer::update_dev(SEXP thnew) {
    static double one[] = {1,0};
    update_Lambda_Ut(thnew);
    M_cholmod_factorize_p(Ut, one, (int*)NULL, (size_t)0, L, &c);
    *ldL2 = M_chm_factor_ldetL2(L);
    *ldRX2 = 0;			// in case p == 0
    cu = solvePL(N_AS_CHM_DN(Zty, q, 1));
    if (p) {
	CHM_r *tmp1 = ZtXp->crossprod_SP(Lambda);
	CHM_r *tmp2 = tmp1->solveCHM_FR(L, CHOLMOD_P);
	tmp1->freeA(); delete tmp1;
	tmp1 = tmp2->solveCHM_FR(L, CHOLMOD_L);
	tmp2->freeA(); delete tmp2;
	RZXp->copy_contents(tmp1);
	tmp1->freeA(); delete tmp1;
	RXp->downdate(RZXp, -1.0, XtXp, 1.0);
	*ldRX2 = RXp->ldet2();
	CHM_DN tt = M_cholmod_copy_dense(N_AS_CHM_DN(Xty, p, 1), &c);
	RZXp->drmult(1/*transpose*/, -1.0, 1.0, cu, tt);
	CHM_DN ff = RXp->solveA(tt);
	M_cholmod_free_dense(&tt, &c);
	dble_cpy(fixef, (double*)ff->x, p);
	RZXp->drmult(0/*no trans*/, -1.0, 1.0, ff, cu);
	M_cholmod_free_dense(&ff, &c);
    }
    CHM_r *td1 = new CHM_rd(cu);
    CHM_r *td2 = td1->solveCHM_FR(L, CHOLMOD_Lt);
    td1->freeA(); delete td1;
    td1 = td2->solveCHM_FR(L, CHOLMOD_Pt);
    td2->freeA(); delete td2;
    dble_cpy(u, (double*)(dynamic_cast<CHM_rd*>(td1)->A)->x, q);
    td1->freeA(); delete td1;
    update_gamma();
    dble_cpy(mu, gam, n);
    update_pwrss();
    if (REML) {
	double nmp = (double)(n - p);
	return *ldL2 + *ldRX2 + nmp * (1 + log(2 * PI * (*pwrss)/nmp));
    }				       
    return *ldL2 + n * (1 + log(2 * PI * (*pwrss)/((double)n)));
}    

glmer::glmer(SEXP rho) : merenv(rho) {
    muEta = VAR_REAL_NULL(rho, install("muEta"), n);
    var = VAR_REAL_NULL(rho, install("var"), n);
    wtres = VAR_REAL_NULL(rho, install("wtres"), n);
    fam.initGL(rho);
}

void glmer::update_sqrtrwt() {
    varFunc();
    for (int j = 0; j < n; j++)	{
	sqrtrwt[j] = sqrt((weights ? weights[j] : 1) / var[j]);
	sqrtXwt[j] = muEta[j] * sqrtrwt[j];
    }
}

double glmer::update_wtres() {
    double ans = 0;
    for (int j = 0; j < n; j++) {
	wtres[j] = sqrtrwt[j] * (y[j] - mu[j]);
	ans += wtres[j] * wtres[j];
    }
    return ans;
}

double glmer::PIRLS() {
    return 1;
}

double glmer::PIRLSbeta() {
    int cvg = 1;//, info;
    double *betaold = new double[p], *varold = (double*)NULL,
	*cbeta = new double[p],
	*tmp = new double[q], *uold = new double[q];
//	cfac = ((double)n)/((double)(q+p));//, crit, pwrss_old,
//	step;
//    static double one[] = {1,0}, zero[] = {0,0}, mone[] = {-1,0};
//    CHM_DN SOL, cRZX = N_AS_CHM_DN(RZX, q, p), cV,
//	cwtres = N_AS_CHM_DN(wtres, n, 1), ctmp = N_AS_CHM_DN(tmp, q, 1);
//    CHM_SP U;
//    R_CheckStack();
    if (!cvg) error(_("Convergence failure in PIRLS"));
    
//    delete[] V;
    delete[] betaold; delete[] cbeta;
    delete[] tmp; delete[] uold; 
    if (var) delete[] varold; 

    return 1;
}

static const int CM_MAXITER = 300;
static const double CM_SMIN = 0.001;

double glmer::IRLS() {
    CHM_r *V, *VtV;    
    CHM_DN cwtres = N_AS_CHM_DN(wtres, n, 1), deltaf,
	Vtwr = M_cholmod_allocate_dense((size_t)p, (size_t)1, (size_t)p,
				       CHOLMOD_REAL, &c);
    double *betaold = new double[p], step, wrss0, wrss1;
    
    update_sqrtrwt();
    for (int i = 0; i < CM_MAXITER; i++) { // iterate until weights stabilize
	dble_cpy(betaold, fixef, p);
	V = Vp();		// derive V from X
	update_gamma();		// linear predictor
	linkinv();		// mu and muEta
	VtV = V->AtA();		// crossproduct
	RXp->update(VtV);	// factor
	
	wrss0 = update_wtres(); 
	V->drmult(1/*transpose*/, 1, 0, cwtres, Vtwr);
	deltaf = RXp->solveA(Vtwr);
	wrss1 = wrss0;		// force one evaluation of the loop
	for (step = 1.; wrss0 <= wrss1 && step > CM_SMIN; step /= 2.) {
	    for (int k = 0; k < p; k++)
		fixef[k] = betaold[k] + step * ((double*)deltaf->x)[k];
	    update_gamma();
	    linkinv();
	    wrss1 = update_wtres();
	    Rprintf("step = %g, wrss0 = %g; wrss1 = %g\nfixef: %g",
		    step, wrss0, wrss1, fixef[0]);
	    for (int k = 1; k < p; k++) Rprintf(", %g", fixef[k]);
	    Rprintf("\n");
	}
	VtV->freeA(); delete VtV;
	V->freeA(); delete V;
	M_cholmod_free_dense(&deltaf, &c);
	break;
    }
    delete[] betaold;
    M_cholmod_free_dense(&Vtwr, &c);
    return 1;
}

#if 0
double glmerold::PIRLSbeta() {
    int cvg, info;
    double *betaold = new double[p], *varold = (double*)NULL,
	*cbeta = new double[p],
	*tmp = new double[q], *uold = new double[q],
	*wtres = new double[n],
	cfac = ((double)n)/((double)(q+p)), crit, pwrss_old,
	step;
    static double d1[2] = {1,0}, d0[2] = {0,0}, dm1[2] = {-1,0};
    CHM_DN SOL, cRZX = N_AS_CHM_DN(RZX, q, p), cV,
	cwtres = N_AS_CHM_DN(wtres, n, 1), ctmp = N_AS_CHM_DN(tmp, q, 1);
    CHM_SP U;
    R_CheckStack();

    // reset u to initial values.  This can result in more
    // iterations but it gives a repeatable function evaluation for
    // the optimizer.
    dble_zero(u, q); 
    dble_zero(fixef, p);

    V = new double[n * p];
    cV = N_AS_CHM_DN(V, n, p);

    cvg = FALSE;

    for (int ii = 0; ii < CM_MAXITER; ii++)
    {
	for (int i = 0; i < CM_MAXITER; i++) {
	    dble_cpy(uold, u, q); // record current coefficients
	    dble_cpy(betaold, fixef, p);

	    pwrss_old = update_mu();
	    U = A_to_U();	// create U
	    if (!M_cholmod_factorize_p(U, d1, (int*)NULL, 0 /*fsize*/, L, &c))
		error(_("cholmod_factorize_p failed: status %d, minor %d from ncol %d"),
		      c.status, L->minor, L->n);
	    X_to_V();			// update V
	    if (!M_cholmod_sdmult(U, 0,	// no transpose
				  d1, d0, cV, cRZX, &c))
		error(_("cholmod_sdmult failed: status %d"), c.status);
	    if (!(SOL = M_cholmod_solve(CHOLMOD_L, L, cRZX, &c)))
		error(_("cholmod_solve (CHOLMOD_L) failed: status %d"), c.status);
	    dble_cpy(RZX, (double*)(SOL->x), q * p);
	    M_cholmod_free_dense(&SOL, &c);
				// solve for RX in downdated V'V
	    F77_CALL(dsyrk)("U", "T", &p, &n, d1, V, &n, d0, RX, &p); //V'V
	    F77_CALL(dsyrk)("U", "T", &p, &q, dm1, RZX, &q, d1, RX, &p);
	    F77_CALL(dpotrf)("U", &p, RX, &p, &info);
	    if (info)
		error(_("Downdated V'V is not positive definite, %d."), info);
				// tmp := U %*% wtdResid 
	    M_cholmod_sdmult(U, 0, // no transpose
			     d1, d0, cwtres, ctmp, &c);
	    for (int j = 0; j < q; j++) tmp[j] -= u[j]; // tmp := tmp - u
	    M_cholmod_free_sparse(&U, &c);
				// solve L %*% tmp = tmp
	    if (!(SOL = M_cholmod_solve(CHOLMOD_L, L, ctmp, &c)))
		error(_("cholmod_solve (CHOLMOD_L) failed"));
	    dble_cpy(tmp, (double*)(SOL->x), q);
	    M_cholmod_free_dense(&SOL, &c);
				// solve RX'cbeta = V'wtres - RZX'cu
	    F77_CALL(dgemv)("T", &n, &p, d1, V, &n, wtres, &i1,
			    d0, cbeta, &i1);
	    F77_CALL(dgemv)("T", &q, &p, dm1, RZX, &q, tmp, &i1,
			    d1, cbeta, &i1);
	    F77_CALL(dtrsv)("U", "T", "N", &p, RX, &p, cbeta, &i1);
				// evaluate convergence criterion
	    double cul2 = sqr_length(tmp, q), cbetal2 = sqr_length(cbeta, p);
	    crit = cfac * (cul2 + cbetal2)/ pwrss_old;
//	    if (verb < 0) Rprintf("cul2 = %g, cbetal2 = %g\n", cul2, cbetal2);
	    if (crit < CM_TOL) {	// don't do needless evaluations 
		cvg = TRUE;
		break;
	    }
				// solve for delta-beta
	    F77_CALL(dtrsv)("U", "N", "N", &p, RX, &p, cbeta, &i1);
				// solve t(L) %*% SOL = tmp - RZX cbeta
	    F77_CALL(dgemv)("N", &q, &p, dm1, RZX, &q, cbeta, &i1,
			    d1, tmp, &i1);
	    if (!(SOL = M_cholmod_solve(CHOLMOD_Lt, L, ctmp, &c)))
		error(_("cholmod_solve (CHOLMOD_Lt) failed"));
	    dble_cpy(tmp, (double*)(SOL->x), q);
	    M_cholmod_free_dense(&SOL, &c);
	    
	    for (step = 1; step > CM_SMIN; step /= 2) { // step halving 
		for (int j = 0; j < q; j++)
		    u[j] = uold[j] + step * tmp[j];
		for (int j = 0; j < p; j++)
		    fixef[j] = betaold[j] + step * cbeta[j];
		d[pwrss_POS] = update_mu();
		if (verb < 0)
		    Rprintf("%2d,%8.6f,%12.4g: %15.6g %15.6g %15.6g %15.6g %15.6g\n",
			    i, step, crit, d[pwrss_POS], pwrss_old, fixef[0], u[0], u[1]);
		if (d[pwrss_POS] < pwrss_old) {
		    pwrss_old = d[pwrss_POS];
		    break;
		}
	    }
	    if (step <= CM_SMIN) break;
	    if (!(muEta || etaGamma)) { // Gaussian linear mixed models require
		cvg = TRUE;		// only 1 iteration
		break;
	    }
	}
	if (!cvg) break;
	if (!var) break;
	eval_varFunc();
	crit = 0;
	for (int j = 0; j < n; j++) {
	    double diff = varold[j] - var[j];
	    crit += (diff * diff)/(var[j]*var[j] + varold[j]*varold[j]);
	}
	if (verb < 0) Rprintf("%2d,%8.6f %15.6g %15.6g %15.6g %15.6g\n",
			      ii, crit, var[0], var[1], varold[0], varold[1]);
	if (crit < CM_TOL) break;
    }
    
    if (!cvg) error(_("Convergence failure in PIRLS"));
    
    delete[] V; delete[] betaold; delete[] cbeta;
    delete[] tmp; delete[] uold; delete[]wtres;
    if (var) delete[] varold; 
    
    d[ldRX2_POS] = 0;
    for (int j = 0; j < p; j++) d[ldRX2_POS] += 2 * log(RX[j * (p + 1)]);
    d[ldL2_POS] = M_chm_factor_ldetL2(L);
    d[sigmaML_POS] = sqrt(d[pwrss_POS]/(srwt ? sqr_length(srwt, n) : (double) n));
    d[sigmaREML_POS] = (etaGamma || muEta) ? NA_REAL :
	d[sigmaML_POS] * sqrt((((double) n)/((double)(n - p))));
    
    return update_dev();
    return 0;
}

double glmerold::PIRLS() {
    update_sqrtrwt();		// update variance and sqrtrwt
    return 0;
}

double glmerold::IRLS() {
    double *betaold = new double[p], *varold = new double[n];
    update_sqrtrwt();		// update variance and sqrtrwt
    for (int ii = 0; ii < CM_MAXITER; ii++) {
	dble_cpy(varold, var, n);
	for (int i = 0; i < CM_MAXITER; i++) {
	    dble_cpy(betaold, fixef, p);
	}
    }	
    delete[] betaold;
    delete[] varold;
    return 0;
}
#endif

extern "C" {

    SEXP glmer_PIRLS(SEXP rho) {
	return ScalarReal(glmer(rho).PIRLS());
    }

    SEXP glmer_IRLS(SEXP rho) {
	return ScalarReal(glmer(rho).IRLS());
    }

    SEXP glmer_PIRLSbeta(SEXP rho) {
	return ScalarReal(glmer(rho).PIRLSbeta());
    }
	
    SEXP glmer_update_sqrtrwt(SEXP rho) {
	glmer(rho).update_sqrtrwt();
	return R_NilValue;
    }


    SEXP lmer_validate(SEXP rho) {
	return ScalarLogical(lmer(rho).validate());
    }

    SEXP lmer_deviance(SEXP rho, SEXP thnew) {
	return ScalarReal(lmer(rho).update_dev(thnew));
    }

/**
 * Duplicate selected names from an environment
 *
 * @param dest environment in which to write the copies
 * @param src source environment
 * @param nms character vector of names of objects to copy
 * @return dest
 */
    SEXP lme4_dup_env_contents(SEXP dest, SEXP src, SEXP nms) {
	if (!isEnvironment(dest))
	    error(_("dest must be an environment"));
	if (!isEnvironment(src))
	    error(_("src must be an environment"));
	if (!isString(nms))
	    error(_("nms must be a character variable"));
	for (int i = 0; i < LENGTH(nms); i++) {
	    SEXP nmsym = install(CHAR(STRING_ELT(nms, i)));
	    defineVar(nmsym, duplicate(findVarBound(src, nmsym)), dest);
	}
	return dest;
    }

/// Check validity of an merenvtrms environment
    SEXP merenvtrms_validate(SEXP rho) {
	return ScalarLogical(merenvtrms(rho).validate());
    }

/**
 * Create the conditional variance arrays
 *
 * @param rho pointer to an merenvtrms environment
 * @param scale a numeric scalar -- the non-negative scale factor
 *
 * @return a list of 3D arrays corresponding to the factors in flist
 */
    SEXP merenvtrms_condVar(SEXP rho, SEXP scale) {
	return merenvtrms(rho).condVar(asReal(scale));
    }
    SEXP merenvtrms_show(SEXP rho) {
	merenvtrms(rho).show();
	return R_NilValue;
    }
	
}

