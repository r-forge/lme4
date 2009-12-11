#include "merenv.h"
#include "merenv.hpp"

#include "lme4utils.hpp"

// Definition of methods for the merenv class

//merenv::merenv(SEXP rho)
void merenv::initMer(SEXP rho)
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

    sl = findVarBound(rho, install("Lind"));
    if (!(nLind = LENGTH(sl)) || !isInteger(sl))
	error(_("Lind vector must be integer"));
    Lind = INTEGER(sl);
    
    // allow for Z and X to have a multiple of n rows.
    Zt = VAR_CHM_SP(rho, lme4_ZtSym, q, 0);
    N = (int)Zt->ncol;
    Ut = VAR_CHM_SP(rho, lme4_UtSym, q, N);
    L = new cholmod_factor;
    M_as_cholmod_factor(L, findVarInFrame(rho, lme4_LSym));
				// versions of Lambda
//FIXME: Use S4 class to determine diagonalLambda
    if (asLogical(findVarBound(rho, install("diagonalLambda")))) {
	Lambdax = VAR_dMatrix_x(rho, lme4_LambdaSym, q, q);
	if (nLind != q)
	    error(_("Lind should be of length q = %d for diagonal Lambda"),
		  q);
	Lambda = (CHM_SP) NULL;
    } else {
	Lambda = VAR_CHM_SP(rho, lme4_LambdaSym, q, q);
	int nnz = M_cholmod_nnz(Lambda, &c);
	if (nnz != nLind)
	    error(_("length(Lind) = %d should be  %d"), nLind, nnz);
	Lambdax = (double*)(Lambda->x);
    }
    eta = VAR_REAL_NULL(rho, install("fitted"), N);
    ldL2 = VAR_REAL_NULL(rho, install("ldL2"), 1);
    prss = VAR_REAL_NULL(rho, install("prss"), 1);
    weights = VAR_REAL_NULL(rho, install("weights"), n, TRUE);
    offset = VAR_REAL_NULL(rho, lme4_offsetSym, N, TRUE);
//    Rprintf(
//	"In merenv(SEXP), dimensions are N = %d, n = %d, p = %d, q = %d\n",
//	N, n, p, q);
}

void merenv::update_Lambda_Ut(SEXP thnew) {
				// check and copy thnew
    if (!isReal(thnew) || LENGTH(thnew) != nth)
	error(_("theta must be numeric and of length %d"), nth);
    dble_cpy(theta, REAL(thnew), nth);
				// update Lambda
    for (int i = 0; i < nLind; i++) Lambdax[i] = theta[Lind[i] - 1];
				// update Ut from Lambda and Zt
    if (Lambda) {
	CHM_SP Lamtr = M_cholmod_transpose(Lambda, TRUE/*values*/, &c),
	    tmp = M_cholmod_ssmult(Lamtr, Zt, 0/*stype*/,
				   TRUE/*values*/, TRUE/*sorted*/, &c);
	M_cholmod_free_sparse(&Lamtr, &c);
	// Should we store t(Lambda), instead of Lambda?
	int *it = (int*)(tmp->i), *iu = (int*)(Ut->i),
	    *pt = (int*)(tmp->p), *pu = (int*)(Ut->p);
	double *xt = (double*)(tmp->x), *xu = (double*)(Ut->x); 

	for (int j = 0; j <= n; j++)
	    if (pt[j] != pu[j])
		error(_("Ut is not consistent with Lambda %*% Zt"));
	for (int j = 0; j < pu[n]; j++) {
	    if (it[j] != iu[j])
		error(_("Ut is not consistent with Lambda %*% Zt"));
	    xu[j] = xt[j];
	}
	M_cholmod_free_sparse(&tmp, &c);
    } else {			// special case for diagonal Lambda
	int *iz = (int*)(Zt->i), nnz = ((int*)(Zt->p))[n];
	double *xu = (double*)(Ut->x), *xz = (double*)(Zt->x);
	for (int i = 0; i < nnz; i++) xu[i] = xz[i] * Lambdax[iz[i]];
    }
}

double merenv::update_prss() {
    *prss = sqr_length(u, q);
    for (int i = 0; i < n; i++) {
	double resi = y[i] - eta[i];
	*prss += resi * resi;
    }
    return *prss;
}

void merenv::update_eta_Ut() {    
    CHM_DN cu = N_AS_CHM_DN(u, q, 1),
	ceta = N_AS_CHM_DN(eta, N, 1);
    if (offset)			// initialize to offset if used
	dble_cpy(eta, offset, N);
    else			// otherwise initialize to zero
	dble_zero(eta, N);
    M_cholmod_sdmult(Ut, 1/*transpose*/, &one, &one, cu, ceta, &c);
}

CHM_DN merenv::crossprod_Lambda(CHM_DN rhs, CHM_DN ans) {
    if (((int)(rhs->nrow)) != q)
	error(_("in crossprod_Lambda, rhs->nrow = %d, should be %d"),
	      rhs->nrow, q);
    if (Lambda) {
	M_cholmod_sdmult(Lambda, 1/*transpose*/, &one, &zero,
			 rhs, ans, &c);
    } else {
	double *ax = (double*)(ans->x), *rx = (double*)(rhs->x);
	int nc = rhs->ncol;
    	for (int j = 0; j < nc * q; j++) ax[j] = rx[j] * Lambdax[j % q];
    }
    return ans;
}

CHM_SP merenv::spcrossprod_Lambda(CHM_SP src) {
    if (((int)(src->nrow)) != q)
	error(_("in spcrossprod_Lambda, src->nrow = %d, should be %d"),
	      src->nrow, q);
    if (Lambda) {
	CHM_SP Lamtr = M_cholmod_transpose(Lambda, TRUE/*values*/, &c),
	    tmp = M_cholmod_ssmult(Lamtr, Zt, 0/*stype*/,
				   TRUE/*values*/, TRUE/*sorted*/, &c);
	M_cholmod_free_sparse(&Lamtr, &c);
	return tmp;
    }
				// special case for diagonal Lambda
    CHM_SP ans = M_cholmod_copy_sparse(src, &c);
    CHM_DN lambda = N_AS_CHM_DN(Lambdax, q, 1);
    if (!M_cholmod_scale(lambda, CHOLMOD_ROW, ans, &c))
	error(_("Error return from cholmod_scale"));
    return ans;
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

// Definition of methods for the mersparse and merdense classes

void mersparse::update_eta() {
    update_eta_Ut();
    M_cholmod_sdmult(X, 0/*no transpose*/, &one, &one,
		     N_AS_CHM_DN(fixef, p, 1),
		     N_AS_CHM_DN(eta, N, 1), &c);
}

void merdense::update_eta() {
    update_eta_Ut();
    if (p) 
	F77_CALL(dgemv)("N", &n, &p, &one, X, &n, fixef,
			&i1, &one, eta, &i1);
}

merdense::merdense(SEXP rho) {
//    Rprintf(
//	"In merdense(SEXP), dimensions are N = %d, n = %d, p = %d, q = %d\n",
//	N, n, p, q);
    X = VAR_dMatrix_x(rho, lme4_XSym, N, p);
    RX = VAR_dMatrix_x(rho, lme4_RXSym, p, p);
    RZX = VAR_dMatrix_x(rho, lme4_RZXSym, q, p);
}

mersparse::mersparse(SEXP rho) {
//    Rprintf(
//	"In mersparse(SEXP), dimensions are N = %d, n = %d, p = %d, q = %d\n",
//	N, n, p, q);
    X = VAR_CHM_SP(rho, lme4_XSym, N, p);
    RX = VAR_CHM_SP(rho, lme4_RXSym, p, p);
    RZX = VAR_CHM_SP(rho, lme4_RZXSym, q, p);
}

// definition of methods for the lmer, lmerdense and lmersparse classes

lmer::lmer(SEXP rho) {
    initMer(rho);
//    Rprintf("In lmer(SEXP), dimensions are N = %d, n = %d, p = %d, q = %d\n",
//	    N, n, p, q);
    if (N != n)
	error(_("nrow(X) = %d must match length(y) = %d for lmer"),
	      N, n);
    REML = asLogical(findVarBound(rho, install("REML")));
    ldRX2 = VAR_REAL_NULL(rho, install("ldRX2"), 1, TRUE);
    Xty = VAR_REAL_NULL(rho, install("Xty"), p);
    Zty = VAR_dMatrix_x(rho, install("Zty"), q, 1);
}

lmerdense::lmerdense(SEXP rho) : lmer(rho), merdense(rho) {
//    Rprintf(
//	"In lmerdense(SEXP), dimensions are N = %d, n = %d, p = %d, q = %d\n",
//	    N, n, p, q);
    ZtX = VAR_dMatrix_x(rho, install("ZtX"), q, p);
    XtX = VAR_dMatrix_x(rho, install("XtX"), p, p);
}

lmersparse::lmersparse(SEXP rho) : lmer(rho), mersparse(rho) {
//    Rprintf(
//	"In lmersparse(SEXP), dimensions are N = %d, n = %d, p = %d, q = %d\n",
//	    N, n, p, q);
    ZtX = VAR_CHM_SP(rho, install("ZtX"), q, p);
    XtX = VAR_CHM_SP(rho, install("XtX"), p, p);
}

void lmer::LMMdev1() {		// update L, create cu
    CHM_DN cZty = N_AS_CHM_DN(Zty, q, 1);
    M_cholmod_factorize_p(Ut, &one, (int*)NULL, (size_t)0, L, &c);
    *ldL2 = M_chm_factor_ldetL2(L);
    cu = solvePL(cZty);
}

void lmer::LMMdev2() {		// solve for u
    CHM_DN tmp1, tmp2;
    tmp1 = M_cholmod_solve(CHOLMOD_Lt, L, cu, &c);
    M_cholmod_free_dense(&cu, &c);
    tmp2 = M_cholmod_solve(CHOLMOD_Pt, L, tmp1, &c);
    M_cholmod_free_dense(&tmp1, &c);
    dble_cpy(u, (double*)(tmp2->x), q);
    M_cholmod_free_dense(&tmp2, &c);
}    

double lmer::LMMdev3() {
    update_prss();
    if (REML) {
	double nmp = (double)(n - p);
	return *ldL2 + *ldRX2 + nmp * (1 + log(2 * PI * (*prss)/nmp));
    }				       
    return *ldL2 + n * (1 + log(2 * PI * (*prss)/((double)n)));
}

/**
 * Update to new value of theta and evaluate the deviance or REML
 * criterion.
 *
 * @param thnew pointer to an numeric vector theta
 *
 * @return deviance or REML criterion according to the value of REML
 */
double lmerdense::update_dev(SEXP thnew) {
    int info;
    update_Lambda_Ut(thnew);
    LMMdev1();
    
    CHM_DN PLZtX = solvePL(N_AS_CHM_DN(ZtX, q, p));
    dble_cpy(RZX, (double*)(PLZtX->x), q * p);
    M_cholmod_free_dense(&PLZtX, &c);
    *ldRX2 = 0;
				// downdate and factor XtX, solve for fixef
    if (p) {
	dble_cpy(RX, XtX, p * p);
	F77_CALL(dsyrk)("U", "T", &p, &q, &mone, RZX, &q, &one, RX, &p);
	dble_cpy(fixef, Xty, p);
	F77_CALL(dgemv)("T", &q, &p, &mone, RZX, &q, (double*)(cu->x),
			&i1, &one, fixef, &i1);
	F77_CALL(dposv)("U", &p, &i1, RX, &p, fixef, &p, &info);
	if (info)
	    error(_("Downdated X'X is not positive definite, %d."), info);
				// evaluate ldRX2
	for (int i = 0; i < p; i++) *ldRX2 += 2 * log(RX[i * (p + 1)]);
	F77_CALL(dgemv)("N", &q, &p, &mone, RZX, &q, fixef, &i1, &one,
			(double*)(cu->x), &i1);
    }
    LMMdev2();
    update_eta();
    return LMMdev3();
}

double lmersparse::update_dev(SEXP thnew) {
    update_Lambda_Ut(thnew);
    LMMdev1();
				// update RZX
    CHM_SP PLZtX = solvePL(ZtX);
    CHM_SP_copy_in_place(RZX, PLZtX);
    M_cholmod_free_sparse(&PLZtX, &c);
				// downdate and factor XtX
    CHM_SP RZXt = M_cholmod_transpose(RZX, 1/*values*/, &c);
    CHM_SP RZXtRZX = M_cholmod_aat(RZXt, (int*)NULL/*fset*/,
				   (size_t)0/*fsize*/,
				   1/*mode = numerical*/, &c);
    M_cholmod_free_sparse(&RZXt, &c);
    CHM_SP DD = M_cholmod_add(XtX, RZXtRZX, &one, &mone, 1/*values*/,
			      1/*sorted*/, &c);
    M_cholmod_free_sparse(&RZXtRZX, &c);
    int *Perm = new int[p];
    for (int i = 0; i < p; i++) Perm[i] = i;
    CHM_FR LL = M_cholmod_analyze_p(DD, Perm, (int*)NULL, (size_t) 0, &c);
    if (!M_cholmod_factorize(DD, LL, &c))
	error(_("Downdated X'X is not positive definite"));
    delete[] Perm;
    M_cholmod_free_sparse(&DD, &c);
    *ldRX2 = M_chm_factor_ldetL2(LL);
				// evaluate fixef
    CHM_DN XtymRZXtu = M_cholmod_copy_dense(N_AS_CHM_DN(Xty, p, 1), &c);
    M_cholmod_sdmult(RZX, 1/*transpose*/, &mone, &one, cu, XtymRZXtu, &c);
    CHM_DN cfixef = M_cholmod_solve(CHOLMOD_A, LL, XtymRZXtu, &c);
    M_cholmod_free_dense(&XtymRZXtu, &c);
    dble_cpy(fixef, (double*)cfixef->x, p);
    M_cholmod_sdmult(RZX, 0/*no transpose*/, &mone, &one, cfixef, cu, &c);
    M_cholmod_free_dense(&cfixef, &c);
				// evaluate and store sparse RX
    CHM_SP Lsp = M_cholmod_factor_to_sparse(LL, &c);
    M_cholmod_free_factor(&LL, &c);
    CHM_SP Lspt = M_cholmod_transpose(Lsp, 1/*values*/, &c);
    M_cholmod_free_sparse(&Lsp, &c);
    CHM_SP_copy_in_place(RX, Lspt);
    M_cholmod_free_sparse(&Lspt, &c);
    LMMdev2();
    update_eta();
    return LMMdev3();
}

// definitions of methods for the merenvtrms classes

merenvtrms::merenvtrms(SEXP rho) {
    initMer(rho);
    flist = findVarBound(rho, install("flist"));
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

    if (!Lambda && ntrm == nfac) { // simple scalar terms only
	int offset = 0;
	double *cc = new double[q];
	CHM_DN col = N_AS_CHM_DN(cc, q, 1);
	
	for (int i = 0; i < nfac; i++) {
	    int nli = nl[i];
	    SET_VECTOR_ELT(ans, i, alloc3DArray(REALSXP, 1, 1, nli));
	    double *ai = REAL(VECTOR_ELT(ans, i));
//FIXME: Check with Tim to see if this is better done with dense or
//sparse columns. 
	    for (int j = 0; j < nli; j++) {
		for (int jj = 0; jj < q; jj++) cc[jj] = 0;
		cc[offset] = Lambdax[offset] * scale;
		CHM_DN sol =
		    M_cholmod_solve(CHOLMOD_A, L, col, &c);
		ai[j] = cc[offset] * ((double*)sol->x)[offset];
		M_cholmod_free_dense(&sol, &c);
		offset++;
	    }
	}
	delete[] cc;
	UNPROTECT(1);
	return ans;
    }

    CHM_SP lam = Lambda;
    if (!Lambda) {		// expand the diagonal Lambda
	lam = M_cholmod_speye(q, q, CHOLMOD_REAL, &c);
	Memcpy((double*)lam->x, Lambdax, q);
    }
    
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
		M_cholmod_submatrix(lam, (int*)NULL, -1, cset, nct,
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

    if (!Lambda)
	M_cholmod_free_sparse(&lam, &c);
    UNPROTECT(1);
    return ans;
}

glmer::glmer(SEXP rho) {
    initMer(rho);
    mu = VAR_REAL_NULL(rho, install("mu"), n);
    muEta = VAR_REAL_NULL(rho, install("muEta"), n);
    var = VAR_REAL_NULL(rho, install("var"), n);
    sqrtrwt = VAR_REAL_NULL(rho, install("sqrtrwt"), n);
    family = findVarBound(rho, install("family"));
    SEXP nms = getAttrib(family, R_NamesSymbol);
    if (!isNewList(family) ||!isString(nms) ||
	LENGTH(nms) != LENGTH(family))
	error(_("Object \"%s\" must be a named list"), "family");
    fam = CHAR(STRING_ELT(getListElement(family, nms, "family"), 1));
    lnk = CHAR(STRING_ELT(getListElement(family, nms, "link"), 1));
    ll = &identitylnk;
    vv = &constvr;
    if (!strcmp(lnk, logitlnk.nm())) ll = &logitlnk;
    if (!strcmp(fam, mu1muvr.distnm())) vv = &mu1muvr;
    Rprintf("variance name: %s\n", vv->distnm());
    Rprintf("link name: %s\n", ll->nm());
}

glmerdense::glmerdense(SEXP rho) : glmer(rho), merdense(rho) {
}

glmersparse::glmersparse(SEXP rho) : glmer(rho), mersparse(rho) {
}

extern "C" {

    SEXP glmer_linkinv(SEXP rho) {
	if (asLogical(findVarBound(rho, install("sparseX"))))
	    glmersparse(rho).linkinv();
	else
	    glmerdense(rho).linkinv();
	return R_NilValue;
    }
	
/**
 * Evaluate the deviance or REML criterion
 *
 * @param rho pointer to an lmerenv environment
 * @param thnew pointer to an numeric vector theta
 *
 * @return deviance value
 */
    SEXP lmerenv_deviance(SEXP rho, SEXP thnew) {
	if (asLogical(findVarBound(rho, install("sparseX"))))
	    return ScalarReal(lmersparse(rho).update_dev(thnew));
	else
	    return ScalarReal(lmerdense(rho).update_dev(thnew));
    }

/**
 * Check validity of an merenv environment
 *
 * @param x pointer to an merenv environment
 * @return TRUE if successful, otherwise it throws an error
 */
    SEXP lmerenv_validate(SEXP rho) {
	if (asLogical(findVarBound(rho, install("sparseX"))))
	    return ScalarLogical(lmersparse(rho).validate());
	else
	    return ScalarLogical(lmerdense(rho).validate());
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

/**
 * Check validity of an merenvtrms environment
 *
 * @param rho pointer to an merenvtrms environment
 * @return TRUE if successful, otherwise it throws an error
 */
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

