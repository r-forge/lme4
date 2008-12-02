// The ST representation of the variance factor in mixed models.

// The variance factor in a mixed model is a q by q matrix
// Lambda(theta) such that

// Var(B) = sigma^2 * Lambda(theta) %*% t(Lambda(theta))

// This particular type of variance factor is created from a q by q
// diagonal scale matrix S and a unit lower-triangular matrix T such
// that Lambda(theta) = T(theta) %*% S(theta) where T and S are
// generated in I blocks, of sizes n_i*q_i by n_i*q_i, i = 1,...,I

// The ith block of S is determined by q_i non-negative parameters,
// each repeated n_i times on the diagonal of S.  The ith block of T
// is determined by (q_i*(q_i - 1))/2 parameters, each creating a
// diagonal block of size n_i by n_i.  These blocks are filled in
// below the q_i diagonal blocks of size n_i by n_i.  Note that when
// q_i == 1 the corresponding block of T is an identity matrix of size
// n_i by n_i.

// The representation is defined by ST, a list of matrices, and Gp,
// the group pointers.

#include "ST.h"

void STinternal::assign_vals(SEXP ST, SEXP Gpp)
{
    nt = LENGTH(ST);
    if (!isInteger(Gpp) || !isNewList(ST) ||
	LENGTH(Gpp) != (nt + 1))
	error(_("Incompatible ST and Gp slots"));
    Gp = INTEGER(Gpp);
    if (Gp[0]) error(_("Gp[1] = %d != 0"), Gp[0]);
    st = new double*[nt];
    nc = new int[nt];
    nlev = new int[nt];
    maxnc = -1;
    np = 0;
    for (int i = 0; i < nt; i++) {
	SEXP STi = VECTOR_ELT(ST, i);
	int *dd = INTEGER(getAttrib(STi, R_DimSymbol));
	if (!(isReal(STi) && isMatrix(STi) && dd[0] && dd[0] == dd[1]))
	    error(_("ST[[%d]] is not a non-empty, square numeric matrix"), i + 1);
	int nci = dd[0];
	int Gpd = Gp[i + 1] - Gp[i];

	if (nci > maxnc) maxnc = nci;
	st[i] = REAL(STi);
	nc[i] = nci;
	if (Gpd <= 0 || Gpd % nci)
	    error(_("diff(Gp)[%d] is not a positive multiple of nc[%d]"),
		  i + 1, i + 1);
	nlev[i] = (Gp[i + 1] - Gp[i])/nci;
	np += (nci * (nci + 1))/2;
    }
}

double* STinternal::Sdiag(double *d)
{
    for (int i = 0, pos = 0; i < nt; i++)
	for (int c = 0; c < nc[i]; c++) {
	    double dd = st[i][c * (nc[i] + 1)];
	    for (int j = 0; j < nlev[i]; j++)
		d[pos++] = dd;
	}
    return d;
}

CHM_SP STinternal::Tmatrix()
{
    if (maxnc < 2)
	return M_cholmod_speye((size_t)Gp[nt], (size_t)Gp[nt],
			       CHOLMOD_REAL, &c);
    int nnz = 0;
    for (int i = 0; i < nt; i++)
	nnz += nlev[i] * (nc[i] * (nc[i] + 1)) / 2;
    CHM_SP A = M_cholmod_allocate_sparse((size_t) Gp[nt], // nrow
					 (size_t) Gp[nt], // ncol
					 (size_t) nnz,
					 1, // sorted
					 1, // packed
					 0, // not symmetric
					 CHOLMOD_REAL, &c);
    double *ax = (double*)(A->x);
    int *ai = (int*)(A->i), *ap = (int*)(A->p);
    for (int k = 0, p = 0; k < nt; k++) { // kth term
	for (int j = 0; j < nc[k]; j++) { // jth col of st[k]
	    for (int jj = 0; jj < nlev[k]; jj++) {
		for (int i = j; i < nc[k]; i++) { // ith row
		    if (i == j) { // diagonal block
			int Tcol = Gp[k] + j * nlev[k] + jj;
			ap[Tcol + 1] = ap[Tcol] + nc[k] - j;
			ax[p] = 1.;
		    } else ax[p] = st[k][i + j * nc[k]];
		    ai[p++] = Gp[k] + i * nlev[k] + jj;
		}
	    }
	}
    }
    return A;
}

CHM_SP STinternal::Lambda()
{
    CHM_SP A = Tmatrix();
    double *ax = (double*)(A->x), *dd = Sdiag(new double[Gp[nt]]);
    int *ap = (int*)(A->p), nc = (int)(A->ncol);
    for (int j = 0; j < nc; j++) /* scale column j */
	for (int p = ap[j]; p < ap[j+1]; p++) ax[p] *= dd[j];
    delete[] dd;
    return A;
}

CHM_SP STinternal::create_A(CHM_SP Zt)
{
    if (((int)(Zt->nrow)) != Gp[nt])
	error(_("nrow(Zt) = %d != ncol(Lambda) = %d"),
	      Zt->nrow, Gp[nt]);
    if (maxnc < 2) { // scale the rows
	CHM_SP A = M_cholmod_copy_sparse(Zt, &c);
	double *ax = (double*)(A->x), *dd = Sdiag(new double[Gp[nt]]);
	int *ap = (int*)(A->p), *ai = (int*)(A->i);
	for (int j = 0; j < ((int)(A->ncol)); j++)
	    for (int p = ap[j]; p < ap[j+1]; p++)
		ax[p] *= dd[ai[p]];
	delete[] dd;
	return A;
    }
    CHM_SP Lmb = Lambda();
    CHM_SP Lmbt = M_cholmod_transpose(Lmb, 1, &c);
    M_cholmod_free_sparse(&Lmb, &c);
    CHM_SP A = M_cholmod_ssmult(Lmbt, Zt, 0, 1, 1, &c);
    M_cholmod_free_sparse(&Lmbt, &c);
    return A;
}

double* STinternal::getPars(double *pars)
{
    for (int i = 0, pos = 0; i < nt; i++) {
	int nci = nc[i], ncp1 = nc[i] + 1;
	for (int j = 0; j < nci; j++)
	    pars[pos++] = st[i][j * ncp1];
	for (int j = 0; j < (nci - 1); j++)
	    for (int k = j + 1; k < nci; k++)
		pars[pos++] = st[i][k + j * nci];
    }
    return pars;
}

void STinternal::setPars(const double *pars)
{
    double *lower = new double[np], *upper = new double[np];
    bounds(lower, upper);
    for (int i = 0; i < np; i++)
	if (pars[i] < lower[i] || pars[i] > upper[i])
	    error(_("pars[%d] = %g is not in [%g,%g]"),
		  i + 1, pars[i], lower[i], upper[i]);
    delete[] lower; delete[] upper;

    for (int i = 0, pos = 0; i < nt; i++) {
	int nci = nc[i], ncp1 = nc[i] + 1;
	double *sti = st[i];

	for (int j = 0; j < nci; j++)
	    sti[j * ncp1] = pars[pos++];
	for (int j = 0; j < (nci - 1); j++)
	    for (int k = j + 1; k < nci; k++)
		sti[k + j * nci] = pars[pos++];
    }
}

void STinternal::initialize(SEXP Zt)
{
    int *Zdims = INTEGER(GET_SLOT(Zt, lme4_DimSym)),
	*zi = INTEGER(GET_SLOT(Zt, lme4_iSym));
    int nnz = INTEGER(GET_SLOT(Zt, lme4_pSym))[Zdims[1]];
    double *rowsqr = Alloca(Zdims[0], double),
	*zx = REAL(GET_SLOT(Zt, lme4_xSym));
    R_CheckStack();
	
    AZERO(rowsqr, Zdims[0]);
    for (int i = 0; i < nnz; i++) 
	rowsqr[zi[i]] += zx[i] * zx[i];
    for (int i = 0; i < nt; i++) {
	AZERO(st[i], nc[i] * nc[i]);
	for (int j = 0; j < nc[i]; j++) {
	    double *stij = st[i] + j * (nc[i] + 1);
	    for (int k = 0; k < nlev[i]; k++)
		*stij += rowsqr[Gp[i] + j * nlev[i] + k];
	    *stij = sqrt(nlev[i]/(0.375 * *stij));
	}
    }
}

void STinternal::bounds(double *lower, double *upper)
{
    for (int i = 0; i < np; i++) {
	lower[i] = R_NegInf;
	upper[i] = R_PosInf;
    }
    for (int i = 0, pos = 0; i < nt; i++) { // low = 0 on els of theta_S
	int nci = nc[i];
	for (int j = 0; j < nci; j++) lower[pos + j] = 0.;
	pos += (nci * (nci + 1)) / 2;
    }
}

int STinternal::Gp_grp(int ind)
{
    if (ind < 0 || ind >= Gp[nt])
	error(_("Invalid negative row index, %d, not in [0, %d]"),
	      ind, Gp[nt]);
    for (int i = 0; i < nt; i++)
	if (ind < Gp[i + 1]) return i;
    return -1;                  /* -Wall */
}

void STinternal::update_A(CHM_SP Zt, CHM_SP A)
{
    int *ai = (int*)(A->i), *ap = (int*)(A->p),
	*zi = (int*)(Zt->i), *zp = (int*)(Zt->p);
    int annz = ap[A->ncol], znnz = zp[Zt->ncol], anc = (int)(A->ncol);
    double *ax = (double*)(A->x), *zx = (double*)(Zt->x), one[] = {1, 0};

    if (annz == znnz) {	     // Copy Z' to A unless A has new nonzeros
	Memcpy(ax, zx, znnz);
    } else { // Only for nonlinear models with correlated random effects
	AZERO(ax, annz); 	// Initialize potential nonzeros to 0
	for (int j = 0; j < anc; j++) { // Iterate over columns
	    int pa = ap[j];
	    for (int p = zp[j]; p < zp[j + 1]; p++) { // nonzeros in Z'
		while (ai[pa] < zi[p]) pa++;          // matching pos in A
		if (ai[pa] != zi[p])
		    error(_("nonconforming Zt and A structures, j = %d"),
			  j);
		ax[pa] = zx[p];
	    }
	}
    }
    if (maxnc > 1) { // When T != I multiply A on the left by T'
	for (int j = 0; j < anc; j++) // mult column j by T'
	    for (int p = ap[j]; p < ap[j + 1];) {
		int i = Gp_grp(ai[p]);
		if (nc[i] <= 1) p++;
		else {
		    int nr = p;	// number of rows in `B' in dtrmm call
		    while ((ai[nr] - Gp[i]) < nlev[i]) nr++;
		    nr -= p;	// nr == 1 except in models with carry-over
		    F77_CALL(dtrmm)("R", "L", "N", "U", &nr, nc + i,
				    one, st[i], nc + i, ax + p, &nr);
		    p += (nr * nc[i]);
		}
	    }
    }

    double *dd = Sdiag(new double[Gp[nt]]);
    for (int p = 0; p < annz; p++) ax[p] *= dd[ai[p]]; // scale rows by S
    delete[] dd;
}


/**
 * Create the ranef matrices from u and perm.
 *
 * b = T  %*% S %*% t(P) %*% u
 *
 * @param x an mer object
 */
SEXP STinternal::create_ranef(SEXP uu, SEXP pperm)
{
    int q = Gp[nt];
    if (!isReal(uu) || !isInteger(pperm) || LENGTH(uu) != q ||
	LENGTH(pperm) != q)
	error(_("u must be numeric and perm integer, both of length %d"),
	      q);
    double *b = new double[q], *u = REAL(uu), d1 = 1.;
    int *perm = INTEGER(pperm);

    for (int i = 0; i < q; i++) b[perm[i]] = u[i]; // inverse permutation
    for (int i = 0; i < nt; i++) {
	for (int k = 0; k < nc[i]; k++) { // multiply by \tilde{S}_i
	    double dd = st[i][k * (nc[i] + 1)];
	    int base = Gp[i] + k * nlev[i];
	    for (int kk = 0; kk < nlev[i]; kk++) b[base + kk] *= dd;
	}
	if (nc[i] > 1) {	// multiply by \tilde{T}_i
	    F77_CALL(dtrmm)("R", "L", "T", "U", nlev + i, nc + i, &d1,
			    st[i], nc + i, b + Gp[i], nlev + i);
	}
    }
    SEXP ans = PROTECT(allocVector(VECSXP, nt));
    for (int i = 0; i < nt; i++) {
	SET_VECTOR_ELT(ans, i, allocMatrix(REALSXP, nlev[i], nc[i]));
	Memcpy(REAL(VECTOR_ELT(ans, i)), b + Gp[i], nlev[i] * nc[i]);
    }
    delete[] b;
    UNPROTECT(1);
    return ans;
}

void STinternal::chol(SEXP ans)
{
    for (int k = 0; k < nt; k++) {
	if (nc[k] > 1) {	// nothing to do for nc[k] == 1
	    int nck = nc[k], nckp1 = nc[k] + 1;
	    double *ak = REAL(VECTOR_ELT(ans, k)), *stk = st[k];
	    
	    for (int j = 0; j < nck; j++) {
		double dd = stk[j * nckp1]; // diagonal el
		for (int i = j + 1; i < nck; i++) {
		    ak[j + i * nck] = dd * stk[i + j * nck];
		    ak[i + j * nck] = 0;
		}
	    }
	}
    }
}


/**
 * Populate the st, nc and nlev arrays.  Return the maximum element of nc.
 *
 * @param ST pointer to a list (length nt) of matrices
 * @param Gp group pointers (length nt + 1)
 * @param st length nt array of (double*) pointers to be filled with
 * pointers to the contents of the matrices in ST.  Not used if NULL.
 * @param nc length nt array to be filled with the number of columns
 * @param nlev length nt array to be filled with the number of
 *        levels of the grouping factor for each term
 * 
 * @return maximum element of nc
 */
static int			/* populate the st, nc and nlev arrays */
ST_nc_nlev(const SEXP ST, const int *Gp, double **st, int *nc, int *nlev)
{
    int ans = 0, nt = LENGTH(ST);

    for (int i = 0; i < nt; i++) {
	SEXP STi = VECTOR_ELT(ST, i);
	int nci = *INTEGER(getAttrib(STi, R_DimSymbol));

	if (nci > ans) ans = nci;
	if (st) st[i] = REAL(STi);
	nc[i] = nci;
	nlev[i] = (Gp[i + 1] - Gp[i])/nci;
    }
    return ans;
}

extern "C" {

/**
 * Return a list of (upper) Cholesky factors from the ST list
 *
 * @param x an ST object
 *
 * @return a list of upper Cholesky factors
 */
SEXP ST_chol(SEXP x)
{
    SEXP ans = PROTECT(duplicate(GET_SLOT(x, lme4_STSym)));
    STinternal(x).chol(ans);
    UNPROTECT(1);
    return ans;
}

/**
 * Generate the A matrix from Zt.
 *
 * @param x an ST object
 * @param Zt the dgCMatrix representation of Z'
 *
 * @return the A matrix as an dgCMatrix object
 */
    SEXP ST_create_A(SEXP x, SEXP rho)
    {
	CHM_SP A = STinternal(x).create_A(AS_CHM_SP(findVarInFrame(rho, lme4_ZtSym)));
	SEXP ans = CHM_SP2SEXP(A, "dgCMatrix");
	M_cholmod_free_sparse(&A, &c);
	return ans;
    }

/**
 * Return the T matrix as a dtCMatrix object
 *
 * @param x an ST object
 
 * @return the T matrix as an dgCMatrix object
 */
    SEXP ST_Tmatrix(SEXP x)
    {
	CHM_SP A = STinternal(x).Tmatrix();
	SEXP ans = CHM_SP2SEXP(A, "dtCMatrix", "L", "N");
	M_cholmod_free_sparse(&A, &c);
	return ans;
    }

/**
 * Return the Lambda matrix as a dgCMatrix object
 *
 * @param x an ST object
 
 * @return the Lambda matrix as an dgCMatrix object
 */
    SEXP ST_Lambda(SEXP x)
    {
	CHM_SP A = STinternal(x).Lambda();
	SEXP ans = CHM_SP2SEXP(A, "dtCMatrix", "L", "N");
	M_cholmod_free_sparse(&A, &c);
	return ans;
    }
    

/**
 * Return the bounds on the parameter vector
 *
 * @param x an ST object
 *
 * @return numeric matrix with 2 columns and np rows
 */
    SEXP ST_bounds(SEXP x)
    {
	STinternal ST = STinternal(x);
	int np = ST.npars();
	SEXP ans = allocMatrix(REALSXP, np, 2);
	double *low = REAL(ans);
	ST.bounds(low, low + np);
	return ans;
    }
    

/**
 * Update the A array from the Zt array
 *
 * @param x an ST object
 *
 * @return numeric vector
 */
    SEXP ST_update_A(SEXP ST, SEXP rho)
    {
	STinternal(ST).update_A(AS_CHM_SP(findVarInFrame(rho, lme4_ZtSym)),
				AS_CHM_SP(findVarInFrame(rho, lme4_ASym)));
	return R_NilValue;
    }

/**
 * Extract the parameters from the ST slot of an mer object
 *
 * @param x an mer object
 *
 * @return pointer to a REAL vector
 */
    SEXP ST_getPars(SEXP x)
    {
	STinternal ST(x);
	SEXP ans = PROTECT(allocVector(REALSXP, ST.npars()));
	ST.getPars(REAL(ans));
	UNPROTECT(1); 
	return ans;
    }


/**
 * Evaluate starting estimates for the elements of ST
 *
 * @param ST pointers to the nt ST factorizations of the diagonal
 *     elements of Sigma 
 * @param Gpp length nt+1 vector of group pointers for the rows of Zt
 * @param Zt transpose of Z matrix
 *
 */
    SEXP ST_initialize(SEXP x, SEXP rho)
    {
	STinternal(x).initialize(findVarInFrame(rho, lme4_ZtSym));
	return R_NilValue;
    }

/**
 * Update the ST slot of an mer object from a REAL vector of
 * parameters and update the sparse model matrix A
 *
 * @param x an mer object
 * @param pars a REAL vector of the appropriate length
 *
 * @return R_NilValue
 */
    SEXP ST_setPars(SEXP x, SEXP pars, SEXP rho)
    {
	SEXP rpars = PROTECT(coerceVector(pars, REALSXP));
	STinternal ST(x);
	int np = ST.npars();
	
	if (LENGTH(pars) != np)
	    error(_("pars must be a real vector of length %d"), np);
	ST.setPars(REAL(rpars));
	UNPROTECT(1);
	ST.update_A(AS_CHM_SP(findVarInFrame(rho, lme4_ZtSym)),
		    AS_CHM_SP(findVarInFrame(rho, lme4_ASym)));
	return R_NilValue;
    }

/**
 * Create the random effects in the original scale as a list of matrices
 *
 * @param x an ST object
 * @param u the vector of orthogonal random effects
 * @param u the permutation vector
 *
 * @return a list of matrices
 */
    SEXP ST_create_ranef(SEXP x, SEXP u, SEXP perm)
    {
	return STinternal(x).create_ranef(u, perm);
    }

    SEXP ST_validate(SEXP x)
    {
	return STinternal(x).validate();
    }
/**
 * Extract the conditional variances of the random effects in an mer
 * object.  Some people called these posterior variances, hence the name.
 *
 * @param x pointer to an mer object
 * @param which pointer to a logical vector
 *
 * @return pointer to a list of arrays
 */
    SEXP mer_postVar(SEXP x, SEXP which)
    {
	int *Gp = Gp_SLOT(x), *dims = DIMS_SLOT(x), *ww;
	SEXP ans, flistP = GET_SLOT(x, lme4_flistSym);
	const int nf = LENGTH(flistP), nt = dims[nt_POS], q = dims[q_POS];
	int nr = 0, pos = 0;
//	int *asgn = INTEGER(getAttrib(flistP, install("assign")));
	double *vv, one[] = {1,0}, sc;
	CHM_SP sm1, sm2;
	CHM_DN dm1;
	CHM_FR L = L_SLOT(x);
	int *Perm = (int*)(L->Perm), *iperm = Alloca(q, int),
	    *nc = Alloca(nt, int), *nlev = Alloca(nt, int);
	double **st = Alloca(nt, double*);
	R_CheckStack();
	
/* FIXME: Write the code for nt != nf */
	if (nt != nf) error(_("Code not written yet"));
	/* determine length of list to return */
	if (!isLogical(which) || LENGTH(which) != nf)
	    error(_("which must be a logical vector of length %d"), nf);
	ww = LOGICAL(which);
	for (int i = 0; i < nt; i++) if (ww[i]) nr++;
	if (!nr) return(allocVector(VECSXP, 0));
	ans = PROTECT(allocVector(VECSXP, nr));
	
	ST_nc_nlev(GET_SLOT(x, lme4_STSym), Gp, st, nc, nlev);
	for (int j = 0; j < q; j++) iperm[Perm[j]] = j; /* inverse permutation */
	sc = dims[useSc_POS] ?
	    (DEV_SLOT(x)[dims[isREML_POS] ? sigmaREML_POS : sigmaML_POS]) : 1;
	for (int i = 0; i < nt; i++) {
	    if (ww[i]) {
		const int nci = nc[i];
		const int ncisqr = nci * nci;
		CHM_SP rhs = M_cholmod_allocate_sparse(q, nci, nci,
						       1/*sorted*/, 1/*packed*/,
						       0/*stype*/, CHOLMOD_REAL, &c);
		
		SET_VECTOR_ELT(ans, pos, alloc3DArray(REALSXP, nci, nci, nlev[i]));
		vv = REAL(VECTOR_ELT(ans, pos));
		pos++;
		for (int j = 0; j <= nci; j++) ((int *)(rhs->p))[j] = j;
		for (int j = 0; j < nci; j++)
		    ((double *)(rhs->x))[j] = st[i][j * (nci + 1)] * sc;
		for (int k = 0; k < nlev[i]; k++) {
		    double *vvk = vv + k * ncisqr;
		    for (int j = 0; j < nci; j++)
			((int*)(rhs->i))[j] = iperm[Gp[i] + k + j * nlev[i]];
		    sm1 = M_cholmod_spsolve(CHOLMOD_L, L, rhs, &c);
		    sm2 = M_cholmod_transpose(sm1, 1 /*values*/, &c);
		    M_cholmod_free_sparse(&sm1, &c);
		    sm1 = M_cholmod_aat(sm2, (int*)NULL, (size_t)0, 1 /*mode*/, &c);
		    dm1 = M_cholmod_sparse_to_dense(sm1, &c);
		    M_cholmod_free_sparse(&sm1, &c); M_cholmod_free_sparse(&sm2, &c);
		    Memcpy(vvk, (double*)(dm1->x), ncisqr);
		    M_cholmod_free_dense(&dm1, &c);
		    if (nci > 1) {
			F77_CALL(dtrmm)("L", "L", "N", "U", nc + i, nc + i,
					one, st[i], nc + i, vvk, nc + i);
			F77_CALL(dtrmm)("R", "L", "T", "U", nc + i, nc + i,
					one, st[i], nc + i, vvk, nc + i);
		    }
		}
		M_cholmod_free_sparse(&rhs, &c);
	    }
	}
	UNPROTECT(1);
	return ans;
    }
    
}
