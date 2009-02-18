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

// For improved convergence we return and set the parameters
// determining S on the square root scale.  That is, the diagonal
// element of S is the square of the corresponding parameter.

// The representation is defined by ST, a list of matrices, and Gp,
// the group pointers.

#include "ST.h"
#include "lme4utils.hpp"

class STinternal {
public:
    STinternal(SEXP x);		  //< external, R-level object
    ~STinternal() {delete[] nlev; delete[] nc; delete[] st;}

    void initialize(SEXP Zt);
    int Gp_grp(int ind);
    void bounds(double *lower, double *upper);
    double *Sdiag(double *d);
    CHM_SP Tmatrix();
    CHM_SP Lambda();
    CHM_SP create_A(CHM_SP Zt);
    void update_A(CHM_SP Zt, CHM_SP A);
    double *getPars(double *pars);
    int npars() {return np;}
    void chol(SEXP ans);
    SEXP create_ranef(SEXP u, SEXP perm);
    SEXP condVar(CHM_FR L, SEXP pperm, SEXP flistP, SEXP which);
    void setPars(const double *pars);
/**
 * Validate method.  This part is a no-op because validation is
 * done in the constructor.
 */
    SEXP validate() {return ScalarLogical(1);}

private:
    double **st;
    int *Gp, *nc, *nlev, nt, maxnc, np;
};

STinternal::STinternal(SEXP x)
{
    SEXP ST = GET_SLOT(x, lme4_STSym), Gpp = GET_SLOT(x, lme4_GpSym);
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

/**
 * Assign values to the diagonal of S
 *
 * @param d pointer to a vector of Gp[nt] values
 * @return d
 */
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

/**
 * Create the T matrix as a CHM_SP object
 *
 * @return T as a CHM_SP object
 */
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

/**
 * Create the Lambda matrix as a CHM_SP object
 *
 * @return Lambda as a CHM_SP object
 */
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

/**
 * Create A from Zt
 *
 * @return A as a CHM_SP object
 */
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
	    pars[pos++] = sqrt(st[i][j * ncp1]);
	for (int j = 0; j < (nci - 1); j++)
	    for (int k = j + 1; k < nci; k++)
		pars[pos++] = st[i][k + j * nci];
    }
    return pars;
}

/**
 * Install new parameters in the ST slot.
 *
 * @param pars double vector of the appropriate length
 *
 */
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

	for (int j = 0; j < nci; j++) {
	    sti[j * ncp1] = pars[pos] * pars[pos];
	    pos++;
	}
	for (int j = 0; j < (nci - 1); j++)
	    for (int k = j + 1; k < nci; k++)
		sti[k + j * nci] = pars[pos++];
    }
}

/**
 * Initialize the parameters in the ST slot from the Zt matrix
 *
 * @param Zt sparse transposed random effects model matrix
 *
 */
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

/**
 * Fill numeric vectors lower and upper, of length np with parameter bounds
 *
 * @param lower pointer to an numeric vector of length np
 * @param upper pointer to an numeric vector of length np
 */
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

/**
 * Utility that returns the index of the term corresponding to a row
 * or column index in Lambda.
 *
 * @param ind index in Lambda - must be in the range [0, Gp[nt]]
 */
int STinternal::Gp_grp(int ind)
{
    if (ind < 0 || ind >= Gp[nt])
	error(_("Invalid negative row index, %d, not in [0, %d]"),
	      ind, Gp[nt]);
    for (int i = 0; i < nt; i++)
	if (ind < Gp[i + 1]) return i;
    return -1;                  /* -Wall */
}


/**
 * Update A from Zt
 *
 * @param Zt original model matrix
 * @param A scaled model matrix
 */
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
 * @param uu pointer to the spherical random effects
 * @param pperm pointer to the 0-based permutation vector
 */
SEXP STinternal::create_ranef(SEXP uu, SEXP pperm)
{
    const int q = Gp[nt];
    if (!isReal(uu) || !isInteger(pperm) || LENGTH(uu) != q ||
	LENGTH(pperm) != q)
	error(_("u must be numeric and perm integer, both of length %d"),
	      q);
    double *b = new double[q], *u = REAL(uu);
    const double d1 = 1.;
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

/**
 * Evaluate the conditional variances, up to the common scale parameter.
 *
 * b = T  %*% S %*% t(P) %*% u
 *
 * @param L current Cholesky factor
 * @param flistP pointer to the factor list
 * @param which pointer to a logical vector
 * @param pperm pointer to the 0-based permutation vector
 */
SEXP STinternal::condVar(CHM_FR L, SEXP pperm, SEXP flistP, SEXP which)
{
    SEXP ans;
    const int nf = LENGTH(flistP), q = Gp[nt];
    if (!isInteger(pperm) || LENGTH(pperm) != q)
	error(_("perm must be an integer vector of length %d"), q);
    int nr = 0, pos = 0;
    int *asgn = INTEGER(getAttrib(flistP, install("assign")));
    double *vv, one[] = {1,0};
    CHM_SP sm1, sm2;
    CHM_DN dm1;
    int *perm = INTEGER(pperm), *iperm = new int[q];
	
    if (nt != nf) {
	for (int i = 0; i < nt; i++) {
	    if (asgn[i] < 1 || asgn[i] > nf)
		error(_("asgn[%d] is not in [1,%d]"), i + 1, nf);
	}
	// FIXME: Write the code for nt != nf
	error(_("Code for more terms than factors not yet written"));
    }
				// determine length of list to return 
    if (!isLogical(which) || LENGTH(which) != nf)
	error(_("which must be a logical vector of length %d"), nf);
    int *ww = LOGICAL(which);
    for (int i = 0; i < nt; i++) if (ww[i]) nr++;
    if (!nr) return(allocVector(VECSXP, 0));
    ans = PROTECT(allocVector(VECSXP, nr));
	
    for (int j = 0; j < q; j++) iperm[perm[j]] = j; // inverse permutation
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
		((double *)(rhs->x))[j] = st[i][j * (nci + 1)];
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
 * @param rho an evironment with Zt matrix
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
 * @param ST an ST object
 * @param rho an evironment that contains Zt and A
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
 * @param x an ST object
 * @param rho environment that contains Zt
 *
 * @return R_NilValue
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
 * @param rho environment that contains Zt and A
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
 * @param perm the permutation vector
 *
 * @return a list of matrices
 */
    SEXP ST_create_ranef(SEXP x, SEXP u, SEXP perm)
    {
	return STinternal(x).create_ranef(u, perm);
    }

/**
 * Evaluate the conditional variances of the random effects, up to the
 * common scale parameter. Some people called these posterior
 * variances, hence the name.
 *
 * @param x pointer to an mer object
 * @param L pointer to the current Cholesky factor
 * @param perm pointer to the permutaion vector
 * @param flist pointer to the factor list
 * @param which pointer to a logical vector
 *
 * @return pointer to a list of arrays
 */
    SEXP ST_postVar(SEXP x, SEXP L, SEXP perm, SEXP flist, SEXP which)
    {
	return STinternal(x).condVar(AS_CHM_FR(L), perm, flist, which);
    }

    SEXP ST_validate(SEXP x)
    {
	return STinternal(x).validate();
    }

}

