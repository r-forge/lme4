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
#include <cmath>		 /* for sqrt, etc. */
#include <R_ext/Lapack.h>        /* for Lapack (dpotrf, etc.) and BLAS */
#include "slotdefs.hpp"

extern "C" {


/**
 * Return the index of the term associated with parameter index ind
 *
 * @param ind an index in [0, Gp[nt] - 1]
 * @param nt total number of terms
 * @param Gp group pointers, a vector of length nt+1 with Gp[0] = 0
 *
 * @return sum of squares
 */
static inline int Gp_grp(int ind, int nt, const int *Gp)
{
    for (int i = 0; i < nt; i++) if (ind < Gp[i + 1]) return i;
    error(_("invalid row index %d (max is %d)"), ind, Gp[nt]);
    return -1;                  /* -Wall */
}

/**
 * Determine the index in theta_S corresponding to index i in the u vector
 *
 * @param i index in u vector (0 <= i < q)
 * @param nt number of random effects terms in the model
 * @param Gp vector of group pointers into u (length nt + 1, Gp[0] == 0)
 * @param nlev vector of length nt giving the number of levels per term
 * @param spt vector of group pointers into theta_S (length nt + 1,
 *            spt[0] == 0)
 */

static inline int theta_S_ind(int i, int nt, int *Gp, int *nlev, int *spt)
{
	int trm = Gp_grp(i, nt, Gp);
	return (spt[trm] + (i - Gp[trm]) / nlev[trm]);
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

/**
 * Multiply A on the left by T'
 *
 * @param A sparse model matrix
 * @param Gp group pointers
 * @param nc number of columns per term
 * @param nlev number of levels per term
 * @param st ST arrays for each term
 * @param nt number of terms
 * 
 */
static void Tt_Zt(CHM_SP A, int *Gp, int *nc, int *nlev, double **st, int nt)
{
    int *ai = (int*)(A->i), *ap = (int *)(A->p);
    double *ax = (double*)(A->x), one[] = {1,0};

    for (int j = 0; j < ((int)(A->ncol)); j++) /* multiply column j by T' */
	for (int p = ap[j]; p < ap[j + 1];) {
	    int i = Gp_grp(ai[p], nt, Gp);
	    
	    if (nc[i] <= 1) p++;
	    else {
		int nr = p;	/* number of rows in `B' in dtrmm call */
		while ((ai[nr] - Gp[i]) < nlev[i]) nr++;
		nr -= p;	/* nr == 1 except in models with carry-over */
		F77_CALL(dtrmm)("R", "L", "N", "U", &nr, nc + i,
				one, st[i], nc + i, ax + p, &nr);
		p += (nr * nc[i]);
	    }
	}
}

/**
 * Evaluate the sparse model matrix A from the Zt and ST slots
 *
 * @param x an mer object
 */
static void update_A(SEXP x)
{
    CHM_SP A = A_SLOT(x), Zt = Zt_SLOT(x);
    int *Gp = Gp_SLOT(x), *ai = (int*)(A->i), *ap = (int*)(A->p),
	*zi = (int*)(Zt->i), *zp = (int*)(Zt->p), ncmax,
	nt = DIMS_SLOT(x)[nt_POS];
    int annz = ap[A->ncol], znnz = zp[Zt->ncol];
    int *nc = Alloca(nt, int), *nlev = Alloca(nt, int);
    double **st = Alloca(nt, double*), *ax = (double*)(A->x),
	*zx = (double*)(Zt->x);
    R_CheckStack();

    ncmax = ST_nc_nlev(GET_SLOT(x, lme4_STSym), Gp, st, nc, nlev);

    if (annz == znnz) { /* Copy Z' to A unless A has new nonzeros */
	Memcpy(ax, zx, znnz);
    } else { /* Only for nonlinear models with correlated random effects */
	AZERO(ax, annz); 	/* Initialize potential nonzeros to 0 */
	for (int j = 0; j < ((int)(A->ncol)); j++) { /* Iterate over columns */
	    int pa = ap[j];
	    for (int p = zp[j]; p < zp[j + 1]; p++) { /* nonzeros in Z' */
		while (ai[pa] < zi[p]) pa++;          /* matching pos in A */
		if (ai[pa] != zi[p])
		    error(_("nonconforming Zt and A structures, j = %d"), j);
		ax[pa] = zx[p];
	    }
	}
    }
				/* When T != I multiply A on the left by T' */
    if (ncmax > 1) Tt_Zt(A, Gp, nc, nlev, st, nt);
				/* Multiply A on the left by S */
    for (int p = 0; p < annz; p++) {
	int i = Gp_grp(ai[p], nt, Gp);
	ax[p] *= st[i][((ai[p] - Gp[i]) / nlev[i]) * (nc[i] + 1)];
    }
}


/**
 * Update the contents of the ranef slot in an mer object using the
 * current contents of the u and ST slots.
 *
 * b = T  %*% S %*% t(P) %*% u
 *
 * @param x an mer object
 */
static void update_ranef(SEXP x)
{
    int *Gp = Gp_SLOT(x), *dims = DIMS_SLOT(x), *perm = PERM_SLOT(x);
    int nt = dims[nt_POS], q = dims[q_POS];
    double *b = RANEF_SLOT(x), *u = U_SLOT(x), one[] = {1,0};
    int *nc = Alloca(nt, int), *nlev = Alloca(nt, int);
    double **st = Alloca(nt, double*);
    R_CheckStack(); 

    ST_nc_nlev(GET_SLOT(x, lme4_STSym), Gp, st, nc, nlev);
				/* inverse permutation */
    for (int i = 0; i < q; i++) b[perm[i]] = u[i];
    for (int i = 0; i < nt; i++) {
	for (int k = 0; k < nc[i]; k++) { /* multiply by \tilde{S}_i */
	    double dd = st[i][k * (nc[i] + 1)];
	    int base = Gp[i] + k * nlev[i];
	    for (int kk = 0; kk < nlev[i]; kk++) b[base + kk] *= dd;
	}
	if (nc[i] > 1) {	/* multiply by \tilde{T}_i */
	    F77_CALL(dtrmm)("R", "L", "T", "U", nlev + i, nc + i, one,
			    st[i], nc + i, b + Gp[i], nlev + i);
	}
    }
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
    int *asgn = INTEGER(getAttrib(flistP, install("assign")));
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

/**
 * Extract the parameters from ST list
 *
 * @param x an mer object
 * @param pars vector of the appropriate length
 *
 * @return pointer to the parameter vector
 */
static double *ST_getPars(SEXP x, double *pars)
{
    SEXP ST = GET_SLOT(x, lme4_STSym);
    int nt = LENGTH(ST), pos = 0;
    for (int i = 0; i < nt; i++) {
	SEXP STi = VECTOR_ELT(ST, i);
	double *st = REAL(STi);
	int nci = INTEGER(getAttrib(STi, R_DimSymbol))[0];
	int ncp1 = nci + 1;

	for (int j = 0; j < nci; j++)
	    pars[pos++] = st[j * ncp1];
	for (int j = 0; j < (nci - 1); j++)
	    for (int k = j + 1; k < nci; k++)
		pars[pos++] = st[k + j * nci];
    }
    return pars;
}

/**
 * Update the ST, A, L, u, fixef and deviance slots of an mer object.
 *
 * @param x an mer object
 * @param pars double vector of the appropriate length
 * @return updated deviance
 *
 */
static void
ST_setPars(SEXP x, const double *pars)
{
    int *Gp = Gp_SLOT(x), nt = DIMS_SLOT(x)[nt_POS], pos = 0;
    int *nc = Alloca(nt, int), *nlev = Alloca(nt, int);
    double **st = Alloca(nt, double*);
    R_CheckStack();

    ST_nc_nlev(GET_SLOT(x, lme4_STSym), Gp, st, nc, nlev);
				/* install the parameters in the ST slot */
    for (int i = 0; i < nt; i++) {
	int nci = nc[i], ncp1 = nc[i] + 1;
	double *sti = st[i];

	for (int j = 0; j < nci; j++)
	    sti[j * ncp1] = pars[pos++];
	for (int j = 0; j < (nci - 1); j++)
	    for (int k = j + 1; k < nci; k++)
		sti[k + j * nci] = pars[pos++];
    }
    update_A(x);
}

/**
 * Return a list of (upper) Cholesky factors from the ST list
 *
 * @param x an mer object
 *
 * @return a list of upper Cholesky factors
 */
SEXP mer_ST_chol(SEXP x)
{
    SEXP ans = PROTECT(duplicate(GET_SLOT(x, lme4_STSym)));
    int ncmax, nt = DIMS_SLOT(x)[nt_POS];
    int *nc = Alloca(nt, int), *nlev = Alloca(nt, int);
    double **st = Alloca(nt, double*);
    R_CheckStack();

    ncmax = ST_nc_nlev(ans, Gp_SLOT(x), st, nc, nlev);
    for (int k = 0; k < nt; k++) {
	if (nc[k] > 1) {	// nothing to do for nc[k] == 1 
	    int nck = nc[k], nckp1 = nc[k] + 1;
	    double *stk = st[k];

	    for (int j = 0; j < nck; j++) {
		double dd = stk[j * nckp1]; // diagonal el
		for (int i = j + 1; i < nck; i++) {
		    stk[j + i * nck] = dd * stk[i + j * nck];
		    stk[i + j * nck] = 0;
		}
	    }
	}
    }

    UNPROTECT(1);
    return ans;
}

/**
 * Extract the parameters from the ST slot of an mer object
 *
 * @param x an mer object
 *
 * @return pointer to a REAL vector
 */
SEXP mer_ST_getPars(SEXP x)
{
    SEXP ans = PROTECT(allocVector(REALSXP, DIMS_SLOT(x)[np_POS]));
    ST_getPars(x, REAL(ans));

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
SEXP mer_ST_initialize(SEXP ST, SEXP Gpp, SEXP Zt)
{
    int *Gp = INTEGER(Gpp),
	*Zdims = INTEGER(GET_SLOT(Zt, lme4_DimSym)),
	*zi = INTEGER(GET_SLOT(Zt, lme4_iSym)), nt = LENGTH(ST);
    int *nc = Alloca(nt, int), *nlev = Alloca(nt, int),
	nnz = INTEGER(GET_SLOT(Zt, lme4_pSym))[Zdims[1]];
    double *rowsqr = Alloca(Zdims[0], double),
	**st = Alloca(nt, double*),
	*zx = REAL(GET_SLOT(Zt, lme4_xSym));
    R_CheckStack();
    
    ST_nc_nlev(ST, Gp, st, nc, nlev);
    AZERO(rowsqr, Zdims[0]);
    for (int i = 0; i < nnz; i++) rowsqr[zi[i]] += zx[i] * zx[i];
    for (int i = 0; i < nt; i++) {
	AZERO(st[i], nc[i] * nc[i]);
	for (int j = 0; j < nc[i]; j++) {
	    double *stij = st[i] + j * (nc[i] + 1);
	    for (int k = 0; k < nlev[i]; k++)
		*stij += rowsqr[Gp[i] + j * nlev[i] + k];
	    *stij = sqrt(nlev[i]/(0.375 * *stij));
	}
    }
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
SEXP mer_ST_setPars(SEXP x, SEXP pars)
{
  SEXP rpars = PROTECT(coerceVector(pars, REALSXP));
  int npar = DIMS_SLOT(x)[np_POS];

    if (LENGTH(pars) != npar)
	error(_("pars must be a real vector of length %d"), npar);
    ST_setPars(x, REAL(rpars));
    UNPROTECT(1);
    return R_NilValue;
}

/**
 * Externally callable update_ranef.
 * Update the contents of the ranef slot in an mer object.  For a
 * linear mixed model the conditional estimates of the fixed effects
 * and the conditional mode of u are evaluated first.
 *
 * @param x an mer object
 *
 * @return R_NilValue
 */
SEXP mer_update_ranef(SEXP x)
{
    update_ranef(x);
    return R_NilValue;
}

}
