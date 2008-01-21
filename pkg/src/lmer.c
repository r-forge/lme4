#include "lmer.h"
#include <Rmath.h>		/* for dnorm5, etc. */
#include <R_ext/Lapack.h>     /* for Lapack (dpotrf, etc.) and BLAS */
#include <R_ext/stats_package.h> /* for S_nlminb_iterate */
#include "Matrix.h"		 /* for cholmod functions */

extern
#include "Syms.h"
extern	       /** cholmod_common struct initialized in R_init_lme4 */
cholmod_common c;

#ifdef ENABLE_NLS		/** i18n of error messages */
#include <libintl.h>
#define _(String) dgettext ("lme4", String)
#else
#define _(String) (String)
#endif

/* When appropriate, alloca is cleaner than malloc/free.  The storage
 * is freed automatically on return from a function. When using gcc the
 * builtin version is much faster. */
#ifdef __GNUC__
# undef alloca
# define alloca(x) __builtin_alloca((x))
#else
/* this is necessary (and sufficient) for Solaris 10: */
# ifdef __sun
#  include <alloca.h>
# endif
#endif

/** alloca n elements of type t */
#define Alloca(n, t)   (t *) alloca( (size_t) ( (n) * sizeof(t) ) )

/** zero an array */
#define AZERO(x, n) {int _I_, _SZ_ = (n); for(_I_ = 0; _I_ < _SZ_; _I_++) (x)[_I_] = 0;}

/** positions in the deviance vector */
enum devP {
    ML_POS=0,			/**<Maximum likelihood deviance  */
    REML_POS,			/**<REML deviance */
    ldL2_POS,			/**<2*log-determinant of L */
    ldRX2_POS,			/**<2*log-determinant of RX */
    pwrss_POS,			/**<penalized, weighted residual sum of squares */
    disc_POS,			/**<discrepancy */
    usqr_POS,			/**<squared length of u */
    wrss_POS			/**<weighted residual sum of squares  */
};
/** positions in the dims vector */
enum dimP {
    nf_POS=0,			/**<number of terms in random effects */
    n_POS,			/**<number of observations */
    p_POS,			/**<number of fixed-effects parameters */
    q_POS,			/**<number of random effects */
    s_POS,			/**<number of variables in h (1 unless nonlinear) */
    np_POS,			/**<total number of parameters for T and S */
    isREML_POS,			/**<indicator of REML estimation */
    fTyp_POS,			/**<family type for generalized model */
    lTyp_POS,			/**<link type for generalized model */
    nest_POS,			/**<indicator of nested grouping factors */
    useSc_POS,			/**<does the family use a separate scale parameter */
    cvg_POS			/**<convergence indictor from port optimization  */
};

/** Residual degrees of freedom */
#define RDF(dims) (dims[n_POS] - (dims[isREML_POS] ? dims[p_POS] : 0))

/** Extract the L slot, convert it to a cholmod_factor struct and return a pointer. */
#define L_SLOT(x) AS_CHM_FR(GET_SLOT(x, lme4_LSym))

/**
 * Return the sum of squares of the first n elements of x
 *
 * @param n
 * @param x
 *
 * @return sum of squares
 */
static R_INLINE double sqr_length(const double *x, int n)
{
    double ans = 0;
    int i;

    for (i = 0; i < n; i++) ans += x[i] * x[i];
    return ans;
}

/**
 * Permute the vector src according to perm into dest
 *
 * @param dest destination
 * @param src source
 * @param perm NULL or 0-based permutation of length n
 * @param n length of src, dest and perm
 *
 * @return dest
 *
 * \note If perm is NULL the first n elements of src are copied to dest.
 */
static R_INLINE double*
apply_perm(double *dest, const double *src, const int *perm, int n)
{
    for (int i = 0; i < n; i++) dest[i] = src[perm ? perm[i] : i];
    return dest;
}

/**
 * Permute the vector src according to the inverse of perm into dest
 *
 * @param dest destination
 * @param src source
 * @param perm NULL or 0-based permutation of length n
 * @param n length of src, dest and perm
 *
 * @return dest
 *
 * \note If perm is NULL the first n elements of src are copied to dest.
 */
static R_INLINE double*
apply_iperm(double *dest, const double *src, const int *perm, int n)
{
    for (int i = 0; i < n; i++) dest[perm ? perm[i] : i] = src[i];
    return dest;
}

/**
 * Return the group in the (nf, Gp) combination to which ind belongs
 *
 * @param ind a row number
 * @param nf number of groups of rows
 * @param Gp group pointers
 *
 */
static R_INLINE int Gp_grp(int ind, int nf, const int *Gp)
{
    for (int i = 0; i < nf; i++) if (ind < Gp[i + 1]) return i;
    error(_("invalid row index %d (max is %d)"), ind, Gp[nf]);
    return -1;			/* -Wall */
}

/**
 * Extract the slot named nm from the object obj and return a null pointer
 * if the slot has length 0 or a pointer to the REAL contents.
 *
 * @param obj pointer to an S4 object
 * @param nm pointer to a symbol naming the slot to extract
 * 
 * @return pointer to the REAL contents, if nonzero length, otherwise
 * a NULL pointer 
 *
 */
static R_INLINE double *SLOT_REAL_NULL(SEXP obj, SEXP nm)
{
    SEXP pt = GET_SLOT(obj, nm);
    return LENGTH(pt) ? REAL(pt) : (double*) NULL;
}

/* FIXME: Instead of this function add elements sigmaML and sigmaREML
 * to the deviance slot and fill them whenever d[pwrss_POS] changes.
 */
/**
 * Return the REML or ML conditional estimate of sigma, the standard
 * deviation of the per-observation noise term.
 *
 * @param REML non-zero for REML estimate, 0 for ML estimate
 * @param dims vector of dimensions
 * @param d vector of deviance components
 */
static R_INLINE double
get_sigma(int REML, const int* dims, const double* d)
{
    return (dims[useSc_POS] ? 
	    sqrt(d[pwrss_POS] / (dims[n_POS] -
			    (REML ? dims[p_POS] : 0))) : 1);
}

/**
 * Populate the st, nc and nlev arrays.  Return the maximum element of nc.
 *
 * @param ST pointer to a list (length nf) of matrices
 * @param Gp group pointers (length nf + 1)
 * @param st length nf array of (double*) pointers to be filled with
 * pointers to the contents of the matrices in ST.  Not used if NULL.
 * @param nc length nf array to be filled with the number of columns
 * @param nlev length nf array to be filled with the number of
 *        levels of the grouping factor for each term
 * 
 * @return maximum element of nc
 */
static int			/* populate the st, nc and nlev arrays */
ST_nc_nlev(const SEXP ST, const int *Gp, double **st, int *nc, int *nlev)
{
    int ans = 0, i, nf = LENGTH(ST);

    for (i = 0; i < nf; i++) {
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
 * Determine the nonzero positions in the jth column of C
 *
 * @param nz array to hold the answer
 * @param j column index
 * @param nf number of groups of rows
 * @param Gp group pointers
 * @param nc array (length nf) of the number of columns
 * @param nlev array (length nf) of the number of levels
 * @param zi row indices in Zt
 * @param zp column pointers for Zt
 *
 * @return count of nonzeros in the jth column of C
 */
static int
C_nz_col(int *nz, int j, int nf, const int *Gp, const int *nc,
       const int *nlev, const int *zi, const int *zp)
{
    int ans, i, p;

    AZERO(nz, Gp[nf]);
    for (p = zp[j]; p < zp[j + 1]; p++) {
	int zrow = zi[p];
	int k = Gp_grp(zrow, nf, Gp);

	nz[zrow] = 1;		/* T contains the identity */
	if (nc[k] > 1) {
	    int nextra = (zrow - Gp[k]) / nlev[k];
	    for (i = 1; i <= nextra; i++)
		nz[zrow - i * nlev[k]] = 1;
	}
    }
    for (i = 0, ans = 0; i < Gp[nf]; i++) if (nz[i]) ans++;
    return ans;
}

/**
 * Update the contents of the ranef slot in an mer object.  For a
 * linear mixed model the conditional estimates of the fixed effects
 * and the conditional mode of u are evaluated first.
 *
 * b = T  %*% S %*% t(P) %*% u
 *
 * @param x an mer object
 */
static void update_ranef(SEXP x)
{
    SEXP ST = GET_SLOT(x, lme4_STSym);
    int *Gp = INTEGER(GET_SLOT(x, lme4_GpSym)),
	*dims = INTEGER(GET_SLOT(x, lme4_dimsSym)),
	i, i1 = 1, k, nf = LENGTH(ST),
	p = dims[p_POS], q = dims[q_POS];
    double *b = REAL(GET_SLOT(x, lme4_ranefSym)),
	*u = REAL(GET_SLOT(x, lme4_uSym)), one = 1;
    CHM_FR L = L_SLOT(x);
    int *nc = Alloca(nf, int), *nlev = Alloca(nf, int);
    double **st = Alloca(nf, double*);
    R_CheckStack(); 

    if (!SLOT_REAL_NULL(x, lme4_sqrtXWtSym)) { /* LMM */
	double *rx = REAL(GET_SLOT(x, lme4_RXySym)),
	    *rzx = REAL(GET_SLOT(x, lme4_RCXySym)),
	    *d = REAL(GET_SLOT(x, lme4_devianceSym)),
	    *fixef = REAL(GET_SLOT(x, lme4_fixefSym)),
	    mone = -1, one = 1;
	CHM_DN cu = N_AS_CHM_DN(u, q, 1), sol;
	R_CheckStack();
				/* solve RX beta-hat = del2 */
	F77_CALL(dtrsv)("U", "N", "N", &p, rx, &p, fixef, &i1);
	F77_CALL(dgemv)("N", &q, &p, &mone, rzx, &q,
			fixef, &i1, &one, u, &i1);
	if (!(sol = M_cholmod_solve(CHOLMOD_Lt, L, cu, &c)))
	    error(_("cholmod_solve (CHOLMOD_Lt) failed:"));
	Memcpy(u, (double*)(sol->x), q);
	M_cholmod_free_dense(&sol, &c);
	d[usqr_POS] = sqr_length(u, q);
	d[wrss_POS] = d[disc_POS] = d[pwrss_POS] - d[usqr_POS];
    }
    ST_nc_nlev(ST, Gp, st, nc, nlev);
				/* inverse permutation */
    for (i = 0; i < q; i++) b[((int*)L->Perm)[i]] = u[i];
    for (i = 0; i < nf; i++) {
	for (k = 0; k < nc[i]; k++) { /* multiply by \tilde{S}_i */
	    double dd = st[i][k * (nc[i] + 1)];
	    int base = Gp[i] + k * nlev[i], kk;
	    for (kk = 0; kk < nlev[i]; kk++) b[base + kk] *= dd;
	}
	if (nc[i] > 1) {	/* multiply by \tilde{T}_i */
	    F77_CALL(dtrmm)("R", "L", "T", "U", nlev + i, nc + i, &one,
			    st[i], nc + i, b + Gp[i], nlev + i);
	}
    }
}

/**
 * Evaluate starting estimates for the elements of ST
 *
 * @param ST pointers to the nf ST factorizations of the diagonal
 *     elements of Sigma 
 * @param Gpp length nf+1 vector of group pointers for the rows of Zt
 * @param Zt transpose of Z matrix
 *
 */
SEXP mer_ST_initialize(SEXP ST, SEXP Gpp, SEXP Zt)
{
    int *Gp = INTEGER(Gpp),
	*Zdims = INTEGER(GET_SLOT(Zt, lme4_DimSym)),
	*zi = INTEGER(GET_SLOT(Zt, lme4_iSym)),
	i, j, k, nf = LENGTH(ST);
    int *nc = Alloca(nf, int), *nlev = Alloca(nf, int),
	nnz = INTEGER(GET_SLOT(Zt, lme4_pSym))[Zdims[1]];
    double *rowsqr = Alloca(Zdims[0], double),
	**st = Alloca(nf, double*),
	*zx = REAL(GET_SLOT(Zt, lme4_xSym));
    R_CheckStack();
    
    ST_nc_nlev(ST, Gp, st, nc, nlev);
    AZERO(rowsqr, Zdims[0]);
    for (i = 0; i < nnz; i++) rowsqr[zi[i]] += zx[i] * zx[i];
    for (i = 0; i < nf; i++) {
	AZERO(st[i], nc[i] * nc[i]);
	for (j = 0; j < nc[i]; j++) {
	    double *stij = st[i] + j * (nc[i] + 1);
	    for (k = 0; k < nlev[i]; k++)
		*stij += rowsqr[Gp[i] + j * nlev[i] + k];
	    *stij = sqrt(nlev[i]/(0.375 * *stij));
	}
    }
    return R_NilValue;
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
    int i, nf = LENGTH(ST), pos = 0;
    for (i = 0; i < nf; i++) {
	SEXP STi = VECTOR_ELT(ST, i);
	double *st = REAL(STi);
	int nci = INTEGER(getAttrib(STi, R_DimSymbol))[0];
	int j, k, ncp1 = nci + 1;

	for (j = 0; j < nci; j++)
	    pars[pos++] = st[j * ncp1];
	for (j = 0; j < (nci - 1); j++)
	    for (k = j + 1; k < nci; k++)
		pars[pos++] = st[k + j * nci];
    }
    return pars;
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
    SEXP ans =
	PROTECT(allocVector(REALSXP, 
			    INTEGER(GET_SLOT(x, lme4_dimsSym))[np_POS]));
    ST_getPars(x, REAL(ans));
    UNPROTECT(1); 
    return ans;
}

/**
 * Update the ST and C slots of an mer object.
 *
 * @param x an mer object
 * @param pars double vector of the appropriate length
 *
 */
static void
ST_setPars(SEXP x, const double *pars)
{
    SEXP ST = GET_SLOT(x, lme4_STSym);
    int *Gp = INTEGER(GET_SLOT(x, lme4_GpSym));
    CHM_SP cC = AS_CHM_SP(GET_SLOT(x, lme4_ASym)),
	cZt = AS_CHM_SP(GET_SLOT(x, lme4_ZtSym));
    int *vi = (int*)(cC->i), *vp = (int*)(cC->p),
	*zi = (int*)(cZt->i), *zp = (int*)(cZt->p),
	i, j, ncmax, nf = LENGTH(ST), p, pos = 0;
    int vnnz = vp[cC->ncol], znnz = zp[cZt->ncol];
    int *nc = Alloca(nf, int), *nlev = Alloca(nf, int);
    double **st = Alloca(nf, double*),
	*vx = (double*)(cC->x), *zx = (double*)(cZt->x),
	one[] = {1,0};
    R_CheckStack();

    ncmax = ST_nc_nlev(ST, Gp, st, nc, nlev);
				/* install the parameters in the ST slot */
    for (i = 0; i < nf; i++) {
	int j, k, nci = nc[i], ncp1 = nc[i] + 1;
	double *sti = st[i];

	for (j = 0; j < nci; j++)
	    sti[j * ncp1] = pars[pos++];
	for (j = 0; j < (nci - 1); j++)
	    for (k = j + 1; k < nci; k++)
		sti[k + j * nci] = pars[pos++];
    }
				/* Copy Z' to V' unless V' has nonzero
				 * positions not in Z' */ 
    if (vnnz == znnz) Memcpy(vx, zx, znnz); 
    else { /* Only occurs in nonlinear models with correlated random effects */
	AZERO(vx, vnnz); 	/* Initialize the potential nonzeros to 0 */
	for (j = 0; j < cC->ncol; j++) { /* Iterate over columns */
	    int pv = vp[j];
	    for (p = zp[j]; p < zp[j + 1]; p++) { /* nonzeros in Z' */
		while (vi[pv] < zi[p]) pv++;	  /* matching pos in V' */
		if (vi[pv] != zi[p])
		    error(_("nonconforming Zt and C structures, j = %d"), j);
		vx[pv] = zx[p];
	    }
	}
    }
				/* When T != I multiply C on the left by T' */
    if (ncmax > 1)		/* T' == I when ncmax == 1 */
	for (j = 0; j < cC->ncol; j++) /* multiply column j by T' */
	    for (p = vp[j]; p < vp[j + 1];) {
		int i = Gp_grp(vi[p], nf, Gp);

		if (nc[i] <= 1) p++;
		else {
		    int nr = p;	/* number of rows in `B' in dtrmm call */
		    while ((vi[nr] - Gp[i]) < nlev[i]) nr++;
		    nr -= p;	/* nr == 1 except in models with carry-over */
		    F77_CALL(dtrmm)("R", "L", "N", "U", &nr, nc + i,
				    one, st[i], nc + i, vx + p, &nr);
		    p += (nr * nc[i]);
		}
	    }
				/* Multiply C on the left by S */
    for (p = 0; p < vnnz; p++) {
	int i = Gp_grp(vi[p], nf, Gp);
	vx[p] *= st[i][((vi[p] - Gp[i]) / nlev[i]) * (nc[i] + 1)];
    }
}

/**
 * Update the ST slot of an mer object from a REAL vector of
 * parameters and update the Cholesky factorization
 *
 * @param x an mer object
 * @param pars a REAL vector of the appropriate length
 *
 * @return R_NilValue
 */
SEXP mer_ST_setPars(SEXP x, SEXP pars)
{
    int npar = INTEGER(GET_SLOT(x, lme4_dimsSym))[np_POS];

    if (!isReal(pars) || LENGTH(pars) != npar)
	error(_("pars must be a real vector of length %d"), npar);
    ST_setPars(x, REAL(pars));
    return R_NilValue;
}

/**
 * Extract the estimate of the common scale parameter from an mer object
 *
 * @param x an lmer object
 * @param which scalar integer (< 0 => REML, 0 => native, > 0 => ML)
 *
 * @return scalar REAL value
 */
SEXP mer_sigma(SEXP x, SEXP which)
{
    int *dims = INTEGER(GET_SLOT(x, lme4_dimsSym)),
	w = LENGTH(which) ? asInteger(which) : 0;
    
    return ScalarReal(get_sigma(w < 0 || (!w && dims[isREML_POS]), dims,
				REAL(GET_SLOT(x, lme4_devianceSym))));
}

/**
 * Extract the posterior variances of the random effects in an mer object
 *
 * @param x pointer to an mer object
 *
 * @return pointer to a list of arrays
 */
SEXP mer_postVar(SEXP x)
{
    int *Gp = INTEGER(GET_SLOT(x, lme4_GpSym)),
	*dims = INTEGER(GET_SLOT(x, lme4_dimsSym));
    int nf = dims[nf_POS], q = dims[q_POS];
    SEXP ans = PROTECT(allocVector(VECSXP, nf));
    int i, j, k;
    double *vv, one[] = {1,0}, sc;
    CHM_SP sm1, sm2;
    CHM_DN dm1;
    CHM_FR L = L_SLOT(x);
    int *Perm = (int*)(L->Perm), *iperm = Alloca(q, int),
	*nc = Alloca(nf, int), *nlev = Alloca(nf, int);
    double **st = Alloca(nf, double*);
    R_CheckStack();

    ST_nc_nlev(GET_SLOT(x, lme4_STSym), Gp, st, nc, nlev);
    for (j = 0; j < q; j++) iperm[Perm[j]] = j; /* inverse permutation */
    sc = get_sigma(dims[isREML_POS], dims, REAL(GET_SLOT(x, lme4_devianceSym)));
    for (i = 0; i < nf; i++) {
	int ncisqr = nc[i] * nc[i];
	CHM_SP rhs = M_cholmod_allocate_sparse(q, nc[i], nc[i],
					       1/*sorted*/, 1/*packed*/,
					       0/*stype*/, CHOLMOD_REAL, &c);

	SET_VECTOR_ELT(ans, i, alloc3DArray(REALSXP, nc[i], nc[i], nlev[i]));
	vv = REAL(VECTOR_ELT(ans, i));
	for (j = 0; j <= nc[i]; j++) ((int *)(rhs->p))[j] = j;
	for (j = 0; j < nc[i]; j++)
	    ((double *)(rhs->x))[j] = st[i][j * (nc[i] + 1)] * sc;
	for (k = 0; k < nlev[i]; k++) {
	    double *vvk = vv + k * ncisqr;
	    for (j = 0; j < nc[i]; j++)
		((int*)(rhs->i))[j] = iperm[Gp[i] + k + j * nlev[i]];
	    sm1 = M_cholmod_spsolve(CHOLMOD_L, L, rhs, &c);
	    sm2 = M_cholmod_transpose(sm1, 1 /*values*/, &c);
	    M_cholmod_free_sparse(&sm1, &c);
	    sm1 = M_cholmod_aat(sm2, (int*)NULL, (size_t)0, 1 /*mode*/, &c);
	    dm1 = M_cholmod_sparse_to_dense(sm1, &c);
	    M_cholmod_free_sparse(&sm1, &c); M_cholmod_free_sparse(&sm2, &c);
	    Memcpy(vvk, (double*)(dm1->x), ncisqr);
	    M_cholmod_free_dense(&dm1, &c);
	    if (nc[i] > 1) {
		F77_CALL(dtrmm)("L", "L", "N", "U", nc + i, nc + i,
				one, st[i], nc + i, vvk, nc + i);
		F77_CALL(dtrmm)("R", "L", "T", "U", nc + i, nc + i,
				one, st[i], nc + i, vvk, nc + i);
	    }
	}
	M_cholmod_free_sparse(&rhs, &c);
    }
    UNPROTECT(1);
    return ans;
}

#if 0
static CHM_SP overlay_cols(CHM_SP A, int s)
{
    int *ai, *ap, *ti, *tj, j, n, nnz, p;
    double *ax, *tx;
    CHM_TR T;
    CHM_SP AA;
    
    if (s < 2) error(_("overlay_cols must have s > 1, s = %d"), s);
    if ((A->ncol) % s)
	error(_("A->ncol = %d must be a multiple of s = %d"), A->ncol, s);
    n = (A->ncol) / s; 
    ai = (int*)A->i; ap = (int*)A->p; ax = (double*)A->x;
    nnz = ap[A->ncol];

    T = M_cholmod_allocate_triplet(A->nrow, n, nnz,
				   0/*not symmetric*/, CHOLMOD_REAL, &c);
    ti = (int*)T->i; tj = (int*)T->j; tx = (double*)T->x;
    for (j = 0; j < n; j++) {
	for (p = ap[j]; p < ap[j + 1]; p++) {
	    ti[p] = ai[p];
	    tj[p] = j % n;
	    tx[p] = ax[p];
	}
    }
    AA = M_cholmod_triplet_to_sparse(T, nnz, &c);
    M_cholmod_free_triplet(&T, &c);
    return AA;
}
#endif

/**
 * Create and initialize L
 *
 * @param CmP pointer to the model matrix for the orthogonal random
 * effects (transposed)
 *
 * @return L
 */
SEXP mer_create_L(SEXP CmP)
{
    double one[] = {1, 0};
    CHM_SP Cm = AS_CHM_SP(CmP);
    CHM_FR L;
    R_CheckStack();

    L = M_cholmod_analyze(Cm, &c);
    if (!M_cholmod_factorize_p(Cm, one, (int*)NULL, 0 /*fsize*/, L, &c))
	error(_("cholmod_factorize_p failed: status %d, minor %d from ncol %d"),
	      c.status, L->minor, L->n);

    return M_chm_factor_to_SEXP(L, 1);
}


/* Utilities and constants for generalized linear models */

static const double LTHRESH = 30.;
static const double MLTHRESH = -30.;
static double MPTHRESH = 0;
static double PTHRESH = 0;
static const double INVEPS = 1/DOUBLE_EPS;


/**
 * Update the eta, v, mu, resid and var slots according to the current
 * values of the parameters and u.  Also evaluate d[wrss_POS] using
 * the current sqrtrWt slot.  The sqrtrWt slot is changed in update_L.
 *
 * @param x pointer to an mer object
 *
 * @return penalized, weighted residual sum of squares
 */
static double update_mu(SEXP x)
{
    SEXP v = GET_SLOT(x, lme4_vSym);
    int *dims = INTEGER(GET_SLOT(x, lme4_dimsSym));
    int i, i1 = 1, n = dims[n_POS], p = dims[p_POS], s = dims[s_POS];
    int ns = n * s;
    double *d = REAL(GET_SLOT(x, lme4_devianceSym)),
	*eta = REAL(GET_SLOT(x, lme4_etaSym)),
	*mu = REAL(GET_SLOT(x, lme4_muSym)),
	*muEta = SLOT_REAL_NULL(x, lme4_muEtaSym),
	*offset = SLOT_REAL_NULL(x, lme4_offsetSym),
	*srwt = SLOT_REAL_NULL(x, lme4_sqrtrWtSym),
	*res = REAL(GET_SLOT(x, lme4_residSym)),
	*var = SLOT_REAL_NULL(x, lme4_varSym),
	*y = REAL(GET_SLOT(x, lme4_ySym)),
	one[] = {1,0};
    CHM_FR L = L_SLOT(x);
    CHM_SP cC = AS_CHM_SP(GET_SLOT(x, lme4_ASym));
    CHM_DN Ptu, ceta = N_AS_CHM_DN(eta, ns, 1),
	cu = AS_CHM_DN(GET_SLOT(x, lme4_uSym));
    R_CheckStack();

				/* eta := offset or eta := 0 */
    for (i = 0; i < ns; i++) eta[i] = offset ? offset[i] : 0;
				/* eta := eta + X \beta */
    F77_CALL(dgemv)("N", &ns, &p, one,
		    REAL(GET_SLOT(x, lme4_XSym)), &ns,
		    REAL(GET_SLOT(x, lme4_fixefSym)), &i1,
		    one, eta, &i1);
				/* eta := eta + C' P' u */
    Ptu = M_cholmod_solve(CHOLMOD_Pt, L, cu, &c);
    if (!M_cholmod_sdmult(cC, 1 /* trans */, one, one, Ptu, ceta, &c))
	error(_("cholmod_sdmult error returned"));
    M_cholmod_free_dense(&Ptu, &c);

    if (LENGTH(v)) {		/* evaluate the nonlinear model */
	SEXP rho = GET_SLOT(x, lme4_envSym);
	SEXP gg, pnames = GET_SLOT(x, lme4_pnamesSym), vv;
	int *gdims, s = dims[s_POS];
	
	for (i = 0; i < s; i++) { /* par. vals. into env. */
	    vv = findVarInFrame(rho,
				install(CHAR(STRING_ELT(pnames, i))));
	    if (!isReal(vv) || LENGTH(vv) != n)
		error(_("Parameter %s in the environment must be a length %d numeric vector"),
		      CHAR(STRING_ELT(pnames, i)), n);
	    Memcpy(REAL(vv), eta + i * n, n);
	}
	vv = PROTECT(eval(GET_SLOT(x, lme4_nlmodelSym), rho));
	if (!isReal(vv) || LENGTH(vv) != n)
	    error(_("evaluated model is not a numeric vector of length %d"), n);
	gg = getAttrib(vv, lme4_gradientSym);
	if (!isReal(gg) || !isMatrix(gg))
	    error(_("gradient attribute of evaluated model must be a numeric matrix"));
	gdims = INTEGER(getAttrib(gg, R_DimSymbol));
	if (gdims[0] != n ||gdims[1] != s)
	    error(_("gradient matrix must be of size %d by %d"), n, s);
				/* colnames of the gradient
				 * corresponding to the order of the
				 * pnames has been checked */
	SET_SLOT(x, lme4_vSym, vv);
	UNPROTECT(1);
	eta = REAL(vv);
    }
				/* eta now points to an n vector  */
    if (muEta) {	/* apply the generalized linear part */
	switch(dims[fTyp_POS]) {
	case 1:			/* binomial with logit link */
	    for (i = 0; i < n; i++) {
		double etai = eta[i];
		double tmp = (etai < MLTHRESH) ? DOUBLE_EPS :
		    ((etai > LTHRESH) ? INVEPS : exp(etai));
		mu[i] = tmp/(1 + tmp);
		muEta[i] = var[i] = mu[i] * (1 - mu[i]);
	    }
	    break;
	case 2:			/* binomial with probit link */
	    if (!MPTHRESH) {
		MPTHRESH = qnorm5(DOUBLE_EPS, 0, 1, 1, 0);
		PTHRESH = -MPTHRESH;
	    }
	    for (i = 0; i < n; i++) {
		double etai = eta[i], tmp;
		mu[i] = (etai < MPTHRESH) ? DOUBLE_EPS :
		    ((etai > PTHRESH) ? 1 - DOUBLE_EPS :
		     pnorm5(etai, 0, 1, 1, 0));
		var[i] = mu[i] * (1 - mu[i]);
		tmp = dnorm4(eta[i], 0, 1, 0);
		muEta[i] =
		    (tmp < DOUBLE_EPS) ? DOUBLE_EPS : tmp;
	    }
	    break;
	case 3:			/* Poisson with log link */
	    for (i = 0; i < n; i++) {
		double tmp = exp(eta[i]);
		muEta[i] = var[i] = mu[i] =
		    (tmp < DOUBLE_EPS) ? DOUBLE_EPS : tmp;
	    }
	    break;
	default:
	    error(_("General form of glmer_linkinv not yet written"));
	} 
    } else {
	Memcpy(mu, eta, n);
    }

    d[wrss_POS] = 0;		/* update resid slot and d[wrss_POS] */
    for (i = 0; i < n; i++) {
	res[i] = (y[i] - mu[i]) * (srwt ? srwt[i] : 1);
	d[wrss_POS] += res[i] * res[i];
    }
				/* store u'u */
    d[usqr_POS] = sqr_length((double*)(cu->x), dims[q_POS]);
    return (d[pwrss_POS] = d[usqr_POS] + d[wrss_POS]);
}

/**
 * Externally callable update_mu.
 * Update the eta, v, mu, resid and var slots according to the current
 * values of the parameters and u.  Also evaluate d[wrss_POS] using
 * the current contents of sqrtrWt.  The sqrtrWt slot is updated in update_L.
 *
 * @param x pointer to an mer object
 *
 * @return penalized, weighted residual sum of squares
 */
SEXP mer_update_mu(SEXP x)
{
    return ScalarReal(update_mu(x));
}

/**
 * Update the A, L, sqrtrWt and resid slots.  It is assumed that
 * update_mu has already been called at the current values of u and
 * the model parameters.
 *
 * @param x pointer to an mer object
 *
 * @return penalized weighted residual sum of squares
 *
 */
static double update_L(SEXP x)
{
    SEXP vP = GET_SLOT(x, lme4_vSym);
    int *dims = INTEGER(GET_SLOT(x, lme4_dimsSym));
    int lv = LENGTH(vP), n = dims[n_POS], s = dims[s_POS];
    double
	*cx = SLOT_REAL_NULL(x, lme4_CxSym),
	*d = REAL(GET_SLOT(x, lme4_devianceSym)),
	*res = REAL(GET_SLOT(x, lme4_residSym)),
	*mu = REAL(GET_SLOT(x, lme4_muSym)),
	*muEta = SLOT_REAL_NULL(x, lme4_muEtaSym),
	*sXwt = SLOT_REAL_NULL(x, lme4_sqrtXWtSym),
	*srwt = SLOT_REAL_NULL(x, lme4_sqrtrWtSym),
	*var =  SLOT_REAL_NULL(x, lme4_varSym),
	*y = REAL(GET_SLOT(x, lme4_ySym)),
	one[] = {1,0};
    CHM_SP A = AS_CHM_SP(GET_SLOT(x, lme4_ASym));
    CHM_FR L = L_SLOT(x);
    R_CheckStack();
    
    if (var) {		    /* Update srwt and res. Reevaluate wrss. */
	double *pwt = SLOT_REAL_NULL(x, lme4_priorWtSym);
	
	d[wrss_POS] = 0;
	for (int j = 0; j < n; j++) {
	    srwt[j] = sqrt((pwt ? pwt[j] : 1) / var[j]);
	    res[j] = srwt[j] * (y[j] - mu[j]);
	    d[wrss_POS] += res[j] * res[j];
	}
    }
    if (sXwt) {			/* Update sXwt and C */
	int *ai = (int*)A->i, *ap = (int*)A->p, i, j;
	double *ax = (double*)(A->x), *grad = (double*) NULL;
	CHM_SP C = A;

	if (lv) grad = REAL(getAttrib(vP, lme4_gradientSym));
	for (j = 0; j < s; j++) {
	    for (i = 0; i < n; i++) {
		int ja = i + j * n;
		sXwt[ja] = (srwt ? srwt[i] : 1) *
		    (muEta ? muEta[i] : 1) *
		    (grad ? grad[ja] : 1);
	    }
	}
	if (s == 1) {		/* C is a scaled version of A */
	    for (j = 0; j < n; j++)
		for (int p = ap[j]; p < ap[j + 1]; p++)
		    cx[p] = ax[p] * sXwt[j];
	    A->x = (void*)cx;
	} else {
	    int *ci, *cp, pa, pc;
	    C = AS_CHM_SP(GET_SLOT(x, lme4_CmSym));
	    R_CheckStack();

	    ci = (int*)C->i; cp = (int*)C->p; cx = (double*)C->x;
	    for (j = 0; j < s; j++)
		for (i = 0; i < n; i++) {
		    int ja = i + j * n;
		    for (pa = ap[ja]; pa < ap[ja + 1]; pa++) {
			for (pc = cp[i]; pc < cp[i + 1]; pc++)
			    if (ci[pc] == ai[pa]) break;
			if (pc >= cp[i + 1])
			    error(_("Structure of Cm and A are not consistent"));
			cx[pc] += ax[pa] * sXwt[ja];
		    }
		}
	    A = C;
	}
    }

    if (!M_cholmod_factorize_p(A, one, (int*)NULL, 0 /*fsize*/, L, &c))
	error(_("cholmod_factorize_p failed: status %d, minor %d from ncol %d"),
	      c.status, L->minor, L->n);

    d[ldL2_POS] = 0;		/* evaluate log(det(L))^2 */
    if (L->is_super) {
	int i;
	for (i = 0; i < L->nsuper; i++) {
	    int j, nrp1 = 1 + ((int *)(L->pi))[i + 1] - ((int *)(L->pi))[i],
		nc = ((int *)(L->super))[i + 1] - ((int *)(L->super))[i];
	    double *x = (double *)(L->x) + ((int *)(L->px))[i];

	    for (j = 0; j < nc; j++) {
		d[ldL2_POS] += 2 * log(fabs(x[j * nrp1]));
	    }
	}
    } else {
	int *li = (int*)(L->i), *lp = (int*)(L->p), j, p;
	double *lx = (double *)(L->x);
	
	for (j = 0; j < L->n; j++) {
	    for (p = lp[j]; li[p] != j && p < lp[j + 1]; p++) {};
	    if (li[p] != j) break; /* what happened to the diagonal element? */
	    d[ldL2_POS] += log(lx[p] * ((L->is_ll) ? lx[p] : 1.));
	}
    }
    return (d[pwrss_POS] = d[usqr_POS] + d[wrss_POS]);
}

/**
 * Externally callable update_L.
 * Update the A, L, sqrtrWt and resid slots.  It is assumed that
 * update_mu has already been called at the current values of u and
 * the model parameters.
 *
 * @param x an mer object
 *
 * @return penalized weighted residual sum of squares
 *
 */
SEXP mer_update_L(SEXP x)
{
    return ScalarReal(update_L(x));
}

/** Maximum number of iterations in update_u */
#define CM_MAXITER  300
/** Tolerance level for convergence criterion in update_u */
#define CM_TOL      1e-10
/** Minimum step factor in update_u */
#define CM_SMIN     1e-5

#if 0
static void
PAHGv(double *dest, const int *perm, const CHM_SP A,
       const double *sXwt, double *v, int s)
{
    int *ai = (int*)(A->i), *ap = (int*)(A->p), j, k,
	ns = A->ncol, p, q = A->nrow;
    int n = ns/s;
    double *ax = (double*)(A->x),
	*td = Alloca(q, double),
	*ts = Alloca(ns, double);
    R_CheckStack();

    for (k = 0; k < s; k++)
	for (j = 0; j < n; j++)
	    ts[j + k * n] = sXwt[j + k * n] * v[j];
    AZERO(td, q);
    for (j = 0; j < ns; j++) {
	for (p = ap[j]; p < ap[j + 1]; p++)
	    td[ai[p]] += ts[j] * ax[p];
    }
    if (perm) apply_perm(dest, td, perm, q);
}
#endif
    
/**
 * Iterate to determine the conditional modes of the random effects.
 *
 * @param x pointer to an mer object
 * @param verb indicator of verbose iterations 
 *             (negative values produce a lot of output)
 *
 * @return number of iterations to convergence (0 for non-convergence) 
 */
static int update_u(SEXP x, int verb)
{
    int *dims = INTEGER(GET_SLOT(x, lme4_dimsSym));
    int i, j, n = dims[n_POS], q = dims[q_POS];
    double *Cx = SLOT_REAL_NULL(x, lme4_CxSym),
	*sXwt = SLOT_REAL_NULL(x, lme4_sqrtXWtSym),
	*res = REAL(GET_SLOT(x, lme4_residSym)), 
	*u = REAL(GET_SLOT(x, lme4_uSym)),
	cfac = ((double)n) / ((double)q), 
	crit, pwrss, pwrss_old, step;
    double *tmp = Alloca(q, double), *tmp1 = Alloca(q, double),
	*uold = Alloca(q, double), one[] = {1,0}, zero[] = {0,0};
    CHM_FR L = L_SLOT(x);
    CHM_DN cres = N_AS_CHM_DN(res, n, 1),
	ctmp = N_AS_CHM_DN(tmp, q, 1), sol;
    CHM_SP C = AS_CHM_SP(GET_SLOT(x, lme4_CmSym));
    R_CheckStack();
    
    if (!(L->is_ll)) error(_("L must be LL', not LDL'"));
    if (q > n) error(_("q = %d > n = %d"), q, n);
    if (!sXwt) error(_("sqrtXWt slot is empty in update_u"));
    if (Cx) {		    /* A and C have the same structure */
	C = AS_CHM_SP(GET_SLOT(x, lme4_ASym));
	R_CheckStack();
	C->x = (void*)Cx;
    }
    
    /* resetting u to zero at the beginning of each evaluation
     * requires more iterations but is necessary to obtain a
     * repeatable evaluation.  If this is not done the optimization
     * algorithm can take wild steps. */
    AZERO(u, q);
    update_mu(x);
    for (i = 0; ; i++) {
	
	Memcpy(uold, u, q);
	pwrss_old = update_L(x);
				/* tmp := PC %*% wtdResid */
	M_cholmod_sdmult(C, 0 /* notrans */, one, zero, cres, ctmp, &c);
	Memcpy(tmp1, tmp, q);
	apply_perm(tmp, tmp1, (int*)L->Perm, q);
				/* tmp := tmp - u */
	for (j = 0; j < q; j++) tmp[j] -= u[j];
				/* solve L %*% sol = tmp */
	if (!(sol = M_cholmod_solve(CHOLMOD_L, L, ctmp, &c)))
	    error(_("cholmod_solve (CHOLMOD_L) failed"));
	Memcpy(tmp, (double*)(sol->x), q);
	M_cholmod_free_dense(&sol, &c);
				/* evaluate convergence criterion */
	crit = cfac * sqr_length(tmp, q) / pwrss_old;
	if (crit < CM_TOL) break; /* don't do needless evaluations */
				/* solve t(L) %*% sol = tmp */
	if (!(sol = M_cholmod_solve(CHOLMOD_Lt, L, ctmp, &c)))
	    error(_("cholmod_solve (CHOLMOD_Lt) failed"));
	Memcpy(tmp, (double*)(sol->x), q);
	M_cholmod_free_dense(&sol, &c);

	for (step = 1; step > CM_SMIN; step /= 2) { /* step halving */
	    for (j = 0; j < q; j++) u[j] = uold[j] + step * tmp[j];
	    pwrss = update_mu(x);
	    if (verb < 0)
		Rprintf("%2d,%8.6f,%12.4g: %15.6g %15.6g %15.6g %15.6g\n",
			i, step, crit, pwrss, pwrss_old, u[1], u[2]);
	    if (pwrss < pwrss_old) {
		pwrss_old = pwrss;
		break;
	    }
	}
	if (step <= CM_SMIN || i > CM_MAXITER) return 0;
    }
    return i;
}

/**
 * Externally callable update_u.
 *
 * Iterate to determine the conditional modes of the random effects.
 *
 * @param x pointer to an mer object
 * @param verbP scalar integer indicator of verbose output
 *             (negative values produce a lot of output)
 *
 * @return number of iterations to convergence (0 for non-convergence) 
 */
SEXP mer_update_u(SEXP x, SEXP verbP)
{
    return ScalarInteger(update_u(x, asInteger(verbP)));
}

/**
 * Evaluate y * log(y/mu) with the correct limiting value at y = 0.
 *
 * @param y 
 * @param mu
 *
 * @return y * log(y/mu) for y > 0, 0 for y == 0.
 */
static R_INLINE double y_log_y(double y, double mu)
{
    return (y) ? (y * log(y/mu)) : 0;
}

/* FIXME: change this to use the variance type slot.  If
 * dims[useSc_POS] is nonzero then use the formula for the profiled
 * deviance.  (Does everything work out nicely if sigma is identically
 * 1 when useSc is false?)
 */

/**
 * Evaluate the discrepancy and log of the penalized discrepancy.
 * update_mu must be called first.
 *
 * @param x pointer to an mer object
 *
 */
SEXP mer_update_dev(SEXP x)
{
    double *d = REAL(GET_SLOT(x, lme4_devianceSym));
    int *dims = INTEGER(GET_SLOT(x, lme4_dimsSym));
    int i, n = dims[n_POS];

    d[disc_POS] = d[wrss_POS];
    if (LENGTH(GET_SLOT(x, lme4_muEtaSym))) { /* generalized */
	SEXP pwtp = GET_SLOT(x, lme4_priorWtSym);
	int pw = LENGTH(pwtp);
	double *d = REAL(GET_SLOT(x, lme4_devianceSym)),
	    *mu = REAL(GET_SLOT(x, lme4_muSym)),
	    *wts = REAL(pwtp),
	    *y = REAL(GET_SLOT(x, lme4_ySym));

	d[disc_POS] = 0;
	switch(dims[fTyp_POS]) {
	case 1:		      /* binomial with logit or probit link */
	case 2:
	    for (i = 0; i < n; i++) {
		double mui = mu[i], yi = y[i];
		d[disc_POS] += 2 * (pw ? wts[i] : 1) *
		    (y_log_y(yi, mui) + y_log_y(1 - yi, 1 - mui));
	    }
	    break;
	case 3:			/* Poisson with log link */
	    for (i = 0; i < n; i++) {
		double mui = mu[i], yi = y[i];
		d[disc_POS] += 2 * (pw ? wts[i] : 1) *
		    (y_log_y(yi, mui) - (yi - mui));
	    }
	    break;
	default:
	    error(_("General form of glmer_dev_resids not yet written"));
	}
    }

    d[ML_POS] = d[disc_POS] + d[ldL2_POS] + d[usqr_POS];
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

/* FIXME: mer_create_A is probably better done in R code using Ztl,
 * not Zt */

/**
 * Create the A matrix pattern from Zt, ST and Gp.
 *
 * @param Zt pointer to the sparse transposed model matrix
 * @param ST pointer to the list of condensed ST factors of the
 * variance-covariance terms
 * @param GpP pointer to the Gp array
 *
 * @return A
 */
SEXP mer_create_A(SEXP Zt, SEXP ST, SEXP GpP)
{
    SEXP ans;
    int *Gp = INTEGER(GpP), *nnz, *nz, *vi, *vp, *zdims, *zi, *zp,
	Vnnz, ZtOK, j, nf = LENGTH(ST);
    int *nc = Alloca(nf, int), *nlev = Alloca(nf, int);
    R_CheckStack();

				/* Trivial case, all nc == 1 */
    if (ST_nc_nlev(ST, Gp, (double**)NULL, nc, nlev) <= 1)
	return duplicate(Zt);
				/* Check the nonzero pattern in Zt */
    zdims = INTEGER(GET_SLOT(Zt, lme4_DimSym));
    nnz = Alloca(zdims[1], int);
    nz = Alloca(zdims[0], int);
    R_CheckStack();

    zi = INTEGER(GET_SLOT(Zt, lme4_iSym));
    zp = INTEGER(GET_SLOT(Zt, lme4_pSym));
    for (j = 0, ZtOK = 1; j < zdims[1]; j++) {
	nnz[j] = C_nz_col(nz, j, nf, Gp, nc, nlev, zi, zp);
	if (nnz[j] != (zp[j + 1] - zp[j])) ZtOK = 0;
    }
    if (ZtOK) return duplicate(Zt);
				/* Must create a new dgCMatrix object */ 
    ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dgCMatrix")));
    SET_SLOT(ans, lme4_DimSym, duplicate(GET_SLOT(Zt, lme4_DimSym)));
    SET_SLOT(ans, lme4_DimNamesSym, allocVector(VECSXP, 2));
				/* create and populate the p slot */
    SET_SLOT(ans, lme4_pSym, allocVector(INTSXP, zdims[1] + 1));
    vp = INTEGER(GET_SLOT(ans, lme4_pSym));
    vp[0] = 0;
    for (j = 0; j < zdims[1]; j++) vp[j + 1] = vp[j] + nnz[j];
    Vnnz = vp[zdims[1]];
				/* create and populate the i slot */
    SET_SLOT(ans, lme4_iSym, allocVector(INTSXP, Vnnz));
    vi = INTEGER(GET_SLOT(ans, lme4_iSym));
    for (j = 0; j < zdims[1]; j++) {
	int i, pos = vp[j];
	C_nz_col(nz, j, nf, Gp, nc, nlev, zi, zp);
	for (i = 0; i < zdims[0]; i++) if (nz[i]) vi[pos++] = i;
    }
				/* create and zero the x slot */
    SET_SLOT(ans, lme4_xSym, allocVector(REALSXP, Vnnz));
    AZERO(REAL(GET_SLOT(ans, lme4_xSym)), Vnnz);

    UNPROTECT(1); 
    return ans;
}

#if 0
/**
 * Created a weighted, possibly folded, copy of src in dest.
 *
 * @param dest destination pointer
 * @param src source pointer
 * @param wt weights
 * @param n length of dest vector
 * @param s fold factor
 *
 * @return squared length of weighted copy
 */
static R_INLINE double
weight_x(double* dest, const double* src, const double *wt,
	 int n, int s)
{
    int i, k;
    double ans = 0;
    
    for (i = 0; i < n; i++) {
	double accum = 0;
	for (k = 0; k < s; k++) {
	    int ii = i + k * n;
	    accum += wt[ii] * src[ii];
	}
	ans += accum * accum;
	dest[i] = accum;
    }
    return ans;
}
#endif

/**
 * Create PAX in dest.
 *
 * @param dest values to be calculated
 * @param perm NULL or a 0-based permutation vector defining P
 * @param A sparse matrix
 * @param X dense matrix
 * @param nc number of columns in X
 *
 */
static void
P_sdmult(double *dest, const int *perm, const CHM_SP A,
	 const double *X, int nc)
{
    int *ai = (int*)(A->i), *ap = (int*)(A->p), m = A->nrow, n = A->ncol;
    double *ax = (double*)(A->x), *tmp = Alloca(m, double);
    R_CheckStack();

    for (int k = 0; k < nc; k++) {
	AZERO(tmp, m);
	for (int j = 0; j < n; j++) {
	    for (int p = ap[j]; p < ap[j + 1]; p++)
		tmp[ai[p]] += X[j + k * n] * ax[p];
	}
	apply_perm(dest + k * m, tmp, perm, m);
    }
}
	
/**
 * Update the RCX and RX slots in an mer object. update_L should be
 * called before update_RX
 *
 * @param x pointer to an mer object
 *
 * @return profiled deviance or REML deviance
 */
static double update_RX(SEXP x)
{
    int *dims = INTEGER(GET_SLOT(x, lme4_dimsSym)), *perm;
    int i, i1 = 1, j, k, n = dims[n_POS];
    int p = dims[p_POS], q = dims[q_POS], s = dims[s_POS];
    double *X = REAL(GET_SLOT(x, lme4_XSym)),
	*cx = SLOT_REAL_NULL(x, lme4_CxSym),
	*d = REAL(GET_SLOT(x, lme4_devianceSym)), *fixef,
	*rzxy, *rxy = REAL(GET_SLOT(x, lme4_RXySym)),
	*sXwt = SLOT_REAL_NULL(x, lme4_sqrtXWtSym), *u,
	*y = REAL(GET_SLOT(x, lme4_ySym)),
    	dn = (double) n, dnmp = (double) (n - p),
	mone[] = {-1,0}, one[] = {1,0}, zero[] = {0,0};
    CHM_SP A = AS_CHM_SP(GET_SLOT(x, lme4_ASym));
    CHM_DN cRCXy = AS_CHM_DN(GET_SLOT(x, lme4_RCXySym)), ans, cu;
    CHM_FR L = L_SLOT(x);
    R_CheckStack();

    rzxy = (double*)cRCXy->x;
    perm = (int*)L->Perm;
    if (sXwt) {			/* Create W^{1/2}GHX in WX */
	double *WX = Calloc(n * p, double);

	AZERO(WX, n * p);
	for (j = 0; j < p; j++)
	    for (k = 0; k < s; k++)
		for (i = 0; i < n; i++)
		    WX[i + j * n] +=
			sXwt[i + k * n] * X[i + n * (k + j * s)];
	if (cx) A->x = (void*)cx;
	else {			/* Replace A by C */
	    A = AS_CHM_SP(GET_SLOT(x, lme4_CmSym));
	    R_CheckStack();
	}
	X = WX;
    }
				/* initialize RZX to PAX */
    P_sdmult(rzxy, perm, A, X, p);
				/* initialize RX to X'X */
    F77_CALL(dsyrk)("U", "T", &p, &n, one, X, &n, zero, rxy, &p);
				/* solve L %*% ans = PAW^{1/2}[X:y] */
    ans = M_cholmod_solve(CHOLMOD_L, L, cRCXy, &c);
    Memcpy(rzxy, (double*)(ans->x), q * p);
    M_cholmod_free_dense(&ans, &c);
    				/* downdate X'X and factor  */
    F77_CALL(dsyrk)("U", "T", &p, &q, mone, rzxy, &q, one, rxy, &p);
    F77_CALL(dpotrf)("U", &p, rxy, &p, &j);
    if (j)
	error(_("Downdated X'X is not positive definite, %d."), j);
				/* accumulate log(det(RX)^2)  */
    for (j = 0, d[ldRX2_POS] = 0; j < p; j++)
	d[ldRX2_POS] += 2 * log(rxy[j * (p + 1)]);

    if (SLOT_REAL_NULL(x, lme4_vSym) || SLOT_REAL_NULL(x, lme4_muEtaSym)) {
	Free(X);
	return d[ML_POS];
    }
				/* linear mixed models only */
    fixef = REAL(GET_SLOT(x, lme4_fixefSym));
    u = REAL(GET_SLOT(x, lme4_uSym));
    cu = N_AS_CHM_DN(u, q, 1);
    R_CheckStack();
				/* solve L del1 = PAy */
    if (sXwt) {
	double *tmp = Calloc(n, double);
	for (i = 0; i < n; i++) tmp[i] = sXwt[i] * y[i];
	y = tmp;
    }
    P_sdmult(u, perm, A, y, 1);
    ans = M_cholmod_solve(CHOLMOD_L, L, cu, &c);
    Memcpy(u, (double*)ans->x, q);
    M_cholmod_free_dense(&ans, &c);
				/* solve RX' del2 = X'y - RZX'del1 */
    F77_CALL(dgemv)("T", &n, &p, one, X, &n, y, &i1, zero, fixef, &i1);
    F77_CALL(dgemv)("T", &q, &p, mone, rzxy, &q,
		    u, &i1, one, fixef, &i1);
    F77_CALL(dtrsv)("U", "T", "N", &p, rxy, &p, fixef, &i1);
    d[pwrss_POS] = sqr_length(y, n)
	- (sqr_length(fixef, p) + sqr_length(u, q));
    if (d[pwrss_POS] < 0)
	error(_("Calculated PWRSS for a LMM is negative"));
    d[ML_POS] = d[ldL2_POS] +
	dn * (1 + log(d[pwrss_POS]) + log(2 * PI / dn));
    d[REML_POS] = d[ldL2_POS] + d[ldRX2_POS] +
	dnmp * (1. + log(d[pwrss_POS]) + log(2. * PI / dnmp));

    if (sXwt) Free(y);
    return d[dims[isREML_POS] ? REML_POS : ML_POS];
}

/**
 * Externally callable update_RX
 *
 * @param x pointer to an mer object
 *
 * @return profiled deviance or REML deviance
 */
SEXP
mer_update_RX(SEXP x)
{
    return ScalarReal(update_RX(x));
}

/**
 * Optimize the profiled deviance of an lmer object or the Laplace
 * approximation to the deviance of a nlmer or glmer object.
 *
 * @param x pointer to an mer object
 * @param verbp pointer to indicator of verbose output
 *
 * @return R_NilValue
 */
SEXP 
mer_optimize(SEXP x, SEXP verbp)
{
    SEXP ST = GET_SLOT(x, lme4_STSym);
    int *dims = INTEGER(GET_SLOT(x, lme4_dimsSym)),
	i, j, pos, verb = asInteger(verbp);
    int lmm = !(LENGTH(GET_SLOT(x, lme4_muEtaSym)) ||
		LENGTH(GET_SLOT(x, lme4_vSym))),
	nf = dims[nf_POS];
/* FIXME: need to add 1 to nv for GLMM with a scale par (maybe).
 * Not really but we do need to profile out the scale parameter when
 * useSc is TRUE. */
    int nv = dims[np_POS] + (lmm ? 0 : dims[p_POS]);
    int liv = S_iv_length(OPT, nv), lv = S_v_length(OPT, nv);
    double *g = (double*)NULL, *h = (double*)NULL, fx = R_PosInf;
    double *fixef = REAL(GET_SLOT(x, lme4_fixefSym));
    int *iv = Alloca(liv, int);
    double *b = Alloca(2 * nv, double), *d = Alloca(nv, double),
	*v = Alloca(lv, double), *xv = Alloca(nv, double);
    R_CheckStack();

    ST_getPars(x, xv);
    if (!lmm) {
	double eta = 1.e-5; /* estimated rel. error on computed lpdisc */

	Memcpy(xv + dims[np_POS], fixef, dims[p_POS]);
	v[31] = eta;		/* RFCTOL */
	v[36] = eta;		/* SCTOL */
	v[41] = eta;		/* ETA0 */
    }
				/* initialize the state vectors v and iv */
    S_Rf_divset(OPT, iv, liv, lv, v);
    if (verb) iv[OUTLEV] = 1;
				/* set the bounds to plus/minus Infty  */
    for (i = 0; i < nv; i++) {
	b[2*i] = R_NegInf; b[2*i+1] = R_PosInf; d[i] = 1;
    }
				/* reset lower bounds on elements of theta_S */
    for (i = 0, pos = 0; i < nf; i++) {
	int nc = *INTEGER(getAttrib(VECTOR_ELT(ST, i), R_DimSymbol));
	for (j = 0; j < nc; j++) b[pos + 2*j] = 0;
	pos += nc * (nc + 1);
    }
    do {
	ST_setPars(x, xv);		/* update ST and C */
/* FIXME: Change this so that update_dev is always called and that
 * function does the selection of methods based on lmm. The Memcpy
 * call should always be used but the number of values to copy can be
 * zero. */
	if (lmm) {
	    update_L(x);
	    fx = update_RX(x);
	} else {
	    
	    Memcpy(fixef, xv + dims[np_POS], dims[p_POS]);
	    update_u(x, verb);
	    mer_update_dev(x);
	    fx = REAL(GET_SLOT(x, lme4_devianceSym))[ML_POS];
	}
	S_nlminb_iterate(b, d, fx, g, h, iv, liv, lv, nv, v, xv);
    } while (iv[0] == 1 || iv[0] == 2);
    update_RX(x);
    dims[cvg_POS] = iv[0];
    return R_NilValue;
}

/* Functions for sampling from a Wishart distribution */
/* FIXME: Move these to the R sources */

/**
 * Simulate the Cholesky factor of a standardized Wishart variate with
 * dimension p and df degrees of freedom.
 *
 * @param df degrees of freedom
 * @param p dimension of the Wishart distribution
 * @param upper if 0 the result is lower triangular, otherwise upper
                triangular
 * @param ans array of size p * p to hold the result
 *
 * @return ans
 */
static double
*std_rWishart_factor(double df, int p, int upper, double ans[])
{
    int i, j, pp1 = p + 1;

    if (df < (double) p || p <= 0)
	error("inconsistent degrees of freedom and dimension");
    for (j = 0; j < p; j++) {	/* jth column */
	ans[j * pp1] = sqrt(rchisq(df - (double) j));
	for (i = 0; i < j; i++) {
	    int uind = i + j * p, /* upper triangle index */
		lind = j + i * p; /* lower triangle index */
	    ans[(upper ? uind : lind)] = norm_rand();
	    ans[(upper ? lind : uind)] = 0;
	}
    }
    return ans;
}

/**
 * Simulate a sample of random matrices from a Wishart distribution
 *
 * @param ns Number of samples to generate
 * @param dfp Degrees of freedom
 * @param scal Positive-definite scale matrix
 *
 * @return
 */
SEXP
lme4_rWishart(SEXP ns, SEXP dfp, SEXP scal)
{
    SEXP ans;
    int *dims = INTEGER(getAttrib(scal, R_DimSymbol)), i, j, k,
	n = asInteger(ns), psqr;
    double *scCp, *ansp, *tmp, df = asReal(dfp), one = 1, zero = 0;

    if (!isMatrix(scal) || !isReal(scal) || dims[0] != dims[1])
	error("scal must be a square, real matrix");
    if (n <= 0) n = 1;
    psqr = dims[0] * dims[0];
    tmp = Alloca(psqr, double);
    scCp = Alloca(psqr, double);
    R_CheckStack();

    Memcpy(scCp, REAL(scal), psqr);
    AZERO(tmp, psqr);
    F77_CALL(dpotrf)("U", &(dims[0]), scCp, &(dims[0]), &j);
    if (j)
	error("scal matrix is not positive-definite");
    PROTECT(ans = alloc3DArray(REALSXP, dims[0], dims[0], n));
    ansp = REAL(ans);
    GetRNGstate();
    for (j = 0; j < n; j++) {
	double *ansj = ansp + j * psqr;
	std_rWishart_factor(df, dims[0], 1, tmp);
	F77_CALL(dtrmm)("R", "U", "N", "N", dims, dims,
			&one, scCp, dims, tmp, dims);
	F77_CALL(dsyrk)("U", "T", &(dims[1]), &(dims[1]),
			&one, tmp, &(dims[1]),
			&zero, ansj, &(dims[1]));

	for (i = 1; i < dims[0]; i++)
	    for (k = 0; k < i; k++)
		ansj[i + k * dims[0]] = ansj[k + i * dims[0]];
    }

    PutRNGstate();
    UNPROTECT(1);
    return ans;
}

/* MCMC sampling functions */

/**
 * Update the fixed effects in an MCMC sample from an mer object.
 *
 * @param x an lmer model with mer_update_b
 * @param sigma current standard deviation of the per-observation
 *        noise terms.
 */
static void MCMC_beta(SEXP x, double sigma)
{
    SEXP fixefP = GET_SLOT(x, lme4_fixefSym),
	uvecP = GET_SLOT(x, lme4_uSym);
    int i1 = 1, j, p = LENGTH(fixefP), q = LENGTH(uvecP);
    int pp1 = p + 1;
    double *fix = REAL(fixefP), *uv = REAL(uvecP), mone = -1, one = 1;
    CHM_FR L = L_SLOT(x);
    double *del1 = Alloca(q, double), *del2 = Alloca(p, double);
    CHM_DN sol, rhs = N_AS_CHM_DN(del1, q, 1);
    R_CheckStack();

    update_ranef(x);
				/* Update the fixed effects */
    for (j = 0; j < p; j++) del2[j] = sigma * norm_rand();
    F77_CALL(dtrsv)("U", "N", "N", &p, REAL(GET_SLOT(x, lme4_RXySym)),
		    &pp1, del2, &i1);
    for (j = 0; j < p; j++) fix[j] += del2[j];
				/* Update the orthogonal random effects */
    for (j = 0; j < q; j++) del1[j] = sigma * norm_rand();
    F77_CALL(dgemv)("N", &q, &p, &mone, REAL(GET_SLOT(x, lme4_RCXySym)),
		    &q, del2, &i1, &one, del1, &i1);
    sol = M_cholmod_solve(CHOLMOD_Lt, L, rhs, &c);
    for (j = 0; j < q; j++) uv[j] += ((double*)(sol->x))[j];
    M_cholmod_free_dense(&sol, &c);
				/* Update the random effects */
    update_ranef(x);
}

/**
 * Update the orthogonal random effects in an MCMC sample
 * from an mer object.
 *
 * @param x an lmer model with mer_update_b
 * @param sigma current standard deviation of the per-observation
 *        noise terms.
 */
static void MCMC_u(SEXP x, double sigma)
{
}

static void
safe_pd_matrix(double x[], const char uplo[], int n, double thresh)
{
    int info, lwork = 3 * n, nm1 = n - 1;
    double *work = Alloca(3 * n, double),
	*w = Alloca(n, double),
	*xcp = Memcpy(Alloca(n * n, double), x, n * n);

    F77_CALL(dsyev)("N", uplo, &n, xcp, &n, w, work, &lwork, &info);
    if (info) error(_("dsyev returned %d"), info);
    if (w[nm1] <= 0) error(_("no positive eigenvalues!"));
    if ((w[0]/w[nm1]) < thresh) {
	int i, np1 = n + 1;
	double incr = w[nm1] * thresh;
	for (i = 0; i < n; i++) x[i * np1] += incr;
    }
}

static void MCMC_ST(SEXP x, double sigma, double *vals, int trans)
{
    SEXP ST = GET_SLOT(x, lme4_STSym);
    int i, info, j, k, ncmax, nf = LENGTH(ST),
	*Gp = INTEGER(GET_SLOT(x, lme4_GpSym));
    double *rr = REAL(GET_SLOT(x, lme4_ranefSym)),
	*scal, *var, *wfac, one = 1, zero = 0;
    double **st = Alloca(nf, double*);
    int *nc = Alloca(nf, int), *nlev = Alloca(nf, int);
    R_CheckStack();

    ncmax = ST_nc_nlev(ST, Gp, st, nc, nlev);
    scal = Alloca(ncmax * ncmax, double);
    wfac = Alloca(ncmax * ncmax, double);
    var = Alloca(ncmax * ncmax, double);
    R_CheckStack();
    for (i = 0; i < nf; i++) {
	int nci = nc[i], nl = nlev[i];
	int ncisq = nci * nci, ncip1 = nci + 1;
	
	AZERO(scal, ncisq);
	F77_CALL(dsyrk)("U", "T", &nci, &nl, &one, rr + Gp[i], &nl,
			&zero, scal, &nci); 
	if (nci == 1) *scal = sqrt(*scal);
	else {
	    safe_pd_matrix(scal, "L", nci, 1e-7);
	    F77_CALL(dpotrf)("U", &nci, scal, &nci, &info);
	    if (info)
		error(_("scal matrix is not pd after safe_pd_matrix!"));
	}
				/* generate random factor from std Wishart */
	AZERO(wfac, ncisq);
	std_rWishart_factor((double) (nl - nci + 1), nci, 1, wfac);
				/* form and store variance-covariance matrix */
	F77_CALL(dtrsm)("L", "U", "T", "N", &nci, &nci,
			&one, wfac, &nci, scal, &nci);
	F77_CALL(dsyrk)("U", "T", &nci, &nci, &one, scal, &nci,
			&zero, var, &nci);
	for (j = 0; j < nci; j++) {
	    *vals++ = (trans ? log(var[j * ncip1]) : var[j * ncip1]);
	}
	for (j = 1; j < nci; j++) {
	    for (k = 0; k < j; k++) {
		*vals++ = (trans ? atanh(var[k + j * nci]/
				     sqrt(var[j * ncip1] * var[k * ncip1]))
		       : var[k + j * nci]);
	    }
	}
	if (nci == 1) {
	    *st[i] = sqrt(*var)/sigma;
	} else {
	    double *sti = st[i];

	    for (j = 0; j < nci; j++) /* copy var to st[i] transposing */
		for (k = j; k < nci; k++)
		    sti[k + j * nci] = var[j + k * nci];
				/* factor the variance-covariance matrix */
	    F77_CALL(dpotrf)("L", &nci, sti, &nci, &info);
	    if (info)
		error(_("new var-cov matrix is not positive definite!"));
	    for (j = 0; j < nci; j++) { /* convert to LDL' then TSST' */
		for (k = j + 1; k < nci; k++)
		    sti[k + j * nci] /= sti[j * ncip1];
		sti[j * ncip1] /= sigma;
	    }
	}
    }
}

/**
 * Generate a new value of sigma, the common scale parameter, for an
 * mer object.
 *
 * @param dims dims slot of the mer object
 * @param d  deviance slot of the mer object
 *
 * @return an updated sigma
 */
static R_INLINE double MCMC_sigma(const int* dims, const double* d)
{
    return (dims[useSc_POS] ? sqrt(d[pwrss_POS]/rchisq(RDF(dims))) : 1);
}

/**
 * Generate a Markov-chain Monte Carlo sample from an mer object
 *
 * @param x pointer to an lmer object
 * @param savebp logical scalar - save random effects?
 * @param nsampp integer scalar - number of samples to generate
 * @param transp logical scalar - transform variance components?
 * @param verbosep logical scalar - verbose output?
 * @param deviancep logical scalar - store the deviance?
 *
 * @return a matrix
 */
SEXP mer_MCMCsamp(SEXP x, SEXP savebp, SEXP nsampp, SEXP transp,
		  SEXP verbosep, SEXP deviancep)
{
    SEXP Pars = PROTECT(mer_ST_getPars(x)), ans;
    int *dims = INTEGER(GET_SLOT(x, lme4_dimsSym));
    int dPOS = dims[REML_POS] ? REML_POS : ML_POS, i, j,
	nsamp = asInteger(nsampp),
	p = dims[p_POS], q = dims[q_POS],
	saveb = asLogical(savebp), trans = asLogical(transp),
	verbose = asLogical(verbosep), dev = asLogical(deviancep);
    double *ansp, *deviance = REAL(GET_SLOT(x, lme4_devianceSym)),
	*d = REAL(GET_SLOT(x, lme4_devianceSym)),
	*fixef = REAL(GET_SLOT(x, lme4_fixefSym)),
	*ranef = (double*)NULL;
    int nrbase = p + 1 + dims[np_POS]; /* rows always included */
    int nrtot = nrbase + dev + (saveb ? q : 0);

    if (nsamp <= 0) nsamp = 1;
    ans = PROTECT(allocMatrix(REALSXP, nrtot, nsamp));
    ansp = REAL(ans);
    for (i = 0; i < nrtot * nsamp; i++) ansp[i] = NA_REAL;
    GetRNGstate();
    if (verbose) Rprintf("%12s %14s\n", "sigma", "deviance");
    if (saveb) ranef = REAL(GET_SLOT(x, lme4_ranefSym));

    for (i = 0; i < nsamp; i++) {
	double *col = ansp + i * nrtot, sigma;
				/* simulate and store new value of sigma */
	sigma = MCMC_sigma(dims, d);
	col[p] = (trans ? 2 * log(sigma) : sigma * sigma);
	update_ranef(x);	/* Conditional estimates of beta, u and b */
	MCMC_beta(x, sigma);	/* Update beta */
	for (j = 0; j < p; j++) col[j] = fixef[j]; /* Store beta */
	MCMC_u(x, sigma);
	if (saveb) {				   /* store b */
	    update_ranef(x);
	    for (j = 0; j < q; j++) col[nrbase + dev + j] = ranef[j];
	}
	MCMC_ST(x, sigma, col + p + 1, trans);
/* 	lmer_update_dev(x);	*/ /* Refactor and evaluate deviance */
	/* FIXME: This deviance should be for sigma, beta, b, ST */
	if (dev) col[nrbase] = deviance[dPOS]; 
	if (verbose) Rprintf("%12.6g %14.8g\n", sigma, deviance[dPOS]);
    }
    PutRNGstate();
				/* Restore pars, refactor, etc. */
    mer_ST_setPars(x, Pars);
/*     lmer_update_dev(x); */
/*     lmer_update_effects(x); */
    UNPROTECT(2);
    return ans;
}

#ifndef BUF_SIZE
/** size of buffer for an error message */
#define BUF_SIZE 127
#endif	

/**
 * Check that slot sym of object x is a numeric matrix of dimension nr
 * by nc.
 *
 * @param buf character buffer of length nb + 1
 * @param nb number of writable positions in the buffer
 * @param x pointer to an mer object
 * @param sym name (symbol, actually) of the slot to check
 * @param nr expected number of rows
 * @param nc expected number of columns
 *
 * @return 0 for success, number of characters written to buf on failure
 */
static int chkDims(char *buf, int nb, SEXP x, SEXP sym, int nr, int nc)
{
    SEXP MP = GET_SLOT(x, sym);
    int *dm = isMatrix(MP) ? INTEGER(getAttrib(MP, R_DimSymbol)) : (int*) NULL;
    if (!dm || !isReal(MP) || dm[0] != nr || dm[1] != nc)
	return snprintf(buf, BUF_SIZE,
			_("Slot %s must be a numeric matrix of size %d by %d"),
			CHAR(PRINTNAME(sym)), nr, nc);
    return 0;
}

/** Check that the length of the sym slot in x is len or, possibly, zero.
 *
 * @param buf character buffer of length nb + 1
 * @param nb number of writable positions in the buffer
 * @param x pointer to an mer object
 * @param sym name (symbol, actually) of the slot to check
 * @param len expected length
 * @param zerok is a length of zero allowed?
 *
 * @return 0 for success, number of characters written to buf on failure

 */
static int chkLen(char *buf, int nb, SEXP x, SEXP sym, int len, int zerok)
{
    int ll;

    if (!(ll = LENGTH(GET_SLOT(x, sym))) == len && zerok && !ll)
	return snprintf(buf, BUF_SIZE, _("Slot %s must have length %d."),
			CHAR(PRINTNAME(sym)), len);
    return 0;
}

/**
 * Check validity of an object that inherits from the mer class.
 *
 * @param x Pointer to an lmer or glmer or nlmer object
 *
 * @return TRUE if the object is a valid mer object, otherwise a string
 *         that describes the violation.
 */
SEXP mer_validate(SEXP x)
{
    SEXP GpP = GET_SLOT(x, lme4_GpSym),
	ST = GET_SLOT(x, lme4_STSym),
	devianceP = GET_SLOT(x, lme4_devianceSym),
	dimsP = GET_SLOT(x, lme4_dimsSym),
	flistP = GET_SLOT(x, lme4_flistSym);
    int *Gp = INTEGER(GpP), *dd = INTEGER(dimsP);
    int i, n = dd[n_POS], nf = dd[nf_POS], nq, p = dd[p_POS],
	q = dd[q_POS], s = dd[s_POS];
    int nv = n * s;
    CHM_SP Zt = AS_CHM_SP(GET_SLOT(x, lme4_ZtSym)),
	A =  AS_CHM_SP(GET_SLOT(x, lme4_ASym));
    CHM_FR L = L_SLOT(x);
    char *buf = Alloca(BUF_SIZE + 1, char);
    R_CheckStack();
				/* check lengths */
    if (nf < 1 || LENGTH(flistP) != nf || LENGTH(ST) != nf)
	return mkString(_("Slots ST, and flist must have length dims['nf']"));
    if (LENGTH(GpP) != nf + 1)
	return mkString(_("Slot Gp must have length dims['nf'] + 1"));
    if (Gp[0] != 0 || Gp[nf] != q)
	return mkString(_("Gp[1] != 0 or Gp[dims['nf'] + 1] != dims['q']"));
    if (LENGTH(devianceP) != (wrss_POS + 1) ||
	LENGTH(getAttrib(devianceP, R_NamesSymbol)) != (wrss_POS + 1))
	return mkString(_("deviance slot not named or incorrect length"));
    if (LENGTH(dimsP) != (cvg_POS + 1) ||
	LENGTH(getAttrib(dimsP, R_NamesSymbol)) != (cvg_POS + 1))
	return mkString(_("dims slot not named or incorrect length"));
    if (L->n != q || !L->is_ll || !L->is_monotonic)
	return mkString(_("Slot L must be a monotonic LL' factorization of size dims['q']"));
    if (Zt->nrow != q || Zt->ncol != nv)
	return mkString(_("Slot Zt must by dims['q']  by dims['n']*dims['s']"));
    if (A->nrow != q || A->ncol != nv)
	return mkString(_("Slot A must be dims['q']  by dims['n']*dims['s']"));
    if (chkLen(buf, BUF_SIZE, x, lme4_etaSym, n, 0)) return(mkString(buf));
    if (chkLen(buf, BUF_SIZE, x, lme4_fixefSym, p, 0)) return(mkString(buf));
    if (chkLen(buf, BUF_SIZE, x, lme4_muEtaSym, n, 1)) return(mkString(buf));
    if (chkLen(buf, BUF_SIZE, x, lme4_muSym, n, 0)) return(mkString(buf));
    if (chkLen(buf, BUF_SIZE, x, lme4_offsetSym, n, 1)) return(mkString(buf));
    if (chkLen(buf, BUF_SIZE, x, lme4_priorWtSym, n, 1)) return(mkString(buf));
    if (chkLen(buf, BUF_SIZE, x, lme4_ranefSym, q, 0)) return(mkString(buf));
    if (chkLen(buf, BUF_SIZE, x, lme4_residSym, n, 0)) return(mkString(buf));
    if (chkLen(buf, BUF_SIZE, x, lme4_sqrtrWtSym, n, 1)) return(mkString(buf));
    if (chkLen(buf, BUF_SIZE, x, lme4_uSym, q, 0)) return(mkString(buf));
    if (chkLen(buf, BUF_SIZE, x, lme4_vSym, nv, 1)) return(mkString(buf));
    if (chkLen(buf, BUF_SIZE, x, lme4_varSym, n, 1)) return(mkString(buf));
    if (chkLen(buf, BUF_SIZE, x, lme4_ySym, n, 0)) return(mkString(buf));
    if (chkDims(buf, BUF_SIZE, x, lme4_XSym, nv, p)) return(mkString(buf));
    if (chkDims(buf, BUF_SIZE, x, lme4_RCXySym, q, p)) return(mkString(buf));
    if (chkDims(buf, BUF_SIZE, x, lme4_RXySym, p, p)) return(mkString(buf));

    nq = 0;
    for (i = 0; i < nf; i++) {
	SEXP STi = VECTOR_ELT(ST, i), fli = VECTOR_ELT(flistP, i);
	int *dm = INTEGER(getAttrib(STi, R_DimSymbol));
	if (!isMatrix(STi) || !isReal(STi) || dm[0] != dm[1])
	    return
		mkString(_("Slot ST must be a list of square numeric matrices"));
	if (Gp[i] > Gp[i + 1])
	    return mkString(_("Gp must be non-decreasing"));
	if (!isFactor(fli))
	    return mkString(_("flist must be a list of factors"));
	nq += dm[0] * LENGTH(getAttrib(fli, R_LevelsSymbol));
    }
    if (q != nq)
	return mkString(_("q is not sum of columns by levels"));
    return ScalarLogical(1);
}

/* Functions for the pedigree class */

/**
 * Create the left Cholesky factor of the numerator relationship
 * matrix from a pedigree.
 *
 * @param x a pedigree object
 * @param ans T stored as a non-unit, lower dtCMatrix
 *
 * @return ans with elements modified to incorporate D
 */
SEXP pedigree_chol(SEXP x, SEXP ans)
{
    SEXP Sire = GET_SLOT(x, install("sire"));
    int *ai = INTEGER(GET_SLOT(ans, lme4_iSym)),
	*ap = INTEGER(GET_SLOT(ans, lme4_pSym)),
	*dam = INTEGER(GET_SLOT(x, install("dam"))),
	*sire = INTEGER(Sire), 
	i, j, n = LENGTH(Sire);
    double *ax = REAL(GET_SLOT(ans, lme4_xSym)), *F, Di, tmp;

    setAttrib(ans, install("F"), allocVector(REALSXP, n));
    F = REAL(getAttrib(ans, install("F")));
    for (i = 0; i < n; i++) {
	int k, p = sire[i] - 1, q = dam[i] - 1;
	if (sire[i] == NA_INTEGER) {
	    F[i] = 0;
	    Di = (dam[i] == NA_INTEGER) ? 1 : sqrt(0.75 - 0.25 * F[q]);
	} else {
	    if (dam[i] == NA_INTEGER) { /* sire only */
		F[i] = 0;
		Di = sqrt(0.75 - 0.25 * F[p]);
	    } else {		/* both parents in pedigree */
		Di = sqrt(0.5 - 0.25 * (F[p] + F[q]));
		F[i] = NA_REAL;
		if ((ap[i + 1] - ap[i]) > 1) {	  /* skip if no progeny */
		    if (p > q) {j = p; p = q; q = j;} /* ensure p <= q */
		    for (j = 0, F[i] = 0; j <= p; j++) {
			for (k = ap[j], tmp = 0;
			     k < ap[j + 1] && ai[k] <= q; k++) {
			    int ii = ai[k];
			    if (ii == p) tmp = ax[k];
			    if (ii == q) F[i] += tmp * ax[k]/2;
			}
		    }
		}
	    }
	}
	for (j = ap[i]; j < ap[i + 1]; j++) ax[j] *= Di;
    }
    return ans;
}

/* NOTE: This function requires that missing parents be coded as zero */
/**
 * Create the inbreeding coefficients according to the algorithm given
 * in "Comparison of four direct algorithms for computing inbreeding
 * coefficients" by Mehdi Sargolzaei and Hiroaki Iwaisaki, Animal
 * Science Journal (2005) 76, 401--406.  This function is a modified
 * version of the code published in an appendix to that paper.
 *
 * @param x a pedigree object
 *
 * @return a list of the inbreeding coefficients
 */

SEXP pedigree_inbreeding(SEXP x)
{
    SEXP ans, sp = GET_SLOT(x, install("sire"));
    int i, j, t, n = LENGTH(sp), S, D;
    int *SI, *MI,		/* start and minor */
	*sire = INTEGER(sp),
	*dam = INTEGER(GET_SLOT(x, install("dam")));
    double *F = Alloca(n + 1, double), /* inbreeding coefficients */
      *L = Alloca(n + 1, double), *B = Alloca(n + 1, double);
    int *Anc = Alloca(n + 1, int),	/* ancestor */
	*LAP = Alloca(n + 1, int); 	/* longest ancestoral path */
    R_CheckStack();
    
    F[0] =-1; LAP[0] =-1; /* set F and lap for unknown parents */
    for(i = 1, t = -1; i <= n; i++) { 	/* evaluate LAP and its maximum */
	S = sire[i]; D = dam[i]; /* parents of animal i */
	LAP[i] = ((LAP[S] < LAP[D]) ? LAP[D] : LAP[S]) + 1;
	if (LAP[i] > t) t = LAP[i];
    }
    SI = Alloca(t + 1, int);
    MI = Alloca(t + 1, int);
    for(i = 0; i <= t ; ++i) SI[i] = MI[i] = 0; /* initialize start and minor */
    for(i = 1; i <= n; i++) { 	/* evaluate F */
	S = sire[i]; D = dam[i]; /* parents of animal i */
	B[i] = 0.5 - 0.25 * (F[S] + F[D]); 
				/* adjust start and minor */
	for (j = 0; j < LAP[i]; j++) {++SI[j]; ++MI[j];} 
	if (S == 0 || D == 0) { /* both parents unknown */
	    F[i] = L[i] = 0; continue;
	}
	if(S == sire[i-1] && D == dam[i-1]) { /* full-sib with last animal */
	    F[i] = F[i-1]; L[i] = L[i-1]; continue;
	}
    
	F[i] = -1; L[i] = 1; 
	t = LAP[i]; /* largest lap group number in the animal's pedigree */
	Anc[MI[t]++] = i; /* initialize Anc and increment MI[t] */
	while(t > -1) { /* from the largest lap group to zero */
	    j = Anc[--MI[t]]; /* next ancestor */
	    S = sire[j]; D = dam[j]; /* parents of the ancestor */
	    if (S) {
		if (!L[S]) Anc[MI[LAP[S]]++] = S; 
				/* add sire in its lap group in Anc
				 * array if it is not added yet and
				 * increment the minor index for the group */ 
		L[S] += 0.5 * L[j]; /* contribution to sire */
	    }
	    if (D) {
		if (!L[D]) Anc[MI[LAP[D]]++] = D;
		L[D] += 0.5 * L[j]; /* contribution to dam */
	    }
	    F[i] += L[j] * L[j] * B[j];
	    L[j] = 0; /*clear L[j] for the evaluation of the next animal */
	    if (MI[t] == SI[t]) --t; /* move to the next lap group when
				      * all ancestors in group t have been
				      * evaluated */
	} 
    }
    ans = PROTECT(allocVector(REALSXP, n));
    Memcpy(REAL(ans), F + 1, n);
    UNPROTECT(1);
    return ans;
}


#if 0

/* Gauss-Hermite Quadrature x positions and weights */
static const double
    GHQ_x1[1] = {0},
    GHQ_w1[1] = {1},
    GHQ_x2[1] = {1},
    GHQ_w2[1] = {0.5},
    GHQ_x3[2] = {1.7320507779261, 0},
    GHQ_w3[2] = {0.166666666666667, 0.666666666666667},
    GHQ_x4[2] = {2.3344141783872, 0.74196377160456},
    GHQ_w4[2] = {0.0458758533899086, 0.454124131589555},
    GHQ_x5[3] = {2.85696996497785, 1.35562615677371, 0},
    GHQ_w5[3] = {0.0112574109895360, 0.222075915334214,
		 0.533333317311434},
    GHQ_x6[3] = {3.32425737665988, 1.88917584542184,
		 0.61670657963811},
    GHQ_w6[3] = {0.00255578432527774, 0.0886157433798025,
		 0.408828457274383},
    GHQ_x7[4] = {3.7504396535397, 2.36675937022918,
		 1.15440537498316, 0},
    GHQ_w7[4] = {0.000548268839501628, 0.0307571230436095,
		 0.240123171391455, 0.457142843409801},
    GHQ_x8[4] = {4.14454711519499, 2.80248581332504,
		 1.63651901442728, 0.539079802125417},
    GHQ_w8[4] = {0.000112614534992306, 0.0096352198313359,
		 0.117239904139746, 0.373012246473389},
    GHQ_x9[5] = {4.51274578616743, 3.20542894799789,
		 2.07684794313409, 1.02325564627686, 0},
    GHQ_w9[5] = {2.23458433364535e-05, 0.00278914123744297,
		 0.0499164052656755, 0.244097495561989,
		 0.406349194142045},
    GHQ_x10[5] = {4.85946274516615, 3.58182342225163,
		  2.48432579912153, 1.46598906930182,
		  0.484935699216176},
    GHQ_w10[5] = {4.31065250122166e-06, 0.000758070911538954,
		  0.0191115799266379, 0.135483698910192,
		  0.344642324578594},
    GHQ_x11[6] = {5.18800113558601, 3.93616653976536,
		  2.86512311160915, 1.87603498804787,
		  0.928868981484148, 0},
    GHQ_w11[6] = {8.12184954622583e-07, 0.000195671924393029,
		  0.0067202850336527, 0.066138744084179,
		  0.242240292596812, 0.36940835831095};

static const double
    *GHQ_x[12] = {(double *) NULL, GHQ_x1, GHQ_x2, GHQ_x3, GHQ_x4,
		  GHQ_x5, GHQ_x6, GHQ_x7, GHQ_x8, GHQ_x9, GHQ_x10,
		  GHQ_x11},
    *GHQ_w[12] = {(double *) NULL, GHQ_w1, GHQ_w2, GHQ_w3, GHQ_w4,
		  GHQ_w5, GHQ_w6, GHQ_w7, GHQ_w8, GHQ_w9, GHQ_w10,
		  GHQ_w11};


#endif
