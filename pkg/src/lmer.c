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
    lpdisc_POS,			/**<log penalized discrepancy */
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
    nest_POS,			/**<indicator of nested grouping factors */
    useSc_POS,			/**<does the family use a separate scale parameter */
    cvg_POS			/**<convergence indictor from port optimization  */
};

/** Extract the L slot, convert to a cholmod_factor and return the pointer. */
#define L_SLOT(x) AS_CHM_FR(GET_SLOT(x, lme4_LSym))

/**
 * Return the sum of squares of the first n elements of x
 *
 * @param n
 * @param x
 *
 * @return sum of squares
 */
static double R_INLINE sqr_length(const double *x, int n)
{
    double ans = 0;
    int i;

    for (i = 0; i < n; i++) ans += x[i] * x[i];
    return ans;
}

/**
 * Return the group in the (nf, Gp) combination to which ind belongs
 *
 * @param ind a row number
 * @param nf number of groups of rows
 * @param Gp group pointers
 *
 */
static int R_INLINE Gp_grp(int ind, int nf, const int *Gp)
{
    int i;
    for (i = 0; i < nf; i++) if (ind < Gp[i + 1]) return i;
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
static double R_INLINE
*SLOT_REAL_NULL(SEXP obj, SEXP nm)
{
    SEXP pt = GET_SLOT(obj, nm);
    return LENGTH(pt) ? REAL(pt) : (double*) NULL;
}

/**
 * Return the REML or ML conditional estimate of sigma, the standard
 * deviation of the per-observation noise term.
 *
 * @param REML non-zero for REML estimate, 0 for ML estimate
 * @param dims vector of dimensions
 * @param deviance vector of deviance components
 */
static double R_INLINE
get_sigma(int REML, const int* dims, const double* deviance)
{
    double sigma = 1;

    if (dims[useSc_POS])
	sigma =
	    sqrt((deviance[wrss_POS] + deviance[usqr_POS])/
		 ((double)(dims[n_POS] -
			   (REML ? dims[p_POS] : 0))));
    return sigma;
}

/**
 * Evaluate the logarithm of the square of the determinant of L
 * (i.e. the logarithm of the determinant of LL')
 *
 * @param L
 *
 * @return logarithm of the square of the determinant of L
 */
static double
chm_log_det2(CHM_FR L)
{
    double ans = 0;
    int i;

    if (L->is_super) {
	for (i = 0; i < L->nsuper; i++) {
	    int j, nrp1 = 1 + ((int *)(L->pi))[i + 1] - ((int *)(L->pi))[i],
		nc = ((int *)(L->super))[i + 1] - ((int *)(L->super))[i];
	    double *x = (double *)(L->x) + ((int *)(L->px))[i];

	    for (j = 0; j < nc; j++) {
		ans += 2 * log(fabs(x[j * nrp1]));
	    }
	}
    } else {
	int *li = (int*)(L->i), *lp = (int*)(L->p), j, p;
	double *lx = (double *)(L->x);
	
	for (j = 0; j < L->n; j++) {
	    for (p = lp[j]; li[p] != j && p < lp[j + 1]; p++) {};
	    if (li[p] != j) break; /* what happened to the diagonal element? */
	    ans += log(lx[p] * ((L->is_ll) ? lx[p] : 1.));
	}
    }
    return ans;
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
 * b = T  %*% S %*% t(P) %*% u
 *
 * @param x an mer object
 *
 */
static void update_ranef(SEXP x)
{
    SEXP ST = GET_SLOT(x, lme4_STSym),
	permP = GET_SLOT(GET_SLOT(x, lme4_LSym), lme4_permSym);
    int *Gp = INTEGER(GET_SLOT(x, lme4_GpSym)),
	*perm = INTEGER(permP), i, k,
	m = LENGTH(permP), nf = LENGTH(ST);
    double *b = REAL(GET_SLOT(x, lme4_ranefSym)),
	*u = REAL(GET_SLOT(x, lme4_uvecSym)), one = 1;
    int *nc = Alloca(nf, int), *nlev = Alloca(nf, int);
    double **st = Alloca(nf, double*);
    R_CheckStack(); 

    ST_nc_nlev(ST, Gp, st, nc, nlev);
				/* inverse permutation */
    for (i = 0; i < m; i++) b[perm[i]] = u[i];
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
    CHM_SP cC = AS_CHM_SP(GET_SLOT(x, lme4_CmSym)),
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
    int w = asInteger(which),
	*dims = INTEGER(GET_SLOT(x, lme4_dimsSym));
		
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
		((int*)(rhs->i))[j] = iperm[Gp[i] + j * nlev[i]];
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

/**
 * Create and initialize L
 *
 * @param C pointer to the model matrix for the orthogonal random effects (transposed)
 *
 * @return L
 */
SEXP mer_create_L(SEXP C)
{
    double one[] = {1, 0};
    int fll = c.final_ll;
    CHM_SP cC = AS_CHM_SP(C);
    CHM_FR L;
    R_CheckStack();

    c.final_ll = 1;
    L = M_cholmod_analyze(cC, &c);
    if (!M_cholmod_factorize_p(cC, one, (int*)NULL, 0 /*fsize*/, L, &c))
	error(_("cholmod_factorize_p failed: status %d, minor %d from ncol %d"),
	      c.status, L->minor, L->n);
    c.final_ll = fll;

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
 * the current sqrtWt slot.  The sqrtWt slot is changed in update_L.
 *
 * @param x pointer to an mer object
 *
 * @return penalized, weighted residual sum of squares
 */
static double update_mu(SEXP x)
{
    SEXP X = GET_SLOT(x, lme4_XSym),
	moff = GET_SLOT(x, lme4_offsetSym),
	v = GET_SLOT(x, lme4_vSym);
    int *dims = INTEGER(GET_SLOT(x, lme4_dimsSym));
    int i, ione = 1, n = dims[n_POS], p = dims[p_POS];
    int nans = n * dims[s_POS];
    double *d = REAL(GET_SLOT(x, lme4_devianceSym)),
	*eta = REAL(GET_SLOT(x, lme4_etaSym)),
	*mu = REAL(GET_SLOT(x, lme4_muSym)),
	*muEta = SLOT_REAL_NULL(x, lme4_muEtaSym),
	*swt = SLOT_REAL_NULL(x, lme4_sqrtWtSym),
	*res = REAL(GET_SLOT(x, lme4_residSym)),
	*var = SLOT_REAL_NULL(x, lme4_varSym),
	*y = REAL(GET_SLOT(x, lme4_ySym)),
	one[] = {1,0};
    CHM_FR L = L_SLOT(x);
    CHM_SP cC = AS_CHM_SP(GET_SLOT(x, lme4_CmSym));
    CHM_DN Ptu, ceta = N_AS_CHM_DN(eta, nans, 1),
	cu = AS_CHM_DN(GET_SLOT(x, lme4_uvecSym));
    R_CheckStack();

				/* eta := offset or eta := 0 */
    if (LENGTH(moff)) Memcpy(eta, REAL(moff), nans);
    else AZERO(eta, nans);
				/* eta := eta + X \beta */
    F77_CALL(dgemv)("N", &nans, &p, one, REAL(X), &nans,
		    REAL(GET_SLOT(x, lme4_fixefSym)), &ione,
		    one, eta, &ione);
				/* eta := eta + V P' u */
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
	res[i] = (y[i] - mu[i]) * (swt ? swt[i] : 1);
	d[wrss_POS] += res[i] * res[i];
    }
				/* store u'u */
    d[usqr_POS] = sqr_length((double*)(cu->x), dims[q_POS]);
    
    return d[wrss_POS] + d[usqr_POS];
}

/**
 * Externally callable update_mu.
 * Update the eta, v, mu, resid and var slots according to the current
 * values of the parameters and u.  Also evaluate d[wrss_POS] using
 * the current sqrtWt slot.  The sqrtWt slot is changed in update_L.
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
 * Update the A, L, sqrtWt and resid slots.  It is assumed that
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
    int j, jv, lv = LENGTH(vP),
	n = dims[n_POS], p, s = dims[s_POS];
    double *d = REAL(GET_SLOT(x, lme4_devianceSym)),
	*res = REAL(GET_SLOT(x, lme4_residSym)),
	*mu = REAL(GET_SLOT(x, lme4_muSym)),
	*muEta = SLOT_REAL_NULL(x, lme4_muEtaSym),
	*swt = SLOT_REAL_NULL(x, lme4_sqrtWtSym),
	*y = REAL(GET_SLOT(x, lme4_ySym)),
	one[] = {1,0};
    CHM_SP C = AS_CHM_SP(GET_SLOT(x, lme4_CmSym));
    CHM_SP A = C;
    int *ai = (int*) NULL, *ap = (int*) NULL,
	*vi = (int*) C->i, *vp = (int*) C->p;
    double *ax = (double*) NULL, *vx = (double*) C->x;
    CHM_FR L = L_SLOT(x);
    R_CheckStack();

    if (lv || muEta || swt) {	/* must generate A from C */
	A = AS_CHM_SP(GET_SLOT(x, lme4_ASym));
	R_CheckStack();
	ai = (int*)(A->i); ap = (int*)(A->p); ax = (double*)(A->x);
	
	if (lv) { /* apply the gradient from the nonlinear model */
	    double *grad = REAL(getAttrib(vP, lme4_gradientSym));
	    
	    AZERO(ax, ap[n]);
	    for (jv = 0; jv < n * s; jv++) {
		int iv, ja = jv % n, ia;
		for (iv = vp[jv]; iv < vp[jv + 1]; iv++) {
		    for (ia = ap[ja]; ia < ap[ja + 1]; ia++)
			if (ai[ia] == vi[iv]) break;
		    if (ia == ap[ja + 1])
			error(_("Structure of A incompatible with C, jv = %d, iv = %d"),
			      jv, iv);
		    ax[ia] += grad[jv] * vx[iv];
		}
	    }
	} else {		/* copy contents of C to A */
	    Memcpy(ax, vx, ap[n]);
	}

	if (swt) {		/* update swt and res. Apply swt and,
				 * if used, muEta to A */
	    double *pwt = SLOT_REAL_NULL(x, lme4_priorWtSym),
		*var =  SLOT_REAL_NULL(x, lme4_varSym);
	    d[wrss_POS] = 0;
	    for (j = 0; j < n; j++) {
		swt[j] = sqrt((pwt ? pwt[j] : 1) / (var ? var[j] : 1));
		res[j] = swt[j] * (y[j] - mu[j]);
		d[wrss_POS] += res[j] * res[j];
		for (p = ap[j]; p < ap[j + 1]; p++)
		    ax[p] *= (muEta ? muEta[j] : 1) * swt[j];
	    }
	}
    }

    if (!M_cholmod_factorize_p(A, one, (int*)NULL, 0 /*fsize*/, L, &c))
	error(_("cholmod_factorize_p failed: status %d, minor %d from ncol %d"),
	      c.status, L->minor, L->n);
    d[ldL2_POS] = chm_log_det2(L);
    return d[wrss_POS] + d[usqr_POS];
}

/**
 * Externally callable version of update_L.
 * Update the A, L, sqrtWt and resid slots.  It is assumed that
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

/** Maximum number of iterations in cond_mode */
#define CM_MAXITER  300
/** Tolerance level for convergence criterion in cond_mode */
#define CM_TOL      1e-10
/** Minimum step factor in cond_mode */
#define CM_SMIN     1e-5

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
    double
	*res = REAL(GET_SLOT(x, lme4_residSym)), 
	*u = REAL(GET_SLOT(x, lme4_uvecSym)),
	cfac = ((double)n) / ((double)q), 
	crit, pwrss, pwrss_old, one[] = {1,0},
	step, zero[] = {0,0};
    double *tmp = Alloca(q, double),
	*uold = Alloca(q, double);
    CHM_FR L = L_SLOT(x);
    CHM_DN cres = N_AS_CHM_DN(res, n, 1),
	ctmp = N_AS_CHM_DN(tmp, q, 1), sol;
    CHM_SP A = AS_CHM_SP(GET_SLOT(x, lme4_ASym));
    R_CheckStack();

    if (!(L->is_ll)) error(_("L must be LL', not LDL'"));
    if (q > n) error(_("q = %d > n = %d"), q, n);

    /* resetting u to zero at the beginning of each evaluation
     * requires more iterations but is necessary to obtain a
     * repeatable evaluation.  If this is not done the optimization
     * algorithm can take wild steps. */

    AZERO(u, q);
    update_mu(x);
    for (i = 0; ; i++) {

	Memcpy(uold, u, q);
	pwrss_old = update_L(x);
				/* tmp := A %*% W^{1/2} %*% (y-mu) */
	if (!(M_cholmod_sdmult(A, 0 /* no trans */, one, zero,
			       cres, ctmp, &c))) 
	    error(_("cholmod_sdmult returned error code"));
				/* tmp := P %*% tmp - u */
	sol = M_cholmod_solve(CHOLMOD_P, L, ctmp, &c);
	for (j = 0; j < q; j++) tmp[j] = ((double*)(sol->x))[j] - u[j];
	M_cholmod_free_dense(&sol, &c);
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
 * Externally callable version of update_u.
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
 * Update the contents of the fixef, ranef and uvec slots in an lmer
 * object.
 *
 * @param x an lmer object
 *
 * @return R_NilValue
 */
SEXP mer_update_effects(SEXP x)
{
    SEXP muEta = GET_SLOT(x, lme4_muEtaSym),
	v = GET_SLOT(x, lme4_vSym);
    int *dims = INTEGER(GET_SLOT(x, lme4_dimsSym));
    CHM_FR L = L_SLOT(x);
    R_CheckStack();

    if (!(LENGTH(muEta) || LENGTH(v))) { /* linear mixed model */
	SEXP uvec = GET_SLOT(x, lme4_uvecSym);
	int ione = 1, p = dims[p_POS], pp1 = dims[p_POS] + 1,
	    q = dims[q_POS];
	double *RXy = REAL(GET_SLOT(x, lme4_RXySym)),
	    *RVXy = REAL(GET_SLOT(x, lme4_RVXySym)),
	    *dev = REAL(GET_SLOT(x, lme4_devianceSym)),
	    *fixef = REAL(GET_SLOT(x, lme4_fixefSym)),
	    *u = REAL(uvec), mone = -1, one = 1;
	CHM_DN cu = AS_CHM_DN(uvec), sol;
	R_CheckStack();

	Memcpy(fixef, RXy + p * pp1, p);
	F77_CALL(dtrsv)("U", "N", "N", &p, RXy, &pp1, fixef, &ione);
	Memcpy(u, RVXy + p * q, q);
	F77_CALL(dgemv)("N", &q, &p, &mone, RVXy, &q,
			fixef, &ione, &one, u, &ione);
	if (!(sol = M_cholmod_solve(CHOLMOD_Lt, L, cu, &c)))
	    error(_("cholmod_solve (CHOLMOD_Lt) failed:"));
	Memcpy(u, (double*)(sol->x), q);
	M_cholmod_free_dense(&sol, &c);
	dev[usqr_POS] = sqr_length(u, q);
	dev[disc_POS] = exp(dev[lpdisc_POS]) - dev[usqr_POS];
    }
    update_ranef(x);
    return R_NilValue;
}

/**
 * Create the C matrix pattern from Zt, ST and Gp.  Partition the
 * columns into s groups and overlay them.
 * 
 *
 * @param Zt 
 * @param ST
 * @param GpP pointer to the Gp array
 *
 * @return C
 */
SEXP mer_create_C(SEXP Zt, SEXP ST, SEXP GpP)
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

/**
 * Create PAW^{1/2}v in dest.
 *
 * @param dest values to be calculated
 * @param perm 0-based permutation that defines P
 * @param A 
 * @param swt pointer to diagonal elements of W^{1/2} or NULL
 * @param v
 * @param n length of v and number columns in A
 * @param q length of dest and number of rows in A
 *
 */
static void
PAWv(double *dest, const int *perm, const CHM_SP A,
      const double *swt, const double *v, int n, int q)
{
    int *ai = (int*)(A->i), *ap = (int*)(A->p), j, k, p;
    double *ax = (double*)(A->x), *tmp = Alloca(q, double);
    R_CheckStack();

    AZERO(tmp, q);
    for (j = 0; j < n; j++) {
	double mm = v[j] * (swt ? swt[j] : 1);
	for (p = ap[j]; p < ap[j + 1]; p++)
	    tmp[ai[p]] += mm * ax[p];
    }
    for (k = 0; k < q; k++)	/* apply P */
	dest[k] = tmp[perm ? perm[k] : k];
}

/**
 * Accumulate the weighted cross-product of [X:y] in rxy
 *
 * @param rxy pointer to a (p + 1) by (p + 1) matrix
 * @param swt pointer to diagonal elements of W^{1/2} or NULL
 * @param X pointer to the model matrix for the fixed effects
 * @param y pointer to the response vector
 * @param n number of rows in X
 * @param p number of columns in X
 *
 * \note Although the contents of X and y are not changed, we cannot
 * use the const modifier for the argument because the pointers are
 * copied.
 */
static void
weighted_XytXy(double *rxy, const double *swt, double *X,
	       double *y, int n, int p)
{
    int i, ione = 1, j, pp1 = p + 1;
    double one = 1, zero = 0;
    double *Xcp = X, *ycp = y;
    
    if (swt) {			/* created weighted copies of X and y */
	Xcp = Alloca(n * p, double),
	ycp = Alloca(n, double);
	R_CheckStack();
	
	for (i = 0; i < n; i++) ycp[i] = swt[i] * y[i];
	for (j = 0; j < p; j++)
	    for (i = 0; i < n; i++)
		Xcp[i + j * n] = swt[i] * X[i + j * n];
    }
    F77_CALL(dsyrk)("U", "T", &p, &n, &one, Xcp, &n, &zero, rxy, &pp1);
    F77_CALL(dgemv)("T", &n, &p, &one, Xcp, &n, ycp, &ione, &zero,
		    rxy + p * pp1, &ione);
    rxy[pp1 * pp1 - 1] = F77_CALL(ddot)(&n, ycp, &ione, ycp, &ione);
}
	
/**
 * Update the RCX and RX slots in an mer object
 *
 * @param x pointer to an mer object
 *
 * @return profiled deviance or REML deviance
 */
static double update_RX(SEXP x)
{
    SEXP RVXyP = GET_SLOT(x, lme4_RVXySym);
    int *dims = INTEGER(GET_SLOT(x, lme4_dimsSym));
    int j, n = dims[n_POS], p = dims[p_POS], q = dims[q_POS];
    int pp1 = p + 1;
    double *X = REAL(GET_SLOT(x, lme4_XSym)),
	*d = REAL(GET_SLOT(x, lme4_devianceSym)),
	*rvxy = REAL(RVXyP),
	*rxy = REAL(GET_SLOT(x, lme4_RXySym)),
	*swt = SLOT_REAL_NULL(x, lme4_sqrtWtSym),
	*y = REAL(GET_SLOT(x, lme4_ySym)),
    	dn = (double) n, dnmp = (double) (n - p),
	mone[] = {-1,0}, one[] = {1,0};
    CHM_SP A =
	AS_CHM_SP(GET_SLOT(x, swt ? lme4_ASym : lme4_CmSym));
    CHM_DN cRVXy = AS_CHM_DN(RVXyP), ans;
    CHM_FR L = L_SLOT(x);
    R_CheckStack();

    update_L(x);
				/* Evaluate PAW^{1/2}[X:y] in RVXy */
    for (j = 0; j < p; j++)
	PAWv(rvxy + j * q, (int*)(L->Perm), A, swt, X + j * n, n, q);
    PAWv(rvxy + p * q, (int*)(L->Perm), A, swt, y, n, q);
				/* solve L %*% ans = PAW^{1/2}[X:y] */
    ans = M_cholmod_solve(CHOLMOD_L, L, cRVXy, &c);
    Memcpy(rvxy, (double*)(ans->x), q * pp1);
    M_cholmod_free_dense(&ans, &c);
    				/* downdate XytXy  */
				/* Evaluate [X:y]'W[X:y] in RXy */
    weighted_XytXy(rxy, swt, X, y, n, p);
    F77_CALL(dsyrk)("U", "T", &pp1, &q, mone, rvxy, &q, one, rxy, &pp1);
    F77_CALL(dpotrf)("U", &pp1, rxy, &pp1, &j);
    if (j)
	error(_("Downdated [X:y]'[X:y] is not positive definite, %d."),
	      j);
    for (j = 0, d[ldRX2_POS] = 0; j < (pp1 - 1); j++)
	d[ldRX2_POS] += 2 * log(rxy[j * (pp1 + 1)]);
    /* update other elements of the deviance slot for LMMs */
    if ((!LENGTH(GET_SLOT(x, lme4_muEtaSym))) &&
	(!LENGTH(GET_SLOT(x, lme4_vSym)))) {
	d[lpdisc_POS] = 2 * log(rxy[pp1 * pp1 - 1]);
	d[ML_POS] = d[ldL2_POS] +
	    dn * (1 + d[lpdisc_POS] + log(2 * PI / dn));
	d[REML_POS] = d[ldL2_POS] + d[ldRX2_POS] +
	    dnmp * (1. + d[lpdisc_POS] + log(2. * PI / dnmp));
    }
    return d[dims[isREML_POS] ? REML_POS : ML_POS];
}

/**
 * Externally callable version of update_RX
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
/* FIXME: need to add 1 to nv for GLMM with a scale par (maybe). */
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
	if (lmm)
	    fx = update_RX(x);
	else {
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
 * Update the fixed effects and the random effects in an MCMC sample
 * from an lmer model.
 *
 * @param x an lmer model with mer_update_b
 * @param sigma current standard deviation of the per-observation
 *        noise terms.
 */
static void MCMC_betab(SEXP x, double sigma)
{
    SEXP fixefP = GET_SLOT(x, lme4_fixefSym),
	uvecP = GET_SLOT(x, lme4_uvecSym);
    int ione = 1, j, p = LENGTH(fixefP), q = LENGTH(uvecP);
    int pp1 = p + 1;
    double *fix = REAL(fixefP), *uv = REAL(uvecP), mone = -1, one = 1;
    CHM_FR L = L_SLOT(x);
    double *del1 = Alloca(q, double), *del2 = Alloca(p, double);
    CHM_DN sol, rhs = N_AS_CHM_DN(del1, q, 1);
    R_CheckStack();

    mer_update_effects(x);
				/* Update the fixed effects */
    for (j = 0; j < p; j++) del2[j] = sigma * norm_rand();
    F77_CALL(dtrsv)("U", "N", "N", &p, REAL(GET_SLOT(x, lme4_RXySym)),
		    &pp1, del2, &ione);
    for (j = 0; j < p; j++) fix[j] += del2[j];
				/* Update the orthogonal random effects */
    for (j = 0; j < q; j++) del1[j] = sigma * norm_rand();
    F77_CALL(dgemv)("N", &q, &p, &mone, REAL(GET_SLOT(x, lme4_RVXySym)),
		    &q, del2, &ione, &one, del1, &ione);
    sol = M_cholmod_solve(CHOLMOD_Lt, L, rhs, &c);
    for (j = 0; j < q; j++) uv[j] += ((double*)(sol->x))[j];
    M_cholmod_free_dense(&sol, &c);
				/* Update the random effects */
    update_ranef(x);
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
 * Generate a Markov-Chain Monte Carlo sample from a fitted
 * linear mixed model.
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
	n = dims[n_POS], nsamp = asInteger(nsampp),
	p = dims[p_POS], q = dims[q_POS],
	saveb = asLogical(savebp), trans = asLogical(transp),
	verbose = asLogical(verbosep), dev = asLogical(deviancep);
    double *ansp, *deviance = REAL(GET_SLOT(x, lme4_devianceSym)),
	*fixef = REAL(GET_SLOT(x, lme4_fixefSym)),
	*ranef = (double*)NULL, df = n - (dims[REML_POS] ? p : 0);
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
	sigma = exp(deviance[lpdisc_POS]/2)/sqrt(rchisq(df));
	col[p] = (trans ? 2 * log(sigma) : sigma * sigma);
	mer_update_effects(x);	/* Conditional estimates of beta, u and b */
	MCMC_betab(x, sigma); /* Updated beta, u and b */
	for (j = 0; j < p; j++) col[j] = fixef[j]; /* Store beta */
	if (saveb)				   /* store b */
	    for (j = 0; j < q; j++) col[nrbase + dev + j] = ranef[j];
	MCMC_ST(x, sigma, col + p + 1, trans);
/* 	lmer_update_dev(x);	/\* Refactor and evaluate deviance *\/ */
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

/**
 * Check the dimensions of the matrix pointer MP.
 *
 * @param buf character buffer of length nb + 1
 * @param nb number of writable positions in the buffer
 * @param nm name of matrix - used in the error message
 * @param MP pointer to a matrix
 * @param nr expected number of rows
 * @param nc expected number of columns
 *
 * @return error message or 0-length string in buffer
 */
static char
*chkDims(char *buf, int nb, char *nm, SEXP MP, int nr, int nc)
{
    int *dd = isMatrix(MP) ? INTEGER(getAttrib(MP, R_DimSymbol)) :
	(int *)NULL;

    if (!dd) error(_("Argument MP to chkDims is not a matrix"));
    buf[0] = '\0';
    if (!isReal(MP) || dd[0] != nr || dd[1] != nc)
	snprintf(buf, nb, "Matrix %s must be a %d by %d numeric matrix",
		nm, nr, nc);
    return buf;
}

#ifndef BUF_SIZE
/** size of buffer for an error message */
#define BUF_SIZE 127
#endif	

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
	XP = GET_SLOT(x, lme4_XSym),
	devianceP = GET_SLOT(x, lme4_devianceSym),
	dimsP = GET_SLOT(x, lme4_dimsSym),
	etaP = GET_SLOT(x, lme4_etaSym),
	fixefP = GET_SLOT(x, lme4_fixefSym),
	flistP = GET_SLOT(x, lme4_flistSym),
	muP = GET_SLOT(x, lme4_muSym),
	muEtaP = GET_SLOT(x, lme4_muEtaSym),
	offsetP = GET_SLOT(x, lme4_offsetSym),
	ranefP = GET_SLOT(x, lme4_ranefSym),
	residP = GET_SLOT(x, lme4_residSym),
	uvecP = GET_SLOT(x, lme4_uvecSym),
	vP = GET_SLOT(x, lme4_vSym),
	weightsP = GET_SLOT(x, lme4_priorWtSym),
	y = GET_SLOT(x, lme4_ySym);
    int *Gp = INTEGER(GpP), *dd = INTEGER(dimsP);
    int i, n = dd[n_POS], nf = dd[nf_POS], nq, p = dd[p_POS],
	q = dd[q_POS], s = dd[s_POS];
    int neta = n * s;
    CHM_SP Zt = AS_CHM_SP(GET_SLOT(x, lme4_ZtSym)),
	C =  AS_CHM_SP(GET_SLOT(x, lme4_CmSym));
    CHM_FR L = L_SLOT(x);
    char *buf = Alloca(BUF_SIZE + 1, char);
    R_CheckStack();
				/* check lengths */
    if (nf < 1 || LENGTH(flistP) != nf || LENGTH(ST) != nf)
	return mkString(_("Slots ST, and flist must have length dims['nf']"));
    if (LENGTH(GpP) != nf + 1)
	return mkString(_("Slot Gp must have length dims['nf'] + 1"));
    if (LENGTH(y) != n)
	return mkString(_("Slot y must have length dims['n']"));
    if (LENGTH(muP) != n)
	return mkString(_("Slot mu must have length dims['n']"));
    if (LENGTH(residP) != n)
	return mkString(_("Slot resid must have length dims['n']"));
    if (Gp[0] != 0 || Gp[nf] != q)
	return mkString(_("Gp[1] != 0 or Gp[dims['nf'] + 1] != dims['q']"));
    if (LENGTH(fixefP) != p)
	return mkString(_("Slot fixef must have length ['p']"));
    if (LENGTH(ranefP) != q)
	return mkString(_("Slot ranef must have length dims['q']"));
    if (LENGTH(uvecP) != q)
	return mkString(_("Slot uvec must have length dims['q']"));
    if (LENGTH(etaP) != neta)
	return mkString(_("Slot eta must have length dims['n']*dims['s']"));
    if (LENGTH(weightsP) && LENGTH(weightsP) != n)
	return mkString(_("Slot weights must have length 0 or dims['n']"));
    if (LENGTH(offsetP) && LENGTH(offsetP) != n)
	return mkString(_("Slot offset must have length 0 or dims['n']*dims['s']"));
    if (LENGTH(devianceP) != (wrss_POS + 1) ||
	LENGTH(getAttrib(devianceP, R_NamesSymbol)) != (wrss_POS + 1))
	return mkString(_("deviance slot not named or incorrect length"));
    if (LENGTH(dimsP) != (cvg_POS + 1) ||
	LENGTH(getAttrib(dimsP, R_NamesSymbol)) != (cvg_POS + 1))
	return mkString(_("dims slot not named or incorrect length"));
    if (L->n != q || !L->is_ll || !L->is_monotonic)
	return mkString(_("Slot L must be a monotonic LL' factorization of size dims['q']"));
    if (Zt->nrow != q || Zt->ncol != neta)
	return mkString(_("Slot Zt must by dims['q']  by dims['n']*dims['s']"));
    if (C->nrow != q || C->ncol != neta)
	return mkString(_("Slot Zt must be dims['q']  by dims['n']*dims['s']"));
    dd = INTEGER(getAttrib(XP, R_DimSymbol)); 
    if (!isReal(XP) || dd[1] != p)
	return mkString(_("Slot X must be a numeric matrix with dims['p'] columns"));
    if (dd[0] != 0 && dd[0] != (n * s)) /* special case */
	return mkString(_("Slot X must be have 0 or dims['n']*dims['s'] rows"));

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
    if (strlen(chkDims(buf, (int) BUF_SIZE, "RVXy",
		       GET_SLOT(x, lme4_RVXySym), q, p + 1)))
	return mkString(buf);
    if (strlen(chkDims(buf, (int) BUF_SIZE, "RXy",
		       GET_SLOT(x, lme4_RXySym), p + 1, p + 1)))
	return mkString(buf);
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
