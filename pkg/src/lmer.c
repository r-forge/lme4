#include "lmer.h"
#include <Rmath.h>		/* for dnorm5, etc. */
#include <R_ext/Lapack.h>     /* for Lapack (dpotrf, etc.) and BLAS */
#include <R_ext/stats_package.h> /* for S_nlminb_iterate */
#include "Matrix.h"		 /* for cholmod functions */

extern
#include "Syms.h"
extern	       /** cholmod_common struct initialized in R_init_lme4 */
cholmod_common c;

#ifdef ENABLE_NLS		/** Allow for translation of error messages */
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
    sigmaML_POS,		/**<current ML estimate of sigma */
    sigmaREML_POS,		/**<current REML estimate of sigma */
    pwrss_POS,			/**<penalized weighted residual sum of squares */
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
    vTyp_POS,			/**<variance type for generalized model */
    nest_POS,			/**<indicator of nested grouping factors */
    useSc_POS,			/**<does the family use a separate scale parameter */
    cvg_POS			/**<convergence indictor from port optimization  */
};

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

/** Allocate (alloca) a cholmod_sparse struct, populate it with values
 * from the A slot and return the pointer. */
#define A_SLOT(x) AS_CHM_SP(GET_SLOT(x, lme4_ASym))

/** Allocate (alloca) a cholmod_sparse struct, populate it with values
 * from the Cm slot and return the pointer. */
#define Cm_SLOT(x) AS_CHM_SP(GET_SLOT(x, lme4_CmSym))

/** Return the double pointer to the Cx slot or (double*) NULL if
 * Cx has length 0) */
#define Cx_SLOT(x) SLOT_REAL_NULL(x, lme4_CxSym)

/** Return the double pointer to the deviance slot */
#define DEV_SLOT(x) SLOT_REAL_NULL(x, lme4_devianceSym)

/** Return the integer pointer to the dims slot */
#define DIMS_SLOT(x) INTEGER(GET_SLOT(x, lme4_dimsSym))

/** Return the double pointer to the eta slot */
#define ETA_SLOT(x) SLOT_REAL_NULL(x, lme4_etaSym)

/** Return the double pointer to the fixef slot */
#define FIXEF_SLOT(x) SLOT_REAL_NULL(x, lme4_fixefSym)

/** Return the integer pointer to the Gp slot */
#define Gp_SLOT(x) INTEGER(GET_SLOT(x, lme4_GpSym))

/** Allocate (alloca) a cholmod_factor struct, populate it with values
 * from the L slot and return the pointer. */
#define L_SLOT(x) AS_CHM_FR(GET_SLOT(x, lme4_LSym))

/** Return the integer pointer to method used in likelihood approximation. */
#define METHOD_SLOT(x) INTEGER(GET_SLOT(x, lme4_dimsSym))

/** Return the double pointer to the mu slot */
#define MU_SLOT(x) SLOT_REAL_NULL(x, lme4_muSym)

/** Return the double pointer to the muEta slot or (double*) NULL if
 * muEta has length 0) */
#define MUETA_SLOT(x) SLOT_REAL_NULL(x, lme4_muEtaSym)

/** Return the double pointer to the offset slot or (double*) NULL if
 * offset has length 0) */
#define OFFSET_SLOT(x) SLOT_REAL_NULL(x, lme4_offsetSym)

/** Return the integer pointer to the permutation vector in the L slot */
#define PERM_VEC(x) INTEGER(GET_SLOT(GET_SLOT(x, lme4_LSym), lme4_permSym))

/** Return the double pointer to the pWt slot or (double*) NULL if
 * pWt has length 0) */
#define PWT_SLOT(x) SLOT_REAL_NULL(x, lme4_pWtSym)

/** Return the double pointer to the ranef slot or (double*) NULL if
 *  ranef has length 0) */
#define RANEF_SLOT(x) SLOT_REAL_NULL(x, lme4_ranefSym)

/** Residual degrees of freedom */
#define RDF(dims) (dims[n_POS] - (dims[isREML_POS] ? dims[p_POS] : 0))

/** Return the double pointer to the resid slot */
#define RESID_SLOT(x) SLOT_REAL_NULL(x, lme4_residSym)

/** Return the double pointer to the RX slot */
#define RX_SLOT(x) SLOT_REAL_NULL(x, lme4_RXSym)

/** Return the double pointer to the RZX slot */
#define RZX_SLOT(x) SLOT_REAL_NULL(x, lme4_RZXSym)

/** Return the double pointer to the sqrtrWt slot or (double*) NULL if
 *  sqrtrWt has length 0) */
#define SRWT_SLOT(x) SLOT_REAL_NULL(x, lme4_sqrtrWtSym)

/** Return the double pointer to the sqrtXWt slot or (double*) NULL if
 *  sqrtXWt has length 0) */
#define SXWT_SLOT(x) SLOT_REAL_NULL(x, lme4_sqrtXWtSym)

/** Return the double pointer to the u slot */
#define U_SLOT(x) SLOT_REAL_NULL(x, lme4_uSym)

/** Return the double pointer to the V slot or (double*) NULL if
 * V has length 0) */
#define V_SLOT(x) SLOT_REAL_NULL(x, lme4_VSym)

/** Return the double pointer to the var slot or (double*) NULL if
 * var has length 0) */
#define VAR_SLOT(x) SLOT_REAL_NULL(x, lme4_varSym)

/** Return the double pointer to the X slot */
#define X_SLOT(x) SLOT_REAL_NULL(x, lme4_XSym)

/** Return the double pointer to the y slot */
#define Y_SLOT(x) SLOT_REAL_NULL(x, lme4_ySym)

/** Allocate (alloca) a cholmod_sparse struct, populate it with values
 * from the Zt slot and return the pointer. */
#define Zt_SLOT(x) AS_CHM_SP(GET_SLOT(x, lme4_ZtSym))

/* Constants */

#ifndef BUF_SIZE
/** size of buffer for an error message */
#define BUF_SIZE 127
#endif	

/** Maximum number of iterations in update_u */
#define CM_MAXITER  300
/** Tolerance level for convergence criterion in update_u */
#define CM_TOL      1e-10
/** Minimum step factor in update_u */
#define CM_SMIN     1e-5

#define LTHRESH     30.
#define MLTHRESH   -30.

static double MPTHRESH = 0;
static double PTHRESH = 0;
static const double INVEPS = 1/DOUBLE_EPS;


/* In-line functions */

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
 * Return the index of the term associated with parameter index ind
 *
 * @param ind an index in [0, Gp[nf] - 1]
 * @param nf total number of terms
 * @param Gp group pointers, a vector of length nf+1 with Gp[0] = 0
 *
 * @return sum of squares
 */
static R_INLINE int Gp_grp(int ind, int nf, const int *Gp)
{
    for (int i = 0; i < nf; i++) if (ind < Gp[i + 1]) return i;
    error(_("invalid row index %d (max is %d)"), ind, Gp[nf]);
    return -1;                  /* -Wall */
}

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
    for (int i = 0; i < n; i++) ans += x[i] * x[i];
    return ans;
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

/* Stand-alone utility functions (sorted by function name) */

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
 * Evaluate the sum of the deviance residuals for a GLM family
 *
 * @param mu pointer to the mu vector
 * @param pWt pointer to the vector of prior weights (NULL for
 *            constant weights)
 * @param y pointer to the response vector
 * @param n length of mu and y
 * @param vTyp type of variance function: the 1-based index into
 *        c("constant", "mu(1-mu)", "mu", "mu^2", "mu^3")
 *
 * @return the sum of the deviance residuals
 */
static double
lme4_devResid(const double* mu, const double* pWt, const double* y,
	      int n, int vTyp)
{
    double ans = 0.;
    for (int i = 0; i < n; i++) {
	double mui = mu[i], wi = pWt ? pWt[i] : 1, yi = y[i];
	double ri = yi - mui;
	switch(vTyp) {
	case 1:			/* constant variance */
	    ans += wi * ri * ri;
	    break;
	case 2:			/* mu(1-mu) variance */
	    ans += 2 * wi *
		(y_log_y(yi, mui) + y_log_y(1 - yi, 1 - mui));
	    break;
	case 3:			/* mu variance */
	    ans += 2 * wi * (y_log_y(yi, mui) - (yi - mui));
	    break;
	case 4:			/* mu^2 variance */
	    ans += 2 * wi * (y_log_y(yi, mui) - (yi - mui)/mui);
	    break;
	case 5:			/* mu^3 variance */
	    ans += wi * (ri * ri)/(yi * mui * mui);
	    break;
	default:
	    error(_("Unknown vTyp value %d"), vTyp);
	}
    }
    return ans;
}

/**
 * Evaluate the derivative d mu/d eta for the GLM link function of
 * type lTyp
 *
 * @param mu pointer to the mu vector
 * @param muEta pointer to the muEta vector
 * @param eta pointer to the eta vector
 * @param n length of mu, muEta and eta
 * @param lTyp type of link: the 1-based index into
 *        c("logit", "probit", "cauchit", "cloglog", "identity",
 *          "log", "sqrt", "1/mu^2", "inverse")
 */
static void
lme4_muEta(double* mu, double* muEta, const double* eta, int n, int lTyp)
{
    for (int i = 0; i < n; i++) { /* apply the generalized linear part */
	double etai = eta[i], tmp;
	switch(lTyp) {
	case 1:		/* logit */
	    tmp = (etai < MLTHRESH) ? DOUBLE_EPS :
	    ((etai > LTHRESH) ? INVEPS : exp(etai));
	    mu[i] = tmp/(1 + tmp);
	    muEta[i] = mu[i] * (1 - mu[i]);
	    break;
	case 2:		/* probit */
	    if (!MPTHRESH) {
		MPTHRESH = qnorm5(DOUBLE_EPS, 0, 1, 1, 0);
		PTHRESH = -MPTHRESH;
	    }
	    mu[i] = (etai < MPTHRESH) ? DOUBLE_EPS :
		((etai > PTHRESH) ? 1 - DOUBLE_EPS :
		 pnorm5(etai, 0, 1, 1, 0));
	    tmp = dnorm4(eta[i], 0, 1, 0);
	    muEta[i] = (tmp < DOUBLE_EPS) ? DOUBLE_EPS : tmp;
	    break;
	case 3:		/* cauchit */
	    error(_("cauchit link not yet coded"));
	    break;
	case 4:		/* cloglog */
	    error(_("cloglog link not yet coded"));
	    break;
	case 5:		/* identity */
	    mu[i] = eta[i];
	    muEta[i] = 1.;
	    break;
	case 6:		/* log */
	    tmp = exp(eta[i]);
	    muEta[i] = mu[i] =
		(tmp < DOUBLE_EPS) ? DOUBLE_EPS : tmp;
	    break;
	case 7:		/* sqrt */
	    mu[i] = etai * etai;
	    muEta[i] = 2 * etai;
	    break;
	default:
	    error(_("General form of glmer_linkinv not yet written"));
	}
    }
}

/**
 * Evaluate the GLM variance function of type vTyp given mu
 *
 * @param var pointer to the variance vector
 * @param mu pointer to the mu vector
 * @param n length of var and mu
 * @param vTyp type of variance function: the 1-based index into
 *        c("constant", "mu(1-mu)", "mu", "mu^2", "mu^3")
 *
 */
static void
lme4_varFunc(double* var, const double* mu, int n, int vTyp)
{
    for (int i = 0; i < n; i++) {
	double mui = mu[i];
	switch(vTyp) {
	case 1:			/* constant variance */
	    var[i] = 1.;
	    break;
	case 2:			/* mu(1-mu) variance */
	    if (mui <= 0 || mui >= 1)
		error(_("mu[i] must be in the range (0,1): mu = %g, i = %d"),
		      mu, i);
	    var[i] = mui * (1 - mui);
	    break;
	case 3:			/* mu variance */
	    if (mui <= 0)
		error(_("mu[i] must be positive: mu = %g, i = %d"), mu, i);
	    var[i] = mui;
	    break;
	case 4:			/* mu^2 variance */
	    if (mui <= 0)
		error(_("mu[i] must be positive: mu = %g, i = %d"), mu, i);
	    var[i] = mui * mui;
	    break;
	case 5:			/* mu^3 variance */
	    if (mui <= 0)
		error(_("mu[i] must be positive: mu = %g, i = %d"), mu, i);
	    var[i] = mui * mui * mui;
	    break;
	default:
	    error(_("Unknown vTyp value %d"), vTyp);
	}
    }
}
	
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
    int nf = LENGTH(ST), pos = 0;
    for (int i = 0; i < nf; i++) {
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
    int ans = 0, nf = LENGTH(ST);

    for (int i = 0; i < nf; i++) {
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
 * Determine the index in theta_S corresponding to index i in the u vector
 *
 * @param i index in u vector (0 <= i < q)
 * @param nf number of random effects terms in the model
 * @param Gp vector of group pointers into u (length nf + 1, Gp[0] == 0)
 * @param nlev vector of length nf giving the number of levels per term
 * @param spt vector of group pointers into theta_S (length nf + 1,
 *            spt[0] == 0)
 */

static R_INLINE int
theta_S_ind(int i, int nf, int *Gp, int *nlev, int *spt)
{
	int trm = Gp_grp(i, nf, Gp);
	return (spt[trm] + (i - Gp[trm]) / nlev[trm]);
}

/**
 * Multiply A on the left by T'
 *
 * @param A sparse model matrix
 * @param Gp group pointers
 * @param nc number of columns per term
 * @param nlev number of levels per term
 * @param st ST arrays for each term
 * @param nf number of terms
 * 
 */
static void Tt_Zt(CHM_SP A, int *Gp, int *nc, int *nlev, double **st, int nf)
{
    int *ai = (int*)(A->i), *ap = (int *)(A->p);
    double *ax = (double*)(A->x), one[] = {1,0};

    for (int j = 0; j < A->ncol; j++) /* multiply column j by T' */
	for (int p = ap[j]; p < ap[j + 1];) {
	    int i = Gp_grp(ai[p], nf, Gp);
	    
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

/* Level-1 utilties that call at least one of the stand-alone utilities */

/**
 * Determine the conditional estimates of the fixed effects
 * and the conditional mode of u for a linear mixed model
 *
 * @param x an mer object
 */
static double lmm_update_fixef_u(SEXP x)
{
    double ans = NA_REAL;
    if (!V_SLOT(x) && !MUETA_SLOT(x)) { /* linear mixed model */
	int *dims = DIMS_SLOT(x), i1 = 1;
	int n = dims[n_POS], p = dims[p_POS], q = dims[q_POS];
	double *WX = (double*) NULL, *X = X_SLOT(x),
	    *d = DEV_SLOT(x), *fixef = FIXEF_SLOT(x),
	    *RZX = RZX_SLOT(x), *RX = RX_SLOT(x),
	    *sXwt = SXWT_SLOT(x), *srwt = SRWT_SLOT(x),
	    *u = U_SLOT(x), *wy = (double*)NULL, *y = Y_SLOT(x),
	    dn = (double) n, dnmp = (double) (n - p),
	    mone[] = {-1,0}, one[] = {1,0}, zero[] = {0,0};
	CHM_SP A = A_SLOT(x);
	CHM_FR L = L_SLOT(x);
	CHM_DN cu = N_AS_CHM_DN(u, q, 1), sol;
	R_CheckStack();
	
	if (sXwt) {		/* Replace X and y by weighted X and y */
	    WX = Calloc(n * p, double);
	    wy = Calloc(n, double);
	    
	    AZERO(WX, n * p);
	    for (int i = 0; i < n; i++) {
		wy[i] = sXwt[i] * y[i];
		for (int j = 0; j < p; j++)
		    WX[i + j * n] += sXwt[i] * X[i + j * n];
	    }
	    X = WX;
	    y = wy;
	}
				/* solve L del1 = PAy */
	P_sdmult(u, (int*)L->Perm, A, y, 1);
	sol = M_cholmod_solve(CHOLMOD_L, L, cu, &c);
	Memcpy(u, (double*)sol->x, q);
	M_cholmod_free_dense(&sol, &c);
				/* solve RX' del2 = X'y - RZX'del1 */
	F77_CALL(dgemv)("T", &n, &p, one, X, &n,
			y, &i1, zero, fixef, &i1);
	F77_CALL(dgemv)("T", &q, &p, mone, RZX, &q,
			u, &i1, one, fixef, &i1);
	F77_CALL(dtrsv)("U", "T", "N", &p, RX, &p, fixef, &i1);
	d[pwrss_POS] = sqr_length(y, n)
	    - (sqr_length(fixef, p) + sqr_length(u, q));
	if (d[pwrss_POS] < 0)
	    error(_("Calculated PWRSS for a LMM is negative"));
				/* solve RX beta-hat = del2 */
	F77_CALL(dtrsv)("U", "N", "N", &p, RX_SLOT(x), &p, fixef, &i1);
	F77_CALL(dgemv)("N", &q, &p, mone, RZX_SLOT(x), &q,
			fixef, &i1, one, u, &i1);
	if (!(sol = M_cholmod_solve(CHOLMOD_Lt, L, cu, &c)))
	    error(_("cholmod_solve (CHOLMOD_Lt) failed:"));
	Memcpy(u, (double*)(sol->x), q);
	M_cholmod_free_dense(&sol, &c);

	d[usqr_POS] = sqr_length(u, q);
	d[wrss_POS] = d[disc_POS] = d[pwrss_POS] - d[usqr_POS];
	d[ML_POS] = d[ldL2_POS] +
	    dn * (1 + log(d[pwrss_POS]) + log(2 * PI / dn));
	d[REML_POS] = d[ldL2_POS] + d[ldRX2_POS] +
	    dnmp * (1. + log(d[pwrss_POS]) + log(2. * PI / dnmp));
	d[sigmaML_POS] = sqrt(d[pwrss_POS]/
			      (srwt ? sqr_length(srwt, n) : dn));
	d[sigmaREML_POS] = d[sigmaML_POS] * sqrt(dn/dnmp);

	if (wy) Free(wy);
	if (WX) Free(WX);
	ans = d[dims[isREML_POS] ? REML_POS : ML_POS];
    }
    return ans;
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
	nf = DIMS_SLOT(x)[nf_POS];
    int annz = ap[A->ncol], znnz = zp[Zt->ncol];
    int *nc = Alloca(nf, int), *nlev = Alloca(nf, int);
    double **st = Alloca(nf, double*), *ax = (double*)(A->x),
	*zx = (double*)(Zt->x);
    R_CheckStack();

    ncmax = ST_nc_nlev(GET_SLOT(x, lme4_STSym), Gp, st, nc, nlev);

    if (annz == znnz) { /* Copy Z' to A unless A has new nonzeros */
	Memcpy(ax, zx, znnz);
    } else { /* Only for nonlinear models with correlated random effects */
	AZERO(ax, annz); 	/* Initialize potential nonzeros to 0 */
	for (int j = 0; j < A->ncol; j++) { /* Iterate over columns */
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
    if (ncmax > 1) Tt_Zt(A, Gp, nc, nlev, st, nf);
				/* Multiply A on the left by S */
    for (int p = 0; p < annz; p++) {
	int i = Gp_grp(ai[p], nf, Gp);
	ax[p] *= st[i][((ai[p] - Gp[i]) / nlev[i]) * (nc[i] + 1)];
    }
}

/**
 * Update the L, sqrtrWt and resid slots.  It is assumed that
 * update_mu has already been called at the current values of u and
 * the model parameters and that A has been updated.
 *
 * @param x pointer to an mer object
 *
 * @return penalized weighted residual sum of squares
 *
 */
static double update_L(SEXP x)
{
    int *dims = DIMS_SLOT(x);
    int n = dims[n_POS], p = dims[p_POS], s = dims[s_POS];
    double *V = V_SLOT(x), *cx = Cx_SLOT(x), *d = DEV_SLOT(x),
	*res = RESID_SLOT(x), *mu = MU_SLOT(x), *muEta = MUETA_SLOT(x),
	*pwt = PWT_SLOT(x), *sXwt = SXWT_SLOT(x), *srwt = SRWT_SLOT(x),
	*var =  VAR_SLOT(x), *y = Y_SLOT(x), one[] = {1,0};
    CHM_SP A = A_SLOT(x);
    CHM_FR L = L_SLOT(x);
    R_CheckStack();
    
    if (var || pwt) {	   /* Update srwt and res. Reevaluate wrss. */
	d[wrss_POS] = 0;
	for (int j = 0; j < n; j++) {
	    srwt[j] = sqrt((pwt ? pwt[j] : 1.0) / (var ? var[j] : 1.0));
	    res[j] = srwt[j] * (y[j] - mu[j]);
	    d[wrss_POS] += res[j] * res[j];
	}
    }
    if (sXwt) {			/* Update sXwt and C */
	int *ai = (int*)A->i, *ap = (int*)A->p;
	double *ax = (double*)(A->x);
	CHM_SP C = A;

	for (int j = 0; j < s; j++) {
	    for (int i = 0; i < n; i++) {
		int ja = i + j * n;
		sXwt[ja] = (srwt ? srwt[i] : 1) *
		    (muEta ? muEta[i] : 1) * (V ? V[ja] : 1);
	    }
	}
	if (cx) {		/* C is a scaled version of A */
	    for (int j = 0; j < n; j++)
		for (int p = ap[j]; p < ap[j + 1]; p++)
		    cx[p] = ax[p] * sXwt[j];
	    A->x = (void*)cx;
	} else {
	    int *ci, *cp;
	    C = Cm_SLOT(x);
	    R_CheckStack();

	    ci = (int*)C->i; cp = (int*)C->p; cx = (double*)C->x;
	    AZERO(cx, cp[n]);
	    for (int j = 0; j < s; j++)
		for (int i = 0; i < n; i++) {
		    int ja = i + j * n, pc;
		    for (int pa = ap[ja]; pa < ap[ja + 1]; pa++) {
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

    d[ldL2_POS] = M_chm_factor_ldetL2(L);
    d[pwrss_POS] = d[usqr_POS] + d[wrss_POS];
    d[sigmaML_POS] = sqrt(d[pwrss_POS]/
			  (srwt ? sqr_length(srwt, n) : (double) n));
    d[sigmaREML_POS] = (V || muEta) ? NA_REAL :
	d[sigmaML_POS] * sqrt((((double) n)/((double)(n - p))));
    return d[pwrss_POS];
}

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
    int *dims = DIMS_SLOT(x);
    int i1 = 1, n = dims[n_POS], p = dims[p_POS], s = dims[s_POS];
    int ns = n * s;
    double *V = V_SLOT(x), *d = DEV_SLOT(x), *eta = ETA_SLOT(x),
	*etaold = (double*) NULL, *mu = MU_SLOT(x),
	*muEta = MUETA_SLOT(x), *offset = OFFSET_SLOT(x),
	*srwt = SRWT_SLOT(x), *res = RESID_SLOT(x),
	*var = VAR_SLOT(x), *y = Y_SLOT(x), one[] = {1,0};
    CHM_FR L = L_SLOT(x);
    CHM_SP A = A_SLOT(x);
    CHM_DN Ptu, ceta, cu = AS_CHM_DN(GET_SLOT(x, lme4_uSym));
    R_CheckStack();

    if (V) {
	etaold = eta;
	eta = Calloc(ns, double);
    }
				/* eta := offset or eta := 0 */
    for (int i = 0; i < ns; i++) eta[i] = offset ? offset[i] : 0;
				/* eta := eta + X \beta */
    F77_CALL(dgemv)("N", &ns, &p, one, X_SLOT(x), &ns,
		    FIXEF_SLOT(x), &i1, one, eta, &i1);
				/* eta := eta + C' P' u */
    Ptu = M_cholmod_solve(CHOLMOD_Pt, L, cu, &c);
    ceta = N_AS_CHM_DN(eta, ns, 1);
    R_CheckStack();
    if (!M_cholmod_sdmult(A, 1 /* trans */, one, one, Ptu, ceta, &c))
	error(_("cholmod_sdmult error returned"));
    M_cholmod_free_dense(&Ptu, &c);

    if (V) {		/* evaluate the nonlinear model */
	SEXP pnames = VECTOR_ELT(GET_DIMNAMES(GET_SLOT(x, lme4_VSym)), 1),
	    gg, rho = GET_SLOT(x, lme4_envSym), vv;
	int *gdims;
	
	if (!isString(pnames) || LENGTH(pnames) != s)
	    error(_("Slot V must be a matrix with %d named columns"), s);
	for (int i = 0; i < s; i++) { /* par. vals. into env. */
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
	Free(eta);
	eta = etaold;
	Memcpy(eta, REAL(vv), n);
	Memcpy(V, REAL(gg), ns);
	UNPROTECT(1);
    }

    if (muEta) {
	lme4_muEta(mu, muEta, eta, n, dims[lTyp_POS]);
	lme4_varFunc(var, mu, n, dims[vTyp_POS]);
    } else {
	Memcpy(mu, eta, n);
    }

    d[wrss_POS] = 0;		/* update resid slot and d[wrss_POS] */
    for (int i = 0; i < n; i++) {
	res[i] = (y[i] - mu[i]) * (srwt ? srwt[i] : 1);
	d[wrss_POS] += res[i] * res[i];
    }
				/* store u'u */
    d[usqr_POS] = sqr_length((double*)(cu->x), dims[q_POS]);
    d[pwrss_POS] = d[usqr_POS] + d[wrss_POS];
    d[sigmaML_POS] = sqrt(d[pwrss_POS]/
			  (srwt ? sqr_length(srwt, n) : (double) n));
    d[sigmaREML_POS] = (V || muEta) ? NA_REAL :
	d[sigmaML_POS] * sqrt((((double) n)/((double)(n - p))));
    return d[pwrss_POS];
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
    int *Gp = Gp_SLOT(x), *dims = DIMS_SLOT(x), *perm = PERM_VEC(x);
    int nf = dims[nf_POS], q = dims[q_POS];
    double *b = RANEF_SLOT(x), *u = U_SLOT(x), one[] = {1,0};
    int *nc = Alloca(nf, int), *nlev = Alloca(nf, int);
    double **st = Alloca(nf, double*);
    R_CheckStack(); 

    ST_nc_nlev(GET_SLOT(x, lme4_STSym), Gp, st, nc, nlev);
				/* inverse permutation */
    for (int i = 0; i < q; i++) b[perm[i]] = u[i];
    for (int i = 0; i < nf; i++) {
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
 * Update the RZX and RX slots in an mer object. update_L should be
 * called before update_RX
 *
 * @param x pointer to an mer object
 *
 * @return profiled deviance or REML deviance
 */
static double update_RX(SEXP x)
{
    int *dims = DIMS_SLOT(x), info;
    int n = dims[n_POS], p = dims[p_POS], q = dims[q_POS], s = dims[s_POS];
    double *cx = Cx_SLOT(x), *d = DEV_SLOT(x),
	*RZX = RZX_SLOT(x), *RX = RX_SLOT(x), *sXwt = SXWT_SLOT(x),
	*WX = (double*) NULL, *X = X_SLOT(x),
	mone[] = {-1,0}, one[] = {1,0}, zero[] = {0,0};
    CHM_SP A = A_SLOT(x);
    CHM_FR L = L_SLOT(x);
    CHM_DN cRZX = N_AS_CHM_DN(RZX, q, p), ans;
    R_CheckStack();

    if (sXwt) {			/* Create W^{1/2}GHX in WX */
	WX = Calloc(n * p, double);

	AZERO(WX, n * p);
	for (int j = 0; j < p; j++)
	    for (int k = 0; k < s; k++)
		for (int i = 0; i < n; i++)
		    WX[i + j * n] +=
			sXwt[i + k * n] * X[i + n * (k + j * s)];
	X = WX;
	/* Replace A by C, either just the x component or the entire matrix */
	if (cx) A->x = (void*)cx;
	else {
	    A = Cm_SLOT(x);
	    R_CheckStack();
	}
    }
				/* solve L %*% RZX = PAW^{1/2}GHX */
    P_sdmult(RZX, (int*)L->Perm, A, X, p); /* right-hand side */
    ans = M_cholmod_solve(CHOLMOD_L, L, cRZX, &c); /* solution */
    Memcpy(RZX, (double*)(ans->x), q * p);
    M_cholmod_free_dense(&ans, &c);
    				/* downdate X'X and factor  */
    F77_CALL(dsyrk)("U", "T", &p, &n, one, X, &n, zero, RX, &p); /* X'X */
    F77_CALL(dsyrk)("U", "T", &p, &q, mone, RZX, &q, one, RX, &p);
    F77_CALL(dpotrf)("U", &p, RX, &p, &info);
    if (info)
	error(_("Downdated X'X is not positive definite, %d."), info);
				/* accumulate log(det(RX)^2)  */
    d[ldRX2_POS] = 0;
    for (int j = 0; j < p; j++) d[ldRX2_POS] += 2 * log(RX[j * (p + 1)]);

    if (WX) Free(WX);
    return d[ML_POS];
}

/* Level-2 utilities that call at least one of the level-1 utilities */

/**
 * Update the ST and A slots of an mer object.
 *
 * @param x an mer object
 * @param pars double vector of the appropriate length
 *
 */
static void
ST_setPars(SEXP x, const double *pars)
{
    int *Gp = Gp_SLOT(x), nf = DIMS_SLOT(x)[nf_POS], pos = 0;
    int *nc = Alloca(nf, int), *nlev = Alloca(nf, int);
    double **st = Alloca(nf, double*);
    R_CheckStack();

    ST_nc_nlev(GET_SLOT(x, lme4_STSym), Gp, st, nc, nlev);
				/* install the parameters in the ST slot */
    for (int i = 0; i < nf; i++) {
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
    int *dims = DIMS_SLOT(x);
    int i, n = dims[n_POS], q = dims[q_POS];
    double *Cx = Cx_SLOT(x), *sXwt = SXWT_SLOT(x),
	*res = RESID_SLOT(x), *u = U_SLOT(x),
	cfac = ((double)n) / ((double)q), 
	crit, pwrss, pwrss_old, step;
    double *tmp = Alloca(q, double), *tmp1 = Alloca(q, double),
	*uold = Alloca(q, double), one[] = {1,0}, zero[] = {0,0};
    CHM_FR L = L_SLOT(x);
    CHM_DN cres = N_AS_CHM_DN(res, n, 1),
	ctmp = N_AS_CHM_DN(tmp, q, 1), sol;
    CHM_SP C = Cm_SLOT(x);
    R_CheckStack();
    
    if (!sXwt) return(0);	/* nothing to do for LMMs */
    if (!(L->is_ll)) error(_("L must be LL', not LDL'"));
    if (q > n) error(_("q = %d > n = %d"), q, n);
    if (Cx) {		    /* A and C have the same structure */
	C = A_SLOT(x);
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
	for (int j = 0; j < q; j++) tmp[j] -= u[j];
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
	    for (int j = 0; j < q; j++) u[j] = uold[j] + step * tmp[j];
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

/* Level-3 utilities that call at least one level-2 utilities */

/**
 * Update the fixed effects and the orthogonal random effects in an MCMC sample
 * from an mer object.
 *
 * @param x an mer object
 * @param sigma current standard deviation of the per-observation
 *        noise terms.
 * @param fvals pointer to memory in which to store the updated beta
 * @param rvals pointer to memory in which to store the updated b (may
 *              be (double*)NULL)
 */
static void MCMC_beta_u(SEXP x, double sigma, double *fvals, double *rvals)
{
    int *dims = DIMS_SLOT(x);
    int i1 = 1, p = dims[p_POS], q = dims[q_POS];
    double *V = V_SLOT(x), *fixef = FIXEF_SLOT(x), *muEta = MUETA_SLOT(x),
	*u = U_SLOT(x), mone[] = {-1,0}, one[] = {1,0};
    CHM_FR L = L_SLOT(x);
    double *del1 = Alloca(q, double), *del2 = Alloca(p, double);
    CHM_DN sol, rhs = N_AS_CHM_DN(del1, q, 1);
    R_CheckStack();

    if (V || muEta) {
	error(_("Update not yet written"));
    } else {			/* Linear mixed model */
	update_L(x);
	update_RX(x);
	lmm_update_fixef_u(x);
				/* Update beta */
	for (int j = 0; j < p; j++) del2[j] = sigma * norm_rand();
	F77_CALL(dtrsv)("U", "N", "N", &p, RX_SLOT(x), &p, del2, &i1);
	for (int j = 0; j < p; j++) fixef[j] += del2[j];
				/* Update u */
	for (int j = 0; j < q; j++) del1[j] = sigma * norm_rand();
	F77_CALL(dgemv)("N", &q, &p, mone, RZX_SLOT(x), &q,
			del2, &i1, one, del1, &i1);
	sol = M_cholmod_solve(CHOLMOD_Lt, L, rhs, &c);
	for (int j = 0; j < q; j++) u[j] += ((double*)(sol->x))[j];
	M_cholmod_free_dense(&sol, &c);
	update_mu(x);	     /* and parts of the deviance slot */
    }
    Memcpy(fvals, fixef, p);
    if (rvals) {
	update_ranef(x);
	Memcpy(rvals, RANEF_SLOT(x), q);
    }
}

/**
 * Update the theta_T parameters from the ST arrays in place.
 *
 * @param x an mer object
 * @param sigma current standard deviation of the per-observation
 *        noise terms.
 */
/* FIXME: Probably should fold this function into MCMC_S */
static void MCMC_T(SEXP x, double sigma)
{
    int *Gp = Gp_SLOT(x), nf = (DIMS_SLOT(x))[nf_POS];
    double **st = Alloca(nf, double*);
    int *nc = Alloca(nf, int), *nlev = Alloca(nf, int);
    R_CheckStack();

    if (ST_nc_nlev(GET_SLOT(x, lme4_STSym), Gp, st, nc, nlev) < 2) return;
    error("Code for non-trivial theta_T not yet written");
}

/**
 * Update the theta_S parameters from the ST arrays in place.
 *
 * @param x an mer object
 * @param sigma current standard deviation of the per-observation
 *        noise terms.
 */
static void MCMC_S(SEXP x, double sigma)
{
    CHM_SP A = A_SLOT(x), Zt = Zt_SLOT(x);
    int *Gp = Gp_SLOT(x), *ai = (int*)(A->i),
	*ap = (int*)(A->p), *dims = DIMS_SLOT(x), *perm = PERM_VEC(x);
    int annz = ap[A->ncol], info, i1 = 1, n = dims[n_POS],
	nf = dims[nf_POS], ns, p = dims[p_POS], pos,
	q = dims[q_POS], znnz = ((int*)(Zt->p))[Zt->ncol];
    double *R, *ax = (double*)(A->x), *b = RANEF_SLOT(x),
	*eta = ETA_SLOT(x), *offset = OFFSET_SLOT(x),  
	*rr, *ss, one = 1, *u = U_SLOT(x), *y = Y_SLOT(x);
    int *nc = Alloca(nf, int), *nlev = Alloca(nf, int),
	*spt = Alloca(nf + 1, int);
    double **st = Alloca(nf, double*);
    R_CheckStack();

    ST_nc_nlev(GET_SLOT(x, lme4_STSym), Gp, st, nc, nlev);
    ns = 0;			/* ns is length(theta_S) */
    spt[0] = 0;			/* pointers into ss for terms */
    for (int i = 0; i < nf; i++) {
	ns += nc[i];
	spt[i + 1] = spt[i] + nc[i];
    }

    if (annz == znnz) { /* Copy Z' to A unless A has new nonzeros */
	Memcpy(ax, (double*)(Zt->x), znnz);
    } else error("Code not yet written for MCMC_S with NLMMs");
				/* Create T'Zt in A */
    Tt_Zt(A, Gp, nc, nlev, st, nf); 
				/* Create P'u in ranef slot */
    for (int i = 0; i < q; i++) b[perm[i]] = u[i];
				/* Create X\beta + offset in eta slot */
    for (int i = 0; i < n; i++) eta[i] = offset ? offset[i] : 0;
    F77_CALL(dgemv)("N", &n, &p, &one, X_SLOT(x), &n,
		    FIXEF_SLOT(x), &i1, &one, eta, &i1);
				/* Allocate R, rr and ss */
    R = Alloca(ns * ns, double); /* crossproduct matrix then factor */
    rr = Alloca(ns, double);	 /* row of model matrix for theta_S */
    ss = Alloca(ns, double);	 /* right hand side, then theta_S */
    R_CheckStack();
    AZERO(R, ns * ns);
    AZERO(ss, ns);
    /* Accumulate crossproduct from pseudo-data part of model matrix */
    for (int i = 0; i < q; i++) {
	int sj = theta_S_ind(i, nf, Gp, nlev, spt);
	AZERO(rr, ns);
	rr[sj] = b[i];
	F77_CALL(dsyr)("U", &ns, &one, rr, &i1, R, &ns);
    }	
    /* Accumulate crossproduct and residual product of the model matrix. */
    /* This is done one row at a time.  Rows of the model matrix
     * correspond to columns of T'Zt */
    for (int j = 0; j < n; j++) { /* jth column of T'Zt */
	AZERO(rr, ns);
	for (int p = ap[j]; p < ap[j + 1]; p++) {
	    int i = ai[p];	/* row in T'Zt */
	    int sj = theta_S_ind(i, nf, Gp, nlev, spt);

	    rr[sj] += ax[p] * b[i];
	    ss[sj] += rr[sj] * (y[j] - eta[j]);
	}
	F77_CALL(dsyr)("U", &ns, &one, rr, &i1, R, &ns);
    }
    F77_CALL(dposv)("U", &ns, &i1, R, &ns, ss, &ns, &info);
    if (info)
	error(_("Model matrix for theta_S is not positive definite, %d."), info);
    for (int j = 0; j < ns; j++) rr[j] = sigma * norm_rand();
    /* Sample from the conditional Gaussian distribution */
    F77_CALL(dtrsv)("U", "N", "N", &ns, R, &ns, rr, &i1);
    for (int j = 0; j < ns; j++) ss[j] += rr[j];
    /* Copy positive part of solution onto diagonals of ST */
    pos = 0;
    for (int i = 0; i < nf; i++) {
	for (int j = 0; j < nc[i]; j++) {
	    st[i][j * (nc[i] + 1)] = (ss[pos] > 0) ? ss[pos] : 0;
	    pos++;
	}
    }
    update_A(x);
}

/* Externally callable functions */

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

/**
 * Generate a Markov-chain Monte Carlo sample from an mer object
 *
 * @param x pointer to an merMCMC object
 * @param fm pointer to an mer object
 *
 * @return x with samples filled in
 */
SEXP mer_MCMCsamp(SEXP x, SEXP fm)
{
    SEXP devsamp = GET_SLOT(x, lme4_devianceSym);
    int *dims = DIMS_SLOT(x), nsamp = LENGTH(devsamp);
    int n = dims[n_POS], np = dims[np_POS],
	p = dims[p_POS], q = dims[q_POS];
    double
	*STsamp = REAL(GET_SLOT(x, lme4_STSym)),
	*d = DEV_SLOT(fm), *dev = REAL(devsamp),
	*sig = SLOT_REAL_NULL(x, lme4_sigmaSym),
	*fixsamp = FIXEF_SLOT(x), *resamp = RANEF_SLOT(x);

    GetRNGstate();
    /* The first column of storage slots contains the fitted values */
    for (int i = 1; i < nsamp; i++) {
/* FIXME: This is probably wrong for a model with weights. */
	if (sig) 		/* update and store sigma */
	    sig[i] = sqrt(d[pwrss_POS]/rchisq((double)(n + q)));
			/* update L, RX, beta, u, eta, mu, res and d */
	MCMC_beta_u(fm, sig ? sig[i] : 1, fixsamp + i * p,
		    resamp + (resamp ? i : 0) * q); 
	dev[i] = d[ML_POS];
				/* update theta_T, theta_S and A */
 	MCMC_T(fm, sig ? sig[i] : 1);
 	MCMC_S(fm, sig ? sig[i] : 1);
	ST_getPars(fm, STsamp + i * np); /* record theta */
    }
    PutRNGstate();
		/* Restore pars from the first columns of the samples */
    ST_setPars(fm, STsamp);
    Memcpy(FIXEF_SLOT(fm), fixsamp, p);
    update_u(fm, 0);
    update_L(fm);
    update_RX(fm);
    lmm_update_fixef_u(fm);
    update_ranef(fm);

    return x;
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
    int *dims = DIMS_SLOT(x), verb = asInteger(verbp);
    int lmm = !(MUETA_SLOT(x) || V_SLOT(x)), nf = dims[nf_POS];
    int nv = dims[np_POS] + (lmm ? 0 : dims[p_POS]);
    int liv = S_iv_length(OPT, nv), lv = S_v_length(OPT, nv);
    double *g = (double*)NULL, *h = (double*)NULL, fx = R_PosInf;
    double *fixef = FIXEF_SLOT(x);
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
    for (int i = 0; i < nv; i++) {
	b[2*i] = R_NegInf; b[2*i+1] = R_PosInf; d[i] = 1;
    }
				/* reset lower bounds on elements of theta_S */
    for (int i = 0, pos = 0; i < nf; i++) {
	int nc = *INTEGER(getAttrib(VECTOR_ELT(ST, i), R_DimSymbol));
	for (int j = 0; j < nc; j++) b[pos + 2*j] = 0;
	pos += nc * (nc + 1);
    }
    do {
	ST_setPars(x, xv);		/* update ST and A */
/* FIXME: Change this so that update_dev is always called and that
 * function does the selection of methods based on lmm. The Memcpy
 * call should always be used but the number of values to copy can be
 * zero. */
	if (lmm) {
	    update_L(x);
	    update_RX(x);
	    fx = lmm_update_fixef_u(x);
	} else {
	    Memcpy(fixef, xv + dims[np_POS], dims[p_POS]);
	    update_u(x, verb);
	    mer_update_dev(x);
	    fx = DEV_SLOT(x)[ML_POS];
	}
	S_nlminb_iterate(b, d, fx, g, h, iv, liv, lv, nv, v, xv);
    } while (iv[0] == 1 || iv[0] == 2);
    update_RX(x);
    lmm_update_fixef_u(x);
    dims[cvg_POS] = iv[0];
    return R_NilValue;
}

/**
 * Extract the conditional variances of the random effects in an mer
 * object.  Some people called these posterior variances, hence the name.
 *
 * @param x pointer to an mer object
 *
 * @return pointer to a list of arrays
 */
SEXP mer_postVar(SEXP x)
{
    int *Gp = Gp_SLOT(x), *dims = DIMS_SLOT(x);
    int nf = dims[nf_POS], q = dims[q_POS];
    SEXP ans = PROTECT(allocVector(VECSXP, nf));
    double *vv, one[] = {1,0}, sc;
    CHM_SP sm1, sm2;
    CHM_DN dm1;
    CHM_FR L = L_SLOT(x);
    int *Perm = (int*)(L->Perm), *iperm = Alloca(q, int),
	*nc = Alloca(nf, int), *nlev = Alloca(nf, int);
    double **st = Alloca(nf, double*);
    R_CheckStack();

    ST_nc_nlev(GET_SLOT(x, lme4_STSym), Gp, st, nc, nlev);
    for (int j = 0; j < q; j++) iperm[Perm[j]] = j; /* inverse permutation */
    sc = dims[useSc_POS] ?
	(DEV_SLOT(x)[dims[isREML_POS] ? sigmaREML_POS : sigmaML_POS]) : 1;
    for (int i = 0; i < nf; i++) {
	int ncisqr = nc[i] * nc[i];
	CHM_SP rhs = M_cholmod_allocate_sparse(q, nc[i], nc[i],
					       1/*sorted*/, 1/*packed*/,
					       0/*stype*/, CHOLMOD_REAL, &c);

	SET_VECTOR_ELT(ans, i, alloc3DArray(REALSXP, nc[i], nc[i], nlev[i]));
	vv = REAL(VECTOR_ELT(ans, i));
	for (int j = 0; j <= nc[i]; j++) ((int *)(rhs->p))[j] = j;
	for (int j = 0; j < nc[i]; j++)
	    ((double *)(rhs->x))[j] = st[i][j * (nc[i] + 1)] * sc;
	for (int k = 0; k < nlev[i]; k++) {
	    double *vvk = vv + k * ncisqr;
	    for (int j = 0; j < nc[i]; j++)
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
    int ncmax, nf = DIMS_SLOT(x)[nf_POS];
    int *nc = Alloca(nf, int), *nlev = Alloca(nf, int);
    double **st = Alloca(nf, double*);
    R_CheckStack();

    ncmax = ST_nc_nlev(ans, Gp_SLOT(x), st, nc, nlev);
    for (int k = 0; k < nf; k++) {
	if (nc[k] > 1) {	/* nothing to do for nc[k] == 1 */
	    int nck = nc[k], nckp1 = nc[k] + 1;
	    double *stk = st[k];

	    for (int j = 0; j < nck; j++) {
		double dd = stk[j * nckp1]; /* diagonal el */
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
	*zi = INTEGER(GET_SLOT(Zt, lme4_iSym)), nf = LENGTH(ST);
    int *nc = Alloca(nf, int), *nlev = Alloca(nf, int),
	nnz = INTEGER(GET_SLOT(Zt, lme4_pSym))[Zdims[1]];
    double *rowsqr = Alloca(Zdims[0], double),
	**st = Alloca(nf, double*),
	*zx = REAL(GET_SLOT(Zt, lme4_xSym));
    R_CheckStack();
    
    ST_nc_nlev(ST, Gp, st, nc, nlev);
    AZERO(rowsqr, Zdims[0]);
    for (int i = 0; i < nnz; i++) rowsqr[zi[i]] += zx[i] * zx[i];
    for (int i = 0; i < nf; i++) {
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
    int npar = DIMS_SLOT(x)[np_POS];

    if (!isReal(pars) || LENGTH(pars) != npar)
	error(_("pars must be a real vector of length %d"), npar);
    ST_setPars(x, REAL(pars));
    return R_NilValue;
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
    double *d = DEV_SLOT(x);
    int *dims = DIMS_SLOT(x);
    const int method = *METHOD_SLOT(x);
    int n = dims[n_POS];
    double dn = (double)n;

    d[disc_POS] = MUETA_SLOT(x) ?
	lme4_devResid(MU_SLOT(x), PWT_SLOT(x), Y_SLOT(x),
		      dims[n_POS], dims[vTyp_POS]) :
	d[wrss_POS];
    /* evaluate maximum likelihood deviance using AGQ or Laplacian method. */
    if(method)
      {

      }
    else
      d[ML_POS] = dn * ( 1 + log(2*PI / dn) ) + dn * log(d[disc_POS] + d[usqr_POS]) + 2 * d[ldL2_POS];
    return R_NilValue;
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

SEXP merMCMC_VarCorr(SEXP x, SEXP typP)
{
    SEXP ST = GET_SLOT(x, lme4_STSym),
	ncs = GET_SLOT(x, install("nc"));
    int *Gp = Gp_SLOT(x), *Sd = INTEGER(GET_DIM(ST)), *nc = INTEGER(ncs);
    int maxnc = 0, nf = LENGTH(ncs), np = Sd[0], nsamp = Sd[1], pos;
    double *sig = SLOT_REAL_NULL(x, lme4_sigmaSym);
    SEXP ans = PROTECT(allocMatrix(REALSXP, nsamp, np + (sig ? 1 : 0)));
    double *av = REAL(ans), *STx = REAL(ST);
    double *as = av + nsamp * np, *t1, *t2, var;
    int *nlev = Alloca(nf, int);
    R_CheckStack();

    for (int j = 0; j < nf; j++) {
	nlev[j] = (Gp[j + 1] - Gp[j])/nc[j];
	if (maxnc < nc[j]) maxnc = nc[j];
    }
    if (maxnc > 1) {
	t1 = Alloca(maxnc * maxnc, double);
	t2 = Alloca(maxnc * maxnc, double);
	R_CheckStack();
    }
    
    for (int i = 0; i < nsamp; i++) {
	var = 1; pos = 0;
	if (sig) var = as[i] = sig[i] * sig[i];
	for (int k = 0; k < nf; k++) {
	    if (nc[k] < 2) {
		double sd = STx[pos + i * np] * sig[i];
		av[i + nsamp * pos++] = sd * sd;
	    }
	    else error(_("Code not yet written"));
	}
    }
	
    UNPROTECT(1);
    return ans;
}
/**
 * Check validity of an mer object
 *
 * @param x Pointer to an mer object
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
	flistP = GET_SLOT(x, lme4_flistSym), asgnP;
    int *Gp = INTEGER(GpP), *dd = INTEGER(dimsP), *asgn;
    int n = dd[n_POS], nf = dd[nf_POS], nq, nfl,
	p = dd[p_POS], q = dd[q_POS], s = dd[s_POS];
    int nv = n * s;
    CHM_SP Zt = Zt_SLOT(x), A =  A_SLOT(x);
    CHM_FR L = L_SLOT(x);
    char *buf = Alloca(BUF_SIZE + 1, char);
    R_CheckStack();
				/* check lengths */
    if (nf < 1 || LENGTH(ST) != nf)
	return mkString(_("Slot ST must have length dims['nf']"));
    asgnP = getAttrib(flistP, install("assign"));
    if (!isInteger(asgnP) || LENGTH(asgnP) != nf)
	return mkString(_("Slot flist must have integer attribute 'assign' of length dims['nf']"));
    asgn = INTEGER(asgnP);
    nfl = LENGTH(flistP);
    for (int i = 0; i < nf; i++)
	if (asgn[i] <= 0 || asgn[i] > nfl)
	    return mkString(_("All elements of the assign attribute must be in [1,length(ST)]"));
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
    if (chkLen(buf, BUF_SIZE, x, lme4_pWtSym, n, 1)) return(mkString(buf));
    if (chkLen(buf, BUF_SIZE, x, lme4_ranefSym, q, 0)) return(mkString(buf));
    if (chkLen(buf, BUF_SIZE, x, lme4_residSym, n, 0)) return(mkString(buf));
    if (chkLen(buf, BUF_SIZE, x, lme4_sqrtrWtSym, n, 1)) return(mkString(buf));
    if (chkLen(buf, BUF_SIZE, x, lme4_uSym, q, 0)) return(mkString(buf));
    if (chkLen(buf, BUF_SIZE, x, lme4_VSym, nv, 1)) return(mkString(buf));
    if (chkLen(buf, BUF_SIZE, x, lme4_varSym, n, 1)) return(mkString(buf));
    if (chkLen(buf, BUF_SIZE, x, lme4_ySym, n, 0)) return(mkString(buf));
    if (chkDims(buf, BUF_SIZE, x, lme4_XSym, nv, p)) return(mkString(buf));
    if (chkDims(buf, BUF_SIZE, x, lme4_RZXSym, q, p)) return(mkString(buf));
    if (chkDims(buf, BUF_SIZE, x, lme4_RXSym, p, p)) return(mkString(buf));

    nq = 0;
    for (int i = 0; i < LENGTH(flistP); i++) {
	SEXP fli = VECTOR_ELT(flistP, i);
	if (!isFactor(fli))
	    return mkString(_("flist must be a list of factors"));
/* 	nq += dm[0] * LENGTH(getAttrib(fli, R_LevelsSymbol)); */
    }
    for (int i = 0; i < nf; i++) {
	SEXP STi = VECTOR_ELT(ST, i);
	int *dm = INTEGER(getAttrib(STi, R_DimSymbol));
	if (!isMatrix(STi) || !isReal(STi) || dm[0] != dm[1])
	    return
		mkString(_("Slot ST must be a list of square numeric matrices"));
	if (Gp[i] > Gp[i + 1])
	    return mkString(_("Gp must be non-decreasing"));
    }
#if 0
/* FIXME: Need to incorporate the assign attribute in the calculation of nq */
    if (q != nq)
	return mkString(_("q is not sum of columns by levels"));
#endif
    return ScalarLogical(1);
}

/**
 * Check validity of an merMCMC object
 *
 * @param x Pointer to an merMCMC object
 *
 * @return TRUE if the object is a valid merMCMC object, otherwise a string
 *         that describes the violation.
 */
SEXP merMCMC_validate(SEXP x)
{
    SEXP GpP = GET_SLOT(x, lme4_GpSym),
	devianceP = GET_SLOT(x, lme4_devianceSym),
	dimsP = GET_SLOT(x, lme4_dimsSym);
    int *Gp = INTEGER(GpP), *dd = INTEGER(dimsP);
    int nf = dd[nf_POS], np = dd[np_POS], nsamp = LENGTH(devianceP),
	p = dd[p_POS], q = dd[q_POS];
    char *buf = Alloca(BUF_SIZE + 1, char);
    R_CheckStack();
				/* check lengths */
    if (nsamp <= 0)
	return mkString(_("number of samples must be positive"));
    if (LENGTH(dimsP) != (cvg_POS + 1) ||
	LENGTH(getAttrib(dimsP, R_NamesSymbol)) != (cvg_POS + 1))
	return mkString(_("dims slot not named or incorrect length"));
    if (LENGTH(GpP) != nf + 1)
	return mkString(_("Slot Gp must have length dims['nf'] + 1"));
    if (Gp[0] != 0 || Gp[nf] != q)
	return mkString(_("Gp[1] != 0 or Gp[dims['nf'] + 1] != dims['q']"));

    if (chkLen(buf, BUF_SIZE, x, lme4_ncSym, nf, 0))
	return(mkString(buf));
    if (chkLen(buf, BUF_SIZE, x, lme4_sigmaSym, nsamp, !dd[useSc_POS]))
	return(mkString(buf));
    if (chkDims(buf, BUF_SIZE, x, lme4_STSym, np, nsamp))
	return(mkString(buf));
    if (chkDims(buf, BUF_SIZE, x, lme4_fixefSym, p, nsamp))
	return(mkString(buf));
    if (LENGTH(GET_SLOT(x, lme4_ranefSym)))
	if (chkDims(buf, BUF_SIZE, x, lme4_ranefSym, q, nsamp))
	    return(mkString(buf));
    return ScalarLogical(1);
}

/**
 * Update the eta, mu, resid and var slots in a sparseRasch object
 * from the current values of the model parameters in the beta slot.
 *
 * @param x pointer to an sparseRasch object
 *
 * @return R_NilValue
 */
SEXP spR_update_mu(SEXP x)
{
    int *dims = DIMS_SLOT(x);
    int n = dims[n_POS];
    double *d = DEV_SLOT(x), *eta = Calloc(n, double), *mu = MU_SLOT(x),
	*offset = OFFSET_SLOT(x), *srwt = SRWT_SLOT(x),
	*res = RESID_SLOT(x), *y = Y_SLOT(x), one[] = {1,0};
    CHM_SP Zt = Zt_SLOT(x);
    CHM_DN cbeta = AS_CHM_DN(GET_SLOT(x, lme4_fixefSym)),
	ceta = N_AS_CHM_DN(eta, n, 1);
    R_CheckStack();
    
    for (int i = 0; i < n; i++) eta[i] = offset ? offset[i] : 0;
    if (!M_cholmod_sdmult(Zt, 1 /* trans */, one, one, cbeta, ceta, &c))
	error(_("cholmod_sdmult error returned"));
    lme4_muEta(mu, MUETA_SLOT(x), eta, n, dims[lTyp_POS]);
    lme4_varFunc(VAR_SLOT(x), mu, n, dims[vTyp_POS]);

    d[wrss_POS] = 0;		/* update resid slot and d[wrss_POS] */
    for (int i = 0; i < n; i++) {
	res[i] = (y[i] - mu[i]) * srwt[i];
	d[wrss_POS] += res[i] * res[i];
    }

    Free(eta);
    return R_NilValue;
}


/**
 * Optimize the log-likelihood for the sparse representation of a
 * Rasch model.
 *
 * @param x pointer to a sparseRasch object
 * @param verbP pointer to indicator of verbose output
 *
 * @return R_NilValue
 */

SEXP spR_optimize(SEXP x, SEXP verbP)
{
    int *zp, m, n, nnz, verb = asInteger(verbP);
    double *Zcp, *Ztx, *d = DEV_SLOT(x), *fixef = FIXEF_SLOT(x),
	*mu = MU_SLOT(x), *muEta = MUETA_SLOT(x),
	*pWt = PWT_SLOT(x), *res = RESID_SLOT(x),
	*srwt = SRWT_SLOT(x), *tmp, *tmp2,
	*var = VAR_SLOT(x), *y = Y_SLOT(x), cfac, crit,
	one[] = {1, 0}, step, wrss_old, zero[] = {0, 0};
    CHM_SP Zt = Zt_SLOT(x);
    CHM_FR L = L_SLOT(x);
    CHM_DN cres = N_AS_CHM_DN(res, Zt->ncol, 1), ctmp, sol;
    R_CheckStack();
    
    zp = (int*)Zt->p;
    m = Zt->nrow;
    n = Zt->ncol;
    nnz = zp[n];
    Zcp = Calloc(nnz, double);
    Ztx = (double*)Zt->x;
    tmp = Calloc(m, double);
    tmp2 = Calloc(m, double);
    ctmp = N_AS_CHM_DN(tmp, m, 1);
    cfac = ((double)n) / ((double)m);
    R_CheckStack();

    spR_update_mu(x);
    for (int i = 0; ; i++) {
	d[wrss_POS] = 0;
	for (int j = 0; j < n; j++) {
	    srwt[j] = sqrt((pWt ? pWt[j] : 1) / var[j]);
	    res[j] = srwt[j] * (y[j] - mu[j]);
	    d[wrss_POS] += res[j] * res[j];
	    for (int p = zp[j]; p < zp[j + 1]; p++)
		Zcp[p] = srwt[j] * muEta[j] * Ztx[p];
	}
	Zt->x = (void*)Zcp;
	if (!M_cholmod_factorize_p(Zt, zero, (int*)NULL, 0 /*fsize*/, L, &c))
	    error(_("cholmod_factorize_p failed: status %d, minor %d from ncol %d"),
		  c.status, L->minor, L->n);
	wrss_old = d[wrss_POS];
				/* tmp := Zt %*% muEta %*% var %*% wtdResid */
	M_cholmod_sdmult(Zt, 0 /* notrans */, one, zero, cres, ctmp, &c);
	Memcpy(tmp2, tmp, m);
	apply_perm(tmp, tmp2, (int*)L->Perm, m);
				/* solve L %*% sol = tmp */
	if (!(sol = M_cholmod_solve(CHOLMOD_L, L, ctmp, &c)))
	    error(_("cholmod_solve (CHOLMOD_L) failed"));
	Memcpy(tmp, (double*)(sol->x), m);
	M_cholmod_free_dense(&sol, &c);
				/* evaluate convergence criterion */
	crit = cfac * sqr_length(tmp, m) / wrss_old;
	if (crit < CM_TOL) break; /* don't do needless evaluations */
				/* solve t(L) %*% sol = tmp */
	if (!(sol = M_cholmod_solve(CHOLMOD_Lt, L, ctmp, &c)))
	    error(_("cholmod_solve (CHOLMOD_Lt) failed"));
	Memcpy(tmp, (double*)(sol->x), m);
	M_cholmod_free_dense(&sol, &c);

	Memcpy(tmp2, fixef, m);
	for (step = 1; step > CM_SMIN; step /= 2) { /* step halving */
	    for (int j = 0; j < m; j++) fixef[j] = tmp2[j] + step * tmp[j];
	    spR_update_mu(x);
	    if (verb < 0)
		Rprintf("%2d,%8.6f,%12.4g: %15.6g %15.6g %15.6g %15.6g\n",
			i, step, crit, d[wrss_POS], wrss_old, fixef[1], fixef[2]);
	    if (d[wrss_POS] < wrss_old) {
		wrss_old = d[wrss_POS];
		break;
	    }
	}
	if (step <= CM_SMIN || i > CM_MAXITER) return 0;
    }
    Free(tmp2);
    Free(Zcp);
    Free(tmp);
    return R_NilValue;
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
	for (int i = 0; i < n; i++) x[i * np1] += incr;
    }
}

/**
 * Update the ST list of arrays and the sparse model matrix A
 *
 * @param x an mer object
 * @param sigma current standard deviation of the per-observation
 *        noise terms.
 * @param vals pointer to memory in which to store the updated values
 *        of the ST parameters
 */
static void MCMC_ST(SEXP x, double sigma, double *vals)
{
    int *Gp = Gp_SLOT(x), *dims = DIMS_SLOT(x), *perm = PERM_VEC(x);
    int nf = dims[nf_POS], q = dims[q_POS], pos = 0;
    double *Ptu = Calloc(q, double), *u = U_SLOT(x);
    double **st = Alloca(nf, double*);
    int *nc = Alloca(nf, int), *nlev = Alloca(nf, int);
    R_CheckStack();

    ST_nc_nlev(GET_SLOT(x, lme4_STSym), Gp, st, nc, nlev);
				/* inverse permutation of u */
    for (int j = 0; j < q; j++) Ptu[perm[j]] = u[j];

    for (int i = 0; i < nf; i++) {
	int qi = nc[i], nl = nlev[i];
	double *sti = st[i], *ui = Ptu + Gp[i];

	if (qi == 1) {
	    *sti *= sqrt(sqr_length(ui, nl)/rchisq((double)nl))/sigma;
	    vals[pos++] = *sti;
	} else {
	    int info, qisq = qi * qi, qip1 = qi + 1;
	    double invsigsq = 1/(sigma * sigma), one[] = {1,0}, zer[] = {0,0};

	    double *scal = Alloca(qisq, double), *tmp = Alloca(qisq, double);
	    R_CheckStack();
				/* scal := crossprod(U_i)/sigma^2 */
	    F77_CALL(dsyrk)("L", "T", &qi, &nl, &invsigsq, ui, &nl,
			    zer, scal, &qi); 
	    F77_CALL(dpotrf)("L", &qi, scal, &qi, &info);
	    if (info)
		error(_("scal matrix %d is not positive definite"), i + 1);
	    /* solve W_i' B_i = R_i for W_i a random std Wishart factor */
	    F77_CALL(dtrsm)("L", "L", "T", "N", &qi, &qi, one,
			    std_rWishart_factor((double)(nl - qi + 1), qi, 0, tmp),
			    &qi, scal, &qi);
				/* B_i := T(theta) %*% S(theta) %*% B_i */
	    for (int k = 0; k < qi; k++)
		for (int j = 0; j < qi; j++)
		    scal[j * qi + k] *= st[i][k * qip1];
	    F77_CALL(dtrmm)("L", "L", "N", "U", &qi, &qi, one, st[i], &qi, scal, &qi);
				/* Create the lower Cholesky factor */
	    F77_CALL(dsyrk)("L", "T", &qi, &qi, one, scal, &qi, zer, tmp, &qi); 
	    F77_CALL(dpotrf)("L", &qi, tmp, &qi, &info);
	    if (info)
		error(_("crossproduct matrix %d is not positive definite"), i + 1);
				/* Update S_i and T_i */
	    for (int j = 0; j < qi; j++) vals[pos++] = sti[j * qip1] = tmp[j * qip1];
	    for (int j = 0; j < qi; j++)
		for (int k = j + 1; k < qi; k++)
		    vals[pos++] = sti[j * qi + k] = tmp[j * qi + k]/sti[j * qip1];
	}
    }
    update_A(x);
    Free(Ptu);
}

/* Functions for sampling from a Wishart distribution */
/* FIXME: Move these to the R sources */

/**
 * Simulate the Cholesky factor of a standardized Wishart variate with
 * dimension p and nu degrees of freedom.
 *
 * @param nu degrees of freedom
 * @param p dimension of the Wishart distribution
 * @param upper if 0 the result is lower triangular, otherwise upper
                triangular
 * @param ans array of size p * p to hold the result
 *
 * @return ans
 */
static double
*std_rWishart_factor(double nu, int p, int upper, double ans[])
{
    int pp1 = p + 1;

    if (nu < (double) p || p <= 0)
	error("inconsistent degrees of freedom and dimension");

    AZERO(ans, p * p);
    for (int j = 0; j < p; j++) {	/* jth column */
	ans[j * pp1] = sqrt(rchisq(nu - (double) j));
	for (int i = 0; i < j; i++) {
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
 * @param nuP Degrees of freedom
 * @param scal Positive-definite scale matrix
 *
 * @return
 */
SEXP
lme4_rWishart(SEXP ns, SEXP nuP, SEXP scal)
{
    SEXP ans;
    int *dims = INTEGER(getAttrib(scal, R_DimSymbol)), info,
	n = asInteger(ns), psqr;
    double *scCp, *ansp, *tmp, nu = asReal(nuP), one = 1, zero = 0;

    if (!isMatrix(scal) || !isReal(scal) || dims[0] != dims[1])
	error("scal must be a square, real matrix");
    if (n <= 0) n = 1;
    psqr = dims[0] * dims[0];
    tmp = Alloca(psqr, double);
    scCp = Alloca(psqr, double);
    R_CheckStack();

    Memcpy(scCp, REAL(scal), psqr);
    AZERO(tmp, psqr);
    F77_CALL(dpotrf)("U", &(dims[0]), scCp, &(dims[0]), &info);
    if (info)
	error("scal matrix is not positive-definite");
    PROTECT(ans = alloc3DArray(REALSXP, dims[0], dims[0], n));
    ansp = REAL(ans);
    GetRNGstate();
    for (int j = 0; j < n; j++) {
	double *ansj = ansp + j * psqr;
	std_rWishart_factor(nu, dims[0], 1, tmp);
	F77_CALL(dtrmm)("R", "U", "N", "N", dims, dims,
			&one, scCp, dims, tmp, dims);
	F77_CALL(dsyrk)("U", "T", &(dims[1]), &(dims[1]),
			&one, tmp, &(dims[1]),
			&zero, ansj, &(dims[1]));

	for (int i = 1; i < dims[0]; i++)
	    for (int k = 0; k < i; k++)
		ansj[i + k * dims[0]] = ansj[k + i * dims[0]];
    }

    PutRNGstate();
    UNPROTECT(1);
    return ans;
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

#endif
