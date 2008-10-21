#include "Matrix.h"

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
#elif defined(__sun) || defined(_AIX)
/* this is necessary (and sufficient) for Solaris 10 and AIX 6: */
# include <alloca.h>
#endif

#include <R_ext/Boolean.h> 	// ensure that the R version of TRUE and FALSE are used

/** alloca n elements of type t */
#define Alloca(n, t)   (t *) alloca( (size_t) ( (n) * sizeof(t) ) )

extern
#include "Syms.h"
extern	       /** cholmod_common struct initialized in R_init_lme4 */
cholmod_common c;

/** zero an array */
#define AZERO(x, n) {int _I_, _SZ_ = (n); for(_I_ = 0; _I_ < _SZ_; _I_++) (x)[_I_] = 0;}

/** positions in the deviance vector */
enum devP {
    ML_POS=0,			/**<Maximum likelihood estimation criterion  */
    REML_POS,			/**<REML criterion */
    ldL2_POS,			/**<2*log-determinant of L */
    ldRX2_POS,			/**<2*log-determinant of RX */
    sigmaML_POS,		/**<current ML estimate of sigma */
    sigmaREML_POS,		/**<current REML estimate of sigma */
    pwrss_POS,			/**<penalized weighted residual sum of squares */
    disc_POS,			/**<discrepancy */
    usqr_POS,			/**<squared length of u */
    wrss_POS,			/**<weighted residual sum of squares  */
    dev_POS,			/**<deviance - defined for quasi families  */
    llik_POS,			/**<log-likelihood - undefined for quasi families  */
    NULLdev_POS			/**<null deviance */
};

/** positions in the dims vector */
enum dimP {
    nt_POS=0,			/**<number of terms in random effects */
    n_POS,			/**<number of observations */
    p_POS,			/**<number of fixed-effects parameters */
    q_POS,			/**<number of random effects */
    s_POS,			/**<number of variables in h (1 unless nonlinear) */
    np_POS,			/**<total number of parameters for T and S */
    LMM_POS,			/**<is the model a linear mixed model? */
    isREML_POS,			/**<indicator of REML estimation */
    fTyp_POS,			/**<family type for generalized model */
    lTyp_POS,			/**<link type for generalized model */
    vTyp_POS,			/**<variance type for generalized model */
    nest_POS,			/**<indicator of nested grouping factors */
    useSc_POS,			/**<does the family use a separate scale parameter */
    nAGQ_POS,			/**<number of adaptive Gauss-Hermite quadrature pts */
    verb_POS,			/**<verbose output in mer_optimize? */
    mxit_POS,			/**<maximum # of iterations in mer_optimize */
    mxfn_POS,			/**<maximum # of function evaluations in mer_optimize */
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

/** Return the double pointer to the deviance slot */
#define DEV_SLOT(x) SLOT_REAL_NULL(x, lme4_devianceSym)

/** Return the integer pointer to the dims slot */
#define DIMS_SLOT(x) INTEGER(GET_SLOT(x, lme4_dimsSym))

/** Return the double pointer to the eta slot */
#define ETA_SLOT(x) SLOT_REAL_NULL(x, lme4_etaSym)

/** Return the double pointer to the etaGamma slot */
#define ETAGAMMA_SLOT(x) SLOT_REAL_NULL(x, lme4_etaGammaSym)

/** Return the double pointer to the fixef slot */
#define FIXEF_SLOT(x) SLOT_REAL_NULL(x, lme4_fixefSym)

/** Return the double pointer to the ghw slot */
#define GHW_SLOT(x) SLOT_REAL_NULL(x, lme4_ghwSym)

/** Return the double pointer to the ghx slot */
#define GHX_SLOT(x) SLOT_REAL_NULL(x, lme4_ghxSym)

/** Return the integer pointer to the Gp slot */
#define Gp_SLOT(x) INTEGER(GET_SLOT(x, lme4_GpSym))

/** Allocate (alloca) a cholmod_factor struct, populate it with values
 * from the L slot and return the pointer. */
#define L_SLOT(x) AS_CHM_FR(GET_SLOT(x, lme4_LSym))

/** Return the double pointer to the mu slot */
#define MU_SLOT(x) SLOT_REAL_NULL(x, lme4_muSym)

/** Return the double pointer to the muEta slot or (double*) NULL if
 * muEta has length 0) */
#define MUETA_SLOT(x) SLOT_REAL_NULL(x, lme4_muEtaSym)

/** Return the double pointer to the offset slot or (double*) NULL if
 * offset has length 0) */
#define OFFSET_SLOT(x) SLOT_REAL_NULL(x, lme4_offsetSym)

/** Return the integer pointer to the dims slot */
#define PERM_SLOT(x) INTEGER(GET_SLOT(x, lme4_permSym))

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

