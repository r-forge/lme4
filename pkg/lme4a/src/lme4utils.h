#ifndef LME4_LME4UTILS_H
#define LME4_LME4UTILS_H
#ifdef	__cplusplus
#define NO_C_HEADERS
#include <cstdio>
#include <cmath>
#include <cstring>
#include <cstdarg>
extern "C" {
#endif 

#include <R.h>
// Rdefines.h includes Rinternals.h (for SEXP, REAL, etc.) and defines
// GET_SLOT, MAKE_CLASS, NEW_OBJECT, SET_SLOT, etc.
#include <Rdefines.h>
#include "Matrix.h"
extern
#include "Syms.h"
extern	       /** cholmod_common struct initialized in R_init_lme4 */
cholmod_common c;
#include <Rmath.h>	      /* for dnorm5, etc. */
#include <R_ext/Lapack.h>     /* for Lapack (dpotrf, etc.) and BLAS */

#ifdef ENABLE_NLS	/** Allow for translation of error messages */
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

#include <R_ext/Boolean.h> // ensure that the R version of TRUE and FALSE are used

/** alloca n elements of type t */
#define Alloca(n, t)   (t *) alloca( (size_t) ( (n) * sizeof(t) ) )

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
    LMM_POS=0,			/**<is the model a linear mixed model? */
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

// Inlined utilties


/**
 * Copy the first nn elements of src to dest
 *
 * @param src source vector
 * @param dest destination vector
 * @param nn number of elements in src and dest
 *
 * @return dest
 * @fixme this and int_cpy should use a template
 */
static inline double *dble_cpy(double *dest, const double *src, int nn)
{
    for (int i = 0; i < nn; i++)
	dest[i] = src[i];
    return dest;
}

/**
 * Copy the first nn elements of src to dest
 *
 * @param src source vector
 * @param dest destination vector
 * @param nn number of elements in src and dest
 *
 * @return dest
 */
static inline int *int_cpy(int *dest, const int *src, int nn)
{
    for (int i = 0; i < nn; i++)
	dest[i] = src[i];
    return dest;
}

/**
 * Zero the first nn elements of double pointer dest
 *
 * @param dest vector
 * @param nn number of elements in dest
 *
 * @return dest
 */
static inline double *dble_zero(double *dest, int nn)
{
    for (int i = 0; i < nn; i++)
	dest[i] = 0.;
    return dest;
}

/**
 * Zero the first nn elements of int pointer dest
 *
 * @param dest vector
 * @param nn number of elements in dest
 *
 * @return dest
 */
static inline int *int_zero(int *dest, int nn)
{
    for (int i = 0; i < nn; i++)
	dest[i] = 0;
    return dest;
}

/**
 * Evaluate the squared length of the first nn elements of x
 *
 * @param x vector
 * @param nn number of elements in x
 *
 * @return the squared length of x
 */
static inline double sqr_length(const double *x, int nn) {
    double ans = 0;
    for (int i = 0; i < nn; i++)
	ans += x[i] * x[i];
    return ans;
}

static inline SEXP findVarBound(SEXP rho, SEXP nm) {
    SEXP var = findVarInFrame(rho, nm);
    if (var == R_UnboundValue)
	error(_("object named '%s' not found in environment"),
	      CHAR(PRINTNAME(nm)));
    return var;
}

static inline SEXP
ALLOC_SLOT(SEXP obj, SEXP nm, SEXPTYPE type, int length)
{
    SEXP val = allocVector(type, length);
    SET_SLOT(obj, nm, val);
    return val;
}

static inline int compare_int_vecs(int *dest, int *src, int n) {
    for (int i = 0; i < n; i++)
	if (dest[i] != src[i]) return 0;
    return 1;
}

/**
 * Evaluate y * log(y/mu) with the correct limiting value at y = 0.
 *
 * @param y 
 * @param mu
 *
 * @return y * log(y/mu) for y > 0, 0 for y == 0.
 */
inline double y_log_y(double y, double mu)
{
    return (y) ? (y * log(y/mu)) : 0;
}

SEXP lme4_ghq(SEXP np);

#ifdef	__cplusplus
}
#endif 
#endif /* LME4_LME4UTILS_H */
