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
#include <Rmath.h>	   // for dnorm5, etc. 
#include <R_ext/Lapack.h>  // for Lapack (dpotrf, etc.) and BLAS 
#include <R_ext/Boolean.h> // for R version of TRUE and FALSE

#ifdef ENABLE_NLS	/** Allow for translation of error messages */
#include <libintl.h>
#define _(String) dgettext ("lme4", String)
#else
#define _(String) (String)
#endif

// Inlined utilties


/**
 * Copy the first nn elements of src to dest
 *
 * @param src source vector
 * @param dest destination vector
 * @param nn number of elements in src and dest
 *
 * @return dest
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

SEXP lme4_ghq(SEXP np);

#ifdef	__cplusplus
}
#endif 
#endif /* LME4_LME4UTILS_H */
