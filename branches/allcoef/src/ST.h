#ifndef LME4_ST_H
#define LME4_ST_H

#include <R.h>
/* Rdefines.h includes Rinternals.h (for SEXP, REAL, etc.) and defines
 * GET_SLOT, MAKE_CLASS, NEW_OBJECT, SET_SLOT, etc. */
#include <Rdefines.h>
#include <R_ext/Lapack.h>        // for Lapack (dpotrf, etc.) and BLAS
#include "Matrix.h"

#ifdef	__cplusplus
extern "C" {
#endif

    SEXP ST_Lambda(SEXP x);
    SEXP ST_Tmatrix(SEXP x);
    SEXP ST_bounds(SEXP x);
    SEXP ST_chol(SEXP x);
    SEXP ST_create_A(SEXP x, SEXP rho);
    SEXP ST_create_ranef(SEXP ST, SEXP u, SEXP perm);
    SEXP ST_getPars(SEXP x);
    SEXP ST_initialize(SEXP x, SEXP rho);
    SEXP ST_setPars(SEXP x, SEXP pars, SEXP rho);
    SEXP ST_update_A(SEXP ST, SEXP rho);
    SEXP ST_validate(SEXP x);

#ifdef	__cplusplus
}

#include "lme4utils.hpp"

class STinternal {
public:
    STinternal(SEXP x)		  //< external, R-level object
    {assign_vals(GET_SLOT(x, lme4_STSym), GET_SLOT(x, lme4_GpSym));}

    STinternal(SEXP ST,		  //< matrix list 
	       SEXP Gp)		  //< group pointers
    {assign_vals(ST, Gp);}

    ~STinternal() {delete[] nlev; delete[] nc; delete[] st;}

/**
 * Initialize the parameters in the ST slot from the Zt matrix
 *
 * @param Zt sparse transposed random effects model matrix
 *
 */
    void initialize(SEXP Zt);

/**
 * Validate method.  This part is a no-op because validation is
 * done in the constructor.
 */
    SEXP validate() {return ScalarLogical(1);}

/**
 * Utility that returns the index of the term corresponding to a row
 * or column index in Lambda.
 *
 * @param ind index in Lambda - must be in the range [0, Gp[nt]]
 */
    int Gp_grp(int ind);

/**
 * Fill numeric vectors lower and upper, of length np with parameter bounds
 *
 * @param lower pointer to an numeric vector of length np
 * @param upper pointer to an numeric vector of length np
 */
    void bounds(double *lower, double *upper);

/**
 * Assign values to the diagonal of S
 *
 * @param d pointer to a vector of Gp[nt] values
 * @return d
 */
    double *Sdiag(double *d);

/**
 * Create the T matrix as a CHM_SP object
 *
 * @return T as a CHM_SP object
 */
    CHM_SP Tmatrix();
    
/**
 * Create the Lambda matrix as a CHM_SP object
 *
 * @return Lambda as a CHM_SP object
 */
    CHM_SP Lambda();
/**
 * Create A from Zt
 *
 * @return A as a CHM_SP object
 */
    CHM_SP create_A(CHM_SP Zt);

/**
 * Update A from Zt
 *
 * @param Zt original model matrix
 * @param A scaled model matrix
 */
    void update_A(CHM_SP Zt, CHM_SP A);

/**
 * Extract the parameter values
 *
 * @param pars vector of length np
 * @return pars
 */
    double *getPars(double *pars);

    int npars() {return np;}

    SEXP create_ranef(SEXP u, SEXP perm);

    void chol(SEXP ans);

/**
 * Install new parameters in the ST slot.
 *
 * @param pars double vector of the appropriate length
 *
 */
    void setPars(const double *pars);

    int theta_S_ind(int i, int *spt)
    {
	int trm = Gp_grp(i);
	return (spt[trm] + (i - Gp[trm]) / nlev[trm]);
    }
    
private:
    double **st;
    int *Gp, *nc, *nlev, nt, maxnc, np;
    void assign_vals(SEXP ST, SEXP Gpp);
};

#endif /* __cplusplus */

#endif /* LME4_ST_H */
