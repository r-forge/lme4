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
    SEXP ST_create_A(SEXP x, SEXP rho);
    SEXP ST_getPars(SEXP x);
    SEXP ST_initialize(SEXP x, SEXP rho);
    SEXP ST_setPars(SEXP x, SEXP pars, SEXP rho);
    SEXP ST_update_A(SEXP ST, SEXP rho);
    SEXP ST_create_ranef(SEXP ST, SEXP rho);
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
 * Fill the integer vector b, of length 2 * np with parameter bounds
 * in the order low_1, high_1, low_2, ...
 *
 * @param b pointer to an integer vector of length 2 * np
 * @return b with values filled in
 */
    double *bounds(double *b);

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

    void update_ranef(const double *u, const int *perm, double *b);
    
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
    void assign_vals(SEXP ST, SEXP Gpp)
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

};

#endif /* __cplusplus */

#endif /* LME4_ST_H */
