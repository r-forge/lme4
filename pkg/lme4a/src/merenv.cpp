#include "merenv.h"

#include "Matrix.h"
#include "lme4utils.hpp"

class merenv {
public:
    merenv(SEXP rho);		//< instantiate from an environment
    ~merenv(){
	delete L;
	delete Zt;
	delete Ut;
	if (!dLam) delete Lambda;
    }
    double update_dev();

private:
    static int i1;
    int *Lind, *perm, REML, dLam, n, nmp, p, q, sparseX;
    double *RX, *RZX, *X, *ZtX, *Zty, *beta, *dlamb, *mu,
	*offset, *pWt, *theta, *res, *u,  *y;
    CHM_FR L;
    CHM_SP Lambda, Ut, Zt, ZtsX;
};

int merenv::i1 = 1;

/**
 * Copy the first nn elements of src to dest
 *
 * @param src source vector
 * @param dest destination vector
 * @param nn number of elements in src and dest
 *
 * @return dest
 */
inline double *dble_cpy(double *dest, const double *src, int nn)
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
inline double *dble_zero(double *dest, int nn)
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
inline int *int_zero(int *dest, int nn)
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
inline double sqr_length(const double *x, int nn)
{
    double ans = 0;
    for (int i = 0; i < nn; i++)
	ans += x[i] * x[i];
    return ans;
}

CHM_SP VAR_CHM_SP(SEXP rho, SEXP nm, int nrow, int ncol)
{
    CHM_SP ans = new cholmod_sparse;
    const char *pn = CHAR(PRINTNAME(nm));
    SEXP var = findVarInFrame(rho, nm);
    if (var == R_UnboundValue)
	error(_("object named '%s' not found in environment"), pn);
    M_as_cholmod_sparse(ans, var, (Rboolean)TRUE, (Rboolean)FALSE);
    if (((int)(ans->nrow)) != nrow)
	error(_("Number of rows of %s is %d, should be %d"),
	      pn, ans->nrow, nrow);
    if (((int)(ans->ncol)) != ncol)
	error(_("Number of columns of %s is %d, should be %d"),
	      pn, ans->ncol, ncol);
    return ans;
}

double *VAR_dMatrix_x(SEXP rho, SEXP nm, int nrow, int ncol)
{    
    const char *pn = CHAR(PRINTNAME(nm));
    SEXP var = findVarInFrame(rho, nm);
    if (var == R_UnboundValue)
	error(_("object named '%s' not found in environment"), pn);

    // FIXME: Should check here to ensure that the object "is" a dMatrix
    int *dims = INTEGER(GET_SLOT(var, lme4_DimSym));
    if (dims[0] != nrow || dims[1] != ncol)
	error(_("object named '%s' should be %d by %d"), pn, nrow, ncol);
    return REAL(GET_SLOT(var, lme4_xSym));
}


// Definition of methods for the mer class

merenv::merenv(SEXP rho)
{
    // Extract slots that must have positive length.
    // Get dimensions of the problem
    SEXP sl = findVarInFrame(rho, lme4_ySym);
    if (sl == R_UnboundValue)
	error(_("Environment must contain response y"));
    if (!(n = LENGTH(sl)) || !isReal(sl)) // n = length of response
	error(_("Response vector y must be numeric (double)"));
    y = REAL(sl);

    sl = findVarInFrame(rho, lme4_betaSym);
    if (sl == R_UnboundValue)
	error(_("Environment must contain numeric vector beta"));
    if (!(p = LENGTH(sl)) || !isReal(sl)) // p = length of fixef
	error(_("beta vector must be numeric (double)"));
    beta = REAL(sl);

    sl = findVarInFrame(rho, lme4_uSym);
    if (sl == R_UnboundValue)
	error(_("Environment must contain numeric vector u"));
    if (!(q = LENGTH(sl)) || !isReal(sl)) // q = length of u
	error(_("u vector must be numeric (double)"));
    u = REAL(sl);
    
    Zt = VAR_CHM_SP(rho, lme4_ZtSym, q, n);
    Ut = VAR_CHM_SP(rho, lme4_UtSym, q, n);
    L = new cholmod_factor;
    M_as_cholmod_factor(L, findVarInFrame(rho, lme4_LSym));

    if (dLam = asLogical(findVarInFrame(rho,
					install("diagonalLambda")))) {
	dlamb = VAR_dMatrix_x(rho, lme4_LambdaSym, q, q);
    } else {
	Lambda = VAR_CHM_SP(rho, lme4_LambdaSym, q, q);
    }
    if (sparseX = asLogical(findVarInFrame(rho, install("sparseX")))) {
    } else {
	RX = VAR_dMatrix_x(rho, lme4_RXSym, p, p);
	RZX = VAR_dMatrix_x(rho, lme4_RZXSym, q, p);
	X = VAR_dMatrix_x(rho, lme4_XSym, n, p);
	ZtX = VAR_dMatrix_x(rho, install("ZtX"), q, p);
    }
}

/**
 * Evaluate the deviance, possibly using AGQ.
 *
 * @return the updated deviance or approximate deviance
 */
double merenv::update_dev()
{
    return (double)0;
}
/* Externally callable functions */

/**
 * Evaluate the deviance or REML criterion
 *
 * @param x pointer to an mer environment
 *
 * @return deviance value
 */
SEXP lmerenv_deviance(SEXP x, SEXP thnew) {
    return ScalarReal(merenv(x).update_dev());
}
