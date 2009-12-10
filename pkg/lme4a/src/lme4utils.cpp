#include "Matrix.h"

#include "lme4utils.hpp"

SEXP CHM_SP2SEXP(CHM_SP A, const char *cls)
{ 
    return CHM_SP2SEXP(A, cls, (char*)NULL, (char*)NULL);
}

SEXP CHM_SP2SEXP(CHM_SP A, const char *cls, const char *uplo)
{ 
    return CHM_SP2SEXP(A, cls, uplo, (char*)NULL);
}

SEXP CHM_SP2SEXP(CHM_SP A, const char *cls, const char *uplo, const char *diag)
{
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS(cls)));
    int *dims = INTEGER(ALLOC_SLOT(ans, lme4_DimSym, INTSXP, 2));
    dims[0] = (int)(A->nrow);
    dims[1] = (int)(A->ncol);
    int nnz = M_cholmod_nnz(A, &c);
    Memcpy(INTEGER(ALLOC_SLOT(ans, lme4_pSym, INTSXP, A->ncol + 1)),
	   (int*)(A->p), A->ncol + 1);
    Memcpy(INTEGER(ALLOC_SLOT(ans, lme4_iSym, INTSXP, nnz)),
	   (int*)(A->i), nnz);
    Memcpy(REAL(ALLOC_SLOT(ans, lme4_xSym, REALSXP, nnz)),
	   (double*)(A->x), nnz);
    if (uplo) SET_SLOT(ans, install("uplo"), mkString(uplo));
    if (diag) SET_SLOT(ans, install("diag"), mkString(diag));
    UNPROTECT(1);
    return ans;
}

double *VAR_REAL_NULL(SEXP rho, SEXP nm, int len, int nullOK, int absentOK)
{
    const char *pn = CHAR(PRINTNAME(nm));
    SEXP var = findVarInFrame(rho, nm);
    if (var == R_UnboundValue) {
	if (absentOK) {
	    defineVar(nm, allocVector(REALSXP, len), rho);
	    return(REAL(findVarInFrame(rho, nm)));
	} else error(_("object named '%s' not found in environment"), pn);
    }
    int ll = 0;			// -Wall

    if (var == R_NilValue || !(ll  = LENGTH(var))) {
	if (nullOK || !len) return (double*) NULL;
	error(_("numeric object '%s' may not have length 0"), pn);
    }
    if (len && ll != len)
	error(_("Expected numeric object '%s' to be length %d, got %d"),
	      pn, len, ll);
    if (!isReal(var))
	error(_("numeric object '%s' not found in env"), pn);
    return REAL(var);
}

double *VAR_REAL_NULL(SEXP rho, SEXP nm, int len, int nullOK)
{
    return VAR_REAL_NULL(rho, nm, len, nullOK, FALSE);
}

double *VAR_REAL_NULL(SEXP rho, SEXP nm, int len)
{
    return VAR_REAL_NULL(rho, nm, len, FALSE, FALSE);
}

/**
 * Return object with name nm in environment rho in a freshly
 * allocated CHM_SP structure checking on the number of rows and
 * columns.  Values of 0 for nrow or ncol skip the check.
 */
CHM_SP VAR_CHM_SP(SEXP rho, SEXP nm, int nrow, int ncol)
{
    CHM_SP ans = new cholmod_sparse;
    M_as_cholmod_sparse(ans, findVarBound(rho, nm),
			(Rboolean)TRUE, (Rboolean)FALSE);
    if (nrow && ((int)(ans->nrow)) != nrow)
	error(_("Number of rows of %s is %d, should be %d"),
	      CHAR(PRINTNAME(nm)), ans->nrow, nrow);
    if (ncol && ((int)(ans->ncol)) != ncol)
	error(_("Number of columns of %s is %d, should be %d"),
	      CHAR(PRINTNAME(nm)), ans->ncol, ncol);
    return ans;
}

double *VAR_dMatrix_x(SEXP rho, SEXP nm, int nrow, int ncol) {    
    SEXP var = findVarBound(rho, nm);

    // FIXME: Should check here to ensure that the object "is" a dMatrix
    int *dims = INTEGER(GET_SLOT(var, lme4_DimSym));
    if (dims[0] != nrow || dims[1] != ncol)
	error(_("object named '%s' should be %d by %d"),
	      CHAR(PRINTNAME(nm)), nrow, ncol);
    return REAL(GET_SLOT(var, lme4_xSym));
}

CHM_SP CHM_SP_copy_in_place(CHM_SP dest, CHM_SP src) {
    size_t m = dest->nrow, n = dest->ncol;
    int nnzd = M_cholmod_nnz(dest, &c);
    if (m != src->nrow ||
	n != src->ncol ||
	dest->xtype != src->xtype ||
	dest->stype != src->stype ||
	!dest->packed || !src->packed ||
	nnzd != M_cholmod_nnz(src, &c) ||
	!compare_int_vecs((int*)dest->p, (int*)src->p, n+1) ||
	!compare_int_vecs((int*)dest->i, (int*)src->i, nnzd))
	error(_("incompatible CHM_SP objects for copy"));
    dble_cpy((double*)dest->x, (double*)src->x, nnzd);
    return dest;
}

SEXP getListElement(SEXP list, SEXP names, const char *str)
{
    SEXP elmt = (SEXP) NULL;
    const char *tempChar;
    int n = LENGTH(list);

    for (int i = 0; i < n; i++) {
	tempChar = CHAR(STRING_ELT(names, i)); /* ASCII only */
	if( strcmp(tempChar,str) == 0) {
	    elmt = VECTOR_ELT(list, i);
	    break;
	}
    }
    return elmt;
}

