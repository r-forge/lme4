#include <R.h>
#include <Rdefines.h>
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
