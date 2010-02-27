#include "lme4utils.h"
#include "lme4utils.hpp"

/**
 * Generate zeros and weights of Hermite polynomial of order N, for
 * the AGQ method.
 *
 * Derived from Fortran code in package 'glmmML'
 *
 * @param N order of the Hermite polynomial
 * @param x zeros of the polynomial, abscissas for AGQ
 * @param w weights used in AGQ
 */

static void internal_ghq(int N, double *x, double *w)
{
    const double GHQ_EPS = 1e-15;
    const int GHQ_MAXIT = 40;
    int NR, IT, I, K, J;
    double Z = 0, HF = 0, HD = 0;
    double Z0, F0, F1, P, FD, Q, WP, GD, R, R1, R2;
    double HN = 1/(double)N;
    double *X = Calloc(N + 1, double), *W = Calloc(N + 1, double);

    for(NR = 1; NR <= N / 2; NR++){
	if(NR == 1)
	    Z = -1.1611 + 1.46 * sqrt((double)N);
	else
	    Z -= HN * (N/2 + 1 - NR);
	for (IT = 0; IT <= GHQ_MAXIT; IT++) {
	    Z0 = Z;
	    F0 = 1.0;
	    F1 = 2.0 * Z;
	    for(K = 2; K <= N; ++K){
		HF = 2.0 * Z * F1 - 2.0 * (double)(K - 1.0) * F0;
		HD = 2.0 * K * F1;
		F0 = F1;
		F1 = HF;
	    }
	    P = 1.0;
	    for(I = 1; I <= NR-1; ++I){
		P *= (Z - X[I]);
	    }
	    FD = HF / P;
	    Q = 0.0;
	    for(I = 1; I <= NR - 1; ++I){
		WP = 1.0;
		for(J = 1; J <= NR - 1; ++J){
		    if(J != I) WP *= ( Z - X[J] );
		}
		Q += WP;
	    }
	    GD = (HD-Q*FD)/P;
	    Z -= (FD/GD);
	    if (fabs((Z - Z0) / Z) < GHQ_EPS) break;
	}

	X[NR] = Z;
	X[N+1-NR] = -Z;
	R=1.0;
	for(K = 1; K <= N; ++K){
	    R *= (2.0 * (double)K );
	}
	W[N+1-NR] = W[NR] = 3.544907701811 * R / (HD*HD);
    }

    if( N % 2 ){
	R1=1.0;
	R2=1.0;
	for(J = 1; J <= N; ++J){
	    R1=2.0*R1*J;
	    if(J>=(N+1)/2) R2 *= J;
	}
	W[N/2+1]=0.88622692545276*R1/(R2*R2);
	X[N/2+1]=0.0;
    }

    dble_cpy(x, X + 1, N);
    dble_cpy(w, W + 1, N);

    if(X) Free(X);
    if(W) Free(W);
}

/**
 * Return zeros and weights of Hermite polynomial of order n as a list
 *
 * @param np pointer to a scalar integer SEXP
 * @return a list with two components, the abscissas and the weights.
 *
 */
SEXP lme4_ghq(SEXP np)
{
    int n = asInteger(np);
    SEXP ans = PROTECT(allocVector(VECSXP, 2));

    if (n < 1) n = 1;
    SET_VECTOR_ELT(ans, 0, allocVector(REALSXP, n ));
    SET_VECTOR_ELT(ans, 1, allocVector(REALSXP, n ));
    
    internal_ghq(n, REAL(VECTOR_ELT(ans, 0)), REAL(VECTOR_ELT(ans, 1)));
    UNPROTECT(1);
    return ans;
}

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
    CHM_SP ans = (CHM_SP) R_alloc(1, sizeof(cholmod_sparse));
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

/**
 * Return object with name nm in environment rho in a freshly
 * allocated CHM_SP structure checking on the number of rows and
 * columns.  Values of 0 for nrow or ncol skip the check.
 */
CHM_DN VAR_CHM_DN(SEXP rho, SEXP nm, int nrow, int ncol)
{
    CHM_DN ans = (CHM_DN) R_alloc(1, sizeof(cholmod_dense));
    M_as_cholmod_dense(ans, findVarBound(rho, nm));
    if (nrow && ((int)(ans->nrow)) != nrow)
	error(_("Number of rows of %s is %d, should be %d"),
	      CHAR(PRINTNAME(nm)), ans->nrow, nrow);
    if (ncol && ((int)(ans->ncol)) != ncol)
	error(_("Number of columns of %s is %d, should be %d"),
	      CHAR(PRINTNAME(nm)), ans->ncol, ncol);
    return ans;
}

/**
 * Return object with name nm in environment rho in a freshly
 * allocated CHM_SP structure checking on the number of rows and
 * columns.  Values of 0 for nrow or ncol skip the check.
 */
CHM_FR VAR_CHM_FR(SEXP rho, SEXP nm, int n)
{
    CHM_FR ans = (CHM_FR) R_alloc(1, sizeof(cholmod_factor));
    M_as_cholmod_factor(ans, findVarBound(rho, nm));
    if (n && ((int)(ans->n)) != n)
	error(_("Size of %s is %d, should be %d"),
	      CHAR(PRINTNAME(nm)), ans->n, n);
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
