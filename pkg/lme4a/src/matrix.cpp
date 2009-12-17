#include "matrix.hpp"

CHM_r* CHM_rd::solveCHM_FR(CHM_FR L, int sys) {
    return new CHM_rd(M_cholmod_solve(sys, L, A, &c));
}

CHM_rd::CHM_rd(SEXP x) {
    A = new cholmod_dense;
    M_as_cholmod_dense(A, x);
}

void CHM_rd::drmult(int transpose, double alpha, double beta,
		    CHM_DN X, CHM_DN Y){
    int ar = A->nrow, ac = A->ncol, xr = X->nrow, xc = X->ncol,
	yr = Y->nrow, yc = Y->ncol;
    if (transpose) {
	if (!(ac == yr && ar == xr && xc == yc))
	    error(_("Matrices not conformable for ddrmult"));
	F77_CALL(dgemm)("T", "N", &yr, &yc, &xr, &alpha,
			(double*)A->x, &ar, (double*)X->x, &xr,
			&beta, (double*)Y->x, &xc);
    } else {
	if (!(ar == yr && ac == xr && xc == yc))
	    error(_("Matrices not conformable for ddrmult"));
	F77_CALL(dgemm)("N", "N", &yr, &yc, &xr, &alpha,
			(double*)A->x, &yr, (double*)X->x, &xr,
			&beta, (double*)Y->x, &yr);
    }
}

CHM_rs::CHM_rs(SEXP x) {
    A = new cholmod_sparse;
    M_as_cholmod_sparse(A, x, (Rboolean)TRUE, (Rboolean)FALSE);
}

void CHM_rs::drmult(int transpose, double alpha, double beta,
		    CHM_DN X, CHM_DN Y){
    M_cholmod_sdmult(A, transpose, &alpha, &beta, X, Y, &c);
}

CHM_r* CHM_rs::solveCHM_FR(CHM_FR L, int sys) {
    return new CHM_rs(M_cholmod_spsolve(sys, L, A, &c));
}

Cholesky_rd::Cholesky_rd(SEXP x) {
    if (!(IS_S4_OBJECT(x)))
	error(_("S4 object expected but not provided"));
// FIXME: This check should be changed to an S4 "is" check, which
// should be available in Rinternals.h but isn't.
    if (strcmp(CHAR(asChar(getAttrib(x, R_ClassSymbol))),
	       "Cholesky") != 0)
	error(_("Object must be of class \"Cholesky\""));
    uplo = CHAR(asChar(GET_SLOT(x, install("uplo"))));
    int *dims = INTEGER(GET_SLOT(x, lme4_DimSym));
    n = dims[0];
    if (dims[1] != n)
	error(_("Cholesky object must be a square matrix"));
    X = REAL(GET_SLOT(x, lme4_xSym));
}

int Cholesky_rd::update(CHM_r *A) {
    return 0;
}

int Cholesky_rs::update(CHM_r *A) {
    return 0;
}

Cholesky_rs::Cholesky_rs(SEXP x) {
    F = new cholmod_factor;
    M_as_cholmod_factor(F, x);
}

dpoMatrix::dpoMatrix(SEXP x) {
    if (!(IS_S4_OBJECT(x)))
	error(_("S4 object expected but not provided"));
// FIXME: This check should be changed to an S4 inherits check, which
// should be available in Rinternals.h but isn't.
    if (strcmp(CHAR(asChar(getAttrib(x, R_ClassSymbol))),
	       "dpoMatrix") != 0)
	error(_("Object must be of class \"Cholesky\""));
    uplo = CHAR(asChar(GET_SLOT(x, install("uplo"))));
    int *dims = INTEGER(GET_SLOT(x, lme4_DimSym));
    n = dims[0];
    if (dims[1] != n)
	error(_("Cholesky object must be a square matrix"));
    X = REAL(GET_SLOT(x, lme4_xSym));
    factors = GET_SLOT(x, install("factors"));
    if (LENGTH(factors) && !isNewList(factors))
	error(_("\"factors\" slot should be a list"));
}
