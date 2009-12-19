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

///replace current contents by factorization of matrix B
int Cholesky_rd::update(CHM_r *AA) {
    int info = 0;
    CHM_DN A = (dynamic_cast<CHM_rd*>(AA))->A;
    if (A->ncol != (size_t)n || A->nrow != (size_t)n)
	error(_("in Cholesky update A must be %d by %d"), n, n);
    if (n > 0) {
	F77_CALL(dlacpy)(uplo, &n, &n, (double*)A->x, &n, X, &n);
	F77_CALL(dpotrf)(uplo, &n, X, &n, &info);
    }
    return info;
}

int Cholesky_rs::update(CHM_r *AA) {
    CHM_SP A = (dynamic_cast<CHM_rs*>(AA))->A;
    return M_cholmod_factorize(A, F, &c);
}

///replace current contents by factorization of alpha*B - beta*D'D
int Cholesky_rd::downdate(CHM_r *AA, double alpha, CHM_r *CC, double beta) {
    CHM_DN A = (dynamic_cast<CHM_rd*>(AA))->A,
	C = (dynamic_cast<CHM_rd*>(CC))->A;
    int info = 0;
    
    if ((int)C->ncol != n || (int)C->nrow != n || (int)A->nrow != n)
	error(_("In Cholesky downdate C be %d by %d and A must be %d by x"),
	      n, n, n);
    if (n > 0) {
	int k = (int)A->ncol;
	if (k > 0) {
	    F77_CALL(dlacpy)(uplo, &n, &n, (double*)C->x, &n, X, &n);
	    F77_CALL(dsyrk)(uplo, "T", &n, &k, &alpha, (double*)A->x, &n,
			    &beta, X, &n);
	}
	F77_CALL(dpotrf)(uplo, &n, X, &n, &info);
    }
    return info;
}

int Cholesky_rs::downdate(CHM_r *AA, double alpha, CHM_r *CC, double beta) {
    CHM_SP A = (dynamic_cast<CHM_rs*>(AA))->A,
	C = (dynamic_cast<CHM_rs*>(CC))->A;
    size_t n = F->n;

    if (C->ncol != n || C->nrow != n || !C->stype)
	error(_("In Cholesky downdate C be symmetric %d by %d"), n, n);
    if (n > 0) {
	CHM_SP At = M_cholmod_transpose(A, 1/*values*/, &c);
	CHM_SP AtA = M_cholmod_aat(At, (int*)NULL, (size_t)0,
				   -1/*numerical values*/, &c);
	M_cholmod_free_sparse(&At, &c);
	CHM_SP sum = M_cholmod_add(AtA, C, &alpha, &beta, 1/*values*/,
				   1/*sorted*/, &c);
	M_cholmod_free_sparse(&AtA, &c);
	M_cholmod_factorize(sum, F, &c);
	M_cholmod_free_sparse(&sum, &c);
    }
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
