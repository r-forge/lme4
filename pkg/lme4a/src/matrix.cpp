#include "matrix.hpp"
CHM_rd::CHM_rd(SEXP x) {
    A = new cholmod_dense;
    M_as_cholmod_dense(A, x);
}

CHM_rs::CHM_rs(SEXP x) {
    A = new cholmod_sparse;
    M_as_cholmod_sparse(A, x, (Rboolean)TRUE, (Rboolean)FALSE);
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
			(double*)A->x, &ar, (double*)X->x, &xr,
			&beta, (double*)Y->x, &xc);
    }
}

void CHM_rs::drmult(int transpose, double alpha, double beta,
		    CHM_DN X, CHM_DN Y){
    M_cholmod_sdmult(A, transpose, &alpha, &beta, X, Y, &c);
}
