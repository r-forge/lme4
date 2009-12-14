#include "lme4utils.h"		// for SEXP, F77_CALL, etc.
#include "matrix.hpp"

ddiMatrix::ddiMatrix(double *x, int nr, int nc) {
    contents = x;
    nro = nr;
    nco = nc;
}

ddiMatrix::ddiMatrix(double *x, int *dims) {
    contents = x;
    nro = dims[0];
    nco = dims[1];
}

dgeMatrix::dgeMatrix(double *x, int nr, int nc) {
    contents = x;
    nro = nr;
    nco = nc;
}

dgeMatrix::dgeMatrix(double *x, int *dims) {
    contents = x;
    nro = dims[0];
    nco = dims[1];
}

void ddiMatrix::rdprod(double *ans, const double *rhs,
		       int nrhs, double alpha, double beta,
		       int transa) {
    for (int j = 0; j < nrhs; j++) {
	for (int i = 0; i < nco; i++) {
	    int ind = j * nro + i;
	    ans[ind] = alpha * contents[i] * rhs[ind] + beta * rhs[ind];
	}
    }
}

void dgeMatrix::rdprod(double *ans, const double *rhs,
		       int nrhs, double alpha, double beta,
		       int transa) {
    if (nco > 0 && nrhs > 0) {
	const char *tra = transa ? "T" : "N";
	int nrop = transa ? nco : nro;
	int ncop = transa ? nro : nco;	
	F77_CALL(dgemm)(tra, "N", &nrop, &nrhs, &ncop, &alpha,
			contents, &nro, rhs, &nrop, &beta, ans, &ncop);
    }
}

void dgCMatrix::rdprod(double *ans, const double *rhs,
		       int nrhs, double alpha, double beta,
		       int transa) {
}
