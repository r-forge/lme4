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

void CHM_rd::copy_contents(CHM_r *SRC) {
    CHM_DN src = (dynamic_cast<CHM_rd*>(SRC))->A;
    int m = (int)A->nrow, n = (int)A->ncol;
    if (!(m == (int)src->nrow && n == (int)src->ncol))
	error(_("copy_contents mismatch: src is %d by %d, dest is %d by %d"),
	      (int)src->nrow, (int)src->ncol);
    dble_cpy((double*)A->x, (double*)src->x, m * n);
}

CHM_r* CHM_rd::crossprod_SP(CHM_SP Lam) {
    double one = 1, zero = 0;
    CHM_rd *ans = new CHM_rd(M_cholmod_copy_dense(A, &c));
    M_cholmod_sdmult(Lam, 1/*transpose*/, &one, &zero, A, ans->A, &c);
    return ans;
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

///replace current contents by factorization of alpha*crossprod(*AA) + beta*(*CC)
void Cholesky_rd::downdate(CHM_r *AA, double alpha,
			   dMatrix *CC, double beta) {
    if (n > 0) {
	int info;
	dpoMatrix *C = dynamic_cast<dpoMatrix*>(CC);

	if (strcmp(uplo, C->uplo) != 0)
	    error(_("In %s uplo attributes of %s and %s must agree"),
		  "Cholesky_rd::downdate", "X");
	if (C->n != n)
	    error(_("In %s %s must be %d by %d"),
		  "Cholesky_rd::downdate", "C", n, n);
	F77_CALL(dlacpy)(uplo, &n, &n, C->X, &n, X, &n);
	
	CHM_DN A = (dynamic_cast<CHM_rd*>(AA))->A;
	if ((int)A->ncol != n)
	    error(_("In %s %s must be %d by %d"),
		  "Cholesky_rd::downdate", "A", n, A->nrow);
	int k = (int)A->nrow;
	if (k > 0) {
	    F77_CALL(dsyrk)(uplo, "T", &n, &k, &alpha, (double*)A->x, &k,
			    &beta, X, &n);
	}
	F77_CALL(dpotrf)(uplo, &n, X, &n, &info);
	if (info)
	    error(_("In %s, %s returned error code %d on size %d matrix"),
		  "Cholesky_rd::downdate", "dpotrf", info, n);
    }
}

CHM_DN Cholesky_rd::solveA(CHM_DN rhs) {
    int info, nrhs = (int)rhs->ncol;
    CHM_DN ans = M_cholmod_copy_dense(rhs, &c);
    F77_CALL(dpotrs)(uplo, &n, &nrhs, X, &n,
		     (double*)ans->x, &n, &info);
    if (!info)
	error(_("dpotrs in Cholesky_rd::solveA returned error code %d"),
	      info);
    return ans;
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

void CHM_rs::copy_contents(CHM_r *SRC) {
    CHM_SP src = (dynamic_cast<CHM_rs*>(SRC))->A;
    int m = (int)A->nrow, n = (int)A->ncol;
    if (!(m == (int)src->nrow && n == (int)src->ncol))
	error(_("copy_contents mismatch: src is %d by %d, dest is %d by %d"),
	      (int)src->nrow, (int)src->ncol);
    int *di = (int*)(A->i), *si = (int*)(src->i),
	*dp = (int*)(A->p), *sp = (int*)(src->p);
    double *dx = (double*)(A->x), *sx = (double*)(src->x); 
    
    for (int j = 0; j <= n; j++)
	if (dp[j] != sp[j])
	    error(_("copy_contents:dest p[%d] not consistent with src"), j);
    for (int j = 0; j < dp[n]; j++) {
	if (di[j] != si[j])
	    error(_("copy_contents: dest i[%d] not consistent with src"), j);
	dx[j] = sx[j];
    }
}

CHM_r* CHM_rs::crossprod_SP(CHM_SP Lam) {
    CHM_SP Lamtr = M_cholmod_transpose(Lam, TRUE/*values*/, &c);
    CHM_rs *ans = new CHM_rs(
	M_cholmod_ssmult(Lamtr, A, 0/*stype*/, TRUE/*values*/,
			 TRUE/*sorted*/, &c));
    M_cholmod_free_sparse(&Lamtr, &c);
    return ans;
}

Cholesky_rs::Cholesky_rs(SEXP x) {
    F = new cholmod_factor;
    M_as_cholmod_factor(F, x);
}

int Cholesky_rs::update(CHM_r *AA) {
    CHM_SP A = (dynamic_cast<CHM_rs*>(AA))->A;
    return M_cholmod_factorize(A, F, &c);
}

void Cholesky_rs::downdate(CHM_r *AA, double alpha,
			   dMatrix *CC, double beta) {
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
}

CHM_DN Cholesky_rs::solveA(CHM_DN rhs) {
    return M_cholmod_solve(CHOLMOD_A, F, rhs, &c);
}

dpoMatrix::dpoMatrix(SEXP x) {
    if (!(IS_S4_OBJECT(x)))
	error(_("S4 object expected but not provided"));
// FIXME: This check should be changed to an S4 inherits check, which
// should be available in Rinternals.h but isn't.
    if (strcmp(CHAR(asChar(getAttrib(x, R_ClassSymbol))),
	       "dpoMatrix") != 0)
	error(_("Object must be of class \"dpoMatrix\""));
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
