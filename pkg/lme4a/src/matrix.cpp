#include "matrix.hpp"

CHM_r* CHM_rd::solveCHM_FR(CHM_FR L, int sys) {
    return new CHM_rd(M_cholmod_solve(sys, L, A, &c));
}

void CHM_rd::show(const char *nm) {
    Rprintf("CHM_rd %s: nr = %d, nc = %d, xtype = %d\n",
	    nm, A->nrow, A->ncol, A->xtype);
}

void CHM_rd::drmult(int transpose, double alpha, double beta,
		    CHM_DN X, CHM_DN Y){
    int ar = A->nrow, ac = A->ncol, xr = X->nrow, xc = X->ncol,
	yr = Y->nrow, yc = Y->ncol;
    if (xr > 0 && yr > 0) {
	if (transpose) {
	    if (!(ac == yr && ar == xr && xc == yc))
		error(_("Matrices not conformable for ddrmult"));
	    F77_CALL(dgemm)("T", "N", &yr, &yc, &xr, &alpha,
			    (double*)A->x, &ar, (double*)X->x, &xr,
			    &beta, (double*)Y->x, &yr);
	} else {
	    if (!(ar == yr && ac == xr && xc == yc))
		error(_("Matrices not conformable for ddrmult"));
	    F77_CALL(dgemm)("N", "N", &yr, &yc, &xc, &alpha,
			    (double*)A->x, &yr, (double*)X->x, &xr,
			    &beta, (double*)Y->x, &yr);
	}
	return;
    }
    if (beta != 1) {
	int ny = yr * yc;
	double *yp = (double*)Y->x;
	for (int j = 0; j < ny; j++) yp[j] *= beta;
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
    static double one[] = {1,0}, zero[] = {0,0};
    CHM_rd *ans = new CHM_rd(M_cholmod_copy_dense(A, &c));
    M_cholmod_sdmult(Lam, 1/*transpose*/, one, zero, A, ans->A, &c);
    return ans;
}

CHM_r* CHM_rd::AtA() {
    static double one[] = {1,0}, zero[] = {0,0};
    int nc = (int)A->ncol, nr = (int)A->nrow;
    CHM_rd *ans = new CHM_rd(M_cholmod_allocate_dense(nc, nc, nc,
						      CHOLMOD_REAL, &c));
    if (nc > 0 && nr > 0) 
	F77_CALL(dsyrk)("U", "T", &nc, &nr, one, (double*)A->x, &nr,
			zero, (double*)((ans->A)->x), &nc);
    return ans;
}

Cholesky_rd::Cholesky_rd(SEXP x, int nn) {
    if (!(IS_S4_OBJECT(x)))
	error(_("S4 object expected but not provided"));
// FIXME: This check should be changed to an S4 "is" check, which
// should be available in Rinternals.h but isn't.
    if (strcmp(CHAR(asChar(getAttrib(x, R_ClassSymbol))),
	       "Cholesky") != 0)
	error(_("Object must be of class \"Cholesky\""));
    uplo = CHAR(asChar(GET_SLOT(x, install("uplo"))));
    int *dims = INTEGER(GET_SLOT(x, lme4_DimSym));
    n = nn;
    if (dims[0] != n || dims[1] != n)
	error(_("Cholesky object must be a square matrix of size %d"));
    X = REAL(GET_SLOT(x, lme4_xSym));
}

///replace current contents by factorization of matrix B
int Cholesky_rd::update(CHM_r *AA) {
    int info = 0;
    CHM_DN A = (dynamic_cast<CHM_rd*>(AA))->A;
    if (A->ncol != (size_t)n || A->nrow != (size_t)n)
	error(_("%s: matrix to factor must be %d by %d"),
	      "Cholesky_rd::update", n, n);
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
	    error(_("%s: uplo attributes of %s and %s must agree"),
		  "Cholesky_rd::downdate", "X", "C");
	if (C->n != n)
	    error(_("%s: %s must be %d by %d"),
		  "Cholesky_rd::downdate", "C", n, n);
	F77_CALL(dlacpy)(uplo, &n, &n, C->X, &n, X, &n);
	
	CHM_DN A = (dynamic_cast<CHM_rd*>(AA))->A;
	if ((int)A->ncol != n)
	    error(_("%s: %s must be %d by %d"),
		  "Cholesky_rd::downdate", "A", n, A->nrow);
	int k = (int)A->nrow;
	if (k > 0) {
	    F77_CALL(dsyrk)(uplo, "T", &n, &k, &alpha, (double*)A->x, &k,
			    &beta, X, &n);
	} else {
	    for (int i = 0; i < n * n; i++) X[i] *= beta;
	}
	F77_CALL(dpotrf)(uplo, &n, X, &n, &info);
	if (info)
	    error(_("%s: %s returned error code %d on size %d matrix"),
		  "Cholesky_rd::downdate", "dpotrf", info, n);
    }
}

CHM_DN Cholesky_rd::solveA(CHM_DN rhs) {
    int info, nrhs = (int)rhs->ncol;
    CHM_DN ans = M_cholmod_copy_dense(rhs, &c);
    F77_CALL(dpotrs)(uplo, &n, &nrhs, X, &n,
		     (double*)ans->x, &n, &info);
    if (info)
	error(_("dpotrs in Cholesky_rd::solveA returned error code %d"),
	      info);
    return ans;
}

double Cholesky_rd::ldet2() {
    double ans = 0;
    for (int i = 0; i < n; i++) {
	double dxi = X[i * (n + 1)];
	ans += log(dxi * dxi);
    }
    return ans;
}

void CHM_rs::show(const char *nm) {
    int nc = A->ncol, nnz = ((int*)A->p)[A->ncol];
    Rprintf("CHM_SP %s: nr = %d, nc = %d, nzmax = %d, nnz = %d\n",
	    nm, A->nrow, nc, A->nzmax, nnz);
    Rprintf("  xtype = %d, stype = %d, sorted = %d, packed = %d\n  p =",
	    A->xtype, A->stype, A->sorted, A->packed);
    int *pp = (int*)A->p, *ii = (int*)A->i;
    double *xx = (double*)A->x;
    for (int i = 0; i < 9 && i <= nc ; i++) Rprintf(" %d,", pp[i]);
    Rprintf(" %d\n  i =", pp[9]);
    for (int i = 0; i < 9 && i < nnz; i++) Rprintf(" %d,", ii[i]);
    Rprintf(" %d\n  x =", ii[9]);
    for (int i = 0; i < 5 && i < nnz; i++) Rprintf(" %g,", xx[i]);
    Rprintf(" %g\n", xx[5]);
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
	error(_("%s mismatch: src is %d by %d, dest is %d by %d"),
	      "CHM_rs::copy_contents",
	      (int)src->nrow, (int)src->ncol, m, n);
    size_t nnz = ((int*)src->p)[n];
    if (A->nzmax < nnz)
	error(_("%s mismatch: src has %d non-zeros, dest->nzmax = %d"),
	      "CHM_rs::copy_contents", nnz, ((int*)A->p)[n]);
    if (A->itype != src->itype || A->xtype != src->xtype ||
	A->dtype != src->dtype)
	error(_("%s mismatch: on itype, xtype or dtype"),
	      "CHM_rs::copy_contents");
    int_cpy((int*)A->p, (int*)src->p, n + 1);
    int_cpy((int*)A->i, (int*)src->i, nnz);
    dble_cpy((double*)A->x, (double*)src->x, nnz);
    A->stype = src->stype; A->sorted = src->sorted;
    if (!src->packed) error(_("%s: only packed src can be copied"),
			    "CHM_rs::copy_contents");
    A->packed = 1;
}

CHM_r* CHM_rs::crossprod_SP(CHM_SP Lam) {
    CHM_SP Lamtr = M_cholmod_transpose(Lam, TRUE/*values*/, &c);
    CHM_rs *ans = new CHM_rs(
	M_cholmod_ssmult(Lamtr, A, 0/*stype*/, TRUE/*values*/,
			 TRUE/*sorted*/, &c));
    M_cholmod_free_sparse(&Lamtr, &c);
    return ans;
}

CHM_r* CHM_rs::AtA() {
    CHM_SP At = M_cholmod_transpose(A, TRUE/*values*/, &c);
    CHM_rs *ans = new CHM_rs(
	M_cholmod_ssmult(At, A, 1/*stype*/, TRUE/*values*/,
			 TRUE/*sorted*/, &c));
    M_cholmod_free_sparse(&At, &c);
    return ans;
}

int Cholesky_rs::update(CHM_r *AA) {
    CHM_SP A = (dynamic_cast<CHM_rs*>(AA))->A;
    return M_cholmod_factorize(A, F, &c);
}

double Cholesky_rs::ldet2() {
    return M_chm_factor_ldetL2(F);
}

void Cholesky_rs::downdate(CHM_r *AA, double alpha,
			   dMatrix *CC, double beta) {
    CHM_SP A = (dynamic_cast<CHM_rs*>(AA))->A,
	C = (dynamic_cast<CHM_rs*>(CC))->A;
    size_t n = F->n;

    if (C->ncol != n || C->nrow != n || !C->stype)
	error(_("In Cholesky downdate C must be symmetric %d by %d"), n, n);
    if (n > 0) {
//	check_sparse(A, "A");
//	check_sparse(C, "C");
	CHM_SP At = M_cholmod_transpose(A, 1/*values*/, &c);
//	check_sparse(At, "At");
	CHM_SP AtA = M_cholmod_ssmult(At, A, 1/*stype*/, 1/*values*/,
				      0/*sorted*/, &c);
//	check_sparse(AtA, "AtA");
	M_cholmod_free_sparse(&At, &c);
	CHM_SP sum = M_cholmod_add(AtA, C, &alpha, &beta,
				   1/*values*/, 1/*sorted*/, &c);
	M_cholmod_free_sparse(&AtA, &c);
//	check_sparse(sum, "sum");
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
