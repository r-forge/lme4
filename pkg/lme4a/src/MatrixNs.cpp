#include "MatrixNs.h"
#include <R_ext/Lapack.h>

using namespace Rcpp;

namespace MatrixNs{
    Matrix::Matrix(S4 &xp) :
	Dimnames(xp.slot("Dimnames")) {
	IntegerVector Dim(xp.slot("Dim"));
	d_nrow = Dim[0];
	d_ncol = Dim[1];
    }

    int Matrix::nrow() const { return d_nrow; }
    int Matrix::ncol() const { return d_ncol; }

// Check this: Do the dspMatrix and dtpMatrix classes pass this check?
    ddenseMatrix::ddenseMatrix(S4 &xp) : dMatrix(xp) {
	if (!x.size() == d_nrow * d_ncol)
	    ::Rf_error("%s: Dim = (%d, %d) is inconsistent with x.size() = %d",
		       "ddenseMatrix::ddenseMatrix", d_nrow, d_ncol, x.size());
    }

    void dgeMatrix::dgemv(char Tr, double alpha,
			  NumericVector const &X, double beta,
			  NumericVector &Y) const {
	int i1 = 1;
	Trans TR(Tr);
	char tr = TR.TR;
	bool NTR = tr == 'N';
	if (X.size() != (NTR ? d_ncol : d_nrow) ||
	    Y.size() != (NTR ? d_nrow : d_ncol))
	    Rf_error("dgemv \"%c\", dim mismatch (%d, %d), X(%d), Y(%d)",
		     tr, d_nrow, d_ncol, X.size(), Y.size());
	F77_CALL(dgemv)(&tr, &d_nrow, &d_ncol, &alpha, x.begin(), &d_nrow,
			X.begin(), &i1, &beta, Y.begin(), &i1);
    }

    void dgeMatrix::dgemm(char TRA, char TRB,
			  double alpha, const dgeMatrix &B,
			  double beta, dgeMatrix &C) const {
	Trans TrA(TRA), TrB(TRB);
	char trA = TrA.TR, trB = TrB.TR;
	bool NTA = trA == 'N', NTB = trB == 'N';
//	Dimension Bd(B.Dim), Cd(C.Dim);
	int M = NTA ? d_nrow : d_ncol,
	    N = NTB ? B.ncol() : B.nrow(),
	    K = NTA ? d_ncol : d_nrow;
	int Bnr = B.nrow();
	if (NTB ? B.ncol() : B.nrow() != K || C.nrow() != M || C.ncol() != N)
	    Rf_error("dgemm \"%c,%c\", dim mismatch (%d, %d), (%d,%d), (%d,%d)",
		     trA, trB, d_nrow, d_ncol, B.nrow(), B.ncol(),
		     C.nrow(), C.ncol());
	F77_CALL(dgemm)(&trA, &trB, &M, &N, &K, &alpha, x.begin(), &d_nrow,
			B.x.begin(), &Bnr, &beta, C.x.begin(), &M);
    }

    void Cholesky::update(dgeMatrix const &A) {
//	Dimension Ad(A.Dim);
	if (d_nrow != A.ncol())
	    Rf_error("%s dimension mismatch, (%d,%d) vs A(%d,%d)",
		     "Cholesky::update(dgeMatrix)", d_nrow, d_ncol,
		     A.nrow(), A.ncol());
	double alpha = 1., beta = 0.;
	int Anr = A.nrow();
	F77_CALL(dsyrk)(&(uplo.UL), "T", &d_nrow, &Anr, &alpha,
			A.x.begin(), &Anr, &beta, x.begin(), &d_nrow);
	int info;
	F77_CALL(dpotrf)(&(uplo.UL), &d_nrow, x.begin(), &d_nrow, &info);
	if (info)
	    Rf_error("Lapack routine %s returned error code %d",
		     "dpotrf", info);
    }

    void Cholesky::update(dpoMatrix const &A) {
//	Dimension Ad(A.Dim);
	if (d_nrow != A.nrow())
	    Rf_error("%s dimension mismatch, (%d,%d) vs A(%d,%d)",
		     "Cholesky::update(dpoMatrix)", d_nrow, d_ncol,
		     A.nrow(), A.ncol());
	uplo = A.uplo;
	std::copy(A.x.begin(), A.x.end(), x.begin());
	int info;
	F77_CALL(dpotrf)(&(uplo.UL), &d_nrow, x.begin(), &d_nrow, &info);
	if (info)
	    Rf_error("Lapack routine %s returned error code %d",
		     "dpotrf", info);
    }

    void Cholesky::update(Trans Tr, double alpha, const dgeMatrix &A,
			  double beta, const dsyMatrix &C) {
	const char tr = Tr.TR;
	const bool NTR = tr == 'N';
//	Dimension Ad(A.Dim), Cd(C.Dim);
	int Anr = A.nrow(), Anc = A.ncol(), Cnr = C.nrow(), Cnc = C.ncol();
	if (d_nrow != Cnr || NTR ? Anr : Anc != d_nrow)
	    Rf_error("%s(\"%c\") dimension mismatch, (%d,%d), A(%d,%d), C(%d,%d)",
		     "Cholesky::update(dpoMatrix, dgeMatrix)", tr,
		     d_nrow, d_ncol, Anr, Anc, Cnr, Cnc);
	uplo = C.uplo;
	std::copy(C.x.begin(), C.x.end(), x.begin());
	F77_CALL(dsyrk)(&(uplo.UL), &tr, &d_nrow, NTR ? &Anc : &Anr, &alpha,
			A.x.begin(), &Anr, &beta, x.begin(), &d_nrow);
	int info;
	F77_CALL(dpotrf)(&(uplo.UL), &d_nrow, x.begin(), &d_nrow, &info);
	if (info)
	    Rf_error("Lapack routine %s returned error code %d",
		     "dpotrf", info);
    }

    void Cholesky::dpotrs(NumericVector &v) const {
	int info, i1 = 1, vs = v.size();
	if (vs != d_nrow)
	    Rf_error("%s (%d, %d) dimension mismatch (%d, %d)",
		     "Cholesky::dpotrs", d_nrow, d_ncol, vs, 1);
	F77_CALL(dpotrs)(&uplo.UL, &d_nrow, &i1, x.begin(), &d_nrow,
			 v.begin(), &vs, &info);
	if (info)
	    Rf_error("Lapack routine %s returned error code %d",
		     "dpotrs", info);
    }

    double Cholesky::logDet2() {
	int nc = ncol(), stride = nrow() + 1;
	double *rx = x.begin(), ans = 0.;
	for (int i = 0; i < nc; i++, rx += stride)
	    ans += 2. * log(*rx);
	return ans;
    }

    chmFr::chmFr(S4 xp) : cholmod_factor()//, m_sexp(SEXP(xp))
    {
	CharacterVector cl(SEXP(xp.attr("class")));
	char *clnm = cl[0];
	if (!xp.is("CHMfactor"))
	    ::Rf_error("Class %s object passed to %s is not a %s",
		       clnm, "chmFr::chmFr", "CHMfactor");
	Dimension Dim(SEXP(xp.slot("Dim")));
	IntegerVector colcount(SEXP(xp.slot("colcount"))),
	    perm(SEXP(xp.slot("perm"))),
	    type(SEXP(xp.slot("type")));
	NumericVector X(SEXP(xp.slot("x")));

//	pp = (CHM_FR)NULL;
	minor = n = Dim[0];
	Perm = perm.begin();
	ColCount = colcount.begin();
	x = (void*)X.begin();
	ordering = type[0];
	is_ll = (type[1] ? 1 : 0);
	is_super = type[2];
	is_monotonic = (type[3] ? 1 : 0);
	xtype = CHOLMOD_REAL;
	itype = CHOLMOD_LONG;
	dtype = 0;  // CHOLMOD_DOUBLE
	z = (void*)NULL;
	const char* msg = "dCHMfactor with is_super == %s is not %s";
	if (is_super) {
	    IntegerVector
		Pi(SEXP(xp.slot("pi"))),
		SUPER(SEXP(xp.slot("super"))),
		S(SEXP(xp.slot("s")));
	    NumericVector PX(SEXP(xp.slot("px")));

	    if (!xp.is("dCHMsuper"))
		::Rf_error(msg, "TRUE", "dCHMsuper");
	    xsize = X.size();
	    ssize = S.size();
	    maxcsize = type[4];
	    maxesize = type[5];
	    nsuper = SUPER.size() - 1;
	    super = (void*)SUPER.begin();
	    pi = (void*)Pi.begin();
	    s = (void*)S.begin();
	    px = (void*)PX.begin();
	} else {
	    IntegerVector
		I(SEXP(xp.slot("i"))),
		NXT(SEXP(xp.slot("nxt"))),
		NZ(SEXP(xp.slot("nz"))),
		P(SEXP(xp.slot("p"))),
		PRV(SEXP(xp.slot("prv")));

	    if (!xp.is("dCHMsimpl"))
		::Rf_error(msg, "FALSE", "dCHMsimpl");
	    nzmax = X.size();
	    p = (void*)P.begin();
	    i = (void*)I.begin();
	    x = (void*)X.begin();
	    nz = (void*)NZ.begin();
	    next = (void*)NXT.begin();
	    prev = (void*)PRV.begin();
	}
    }

    // chmFr::chmFr(CHM_FR xp) : cholmod_factor(), pp(xp) {
    // 	n = xp->n;
    // 	minor = xp->minor;
    // 	Perm = xp->Perm;
    // 	ColCount = xp->ColCount;
    // 	nzmax = xp->nzmax;
    // 	p = xp->p;
    // 	i = xp->i;
    // 	x = xp->x;
    // 	z = xp->z;
    // 	nz = xp->nz;
    // 	next = xp->next;
    // 	prev = xp->prev;
    // 	nsuper = xp->nsuper;
    // 	ssize = xp->ssize;
    // 	xsize = xp->xsize;
    // 	maxcsize = xp->maxcsize;
    // 	maxesize = xp->maxesize;
    // 	super = xp->super;
    // 	pi = xp->pi;
    // 	px = xp->px;
    // 	s = xp->s;
    // 	ordering = xp->ordering;
    // 	is_ll = xp->is_ll;
    // 	is_super = xp->is_super;
    // 	is_monotonic = xp->is_monotonic;
    // 	itype = xp->itype;
    // 	xtype = xp->xtype;
    // 	dtype = xp->dtype;
    // }
		    
    void chmFr::update(cholmod_sparse const &A, double Imult) {
	M_cholmod_factorize_p((const CHM_SP)&A, &Imult, (int*)NULL,
			      (size_t) 0, this, &c);
    }

    CHM_DN chmFr::solve(int sys, const CHM_DN b) const {
	return M_cholmod_solve(sys, (const CHM_FR)this, b, &c);
    }

    CHM_DN chmFr::solve(int sys, chmDn const &b) const {
	return M_cholmod_solve(sys, (const CHM_FR)this,
			       (const CHM_DN)&b, &c);
    }

    CHM_SP chmFr::spsolve(int sys, const CHM_SP b) const {
	return M_cholmod_spsolve(sys, (const CHM_FR)this, b, &c);
    }
    CHM_SP chmFr::spsolve(int sys, chmSp const &b) const {
	return M_cholmod_spsolve(sys, (const CHM_FR)this,
				 (const CHM_SP)&b, &c);
    }

    chmSp::chmSp(S4 xp) : cholmod_sparse()//, m_sexp(SEXP(xp))
    {
	CharacterVector cl(SEXP(xp.attr("class")));
	char *clnm = cl[0];
	if (!xp.is("CsparseMatrix"))
	    Rf_error("Class %s object passed to %s is not a %s",
		     clnm, "chmSp::chmSp", "CsparseMatrix");
//	pp = (CHM_SP)NULL;
	IntegerVector
	    Dim(SEXP(xp.slot("Dim"))),
	    pp(SEXP(xp.slot("p"))),
	    ii(SEXP(xp.slot("i")));
	nrow = Dim[0];
	ncol = Dim[1];
	p = (void*)pp.begin();
	i = (void*)ii.begin();
	dtype = 0;  // CHOLMOD_DOUBLE
	stype = 0;
	if (!xp.is("generalMatrix")) {
	    if (xp.is("symmetricMatrix")) {
		CharacterVector uplo(SEXP(xp.slot("uplo")));
		char *UL = uplo[0];
		stype = 1;
		if (*UL == 'L' || *UL == 'l') stype = -1;
	    }
	}
	nzmax = ii.size();
	itype = CHOLMOD_LONG;
	packed = (int)true;
	sorted = (int)true;
	xtype = -1;
	if (xp.is("dsparseMatrix")) {
	    NumericVector xx(SEXP(xp.slot("x")));
	    x = xx.begin();
	    xtype = CHOLMOD_REAL;
	}
	if (xp.is("nsparseMatrix")) xtype = CHOLMOD_PATTERN;
	if (xp.is("zsparseMatrix")) {
	    xtype = CHOLMOD_COMPLEX;
	    Rf_error("Not yet defined zsparseMatrix?");
	}
	if (xtype == -1) Rf_error("Unknown (logical?) sparse Matrix type");
    }

/// Wrap the CHM_SP so the chmSp destructor can call M_cholmod_free_sparse
    // chmSp::chmSp(CHM_SP xp) : cholmod_sparse(), pp(xp) {
    // 	nrow = xp->nrow;
    // 	ncol = xp->ncol;
    // 	nzmax = xp->nzmax;
    // 	p = xp->p;
    // 	i = xp->i;
    // 	nz = xp->nz;
    // 	x = xp->x;
    // 	z = xp->z;
    // 	stype = xp->stype;
    // 	itype = xp->itype;
    // 	xtype = xp->xtype;
    // 	dtype = xp->dtype;
    // 	packed = xp->packed;
    // 	sorted = xp->sorted;
    // }

    CHM_SP chmSp::transpose(int values) const {
	return M_cholmod_transpose((const CHM_SP)this, values, &c);
    }

    int chmSp::dmult(char tr, double alpha, double beta,
		     chmDn const &src, chmDn &dest) const {
	return M_cholmod_sdmult((const CHM_SP)this,
				Trans(tr).TR == 'T', &alpha,
				&beta, (const CHM_DN)(&src), &dest, &c);
    }

    void chmSp::update(cholmod_sparse const &nn) {
	size_t nnznn = M_cholmod_nnz((const CHM_SP)&nn, &c);
	if (nn.ncol != ncol || nnznn > nzmax || xtype != nn.xtype ||
	    itype != nn.itype || dtype != nn.dtype || packed != nn.packed)
	    Rf_error("%s: matrices not conformable", "chmSp::update");
	stype = nn.stype;
	nrow = nn.nrow;
	sorted = nn.sorted;
	int *pp = (int*)p, *nnp = (int*)(nn.p);
	std::copy(nnp, nnp + ncol + 1, pp);
	int *ip = (int*)i, *nni = (int*)(nn.i);
	std::copy(nni, nni + nnznn, ip);
	switch(xtype) {
	case CHOLMOD_PATTERN:
	    break;
	case CHOLMOD_REAL: {
	    double *xp = (double*)x, *nnx = (double*)(nn.x);
	    std::copy(nnx, nnx + nnznn, xp);
	    break;
	}
	default:
	    Rf_error("Unsupported CHOLMOD xtype in %s", "chmSp::update");
	}
    }

    CHM_SP chmSp::crossprod() const {
	CHM_SP t1 = this->transpose();
	CHM_SP t2 = ::M_cholmod_aat(t1, (int*)NULL, 0/*fsize*/,
				    xtype/*mode*/, &c);
	::M_cholmod_free_sparse(&t1, &c);
	t1 = ::M_cholmod_copy(t2, 1/*stype*/, xtype/*mode*/, &c);
	::M_cholmod_free_sparse(&t2, &c);
	return t1;
    }

    CHM_SP chmSp::crossprod(const CHM_SP B, int sorted) const {
	CHM_SP t1 = this->transpose();
	CHM_SP t2 = ::M_cholmod_ssmult(t1, B, 0/*stype*/, xtype/*values*/,
				       sorted, &c);
	::M_cholmod_free_sparse(&t1, &c);
	return t2;
    }

    CHM_SP chmSp::crossprod(chmSp const &B, int sorted) const {
	return crossprod((const CHM_SP)&B, sorted);
    }

    CHM_SP chmSp::tcrossprod() const {
	CHM_SP t1 = M_cholmod_aat((const CHM_SP)this, (int*)NULL,
				  0/*fsize*/, xtype/*mode*/, &c);
	CHM_SP t2 = ::M_cholmod_copy(t1, 1/*stype*/, xtype/*mode*/, &c);
	::M_cholmod_free_sparse(&t1, &c);
	return t2;
    }

    CHM_SP chmSp::tcrossprod(const CHM_SP B, int sorted) const {
	CHM_SP t1 = M_cholmod_transpose(B, xtype/*values*/, &c);
	CHM_SP t2 =
	    M_cholmod_ssmult((const CHM_SP)this, t1, 0/*stype*/,
			     xtype/*values*/, sorted, &c);
	M_cholmod_free_sparse(&t1, &c);
	return t2;
    }

    CHM_SP chmSp::tcrossprod(chmSp const &B, int sorted) const {
	return tcrossprod((const CHM_SP)&B, sorted);
    }

    CHM_SP chmSp::smult(chmSp const &B, int stype, int values,
			int sorted) const {
	return M_cholmod_ssmult((const CHM_SP)this, (const CHM_SP)&B,
				stype, values, sorted, &c);
    }

    // chmDn::chmDn(CHM_DN xp) : cholmod_dense(), pp(xp) {
    // 	nrow = pp->nrow;
    // 	ncol = pp->ncol;
    // 	nzmax = pp->nzmax;
    // 	d = pp->d;
    // 	x = pp->x;
    // 	z = pp->z;
    // 	xtype = pp->xtype;
    // 	dtype = pp->dtype;
    // }

    void chmDn::init(double *X, int r, int c) {
	z = (void*)NULL;
	xtype = CHOLMOD_REAL;
	dtype = 0;		// CHOLMOD_DOUBLE
	nrow = (size_t) r;
	ncol = (size_t) c;
	nzmax = (size_t) r * c;
	d = (size_t) r;
	x = (void*) X;
//	pp = (CHM_DN)NULL;
    }
}
