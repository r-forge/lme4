#include "MatrixNs.h"
#include <R_ext/Lapack.h>
#include <cctype>

namespace MatrixNs{
    char chkchar(char x, std::string allowed) {
	char X = toupper(x);
	if (std::find(allowed.begin(), allowed.end(), X) == allowed.end())
	    throw std::range_error("chkchar");
	return X;
    }

    char UpLo(char x) {
	return(chkchar(x, "UL"));
    }

    char UpLo(std::string const& x) {
	return UpLo(x[0]);
    }
    char UpLo(SEXP x) {
	return UpLo(*CHAR(Rf_asChar(x)));
    }

    char Diag(char x) {
	return(chkchar(x, "NU"));
    }
    char Diag(std::string const& x) {
	return Diag(x[0]);
    }
    char Diag(SEXP x) {
	return Diag(*CHAR(Rf_asChar(x)));
    }

    char Trans(char x) {
	return(chkchar(x, "TNC"));
    }
    char Trans(std::string const& x) {
	return Trans(x[0]);
    }
    char Trans(SEXP x) {
	return Trans(*CHAR(Rf_asChar(x)));
    }

    Matrix::Matrix(Rcpp::S4 &xp) :
	Dimnames(xp.slot("Dimnames")) {
	Rcpp::IntegerVector Dim(xp.slot("Dim"));
	d_nrow = Dim[0];
	d_ncol = Dim[1];
    }

    Matrix::Matrix(int nr, int nc)
	: Dimnames(2),
	  d_nrow(nr),
	  d_ncol(nc) {
	if (nr < 0 || nc < 0)
	    throw std::range_error("Matrix(nr,nc)");
    }
	
    int Matrix::nrow() const { return d_nrow; }
    int Matrix::ncol() const { return d_ncol; }

    dMatrix::dMatrix(Rcpp::S4& xp)
	: Matrix(xp),
	  d_x(SEXP(xp.slot("x"))) {
    }

    dMatrix::dMatrix(int nr, int nc, int nx)
	: Matrix(nr, nc),
	  d_x(nx) {
    }

    void dMatrix::setX(Rcpp::NumericVector const& nx) {
	if (nx.size() != d_x.size())
	    throw std::range_error("Size mismatch in setX");
	std::copy(nx.begin(), nx.end(), d_x.begin());
    }

    void dMatrix::setX(Rcpp::NumericMatrix const& mm) {
	setX(Rcpp::NumericVector(SEXP(mm)));
    }

// Check this: Do the dspMatrix and dtpMatrix classes pass this check?
    ddenseMatrix::ddenseMatrix(Rcpp::S4 &xp) : dMatrix(xp) {
	if (!d_x.size() == d_nrow * d_ncol)
	    ::Rf_error("%s: Dim = (%d, %d) is inconsistent with x.size() = %d",
		       "ddenseMatrix::ddenseMatrix", d_nrow, d_ncol, d_x.size());
    }

    ddenseMatrix::ddenseMatrix(int nr, int nc)
	: dMatrix(nr, nc, nr * nc) {
    }

    compMatrix::compMatrix(Rcpp::S4& xp)
	: factors(SEXP(xp.slot("factors"))) {
    }

    generalMatrix::generalMatrix(Rcpp::S4& xp)
	: compMatrix(xp) {
    }

    triangularMatrix::triangularMatrix(Rcpp::S4& xp)
	: d_ul(UpLo(SEXP(xp.slot("uplo")))),
	  d_di(Diag(SEXP(xp.slot("diag")))) {
    }

    triangularMatrix::triangularMatrix(char ul, char di)
	: d_ul(UpLo(ul)),
	  d_di(Diag(ul)) {
    }

    symmetricMatrix::symmetricMatrix(Rcpp::S4& xp)
	: d_ul(UpLo(SEXP(xp.slot("uplo")))) {
    }

    symmetricMatrix::symmetricMatrix(char ul)
	: d_ul(UpLo(ul)) {
    }

// Not sure why dgeMatrix needs Rcpp::S4 but others take Rcpp::S4& in constructor
// Maybe because of compMatrix extracting the List??
    dgeMatrix::dgeMatrix(Rcpp::S4 xp)
	: ddenseMatrix(xp),
	  generalMatrix(xp) {
    }

    dgeMatrix::dgeMatrix(int nr, int nc)
	: ddenseMatrix(nr, nc) {
    }

    int dgeMatrix::dmult(char Tr, double alpha, double beta, 
			 chmDn const &src, chmDn &dest) const {
	int i1 = 1;
	char tr = Trans(Tr);
	if (src.nrow == 1)
	    F77_CALL(dgemv)(&tr, &d_nrow, &d_ncol, &alpha, d_x.begin(),
			    &d_nrow, (double*)src.x, &i1, &beta,
			    (double*)dest.x, &i1);
	else {
	    bool NTR = tr == 'N';
	    int M = NTR ? d_nrow : d_ncol, N = src.ncol,
		K = NTR ? d_ncol : d_nrow;
	    F77_CALL(dgemm)(&tr, "N", &M, &N, &K, &alpha, d_x.begin(),
			    &d_nrow, (double*)src.x, &K, &beta,
			    (double*)dest.x, &M);
	}
	return 0;
    }

    void dgeMatrix::dgemv(char Tr, double alpha, Rcpp::NumericVector const &X,
			  double beta, double *Y) const {
	int i1 = 1;
	char tr = Trans(Tr);
	bool NTR = tr == 'N';
	if (X.size() != (NTR ? d_ncol : d_nrow))
	    Rf_error("dgemv \"%c\", dim mismatch (%d, %d), X(%d)",
		     tr, d_nrow, d_ncol, X.size());
	F77_CALL(dgemv)(&tr, &d_nrow, &d_ncol, &alpha, d_x.begin(), &d_nrow,
			X.begin(), &i1, &beta, Y, &i1);
    }

    void dgeMatrix::dgemv(char Tr, double alpha,
			  Rcpp::NumericVector const &X, double beta,
			  Rcpp::NumericVector &Y) const {
	int i1 = 1;
	char tr = Trans(Tr);
	bool NTR = tr == 'N';
	if (X.size() != (NTR ? d_ncol : d_nrow) ||
	    Y.size() != (NTR ? d_nrow : d_ncol))
	    Rf_error("dgemv \"%c\", dim mismatch (%d, %d), X(%d), Y(%d)",
		     tr, d_nrow, d_ncol, X.size(), Y.size());
	F77_CALL(dgemv)(&tr, &d_nrow, &d_ncol, &alpha, d_x.begin(),
			&d_nrow, X.begin(), &i1, &beta,
			Y.begin(), &i1);
    }

    void dgeMatrix::dgemm(char TRA, char TRB,
			  double alpha, const dgeMatrix &B,
			  double beta, dgeMatrix &C) const {
	char trA = Trans(TRA), trB = Trans(TRB);
	bool NTA = trA == 'N', NTB = trB == 'N';
	int M = NTA ? d_nrow : d_ncol,
	    N = NTB ? B.ncol() : B.nrow(),
	    K = NTA ? d_ncol : d_nrow;
	int Bnr = B.nrow();
	if (NTB ? B.ncol() : B.nrow() != K || C.nrow() != M || C.ncol() != N)
	    Rf_error("dgemm \"%c,%c\", dim mismatch (%d, %d), (%d,%d), (%d,%d)",
		     trA, trB, d_nrow, d_ncol, B.nrow(), B.ncol(),
		     C.nrow(), C.ncol());
	F77_CALL(dgemm)(&trA, &trB, &M, &N, &K, &alpha, d_x.begin(), &d_nrow,
			B.x().begin(), &Bnr, &beta, C.x().begin(), &M);
    }

    dtrMatrix::dtrMatrix(Rcpp::S4& xp)
	: ddenseMatrix(xp),
	  triangularMatrix(xp) {
    }

    dtrMatrix::dtrMatrix(int nr, char ul, char di)
	: ddenseMatrix(nr, nr),
	  triangularMatrix(ul,di) {
    }

    dsyMatrix::dsyMatrix(Rcpp::S4& xp)
	: ddenseMatrix(xp),
	  symmetricMatrix(xp) {
    }

    dsyMatrix::dsyMatrix(int nr, char ul)
	: ddenseMatrix(nr, nr),
	  symmetricMatrix(ul) {
    }

    void dsyMatrix::dsyrk(dgeMatrix const& A, double alpha, double beta) {
	if (d_nrow != A.ncol())
	    Rf_error("%s dimension mismatch, (%d,%d) vs A(%d,%d)",
		     "dsyMatrix::dsyrk", d_nrow, d_ncol,
		     A.nrow(), A.ncol());
	int Anr = A.nrow();
	F77_CALL(dsyrk)(&d_ul, "T", &d_nrow, &Anr, &alpha,
			A.x().begin(), &Anr, &beta, d_x.begin(), &d_nrow);
    }

    dpoMatrix::dpoMatrix(Rcpp::S4& xp)
	: dsyMatrix(xp) {
    }

    dpoMatrix::dpoMatrix(int nr, char ul)
	: dsyMatrix(nr, ul) {
    }

    Cholesky::Cholesky(Rcpp::S4 xp)
	: dtrMatrix(xp) {
    }

    Cholesky::Cholesky(dgeMatrix A, char ul)
	: dtrMatrix(A.ncol(), ul) {
	update(A);
    }

    void Cholesky::update(dgeMatrix const& A) {
	if (d_nrow != A.ncol())
	    Rf_error("%s dimension mismatch, (%d,%d) vs A(%d,%d)",
		     "Cholesky::update(dgeMatrix)", d_nrow, d_ncol,
		     A.nrow(), A.ncol());
	double alpha = 1., beta = 0.;
	int Anr = A.nrow();
	F77_CALL(dsyrk)(&d_ul, "T", &d_nrow, &Anr, &alpha,
			A.x().begin(), &Anr, &beta, d_x.begin(), &d_nrow);
	int info;
	F77_CALL(dpotrf)(&d_ul, &d_nrow, d_x.begin(), &d_nrow, &info);
	if (info)
	    Rf_error("Lapack routine %s returned error code %d",
		     "dpotrf", info);
    }

    void Cholesky::update(dpoMatrix const &A) {
	if (d_nrow != A.nrow())
	    Rf_error("%s dimension mismatch, (%d,%d) vs A(%d,%d)",
		     "Cholesky::update(dpoMatrix)", d_nrow, d_ncol,
		     A.nrow(), A.ncol());
	d_ul = A.uplo();
	std::copy(A.x().begin(), A.x().end(), d_x.begin());
	int info;
	F77_CALL(dpotrf)(&d_ul, &d_nrow, d_x.begin(), &d_nrow, &info);
	if (info)
	    Rf_error("Lapack routine %s returned error code %d",
		     "dpotrf", info);
    }

    void Cholesky::update(char Tr, double alpha, const dgeMatrix &A,
			  double beta, const dsyMatrix &C) {
	const char tr = Trans(Tr);
	const bool NTR = tr == 'N';
	int Anr = A.nrow(), Anc = A.ncol(), Cnr = C.nrow(), Cnc = C.ncol();
	if (d_nrow != Cnr || NTR ? Anr : Anc != d_nrow)
	    Rf_error("%s(\"%c\") dimension mismatch, (%d,%d), A(%d,%d), C(%d,%d)",
		     "Cholesky::update(dpoMatrix, dgeMatrix)", tr,
		     d_nrow, d_ncol, Anr, Anc, Cnr, Cnc);
	d_ul = C.uplo();
	std::copy(C.x().begin(), C.x().end(), d_x.begin());
	F77_CALL(dsyrk)(&d_ul, &tr, &d_nrow, NTR ? &Anc : &Anr, &alpha,
			A.x().begin(), &Anr, &beta, d_x.begin(), &d_nrow);
	int info;
	F77_CALL(dpotrf)(&d_ul, &d_nrow, d_x.begin(), &d_nrow, &info);
	if (info)
	    Rf_error("Lapack routine %s returned error code %d",
		     "dpotrf", info);
    }

    void Cholesky::dpotrs(double *v, int nb) const {
	int info;
	F77_CALL(dpotrs)(&d_ul, &d_nrow, &nb, d_x.begin(), &d_nrow,
			 v, &d_nrow, &info);
	if (info)
	    Rf_error("Lapack routine %s returned error code %d",
		     "dpotrs", info);
    }

    void Cholesky::dpotrs(Rcpp::NumericVector &v) const {
	if (v.size() != d_nrow)
	    Rf_error("%s (%d, %d) dimension mismatch (%d, %d)",
		     "Cholesky::dpotrs", d_nrow, d_ncol, v.size(), 1);
	dpotrs(v.begin());
    }

    void Cholesky::dpotrs(std::vector<double> &v) const {
	if ((int)v.size() != d_nrow)
	    Rf_error("%s (%d, %d) dimension mismatch (%d, %d)",
		     "Cholesky::dpotrs", d_nrow, d_ncol, v.size(), 1);
	dpotrs(&v[0]);
    }

    Rcpp::NumericMatrix Cholesky::solve(int sys, const_CHM_DN B) const {
	if (sys != CHOLMOD_A) Rf_error("Code not yet written");
	if ((int)B->nrow != d_nrow)
	    Rf_error("%s (%d, %d) dimension mismatch (%d, %d)",
		     "Cholesky::solve",
		     d_nrow, d_ncol, B->nrow, B->ncol);
	Rcpp::NumericMatrix ans(B->nrow, B->ncol);
	double *bx = (double*)B->x;
	std::copy(bx, bx + ans.size(), ans.begin());
	dpotrs(ans.begin(), ans.ncol());
	return ans;
    }
	
    Rcpp::NumericMatrix Cholesky::solve(int                 sys,
					Rcpp::NumericMatrix const&  B) const {
	if (sys != CHOLMOD_A) Rf_error("Code not yet written");
	if (B.nrow() != d_nrow)
	    Rf_error("%s (%d, %d) dimension mismatch (%d, %d)",
		     "Cholesky::solve",
		     d_nrow, d_ncol, B.nrow(), B.ncol());
	Rcpp::NumericMatrix ans = clone(B);
	dpotrs(ans.begin(), ans.ncol());
	return ans;
    }
	
    Rcpp::NumericMatrix Cholesky::solve(int                 sys,
					Rcpp::NumericVector const&  B) const {
	if (sys != CHOLMOD_A) Rf_error("Code not yet written");
	if (B.size() != d_nrow)
	    Rf_error("%s (%d, %d) dimension mismatch (%d, %d)",
		     "Cholesky::solve",
		     d_nrow, d_ncol, B.size(), 1);
	Rcpp::NumericMatrix ans(B.size(), 1);
	std::copy(B.begin(), B.end(), ans.begin());
	dpotrs(ans.begin(), ans.ncol());
	return ans;
    }
	
    double Cholesky::logDet2() const {
	int nc = ncol(), stride = nrow() + 1;
	const double *rx = d_x.begin();
	double ans = 0.;
	for (int i = 0; i < nc; i++, rx += stride)
	    ans += 2. * log(*rx);
	return ans;
    }

    chmFr::chmFr(Rcpp::S4 xp)
	: d_xp(xp)
    {
	Rcpp::CharacterVector cl(SEXP(xp.attr("class")));
	char *clnm = cl[0];
	if (!xp.is("CHMfactor"))
	    ::Rf_error("Class %s object passed to %s is not a %s",
		       clnm, "chmFr::chmFr", "CHMfactor");
	Rcpp::Dimension Dim(SEXP(xp.slot("Dim")));
	Rcpp::IntegerVector colcount(SEXP(xp.slot("colcount"))),
	    perm(SEXP(xp.slot("perm"))),
	    type(SEXP(xp.slot("type")));
	Rcpp::NumericVector X(SEXP(xp.slot("x")));

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
	    Rcpp::IntegerVector
		Pi(SEXP(xp.slot("pi"))),
		SUPER(SEXP(xp.slot("super"))),
		S(SEXP(xp.slot("s")));
	    Rcpp::NumericVector PX(SEXP(xp.slot("px")));

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
	    Rcpp::IntegerVector
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
	M_cholmod_factorize_p((const_CHM_SP)&A, &Imult,
			      (int*)NULL, (size_t) 0, this, &c);
    }

    double chmFr::logDet2() const {   // Need Matrix_0.999375-42 or later
//    double chmFr::logDet2() {
	return M_chm_factor_ldetL2(this);
    }

    Rcpp::NumericMatrix chmFr::solve(int sys, const_CHM_DN b) const {
	CHM_DN t1 = M_cholmod_solve(sys, (const_CHM_FR)this, b, &c);
	Rcpp::NumericMatrix ans((int) t1->nrow, (int) t1->ncol);
	double *tx = (double*)t1->x;
	std::copy(tx, tx + ans.size(), ans.begin());
	M_cholmod_free_dense(&t1, &c);
	return ans;
    }

    Rcpp::NumericMatrix chmFr::solve(int sys, Rcpp::NumericMatrix const& b) const {
	const chmDn cb(b);
	return solve(sys, &cb);
    }

    Rcpp::NumericMatrix chmFr::solve(int sys, Rcpp::NumericVector const& b) const {
	const chmDn cb(b);
	return solve(sys, &cb);
    }
    
    CHM_SP chmFr::spsolve(int sys, const_CHM_SP b) const {
	return M_cholmod_spsolve(sys, (const CHM_FR)this, b, &c);
    }
    
    CHM_SP chmFr::spsolve(int sys, chmSp const &b) const {
	return M_cholmod_spsolve(sys, (const CHM_FR)this,
				 (const_CHM_SP)&b, &c);
    }
    
    chmSp::chmSp(Rcpp::S4 xp)
	: d_xp(xp)
    {
	if (!xp.is("CsparseMatrix")) {
	    Rcpp::CharacterVector cls = SEXP(xp.attr("class"));
	    char *clnm = cls[0];
	    Rf_error("Class %s object passed to %s is not a %s",
		     clnm, "chmSp::chmSp", "CsparseMatrix");
	}
	Rcpp::IntegerVector
	    Dim(xp.slot("Dim")), pp(xp.slot("p")), ii(xp.slot("i"));
	nrow = Dim[0];
	ncol = Dim[1];
	p = (void*)pp.begin();
	i = (void*)ii.begin();
	dtype = 0;  // CHOLMOD_DOUBLE
	stype = 0;
	if (!xp.is("generalMatrix")) {
	    if (xp.is("symmetricMatrix")) {
		Rcpp::CharacterVector uplo(SEXP(xp.slot("uplo")));
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
	    Rcpp::NumericVector xx(SEXP(xp.slot("x")));
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
	return M_cholmod_transpose((const_CHM_SP)this, values, &c);
    }

    int chmSp::dmult(char tr, double alpha, double beta,
		     chmDn const &src, chmDn &dest) const {
	return M_cholmod_sdmult((const_CHM_SP)this,
				Trans(tr) == 'T', &alpha,
				&beta, &src, &dest, &c);
    }

    void chmSp::update(cholmod_sparse const &nn) {
	size_t nnznn = M_cholmod_nnz((const_CHM_SP)&nn, &c);
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

    CHM_SP chmSp::crossprod(const_CHM_SP B, int sorted) const {
	CHM_SP t1 = this->transpose();
	CHM_SP t2 = ::M_cholmod_ssmult(t1, B, 0/*stype*/, xtype/*values*/,
				       sorted, &c);
	::M_cholmod_free_sparse(&t1, &c);
	return t2;
    }

    CHM_SP chmSp::crossprod(chmSp const &B, int sorted) const {
	return crossprod((const_CHM_SP)&B, sorted);
    }

    CHM_SP chmSp::tcrossprod() const {
	CHM_SP t1 = M_cholmod_aat((const_CHM_SP)this, (int*)NULL,
				  0/*fsize*/, xtype/*mode*/, &c);
	CHM_SP t2 = ::M_cholmod_copy(t1, 1/*stype*/, xtype/*mode*/, &c);
	::M_cholmod_free_sparse(&t1, &c);
	return t2;
    }

    CHM_SP chmSp::tcrossprod(const_CHM_SP B, int sorted) const {
	CHM_SP t1 = M_cholmod_transpose(B, xtype/*values*/, &c);
	CHM_SP t2 =
	    M_cholmod_ssmult((const_CHM_SP)this, t1, 0/*stype*/,
			     xtype/*values*/, sorted, &c);
	M_cholmod_free_sparse(&t1, &c);
	return t2;
    }

    CHM_SP chmSp::tcrossprod(chmSp const &B, int sorted) const {
	return tcrossprod((const_CHM_SP)&B, sorted);
    }

    CHM_SP chmSp::smult(chmSp const &B, int stype, int values,
			int sorted) const {
	return M_cholmod_ssmult((const_CHM_SP)this, (const_CHM_SP)&B,
				stype, values, sorted, &c);
    }

    void chmSp::scale(int whch, chmDn const& sc) {
	M_cholmod_scale(&sc, whch, this, &c);
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

    void chmDn::init(const double *X, int r, int c) {
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

    chmDn::chmDn(double *xx, int nr, int nc)
	: cholmod_dense() {
	this->init(xx, nr, nc);
    }

    chmDn::chmDn(const double *xx, int nr, int nc)
	: cholmod_dense() {
	    this->init(xx, nr, nc);
    }

    chmDn::chmDn(std::vector<double> &v)
	: cholmod_dense() {
	this->init(&v[0], v.size(), 1);
    }

    chmDn::chmDn(std::vector<double> const &v)
	: cholmod_dense() {
	this->init(&v[0], v.size(), 1);
    }

    chmDn::chmDn(Rcpp::NumericVector &v)
	: cholmod_dense() {
	this->init(v.begin(), v.size(), 1);
    }

    chmDn::chmDn(Rcpp::NumericVector const &v)
	: cholmod_dense() {
	this->init(v.begin(), v.size(), 1);
    }

    chmDn::chmDn(Rcpp::NumericMatrix &m)
	: cholmod_dense() {
	this->init(m.begin(), m.nrow(), m.ncol());
    }

    chmDn::chmDn(Rcpp::NumericMatrix const &m)
	: cholmod_dense() {
	this->init(m.begin(), m.nrow(), m.ncol());
    }

    chmDn::chmDn(ddenseMatrix &m)
	: cholmod_dense() {
	this->init(m.x().begin(), m.nrow(), m.ncol());
    }
    
    chmDn::chmDn(ddenseMatrix const &m)
	: cholmod_dense() {
	this->init(m.x().begin(), m.nrow(), m.ncol());
    }
}
