#include "MatrixNs.h"
#include <R_ext/Lapack.h>
#include <cctype>

using namespace Rcpp;
using namespace std;

namespace MatrixNs{
    char chkchar(char x, std::string allowed) {
	char X = toupper(x);
	if (find(allowed.begin(), allowed.end(), X) == allowed.end())
	    throw range_error("chkchar");
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

    Matrix::Matrix(Rcpp::S4 &xp)
	: d_sexp(SEXP(xp)),
	  d_dimnames(xp.slot("Dimnames")) {
	IntegerVector Dim(xp.slot("Dim"));
	d_nrow = Dim[0];
	d_ncol = Dim[1];
    }

    Matrix::Matrix(int nr, int nc)
	: d_sexp(R_NilValue),
	  d_dimnames(2),
	  d_nrow(nr),
	  d_ncol(nc) {
	if (nr < 0 || nc < 0)
	    throw range_error("Matrix(nr,nc)");
    }
	
    int Matrix::nrow() const { return d_nrow; }
    int Matrix::ncol() const { return d_ncol; }

    modelMatrix::modelMatrix(Rcpp::S4 &xp)
	: d_assign(   xp.slot("assign")),
	  d_contrasts(xp.slot("contrasts")) {
    }


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
	    throw range_error("Size mismatch in setX");
	copy(nx.begin(), nx.end(), d_x.begin());
    }

    void dMatrix::setX(Rcpp::NumericMatrix const& mm) {
	setX(NumericVector(SEXP(mm)));
    }

// Check this: Do the dspMatrix and dtpMatrix classes pass this check?
    ddenseMatrix::ddenseMatrix(Rcpp::S4 &xp) throw (std::runtime_error) : dMatrix(xp) {
	if (!d_x.size() == d_nrow * d_ncol)
	    throw std::runtime_error("dimension mismatch");
	    // ::Rf_error("%s: Dim = (%d, %d) is inconsistent with x.size() = %d",
	    // 	       "ddenseMatrix::ddenseMatrix", d_nrow, d_ncol, d_x.size());
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
			  double beta, double *Y) const throw (std::runtime_error) {
	int i1 = 1;
	char tr = Trans(Tr);
	bool NTR = tr == 'N';
	if (X.size() != (NTR ? d_ncol : d_nrow))
	    throw std::runtime_error("dimension mismatch");
	    // Rf_error("dgemv \"%c\", dim mismatch (%d, %d), X(%d)",
	    // 	     tr, d_nrow, d_ncol, X.size());
	F77_CALL(dgemv)(&tr, &d_nrow, &d_ncol, &alpha, d_x.begin(), &d_nrow,
			X.begin(), &i1, &beta, Y, &i1);
    }

    void dgeMatrix::dgemv(char Tr, double alpha,
			  Rcpp::NumericVector const &X, double beta,
			  Rcpp::NumericVector &Y) const throw (std::runtime_error) {
	int i1 = 1;
	char tr = Trans(Tr);
	bool NTR = tr == 'N';
	if (X.size() != (NTR ? d_ncol : d_nrow) ||
	    Y.size() != (NTR ? d_nrow : d_ncol))
	    throw std::runtime_error("dimension mismatch");
	    // Rf_error("dgemv \"%c\", dim mismatch (%d, %d), X(%d), Y(%d)",
	    // 	     tr, d_nrow, d_ncol, X.size(), Y.size());
	F77_CALL(dgemv)(&tr, &d_nrow, &d_ncol, &alpha, d_x.begin(),
			&d_nrow, X.begin(), &i1, &beta,
			Y.begin(), &i1);
    }

    void dgeMatrix::dgemm(char TRA, char TRB,
			  double alpha, const dgeMatrix &B,
			  double beta, dgeMatrix &C) const throw (std::runtime_error) {
	char trA = Trans(TRA), trB = Trans(TRB);
	bool NTA = trA == 'N', NTB = trB == 'N';
	int M = NTA ? d_nrow : d_ncol,
	    N = NTB ? B.ncol() : B.nrow(),
	    K = NTA ? d_ncol : d_nrow;
	int Bnr = B.nrow();
	if (NTB ? B.ncol() : B.nrow() != K || C.nrow() != M || C.ncol() != N)
	    throw std::runtime_error("dimension mismatch");
	    // Rf_error("dgemm \"%c,%c\", dim mismatch (%d, %d), (%d,%d), (%d,%d)",
	    // 	     trA, trB, d_nrow, d_ncol, B.nrow(), B.ncol(),
	    // 	     C.nrow(), C.ncol());
	F77_CALL(dgemm)(&trA, &trB, &M, &N, &K, &alpha, d_x.begin(), &d_nrow,
			B.x().begin(), &Bnr, &beta, C.x().begin(), &M);
    }

    dgeMatrix::operator SEXP() const {
	S4 ans("dgeMatrix");
	ans.slot("x") = clone(d_x);
	ans.slot("Dimnames") = clone(d_dimnames);
	ans.slot("Dim") = IntegerVector::create(d_nrow, d_ncol);
	return ans;
    }

    ddenseModelMatrix::ddenseModelMatrix(Rcpp::S4 xp)
	: dgeMatrix(  xp),
	  modelMatrix(xp) {
    }

    dtrMatrix::dtrMatrix(Rcpp::S4& xp)
	: ddenseMatrix(xp),
	  triangularMatrix(xp) {
    }

    dtrMatrix::dtrMatrix(int nr, char ul, char di)
	: ddenseMatrix(nr, nr),
	  triangularMatrix(ul,di) {
    }

    void dtrMatrix::dtrtrs(char tr, double* v, int nb) const throw (std::runtime_error) {
	int info;
	F77_CALL(dtrtrs)(&d_ul, &tr, &d_di, &d_nrow, &nb, d_x.begin(),
			 &d_nrow, v, &d_nrow, &info);
	if (info) throw std::runtime_error("Lapack routine dtrtrs returned an error code");
	    // Rf_error("Lapack routine %s returned error code %d",
	    // 	     "dtrtrs", info);
    }

    dtrMatrix::operator SEXP() const {
	S4 ans("dtrMatrix");
	ans.slot("x") = clone(d_x);
	ans.slot("Dimnames") = clone(d_dimnames);
	ans.slot("Dim") = IntegerVector::create(d_nrow, d_ncol);
	ans.slot("uplo") = uplo() == 'U' ? "U" : "L";	
	ans.slot("diag") = diag() == 'U' ? "U" : "N";
	return ans;
    }


    dsyMatrix::dsyMatrix(Rcpp::S4& xp)
	: ddenseMatrix(xp),
	  symmetricMatrix(xp) {
    }

    dsyMatrix::dsyMatrix(int nr, char ul)
	: ddenseMatrix(nr, nr),
	  symmetricMatrix(ul) {
    }

    void dsyMatrix::dsyrk(dgeMatrix const& A, double alpha, double beta) throw (std::runtime_error) {
	if (d_nrow != A.ncol())
	    throw std::runtime_error("dimension mismatch");
	    // Rf_error("%s dimension mismatch, (%d,%d) vs A(%d,%d)",
	    // 	     "dsyMatrix::dsyrk", d_nrow, d_ncol,
	    // 	     A.nrow(), A.ncol());
	int Anr = A.nrow();
	F77_CALL(dsyrk)(&d_ul, "T", &d_nrow, &Anr, &alpha,
			A.x().begin(), &Anr, &beta, d_x.begin(), &d_nrow);
    }

    dsyMatrix::operator SEXP() const {
	S4 ans("dsyMatrix");
	ans.slot("x") = clone(d_x);
	ans.slot("Dimnames") = clone(d_dimnames);
	ans.slot("Dim") = IntegerVector::create(d_nrow, d_ncol);
	ans.slot("uplo") = uplo() == 'U' ? "U" : "L";
	return ans;
    }

    dpoMatrix::dpoMatrix(Rcpp::S4& xp)
	: dsyMatrix(xp) {
    }

    dpoMatrix::dpoMatrix(int nr, char ul)
	: dsyMatrix(nr, ul) {
    }

    dpoMatrix::operator SEXP() const {
	S4 ans("dpoMatrix");
	ans.slot("x") = clone(d_x);
	ans.slot("Dimnames") = clone(d_dimnames);
	ans.slot("Dim") = IntegerVector::create(d_nrow, d_ncol);
	ans.slot("uplo") = uplo() == 'U' ? "U" : "L";
	return ans;
    }

    Cholesky::Cholesky(Rcpp::S4 xp)
	: dtrMatrix(xp) {
    }

    Cholesky::Cholesky(dgeMatrix A, char ul)
	: dtrMatrix(A.ncol(), ul) {
	update(A);
    }

    void Cholesky::update(dgeMatrix const& A) throw (std::runtime_error) {
	if (d_nrow != A.ncol())
	    throw std::runtime_error("dimension mismatch");
	    // Rf_error("%s dimension mismatch, (%d,%d) vs A(%d,%d)",
	    // 	     "Cholesky::update(dgeMatrix)", d_nrow, d_ncol,
	    // 	     A.nrow(), A.ncol());
	double alpha = 1., beta = 0.;
	int Anr = A.nrow();
	F77_CALL(dsyrk)(&d_ul, "T", &d_nrow, &Anr, &alpha,
			A.x().begin(), &Anr, &beta, d_x.begin(), &d_nrow);
	int info;
	F77_CALL(dpotrf)(&d_ul, &d_nrow, d_x.begin(), &d_nrow, &info);
	if (info)
	    throw std::runtime_error("Lapack routine dpotrf returned error code");
	    // Rf_error("Lapack routine %s returned error code %d",
	    // 	     "dpotrf", info);
    }

    void Cholesky::update(dpoMatrix const &A) throw (std::runtime_error) {
	if (d_nrow != A.nrow())
	    throw std::runtime_error("dimension mismatch");
	    // Rf_error("%s dimension mismatch, (%d,%d) vs A(%d,%d)",
	    // 	     "Cholesky::update(dpoMatrix)", d_nrow, d_ncol,
	    // 	     A.nrow(), A.ncol());
	d_ul = A.uplo();
	copy(A.x().begin(), A.x().end(), d_x.begin());
	int info;
	F77_CALL(dpotrf)(&d_ul, &d_nrow, d_x.begin(), &d_nrow, &info);
	if (info)
	    throw std::runtime_error("Lapack routine dpotrf returned error code");
	    // Rf_error("Lapack routine %s returned error code %d",
	    // 	     "dpotrf", info);
    }

    void Cholesky::update(char Tr, double alpha, const dgeMatrix &A,
			  double beta, const dsyMatrix &C) throw (std::runtime_error) {
	const char tr = Trans(Tr);
	const bool NTR = tr == 'N';
	int Anr = A.nrow(), Anc = A.ncol(), Cnr = C.nrow(); //, Cnc = C.ncol();
	if (d_nrow != Cnr || NTR ? Anr : Anc != d_nrow)
	    throw std::runtime_error("dimension mismatch");
	    // Rf_error("%s(\"%c\") dimension mismatch, (%d,%d), A(%d,%d), C(%d,%d)",
	    // 	     "Cholesky::update(dpoMatrix, dgeMatrix)", tr,
	    // 	     d_nrow, d_ncol, Anr, Anc, Cnr, Cnc);
	d_ul = C.uplo();
	copy(C.x().begin(), C.x().end(), d_x.begin());
	F77_CALL(dsyrk)(&d_ul, &tr, &d_nrow, NTR ? &Anc : &Anr, &alpha,
			A.x().begin(), &Anr, &beta, d_x.begin(), &d_nrow);
	int info;
	F77_CALL(dpotrf)(&d_ul, &d_nrow, d_x.begin(), &d_nrow, &info);
	if (info)
	    throw std::runtime_error("Lapack routine dpotrf returned error code");
	    // Rf_error("Lapack routine %s returned error code %d",
	    // 	     "dpotrf", info);
    }

    void Cholesky::dpotrs(double *v, int nb) const throw (std::runtime_error) {
	int info;
	F77_CALL(dpotrs)(&d_ul, &d_nrow, &nb, d_x.begin(), &d_nrow,
			 v, &d_nrow, &info);
	if (info)
	    throw std::runtime_error("Lapack routine dpotrs returned error code");
	    // Rf_error("Lapack routine %s returned error code %d",
	    // 	     "dpotrs", info);
    }

    void Cholesky::dpotrs(Rcpp::NumericVector &v) const throw (std::runtime_error) {
	if (v.size() != d_nrow)
	    throw std::runtime_error("dimension mismatch");
	    // Rf_error("%s (%d, %d) dimension mismatch (%d, %d)",
	    // 	     "Cholesky::dpotrs", d_nrow, d_ncol, v.size(), 1);
	dpotrs(v.begin());
    }

    void Cholesky::dpotrs(std::vector<double> &v) const throw (std::runtime_error) {
	if ((int)v.size() != d_nrow)
	    throw std::runtime_error("dimension mismatch");
	    // Rf_error("%s (%d, %d) dimension mismatch (%d, %d)",
	    // 	     "Cholesky::dpotrs", d_nrow, d_ncol, v.size(), 1);
	dpotrs(&v[0]);
    }

    void Cholesky::inPlaceSolve(int sys, Rcpp::NumericMatrix& B) const throw (std::runtime_error) {
	if (B.nrow() != d_nrow)
	    throw std::runtime_error("dimension mismatch");
	    // Rf_error("%s (%d, %d) dimension mismatch (%d, %d)",
	    // 	     "Cholesky::solve",
	    // 	     d_nrow, d_ncol, B.nrow(), B.ncol());
	if (sys == CHOLMOD_A) return dpotrs(B.begin(), B.ncol());
	if (sys == CHOLMOD_L)
	    return dtrtrs((d_ul == 'U') ? 'T' : 'N', B.begin(), B.ncol());
	if (sys == CHOLMOD_Lt)
	    return dtrtrs((d_ul == 'U') ? 'N' : 'T', B.begin(), B.ncol());
	throw std::runtime_error("Unknown sys argument for Cholesky::inPlaceSolve");
    }

    Rcpp::NumericMatrix Cholesky::solve(int sys, const_CHM_DN B) const {
	NumericMatrix ans(B->nrow, B->ncol);
	double *bx = (double*)B->x;
	int sz = B->nrow * B->ncol;
	copy(bx, bx + sz, ans.begin());
	inPlaceSolve(sys, ans);
	return ans;
    }
	
    Rcpp::NumericMatrix Cholesky::solve(int                      sys,
					const Rcpp::NumericMatrix& B) const {
	NumericMatrix ans = clone(B);
	inPlaceSolve(sys, ans);
	return ans;
    }
	
    Rcpp::NumericMatrix Cholesky::solve(int                      sys,
					const Rcpp::NumericVector& B) const {
	NumericMatrix ans(B.size(), 1);
	copy(B.begin(), B.end(), ans.begin());
	inPlaceSolve(sys, ans);
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

    Cholesky::operator SEXP() const {
	S4 ans("Cholesky");
// FIXME: Don't use cut-and-paste programming here.
	ans.slot("x") = clone(d_x);
	ans.slot("Dimnames") = clone(d_dimnames);
	ans.slot("Dim") = IntegerVector::create(d_nrow, d_ncol);
	ans.slot("uplo") = uplo() == 'U' ? "U" : "L";	
	ans.slot("diag") = diag() == 'U' ? "U" : "N";
	return ans;
    }

    chmFr::chmFr(Rcpp::S4 xp)  throw (wrongS4) 
	: d_xp(xp) {
	const char* clz = as<const char*>(xp.attr("class"));
	if (!xp.is("CHMfactor"))
	    throw wrongS4(clz, "chmFr::chmFr", "CHMfactor");
	Dimension Dim(SEXP(xp.slot("Dim")));
	IntegerVector colcount(SEXP(xp.slot("colcount"))),
	    perm(SEXP(xp.slot("perm"))),
	    type(SEXP(xp.slot("type")));
	NumericVector X(SEXP(xp.slot("x")));

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
	if (is_super) {
	    IntegerVector
		Pi(SEXP(xp.slot("pi"))),
		SUPER(SEXP(xp.slot("super"))),
		S(SEXP(xp.slot("s")));
	    NumericVector PX(SEXP(xp.slot("px")));

	    if (!xp.is("dCHMsuper"))
		throw wrongS4(clz, "is_super = TRUE", "dCHMsuper");
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
		throw wrongS4(clz, "is_super = FALSE", "dCHMsimpl");
	    nzmax = X.size();
	    p = (void*)P.begin();
	    i = (void*)I.begin();
	    x = (void*)X.begin();
	    nz = (void*)NZ.begin();
	    next = (void*)NXT.begin();
	    prev = (void*)PRV.begin();
	}
    }

    chmFr::operator SEXP() const throw (std::runtime_error) {
// Should this clone d_xp before wrapping?
	if (!d_xp) return wrap(d_xp);

	std::string clz = (is_super) ? "dCHMsuper" : "dCHMsimpl";
	S4 ans(clz);
	ans.slot("Dim") = IntegerVector::create(n, n);
	ans.slot("perm") = IntegerVector((int*)Perm, (int*)Perm + n);
	ans.slot("colcount") = IntegerVector((int*)ColCount, (int*)ColCount + n);
	int ntype = is_super ? 6 : 4;
	IntegerVector type(ntype);
	type[0] = ordering;
	type[1] = (is_ll ? 1 : 0);
	type[2] = is_super;
	type[3] = (is_monotonic ? 1 : 0);
	if (is_super) {		// fill in later
	} else {
	    int *ppt = (int*)p;
	    int nx = ppt[n];
	    ans.slot("x") = NumericVector((double*)x, (double*)x + nx);
	    ans.slot("p") = IntegerVector(ppt, ppt + n + 1);
	    ans.slot("i") = IntegerVector((int*)i, (int*)i + nx);
	    ans.slot("nz") = IntegerVector((int*)nz, (int*)nz + n);
	    ans.slot("nxt") = IntegerVector((int*)next, (int*)next + n + 2);
	    ans.slot("prv") = IntegerVector((int*)prev, (int*)prev + n + 2);
	}
	return ans;
    }

    void chmFr::update(cholmod_sparse const &A, double Imult) {
	M_cholmod_factorize_p((const_CHM_SP)&A, &Imult,
			      (int*)NULL, (size_t) 0, this, &c);
    }

    double chmFr::logDet2() const {   // Needs Matrix_0.999375-41 or later
	return M_chm_factor_ldetL2(this);
    }

    Rcpp::NumericMatrix chmFr::solve(int sys, const_CHM_DN b) const {
	CHM_DN t1 = M_cholmod_solve(sys, (const_CHM_FR)this, b, &c);
	NumericMatrix ans((int) t1->nrow, (int) t1->ncol);
	int sz = t1->nrow * t1->ncol;
	double *tx = (double*)t1->x;
	copy(tx, tx + sz, ans.begin());
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

    
    chmSp::chmSp(Rcpp::S4 xp) throw (wrongS4) 
	: d_xp(xp) {
	if (!xp.is("CsparseMatrix"))
	    throw wrongS4(as<const char*>(xp.slot("class")), "chmSp::chmSp",
			  "CsparseMatrix");
	IntegerVector Dim(xp.slot("Dim")), pp(xp.slot("p")), ii(xp.slot("i"));
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
	if (xp.is("zsparseMatrix")) 
	{
	    xtype = CHOLMOD_COMPLEX;
	    throw std::runtime_error("Not yet defined zsparseMatrix?");
	}
	if (xtype == -1) throw std::runtime_error("Unknown (logical?) sparse Matrix type");
    }

    chmSp::operator SEXP() const throw (std::runtime_error) {
	S4 ans;
	if (xtype != CHOLMOD_REAL && xtype != CHOLMOD_PATTERN)
	    throw std::runtime_error("chmSp object must have xtype of REAL or PATTERN");
	if (xtype == CHOLMOD_REAL) {
	    if (stype) ans = S4("dsCMatrix");
	    else ans = S4("dgCMatrix");
	} else {
	    if (stype) ans = S4("nsCMatrix");
	    else ans = S4("ngCMatrix");
	}
	ans.slot("Dim") = IntegerVector::create(nrow, ncol);
	int *ppt = (int*)p, *ipt = (int*) i;
	int nx = ppt[ncol];
	IntegerVector pp(ppt, ppt + ncol + 1), ii(ipt, ipt + nx);
	ans.slot("p") = pp;
	ans.slot("i") = ii;
	if (xtype == CHOLMOD_REAL) {
	    double *xpt = (double*) x;
	    NumericVector xx(xpt, xpt + nx);
	    ans.slot("x") = xx;
	}
	if (stype) {
	    if (stype == 1) ans.slot("uplo") = "U";
	    else if (stype == -1) ans.slot("uplo") = "L";
	    else throw std::runtime_error("unknown stype");
	}
	return ans;
    }

    CHM_SP chmSp::transpose(int values) const {
	return M_cholmod_transpose((const_CHM_SP)this, values, &c);
    }

    int chmSp::dmult(char tr, double alpha, double beta,
		     chmDn const &src, chmDn &dest) const {
	return M_cholmod_sdmult((const_CHM_SP)this,
				Trans(tr) == 'T', &alpha,
				&beta, &src, &dest, &c);
    }

    void chmSp::update(cholmod_sparse const &nn) throw (std::runtime_error) {
	size_t nnznn = M_cholmod_nnz((const_CHM_SP)&nn, &c);
	if (nn.ncol != ncol || nnznn > nzmax || xtype != nn.xtype ||
	    itype != nn.itype || dtype != nn.dtype || packed != nn.packed)
	    throw std::runtime_error("chmSp::update: matrices not conformable");
//	    Rf_error("%s: matrices not conformable", "chmSp::update");
	stype = nn.stype;
	nrow = nn.nrow;
	sorted = nn.sorted;
	int *pp = (int*)p, *nnp = (int*)(nn.p);
	copy(nnp, nnp + ncol + 1, pp);
	int *ip = (int*)i, *nni = (int*)(nn.i);
	copy(nni, nni + nnznn, ip);
	switch(xtype) {
	case CHOLMOD_PATTERN:
	    break;
	case CHOLMOD_REAL: {
	    double *xp = (double*)x, *nnx = (double*)(nn.x);
	    copy(nnx, nnx + nnznn, xp);
	    break;
	}
	default:
	    throw std::runtime_error("Unsupported CHOLMOD xtype in chmSp::update");
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

namespace Rcpp {
    template <> SEXP
    wrap<MatrixNs::dgeMatrix>(const MatrixNs::dgeMatrix& m) {return m;}
    template <> SEXP
    wrap<MatrixNs::dtrMatrix>(const MatrixNs::dtrMatrix& m) {return m;}
    template <> SEXP
    wrap<MatrixNs::dsyMatrix>(const MatrixNs::dsyMatrix& m) {return m;}
    template <> SEXP
    wrap<MatrixNs::dpoMatrix>(const MatrixNs::dpoMatrix& m) {return m;}
    template <> SEXP
    wrap<MatrixNs::Cholesky>(const MatrixNs::Cholesky& m) {return m;}
    template <> SEXP
    wrap<MatrixNs::chmFr>(const MatrixNs::chmFr& m) {return m;}
    template <> SEXP
    wrap<MatrixNs::chmSp>(const MatrixNs::chmSp& m) {return m;}
}
