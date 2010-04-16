#include "MatrixNs.h"
#include <cctype>		// for toupper
#include <R_ext/Lapack.h>

using namespace Rcpp;

extern cholmod_common c;

namespace MatrixNs{

// A function like this should be part of the Rcpp::S4 class
    static bool isClass(S4 x, std::string clname) {
	CharacterVector cl = x.attr("class");
	SEXP pkg = cl.attr("package");
	if (as<std::string>(cl) == clname) return true;
	Function clDef("getClassDef");
	S4 cld = clDef(cl);
	List exts = cld.attr("contains");
	Function ssc(".selectSuperClasses");
	CharacterVector nms = ssc(exts, false, true, false);
	for (int i = 0; i < nms.size(); i++)
	    if (std::string(nms[i]) == clname) return true;
	return false;
    }

    Matrix::Matrix(S4 &xp) :
	Dim(SEXP(xp.slot("Dim"))),
	Dimnames(SEXP(xp.slot("Dimnames"))) {
    }

// These should be templated
    dMatrix::dMatrix(S4 &xp) :
	Matrix(xp),
	x(SEXP(xp.slot("x"))) {
    }

// Check this: Do the dspMatrix and dtpMatrix classes pass this check?
    ddenseMatrix::ddenseMatrix(S4 &xp) : dMatrix(xp) {
	if (!x.size() == Dim[0] * Dim[1])
	    ::Rf_error("%s: Dim = (%d, %d) is inconsistent with x.size() = %d",
		       "ddenseMatrix::ddenseMatrix", Dim[0], Dim[1], x.size());
    }
    
    compMatrix::compMatrix(S4 &xp) : factors(SEXP(xp.slot("factors"))) {
    }

    generalMatrix::generalMatrix(S4 &xp) : compMatrix(xp) {
    }

    triangularMatrix::triangularMatrix(S4 &xp) :
	uplo(SEXP(xp.slot("uplo"))),
	diag(SEXP(xp.slot("diag"))) {
    }

    symmetricMatrix::symmetricMatrix(S4 &xp) :
	compMatrix(xp),
	uplo(SEXP(xp.slot("uplo"))) {
    }

    static inline char trCan(char tt) { // canonical trans character
	char ans = toupper(tt);
	if (ans != 'C' || ans != 'N' || ans != 'T')
	    Rf_error("invalid trans specification %c", ans);
	return ans;
    }

    void dgeMatrix::dgemv(const char* trans,
			  double alpha, const NumericVector &X,
			  double beta, NumericVector &Y) {
	int i1 = 1;
	char tr = trCan(*trans);
	bool NTR = tr == 'N';
	if (X.size() != Dim[NTR ? 1 : 0] || Y.size() != Dim[NTR ? 0 : 1])
	    Rf_error("dgemv \"%c\", dim mismatch (%d, %d), X(%d), Y(%d)",
		     tr, Dim[0], Dim[1], X.size(), Y.size());
	F77_CALL(dgemv)(&tr, &Dim[0], &Dim[1], &alpha, x.begin(), &Dim[0],
			X.begin(), &i1, &beta, Y.begin(), &i1);
    }

    void dgeMatrix::dgemm(const char* transA, const char* transB,
			  double alpha, const dgeMatrix &B,
			  double beta, dgeMatrix &C) {
	int i1 = 1;
	char trA = trCan(*transA), trB = trCan(*transB);
	bool NTA = trA == 'N', NTB = trB == 'N';
	Dimension Bd(B.Dim), Cd(C.Dim);
	int M = Dim[NTA ? 0 : 1], N = Bd[NTB ? 1 : 0], K = Dim[NTA ? 1 : 0];
	if (Bd[NTB ? 0 : 1] != K || Cd[0] != M || Cd[1] != N)
	    Rf_error("dgemm \"%c,%c\", dim mismatch (%d, %d), (%d,%d), (%d,%d)",
		     trA, trB, Dim[0], Dim[1], Bd[0], Bd[1], Cd[0], Cd[1]);
	F77_CALL(dgemm)(&trA, &trB, &M, &N, &K, &alpha, x.begin(), &Dim[0],
			B.x.begin(), &Bd[0], &beta, C.x.begin(), &Cd[0]);
    }

				//! This class is just a wrapper for a CHM_FR struct
    dCHMfactor::dCHMfactor(S4 xp) :
	Dim(SEXP(xp.slot("Dim"))),
	ColCount(SEXP(xp.slot("colcount"))),
	perm(SEXP(xp.slot("perm"))),
	type(SEXP(xp.slot("type"))),
	x(SEXP(xp.slot("x")))
    {
	CharacterVector cls = xp.attr("class");
	std::string cl(cls[0]);
	fa = new cholmod_factor;
	fa->minor = fa->n = (size_t)*(Dim.begin());
	fa->Perm = (void*)perm.begin();
	fa->ColCount = (void*)ColCount.begin();
	fa->x = (void*)x.begin();
	fa->ordering = type[0];
	fa->is_ll = (type[1] ? 1 : 0);
	fa->is_super = 0;
	fa->is_monotonic = (type[3] ? 1 : 0);
	fa->xtype = CHOLMOD_REAL;
	fa->itype = CHOLMOD_LONG;
	fa->dtype = 0;  // CHOLMOD_DOUBLE
	fa->z = (void*)NULL;
	if (cl == "dCHMsimpl") { //FIXME: check is(xp, "dCHMsimpl") instead
	    IntegerVector
		colcount(SEXP(xp.slot("colcount"))),
		i(SEXP(xp.slot("i"))),
		nxt(SEXP(xp.slot("nxt"))),
		nz(SEXP(xp.slot("nz"))),
		p(SEXP(xp.slot("p"))),
		prv(SEXP(xp.slot("prv")));
	    NumericVector x(SEXP(xp.slot("x")));
	    
	    if (fa->is_super)
		::Rf_error("conflict between type[2] = %d and class %s",
			 type[2], cl.c_str());
	    fa->nzmax = (size_t) x.size();
	    fa->p = (void*) p.begin();
	    fa->i = (void*) i.begin();
	    fa->x = (void*) x.begin();
	    fa->nz = (void*) nz.begin();
	    fa->next = (void*) nxt.begin();
	    fa->prev = (void*) prv.begin();
	    return;
	}
	if (cl == "dCHMsuper") { //FIXME: check is(xp, "dCHMsuper") instead
	    IntegerVector
		pi(SEXP(xp.slot("pi"))),
		super(SEXP(xp.slot("super"))),
		s(SEXP(xp.slot("s")));
	    NumericVector px(SEXP(xp.slot("px")));

	    if (!(fa->is_super))
		::Rf_error("conflict between type[2] = %d and class %s",
			 type[2], cl.c_str());
	    fa->xsize = (size_t) x.size();
	    fa->ssize = (size_t) s.size();
	    fa->maxcsize = (size_t) type[4];
	    fa->maxesize = (size_t) type[5];
	    fa->nsuper = super.size() - 1;
	    fa->super = (void*)super.begin();
	    fa->pi = (void*)pi.begin();
	    fa->s = (void*)s.begin();
	    fa->px = (void*)px.begin();
	    return;
	}
	::Rf_error("Class %s is not suitable for dCHMfactor", cl.c_str());
    }

    void dCHMfactor::update(CHM_SP A, double Imult) {
	::M_cholmod_factorize_p(A, &Imult, (int*)NULL, (size_t)0, fa, &c);
    }
    
    chmSp::chmSp(S4 xp) : cholmod_sparse() {
	CharacterVector cl(SEXP(xp.attr("class")));
	char *clnm = cl[0];
	if (!isClass(xp, "CsparseMatrix"))
	    Rf_error("Class %s object passed to %s is not a %s",
		     clnm, "chmSp::chmSp", "CsparseMatrix");
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
	if (!isClass(xp, "generalMatrix")) {
	    if (isClass(xp, "symmetricMatrix")) {
		CharacterVector uplo(SEXP(xp.slot("uplo")));
		char *UL = uplo[0];
		stype = 1;
		if (*UL == 'L' || *UL == 'l') stype = -1;
	    }
	}
	int nzmax = ((int*)p)[ncol];
	if (ii.size() != nzmax)
	    Rf_error("size of i must match p[Dim[2] + 1]");
	itype = CHOLMOD_LONG;
	packed = (int)true;
	sorted = (int)true;
	xtype = -1;
	if (isClass(xp, "dsparseMatrix")) {
	    NumericVector xx(SEXP(xp.slot("x")));
	    x = xx.begin();
	    xtype = CHOLMOD_REAL;
	}
	if (isClass(xp, "nsparseMatrix")) xtype = CHOLMOD_PATTERN;
	if (isClass(xp, "zsparseMatrix")) {
	    xtype = CHOLMOD_COMPLEX;
	    Rf_error("Not yet defined zsparseMatrix?");
	}
	if (xtype == -1) Rf_error("Unknown (logical?) sparse Matrix type");
    }

    static inline
    void chk_mismatch(int a, int b, std::string compNm, std::string methNm) {
	if (a != b)
	    Rf_error("%s: %s mismatch, %d != %d",
		     methNm.c_str(), compNm.c_str(), a, b);
    }

    void chmSp::update(CHM_SP nn) {
	chk_mismatch(nrow, nn->nrow, "nrow", "chmSp::update");
	chk_mismatch(ncol, nn->ncol, "ncol", "chmSp::update");
	chk_mismatch(stype, nn->stype, "stype", "chmSp::update");
	chk_mismatch(itype, nn->itype, "itype", "chmSp::update");
	chk_mismatch(xtype, nn->xtype, "xtype", "chmSp::update");
	chk_mismatch(dtype, nn->dtype, "dtype", "chmSp::update");
	chk_mismatch(packed, nn->packed, "packed", "chmSp::update");
	chk_mismatch(sorted, nn->sorted, "sorted", "chmSp::update");
	int nnz = ::M_cholmod_nnz(this, &c);
	chk_mismatch(nnz, ::M_cholmod_nnz(nn, &c), "nnz", "chmSp::update");
	int *ii = (int*)i, *pp = (int*)p;
	if (!std::equal(pp, pp + ncol + 1, (int*)nn->p))
	    Rf_error("%s: inconsistency in %s", "chmSp::update", "p");
	if (!std::equal(ii, ii + nnz, (int*)nn->i))
	    Rf_error("%s: inconsistency in %s", "chmSp::update", "i");
	double *nnX = (double*)nn->x;
	std::copy(nnX, nnX + nnz, (double*)x);
    }

    void chmDn::init(double *X, int r, int c) {
	z = (void*)NULL;
	xtype = CHOLMOD_REAL;
	dtype = 0;		// CHOLMOD_DOUBLE
	nrow = (size_t) r;
	ncol = (size_t) c;
	nzmax = (size_t) r * c;
	d = (size_t) r;
	x = (void*) X;
    }

    chmDn::chmDn(double *xx, int nr, int nc) : cholmod_dense() {
	this->init(xx, nr, nc);
    }
    
    chmDn::chmDn(std::vector<double> v) : cholmod_dense() {
	this->init(&v[0], v.size(), 1);
    }

    chmDn::chmDn(NumericVector v) : cholmod_dense() {
	this->init(v.begin(), v.size(), 1);
    }

    chmDn::chmDn(NumericMatrix m) : cholmod_dense() {
	this->init(m.begin(), m.nrow(), m.ncol());
    }
    
    chmDn::chmDn(ddenseMatrix m) : cholmod_dense() {
	this->init(m.x.begin(), m.nrow(), m.ncol());
    }
}
