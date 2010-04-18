#include "MatrixNs.h"
#include <R_ext/Lapack.h>

using namespace Rcpp;

namespace MatrixNs{

// A function like this should be part of the Rcpp::S4 class
    static bool isClass(S4 x, std::string clname) {
	CharacterVector cl = x.attr("class");
	if (as<std::string>(cl) == clname) return true;

//	SEXP pkg = cl.attr("package");
	Function clDef("getClassDef");
	S4 cld = clDef(cl);
	List exts = cld.attr("contains");
	Function ssc(".selectSuperClasses");
	CharacterVector nms = ssc(exts, false, true, false);
	for (int i = 0; i < nms.size(); i++)
	    if (std::string(nms[i]) == clname) return true;
	return false;
    }

// Check this: Do the dspMatrix and dtpMatrix classes pass this check?
    ddenseMatrix::ddenseMatrix(S4 &xp) : dMatrix(xp) {
	if (!x.size() == Dim[0] * Dim[1])
	    ::Rf_error("%s: Dim = (%d, %d) is inconsistent with x.size() = %d",
		       "ddenseMatrix::ddenseMatrix", Dim[0], Dim[1], x.size());
    }

    void dgeMatrix::dgemv(Trans Tr,
			  double alpha, const NumericVector &X,
			  double beta, NumericVector &Y) {
	int i1 = 1;
	char tr = Tr.TR;
	bool NTR = tr == 'N';
	if (X.size() != Dim[NTR ? 1 : 0] || Y.size() != Dim[NTR ? 0 : 1])
	    Rf_error("dgemv \"%c\", dim mismatch (%d, %d), X(%d), Y(%d)",
		     tr, Dim[0], Dim[1], X.size(), Y.size());
	F77_CALL(dgemv)(&tr, &Dim[0], &Dim[1], &alpha, x.begin(), &Dim[0],
			X.begin(), &i1, &beta, Y.begin(), &i1);
    }

    void dgeMatrix::dgemm(Trans TrA, Trans TrB,
			  double alpha, const dgeMatrix &B,
			  double beta, dgeMatrix &C) {
//	int i1 = 1;
	char trA = TrA.TR, trB = TrB.TR;
	bool NTA = trA == 'N', NTB = trB == 'N';
	Dimension Bd(B.Dim), Cd(C.Dim);
	int M = Dim[NTA ? 0 : 1], N = Bd[NTB ? 1 : 0], K = Dim[NTA ? 1 : 0];
	if (Bd[NTB ? 0 : 1] != K || Cd[0] != M || Cd[1] != N)
	    Rf_error("dgemm \"%c,%c\", dim mismatch (%d, %d), (%d,%d), (%d,%d)",
		     trA, trB, Dim[0], Dim[1], Bd[0], Bd[1], Cd[0], Cd[1]);
	F77_CALL(dgemm)(&trA, &trB, &M, &N, &K, &alpha, x.begin(), &Dim[0],
			B.x.begin(), &Bd[0], &beta, C.x.begin(), &Cd[0]);
    }

    void Cholesky::update(const dpoMatrix &A) {
	Dimension Ad(A.Dim);
	if (Dim[0] != Ad[0])
	    Rf_error("%s dimension mismatch, (%d,%d) vs A(%d,%d)",
		     "Cholesky::update(dpoMatrix)", Dim[0], Dim[1],
		     Ad[0], Ad[1]);
	uplo = A.uplo;
	std::copy(A.x.begin(), A.x.end(), x.begin());
	int info;
	F77_CALL(dpotrf)(&(uplo.UL), &Dim[0], x.begin(), &Dim[0], &info);
	if (!info)
	    Rf_error("Lapack routine %s returned error code %d",
		     "dpotrf", info);
    }

    void Cholesky::update(Trans Tr, double alpha, const dgeMatrix &A,
			  double beta, const dsyMatrix &C) {
	const char tr = Tr.TR;
	const bool NTR = tr == 'N';
	Dimension Ad(A.Dim), Cd(C.Dim);
	if (Dim[0] != Cd[0] || Ad[NTR ? 0 : 1] != Dim[0])
	    Rf_error("%s(\"%c\") dimension mismatch, (%d,%d), A(%d,%d), C(%d,%d)",
		     "Cholesky::update(dpoMatrix, dgeMatrix)", tr,
		     Dim[0], Dim[1], Ad[0], Ad[1], Cd[0], Cd[1]);
	uplo = C.uplo;
	std::copy(C.x.begin(), C.x.end(), x.begin());
	F77_CALL(dsyrk)(&(uplo.UL), &tr, &Dim[0], &Ad[NTR ? 1 : 0], &alpha,
			A.x.begin(), &Ad[0], &beta, x.begin(), &Dim[0]);
	int info;
	F77_CALL(dpotrf)(&(uplo.UL), &Dim[0], x.begin(), &Dim[0], &info);
	if (info)
	    Rf_error("Lapack routine %s returned error code %d",
		     "dpotrf", info);
    }

    void Cholesky::dpotrs(NumericVector &v) {
	int info, i1 = 1, vs = v.size();
	if (vs != Dim[0])
	    Rf_error("%s (%d, %d) dimension mismatch (%d, %d)",
		     "Cholesky::dpotrs", Dim[0], Dim[1], vs, 1);
	F77_CALL(dpotrs)(&uplo.UL, &Dim[0], &i1, x.begin(), &Dim[0],
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

    chmFa::chmFa(S4 xp) : cholmod_factor() {
	CharacterVector cl(SEXP(xp.attr("class")));
	char *clnm = cl[0];
	if (!isClass(xp, "CHMfactor"))
	    ::Rf_error("Class %s object passed to %s is not a %s",
		       clnm, "chmFa::chmFa", "CHMfactor");
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
	const char* msg = "dCHMfactor with is_super == %s is not %s";
	if (is_super) {
	    IntegerVector
		Pi(SEXP(xp.slot("pi"))),
		SUPER(SEXP(xp.slot("super"))),
		S(SEXP(xp.slot("s")));
	    NumericVector PX(SEXP(xp.slot("px")));

	    if (!isClass(xp, "dCHMsuper"))
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

	    if (!isClass(xp, "dCHMsimpl"))
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
    void chk_mismatch(int a, int b, const std::string compNm, const std::string methNm) {
	if (a != b)
	    Rf_error("%s: %s mismatch, %d != %d",
		     methNm.c_str(), compNm.c_str(), a, b);
    }

    void chmSp::update(cholmod_sparse &nn) {
	chk_mismatch(nrow, nn.nrow, "nrow", "chmSp::update");
	chk_mismatch(ncol, nn.ncol, "ncol", "chmSp::update");
	chk_mismatch(stype, nn.stype, "stype", "chmSp::update");
	chk_mismatch(itype, nn.itype, "itype", "chmSp::update");
	chk_mismatch(xtype, nn.xtype, "xtype", "chmSp::update");
	chk_mismatch(dtype, nn.dtype, "dtype", "chmSp::update");
	chk_mismatch(packed, nn.packed, "packed", "chmSp::update");
	chk_mismatch(sorted, nn.sorted, "sorted", "chmSp::update");
	int nnz = ::M_cholmod_nnz(this, &c), nnznn = ::M_cholmod_nnz(&nn, &c);
	chk_mismatch(nnz, nnznn, "nnz", "chmSp::update");
	int *ii = (int*)i, *pp = (int*)p;
	if (!std::equal(pp, pp + ncol + 1, (int*)nn.p))
	    Rf_error("%s: inconsistency in %s", "chmSp::update", "p");
	if (!std::equal(ii, ii + nnz, (int*)nn.i))
	    Rf_error("%s: inconsistency in %s", "chmSp::update", "i");
	double *nnX = (double*)nn.x;
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
}
