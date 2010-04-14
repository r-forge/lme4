#include "Matrix_ns.h"

using namespace Rcpp;

extern cholmod_common c;

namespace Matrix{
    
    Matrix::Matrix(S4 xp) :
	Dim(SEXP(xp.slot("Dim"))),
	Dimnames(SEXP(xp.slot("Dimnames"))) {
    }

    dMatrix::dMatrix(S4 xp) :
	Matrix(xp),
	x(SEXP(xp.slot("x"))) {
    }

// Check this: Do the dspMatrix and dtpMatrix classes pass this check?
    ddenseMatrix::ddenseMatrix(S4 xp) : dMatrix(xp) {
	if (!x.size() == Dim[0] * Dim[1])
	    ::Rf_error("%s: Dim = (%d, %d) is inconsistent with x.size() = %d",
		       "dgeMatrix::dgeMatrix", Dim[0], Dim[1], x.size());
    }
    
    compMatrix::compMatrix(S4 xp) : factors(SEXP(xp.slot("factors"))) {
    }

    generalMatrix::generalMatrix(S4 xp) : compMatrix(xp) {
    }

    triangularMatrix::triangularMatrix(S4 xp) :
	uplo(SEXP(xp.slot("uplo"))),
	diag(SEXP(xp.slot("diag"))) {
    }

    symmetricMatrix::symmetricMatrix(S4 xp) :
	compMatrix(xp),
	uplo(SEXP(xp.slot("uplo"))) {
    }

    dgeMatrix::dgeMatrix(S4 xp) : ddenseMatrix(xp), generalMatrix(xp) {
    }

    dtrMatrix::dtrMatrix(S4 xp) : ddenseMatrix(xp), triangularMatrix(xp) {
    }
    
    dsyMatrix::dsyMatrix(S4 xp) : ddenseMatrix(xp), symmetricMatrix(xp) {
    }

    Cholesky::Cholesky(S4 xp) : dtrMatrix(xp) {
    }

    dpoMatrix::dpoMatrix(S4 xp) : dsyMatrix(xp) {
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
    
    dgCMatrix::dgCMatrix(S4 xp) :
	i(SEXP(xp.slot("i"))),
	p(SEXP(xp.slot("p"))),
	Dim(SEXP(xp.slot("Dim"))),
	Dimnames(SEXP(xp.slot("Dim"))),
	factors(SEXP(xp.slot("Dim"))),
	x(SEXP(xp.slot("x")))
    {
	int nnz = p[Dim[1]];
	if (i.size() != nnz || x.size() != nnz)
	    Rf_error("size of i and x must match p[Dim[2] + 1]");
	sp = new cholmod_sparse;
	sp->nrow = (size_t)Dim[0];
	sp->ncol = (size_t)Dim[1];
	sp->nzmax = (size_t)nnz;
	sp->p = (void*)(&p[0]);
	sp->i = (void*)(&i[0]);
	sp->x = (void*)(&x[0]);
	sp->stype = 0;
	sp->itype = CHOLMOD_LONG;
	sp->xtype = CHOLMOD_REAL;
	sp->dtype = 0;  // CHOLMOD_DOUBLE
	sp->packed = (int)true;
	sp->sorted = (int)true;
    }

    static inline
    void chk_mismatch(int a, int b, std::string compNm, std::string methNm) {
	if (a != b)
	    Rf_error("%s: %s mismatch, %d != %d",
		     methNm.c_str(), compNm.c_str(), a, b);
    }

    void dgCMatrix::update(CHM_SP nn) {
	chk_mismatch(sp->nrow, nn->nrow, "nrow", "dgCMatrix::update");
	chk_mismatch(sp->ncol, nn->ncol, "ncol", "dgCMatrix::update");
	chk_mismatch(sp->stype, nn->stype, "stype", "dgCMatrix::update");
	chk_mismatch(sp->itype, nn->itype, "itype", "dgCMatrix::update");
	chk_mismatch(sp->xtype, nn->xtype, "itype", "dgCMatrix::update");
	chk_mismatch(sp->dtype, nn->dtype, "dtype", "dgCMatrix::update");
	chk_mismatch(sp->packed, nn->packed, "packed", "dgCMatrix::update");
	chk_mismatch(sp->sorted, nn->sorted, "sorted", "dgCMatrix::update");
	int nnz = ::M_cholmod_nnz(sp, &c);
	chk_mismatch(nnz, ::M_cholmod_nnz(nn, &c), "nnz", "dgCMatrix::update");
	int *spP = (int*)sp->p;
	if (!std::equal(spP, spP + sp->ncol + 1, (int*)nn->p))
	    Rf_error("%s: inconsistency in %s", "dgCMatrix::update", "p");
	spP = (int*)sp->i;
	if (!std::equal(spP, spP + nnz, (int*)nn->i))
	    Rf_error("%s: inconsistency in %s", "dgCMatrix::update", "i");
	double *nnX = (double*)nn->x;
	std::copy(nnX, nnX + nnz, (double*)sp->x);
    }

    chmDn::chmDn(double *x, int nr, int nc) {
	dn = new cholmod_dense;
	dn->nrow = (size_t) nr;
	dn->ncol = (size_t) nc;
	dn->nzmax = (size_t) nr * nc;
	dn->d = (size_t) nr;
	dn->x = (void*) x;
	dn->z = (void*) NULL;
	dn->xtype = CHOLMOD_REAL;
	dn->dtype = 0;  // CHOLMOD_DOUBLE
    }
}
