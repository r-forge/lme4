#include "mer.h"
#include <cstring>

using namespace Rcpp;

extern cholmod_common c;

namespace mer{

    merResp::merResp(S4 xp) :
	Utr(SEXP(xp.slot("Utr"))),
	Vtr(SEXP(xp.slot("Vtr"))),
	cbeta(SEXP(xp.slot("cbeta"))),
	cu(SEXP(xp.slot("cu"))),
	mu(SEXP(xp.slot("mu"))),
	offset(SEXP(xp.slot("offset"))),
	resid(SEXP(xp.slot("resid"))),
	weights(SEXP(xp.slot("weights"))),
	wrss(SEXP(xp.slot("wrss"))),
	y(SEXP(xp.slot("y")))
    {
	int n = y.size(), ws = weights.size();
	
	if (mu.size() != n || resid.size() != n)
	    ::Rf_error("y, mu and resid slots must have equal lengths");
	if (cbeta.size() != Vtr.size())
	    ::Rf_error("cbeta and Vtr slots must have equal lengths");
	if (cu.size() != Utr.size())
	    ::Rf_error("cu and Utr slots must have equal lengths");
	if (ws && ws != n)
	    ::Rf_error("weights slot must have length 0 or n");
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

    SEXP dgCMatrix::dims() const {
	return wrap(clone(Dim));
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

    reModule::reModule(S4 xp) :
	L(S4(SEXP(xp.slot("L")))),
	Lambda(S4(SEXP(xp.slot("Lambda")))),
	Ut(S4(SEXP(xp.slot("Ut")))),
	Zt(S4(SEXP(xp.slot("Zt")))),
	Lind(SEXP(xp.slot("Lind"))),
	lower(SEXP(xp.slot("lower"))),
	theta(SEXP(xp.slot("theta"))),
	u(SEXP(xp.slot("u"))),
	ldL2(SEXP(xp.slot("ldL2"))) {
    }

    void reModule::updateTheta(NumericVector nt) {
	if (nt.size() != theta.size())
	    ::Rf_error("length(theta) = %d != length(newtheta) = %d",
		       theta.size(), nt.size());
	std::copy(nt.begin(), nt.end(), theta.begin());
	double *LamxP = Lambda.x.begin();
	for (int *Li = Lind.begin(); Li < Lind.end(); Li++)
	    *LamxP++ = nt[(*Li) - 1];
	CHM_SP LamTr = ::M_cholmod_transpose(Lambda.sp, 1/*values*/, &c);
	CHM_SP LamTrZt = ::M_cholmod_ssmult(LamTr, Zt.sp, 0/*stype*/,
					    1/*values*/, 1/*sorted*/, &c);
	::M_cholmod_free_sparse(&LamTr, &c);
	Ut.update(LamTrZt);
	::M_cholmod_free_sparse(&LamTrZt, &c);
	//FIXME: Define a CHMfactor class and update L
    }
}

extern "C"
SEXP check_dgCMatrix(SEXP xp) {
    S4 x4(xp);
    CharacterVector cls = x4.attr("class");
    std::string cl(cls[0]);
    if (cl != "dgCMatrix")
	Rf_error("incorrect class, %s, to %s", cl.c_str(), "check_dgCMatrix");
    mer::dgCMatrix rr(x4);
    return rr.dims();
}

extern "C"
SEXP check_resp(SEXP xp) {
    S4 x4(xp);
    mer::merResp rr(x4);
    return Rf_ScalarLogical(1);
}

extern "C"
SEXP update_reModule(SEXP xp, SEXP ntheta) {
    S4 x4(xp);
    mer::reModule rr(x4);
    NumericVector nt(ntheta);
    rr.updateTheta(nt);
}

