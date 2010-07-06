#include "mer.h"
#include "glm.h"
#include <R_ext/BLAS.h>

using namespace std;
using namespace MatrixNs;

namespace glm {
    modelMatrix::modelMatrix(Rcpp::S4 &xp)
	: d_xp(                     xp),
	  d_beta(SEXP(xp.slot("beta"))),
	  d_Vtr(          d_beta.size()) {
    }
    
    void modelMatrix::setBeta(Rcpp::NumericVector const& bbase,
			      Rcpp::NumericVector const&  incr,
			      double                      step) {
	// FIXME: have feModule inherit from modelMatrix so this is not repeated.
	R_len_t p = d_beta.size();
	if (p == 0) return;
	if (bbase.size() != p)
	    throw runtime_error("feModule::setBeta size mismatch of beta and bbase");
	if (step == 0.) {
	    std::copy(bbase.begin(), bbase.end(), d_beta.begin());
	} else {
	    Rcpp::NumericVector res = bbase + incr * step;
	    std::copy(res.begin(), res.end(), d_beta.begin());
	}
    }

    deModMat::deModMat(Rcpp::S4 xp, R_len_t n)
	: modelMatrix(               xp),
	  d_X(   Rcpp::S4(xp.slot("X"))),
	  d_V(          n,   d_X.ncol()),
	  d_R(                      d_X) {
    }

    /** 
     * Update V, VtV and Vtr
     * 
     * @param Ut from the reModule
     * @param Xwt square root of the weights for the model matrices
     * @param wtres weighted residuals
     */
    void deModMat::reweight(Rcpp::NumericMatrix   const&   Xwt,
			    Rcpp::NumericVector   const& wtres) {
	if (d_beta.size() == 0) return;
	chmDn cXwt(Xwt);
	if (Xwt.size() != d_X.nrow())
	    Rf_error("%s: dimension mismatch %s(%d,%d), %s(%d,%d)",
		     "deFeMod::reweight", "X", d_X.nrow(), d_X.ncol(),
		     "Xwt", Xwt.nrow(), Xwt.ncol());
	int Wnc = Xwt.ncol(), Wnr = Xwt.nrow(),
	    Xnc = d_X.ncol(), Xnr = d_X.nrow();
	double *V = d_V.x().begin(), *X = d_X.x().begin();

	if (Wnc == 1) {
	    for (int j = 0; j < Xnc; j++) 
		transform(Xwt.begin(), Xwt.end(), X + j*Xnr,
			  V + j*Xnr, multiplies<double>());
	} else {
	    int i1 = 1;
	    double one = 1., zero = 0.;
	    Rcpp::NumericVector tmp(Xnr), mm(Wnc); 
	    fill(mm.begin(), mm.end(), 1.);
	    for (int j = 0; j < Xnc; j++) {
		transform(Xwt.begin(), Xwt.end(), X + j*Xnr,
			  tmp.begin(), multiplies<double>());
		F77_CALL(dgemv)("N", &Wnr, &Wnc, &one, tmp.begin(),
				&Wnr, mm.begin(), &i1, &zero,
				V + j * Wnr, &i1);
	    }
	}
	d_V.dgemv('T', 1., wtres, 0., d_Vtr);
	d_R.update(d_V);
    }
    
    double deModMat::solveBeta() {
	copy(d_Vtr.begin(), d_Vtr.end(), d_beta.begin());
	d_R.dpotrs(d_beta.begin());
	return 0.;
    }

    spModMat::spModMat(Rcpp::S4 xp, R_len_t n)
	: modelMatrix(               xp),
	  d_X(   Rcpp::S4(xp.slot("X"))),
	  d_F(   Rcpp::S4(xp.slot("RX"))) {
    }

    /** 
     * Update V, VtV and Vtr
     * 
     * @param Ut from the reModule
     * @param Xwt square root of the weights for the model matrices
     * @param wtres weighted residuals
     */
    void spModMat::reweight(Rcpp::NumericMatrix   const&   Xwt,
			    Rcpp::NumericVector   const& wtres) {
	if (d_beta.size() == 0) return;
	double one = 1., zero = 0.;
	int Wnc = Xwt.ncol(), Wnr = Xwt.nrow(),
	    Xnc = d_X.ncol, Xnr = d_X.nrow;
	if (Xwt.size() != (int)d_X.nrow)
	    Rf_error("%s: dimension mismatch %s(%d,%d), %s(%d,%d)",
		     "deFeMod::reweight", "X", Xnr, Xnc,
		     "Xwt", Wnr, Wnc);
	if (Wnc == 1) {
	    if (d_V) M_cholmod_free_sparse(&d_V, &c);
	    d_V = M_cholmod_copy_sparse(&d_X, &c);
	    chmDn csqrtX(Xwt);
	    M_cholmod_scale(&csqrtX, CHOLMOD_ROW, d_V, &c);
	} else throw runtime_error("spModMat::reweight: multiple columns in Xwt");
// FIXME rewrite this using the triplet representation
	
	chmDn cVtr(d_Vtr);
	const chmDn cwtres(wtres);
	M_cholmod_sdmult(d_V, 'T', &one, &zero, &cwtres, &cVtr, &c);

	CHM_SP Vt = M_cholmod_transpose(d_V, 1/*values*/, &c);
	d_F.update(*Vt);
	M_cholmod_free_sparse(&Vt, &c);
    }

    double spModMat::solveBeta() {
	Rcpp::NumericMatrix bb = d_F.solve(CHOLMOD_A, d_Vtr);
	copy(bb.begin(), bb.end(), d_beta.begin());

	return 0.;
    }


}

RCPP_FUNCTION_2(Rcpp::List, glmIRLS, Rcpp::S4 xp, int verb) {
    Rcpp::S4 mm = xp.slot("mm");
    if (mm.is("deFeMod")) {
	glm::mod<glm::deModMat,mer::glmerResp> m(xp);
	return m.IRLS(verb);
    } else if(mm.is("spFeMod")) {
	glm::mod<glm::spModMat,mer::glmerResp> m(xp);
	return m.IRLS(verb);
    } else throw runtime_error("Unknown model matrix type");
}
