#include "mer.h"
#include <R_ext/BLAS.h>    // for F77_CALL(dgemv)

using namespace std;
using namespace MatrixNs;  // for chmDn, chmFr and chmSp

namespace mer {
    deFeMod::deFeMod(Rcpp::S4 xp, int n)
	: feModule(                  xp),
	  d_RZX(Rcpp::clone(Rcpp::S4(xp.slot("RZX")))),
	  d_X(              Rcpp::S4(xp.slot("X"))),
	  d_RX( Rcpp::clone(Rcpp::S4(xp.slot("RX")))),
	  d_UtV(d_RZX.nrow(), d_RZX.ncol()),
	  d_V(             n,   d_X.ncol()),
	  d_VtV(             d_beta.size()) {
    }

    /** 
     * Update V, UtV, VtV and Vtr
     * 
     * @param Ut from the reModule
     * @param Xwt square root of the weights for the model matrices
     * @param wtres weighted residuals
     */
    void deFeMod::reweight(cholmod_sparse        const*    Ut,
			   Rcpp::NumericMatrix   const&   Xwt,
			   Rcpp::NumericVector   const& wtres) {
	if (d_beta.size() == 0) return;
	chmDn cXwt(Xwt);
	double one = 1., zero = 0.;
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
	chmDn cUtV(d_UtV), cV(d_V);
	M_cholmod_sdmult(Ut, 0/*trans*/, &one, &zero, &cV, &cUtV, &c);
	d_VtV.dsyrk(d_V, 1., 0.);
    }

    /** 
     * Solve (V'V)beta = Vtr for beta.
     * 
     */
    void deFeMod::solveBeta() {
	if (d_beta.size() == 0) return;
	MatrixNs::Cholesky chol(d_V);
	copy(d_Vtr.begin(), d_Vtr.end(), d_beta.begin());
	chol.dpotrs(d_beta);
    }

    void deFeMod::updateDcmp(Rcpp::List& ll) {
	Rcpp::List devcomp = ll["devcomp"];
	Rcpp::NumericVector cmp = devcomp["cmp"];
	cmp["ldRX2"] = d_RX.logDet2();
	ll["RZX"] = d_RZX.sexp();
	ll["RX"]  = d_RX. sexp();
	ll["beta"] = d_beta;
    }

    /** 
     * Update RZX and RX
     *   RZX <<- solve(L, solve(L, crossprod(Lambda, ZtX), 
     *                          sys = "P"),
     *                 sys = "L")
     *   RX <<- chol(XtX - crossprod(RZX))
     * 
     * @param Lambda relative covariance factor for the random effects
     * @param L sparse Cholesky factor for the random effects
     */
    void deFeMod::updateRzxRx(MatrixNs::chmSp const& Lambda,
			      MatrixNs::chmFr const&      L) {
	if (d_beta.size() == 0) return;
	chmDn cRZX(d_RZX);
	Lambda.dmult('T', 1., 0., chmDn(d_UtV), cRZX);
	Rcpp::NumericMatrix
	    ans = L.solve(CHOLMOD_L, L.solve(CHOLMOD_P, &cRZX));
	d_RZX.setX(ans);
	d_RX.update('T', -1., d_RZX, 1., d_VtV);
	d_ldRX2 = d_RX.logDet2();
    }

    /** 
     * Update beta
     *	beta <- solve(RX, solve(t(RX), Vtr - crossprod(RZX, cu)))
     * 
     * @param cu intermediate solution of random-effects.
     * 
     * @return cu - RZX %*% beta
     */
    Rcpp::NumericVector deFeMod::updateBeta(Rcpp::NumericVector const& cu) {
	Rcpp::NumericVector ans = clone(cu);
	if (d_beta.size() == 0) return ans;

	copy(d_Vtr.begin(), d_Vtr.end(), d_beta.begin());
	d_RZX.dgemv('T', -1., ans, 1., d_beta);
	d_RX.dpotrs(d_beta);
	d_RZX.dgemv('N', -1., d_beta, 1., ans);
	return ans;
    }
}
