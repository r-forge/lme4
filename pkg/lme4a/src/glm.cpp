#include "mer.h"
#include "glm.h"
#include <R_ext/BLAS.h>

using namespace std;
using namespace MatrixNs;

namespace glm {
    predModule::predModule(Rcpp::S4& xp)
	: d_xp(                     xp),
	  d_coef(SEXP(xp.slot("coef"))),
	  d_Vtr(          d_coef.size()) {
    }
    
    void predModule::setCoef(Rcpp::NumericVector const& cbase,
			     Rcpp::NumericVector const&  incr,
			     double                      step) {
	// FIXME: have feModule inherit from predModule so this is not repeated.
	R_len_t p = d_coef.size();
	if (p == 0) return;
	if (cbase.size() != p)
	    throw runtime_error("feModule::setCoef size mismatch of coef and cbase");
	if (step == 0.) {
	    std::copy(cbase.begin(), cbase.end(), d_coef.begin());
	} else {
//	    Rcpp::NumericVector res = cbase + incr * step;  // needs Rcpp_0.8.3
//	    copy(res.begin(), res.end(), d_coef.begin());
	    double *cb = cbase.begin(), *inc = incr.begin(), *cc = d_coef.begin();
	    for (R_len_t i = 0; i < p; i++) cc[i] = cb[i] + inc[i] * step;
	}
    }

    dPredModule::dPredModule(Rcpp::S4 xp, R_len_t n)
	: predModule(                  xp),
	  d_X(     Rcpp::S4(xp.slot("X"))),
	  d_V(            n,   d_X.ncol()),
	  d_fac( Rcpp::S4(xp.slot("fac"))) {
    }

    /** 
     * Update V, VtV and Vtr
     * 
     * @param Ut from the reModule
     * @param Xwt square root of the weights for the model matrices
     * @param wtres weighted residuals
     */
    void dPredModule::reweight(Rcpp::NumericMatrix   const&   Xwt,
			    Rcpp::NumericVector   const& wtres) {
	if (d_coef.size() == 0) return;
	chmDn cXwt(Xwt);
	if ((Xwt.rows() * Xwt.cols()) != d_X.nrow())
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
	d_fac.update(d_V);
    }
    
    double dPredModule::solveCoef(double wrss) {
	copy(d_Vtr.begin(), d_Vtr.end(), d_coef.begin());
	d_fac.dtrtrs('T', d_coef.begin()); // solve R'c = Vtr
	double ans = sqrt(inner_product(d_coef.begin(), d_coef.end(),
					d_coef.begin(), double())/wrss);
	d_fac.dtrtrs('N', d_coef.begin()); // solve R beta = c;
	return ans;
    }

    sPredModule::sPredModule(Rcpp::S4 xp, R_len_t n)
	: predModule(                  xp),
	  d_X(     Rcpp::S4(xp.slot("X"))),
	  d_fac(  Rcpp::S4(xp.slot("fac"))) {
    }

    /** 
     * Update V, VtV and Vtr
     * 
     * @param Ut from the reModule
     * @param Xwt square root of the weights for the model matrices
     * @param wtres weighted residuals
     */
    void sPredModule::reweight(Rcpp::NumericMatrix   const&   Xwt,
			       Rcpp::NumericVector   const& wtres) {
	if (d_coef.size() == 0) return;
	double one = 1., zero = 0.;
	int Wnc = Xwt.ncol(), Wnr = Xwt.nrow(),
	    Xnc = d_X.ncol, Xnr = d_X.nrow;
	if ((Xwt.rows() * Xwt.cols()) != (int)d_X.nrow)
	    Rf_error("%s: dimension mismatch %s(%d,%d), %s(%d,%d)",
		     "deFeMod::reweight", "X", Xnr, Xnc,
		     "Xwt", Wnr, Wnc);
	if (Wnc == 1) {
	    if (d_V) M_cholmod_free_sparse(&d_V, &c);
	    d_V = M_cholmod_copy_sparse(&d_X, &c);
	    chmDn csqrtX(Xwt);
	    M_cholmod_scale(&csqrtX, CHOLMOD_ROW, d_V, &c);
	} else throw runtime_error("sPredModule::reweight: multiple columns in Xwt");
// FIXME write this combination using the triplet representation
	
	chmDn cVtr(d_Vtr);
	const chmDn cwtres(wtres);
	M_cholmod_sdmult(d_V, 'T', &one, &zero, &cwtres, &cVtr, &c);

	CHM_SP Vt = M_cholmod_transpose(d_V, 1/*values*/, &c);
	d_fac.update(*Vt);
	M_cholmod_free_sparse(&Vt, &c);
    }

    double sPredModule::solveCoef(double wrss) {
	Rcpp::NumericMatrix cc = d_fac.solve(CHOLMOD_L, d_fac.solve(CHOLMOD_P, d_Vtr));
	double ans = sqrt(inner_product(cc.begin(), cc.end(),
					cc.begin(), double())/wrss);
	Rcpp::NumericMatrix bb = d_fac.solve(CHOLMOD_Pt, d_fac.solve(CHOLMOD_Lt, cc));
	copy(bb.begin(), bb.end(), d_coef.begin());
	return ans;
    }


}

RCPP_FUNCTION_2(Rcpp::List, glmIRLS, Rcpp::S4 xp, int verb) {
    Rcpp::S4 pm = xp.slot("pred");
    if (pm.is("dPredModule")) {
	glm::mod<glm::dPredModule,mer::glmerResp> m(xp);
	return m.IRLS(verb);
    } else if(pm.is("sPredModule")) {
	glm::mod<glm::sPredModule,mer::glmerResp> m(xp);
	return m.IRLS(verb);
    } else throw runtime_error("Unknown linear predictor module type");
}
