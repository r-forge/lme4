#include "predModules.h"
#include <R_ext/BLAS.h>
#include "utilities.h"

using namespace std;
using namespace MatrixNs;

namespace matMod {
    predModule::predModule(Rcpp::S4& xp)
	: d_coef(Rcpp::clone(SEXP(xp.slot("coef")))),
	  d_Vtr(                      d_coef.size()) {
    }

    predModule::predModule(int ncoef)
	: d_coef(ncoef), d_Vtr(ncoef) {
    }

    void predModule::setCoef(Rcpp::NumericVector const& cbase,
			     Rcpp::NumericVector const&  incr,
			     double                      step) {
	R_len_t p = d_coef.size();
	if (cbase.size() != p)
	    throw runtime_error("predModule::setCoef size mismatch of coef and cbase");
	if (p == 0) return;
	if (step == 0.) {
	    copy(cbase.begin(), cbase.end(), d_coef.begin());
	} else {
	    Rcpp::NumericVector res = cbase + incr * step;
	    copy(res.begin(), res.end(), d_coef.begin());
	}
    }

    dPredModule::dPredModule(Rcpp::S4 xp, int n)
	: predModule(                  xp),
	  d_X(     Rcpp::S4(xp.slot("X"))),
	  d_V(            n,   d_X.ncol()),
	  d_fac( Rcpp::S4(xp.slot("fac"))) {
    }

    dPredModule::dPredModule(Rcpp::NumericMatrix mm, int n)
        : predModule(mm.ncol()), d_X(mm), d_V(n, mm.ncol()),
	  d_fac(mm.ncol(), 'U') {
    }

    Rcpp::NumericVector dPredModule::linPred() const {
	Rcpp::NumericVector ans(d_X.nrow());
	d_X.dgemv('N', 1., d_coef, 0., ans);
	return ans;
    }

    
    /** 
     * Solve (V'V)coef = Vtr for coef.
     *
     * @param wrss weighted residual sum of squares (defaults to 1.)
     * @return convergence criterion
     */
    double dPredModule::solveCoef(double wrss) {
	if (d_coef.size() == 0) return 0.;
	copy(d_Vtr.begin(), d_Vtr.end(), d_coef.begin());
	d_fac.update(d_V);
	d_fac.dtrtrs('T', d_coef.begin()); // solve R'c = Vtr
	double ans = sqrt(inner_product(d_coef.begin(), d_coef.end(),
					d_coef.begin(), double())/wrss);
	d_fac.dtrtrs('N', d_coef.begin()); // solve R beta = c;
	return ans;
    }

    sPredModule::sPredModule(Rcpp::S4 xp, int n)
	: predModule(                  xp),
	  d_X(     Rcpp::S4(xp.slot("X"))),
	  d_fac(  Rcpp::S4(xp.slot("fac"))) {
	d_V = (CHM_SP)NULL;
    }

    Rcpp::NumericVector sPredModule::linPred() const {
	Rcpp::NumericVector ans(d_X.nr());
	chmDn cans(ans);
	d_X.dmult('N',1.,0.,chmDn(d_coef),cans);
	return ans;
    }

    /** 
     * Update V, Vtr and fac
     *
     * Note: May want to update fac in a separate operation.  For the
     * fixed-effects modules this will update the factor twice because
     * it is separately updated in updateRzxRx.
     * 
     * @param Xwt square root of the weights for the model matrices
     * @param wtres weighted residuals
     */
    void sPredModule::reweight(Rcpp::NumericMatrix   const&   Xwt,
			       Rcpp::NumericVector   const& wtres) throw(std::runtime_error) {
	if (d_coef.size() == 0) return;
	double one = 1., zero = 0.;
	int Wnc = Xwt.ncol();//, Wnr = Xwt.nrow(),
//	    Xnc = d_X.ncol, Xnr = d_X.nrow;
	if ((Xwt.rows() * Xwt.cols()) != (int)d_X.nrow)
	    throw std::runtime_error("dimension mismatch");
	    // Rf_error("%s: dimension mismatch %s(%d,%d), %s(%d,%d)",
	    // 	     "deFeMod::reweight", "X", Xnr, Xnc,
	    // 	     "Xwt", Wnr, Wnc);
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
	Rcpp::NumericMatrix
	    cc = d_fac.solve(CHOLMOD_L, d_fac.solve(CHOLMOD_P, d_Vtr));
	double ans = sqrt(inner_product(cc.begin(), cc.end(),
					cc.begin(), double())/wrss);
	Rcpp::NumericMatrix
	    bb = d_fac.solve(CHOLMOD_Pt, d_fac.solve(CHOLMOD_Lt, cc));
	copy(bb.begin(), bb.end(), d_coef.begin());
	return ans;
    }

}

