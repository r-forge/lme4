// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
#ifndef LME4_GLM_H
#define LME4_GLM_H

#include "mer.h"
namespace glm {
    class modelMatrix {
    protected:
	Rcpp::S4                      d_xp;
	Rcpp::NumericVector  d_beta, d_Vtr;
    public:
	modelMatrix(Rcpp::S4&);
	
	Rcpp::NumericVector const& getBeta() {return d_beta;}

	void setBeta(Rcpp::NumericVector const&,
		     Rcpp::NumericVector const& = Rcpp::NumericVector(),
		     double = 0.);
    };
    
    class deModMat : public modelMatrix {
	MatrixNs::dgeMatrix   d_X, d_V;
	MatrixNs::dpoMatrix      d_VtV;
    public:
	deModMat(Rcpp::S4&,int);
	
	MatrixNs::dgeMatrix const& X() {return d_X;}

	void reweight(Rcpp::NumericMatrix const&,
		      Rcpp::NumericVector const&);
	double solveBeta();
    };
    
    template<typename Tm, typename Tr>
    class mod {
	Tm     mm;
	Tr   resp;
    public:
	mod(Rcpp::S4&);

	double  updateMu();
	double updateWts();
	int N()  const   {return resp.offset().size();}
	int n()  const   {return  resp.wtres().size();}
	int p()  const   {return  mm.getBeta().size();}

	Rcpp::List IRLS(int verb);
	double solveBeta();
    };

    template<typename Tm, typename Tr>
    inline mod<Tm,Tr>::mod(Rcpp::S4& xp)
	: resp (Rcpp::S4(xp.slot("resp"))),
	  mm   (Rcpp::S4(xp.slot("mm")), resp.mu().size()) {
    }

    /** 
     * Update the conditional mean, weighted residuals and wrss in resp
     * from the linear predictor.  Some resp modules also update other
     * information, such as the sqrtXwt matrix. 
     *
     * @return the *penalized*, weighted residual sum of squares
     */
    template<typename Tm, typename Tr> inline
    double mod<Tm,Tr>::updateMu() {
	Rcpp::NumericVector const& offset = resp.offset();
	Rcpp::NumericVector gamma(offset.size());
	std::copy(offset.begin(), offset.end(), gamma.begin());

	MatrixNs::chmDn gg(gamma);
	if (mm.getBeta().size() > 0)
	    mm.X().dmult('N', 1., 1., MatrixNs::chmDn(mm.getBeta()), gg);
	return resp.updateMu(gamma);
    }	

    /**
     * Update the weighted residuals, wrss, sqrtrwt, sqrtXwt, U,
     * cu, UtV, VtV and Vtr.
     *
     * @return penalized, weighted residual sum of squares
     */
    template<typename Tf, typename Tr> inline
    double mod<Tf,Tr>::updateWts() {
	double ans = resp.updateWts();
	mm.reweight(resp.sqrtXwt(), resp.wtres());
	return ans;
    }
    
    template<typename Tf, typename Tr> inline
    double mod<Tf,Tr>::solveBeta() {
	return mm.solveBeta();
    }

    template<typename Tf, typename Tr> inline
    Rcpp::List mod<Tf,Tr>::IRLS(int verb) {
	Rcpp::NumericVector bBase(p()), incB(p()), muBase(n());
	double crit, step, c0, c1;
				// sqrtrwt and sqrtXwt must be set
	updateMu();		// using current beta
	crit = 10. * CM_TOL;
	for (int i = 0; crit >= CM_TOL && i < CM_MAXITER; i++) {
				// store copies of mu and beta
	    std::copy(resp.mu().begin(), resp.mu().end(), muBase.begin());
	    std::copy(mm.getBeta().begin(), mm.getBeta().end(), bBase.begin());
	    c0 = updateWts();
	    solveBeta();
	    std::copy(mm.getBeta().begin(), mm.getBeta().end(), incB.begin());
	    for (c1 = c0, step = 1.; c0 <= c1 && step > CM_SMIN;
		 step /= 2.) {
		c1 = mm.setBeta(bBase, incB, step);
		if (verb > 1) {
		    showincr(step, c0, c1, mm.getBeta(), "beta");
		}
	    }
	    crit = compareVecWt(muBase, resp.mu(), resp.sqrtrwt());
	    if (verb > 1)
		Rprintf("   convergence criterion: %g\n", crit);
	}
	return Rcpp::List::create(Rcpp::_["beta"] = mm.getBeta());
    } // IRLS
}
    
#endif /* LME4_GLM_H */

