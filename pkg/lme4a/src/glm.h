// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
#ifndef LME4_GLM_H
#define LME4_GLM_H

#include "mer.h"
namespace glm {
    class predModule {
    protected:
	Rcpp::S4                      d_xp;
	Rcpp::NumericVector  d_coef, d_Vtr;
    public:
	predModule(Rcpp::S4&);
	
	Rcpp::NumericVector const& coef() const {return d_coef;}
	Rcpp::NumericVector const&  Vtr() const {return  d_Vtr;}

	void setCoef(Rcpp::NumericVector const&,
		     Rcpp::NumericVector const& = Rcpp::NumericVector(),
		     double = 0.);
    };
    
    class dPredModule : public predModule {
	MatrixNs::ddenseModelMatrix   d_X;
	MatrixNs::dgeMatrix           d_V;
	MatrixNs::Cholesky          d_fac;
    public:
	dPredModule(Rcpp::S4,R_len_t); // Maybe use int instead of R_len_t?
	
	MatrixNs::ddenseModelMatrix const& X() const {return d_X;}

	void reweight(Rcpp::NumericMatrix const&,
		      Rcpp::NumericVector const&);
	double solveCoef(double);
    };
    
    class sPredModule : public predModule {
	MatrixNs::chmSp          d_X;
	MatrixNs::chmFr        d_fac;
	CHM_SP                   d_V;
    public:
	sPredModule(Rcpp::S4,R_len_t);
	
	MatrixNs::chmSp      const& X() const {return d_X;}

	void reweight(Rcpp::NumericMatrix const&,
		      Rcpp::NumericVector const&);
	double solveCoef(double);
    };
    
    template<typename Tp, typename Tr>
    class mod {
	Tr   resp;
	Tp   pred;
    public:
	mod(Rcpp::S4&);

	double   setCoef(Rcpp::NumericVector const& base,
			 Rcpp::NumericVector const& incr,
			 double                     step);
	double solveCoef(double);
	double  updateMu();
	double updateWts();
	int N()  const   {return resp.offset().size();}
	int n()  const   {return  resp.wtres().size();}
	int p()  const   {return   pred.coef().size();}

	Rcpp::List IRLS(int verb);
    };

    template<typename Tp, typename Tr>
    inline mod<Tp,Tr>::mod(Rcpp::S4& xp)
	: resp (Rcpp::S4(xp.slot("resp"))),
	  pred (Rcpp::S4(xp.slot("pred")), resp.mu().size()) {
    }

    /** 
     * Update the conditional mean, weighted residuals and wrss in resp
     * from the linear predictor.  Some resp modules also update other
     * information, such as the sqrtXwt matrix. 
     *
     * @return the weighted residual sum of squares
     */
    template<typename Tp, typename Tr> inline
    double mod<Tp,Tr>::updateMu() {
	Rcpp::NumericVector const& offset = resp.offset(), coef = pred.coef();
	Rcpp::NumericVector gamma(offset.size());
	std::copy(offset.begin(), offset.end(), gamma.begin());

	MatrixNs::chmDn gg(gamma);
	if (coef.size() > 0)
	    pred.X().dmult('N', 1., 1., MatrixNs::chmDn(coef), gg);
	return resp.updateMu(gamma);
    }	

    template<typename Tp, typename Tr> inline
    double mod<Tp,Tr>::setCoef(Rcpp::NumericVector const& base,
			       Rcpp::NumericVector const& incr,
			       double                     step) {
	pred.setCoef(base, incr, step);
	return updateMu();
    }	

    /**
     * Update the weighted residuals, wrss, sqrtrwt, sqrtXwt, fac and Vtr
     *
     * @return penalized, weighted residual sum of squares
     */
    template<typename Tp, typename Tr> inline
    double mod<Tp,Tr>::updateWts() {
	double ans = resp.updateWts();
	pred.reweight(resp.sqrtXwt(), resp.wtres());
	return ans;
    }
    
    template<typename Tp, typename Tr> inline
    double mod<Tp,Tr>::solveCoef(double dd) {
	return pred.solveCoef(dd);
    }

    template<typename Tp, typename Tr> inline
    Rcpp::List mod<Tp,Tr>::IRLS(int verb) {
	Rcpp::NumericVector cBase(p()), incr(p());
	double crit, step, c0, c1;
	size_t iter = 0;
				// sqrtrwt and sqrtXwt must be set
// FIXME: Actually sqrtrwt and sqrtXwt will be overwritten.  Consider
// how to handle the situation that these vectors are set from mustart
// or etastart.
	updateMu();		// using current coef
	do {
	    if (++iter > CM_MAXITER)
		throw std::runtime_error("IRLS MAXITER exceeded");
				// store a copy of the coefficients
	    std::copy(pred.coef().begin(), pred.coef().end(), cBase.begin());
	    c0 = updateWts();
	    crit = solveCoef(c0);
	    if (verb > 1) Rprintf("  convergence criterion: %g\n", crit);
	    std::copy(pred.coef().begin(), pred.coef().end(), incr.begin());
	    step = 2.;
	    do {
		if ((step /= 2.) < CM_SMIN)
		    throw std::runtime_error("IRLS step factor reduced beyond SMIN");
		c1 = setCoef(cBase, incr, step);
		if (verb > 1)
		    mer::showincr(step, c0, c1, pred.coef(), "coef");
	    } while (c1 >= c0);
	} while (crit >= CM_TOL);

	return Rcpp::List::create(Rcpp::_["coef"] = pred.coef());
// FIXME Need a templated method to evaluate the deviance  Rcpp::_["deviance"] = );
    }
}

#endif /* LME4_GLM_H */

