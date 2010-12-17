// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
#ifndef LME4_MERMOD_H
#define LME4_MERMOD_H

#include <Rcpp.h>
#include "MatrixNs.h"
#include "utilities.h"
#include "respModule.h"
#include "reModule.h"
#include "feModule.h"
namespace mer {
    enum Alg {Beta, U, BetaU};

    /* Model object template
     * Tf is the type of fixed-effects module (deFeMod or spFeMod)
     * Tr is the type of response module (merResp, glmerResp, nlmerResp or nglmerResp)
     */
    template<typename Tf, typename Tr>  
    class mer {
	reModule re;
	Tr     resp;
	Tf       fe;
    public:
	mer(Rcpp::S4&);

	// These return a scalar numeric value but as an R vector so
	// that it can have attributes attached to it.  See the mkans
	// function.
	Rcpp::NumericVector LMMdeviance(const Rcpp::NumericVector&);
	Rcpp::NumericVector PIRLS      (const Rcpp::NumericVector&,
					const Rcpp::NumericVector&,
					const Rcpp::NumericVector&,
					int,Alg);
	Rcpp::List updateDcmp          (const Rcpp::NumericVector&,
					const Rcpp::NumericVector&,
					const Rcpp::NumericVector&);
	double Laplace                 () const;
	double setBetaU                (const Rcpp::NumericVector&,
					const Rcpp::NumericVector&,
					const Rcpp::NumericVector&,
					const Rcpp::NumericVector&,
					double,Alg);
	double updateMu                ();
	double updateWts               (Alg);

	int N()  const   {return resp.offset().size();}
	int n()  const   {return  resp.wtres().size();}
	int nth()const   {return    re.lower().size();}
	int p()  const   {return     fe.coef().size();}
	int q()  const   {return        re.u().size();}
	int s()  const   {return              N()/n();}

	void solveCoef(Alg);
	void updateRzxRx();
    };
    
    template<typename Tf, typename Tr>
    inline mer<Tf,Tr>::mer(Rcpp::S4& xp)
	: re   (Rcpp::S4(xp.slot("re"))  ),
	  resp (Rcpp::S4(xp.slot("resp"))),
    	  fe   (Rcpp::S4(xp.slot("fe")), resp.mu().size()) {
    }

    /** 
     * Evaluate the profiled deviance or REML criterion for a linear mixed
     * model 
     * 
     * @param nt New value of theta
     * 
     * @return profiled deviance or REML criterion
     */
    template<typename Tf, typename Tr> inline
    Rcpp::NumericVector mer<Tf,Tr>::LMMdeviance(const Rcpp::NumericVector& nt) {
	re.setTheta(nt);
	updateWts(BetaU);
	solveCoef(BetaU);
	updateMu();
	return mkans(Laplace(), fe.coef(), re.u(),
		     re.ldL2(), resp.wrss(), re.sqrLenU());
    }

    template<typename Tf, typename Tr> inline
    double mer<Tf,Tr>::Laplace() const {
	return resp.Laplace(re.ldL2(), fe.ldRX2(), re.sqrLenU());
    }

    /** 
     * Update beta or u or both from base values, increments and a step
     * factor.
     * 
     * @param bBase base of beta vector
     * @param uBase base of u vector
     * @param incB increment for beta vector
     * @param incU increment for u vector
     * @param step step factor
     * @param alg  algorithm
     * 
     * @return penalized, weighted residual sum of squares
     */
    template<typename Tf, typename Tr> inline
    double mer<Tf,Tr>::setBetaU(Rcpp::NumericVector const& bBase,
				Rcpp::NumericVector const& uBase,
				Rcpp::NumericVector const&  incB,
				Rcpp::NumericVector const&  incU,
				double                      step,
				Alg                          alg) {
	if (alg !=    U) fe.setCoef(bBase, incB, step);
	if (alg != Beta) re.setU(uBase, incU, step);
	return updateMu();
    }
    
    template<typename Tf, typename Tr> inline
    void mer<Tf,Tr>::solveCoef(Alg alg) {
	switch(alg) {
	case Beta:
	    fe.solveCoef();
	    break;
	case U:
	    re.solveU();
	    break;
	case BetaU:
	    fe.updateRzxRx(re.Lambda(), re.L());
	    re.updateU(fe.updateBeta(re.cu()));
	}
    }

    /** 
     * Update the conditional mean, weighted residuals and wrss in resp
     * from the linear predictor.  Some resp modules also update other
     * information, such as the sqrtXwt matrix. 
     *
     * @return the *penalized*, weighted residual sum of squares
     */
    template<typename Tf, typename Tr> inline
    double mer<Tf,Tr>::updateMu() {
				// needs Rcpp_0.8.3
	Rcpp::NumericVector gamma = re.linPred() + fe.linPred(); 
	return resp.updateMu(gamma) + re.sqrLenU();
    }	

    /**
     * Update the RZX and RX slots in fe
     * cu, UtV, VtV and Vtr.
     *
     * @return penalized, weighted residual sum of squares
     */
    template<typename Tf, typename Tr> inline
    void mer<Tf,Tr>::updateRzxRx() {
	re.reweight(resp.sqrtXwt(), resp.wtres());
	fe.reweight(resp.sqrtXwt(), resp.wtres());
	fe.updateUtV(re.Ut());
	fe.updateRzxRx(re.Lambda(), re.L());
    }

    /**
     * Update the weighted residuals, wrss, sqrtrwt, sqrtXwt, U,
     * cu, UtV, VtV and Vtr.
     *
     * @return penalized, weighted residual sum of squares
     */
    template<typename Tf, typename Tr> inline
    double mer<Tf,Tr>::updateWts(Alg alg) {
	double ans = resp.updateWts();
	ans += re.sqrLenU();
	if (alg != Beta) re.reweight(resp.sqrtXwt(), resp.wtres());
	if (alg != U) fe.reweight(resp.sqrtXwt(), resp.wtres());
	if (alg == BetaU) fe.updateUtV(re.Ut());
	return ans;
    }

#define CM_TOL 1.e-4
#define CM_MAXITER 200
#define CM_SMIN 1.e-4

    template<typename Tf, typename Tr> inline
    Rcpp::NumericVector mer<Tf,Tr>::PIRLS(const Rcpp::NumericVector& theta,
					  const Rcpp::NumericVector&  beta,
					  const Rcpp::NumericVector&    u0, int verb, Alg alg) {
	Rcpp::NumericVector bBase(p()), incB(p()), incU(q()), muBase(n()), uBase(q());
	double crit, step, c0, c1;
				// sqrtrwt and sqrtXwt must be set
	re.setTheta(theta);
	if (alg == U) fe.setCoef(beta);
	re.setU(u0);

	updateMu();		// using current beta and u
	crit = 10. * CM_TOL;
	for (int i = 0; crit >= CM_TOL && i < CM_MAXITER; i++) {
				// store copies of mu, u and beta
	    std::copy(resp.mu().begin(), resp.mu().end(), muBase.begin());
	    std::copy(re.u().begin(), re.u().end(), uBase.begin());
	    std::copy(fe.coef().begin(), fe.coef().end(), bBase.begin());
	    c0 = updateWts(alg);
	    solveCoef(alg);
	    std::copy(fe.coef().begin(), fe.coef().end(), incB.begin());
	    std::copy(re.u().begin(), re.u().end(), incU.begin());
	    for (c1 = c0, step = 1.; c0 <= c1 && step > CM_SMIN;
		 step /= 2.) {
		c1 = setBetaU(bBase, uBase, incB, incU, step, alg);
		if (verb > 1) {
		    showincr(step, c0, c1, fe.coef(), "beta");
		    showdbl(re.u().begin(), "u", re.u().size());
		}
	    }
	    crit = compareVecWt(muBase, resp.mu(), resp.sqrtrwt());
	    if (verb > 1)
		Rprintf("   convergence criterion: %g\n", crit);
	}
	return mkans(Laplace(), fe.coef(), re.u(),
		     re.ldL2(), resp.wrss(), re.sqrLenU());
    } // PIRLS

    template<typename Tf, typename Tr> inline
    Rcpp::List mer<Tf,Tr>::updateDcmp(const Rcpp::NumericVector&   th,
				      const Rcpp::NumericVector& beta,
				      const Rcpp::NumericVector&    u) {
	re.setTheta(th);
	re.setU(u);
	fe.setCoef(beta);
	updateMu();
	updateWts(BetaU);
	re.reweight(resp.sqrtXwt(), resp.wtres());
	fe.reweight(resp.sqrtXwt(), resp.wtres());
	fe.updateUtV(re.Ut());
	fe.updateRzxRx(re.Lambda(), re.L());
	return
	    Rcpp::List::create(Rcpp::_["ldRX2"]   = Rcpp::wrap(fe.ldRX2())
			       ,Rcpp::_["mu"]     = Rcpp::clone(resp.mu())
			       ,Rcpp::_["devRes"] = resp.devResid()
			       ,Rcpp::_["RX"]     = fe.RX()
			       ,Rcpp::_["RZX"]    = fe.RZX()
			       ,Rcpp::_["L"]      = re.L()
			       ,Rcpp::_["Lambda"] = re.Lambda()
		);
    }
}

namespace glm {
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
//	pred.reweight(resp.sqrtXwt(), resp.wtres());
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

#endif
