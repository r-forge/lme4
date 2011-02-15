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
	double updateMu                (double);
	double updateWts               (Alg);

	int N()  const   {return resp.offset().size();}
	int n()  const   {return  resp.wtres().size();}
	int nth()const   {return    re.lower().size();}
	int p()  const   {return     fe.coef().size();}
	int q()  const   {return        re.u().size();}
	int s()  const   {return              N()/n();}

	void solveIncr   (Alg);
	void installCoef0(Alg);
	void updateRzxRx ();
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
	solveIncr(BetaU);
	updateMu(1.);
	return mkans(Laplace(), fe.coef(), re.u(),
		     re.ldL2(), resp.wrss(), re.sqrLenU());
    }

    template<typename Tf, typename Tr> inline
    double mer<Tf,Tr>::Laplace() const {
	return resp.Laplace(re.ldL2(), fe.ldRX2(), re.sqrLenU());
    }

    template<typename Tf, typename Tr> inline
    void mer<Tf,Tr>::solveIncr(Alg alg) {
	switch(alg) {
	case Beta:
	    fe.solveIncr();
	    break;
	case U:
	    re.solveIncr();
	    break;
	case BetaU:
	    fe.updateRzxRx(re.Lambda(), re.L());
	    re.updateIncr(fe.updateIncr(re.cu()));
	}
    }

    template<typename Tf, typename Tr> inline
    void mer<Tf,Tr>::installCoef0(Alg alg) {
	switch(alg) {
	case Beta:
	    fe.installCoef0();
	    break;
	case U:
	    re.installU0();
	    break;
	case BetaU:
	    fe.installCoef0();	
	    re.installU0();
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
    double mer<Tf,Tr>::updateMu(double step) {
#ifdef USE_RCPP_SUGAR
	Rcpp::NumericVector gamma = re.linPred1(step) + fe.linPred1(step); 
#else
	Rcpp::NumericVector feP = fe.linPred1(step);
	Rcpp::NumericVector gamma = re.linPred1(step);
	std::transform(feP.begin(), feP.end(), gamma.begin(),
		       gamma.begin(), std::plus<double>());
#endif
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
	Rcpp::NumericVector muBase(n());
	double crit, step, c0, c1;
				// sqrtrwt and sqrtXwt must be set
	re.setTheta(theta);
	fe.setCoef0(beta);
	re.setU0(u0);

	updateMu(1.);		// using current beta and u
	crit = 10. * CM_TOL;
	for (int i = 0; crit >= CM_TOL && i < CM_MAXITER; i++) {
				// store a copy of mu
	    std::copy(resp.mu().begin(), resp.mu().end(), muBase.begin());
	    c0 = updateWts(alg);
	    solveIncr(alg);
	    for (c1 = c0, step = 1.; c0 <= c1 && step > CM_SMIN;
		 step /= 2.) {
		c1 = updateMu(step);
		if (verb > 1) {
		    showincr(step, c0, c1, fe.coef(), "beta");
		    showdbl(re.u().begin(), "u", re.u().size());
		}
	    }
	    installCoef0(alg);
	    crit = compareVecWt(muBase, resp.mu(), resp.sqrtrwt());
	    if (verb > 1)
		Rprintf("   convergence criterion: %g\n", crit);
	}
	return mkans(Laplace(), fe.coef(), re.u(),
		     re.ldL2(), resp.wrss(), re.sqrLenU());
    } // PIRLS

//FIXME: This should be avoided eventually
    template<typename Tf, typename Tr> inline
    Rcpp::List mer<Tf,Tr>::updateDcmp(const Rcpp::NumericVector&   th,
				      const Rcpp::NumericVector& beta,
				      const Rcpp::NumericVector&    u) {
	re.setTheta(th);
	re.setU0(u);
	re.setIncr(Rcpp::NumericVector(u.size()));
	fe.setCoef0(beta);
	fe.setIncr(Rcpp::NumericVector(beta.size()));
	updateMu(1.);
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

#if 0
	double   setCoef(Rcpp::NumericVector const& base,
			 Rcpp::NumericVector const& incr,
			 double                     step);
#endif
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

    /**
     * Update the weighted residuals, wrss, sqrtrwt, sqrtXwt, fac and Vtr
     *
     * @return penalized, weighted residual sum of squares
     */
    template<typename Tp, typename Tr> inline
    double mod<Tp,Tr>::updateWts() {
	double ans = resp.updateWts();
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
// FIXME: Replace this.
//		c1 = setCoef(cBase, incr, step);
		c1 = 0.;
		if (verb > 1)
		    mer::showincr(step, c0, c1, pred.coef(), "coef");
	    } while (c1 >= c0);
	} while (crit >= CM_TOL);

	return Rcpp::List::create(Rcpp::_["coef"] = pred.coef());
// FIXME Need a templated method to evaluate the deviance  Rcpp::_["deviance"] = );
    }
}

#endif
