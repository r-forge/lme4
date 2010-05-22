// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
#ifndef LME4_MER_H
#define LME4_MER_H

#include <Rcpp.h>
#include "MatrixNs.h"
#include "glmFamily.h"

namespace mer {

    double compare_vec(const Rcpp::NumericVector v1, const Rcpp::NumericVector v2);
    void showdbl(const double* x, const char* nm, int n);
    void showint(const int* x, const char* nm, int n);
    void showCHM_DN(const CHM_DN x, const std::string &nm);
    void showCHM_FR(const CHM_FR x, const std::string &nm);
    void showCHM_SP(const CHM_SP x, const std::string &nm);

    class reModule {
	MatrixNs::chmFr     d_L;
	MatrixNs::chmSp     d_Lambda, d_Ut, d_Zt;
	Rcpp::IntegerVector Lind;
	Rcpp::NumericVector lower, theta;
	Rcpp::NumericVector d_Utr, d_cu, d_u;
	double             *d_ldL2, d_sqlLenU;
    public:
	reModule(Rcpp::S4);

	const Rcpp::NumericVector  &cu() const {return d_cu;}
	const Rcpp::NumericVector   &u() const {return d_u;}
	const Rcpp::NumericVector &Utr() const {return d_Utr;}
	const MatrixNs::chmFr       &L() const {return d_L;}
	const MatrixNs::chmSp  &Lambda() const {return d_Lambda;}
	const MatrixNs::chmSp      &Ut() const {return d_Ut;}
	const MatrixNs::chmSp      &Zt() const {return d_Zt;}
	double                    ldL2() const {return *d_ldL2;}
	double                 sqrLenU() const {return d_sqlLenU;}

	Rcpp::NumericVector rwUpdateL(Rcpp::NumericMatrix const&,
				      Rcpp::NumericVector const&);
	void setU(const cholmod_dense*);
	void setU(Rcpp::NumericVector const&);
	void updateSqrLenU();
	void updateTheta(Rcpp::NumericVector const&);
	void updateU(MatrixNs::chmDn const&);
    };

    class feModule {
    protected:
	Rcpp::NumericVector d_beta, d_Vtr;
	double             *d_ldRX2;
    public:
	feModule(Rcpp::S4 xp);

	void setBeta(Rcpp::NumericVector const&);

	const Rcpp::NumericVector &beta() const {return d_beta;}
	const Rcpp::NumericVector  &Vtr() const {return d_Vtr;}
	double                    ldRX2() const {return *d_ldRX2;}
    };

    class deFeMod : public feModule {
    protected:
	MatrixNs::dgeMatrix d_X, d_RZX, d_UtV, d_V;
	MatrixNs::dpoMatrix d_VtV;
	MatrixNs::Cholesky  d_RX;
    public:
	deFeMod(Rcpp::S4 xp);

	const MatrixNs::Cholesky   &RX() const{return d_RX;}
	const MatrixNs::dgeMatrix   &X() const{return d_X;}
	const MatrixNs::dgeMatrix &RZX() const{return d_RZX;}
	const MatrixNs::dgeMatrix &UtV() const{return d_UtV;}
	const MatrixNs::dpoMatrix &VtV() const{return d_VtV;}
	const MatrixNs::dgeMatrix   &V() const{return d_V;}

	void updateRzxRx(MatrixNs::chmSp const&,
			 MatrixNs::chmFr const&);
	void updateUtV(MatrixNs::chmSp const&);
	void updateV(Rcpp::NumericMatrix const&);
	void updateVtr(Rcpp::NumericVector const&);
	void updateRX(bool);
	void solveA(Rcpp::NumericVector&) const;
	Rcpp::NumericVector updateBeta(Rcpp::NumericVector const&);
    };

    class spFeMod : public feModule {
    protected:
	MatrixNs::chmSp    d_RZX, d_UtV, d_V, d_VtV, d_X;
	MatrixNs::chmFr    d_RX;
    public:
	spFeMod(Rcpp::S4 xp);

	const MatrixNs::chmSp      &X() const{return d_X;}
	const MatrixNs::chmSp    &RZX() const{return d_RZX;}
	const MatrixNs::chmSp    &UtV() const{return d_UtV;}
	const MatrixNs::chmSp      &V() const{return d_V;}
	const MatrixNs::chmSp    &VtV() const{return d_VtV;}
	const MatrixNs::chmFr     &RX() const{return d_RX;}

	void updateRzxRx(MatrixNs::chmSp const&,
			 MatrixNs::chmFr const&);
	void updateUtV(MatrixNs::chmSp const&,
		       Rcpp::NumericMatrix const&);
	void updateVtr(Rcpp::NumericMatrix const&,
		       Rcpp::NumericVector const&);
	Rcpp::NumericVector updateBeta(Rcpp::NumericVector const&);
    };

    class merResp {
    protected:
	double              *d_wrss;
	Rcpp::NumericVector  d_offset, d_sqrtrwt, d_wtres, mu, weights, y;
	Rcpp::NumericMatrix  d_sqrtXwt;
    public:
	merResp(Rcpp::S4);

	      bool                    reml() const{return false;}
	const Rcpp::NumericMatrix &sqrtXwt() const{return d_sqrtXwt;}
	const Rcpp::NumericVector &sqrtrwt() const{return d_sqrtrwt;}
	const Rcpp::NumericVector   &wtres() const{return d_wtres;}
	const Rcpp::NumericVector  &offset() const{return d_offset;}
	double                        wrss() const{return *d_wrss;}

	double updateMu(Rcpp::NumericVector const&);
	double updateWrss();
// FIXME: This should be const but runs into trouble in glmerResp with
// family.devResid using an extractor of a component of a list,
// lst["dev.resids"]
	double unscaledDeviance() {return *d_wrss;} 
	void updateL(MatrixNs::chmSp const&, MatrixNs::chmFr);
	void updateSqrtRWt(){}	// overridden by glmerResp
	void updateSqrtXWt(){}	// overridden by glmerResp
    };

    class lmerResp : public merResp {
	bool d_reml;
    public:
	lmerResp(Rcpp::S4);
	bool reml() const{return d_reml;}
    };

    class glmerResp : public merResp {
    protected:
	glmFamily family;
	Rcpp::NumericVector eta, muEta, n, d_var;
	void   linkFun(){family.linkFun(eta, mu);}
	void   linkInv(){family.linkInv(mu, eta);}
	void     MuEta(){family.muEta(muEta, eta);}
	void  variance(){family.variance(d_var, mu);}
    public:
	glmerResp(Rcpp::S4 xp);

	const Rcpp::NumericVector &var() const {return d_var;}
	double unscaledDeviance();
	double updateMu(Rcpp::NumericVector const &gamma);
	Rcpp::NumericVector devResid();
	void updateSqrtRWt();
	void updateSqrtXWt();
    };
    
    class nlmerResp : public merResp {
	Rcpp::Environment nlenv;
	Rcpp::Language nlmod;
	Rcpp::CharacterVector pnames;
    public:
	nlmerResp(Rcpp::S4 xp);
	double updateMu(Rcpp::NumericVector const &gamma);
    };
    
/* Model object template
 * Tf is the type of fixed-effects module (deFeMod or spFeMod)
 * Tr is the type of response module (merResp, glmerResp, nlmerResp or nglmerResp)
 */
    template<typename Tf, typename Tr>  
    class mer {
	reModule re;
	Tf fe;
	Tr resp;
	static double l2PI,;
    public:
	mer(Rcpp::S4 xp);

	bool reml() const{return resp.reml();}
	double IRLS(int verbose);
	double Laplace() {
	    return resp.unscaledDeviance() + re.sqrLenU() + re.ldL2();
	}
	double PIRLS(int verbose);
	double pwrss() const{return resp.wrss() + re.sqrLenU();}
	double profDev() const;
	double profREML() const;
	void solveBeta();
	void solveBetaU();
	void solveU();
	void updateTheta(const Rcpp::NumericVector&nt){re.updateTheta(nt);}
	double updateMu();
    };
    
    template<typename Tf, typename Tr>
    inline mer<Tf,Tr>::mer(Rcpp::S4 xp)
	: re(Rcpp::S4(xp.slot("re"))),
	  fe(Rcpp::S4(xp.slot("fe"))),
	  resp(Rcpp::S4(xp.slot("resp"))) {
    }

    template<typename Tf, typename Tr>
    double mer<Tf,Tr>::l2PI = log(2. * PI);

#define CM_TOL 1.e-5
#define CM_MAXITER 30
#define CM_SMIN 1.e-4

/**
 * Determine the Euclidean distance between two vectors relative to the
 * square root of the product of their lengths
 */

/** 
 * Evaluate the profiled deviance (only meaningful for lmer and nlmer)
 * 
 * @return ldL2 + n *(1 + log(2 * pi * pwrss/n))
 */
    template<typename Tf, typename Tr> inline
    double mer<Tf,Tr>::profDev() const {
	double nn = (double)resp.wtres().size();
	return re.ldL2() + nn * (1 + l2PI + log(pwrss()/nn));
    }

/** 
 * Evaluate the profiled REML criterion (only meaningful for lmer)
 * 
 * @return ldL2 + ldRX2 + (n - p) * (1 + log(2 * pi * pwrss/(n - p)))
 */
    template<typename Tf, typename Tr> inline
    double mer<Tf, Tr>::profREML() const {
	double nmp = (double)(resp.wtres().size() - fe.beta().size());
	return re.ldL2()+fe.ldRX2()+nmp*(1 + l2PI + log(pwrss()/nmp));
    }

/** 
 * Update the conditional mean, weighted residuals and wrss in resp
 * from the linear predictor.  Some resp modules also update other
 * information, such as the sqrtXwt matrix. 
 */
    template<typename Tf, typename Tr> inline
    double mer<Tf,Tr>::updateMu() {
	const Rcpp::NumericVector u = re.u(), offset = resp.offset();
	Rcpp::NumericVector b(u.size()), gamma(offset.size());
	std::copy(offset.begin(), offset.end(), gamma.begin());
	MatrixNs::chmDn bb(b), gg(gamma);

	re.Lambda().dmult('N', 1., 0., MatrixNs::chmDn(u), bb);
	re.Zt().dmult('T', 1., 1., bb, gg);
	fe.X().dmult('N', 1., 1., MatrixNs::chmDn(fe.beta()), gg);
	return resp.updateMu(gamma);
    }	

    template<typename Tf, typename Tr> inline
    void mer<Tf,Tr>::solveBeta() {
    }

    template<typename Tf, typename Tr> inline
    void mer<Tf,Tr>::solveBetaU() {
	fe.updateRzxRx(re.Lambda(), re.L());
	Rcpp::NumericVector cu = fe.updateBeta(re.cu());
	re.updateU(MatrixNs::chmDn(cu));
    }
    
    template<typename Tf, typename Tr> inline
    void mer<Tf,Tr>::solveU() {
    }

    template<typename Tf, typename Tr> inline
    double mer<Tf,Tr>::IRLS(int verbose) {
	const Rcpp::NumericVector var = resp.var(), bb = fe.beta();
	Rcpp::NumericVector varold(var.size());
	int p = bb.size();
	Rcpp::NumericVector betabase(p), newbeta(p), incr(p); 
	double crit, step, pwrss0, pwrss1;

	crit = 10. * CM_TOL;
	    // weights, var and sqrtrwt are assumed established.
	for (int i = 0; crit >= CM_TOL && i < CM_MAXITER; i++) {
				// store copies of var and beta
	    std::copy(bb.begin(), bb.end(), betabase.begin());
	    std::copy(var.begin(), var.end(), varold.begin());
				// using current beta and u
	    pwrss0 = updateMu() + re.sqrLenU();
	    resp.updateSqrtXWt();	// muEta and sqrtXwt
	    fe.updateV(resp.sqrtXwt());	// derive V from X
	    fe.updateVtr(resp.wtres());	// crossprod(V, wtres)
	    fe.updateRX(false);		// factor V'V
	    std::copy(fe.Vtr().begin(), fe.Vtr().end(), incr.begin());
	    fe.RX().dpotrs(incr); // evaluate the increment
	    if (verbose > 1) showdbl(&incr[0], "deltab", p);
	    pwrss1 = pwrss0;	// force one evaluation of the loop
	    for (step = 1.; pwrss0 <= pwrss1 && step > CM_SMIN; step /= 2.) {
				// newbeta <- betabase + incr * mult
		std::transform(incr.begin(), incr.end(), newbeta.begin(),
		       std::bind2nd(std::multiplies<double>(), step));
		std::transform(betabase.begin(), betabase.end(),
			       newbeta.begin(), newbeta.begin(),
			       std::plus<double>());
		fe.setBeta(newbeta);
		pwrss1 = updateMu() + re.sqrLenU();
		if (verbose > 1) {
		    Rcpp::NumericVector bb = fe.beta();
		    Rprintf("step = %g, pwrss0 = %g; pwrss1 = %g\n",
			    step, pwrss0, pwrss1, *bb.begin());	
		    showdbl(bb.begin(), "beta", p);
		}
	    }
	    resp.updateSqrtRWt();
	    crit = compare_vec(varold, var);
	    if (verbose > 1) Rprintf("convergence criterion: %g\n", crit);
	}
	return Laplace();
    } // IRLS

    template<typename Tf, typename Tr> inline
    double mer<Tf,Tr>::PIRLS(int verbose) {
	const Rcpp::NumericVector u = re.u(), var = resp.var();
	int n = var.size(), q = u.size();
	Rcpp::NumericVector newU(q), ubase(q), varold(n);
	double crit, step, pwrss0, pwrss1;

	crit = 10. * CM_TOL;
	// weights, var and sqrtrwt are assumed established.
	for (int i = 0; crit >= CM_TOL && i < CM_MAXITER; i++) {
				// store copies of var and u
	    std::copy(u.begin(), u.end(), ubase.begin());
	    std::copy(var.begin(), var.end(), varold.begin());
	    updateMu();		// linear predictor
	    pwrss0 = resp.updateWrss() + re.sqrLenU();
	    resp.updateSqrtXWt(); // muEta and sqrtXwt
				// solve for incr
	    Rcpp::NumericVector incr = re.rwUpdateL(resp.sqrtXwt(), resp.wtres());
	    pwrss1 = pwrss0;	// force one evaluation of the loop
	    for (step = 1.; pwrss0 <= pwrss1 && step > CM_SMIN; step /= 2.) {
		std::transform(incr.begin(), incr.end(), newU.begin(),
			       std::bind2nd(std::multiplies<double>(),
					    step));
		std::transform(ubase.begin(), ubase.end(),
			       newU.begin(), newU.begin(),
			       std::plus<double>());
		re.setU(newU);
		updateMu();
		pwrss1 = resp.updateWrss() + re.sqrLenU();
		if (verbose > 1) {
		    Rprintf("step = %g, pwrss0 = %g; pwrss1 = %g\n",
			    step, pwrss0, pwrss1);
		    showdbl(newU.begin(), "u", newU.size());
		}
	    }
	    resp.updateSqrtRWt();
	    crit = compare_vec(varold, var);
	    if (verbose > 1) Rprintf("convergence criterion: %g\n", crit);
	}
	return Laplace();
    } // PIRLS
}

#endif
