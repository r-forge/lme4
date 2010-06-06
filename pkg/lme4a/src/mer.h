// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
#ifndef LME4_MER_H
#define LME4_MER_H

#include <Rcpp.h>
#include "MatrixNs.h"
#include "glmFamily.h"

namespace mer {
				// utilities
    double compareVecWt(  Rcpp::NumericVector const&,
			  Rcpp::NumericVector const&,
			  Rcpp::NumericVector const&);
    void showCHM_DN(const_CHM_DN, std::string const&);
    void showCHM_FR(const_CHM_FR, std::string const&);
    void showCHM_SP(const_CHM_SP, std::string const&);
    void showdbl(const double*, const char*, int);
    void showincr(double,double,double,
		  Rcpp::NumericVector const&, const char*);
    void showint(const int*, const char*, int);

    class reModule {
	MatrixNs::chmFr     d_L;
	MatrixNs::chmSp     d_Lambda, d_Ut, d_Zt;
	Rcpp::IntegerVector d_Lind;
	Rcpp::NumericVector d_Utr, d_cu, d_lower, d_theta, d_u;
	double             *d_ldL2, d_sqrLenU;
    public:
	reModule(Rcpp::S4);

	const Rcpp::NumericVector  &cu() const {return  d_cu;}
 	const Rcpp::NumericVector   &u() const {return  d_u;}
	const Rcpp::NumericVector &Utr() const {return  d_Utr;}
	const MatrixNs::chmFr       &L() const {return  d_L;}
	const MatrixNs::chmSp  &Lambda() const {return  d_Lambda;}
	const MatrixNs::chmSp      &Ut() const {return  d_Ut;}
	const MatrixNs::chmSp      &Zt() const {return  d_Zt;}
	double                    ldL2() const {return *d_ldL2;}
	double                 sqrLenU() const {return  d_sqrLenU;}

	void reweight(Rcpp::NumericMatrix const&,
		      Rcpp::NumericVector const&);
	void setU(Rcpp::NumericVector const&,
		  Rcpp::NumericVector const& = Rcpp::NumericVector(),
		  double = 0.);
	void solveU();
	void updateLcu();
	void updateLambda(Rcpp::NumericVector const&);
	void updateU(Rcpp::NumericVector const&);
	void zeroU();
    };

    class feModule {
    protected:
	Rcpp::NumericVector d_beta, d_Vtr;
	double             *d_ldRX2;
    public:
	feModule(Rcpp::S4 xp);

	void setBeta(Rcpp::NumericVector const&,
		     Rcpp::NumericVector const& = Rcpp::NumericVector(),
		     double = 0.);
	const Rcpp::NumericVector &beta() const {return  d_beta;}
	const Rcpp::NumericVector  &Vtr() const {return  d_Vtr;}
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

	Rcpp::NumericVector updateBeta(Rcpp::NumericVector const&);

	void reweight(MatrixNs::chmSp     const&,
		      Rcpp::NumericMatrix const&,
		      Rcpp::NumericVector const&);
	void solveBeta();
	void updateRzxRx(MatrixNs::chmSp const&,
			 MatrixNs::chmFr const&);
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

	Rcpp::NumericVector updateBeta(Rcpp::NumericVector const&);

	void reweight(MatrixNs::chmSp     const&,
		      Rcpp::NumericMatrix const&,
		      Rcpp::NumericVector const&);
	void solveBeta();
	void updateRzxRx(MatrixNs::chmSp const&,
			 MatrixNs::chmFr const&);
    };

    class merResp {
    protected:
	double              *d_wrss;
	Rcpp::NumericVector  d_offset, d_sqrtrwt, d_wtres, d_mu, d_weights, d_y;
	Rcpp::NumericMatrix  d_sqrtXwt;
    public:
	merResp(Rcpp::S4);

	const Rcpp::NumericVector      &mu() const{return  d_mu;}
	const Rcpp::NumericVector  &offset() const{return  d_offset;}
	const Rcpp::NumericMatrix &sqrtXwt() const{return  d_sqrtXwt;}
	const Rcpp::NumericVector &sqrtrwt() const{return  d_sqrtrwt;}
	const Rcpp::NumericVector   &wtres() const{return  d_wtres;}
	double                        wrss() const{return *d_wrss;}

	double                   updateWts(){return updateWrss();}
	double                  updateWrss();
    };

    class lmerResp : public merResp {
	int d_reml;
    public:
	lmerResp(Rcpp::S4);
	double Laplace(double,double,double)const;
	double updateMu(Rcpp::NumericVector const&);
    };

    class glmerResp : public merResp {
    protected:
	double              *d_devres;
	glmFamily            family;
	Rcpp::NumericVector  d_eta, d_muEta, d_n, d_var;
    public:
	glmerResp(Rcpp::S4 xp);

	const Rcpp::NumericVector &var() const{return d_var;}
	const Rcpp::NumericVector &eta() const{return d_var;}
	double                 Laplace(double,double,double) const;
	double                  devres() const{return *d_devres;}
	double                updateMu(Rcpp::NumericVector const&);
	double               updateWts();
	Rcpp::NumericVector   devResid();
    };
    
    class nlmerResp : public merResp {
	Rcpp::Environment nlenv;
	Rcpp::Language nlmod;
	Rcpp::CharacterVector pnames;
    public:
	nlmerResp(Rcpp::S4 xp);
	double updateMu(Rcpp::NumericVector const &gamma);
	double Laplace(double,double,double) const;
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
    public:
	mer(Rcpp::S4 xp);

	double IRLS     (int);
	double Laplace  () const {
	    return resp.Laplace(re.ldL2(), fe.ldRX2(), re.sqrLenU());
	}
	double PIRLS    (int);
	double PIRLSBeta(int);
	double setBeta  (Rcpp::NumericVector const&,
			 Rcpp::NumericVector const&, double);
	double setU     (Rcpp::NumericVector const&,
			 Rcpp::NumericVector const&, double);
	double setBetaU (Rcpp::NumericVector const&,
			 Rcpp::NumericVector const&,
			 Rcpp::NumericVector const&,
			 Rcpp::NumericVector const&, double);
	double updateMu ();
	double updateWts();

	int N() const      {return resp.offset().size();}
	int n() const      {return  resp.wtres().size();}
	int p() const      {return     fe.beta().size();}
	int q() const      {return        re.u().size();}
	int s() const      {return              N()/n();}

	void solveBeta()   {fe.solveBeta();}
	void solveBetaU();
	void solveU();
	void updateLambda(Rcpp::NumericVector const&);
	void updateRzxRx() {fe.updateRzxRx(re.Lambda(), re.L());}
	void zeroU()       {re.zeroU();}
    };
    
    template<typename Tf, typename Tr>
    inline mer<Tf,Tr>::mer(Rcpp::S4 xp)
	: re    (Rcpp::S4(xp.slot("re"))),
	  fe    (Rcpp::S4(xp.slot("fe"))),
	  resp (Rcpp::S4(xp.slot("resp"))) {
    }

#define CM_TOL 1.e-4
#define CM_MAXITER 30
#define CM_SMIN 1.e-4

    template<typename Tf, typename Tr> inline
    double mer<Tf,Tr>::setBeta(Rcpp::NumericVector const&  bb,
			       Rcpp::NumericVector const& inc,
			       double                    step) {
	fe.setBeta(bb, inc, step);
	return updateMu();
    }
    
    template<typename Tf, typename Tr> inline
    double mer<Tf,Tr>::setBetaU(Rcpp::NumericVector const& bBase,
				Rcpp::NumericVector const& uBase,
				Rcpp::NumericVector const&  incB,
				Rcpp::NumericVector const&  incU,
				double                      step) {
	fe.setBeta(bBase, incB, step);
	re.setU(uBase, incU, step);
	return updateMu();
    }
    
    template<typename Tf, typename Tr> inline
    double mer<Tf,Tr>::setU(Rcpp::NumericVector const& uBase,
			    Rcpp::NumericVector const&  incU,
			    double                      step) {
	re.setU(uBase, incU, step);
	return updateMu();
    }
    
    template<typename Tf, typename Tr> inline
    void mer<Tf,Tr>::solveBetaU() {
	fe.updateRzxRx(re.Lambda(), re.L());
	re.updateU(fe.updateBeta(re.cu()));
    }

    template<typename Tf, typename Tr> inline
    void mer<Tf,Tr>::solveU() {
	re.solveU();
    }

    template<typename Tf, typename Tr> inline
    void mer<Tf,Tr>::updateLambda(Rcpp::NumericVector const& nt) {
	re.updateLambda(nt);
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
	const Rcpp::NumericVector &u = re.u(), &offset = resp.offset();
	Rcpp::NumericVector b(q()), gamma(offset.size());
	std::copy(offset.begin(), offset.end(), gamma.begin());
	MatrixNs::chmDn bb(b), gg(gamma);

	re.Lambda().dmult('N', 1., 0., MatrixNs::chmDn(u), bb);
	re.Zt().dmult('T', 1., 1., bb, gg);
	if (fe.beta().size() > 0)
	    fe.X().dmult('N', 1., 1., MatrixNs::chmDn(fe.beta()), gg);
	return resp.updateMu(gamma) + re.sqrLenU();
    }	

    /**
     * Update the weighted residuals, wrss, sqrtrwt, sqrtXwt, U, Utr,
     * cu, UtV, VtV and Vtr.
     *
     * @return penalized, weighted residual sum of squares
     */
    template<typename Tf, typename Tr> inline
    double mer<Tf,Tr>::updateWts() {
	double ans = resp.updateWts();
	ans += re.sqrLenU();
	re.reweight(resp.sqrtXwt(), resp.wtres());
	fe.reweight(re.Ut(), resp.sqrtXwt(), resp.wtres());
	return ans;
    }

    template<typename Tf, typename Tr> inline
    double mer<Tf,Tr>::IRLS(int verb) {
	Rcpp::NumericVector bBase(p()), incB(p()), muBase(n());
	double crit, step, c0, c1;
				// sqrtrwt and sqrtXwt must be set
	crit = 10. * CM_TOL;
	updateMu();		// using current beta and u
	for (int i = 0; crit >= CM_TOL && i < CM_MAXITER; i++) {
				// store copies of var and beta
	    std::copy(fe.beta().begin(), fe.beta().end(), bBase.begin());
	    std::copy(resp.mu().begin(), resp.mu().end(), muBase.begin());
	    c0 = updateWts();
	    solveBeta();	// calculate increment in fe.d_beta
	    std::copy(fe.beta().begin(), fe.beta().end(), incB.begin());
	    for (c1 = c0, step = 1.; c0 <= c1 && step > CM_SMIN;
		 step /= 2.) {
		c1 = setBeta(bBase, incB, step);
		if (verb > 1) showincr(step, c0, c1, fe.beta(), "beta");
	    }
	    crit = compareVecWt(muBase, resp.mu(), resp.sqrtrwt());
	    if (verb > 1)
		Rprintf("   convergence criterion: %g\n", crit);
	}
	return Laplace();
    } // IRLS

    template<typename Tf, typename Tr> inline
    double mer<Tf,Tr>::PIRLS(int verb) {
	Rcpp::NumericVector incU(q()), muBase(n()), uBase(q());
	double crit, step, c0, c1;

	crit = 10. * CM_TOL;
	zeroU();		// start from a known position
	updateMu();
	for (int i = 0; crit >= CM_TOL && i < CM_MAXITER; i++) {
				// store copies of mu and u
	    std::copy(re.u().begin(), re.u().end(), uBase.begin());
	    std::copy(resp.mu().begin(), resp.mu().end(), muBase.begin());
	    c0 = updateWts();	// sqrtrwt, sqrtXwt, wtres, pwrss
	    solveU();
	    std::copy(re.u().begin(), re.u().end(), incU.begin());
	    for (c1 = c0, step = 1.; c0 <= c1 && step > CM_SMIN;
		 step /= 2.) {
		c1 = setU(uBase, incU, step);
		if (verb > 1) showincr(step, c0, c1, re.u(), "u");
	    }
	    crit = compareVecWt(muBase, resp.mu(), resp.sqrtrwt());
	    if (verb > 1)
		Rprintf("   convergence criterion: %g\n", crit);
	}
	return Laplace();
    } // PIRLS

    template<typename Tf, typename Tr> inline
    double mer<Tf,Tr>::PIRLSBeta(int verb) {
	Rcpp::NumericVector bBase(p()), incB(p()), incU(q()), muBase(n()), uBase(q());
	double crit, step, c0, c1;
	
	crit = 10. * CM_TOL;
	updateMu();
	for (int i = 0; crit >= CM_TOL && i < CM_MAXITER; i++) {
				// store copies of mu, beta and u
	    std::copy(resp.mu().begin(), resp.mu().end(), muBase.begin());
	    std::copy(re.u().begin(), re.u().end(), uBase.begin());
	    std::copy(fe.beta().begin(), fe.beta().end(), bBase.begin());
	    c0 = updateWts();	// sqrtrwt, sqrtXwt, wtres, wrss
	    solveBetaU();
	    std::copy(re.u().begin(), re.u().end(), incU.begin());
	    std::copy(fe.beta().begin(), fe.beta().end(), incB.begin());
	    for (c1 = c0, step = 1.; c0 <= c1 && step > CM_SMIN;
		 step /= 2.) {
		c1 = setBetaU(bBase, uBase, incB, incU, step);
		if (verb > 1) {
		    showincr(step, c0, c1, fe.beta(), "beta");
		    showdbl(re.u().begin(), "u", re.u().size());
		}
	    }
	    crit = compareVecWt(muBase, resp.mu(), resp.sqrtrwt());
	    if (verb > 1)
		Rprintf("   convergence criterion: %g\n", crit);
	}
	return Laplace();
    } // PIRLSBeta
}

#endif
