// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
#ifndef LME4_MER_H
#define LME4_MER_H

#include <Rcpp.h>
#include "MatrixNs.h"
#include "glmFamily.h"

namespace mer {

    double compare_vec(Rcpp::NumericVector const&,Rcpp::NumericVector const&);
    void showCHM_DN(const_CHM_DN, std::string const&);
    void showCHM_FR(const_CHM_FR, std::string const&);
    void showCHM_SP(const_CHM_SP, std::string const&);
    void showdbl(const double*, const char*, int);
    void showincr(double,double,double, Rcpp::NumericVector const&, const char*);
    void showint(const int*, const char*, int);

    class reModule {
	MatrixNs::chmFr     d_L;
	MatrixNs::chmSp     d_Lambda, d_Ut, d_Zt;
	Rcpp::IntegerVector Lind;
	Rcpp::NumericVector lower, theta;
	Rcpp::NumericVector d_Utr, d_cu, d_u;
	double             *d_ldL2, d_sqlLenU;
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
	double                 sqrLenU() const {return  d_sqlLenU;}

	Rcpp::NumericVector rwUIncr(Rcpp::NumericMatrix const&,
				    Rcpp::NumericVector const&);
	void setU(Rcpp::NumericVector const&,
		  Rcpp::NumericVector const& = Rcpp::NumericVector(),
		  double = 0.);
	void updateTheta(Rcpp::NumericVector const&);
	void updateU(MatrixNs::chmDn const&);
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

	void updateRzxRx(MatrixNs::chmSp const&,
			 MatrixNs::chmFr const&);
	void updateUtV(MatrixNs::chmSp const&);
	void updateV(Rcpp::NumericMatrix const&);
	void updateVtr(Rcpp::NumericVector const&);
	void updateRX(bool);
	void solveA(Rcpp::NumericVector&) const;

	Rcpp::NumericVector rwBetaIncr(Rcpp::NumericMatrix const&,
				       Rcpp::NumericVector const&);
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
	Rcpp::NumericVector rwBetaIncr(Rcpp::NumericMatrix const&,
				       Rcpp::NumericVector const&);
	Rcpp::NumericVector updateBeta(Rcpp::NumericVector const&);
    };

    class merResp {
    protected:
	double              *d_wrss;
	Rcpp::NumericVector  d_offset, d_sqrtrwt, d_wtres, d_mu, d_weights, d_y;
	Rcpp::NumericMatrix  d_sqrtXwt;
    public:
	merResp(Rcpp::S4);

	const Rcpp::NumericMatrix &sqrtXwt() const{return d_sqrtXwt;}
	const Rcpp::NumericVector &sqrtrwt() const{return d_sqrtrwt;}
	const Rcpp::NumericVector   &wtres() const{return d_wtres;}
	const Rcpp::NumericVector  &offset() const{return d_offset;}
	double                        wrss() const{return *d_wrss;}

	double updateWrss();
	void updateL(MatrixNs::chmSp const&, MatrixNs::chmFr);
	void updateWts(){}	// overridden by glmerResp
    };

    class lmerResp : public merResp {
	bool d_reml;
    public:
	lmerResp(Rcpp::S4);
	double Laplace(double,double,double)const;
	double updateMu(Rcpp::NumericVector const&);
    };

    class glmerResp : public merResp {
    protected:
	glmFamily family;
	Rcpp::NumericVector d_eta, d_muEta, d_n, d_var;
    public:
	glmerResp(Rcpp::S4 xp);

	const Rcpp::NumericVector &var() const {return d_var;}
	const Rcpp::NumericVector &eta() const {return d_var;}
	double Laplace(double,double,double) const;
	double updateMu(Rcpp::NumericVector const &gamma);
	Rcpp::NumericVector devResid();
	void updateWts();
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

	double IRLS(int);
	double Laplace() const {
	    return resp.Laplace(re.ldL2(), fe.ldRX2(), re.sqrLenU());
	}
	double PIRLS(int);
	double setBeta(Rcpp::NumericVector const&,
		       Rcpp::NumericVector const&, double);
	double setU(Rcpp::NumericVector const&,
		    Rcpp::NumericVector const&, double);
	double updateMu();

	int N()const{return resp.offset().size();}
	int n()const{return  resp.wtres().size();}
	int p()const{return     fe.beta().size();}
	int q()const{return        re.u().size();}

	void solveBeta();
	void solveBetaU();
	void solveU();
	void updateTheta(const Rcpp::NumericVector&nt){re.updateTheta(nt);}
	void updateRzxRx(){fe.updateRzxRx(re.Lambda(), re.L());}
    };
    
    template<typename Tf, typename Tr>
    inline mer<Tf,Tr>::mer(Rcpp::S4 xp)
	: re(Rcpp::S4(xp.slot("re"))),
	  fe(Rcpp::S4(xp.slot("fe"))),
	  resp(Rcpp::S4(xp.slot("resp"))) {
    }

#define CM_TOL 1.e-4
#define CM_MAXITER 30
#define CM_SMIN 1.e-4

/** 
 * Update the conditional mean, weighted residuals and wrss in resp
 * from the linear predictor.  Some resp modules also update other
 * information, such as the sqrtXwt matrix. 
 *
 * @return the *penalized*, weighted residual sum of squares
 */
    template<typename Tf, typename Tr> inline
    double mer<Tf,Tr>::updateMu() {
	const Rcpp::NumericVector u = re.u(), offset = resp.offset();
	Rcpp::NumericVector b(q()), gamma(offset.size());
	std::copy(offset.begin(), offset.end(), gamma.begin());
	MatrixNs::chmDn bb(b), gg(gamma);

	re.Lambda().dmult('N', 1., 0., MatrixNs::chmDn(u), bb);
	re.Zt().dmult('T', 1., 1., bb, gg);
	fe.X().dmult('N', 1., 1., MatrixNs::chmDn(fe.beta()), gg);
	return resp.updateMu(gamma) + re.sqrLenU();
    }	

    template<typename Tf, typename Tr> inline
    double mer<Tf,Tr>::setBeta(Rcpp::NumericVector const&  bb,
			       Rcpp::NumericVector const& inc,
			       double                    step) {
	fe.setBeta(bb, inc, step);
	return updateMu();
    }
    
    template<typename Tf, typename Tr> inline
    double mer<Tf,Tr>::setU(Rcpp::NumericVector const&  uu,
			    Rcpp::NumericVector const& inc,
			    double                    step) {
	re.setU(uu, inc, step);
	return updateMu();
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
    double mer<Tf,Tr>::IRLS(int verb) {
	const Rcpp::NumericVector &var = resp.var(), &beta = fe.beta();
	Rcpp::NumericVector betabase(p()), varold(n());
	double crit, step, c0, c1;
				// sqrtrwt and sqrtXwt must be set
	crit = 10. * CM_TOL;
	for (int i = 0; crit >= CM_TOL && i < CM_MAXITER; i++) {
				// store copies of var and beta
	    std::copy(beta.begin(), beta.end(), betabase.begin());
	    std::copy(var.begin(), var.end(), varold.begin());
	    c0 = updateMu();	// using current beta and u
	    Rcpp::NumericVector
		incr = fe.rwBetaIncr(resp.sqrtXwt(), resp.wtres());
	    for (c1 = c0, step = 1.; c0 <= c1 && step > CM_SMIN;
		 step /= 2.) {
		c1 = setBeta(betabase, incr, step);
		if (verb > 1) showincr(step, c0, c1, beta, "beta");
	    }
	    resp.updateWts(); // sqrtrwt, muEta and sqrtXwt
	    crit = compare_vec(varold, var);
	    if (verb > 1) Rprintf("convergence criterion: %g\n", crit);
	}
	return Laplace();
    } // IRLS

    template<typename Tf, typename Tr> inline
    double mer<Tf,Tr>::PIRLS(int verb) {
	const Rcpp::NumericVector u = re.u(), var = resp.var();
	Rcpp::NumericVector ubase(u.size()), varold(var.size());
	double crit, step, c0, c1;

	crit = 10. * CM_TOL;
	for (int i = 0; crit >= CM_TOL && i < CM_MAXITER; i++) {
				// store copies of var and u
	    std::copy(u.begin(), u.end(), ubase.begin());
	    std::copy(var.begin(), var.end(), varold.begin());
	    c0 = updateMu();
	    Rcpp::NumericVector
		incr = re.rwUIncr(resp.sqrtXwt(), resp.wtres());
	    for (c1 = c0, step = 1.; c0 <= c1 && step > CM_SMIN;
		 step /= 2.) {
		c1 = setU(ubase, incr, step);
		if (verb > 1) showincr(step, c0, c1, u, "u");
	    }
	    resp.updateWts();	// sqrtrwt, muEta and sqrtXwt
	    crit = compare_vec(varold, var);
	    if (verb > 1) Rprintf("convergence criterion: %g\n", crit);
	}
	return Laplace();
    } // PIRLS
}

#endif
