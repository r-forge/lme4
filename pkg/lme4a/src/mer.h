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
    void showCHM_SP(const MatrixNs::chmSp&, std::string const&);
    void showdMat(MatrixNs::ddenseMatrix const&, std::string const&);
    void showdbl(const double*, const char*, int);
    void showdbl(const Rcpp::NumericVector& vv, const char* nm);
    void showincr(double,double,double,
		  Rcpp::NumericVector const&, const char*);
    void showint(const int*, const char*, int);

    struct lengthFun : std::unary_function<Rcpp::RObject, R_len_t> {
	inline R_len_t operator() (Rcpp::RObject const& x) {return Rf_length(SEXP(x));}
    };

    struct nlevsFun : std::unary_function<Rcpp::RObject, R_len_t> {
	inline R_len_t operator() (Rcpp::RObject const& x) {
	    return Rf_length(Rf_getAttrib(SEXP(x), R_LevelsSymbol));
	}
    };

    struct sqrtquotFun : std::binary_function<double,double,double> {
	inline double operator() (double x, double y) {
	    return sqrt(x / y);
	}
    };

    struct sqrtFun : std::unary_function<double,double> {
	inline double operator() (double x) {
	    return sqrt(x);
	}
    };

    Rcpp::NumericVector mkans(double, const Rcpp::NumericVector&, const Rcpp::NumericVector&);

    class reModule {
    protected:
	Rcpp::S4            d_xp;
	MatrixNs::chmFr     d_L;
	MatrixNs::chmSp     d_Lambda, d_Zt;
	Rcpp::IntegerVector d_Lind;
	Rcpp::NumericVector d_lower, d_u, d_cu;
	double              d_ldL2, d_sqrLenU;
	CHM_SP              d_Ut;
    public:
	reModule(Rcpp::S4);

	const Rcpp::NumericVector    &cu() const {return d_cu;}
 	const Rcpp::NumericVector     &u() const {return d_u;}
	const Rcpp::S4               &xp() const {return d_xp;}
	const Rcpp::NumericVector &lower() const {return d_lower;}
	const MatrixNs::chmFr         &L() const {return d_L;}
	const MatrixNs::chmSp    &Lambda() const {return d_Lambda;}
	const cholmod_sparse         *Ut() const {return d_Ut;}
	const MatrixNs::chmSp        &Zt() const {return d_Zt;}
	double                      ldL2() const {return d_ldL2;}
	double                   sqrLenU() const {return d_sqrLenU;}

	void reweight(Rcpp::NumericMatrix const&,
		      Rcpp::NumericVector const&);
	void setU(Rcpp::NumericVector const&,
		  Rcpp::NumericVector const& = Rcpp::NumericVector(),
		  double = 0.);
	void solveU();
//	void updateDcmp(Rcpp::List&) const;  //needs Matrix_0.999375-42 or later
	void updateDcmp(Rcpp::List&);
	void updateLambda(Rcpp::NumericVector const&);
	void updateU(Rcpp::NumericVector const&);
	void zeroU();
    };

    class reTrms : public reModule {
	Rcpp::List          d_flist, d_cnms;
	Rcpp::IntegerVector d_assign;
    public:
	reTrms(Rcpp::S4);

	Rcpp::IntegerVector const& assign() const {return d_assign;}
	Rcpp::IntegerVector         nlevs() const; // number of levels per factor
	Rcpp::IntegerVector         ncols() const; // number of columns per term
	Rcpp::IntegerVector         nctot() const; // total number of columns per factor
	Rcpp::IntegerVector       offsets() const; // offsets into b vector for each term
	Rcpp::IntegerVector         terms(int) const; // 0-based indices of terms for a factor
	Rcpp::List                condVar(double);
	Rcpp::List const&            cnms() const {return d_cnms;}
	Rcpp::List const&           flist() const {return d_flist;}
    };
	
    class feModule {
    protected:
	Rcpp::NumericVector d_beta, d_Vtr;
	double              d_ldRX2;
    public:
	feModule(Rcpp::S4 xp);

	void setBeta(Rcpp::NumericVector const&,
		     Rcpp::NumericVector const& = Rcpp::NumericVector(),
		     double = 0.);
// Don't want to use the name beta because the R include files remap it to Rf_beta
	const Rcpp::NumericVector& getBeta() const {return  d_beta;}
	double                       ldRX2() const {return  d_ldRX2;}
	virtual void updateDcmp(Rcpp::List&) = 0;
    };

    class deFeMod : public feModule {
    protected:
	MatrixNs::dgeMatrix d_RZX, d_X;
	MatrixNs::Cholesky        d_RX;
	MatrixNs::dgeMatrix d_UtV, d_V;
	MatrixNs::dpoMatrix      d_VtV;

    public:
	deFeMod(Rcpp::S4 xp,int);

	const MatrixNs::Cholesky   &RX() const{return d_RX;}
	const MatrixNs::dgeMatrix   &X() const{return d_X;}
	const MatrixNs::dgeMatrix &RZX() const{return d_RZX;}

	Rcpp::NumericVector updateBeta(Rcpp::NumericVector const&);

	void reweight(cholmod_sparse      const*,
		      Rcpp::NumericMatrix const&,
		      Rcpp::NumericVector const&);
	void solveBeta();
	void updateDcmp(Rcpp::List&) ;
	void updateRzxRx(MatrixNs::chmSp const&,
			 MatrixNs::chmFr const&);
    };

    class spFeMod : public feModule {
    protected:
	MatrixNs::chmSp    d_RZX, d_X;
	MatrixNs::chmFr    d_RX;
	CHM_SP             d_UtV, d_V, d_VtV;
    public:
	spFeMod(Rcpp::S4 xp,int);
	~spFeMod();

	const MatrixNs::chmSp      &X() const{return d_X;}
	const MatrixNs::chmSp    &RZX() const{return d_RZX;}
	const MatrixNs::chmFr     &RX() const{return d_RX;}

	Rcpp::NumericVector updateBeta(Rcpp::NumericVector const&);

	void reweight(cholmod_sparse      const*,
		      Rcpp::NumericMatrix const&,
		      Rcpp::NumericVector const&);
	void solveBeta();
//	void updateDcmp(Rcpp::List&) const; // needs Matrix_0.999375-42 or later
	void updateDcmp(Rcpp::List&);
	void updateRzxRx(MatrixNs::chmSp const&,
			 MatrixNs::chmFr const&);
    };

    class merResp {
    protected:
	Rcpp::S4             d_xp;
	double               d_wrss;
	Rcpp::NumericVector  d_offset, d_weights, d_y; // read-only slots
	Rcpp::NumericVector  d_sqrtrwt, d_wtres, d_mu; // writeable slots
	Rcpp::NumericMatrix  d_sqrtXwt;
    public:
	merResp(Rcpp::S4);

	const Rcpp::NumericVector      &mu() const{return d_mu;}
	const Rcpp::NumericVector  &offset() const{return d_offset;}
	const Rcpp::NumericMatrix &sqrtXwt() const{return d_sqrtXwt;}
	const Rcpp::NumericVector &sqrtrwt() const{return d_sqrtrwt;}
	const Rcpp::NumericVector   &wtres() const{return d_wtres;}
	const Rcpp::S4                 &xp() const{return d_xp;}

	double                        wrss() const{return d_wrss;}
	double                   updateWts(){return updateWrss();}
	double                  updateWrss();

	virtual void            updateDcmp(Rcpp::List&) const = 0;
    };

    class lmerResp : public merResp {
	int d_reml;
    public:
	lmerResp(Rcpp::S4);
	double Laplace(double,double,double)const;
	double updateMu(Rcpp::NumericVector const&);
	void updateDcmp(Rcpp::List&) const;
    };

    class glmerResp : public merResp {
    protected:
	glm::glmFamily           family;
	Rcpp::NumericVector  d_eta, d_n;
    public:
	glmerResp(Rcpp::S4 xp);

	Rcpp::NumericVector   devResid() const;

	const Rcpp::NumericVector &eta() const{return d_eta;}

	double                 Laplace(double,double,double) const;
	double                updateMu(Rcpp::NumericVector const&);
	double               updateWts();

	void  updateDcmp(Rcpp::List&) const;
    };
    
    class nlmerResp : public merResp {
	Rcpp::Environment nlenv;
	Rcpp::Language nlmod;
	Rcpp::CharacterVector pnames;
    public:
	nlmerResp(Rcpp::S4 xp);
	double updateMu(Rcpp::NumericVector const &gamma);
	double Laplace(double,double,double) const;
	void  updateDcmp(Rcpp::List&) const;
    };
    
    enum Alg {Beta, U, BetaU};

    /* Model object template
     * Tf is the type of fixed-effects module (deFeMod or spFeMod)
     * Tr is the type of response module (merResp, glmerResp, nlmerResp or nglmerResp)
     */
    template<typename Tf, typename Tr>  
    class mer {
	reModule re;
	Tr resp;
	Tf fe;
    public:
	mer(Rcpp::S4&);

	double Laplace  () const {
	    return resp.Laplace(re.ldL2(), fe.ldRX2(), re.sqrLenU());
	}
	Rcpp::NumericVector LMMdeviance(const Rcpp::NumericVector&,
					const Rcpp::NumericVector&);
	double PIRLS      (Rcpp::NumericVector const&,int,Alg);
	double setBetaU   (Rcpp::NumericVector const&,
			   Rcpp::NumericVector const&,
			   Rcpp::NumericVector const&,
			   Rcpp::NumericVector const&,
			   double,Alg);
	double updateMu   ();
	double updateWts  ();

	int N()  const   {return resp.offset().size();}
	int n()  const   {return  resp.wtres().size();}
	int nth()const   {return    re.lower().size();}
	int p()  const   {return  fe.getBeta().size();}
	int q()  const   {return        re.u().size();}
	int s()  const   {return              N()/n();}

	void solveCoef(Alg);
	void updateRzxRx();
	void updateDcmp(Rcpp::List&,Rcpp::NumericVector const&);
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
    Rcpp::NumericVector mer<Tf,Tr>::LMMdeviance(const Rcpp::NumericVector& nt, const Rcpp::NumericVector& u0) {
	re.updateLambda(nt);
	re.setU(u0);
	updateWts();
	solveCoef(BetaU);
	updateMu();
	return mkans(Laplace(), fe.getBeta(), re.u());
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
	if (alg !=    U) fe.setBeta(bBase, incB, step);
	if (alg != Beta) re.setU(uBase, incU, step);
	return updateMu();
    }
    
    template<typename Tf, typename Tr> inline
    void mer<Tf,Tr>::solveCoef(Alg alg) {
	switch(alg) {
	case Beta:
	    fe.solveBeta();
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
	const Rcpp::NumericVector &u = re.u(), &offset = resp.offset();
	Rcpp::NumericVector b(q()), gamma(offset.size());
	std::copy(offset.begin(), offset.end(), gamma.begin());
	MatrixNs::chmDn bb(b), gg(gamma);

	re.Lambda().dmult('N', 1., 0., MatrixNs::chmDn(u), bb);
	re.Zt().dmult('T', 1., 1., bb, gg);
	if (fe.getBeta().size() > 0)
	    fe.X().dmult('N', 1., 1., MatrixNs::chmDn(fe.getBeta()), gg);
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
	fe.reweight(re.Ut(), resp.sqrtXwt(), resp.wtres());
	fe.updateRzxRx(re.Lambda(), re.L());
    }

    /**
     * Update the weighted residuals, wrss, sqrtrwt, sqrtXwt, U,
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

#define CM_TOL 1.e-4
#define CM_MAXITER 200
#define CM_SMIN 1.e-4

    template<typename Tf, typename Tr> inline
    double mer<Tf,Tr>::PIRLS(Rcpp::NumericVector const& u0, int verb, Alg alg) {
	Rcpp::NumericVector bBase(p()), incB(p()), incU(q()), muBase(n()), uBase(q());
	double crit, step, c0, c1;
				// sqrtrwt and sqrtXwt must be set
	crit = 10. * CM_TOL;
	if (u0.size() != re.u().size())
	    throw std::runtime_error("lengths of u0 and re.u() must match");
	re.setU(u0);
	updateMu();		// using current beta and u
	for (int i = 0; crit >= CM_TOL && i < CM_MAXITER; i++) {
				// store copies of mu, u and beta
	    std::copy(resp.mu().begin(), resp.mu().end(), muBase.begin());
	    std::copy(re.u().begin(), re.u().end(), uBase.begin());
	    std::copy(fe.getBeta().begin(), fe.getBeta().end(), bBase.begin());
	    c0 = updateWts();
	    solveCoef(alg);
	    std::copy(fe.getBeta().begin(), fe.getBeta().end(), incB.begin());
	    std::copy(re.u().begin(), re.u().end(), incU.begin());
	    for (c1 = c0, step = 1.; c0 <= c1 && step > CM_SMIN;
		 step /= 2.) {
		c1 = setBetaU(bBase, uBase, incB, incU, step, alg);
		if (verb > 1) {
		    showincr(step, c0, c1, fe.getBeta(), "beta");
		    showdbl(re.u().begin(), "u", re.u().size());
		}
	    }
	    crit = compareVecWt(muBase, resp.mu(), resp.sqrtrwt());
	    if (verb > 1)
		Rprintf("   convergence criterion: %g\n", crit);
	}
	return Laplace();
    } // PIRLS

    template<typename Tf, typename Tr> inline
    void mer<Tf,Tr>::updateDcmp(Rcpp::List& ans, Rcpp::NumericVector const& pars) {
	Rcpp::NumericVector th(nth()), beta(p()), u(q());
	double *pp = pars.begin();
	std::copy(pp, pp + nth(), th.begin());
	std::copy(pp + nth(), pp + nth() + p(), beta.begin());
	std::copy(pp + nth() + p(), pp + nth() + p() + q(), u.begin());
	LMMdeviance(pars, u);
	fe.setBeta(beta);
	(&resp)->updateDcmp(ans);
	(&re)->updateDcmp(ans);
	(&fe)->updateDcmp(ans);
    }
}

#endif
