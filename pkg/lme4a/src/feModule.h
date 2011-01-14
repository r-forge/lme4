// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
#ifndef LME4_FEMODULE_H
#define LME4_FEMODULE_H

#include <Rcpp.h>
#include "MatrixNs.h"
#include "predModules.h"

namespace mer {
    class feModule {
    protected:
	double   d_ldRX2;
    public:
	double     ldRX2() const {return  d_ldRX2;}
    };

    class deFeMod : public feModule, public matMod::dPredModule {
    protected:
	MatrixNs::dgeMatrix      d_RZX;
	MatrixNs::dgeMatrix      d_UtV;
	MatrixNs::dpoMatrix      d_VtV;
	Rcpp::NumericVector d_coef0, d_incr;

    public:
	deFeMod(Rcpp::S4,int);
	deFeMod(Rcpp::S4,int,int,int);

	Rcpp::NumericVector linPred() const;
	Rcpp::NumericVector linPred1(double) const;
	Rcpp::NumericVector updateBeta(Rcpp::NumericVector const&);
	Rcpp::NumericVector updateIncr(Rcpp::NumericVector const&);

	void updateRzxRxp(Rcpp::S4, Rcpp::S4);
	void updateRzxpRxpp(Rcpp::XPtr<MatrixNs::chmSp>,
			    Rcpp::XPtr<MatrixNs::chmFr>);
	void updateRzxRx(MatrixNs::chmSp const&, MatrixNs::chmFr const&);
	void updateUtV(   cholmod_sparse const*);
	void updateUtVp(Rcpp::XPtr<cholmod_sparse>);
	void reweight(Rcpp::NumericMatrix const&,
		      Rcpp::NumericVector const&) throw(std::runtime_error);
	void installCoef0();
	void setCoef0(const Rcpp::NumericVector&) throw (std::runtime_error);
	void setIncr(const Rcpp::NumericVector&) throw (std::runtime_error);
	
	double                          ldRX2() const {return  d_ldRX2;}
	const Rcpp::NumericVector&       coef() const {return d_coef;}
	const Rcpp::NumericVector&       incr() const {return d_incr;}
	const Rcpp::NumericVector&      coef0() const {return d_coef0;}
	const Rcpp::NumericVector&        Vtr() const {return d_Vtr;}
	const MatrixNs::Cholesky&          RX() const {return d_fac;}
	const MatrixNs::dgeMatrix&        RZX() const {return d_RZX;}
	const MatrixNs::dgeMatrix&        UtV() const {return d_UtV;}
	const MatrixNs::dpoMatrix&        VtV() const {return d_VtV;}
	const MatrixNs::ddenseModelMatrix&  X() const {return d_X;}
	const MatrixNs::dgeMatrix&          V() const {return d_V;}
    };

    class spFeMod : public feModule, public matMod::sPredModule {
    protected:
	MatrixNs::chmSp    d_RZX;
	CHM_SP             d_UtV, d_VtV;
    public:
	spFeMod(Rcpp::S4 xp,int);
	~spFeMod();

	MatrixNs::chmFr const&     RX() const {return d_fac;}
	MatrixNs::chmSp const&    RZX() const {return d_RZX;}

	Rcpp::NumericVector updateBeta(Rcpp::NumericVector const&);

	void updateRzxRx(const MatrixNs::chmSp&,
			 const MatrixNs::chmFr&);
	void updateUtV(  cholmod_sparse const*);
    };

}
#endif
