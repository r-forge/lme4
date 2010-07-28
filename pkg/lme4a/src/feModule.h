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

    public:
	deFeMod(Rcpp::S4 xp,int);

	MatrixNs::Cholesky  const&  RX() const {return d_fac;}
	MatrixNs::dgeMatrix const& RZX() const {return d_RZX;}

	Rcpp::NumericVector updateBeta(Rcpp::NumericVector const&);

	void updateRzxRx(MatrixNs::chmSp const&,
			 MatrixNs::chmFr const&);
	void updateUtV(   cholmod_sparse const*);
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
