// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
#ifndef LME4_PREDMODULES_H
#define LME4_PREDMODULES_H

#include <Rcpp.h>
#include "MatrixNs.h"

namespace matMod {	   // matMod is the namespace for matrixModels
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
    protected:
	MatrixNs::ddenseModelMatrix   d_X;
	MatrixNs::dgeMatrix           d_V;
	MatrixNs::Cholesky          d_fac;
    public:
	dPredModule(Rcpp::S4,int);
	
	MatrixNs::ddenseModelMatrix const& X() const {return d_X;}
	Rcpp::NumericVector linPred() const;

	void reweight(Rcpp::NumericMatrix const&,
		      Rcpp::NumericVector const&) throw(std::runtime_error);
	double solveCoef(double = 1.);
    };
    
    class sPredModule : public predModule {
    protected:
	MatrixNs::chmSp          d_X;
	MatrixNs::chmFr        d_fac;
	CHM_SP                   d_V;
    public:
	sPredModule(Rcpp::S4,int);
	
	MatrixNs::chmSp      const& X() const {return d_X;}
	const_CHM_SP                V() const {return d_V;}

	Rcpp::NumericVector linPred() const;
	void reweight(Rcpp::NumericMatrix const&,
		      Rcpp::NumericVector const&) throw(std::runtime_error);
	double solveCoef(double = 1.);
    };
    
}

#endif /* LME4_GLM_H */

