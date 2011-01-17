// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
#ifndef LME4_REMODULE_H
#define LME4_REMODULE_H

#include <Rcpp.h>
#include "MatrixNs.h"

namespace mer {
// Create polymorphic functions updateRzxpRxpp in deFeMod to
// accomodate this class.
    class ssre {
	MatrixNs::chmSp     d_Z;
	Rcpp::IntegerVector d_ff;
	Rcpp::NumericVector d_L, d_u, d_cu;
	double              d_ldL2;
	CHM_SP              d_Ut;
    public:
	ssre(Rcpp::S4);

	const Rcpp::NumericVector    &cu() const {return d_cu;}
 	const Rcpp::NumericVector     &u() const {return d_u;}
	const cholmod_sparse         *Ut() const {return d_Ut;}
	const MatrixNs::chmSp         &Z() const {return d_Z;}
	double                      ldL2() const {return d_ldL2;}
	// used in PIRLS when beta is not part of the solution
	double                    solveU();
	double                   sqrLenU() const {return sum(d_u * d_u);}

	void reweight    (const Rcpp::NumericMatrix&,
			  const Rcpp::NumericVector&);
	void setU        (const Rcpp::NumericVector&,
			  const Rcpp::NumericVector& = Rcpp::NumericVector(),
			  double = 0.);
	void setTheta(const Rcpp::NumericVector&);
	// used in PIRLSbeta for solving for both u and beta
	void updateU     (const Rcpp::NumericVector&);
    };
	    
    class reModule {
    protected:
	MatrixNs::chmFr     d_L;
	MatrixNs::chmSp     d_Lambda, d_Zt;
	Rcpp::IntegerVector d_Lind;
	Rcpp::NumericVector d_lower, d_theta, d_u0, d_incr, d_u, d_cu;
	double              d_ldL2, d_CcNumer; // numerator of conv. crit. 
	CHM_SP              d_Ut;
    public:
	reModule(Rcpp::S4);
	reModule(Rcpp::S4, Rcpp::S4, Rcpp::S4,
		 Rcpp::IntegerVector, Rcpp::NumericVector)
	    throw(MatrixNs::wrongS4);

	cholmod_sparse      const*     Ut() const {return d_Ut;}

	double                       ldL2() const {return d_ldL2;}
	// used in PIRLS when beta is not part of the solution
	double                     solveU();
	double                  solveIncr();
	double                    sqrLenU() const {return sum(d_u * d_u);}
	double                    CcNumer() const {return d_CcNumer;}

	MatrixNs::chmFr     const&      L() const {return d_L;}
	MatrixNs::chmSp     const& Lambda() const {return d_Lambda;}
	MatrixNs::chmSp     const&     Zt() const {return d_Zt;}
	Rcpp::XPtr<MatrixNs::chmSp>Lambdap() {
	    Rcpp::XPtr<MatrixNs::chmSp> pt(&d_Lambda, false);
	    return pt;
	}
	Rcpp::XPtr<MatrixNs::chmFr>    Lp() {
	    Rcpp::XPtr<MatrixNs::chmFr> pt(&d_L, false);
	    return pt;
	}

	Rcpp::IntegerVector const&   Lind() const {return d_Lind;}

	Rcpp::NumericVector const&     cu() const {return d_cu;}
	Rcpp::NumericVector const&  lower() const {return d_lower;}
	Rcpp::NumericVector const&  theta() const {return d_theta;}
 	Rcpp::NumericVector const&      u() const {return d_u;}
 	Rcpp::NumericVector const&     u0() const {return d_u0;}
 	Rcpp::NumericVector const&   incr() const {return d_incr;}
 	Rcpp::NumericVector             b() const;

	Rcpp::NumericVector       linPred() const;
	Rcpp::NumericVector       linPred1(double);
	Rcpp::XPtr<cholmod_sparse>    Utp() const;

	void reweight    (const Rcpp::NumericMatrix&,
			  const Rcpp::NumericVector&, bool useU0=false);
	void installU0   ();
	void setU0       (const Rcpp::NumericVector&) throw (std::runtime_error);
	void setU        (const Rcpp::NumericVector&,
			  const Rcpp::NumericVector& = Rcpp::NumericVector(),
			  double = 0.) throw (std::runtime_error);
	void setTheta    (const Rcpp::NumericVector&)
	    throw (std::runtime_error);
	// used in PIRLSbeta to solve for u only
	void updateU     (const Rcpp::NumericVector&);
	void updateIncr  (const Rcpp::NumericVector&);
    };

    class reTrms : public reModule {
	Rcpp::List          d_flist, d_cnms;
	Rcpp::IntegerVector d_assign;
    public:
	reTrms(Rcpp::S4);

	const Rcpp::IntegerVector& assign() const {return d_assign;}
	const Rcpp::List&            cnms() const {return d_cnms;}
	const Rcpp::List&           flist() const {return d_flist;}

	Rcpp::IntegerVector         nlevs() const; // number of levels per factor
	Rcpp::IntegerVector         ncols() const; // number of columns per term
	Rcpp::IntegerVector         nctot() const; // total number of columns per factor
	Rcpp::IntegerVector       offsets() const; // offsets into b vector for each term
	Rcpp::IntegerVector         terms(R_len_t) const; // 0-based indices of terms for a factor
	Rcpp::List                condVar(double);
    };
}

#endif
