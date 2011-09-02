//
// predModule.cpp: implementation of predictor module using Eigen
//
// Copyright (C)       2011 Douglas Bates, Martin Maechler and Ben Bolker
//
// This file is part of lme4Eigen.

#include "predModule.h"

namespace lme4Eigen {
    using std::copy;
    using std::cout;
    using std::endl;
    using std::fill;

#if 0
    void cholmod_dump_common(cholmod_common& c) {
	cout << "Cholmod common structure" << endl;
	cout << "status = " << c.status
	     << ", dbound = " << c.dbound
	     << ", grow0 = " << c.grow0
	     << ", grow1 = " << c.grow1
	     << ", grow2 = " << c.grow2
	     << ", maxrank = " << c.maxrank << endl;
	cout << "supernodal_switch = " << c.supernodal_switch
	     << ", final_asis = " << c.final_asis
	     << ", final_super = " << c.final_super
	     << ", final_ll = " << c.final_ll
	     << ", final_pack = " << c.final_pack
	     << ", final_monotonic = " << c.final_monotonic
	     << ", final_resymbol = " << c.final_resymbol << endl;
	cout << "prefer_upper = " << c.prefer_upper
	     << ", print = " << c.print
	     << ", precise = " << c.precise << endl;
	cout << "nmethods = " << c.nmethods
	     << ", current = " << c.current
	     << ", selected = " << c.selected
	     << ", postorder = " << c.postorder << endl;
	cout << "method numbers: " << c.method[0].ordering;
	for (int i = 1; i < c.nmethods; ++i) cout << ", " << c.method[i].ordering;
	cout << endl;
    }
#endif

    merPredD::merPredD(S4 X, S4 Zt, S4 Lambdat, IntegerVector Lind,
		       NumericVector theta)
	: d_X(X),
          d_Zt(Zt),
	  d_theta(clone(theta)),
	  d_Lind(Lind),
	  d_n(d_X.rows()),
	  d_p(d_X.cols()),
	  d_q(d_Zt.rows()),
	  d_nnz(-1),
	  d_Lambdat(Lambdat),
	  d_RZX(d_q, d_p),
	  d_V(d_X),
	  d_VtV(d_p,d_p),
	  d_Vtr(d_p),
	  d_Utr(d_q),
	  d_delb(d_p),
	  d_delu(d_q),
	  d_beta0(d_p),
	  d_u0(d_q),
	  d_Ut(d_Zt),
	  d_RX(d_p)
    {				// Check consistency of dimensions
	if (d_n != d_Zt.cols())
	    throw invalid_argument("Z dimension mismatch");
	if (d_Lind.size() != d_Lambdat.nonZeros())
	    throw invalid_argument("size of Lind does not match nonzeros in Lambda");
	// checking of the range of Lind is now done in R code for reference class
				// initialize beta0, u0, delb, delu and VtV
	d_beta0.setZero(); d_u0.setZero(); d_delu.setZero(); d_delb.setZero();
	d_VtV.setZero().selfadjointView<Eigen::Upper>().rankUpdate(d_V.adjoint());

	setTheta(d_theta);		// starting values into Lambda
        d_L.cholmod().final_ll = 1;	// force an LL' decomposition
// FIXME: updating d_Ut should be done in a method to accomodate the nlmer case
	d_LamtUt = d_Lambdat * d_Ut;
	d_L.analyzePattern(d_LamtUt); // perform symbolic analysis
        if (d_L.info() != Eigen::Success)
	    throw runtime_error("CholeskyDecomposition.analyzePattern failed");
    }

    void merPredD::updateLamtUt() {
	fill(d_LamtUt._valuePtr(), d_LamtUt._valuePtr() + d_LamtUt.nonZeros(), Scalar());
	for (Index j = 0; j < d_Ut.cols(); ++j) {
	    for (SpMatrixd::InnerIterator rhsIt(d_Ut, j); rhsIt; ++rhsIt) {
		SpMatrixd::InnerIterator prdIt(d_LamtUt, j);
		Scalar                       y = rhsIt.value();
		for (SpMatrixd::InnerIterator lhsIt(d_Lambdat, rhsIt.index()); lhsIt; ++lhsIt) {
		    Index      i = lhsIt.index();
		    while (prdIt.index() != i && prdIt) ++prdIt;
		    if (!prdIt) throw runtime_error("logic error in updateLamtUt");
		    prdIt.valueRef() += lhsIt.value() * y;
		}
	    }
	}
    }

    VectorXd merPredD::b(const double& f) const {return d_Lambdat.adjoint() * u(f);}

    VectorXd merPredD::beta(const double& f) const {return d_beta0 + f * d_delb;}

    VectorXd merPredD::linPred(const double& f) const {
	return d_X * beta(f) + d_Zt.adjoint() * b(f);
    }

    VectorXd merPredD::u(const double& f) const {return d_u0 + f * d_delu;}

    merPredD::Scalar merPredD::sqrL(const double& f) const {return u(f).squaredNorm();}

    void merPredD::updateL() {
	updateLamtUt();
	d_L.factorize_p(d_LamtUt, ArrayXi(), 1.);
	d_ldL2 = ::M_chm_factor_ldetL2(d_L.factor());
    }

    void merPredD::setTheta(const NumericVector& theta) {
	// cout << theta[0];
	// for (int i = 1; i < theta.size(); ++i) cout << ", " << theta[i];
	// cout << endl;
	if (theta.size() != d_theta.size())
	    throw invalid_argument("theta size mismatch");
				// update theta
	copy(theta.begin(), theta.end(), d_theta.begin());
				// update Lambdat
	int    *lipt = d_Lind.begin();
	double *LamX = d_Lambdat._valuePtr(), *thpt = d_theta.begin();
	for (int i = 0; i < d_Lind.size(); ++i) {
	    LamX[i] = thpt[lipt[i] - 1];
	}
    }

    void merPredD::solve() {
	d_delu          = d_Utr;
	d_L.solveInPlace(d_delu, CHOLMOD_P);
	d_L.solveInPlace(d_delu, CHOLMOD_L);	
				// d_delu now contains cu
	d_delb          = d_RX.solve(d_Vtr - d_RZX.adjoint() * d_delu);
	d_delu         -= d_RZX * d_delb;
	d_L.solveInPlace(d_delu, CHOLMOD_Lt);
	d_L.solveInPlace(d_delu, CHOLMOD_Pt);
    }

    void merPredD::updateXwts(const VectorXd& sqrtXwt) {
	if (d_X.rows() != sqrtXwt.size())
	    throw invalid_argument("updateXwts: dimension mismatch");
	if (d_V.rows() == d_X.rows()) {  //FIXME: Generalize this for nlmer
	    DiagonalMatrix<double, Dynamic> W(sqrtXwt.asDiagonal());
	    d_V         = W * d_X;
	    d_Ut        = d_Zt * W;
	} else throw invalid_argument("updateRes: no provision for nlmer yet");
	d_VtV.setZero().selfadjointView<Upper>().rankUpdate(d_V.adjoint());
    }

    void merPredD::updateDecomp() { // update L, RZX and RX
	updateL();
	d_RZX           = d_Lambdat * (d_Ut * d_V);
	d_L.solveInPlace(d_RZX, CHOLMOD_P);
	d_L.solveInPlace(d_RZX, CHOLMOD_L);

	MatrixXd      VtVdown(d_VtV);
	d_RX.compute(VtVdown.selfadjointView<Eigen::Upper>().rankUpdate(d_RZX.adjoint(), -1));
	if (d_RX.info() != Eigen::Success)
	    ::Rf_error("Downdated VtV is not positive definite");
	d_ldRX2         = 2. * d_RX.matrixLLT().diagonal().array().abs().log().sum();
    }

    void merPredD::updateRes(const VectorXd& wtres) {
	if (d_V.rows() != wtres.size())
	    throw invalid_argument("updateRes: dimension mismatch");
	d_Vtr           = d_V.adjoint() * wtres;
	d_Utr           = d_Lambdat * (d_Ut * wtres);
    }

    void merPredD::setBeta0(const VectorXd& nBeta) {
	if (nBeta.size() != d_p) throw invalid_argument("setBeta0: dimension mismatch");
	copy(nBeta.data(), nBeta.data() + d_p, d_beta0.data());
    }

    void merPredD::setU0(const VectorXd& newU0) {
	if (newU0.size() != d_q) throw invalid_argument("setU0: dimension mismatch");
	copy(newU0.data(), newU0.data() + d_q, d_u0.data());
    }

    IntegerVector merPredD::Pvec() const {
	const cholmod_factor* cf = d_L.factor();
	int*                 ppt = (int*)cf->Perm;
	return IntegerVector(ppt, ppt + cf->n);
    }
}

extern "C" {
    using Rcpp::IntegerVector;
    using Rcpp::NumericVector;
    using Rcpp::S4;
    using Rcpp::XPtr;
    using Rcpp::as;
    using Rcpp::wrap;

    SEXP merPredDCreate(SEXP Xs, SEXP Zts, SEXP Lambdats, SEXP Linds, SEXP thetas) {
	BEGIN_RCPP;
	S4 X(Xs), Zt(Zts), Lambdat(Lambdats);
	IntegerVector Lind(Linds);
	NumericVector theta(thetas);
	lme4Eigen::merPredD *ans = new lme4Eigen::merPredD(X, Zt, Lambdat, Lind, theta);
	return wrap(XPtr<lme4Eigen::merPredD>(ans, true));
	END_RCPP;
    }

    SEXP merPredDsetBeta0(SEXP ptr, SEXP beta0) {
	BEGIN_RCPP;
	XPtr<lme4Eigen::merPredD>(ptr)->setBeta0(as<Eigen::VectorXd>(beta0));
	END_RCPP;
    }
    
    SEXP merPredDsetTheta(SEXP ptr, SEXP theta) {
	BEGIN_RCPP;
	XPtr<lme4Eigen::merPredD>(ptr)->setTheta(NumericVector(theta));
	END_RCPP;
    }
    
    SEXP merPredDsetU0(SEXP ptr, SEXP u0) {
	BEGIN_RCPP;
	XPtr<lme4Eigen::merPredD>(ptr)->setU0(as<Eigen::VectorXd>(u0));
	END_RCPP;
    }
    
    SEXP merPredDLambdat(SEXP ptr) {
	BEGIN_RCPP;
	return wrap(XPtr<lme4Eigen::merPredD>(ptr)->Lambdat());
	END_RCPP;
    }
    
    SEXP merPredDLamtUt(SEXP ptr) {
	BEGIN_RCPP;
	return wrap(XPtr<lme4Eigen::merPredD>(ptr)->LamtUt());
	END_RCPP;
    }
    
    SEXP merPredDPvec(SEXP ptr) {
	BEGIN_RCPP;
	return wrap(XPtr<lme4Eigen::merPredD>(ptr)->Pvec());
	END_RCPP;
    }
    
    SEXP merPredDRX(SEXP ptr) {
	BEGIN_RCPP;
	return wrap(XPtr<lme4Eigen::merPredD>(ptr)->RX());
	END_RCPP;
    }
    
    SEXP merPredDRXdiag(SEXP ptr) {
	BEGIN_RCPP;
	return wrap(XPtr<lme4Eigen::merPredD>(ptr)->RXdiag());
	END_RCPP;
    }
    
    SEXP merPredDRZX(SEXP ptr) {
	BEGIN_RCPP;
	return wrap(XPtr<lme4Eigen::merPredD>(ptr)->RZX());
	END_RCPP;
    }
    
    SEXP merPredDZt(SEXP ptr) {
	BEGIN_RCPP;
	return wrap(XPtr<lme4Eigen::merPredD>(ptr)->Zt());
	END_RCPP;
    }
    
    SEXP merPredDVtV(SEXP ptr) {
	BEGIN_RCPP;
	return wrap(XPtr<lme4Eigen::merPredD>(ptr)->VtV());
	END_RCPP;
    }
    
    SEXP merPredDbeta0(SEXP ptr) {
	BEGIN_RCPP;
	return wrap(XPtr<lme4Eigen::merPredD>(ptr)->beta0());
	END_RCPP;
    }
    
    SEXP merPredDdelb(SEXP ptr) {
	BEGIN_RCPP;
	return wrap(XPtr<lme4Eigen::merPredD>(ptr)->delb());
	END_RCPP;
    }
    
    SEXP merPredDdelu(SEXP ptr) {
	BEGIN_RCPP;
	return wrap(XPtr<lme4Eigen::merPredD>(ptr)->delu());
	END_RCPP;
    }

    SEXP merPredDu0(SEXP ptr) {
	BEGIN_RCPP;
	return wrap(XPtr<lme4Eigen::merPredD>(ptr)->u0());
	END_RCPP;
    }
    
    SEXP merPredDunsc(SEXP ptr) {
	BEGIN_RCPP;
	return wrap(XPtr<lme4Eigen::merPredD>(ptr)->unsc());
	END_RCPP;
    }
    
    SEXP merPredDldL2(SEXP ptr) {
	BEGIN_RCPP;
	return ::Rf_ScalarReal(XPtr<lme4Eigen::merPredD>(ptr)->ldL2());
	END_RCPP;
    }
    
    SEXP merPredDldRX2(SEXP ptr) {
	BEGIN_RCPP;
	return ::Rf_ScalarReal(XPtr<lme4Eigen::merPredD>(ptr)->ldRX2());
	END_RCPP;
    }
    
    SEXP merPredDsqrL(SEXP ptr, SEXP fac) {
	BEGIN_RCPP;
	return ::Rf_ScalarReal(XPtr<lme4Eigen::merPredD>(ptr)->sqrL(::Rf_asReal(fac)));
	END_RCPP;
    }
    
    SEXP merPredDtheta(SEXP ptr) {
	BEGIN_RCPP;
	return wrap(XPtr<lme4Eigen::merPredD>(ptr)->theta());
	END_RCPP;
    }
    
    SEXP merPredDinstallPars(SEXP ptr, SEXP fac) {
	BEGIN_RCPP;
	XPtr<lme4Eigen::merPredD>(ptr)->installPars(::Rf_asReal(fac));
	END_RCPP;
    }
    
    SEXP merPredDlinPred(SEXP ptr, SEXP fac) {
	BEGIN_RCPP;
	return wrap(XPtr<lme4Eigen::merPredD>(ptr)->linPred(::Rf_asReal(fac)));
	END_RCPP;
    }
	    
}
