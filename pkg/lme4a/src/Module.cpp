#include "glmFamily.h"
#include "respModule.h"
#include "feModule.h"
#include "reModule.h"

using namespace Rcpp;

RCPP_MODULE(lme4) {

class_<glm::glmFamily>("glmFamily")

    .constructor<List>()

    .property("family", &glm::glmFamily::fam,
	      "family name")
    .property("link",   &glm::glmFamily::lnk,
	      "name of link function")

    .method("linkFun",  &glm::glmFamily::linkFun,
	    "apply the link function (eta from mu)")
    .method("linkInv",  &glm::glmFamily::linkInv,
	    "apply the inverse link function (mu from eta)")
    .method("muEta",    &glm::glmFamily::muEta,
	    "evaluate the gradient, d mu/d eta given eta")
    .method("devResid", &glm::glmFamily::devResid,
	    "evaluate the deviance residuals given mu, the prior weights and y")
    .method("variance", &glm::glmFamily::variance,
	    "evaluate the variance function given mu")
    ;

class_<MatrixNs::chmSp>("chmSp")

    .constructor<S4>()

    .property("nrow",  &MatrixNs::chmSp::nr,
	      "number of rows")
    .property("ncol",  &MatrixNs::chmSp::nc,
	      "number of columns")
    .property("nnz",   &MatrixNs::chmSp::nnz,
	      "number of non-zeros")
    .property("nzmax", &MatrixNs::chmSp::nzMax,
	      "maximum number of non-zeros")
    ;

using namespace mer;

class_<lmerResp>("lmerResp")

    .constructor<S4>()
    .constructor<int,NumericVector>()
    .constructor<int,NumericVector,NumericVector>()
    .constructor<int,NumericVector,NumericVector,NumericVector>()


    .property("REML",     &lmerResp::REML, &lmerResp::setReml,
	      "integer which is 0 for ML estimation and p for REML")
    .property("mu",       &lmerResp::mu, "mean vector")
    .property("offset",   &lmerResp::offset, &lmerResp::setOffset,
	      "offset vector (present even if it is all zeros")
    .property("sqrtXwt",  &lmerResp::sqrtXwt,
	      "Matrix of square roots of weights applied to X")
    .property("sqrtrwt",  &lmerResp::sqrtrwt,
	      "vector of square roots of weights applied to residuals")
    .property("weights",  &lmerResp::weights, &lmerResp::setWeights,
	      "prior weights vector (present even if it is all ones")
    .property("wrss",     &lmerResp::wrss,
	      "weighted residual sum of squares")
    .property("wtres",    &lmerResp::wtres,
	      "weighted residual vector")
    .property("y",        &lmerResp::y, "response vector")

    .method("Laplace",    &lmerResp::Laplace)
    .method("updateMu",   &lmerResp::updateMu)
    .method("updateWts",  &lmerResp::updateWts)
    ;

class_<glmerResp>("glmerResp")

    .constructor<S4>()
    .constructor<List,NumericVector>()
    .constructor<List,NumericVector,NumericVector>()
    .constructor<List,NumericVector,NumericVector,NumericVector>()
    .constructor<List,NumericVector,NumericVector,NumericVector,NumericVector>()
 
    .property("devResid",      &glmerResp::devResid,
	    "deviance residuals")
    .property("eta",           &glmerResp::eta,
	      "current value of the linear predictor")
    .property("family",        &glmerResp::family,
	      "name of the glm family")
    .property("link",          &glmerResp::link,
	      "name of the glm link")
    .property("mu",            &glmerResp::mu,
	      "current value of the mean vector")
    .property("muEta",         &glmerResp::muEta,
	      "current value of d mu/d eta")
    .property("offset",        &glmerResp::offset,
	      "offset vector - const")
    .property("pwrss",         &glmerResp::pwrss, &glmerResp::setPwrss,
	      "penalized, weighted residual sum of squares at last evaluation")
    .property("residDeviance", &glmerResp::residDeviance,
	      "sum of the deviance residuals")
    .property("sqrtXwt",       &glmerResp::sqrtXwt,
	      "square root of the weights applied to the model matrix")
    .property("sqrtrwt",       &glmerResp::sqrtrwt,
	      "square root of the weights applied to the residuals")
    .property("sqrtWrkWt",     &glmerResp::sqrtWrkWt,
	      "square root of the weights applied to the working response")
    .property("weights",       &glmerResp::weights,
	      "prior weights - const")
    .property("variance",      &glmerResp::variance,
	      "current (unscaled) variances")
    .property("wtres",         &glmerResp::wtres,
	      "current value of the weighted residuals")
    .property("wrss",          &glmerResp::wrss,
	      "weighted residual sum of squares")
    .property("wrkResids",     &glmerResp::wrkResids,
	      "working residuals - on the eta scale")
    .property("wrkResp",       &glmerResp::wrkResp,
	      "working response - on the eta scale")
    .property("y",             &glmerResp::y,
	      "numeric response vector - const")

    .method("Laplace",         &glmerResp::Laplace,
	    "Evaluate the Laplace approximation to the deviance, given ldL2, ldRX2 (not used) and sqrLenU")
    .method("updateMu",        &glmerResp::updateMu,
	    "update mu and derived quantities from a new value of eta")
    .method("updateWts",       &glmerResp::updateWts,
	    "update the residual and X weights.")
    ;

class_<deFeMod>("deFeMod")

    .constructor<S4,int>()
    .constructor<S4,int,int,int>()

    .property("CcNumer",       &deFeMod::CcNumer,
	      "contribution to numerator of the convergence criterion")
    .property("coef",          &deFeMod::coef,
	      "coefficient vector")
    .property("coef0",         &deFeMod::coef0, &deFeMod::setCoef0,
	      "base coefficient vector")
    .property("incr",          &deFeMod::incr, &deFeMod::setIncr,
	      "full increment of coefficient vector")
    .property("ldRX2",         &deFeMod::ldRX2,
	      "log of the square of the determinant of RX")
    .property("RX",            &deFeMod::RX,
	      "Cholesky factor")
    .property("RZX",           &deFeMod::RZX,
	      "Off-diagonal block in large Cholesky factor")
    .property("UtV",           &deFeMod::UtV,
	      "Weighted crossproduct")
    .property("V",             &deFeMod::V,
	      "scaled model matrix")
    .property("Vtr",           &deFeMod::Vtr,
	      "weighted cross product of model matrix and residual")
    .property("VtV",           &deFeMod::VtV,
	      "Weighted crossproduct")
    .property("X",             &deFeMod::X,
	      "original model matrix")

    .method("installCoef0",    &deFeMod::installCoef0,
	    "install coef as coef0")
    .method("linPred1",        &deFeMod::linPred1,
	    "update coef = coef0 + fac * incr and return linear predictor contribution")
    .method("reweight",        &deFeMod::reweight,
	    "update V and Vtr for new weights")
    .method("solveIncr",       &deFeMod::solveIncr,
	    "update fac from VtV and solve for increment")
    .method("updateIncr",      &deFeMod::updateIncr,
	    "update the increment vector given cu")
    .method("updateRzxpRxpp",  &deFeMod::updateRzxpRxpp,
	    "update the triangular factor sections given external pointers to Lambda and L")
    .method("updateUtV",      &deFeMod::updateUtVp,
	    "update UtV given a pointer to Ut")
    ;

class_<reModule>("reModule")

    .constructor<S4>()
    .constructor<S4,S4,S4,IntegerVector,NumericVector>()

    .property("b",         &reModule::b,
	      "random effects on original scale")
    .property("CcNumer",   &reModule::CcNumer,
	      "contribution to numerator of the convergence criterion")
    .property("cu",        &reModule::cu,
	      "intermediate solution for u")
    .property("incr",      &reModule::incr, &reModule::setIncr,
	      "full increment for u")
    .property("L",         &reModule::L,
	      "sparse Cholesky factor")
    .property("Lambda",    &reModule::Lambda,
	      "relative covariance factor")
    .property("Lambdap",   &reModule::Lambdap,
	      "external pointer to the relative covariance factor")
    .property("ldL2",      &reModule::ldL2,
	      "logarithm of the square of the determinant of L")
    .property("Lind",      &reModule::Lind,
	      "1-based index vector into theta for Lambda@x")
    .property("linPred",   &reModule::linPred,
	      "linear predictor contribution")
    .property("lower",     &reModule::lower,
	      "lower bounds on the theta parameters")
    .property("Lp",        &reModule::Lp,
	      "external pointer to the sparseCholesky factor")
    .property("sqrLenU",   &reModule::sqrLenU,
	      "squared length of u, the orthogonal random effects")
    .property("theta",     &reModule::theta,   &reModule::setTheta,
	      "current value of variance component parameters")
    .property("u",         &reModule::u,
	      "orthonormal random effects vector")
    .property("u0",        &reModule::u0, &reModule::setU0,
	      "base orthonormal random effects vector")
    .property("Ut",        &reModule::Utp,
	      "pointer to the weighted, orthogonal design matrix")
    .property("Zt",        &reModule::Zt,
	      "transpose of the model matrix for the random effects")

    .method("installU0",   &reModule::installU0,
	    "install u as u0")
    .method("reweight",    &reModule::reweight,
	    "update L, Ut and cu for new weights")
    .method("solveU",      &reModule::solveU,
	    "solve for u (or the increment for u) only.  Returns squared length of c1")
    .method("solveIncr",   &reModule::solveIncr,
	    "solve for the increment to u only.  Returns squared length of c1")
    .method("updateIncr",  &reModule::updateIncr,
	    "update the increment given the updated cu from the feModule's updateIncr method")
    .method("linPred1",    &reModule::linPred1,
	    "update u = u0 + fac * incr and return linear predictor contribution")
    ;

class_<reTrms>("reTrms")
    .constructor<S4>()

    .property("assign",        &reTrms::assign,
	      "assignment of terms to grouping factors")
    .property("cnms",          &reTrms::cnms,
	      "list of column names per term")
    .property("flist",         &reTrms::flist,
	      "list of grouping factors")
    .property("ncols",         &reTrms::ncols,
	      "number of columns per term")
    .property("nctot",         &reTrms::nctot,
	      "total number of columns per grouping factor")
    .property("nlevs",         &reTrms::nlevs,
	      "number of levels per grouping factor (after removing unused levels)")
    .property("offsets",       &reTrms::offsets,
	      "offsets per term into random effects")

    .method("terms",         &reTrms::terms,
	    "returns the terms associated with a grouping factor (argument is a 0-based index)")
    .method("condVar",       &reTrms::condVar,
	    "returns a list of 3D arrays, argument is scalar scale factor")
    ;
}
