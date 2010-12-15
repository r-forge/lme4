#include "glmFamily.h"
#include "respModule.h"
#include "feModule.h"
#include "reModule.h"

using namespace Rcpp;

RCPP_MODULE(lme4) {

class_<glm::glmFamily>("glmFamily")

    .constructor<List>()

    .property("family", &glm::glmFamily::fam)
    .property("link",   &glm::glmFamily::lnk)

    .method("linkFun",  &glm::glmFamily::linkFun)
    .method("linkInv",  &glm::glmFamily::linkInv)
    .method("muEta",    &glm::glmFamily::muEta)
    .method("devResid", &glm::glmFamily::devResid)
    .method("variance", &glm::glmFamily::variance)
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

//    .property("devResid", &lmerResp::devResid)
    .property("mu",       &lmerResp::mu)
    .property("offset",   &lmerResp::offset)
    .property("sqrtXwt",  &lmerResp::sqrtXwt)
    .property("sqrtrwt",  &lmerResp::sqrtrwt)
    .property("weights",  &lmerResp::weights)
    .property("wrss",     &lmerResp::wrss)
    .property("wtres",    &lmerResp::wtres)
    .property("y",        &lmerResp::y)

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
    .property("residDeviance", &glmerResp::residDeviance,
	      "sum of the deviance residuals")
    .property("sqrtXwt",       &glmerResp::sqrtXwt,
	      "square root of the weights applied to the model matrix")
    .property("sqrtrwt",       &glmerResp::sqrtrwt,
	      "square root of the weights applied to the residuals")
    .property("sqrtWrkWt",     &glmerResp::sqrtrwt,
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

    .method("updateMu",        &glmerResp::updateMu,
	    "update mu and derived quantities from a new value of eta")
    .method("updateWts",       &glmerResp::updateWts,
	    "update the residual and X weights.")
    ;

class_<deFeMod>("deFeMod")

    .constructor<S4,int>()

    .property("coef",          &deFeMod::coef,
	      "coefficient vector")
    .property("Vtr",           &deFeMod::Vtr,
	      "weighted cross product of model matrix and residual")
    .property("ldRX2",         &deFeMod::ldRX2,
	      "log of the square of the determinant of RX")
    .property("V",             &deFeMod::V,
	      "scaled model matrix")
    .property("RX",            &deFeMod::RX,
	      "Cholesky factor")
    .property("RZX",           &deFeMod::RZX,
	      "Off-diagonal block in large Cholesky factor")
    .property("UtV",           &deFeMod::UtV,
	      "Weighted crossproduct")
    .property("VtV",           &deFeMod::VtV,
	      "Weighted crossproduct")
    .method("updateBeta",      &deFeMod::updateBeta,
	    "update the coefficient vector given cu")
    ;

class_<reModule>("reModule")

    .constructor<S4>()
    .constructor<S4,S4,S4,IntegerVector,NumericVector>()

    .property("sqrLenU",       &reModule::sqrLenU,
	      "squared length of u, the orthogonal random effects")
    .property("linPred",       &reModule::linPred,
	      "linear predictor contribution")
    .property("ldL2",          &reModule::ldL2,
	      "logarithm of the square of the determinant of L")
    .property("cu",            &reModule::cu,
	      "intermediate solution for u")
    .property("Lind",          &reModule::Lind,
	      "1-based index vector into theta for Lambda@x")
    .property("lower",         &reModule::lower,
	      "lower bounds on the theta parameters")
    .property("theta",         &reModule::theta,   &reModule::setTheta,
	      "current value of variance component parameters")
    .property("u",             &reModule::u,
	      "orthonormal random effects vector")
    .property("Zt",            &reModule::Zt,
	      "transpose of the model matrix for the random effects")
    .property("Lambda",        &reModule::Lambda,
	      "relative covariance factor")
    .property("L",             &reModule::L,
	      "sparse Cholesky factor")

    .method("reweight",        &reModule::reweight,
	    "update L, Ut and cu for new weights")
    .method("setU",            &reModule::setU,
	    "set a new value of u, possibly with an increment and step")
    .method("solveU",          &reModule::solveU,
	    "solve for u (or the increment for u) only.  Returns squared length of c1")
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
