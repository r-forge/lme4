#include "mer.h"
#include <R_ext/BLAS.h>

using namespace Rcpp;
using namespace MatrixNs;
using namespace std;

namespace mer{ // utilities defined here, class constructors and
	       // member methods in class-specific source files

    void showdbl(const double* x, const char* nm, int n) {
	if (n < 1) {
	    Rprintf("%s[%d]:\n", nm, n);
	    return;
	}
	int n5 = (n < 5) ? n : 5;
	Rprintf("%s[1:%d]: %g", nm, n, x[0]);
	for (int i = 1; i < n5; i++) Rprintf(", %g", x[i]);
	if (n > 5) Rprintf(", ...");
	Rprintf("\n");
    }

    void showdbl(const Rcpp::NumericVector& vv, const char* nm) {
	showdbl(vv.begin(), nm, vv.size());
    }

    void showincr(double step, double c0, double c1,
		  NumericVector const& incr, const char* nm) {
	Rprintf("step = %8.5f, pwrss0 = %12g, pwrss1 = %12g\n",
		step, c0, c1);
	showdbl(incr.begin(), nm, incr.size());
    }

    void showint(const int* x, const char* nm, int n) {
	if (n < 1) {
	    Rprintf("%s[%d]:\n", nm, n);
	    return;
	}
	int n20 = (n < 20) ? n : 20;
	Rprintf("%s[1:%d]: %d", nm, n, x[0]);
	for (int i = 1; i < n20; i++) Rprintf(", %d", x[i]);
	if (n > 20) Rprintf(", ...");
	Rprintf("\n");
    }

    void showCHM_DN(const cholmod_dense* x, const string &nm) {
	Rprintf("%s: nrow = %d, ncol = %d, nzmax = %d, d = %d, xtype = %d, dtype = %d\n",
		nm.c_str(), x->nrow, x->ncol, x->nzmax, x->d, x->xtype, x->dtype);
	showdbl((double*)x->x, "x", x->nzmax);
    }

    void showCHM_FR(const cholmod_factor* x, const string &nm) {
	Rprintf("%s: n = %d, minor = %d, xtype = %d, itype = %d, dtype = %d\n",
		nm.c_str(), x->n, x->minor, x->xtype, x->itype, x->dtype);
	Rprintf("ordering = %d, is_ll = %d, is_super = %d, is_monotonic = %d, nzmax = %d\n",
		x->ordering, x->is_ll, x->is_super, x->is_monotonic, x->nzmax);
	int *pp = (int*)x->p;
	showint(pp, "p", (x->n) + 1);
	showint((int*)x->Perm, "Perm", x->n);
	if (pp[x->n] > 0) {
	    showint((int*)x->i, "i", x->n);
	    showdbl((double*)x->x, "x", x->n);
	} else Rprintf("nnz = 0\n");
    }

    void showCHM_SP(const cholmod_sparse* x, const string &nm) {
	Rprintf("%s: nrow = %d, ncol = %d, xtype = %d, stype = %d, itype = %d, dtype = %d\n",
		nm.c_str(), x->nrow, x->ncol, x->xtype, x->stype,
		x->itype, x->dtype);
	int nc = x->ncol, nnz = M_cholmod_nnz(x, &c);
	showint((int*)x->p, "p", nc + 1);
	if (nnz > 0) {
	    showint((int*)x->i, "i", nnz);
	    showdbl((double*)x->x, "x", nnz);
	} else Rprintf("nnz = 0\n");
    }

    void showCHM_SP(const MatrixNs::chmSp& x, const std::string& nm) {
	showCHM_SP(&x, nm);
    }

    void showdMat(MatrixNs::ddenseMatrix const& x, string const& nm) {
	int m = x.ncol(), n = x.nrow();
	int mn = m * n;
	double *xx = x.x().begin();
	if (mn < 1) {
	    Rprintf("%s[%d,%d]:\n", nm.c_str(), m, n);
	    return;
	}
	int mn7 = (mn < 7) ? mn : 5;
	Rprintf("%s[1:%d,1:%d]: %g", nm.c_str(), m, n, xx[0]);
	for (int i = 1; i < mn7; i++) Rprintf(", %g", xx[i]);
	if (mn > 7) Rprintf(", ...");
	Rprintf("\n");
    }

    /** 
     * Determine the weighted Euclidean distance between two vectors,
     * relative to the square root of the product of their lengths
     * 
     * @param v1 First vector to compare
     * @param v2 Second vector to compare
     * @param wt square root of the weights
     * 
     * @return relative difference between the matrices
     */
    double compareVecWt(NumericVector const& v1,
			NumericVector const& v2,
			NumericVector const& wt) {
	int n = v1.size();
	double num, d1, d2;
	if (v2.size() != n || wt.size() != n)
	    Rf_error("%s: size mismatch, %d != %d or != %d\n",
		     "compareVecWt", n, v2.size(), wt.size());
	vector<double> a(n);

	transform(v1.begin(), v1.end(), v2.begin(),
		       a.begin(), minus<double>());
	transform(a.begin(), a.end(), wt.begin(), a.begin(),
		       multiplies<double>());
	num = inner_product(a.begin(), a.end(), a.begin(), double());

	transform(v1.begin(), v1.end(), wt.begin(), a.begin(),
		       multiplies<double>());
	d1 = inner_product(a.begin(), a.end(), a.begin(), double());
	transform(v2.begin(), v2.end(), wt.begin(), a.begin(),
		       multiplies<double>());
	d2 = inner_product(a.begin(), a.end(), a.begin(), double());

	if (d1 == 0) {
	    return (d2 == 0) ? 0. : sqrt(num/d2);
	}
	return (d2 == 0) ? sqrt(num/d1) : sqrt(num/sqrt(d1 * d2));
    }

    NumericVector mkans(double x, const NumericVector& beta,
			const NumericVector& u) {
	NumericVector ans(1);
	*ans.begin() = x;
	ans.attr("beta") = beta;
	ans.attr("u") = u;
	return ans;
    }
} // namespace mer

/** 
 * Evaluate the profiled deviance or the profiled REML criterion from
 * a theta parameter value for an lmer model.
 * 
 * @param xp merMod object
 * @param nt new value of parameters - either theta or c(theta,beta)
 * 
 * @return a deviance evaluation
 */

RCPP_FUNCTION_5(NumericVector, merDeviance, S4 xp, NumericVector pars, NumericVector u0, int verb, int alg) {
    S4 fe(xp.slot("fe")), resp(xp.slot("resp"));
    bool de(fe.is("deFeMod"));
    if (!de && !fe.is("spFeMod"))
	throw runtime_error("fe slot is neither deFeMod nor spFeMod");

    if (resp.is("lmerResp")) {
	if (de) {
	    mer::mer<mer::deFeMod,mer::lmerResp> lm(xp);
	    return lm.LMMdeviance(pars, u0);
	} else {
	    mer::mer<mer::spFeMod,mer::lmerResp> lm(xp);
	    return lm.LMMdeviance(pars, u0);
	}
    }

    if (alg < 1 || alg > 3) throw range_error("alg must be 1, 2 or 3");
    mer::Alg aa = (alg == 1) ? mer::Beta : ((alg == 2) ? mer::U : mer::BetaU);
// FIXME: Need to set beta for alg == 3
    if (resp.is("glmerResp")) {
	if (de) {
	    mer::mer<mer::deFeMod,mer::glmerResp> glmr(xp);
	    return glmr.PIRLS(u0, verb, aa);
	} else {
	    mer::mer<mer::spFeMod,mer::glmerResp> glmr(xp);
	    return glmr.PIRLS(u0, verb, aa);
	}
    } else if (resp.is("nlmerResp")) {
	if (de) {
	    mer::mer<mer::deFeMod,mer::nlmerResp> nlmr(xp);
	    return nlmr.PIRLS(u0, verb, aa);
	} else {
	    mer::mer<mer::spFeMod,mer::nlmerResp> nlmr(xp);
	    return nlmr.PIRLS(u0, verb, aa);
	}
    }

    throw runtime_error("resp slot is not lmerResp or glmerResp or nlmerResp");
}
#if 0
RCPP_FUNCTION_VOID_2(feSetBeta, S4 xp, NumericVector nbeta) {
    mer::feModule fe(xp);
    fe.setBeta(nbeta);
}

RCPP_FUNCTION_VOID_2(reUpdateLambda, S4 xp, NumericVector nth) {
    mer::reModule re(xp);
    re.updateLambda(nth);
}
#endif
RCPP_FUNCTION_VOID_1(updateRzxRx, S4 xp) {
    S4 fe(xp.slot("fe")), resp(xp.slot("resp"));
    bool de = fe.is("deFeMod");
    if (!de && !fe.is("spFeMod"))
	throw runtime_error("fe slot is neither deFeMod nor spFeMod");
    if (resp.is("glmerResp")) {
	if (de) {
	    mer::mer<mer::deFeMod,mer::glmerResp> glmr(xp);
	    glmr.updateRzxRx();
	} else {
	    mer::mer<mer::spFeMod,mer::glmerResp> glmr(xp);
	    glmr.updateRzxRx();
	}
    } else if (resp.is("nlmerResp")) {
	if (de) {
	    mer::mer<mer::deFeMod,mer::nlmerResp> nlmr(xp);
	    nlmr.updateRzxRx();
	} else {
	    mer::mer<mer::spFeMod,mer::nlmerResp> nlmr(xp);
	    nlmr.updateRzxRx();
	}
    } else 
	throw runtime_error("resp slot is not glmerResp or nlmerResp in updateRzxRx");
}

RCPP_FUNCTION_2(List, updateDc, S4 xp, NumericVector pars) {
    List ans = List::create(_["devcomp"] = clone(SEXP(xp.slot("devcomp"))));

    S4 fe(xp.slot("fe")), resp(xp.slot("resp"));
    bool de(fe.is("deFeMod"));
    if (!de && !fe.is("spFeMod"))
	throw runtime_error("fe slot is neither deFeMod nor spFeMod");

    if (resp.is("lmerResp")) {
	if (de) {
	    mer::mer<mer::deFeMod,mer::lmerResp> lm(xp);
	    lm.updateDcmp(ans, pars);
	    return ans;
	} else {
	    mer::mer<mer::spFeMod,mer::lmerResp> lm(xp);
	    lm.updateDcmp(ans, pars);
	    return ans;
	}
    }
    if (resp.is("glmerResp")) {
	if (de) {
	    mer::mer<mer::deFeMod,mer::glmerResp> glmr(xp);
	    glmr.updateDcmp(ans, pars);
	    return ans;
	} else {
	    mer::mer<mer::spFeMod,mer::glmerResp> glmr(xp);
	    glmr.updateDcmp(ans, pars);
	    return ans;
	}
    } else if (resp.is("nlmerResp")) {
	if (de) {
	    mer::mer<mer::deFeMod,mer::nlmerResp> nlmr(xp);
	    nlmr.updateDcmp(ans, pars);
	    return ans;
	} else {
	    mer::mer<mer::spFeMod,mer::nlmerResp> nlmr(xp);
	    nlmr.updateDcmp(ans, pars);
	    return ans;
	}
    }

    throw runtime_error("resp slot is not lmerResp or glmerResp or nlmerResp");
}    
RCPP_FUNCTION_2(List, reTrmsCondVar, S4 xp, double scale) {
    mer::reTrms trms(xp);
    return trms.condVar(scale);
}
