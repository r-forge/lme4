#include "utilities.h"

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

    NumericVector mkans(double x,
			const NumericVector& beta,
			const NumericVector&    u,
			double               ldL2,
			double               wrss,
			double               ussq) {
	NumericVector ans(1, x);
	ans.attr("beta") = beta;
	ans.attr("u")    = u;
	ans.attr("ldL2") = ldL2;
	ans.attr("wrss") = wrss;
	ans.attr("ussq") = ussq; // is this needed in addition to u?
	return ans;
    }
} // namespace mer
