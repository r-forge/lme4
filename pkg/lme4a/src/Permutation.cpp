#include "MatrixNs.h"

using namespace Rcpp;
using namespace std;

namespace MatrixNs {

    Permutation::Permutation(Rcpp::IntegerVector &vv)
	: d_perm(vv),
	  n(vv.size()) {
	int *vpt = vv.begin();
	std::vector<bool> chk(n);
	std::fill(chk.begin(), chk.end(), false);
	for (int i = 0; i < n; i++) {
	    int vi = vpt[i];
	    if (vi < 0 || n <= vi)
		throw std::runtime_error("permutation elements must be in [0,n)");
	    if (chk[vi])
		throw std::runtime_error("permutation is not a permutation");
	    chk[vi] = true;
	}
    }
    
    Rcpp::NumericVector
    Permutation::forward(Rcpp::NumericVector const& vv) const {
	if (vv.size() != n)
	    throw runtime_error("size mismatch in permutation");
	NumericVector ans(n);
	double *vpt = vv.begin(), *apt = ans.begin();
	int *ppt = d_perm.begin();
	for (int i = 0; i < n; ++i) apt[i] = vpt[ppt[i]];
	return ans;
    }
}

					       
