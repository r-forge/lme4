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
    
}

					       
