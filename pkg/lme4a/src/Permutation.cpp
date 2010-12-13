#if 0
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
		throw runtime_error("permutation elements must be in [0,n)");
	    if (chk[vi])
		throw runtime_error("permutation is not a permutation");
	    chk[vi] = true;
	}
    }
    
}

RCPP_FUNCTION_2(List, lme4_PermChk, IntegerVector perm, IntegerVector x) {
    IntegerVector zerob = clone(perm); // modifiable copy
    int bb = *(std::min_element(zerob.begin(), zerob.end()));
    if (bb != 0) zerob = zerob - bb;
    MatrixNs::Permutation pp(zerob);
    return List::create(_["forw"] = pp.forward(x),
			_["inv"] = pp.inverse(x));
}

#endif

