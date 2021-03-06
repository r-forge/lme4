#include "reModule.h"

using namespace Rcpp;
using namespace MatrixNs;
using namespace std;

namespace mer{

    reTrms::reTrms(Rcpp::S4 xp)
	: reModule(xp),
	  d_flist(xp.slot("flist")),
	  d_cnms(xp.slot("cnms")),
	  d_assign(d_flist.attr("assign")) {
    }

    Rcpp::IntegerVector reTrms::ncols() const {
	int nt = d_cnms.size();
	IntegerVector ans(nt);
	int *ap = ans.begin();
	for (int i = 0; i < nt; i++) {
	    CharacterVector nmi = d_cnms(i);
	    ap[i] = nmi.size();
	}
	return ans;
    }

    Rcpp::IntegerVector reTrms::nctot() const {
	IntegerVector ans(d_flist.size()), nc = ncols();
	int *asg = d_assign.begin(), *ansp = ans.begin(), *ncp = nc.begin(), nnc = nc.size();
	fill(ans.begin(), ans.end(), int());
	for (int i = 0; i < nnc; i++) ansp[asg[i] - 1] += ncp[i];
	return ans;
    }

    Rcpp::IntegerVector reTrms::nlevs() const {
	int nfac = d_flist.size();
	IntegerVector ans(nfac);
	int *ap = ans.begin();
	for (int i = 0; i < nfac; i++) {
	    IntegerVector fi = d_flist(i);
	    CharacterVector ll = fi.attr("levels");
	    ap[i] = ll.size();
	}
	return ans;
    }

    Rcpp::IntegerVector reTrms::offsets() const {
	int *asgn = d_assign.begin(), ntrm = d_cnms.size();
	IntegerVector ans(ntrm), nc = ncols(), nl = nlevs();
	int *ansp = ans.begin(), *ncp = nc.begin(),  *nlp = nl.begin(), offset = 0;
	for (R_len_t i = 0; i < ntrm; ++i) {
	    ansp[i] = offset;
	    offset += ncp[i] * nlp[asgn[i] - 1];
	}
	return ans;
    }

    Rcpp::IntegerVector reTrms::terms(R_len_t ii) const {
	R_len_t ip1 = ii + 1;   // for matching 1-based indices in d_assign
	IntegerVector ans(count(d_assign.begin(), d_assign.end(), ip1));
	int *asgn = d_assign.begin(), ai = 0;
	for (R_len_t i = 0; i < d_assign.size(); ++i)
	    if (asgn[i] == ip1) ans[ai++] = i;
	return ans;
    }

    Rcpp::List reTrms::condVar(double scale) {
	if (scale < 0 || !R_finite(scale))
	    throw runtime_error("scale must be non-negative and finite");
	int nf = d_flist.size();
	IntegerVector nc = ncols(), nl = nlevs(), nct = nctot(), off = offsets();
	
    	List ans(nf);
	CharacterVector nms = d_flist.names();
	ans.names() = clone(nms);
	
	for (int i = 0; i < nf; i++) {
	    int ncti = nct[i], nli = nl[i];
	    IntegerVector trms = terms(i);
	    int *cset = new int[ncti], nct2 = ncti * ncti;
	    
	    NumericVector ansi(Dimension(ncti, ncti, nli));
	    ans[i] = ansi;
	    double *ai = ansi.begin();
	    
	    for (int j = 0; j < nli; j++) {
		int kk = 0;
		for (int jj = 0; jj < trms.size(); jj++) {
		    int tjj = trms[jj];
		    for (int k = 0; k < nc[tjj]; k++)
			cset[kk++] = off[tjj] + j * nc[tjj] + k;
		}
		CHM_SP cols =
		    M_cholmod_submatrix(&d_Lambda, (int*)NULL, -1, cset, ncti,
					1/*values*/, 1/*sorted*/, &c);
		CHM_SP sol = d_L.spsolve(CHOLMOD_A, cols);
		CHM_SP tcols = M_cholmod_transpose(cols, 1/*values*/, &c);
		M_cholmod_free_sparse(&cols, &c);
		CHM_SP var = M_cholmod_ssmult(tcols, sol, 0/*stype*/,
					      1/*values*/, 1/*sorted*/, &c);
		M_cholmod_free_sparse(&sol, &c);
		M_cholmod_free_sparse(&tcols, &c);
		CHM_DN dvar = M_cholmod_sparse_to_dense(var, &c);
		M_cholmod_free_sparse(&var, &c);
		Memcpy(ai + j * nct2, (double*)dvar->x, nct2);
		M_cholmod_free_dense(&dvar, &c);
	    }
	    delete[] cset;
	    transform(ansi.begin(), ansi.end(), ansi.begin(),
		      bind2nd(multiplies<double>(), scale * scale));
	}
	return ans;
    }
}

