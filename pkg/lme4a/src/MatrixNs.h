// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-

// Can't name the file Matrix.h because of the Matrix.h in Matrix/include

#ifndef LME4_MATRIX_H
#define LME4_MATRIX_H

//#include <RcppCommon.h>

// namespace Rcpp {
//     template <int RTYPE> class Vector;
//     template <int RTYPE> class Matrix;
//     class S4;
//     typedef Vector<REALSXP> NumericVector;
//     typedef Vector<VECSXP>  List;
//     typedef Matrix<REALSXP> NumericMatrix;
// }

#include <Rcpp.h>
#include <Matrix.h>

extern cholmod_common c;

namespace MatrixNs {
				// utilities
    char UpLo(char);		// checks for and returns 'U' or 'L' 
    char UpLo(std::string const&);
    char UpLo(const SEXPREC*);

    char Diag(char);		// checks for and returns 'U' or 'N' 
    char Diag(std::string const&);
    char Diag(const SEXPREC*);

    char Trans(char);		// checks for and returns 'T' or 'N' 
    char Trans(std::string const&);
    char Trans(const SEXPREC*);

    class Matrix {
    protected:
	SEXP d_sexp;
	Rcpp::List d_dimnames;
	int d_nrow, d_ncol;
    public:
	Matrix(Rcpp::S4&);
	Matrix(int,int);
	int nrow() const;
	int ncol() const;
	const SEXPREC* sexp() const {return d_sexp ? d_sexp : R_NilValue;}
//	operator SEXP() const {return d_sexp ? d_sexp : R_NilValue;}
    };

    class dMatrix : public Matrix {
    protected:
	Rcpp::NumericVector d_x;
    public:
	dMatrix(Rcpp::S4&);
	dMatrix(int,int,int=0);

	const Rcpp::NumericVector& x() const {return d_x;}
	Rcpp::NumericVector& X() {return d_x;}
	void setX(const Rcpp::NumericVector&);
	void setX(const Rcpp::NumericMatrix&);
    };

    class ddenseMatrix : public dMatrix {
    public:
	ddenseMatrix(Rcpp::S4&);
	ddenseMatrix(int,int);
    };

// C++ classes mirroring virtual S4 structure classes do not inherit
// from Matrix, so as to avoid multiple definitions of Dim and
// Dimnames slots.  You can get around this for concrete classes but
// these are pure virtual classes.

    class compMatrix {		// composite (factorizable) Matrix
    protected:
	Rcpp::List factors;
    public:
	compMatrix(Rcpp::S4&);
	compMatrix(){};		// factor list is empty
    };


    class generalMatrix : public compMatrix { //< general structure
    public:
	generalMatrix(Rcpp::S4&);
	generalMatrix(){};
    };

    class triangularMatrix {
    protected:
	char d_ul, d_di;
    public:
	triangularMatrix(Rcpp::S4&);
	triangularMatrix(char='U',char='N');

	char diag() const {return d_di;} 
	char uplo() const {return d_ul;} 

    };
    
    class symmetricMatrix : public compMatrix {
    protected:
	char d_ul;
    public:
	symmetricMatrix(Rcpp::S4&);
	symmetricMatrix(char='U');

	char uplo() const {return d_ul;} 
    };

// Concrete classes are initialized from Rcpp::S4 not Rcpp::S4& so
// they can be constructed from slots extracted during the
// construction of composite objects.

    class chmDn; 		// forward declaration

    class dgeMatrix : public ddenseMatrix, public generalMatrix {
    public:
	dgeMatrix(Rcpp::S4);
	dgeMatrix(int,int);

	int  dmult(char,double,double, const chmDn&,chmDn&) const;
	void dgemv(char,double,const Rcpp::NumericVector&,
		   double,Rcpp::NumericVector&) const;
	void dgemv(char,double,const Rcpp::NumericVector&,
		   double,double*) const;
	void dgemm(char,char,double,const dgeMatrix&,
		   double,dgeMatrix&) const;
    };

    class dtrMatrix : public ddenseMatrix, public triangularMatrix {
    public:
	dtrMatrix(Rcpp::S4&);
	dtrMatrix(int,char='U',char='N');

//	void dtrtrs(char,Rcpp::NumericVector&) const;
//	void dtrtrs(char,std::vector<double>&) const;
	void dtrtrs(char,double*,int = 1) const;
    };

    class dsyMatrix : public ddenseMatrix, public symmetricMatrix {
    public:
	dsyMatrix(Rcpp::S4&);
	dsyMatrix(int,char='U');

	void dsyrk(dgeMatrix const&,double,double);
    };

    class dpoMatrix : public dsyMatrix {
    public:
	dpoMatrix(Rcpp::S4&);
	dpoMatrix(int,char='U');
    };

    class Cholesky : public dtrMatrix {
    public:
	Cholesky(Rcpp::S4);
	Cholesky(dgeMatrix, char = 'U');

	Rcpp::NumericMatrix solve(int,const_CHM_DN) const;
	Rcpp::NumericMatrix solve(int,const Rcpp::NumericMatrix&) const;
	Rcpp::NumericMatrix solve(int,const Rcpp::NumericVector&) const;

	double logDet2() const;

	void dpotrs(Rcpp::NumericVector&) const;
	void dpotrs(std::vector<double>&) const;
	void dpotrs(double*, int = 1) const;

	void inPlaceSolve(int, Rcpp::NumericMatrix&) const;

	void update(const dpoMatrix&); // chol(A)
	void update(const dgeMatrix&); // chol(crossprod(X))
	void update(char,double,const dgeMatrix&,double,const dsyMatrix&);
    };

    class chmDn : public cholmod_dense {
	void init(const double*, int, int);
    public:
	chmDn(double*, int, int);
	chmDn(const double*, int, int);
	chmDn(std::vector<double> &);
	chmDn(std::vector<double> const&);
	chmDn(Rcpp::NumericVector&);
	chmDn(Rcpp::NumericVector const&);
	chmDn(Rcpp::NumericMatrix&);
	chmDn(Rcpp::NumericMatrix const&);
    	chmDn(ddenseMatrix&);
    	chmDn(ddenseMatrix const&);
//	chmDn(CHM_DN);
//	~chmDn() {if (pp) ::M_cholmod_free_dense(&pp, &c);}

	int nr() const { return nrow; }
	int nc() const { return ncol; }
	double* begin() {return (double*)x;} // template this
//	const double* begin() {return const (double*)x;}
	double* end() {return begin() + nrow * ncol;}
//	const double* end() {
//	    const double *ee = begin() + nrow * ncol;
//	    return ee;
//	}
//    protected:
//	CHM_DN pp;
    };

    class chmSp : public cholmod_sparse { 
	Rcpp::S4        d_xp;
    public:
	chmSp(Rcpp::S4);

	const Rcpp::S4&  S4() const {return d_xp;}
	const SEXPREC* sexp() const {return SEXP(d_xp);}
	CHM_SP crossprod() const;
	CHM_SP crossprod(const_CHM_SP, int = 1) const;
	CHM_SP crossprod(const chmSp&, int = 1) const;

	CHM_SP tcrossprod() const;
	CHM_SP tcrossprod(const_CHM_SP, int = 1) const;
	CHM_SP tcrossprod(const chmSp&, int = 1) const;

	CHM_SP transpose(int values = 1) const;

	CHM_SP smult(const chmSp&, int, int, int) const;
	int dmult(char,double,double,const chmDn&,chmDn&) const;
	
	void scale(int,const chmDn&);
	void update(const cholmod_sparse&);
    };

    class chmFr : public cholmod_factor {
	Rcpp::S4        d_xp;
    public:
	chmFr(Rcpp::S4);

	const Rcpp::S4&  S4() const {return d_xp;}
	const SEXPREC* sexp() const {return SEXP(d_xp);}
//      double logDet2() const {   // Need Matrix_0.999375-42 or later
	double logDet2();

	void update(const cholmod_sparse&, double Imult = 0.); 

	Rcpp::NumericMatrix solve(int, const_CHM_DN) const;
	Rcpp::NumericMatrix solve(int, const Rcpp::NumericMatrix&) const;
	Rcpp::NumericMatrix solve(int, const Rcpp::NumericVector&) const;

	CHM_SP spsolve(int sys, const_CHM_SP b) const;
	CHM_SP spsolve(int sys, const chmSp&b) const;
    };

    class Permutation {
	Rcpp::IntegerVector d_perm;
	int n;
    public:
	Permutation(Rcpp::IntegerVector&);

	template<int Rtype>
	Rcpp::Vector<Rtype> forward(const Rcpp::Vector<Rtype>&) const;
	template<int Rtype>
	Rcpp::Vector<Rtype> inverse(const Rcpp::Vector<Rtype>&) const;
	template<int Rtype>
	void inPlaceFwd(Rcpp::Vector<Rtype>&) const;
	template<int Rtype>
	void inPlaceInv(Rcpp::Vector<Rtype>&) const;
    };

    template<int Rtype> inline
    Rcpp::Vector<Rtype> Permutation::forward(const Rcpp::Vector<Rtype>& vv) const {
    	if (vv.size() != n)
	    throw std::runtime_error("size mismatch in permutation");
	Rcpp::Vector<Rtype> ans(n);
	int *ppt = d_perm.begin();
	for (R_len_t i = 0; i < n; ++i) ans[i] = vv[ppt[i]];
	return ans;
    }

    template<int Rtype> inline
    Rcpp::Vector<Rtype> Permutation::inverse(const Rcpp::Vector<Rtype>& vv) const {
    	if (vv.size() != n)
	    throw std::runtime_error("size mismatch in permutation");
	Rcpp::Vector<Rtype> ans(n);
	int *ppt = d_perm.begin();
	for (R_len_t i = 0; i < n; ++i) ans[ppt[i]] = vv[i];
	return ans;
    }
    
    template<int Rtype> inline
    void Permutation::inPlaceFwd(Rcpp::Vector<Rtype>& vv) const {
	Rcpp::Vector<Rtype> ans = forward(vv);
	std::copy(ans.begin(), ans.end(), vv.begin());
    }

    template<int Rtype> inline
    void Permutation::inPlaceInv(Rcpp::Vector<Rtype>& vv) const {
	Rcpp::Vector<Rtype> ans = inverse(vv);
	std::copy(ans.begin(), ans.end(), vv.begin());
    }
}

#endif
