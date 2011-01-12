// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-

// Can't name the file Matrix.h because of the Matrix.h in Matrix/include

#ifndef LME4_MATRIX_H
#define LME4_MATRIX_H

#include <RcppCommon.h>
// forward declarations
namespace MatrixNs {
    class Cholesky;
    class chmFr;
    class chmSp;
    class ddenseModelMatrix;
    class dgeMatrix;
    class dpoMatrix;
    class dsyMatrix;
    class dtrMatrix;
}

namespace Rcpp {
    namespace traits {
	template <> class
	is_convertible<SEXP,MatrixNs::Cholesky> : public true_type{};
	template <> class
	is_convertible<SEXP,MatrixNs::chmFr> : public true_type{};
	template <> class
	is_convertible<SEXP,MatrixNs::chmSp> : public true_type{};
	template <> class
	is_convertible<SEXP,MatrixNs::dgeMatrix> : public true_type{};
	template <> class
	is_convertible<SEXP,MatrixNs::dpoMatrix> : public true_type{};
	template <> class
	is_convertible<SEXP,MatrixNs::dsyMatrix> : public true_type{};
	template <> class
	is_convertible<SEXP,MatrixNs::dtrMatrix> : public true_type{};
    }
}

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

    class wrongS4 : public std::exception {
    public:
	wrongS4(const char* S4class, const char* meth,
		const char* desired)
	    : d_clz(strdup(S4class)), d_methNm(strdup(meth)), d_desired(strdup(desired)) {
	}
	inline const char* what() const throw(){
	    return Rcpp::sprintf<100>("Class %s object passed to %s is not a %s",
				      d_clz, d_methNm, d_desired).c_str();
	}
	~wrongS4() throw(){free(d_clz); free(d_methNm); free(d_desired);}
    private:
	char *d_clz, *d_methNm, *d_desired;
    };

    class Matrix {
    protected:
	SEXP           d_sexp;
	Rcpp::List d_dimnames;
	int    d_nrow, d_ncol;
    public:
	Matrix(Rcpp::S4&);
	Matrix(int,int);

	int nrow() const {return d_nrow;}
	int ncol() const {return d_ncol;}
	const SEXPREC* sexp() const {return d_sexp ? d_sexp : R_NilValue;}
    };

    class modelMatrix {
    protected:
	Rcpp::IntegerVector    d_assign;
	Rcpp::List          d_contrasts;
    public:
	modelMatrix(Rcpp::S4&);

	Rcpp::IntegerVector const& assign() const {return d_assign;}
	Rcpp::List const&       contrasts() const {return d_contrasts;}
    };

    class dMatrix : public Matrix {
    protected:
	Rcpp::NumericVector d_x;
    public:
	dMatrix(Rcpp::S4&);
	dMatrix(int,int,int);
	dMatrix(int,int,Rcpp::NumericVector);

	const Rcpp::NumericVector& x() const {return d_x;}
	Rcpp::NumericVector& X() {return d_x;}
	void setX(const Rcpp::NumericVector&);
	void setX(const Rcpp::NumericMatrix&);
    };

    class ddenseMatrix : public dMatrix {
    public:
	ddenseMatrix(Rcpp::S4&) throw (std::runtime_error);
	ddenseMatrix(int,int);
	ddenseMatrix(int,int,Rcpp::NumericVector);
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
	dgeMatrix(Rcpp::NumericMatrix);

	int  dmult(char,double,double,chmDn const&,chmDn&) const;
	void dgemv(char,double,Rcpp::NumericVector const&,
		   double,Rcpp::NumericVector&) const throw (std::runtime_error);
	void dgemv(char,double,Rcpp::NumericVector const&,
		   double,double*) const throw (std::runtime_error);
	void dgemm(char,char,double,dgeMatrix const&,
		   double,dgeMatrix&) const throw (std::runtime_error);
	operator SEXP() const;
   };

// Inherits from dgeMatrix, not ddenseMatrix as in the R class
// definition so that the dmult method is available.  Note: we need to
// use a method that is defined for both dense and sparse matrices,
// which, at present, means using the chmDn class for the vector arguments.

    class ddenseModelMatrix : public dgeMatrix, public modelMatrix {
    public:
	ddenseModelMatrix(Rcpp::S4);

	operator SEXP() const;
    };



    class dtrMatrix : public ddenseMatrix, public triangularMatrix {
    public:
	dtrMatrix(Rcpp::S4&);
	dtrMatrix(int,char='U',char='N');

	operator SEXP() const;
	void dtrtrs(char,double*,int = 1) const throw(std::runtime_error);
    };

    class dsyMatrix : public ddenseMatrix, public symmetricMatrix {
    public:
	dsyMatrix(Rcpp::S4&);
	dsyMatrix(int,char='U');

	operator SEXP() const;
	void dsyrk(dgeMatrix const&,double,double)
	    throw (std::runtime_error);
    };

    class dpoMatrix : public dsyMatrix {
    public:
	dpoMatrix(Rcpp::S4&);
	dpoMatrix(int,char='U');
	operator SEXP() const;
    };

    class Cholesky : public dtrMatrix {
    public:
	Cholesky(Rcpp::S4);
	Cholesky(dgeMatrix, char = 'U');
	Cholesky(int, char = 'U');

	Rcpp::NumericMatrix solve(int,const_CHM_DN) const;
	Rcpp::NumericMatrix solve(int,const Rcpp::NumericMatrix&) const;
	Rcpp::NumericMatrix solve(int,const Rcpp::NumericVector&) const;

	double logDet2() const;

	void dpotrs(Rcpp::NumericVector&) const throw (std::runtime_error);
	void dpotrs(std::vector<double>&) const throw (std::runtime_error);
	void dpotrs(double*, int = 1) const throw (std::runtime_error);

	void inPlaceSolve(int, Rcpp::NumericMatrix&) const throw (std::runtime_error);

	void update(const dpoMatrix&) throw (std::runtime_error); // chol(A)
	void update(const dgeMatrix&) throw (std::runtime_error); // chol(crossprod(X))
	void update(char,double,const dgeMatrix&,double,const dsyMatrix&)
	    throw (std::runtime_error);
	operator SEXP() const;
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

	int nr() const { return nrow; }
	int nc() const { return ncol; }
	double* begin() {return (double*)x;} // template this
	double* end() {return begin() + nrow * ncol;}
	operator SEXP() const;
    };

    class chmSp : public cholmod_sparse { 
	Rcpp::S4    d_xp;
    public:
	chmSp(Rcpp::S4) throw(wrongS4);

	int              nr() const {return nrow;}
	int              nc() const {return ncol;}
	int             nnz() const {return ((int*)p)[ncol];}
	int           nzMax() const {return nzmax;}

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
	void update(const cholmod_sparse&) throw(std::runtime_error);
	operator SEXP() const throw (std::runtime_error);
    };

    class chmFr : public cholmod_factor {
	Rcpp::S4     d_xp;
    public:
	chmFr(Rcpp::S4) throw (wrongS4);

	double      logDet2() const;

	void update(const cholmod_sparse&, double Imult = 0.); 

	Rcpp::NumericMatrix solve(int, const_CHM_DN) const;
	Rcpp::NumericMatrix solve(int, const Rcpp::NumericMatrix&) const;
	Rcpp::NumericMatrix solve(int, const Rcpp::NumericVector&) const;

	CHM_SP spsolve(int sys, const_CHM_SP b) const;
	CHM_SP spsolve(int sys, const chmSp&b) const;
	operator SEXP() const throw (std::runtime_error);
    };
#if 0
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
#endif
}

namespace Rcpp {
    template <>
    SEXP wrap<MatrixNs::dgeMatrix>(const MatrixNs::dgeMatrix& m);
    template <>
    SEXP wrap<MatrixNs::dtrMatrix>(const MatrixNs::dtrMatrix& m);
    template <>
    SEXP wrap<MatrixNs::dsyMatrix>(const MatrixNs::dsyMatrix& m);
    template <>
    SEXP wrap<MatrixNs::dpoMatrix>(const MatrixNs::dpoMatrix& m);
    template <>
    SEXP wrap<MatrixNs::Cholesky>(const MatrixNs::Cholesky& m);
    template <>
    SEXP wrap<MatrixNs::chmFr>(const MatrixNs::chmFr& m);
    template <>
    SEXP wrap<MatrixNs::chmSp>(const MatrixNs::chmSp& m);
}

#endif
