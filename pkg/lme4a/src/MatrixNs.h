// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-

// Can't name the file Matrix.h because of the Matrix.h in Matrix/include

#ifndef LME4_MATRIX_H
#define LME4_MATRIX_H

#include <Rcpp.h>
#include <Matrix.h>

extern cholmod_common c;

namespace MatrixNs {

    class Matrix {
    protected:
	Rcpp::List Dimnames;
	int d_nrow, d_ncol;
    public:
	Matrix(Rcpp::S4 &xp);
	int nrow() const;
	int ncol() const;
    };

    class dMatrix : public Matrix {
    public:
	Rcpp::NumericVector x;
	dMatrix(Rcpp::S4 &xp) : Matrix(xp), x(SEXP(xp.slot("x"))) { }
    };

    class ddenseMatrix : public dMatrix {
    public:
	ddenseMatrix(Rcpp::S4&);
    };

// C++ classes mirroring virtual S4 structure classes do not inherit
// from Matrix, so as to avoid multiple definitions of Dim and
// Dimnames slots.  You can get around this for concrete classes but
// these are pure virtual classes.

    class compMatrix {		//< composite (factorizable) Matrix
    protected:
	Rcpp::List factors;
    public:
	compMatrix(Rcpp::S4 &xp) :
	    factors(SEXP(xp.slot("factors"))) {
	}
    };

    class generalMatrix : public compMatrix { //< general structure
    public:
	generalMatrix(Rcpp::S4 &xp) : compMatrix(xp) { }
    };

    class UpLo {		//< validated uplo char
    private: 
	char chk (char x) {
	    if (x == 'L' || x == 'l') return 'L';
	    if (x == 'U' || x == 'u') return 'U';
	    throw std::range_error("uplo");
	}
    public:
	char UL;
	UpLo(char &x) : UL(chk(x)) {}
	UpLo(const std::string &x) : UL(chk(x[0])) {}
	UpLo(const SEXP x) : UL(chk(*CHAR(Rf_asChar(x)))) {}
    };

    class Diag {		//< validated diag char
    private: 
	char chk (char x) {
	    if (x == 'N' || x == 'n') return 'N';
	    if (x == 'U' || x == 'u') return 'U';
	    throw std::range_error("diag");
	}
    public:
	char DD;
	Diag(char x) : DD(chk(x)) {}
	Diag(const std::string &x) : DD(chk(x[0])) {}
	Diag(const SEXP x) : DD(chk(*CHAR(Rf_asChar(x)))) {};
    };

    class Trans {		//< validated trans char
    private: 
	char chk (char x) {
	    if (x == 'C' || x == 'c') return 'C';
	    if (x == 'N' || x == 'n') return 'N';
	    if (x == 'T' || x == 't') return 'T';
	    throw std::range_error("trans");
	}
    public:
	char TR;
	Trans(char x) : TR(chk(x)) {}
	Trans(const std::string &x) : TR(chk(x[0])) {}
	Trans(const SEXP x) : TR(chk(*CHAR(Rf_asChar(x)))) {};
    };

    class triangularMatrix {
    public:
	UpLo uplo;
	Diag diag;
	triangularMatrix(Rcpp::S4 &xp) :
	    uplo(SEXP(xp.slot("uplo"))),
	    diag(SEXP(xp.slot("diag"))) { }
    };
    
    class symmetricMatrix : public compMatrix {
    public:
	UpLo uplo;
	symmetricMatrix(Rcpp::S4 &xp) :
	    compMatrix(xp),
	    uplo(SEXP(xp.slot("uplo"))) { }
    };

// Concrete classes are initialized from Rcpp::S4 not Rcpp::S4& so
// they can be constructed from slots extracted during the
// construction of composite objects.

    class chmDn; 		// forward declaration

    class dgeMatrix : public ddenseMatrix, public generalMatrix {
    public:
	dgeMatrix(Rcpp::S4 xp) : ddenseMatrix(xp), generalMatrix(xp) {}
	int dmult(char tr, double alpha, double beta, 
	 	  chmDn const &src, chmDn &dest) const;
	void dgemv(char,double,Rcpp::NumericVector const&,
		   double, Rcpp::NumericVector&) const;
	void dgemv(char,double,Rcpp::NumericVector const&,
		   double,double*) const;
	void dgemm(char,char,double,dgeMatrix const&,
		   double,dgeMatrix&) const;
    };

    class dtrMatrix : public ddenseMatrix, public triangularMatrix {
    public:
	dtrMatrix(Rcpp::S4 xp) : ddenseMatrix(xp), triangularMatrix(xp) {}
    };

    class dsyMatrix : public ddenseMatrix, public symmetricMatrix {
    public:
	dsyMatrix(Rcpp::S4 xp) : ddenseMatrix(xp), symmetricMatrix(xp) {}

	void dsyrk(dgeMatrix&,double,double);
    };

    class dpoMatrix : public dsyMatrix {
    public:
	dpoMatrix(Rcpp::S4 xp) : dsyMatrix(xp) {}
    };

    class Cholesky : public dtrMatrix {
    public:
	Cholesky(Rcpp::S4 xp) : dtrMatrix(xp) {diag = Diag('N');}

	Rcpp::NumericMatrix solve(int, const_CHM_DN) const;
	Rcpp::NumericMatrix solve(int,Rcpp::NumericMatrix const&) const;
	Rcpp::NumericMatrix solve(int,Rcpp::NumericVector const&) const;

	void update(dpoMatrix const&); // chol(A)
	void update(dgeMatrix const&); // chol(crossprod(X))
	void update(Trans,double,dgeMatrix const&,double,
		    const dsyMatrix&);
	void update(char Tr, double alpha, dgeMatrix const& A,
		    double beta, dsyMatrix const& C) {
	    update(Trans(Tr), alpha, A, beta, C);
	}

	void dpotrs(Rcpp::NumericVector&) const;
	void dpotrs(std::vector<double>&) const;
	void dpotrs(double*, int = 1) const;

	double logDet2();
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
    public:
	chmSp(Rcpp::S4);

	CHM_SP crossprod() const;
	CHM_SP crossprod(const cholmod_sparse*, int sorted = 1) const;
	CHM_SP crossprod(chmSp const &B, int sorted = 1) const;

	CHM_SP tcrossprod() const;
	CHM_SP tcrossprod(const_CHM_SP, int sorted = 1) const;
	CHM_SP tcrossprod(chmSp const &B, int sorted = 1) const;

	CHM_SP transpose(int values = 1) const;

	CHM_SP smult(chmSp const&, int, int, int) const;
	int dmult(char, double, double, chmDn const&, chmDn&) const;
	
	void scale(int,chmDn const&);
	void update(cholmod_sparse const&);
    };

    class chmFr : public cholmod_factor {
//	CHM_FR pp;
//	SEXP m_sexp;
    public:
	chmFr(Rcpp::S4);
//	chmFr(CHM_FR);
//	~chmFr() {if (pp) ::M_cholmod_free_factor(&pp, &c);}

	void update(cholmod_sparse const &A, double Imult = 0.); 

	Rcpp::NumericMatrix solve(int, const_CHM_DN) const;
	Rcpp::NumericMatrix solve(int,Rcpp::NumericMatrix const&) const;
	Rcpp::NumericMatrix solve(int,Rcpp::NumericVector const&) const;

	CHM_SP spsolve(int sys, const_CHM_SP b) const;
	CHM_SP spsolve(int sys, chmSp const &b) const;

    };
}

#endif
