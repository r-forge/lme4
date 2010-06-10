// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-

// Can't name the file Matrix.h because of the Matrix.h in Matrix/include

#ifndef LME4_MATRIX_H
#define LME4_MATRIX_H

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
	Rcpp::List Dimnames;
	int d_nrow, d_ncol;
    public:
	Matrix(Rcpp::S4&);
	Matrix(int,int);
	int nrow() const;
	int ncol() const;
    };

    class dMatrix : public Matrix {
    protected:
	Rcpp::NumericVector d_x;
    public:
	dMatrix(Rcpp::S4&);
	dMatrix(int,int,int=0);

	Rcpp::NumericVector const& x() const{return d_x;}
	Rcpp::NumericVector& X() {return d_x;}
	void setX(Rcpp::NumericVector const&);
	void setX(Rcpp::NumericMatrix const&);
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
	dtrMatrix(Rcpp::S4&);
	dtrMatrix(int,char='U',char='N');

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

	Rcpp::NumericMatrix solve(int, const_CHM_DN) const;
	Rcpp::NumericMatrix solve(int,Rcpp::NumericMatrix const&) const;
	Rcpp::NumericMatrix solve(int,Rcpp::NumericVector const&) const;

	void update(dpoMatrix const&); // chol(A)
	void update(dgeMatrix const&); // chol(crossprod(X))
	void update(char,double,dgeMatrix const&,
		    double, dsyMatrix const&);

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
