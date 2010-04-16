// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-

// Can't name the file Matrix.h because of the Matrix.h in Matrix/include

#ifndef LME4_MATRIX_H
#define LME4_MATRIX_H

#include <Rcpp.h>
#include <Matrix.h>

namespace MatrixNs {

    class Matrix {
    public:
	Matrix(Rcpp::S4&);
	int nrow(){return Dim[0];}
	int ncol(){return Dim[1];}
	
	Rcpp::Dimension Dim;
	Rcpp::List Dimnames;
    };

    class dMatrix : public Matrix {
    public:
	dMatrix(Rcpp::S4&);
	
	Rcpp::NumericVector x;
    };

    class ddenseMatrix : public dMatrix {
    public:
	ddenseMatrix(Rcpp::S4&);
    };

// C++ classes mirroring virtual S4 structure classes do not inherit
// from Matrix, so as to avoid multiple definitions of Dim and
// Dimnames slots.  You can get around this for concrete but classes
// these are pure virtual classes.

    class compMatrix {		// composite (factorizable) Matrix
    public:
	compMatrix(Rcpp::S4&);
	
	Rcpp::List factors;
    };

    class generalMatrix : public compMatrix {
    public:
	generalMatrix(Rcpp::S4&);
    };

    class triangularMatrix {
    public:
	triangularMatrix(Rcpp::S4&);

	Rcpp::CharacterVector uplo, diag;
    };
    
    class symmetricMatrix : public compMatrix {
    public:
	symmetricMatrix(Rcpp::S4&);

	Rcpp::CharacterVector uplo;
    };
	
    class dgeMatrix : public ddenseMatrix, public generalMatrix {
    public:
	dgeMatrix(Rcpp::S4 xp) : ddenseMatrix(xp), generalMatrix(xp) {}

	void dgemv(const char*,double,const Rcpp::NumericVector&,
		   double,Rcpp::NumericVector&);
	void dgemm(const char*,const char*,double,const dgeMatrix&,
		   double,dgeMatrix&);
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
	Cholesky(Rcpp::S4 xp) : dtrMatrix(xp) {}
    };

    class chmDn : public cholmod_dense {
    public:
	chmDn(double*, int, int);
	chmDn(std::vector<double>);
	chmDn(Rcpp::NumericVector);
	chmDn(Rcpp::NumericMatrix);
	chmDn(ddenseMatrix);

	int nr() const { return nrow; }
	int nc() const { return ncol; }

    private:
	void init(double*, int, int);
    };

    class chmSp : public cholmod_sparse { 
    public:
	chmSp(Rcpp::S4);
	void update(CHM_SP);
    };

    class chmFa : public cholmod_factor {
    public:
	chmFa(Rcpp::S4);
	void update(chmSp, double);
    };

    class dCHMfactor {		// wrapper for CHM_FR struct
    public:
	dCHMfactor(Rcpp::S4);
	~dCHMfactor(){delete fa;}
	void update(CHM_SP,double);
	
	CHM_FR fa;
                               // slots common to dCHMsimpl and dCHMsuper
	Rcpp::IntegerVector Dim, ColCount, perm, type; 
	Rcpp::NumericVector x;
    };

}




#endif
