// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
// Can't call this Matrix.h because of the Matrix.h in Matrix/include
#ifndef LME4_MATRIX_H
#define LME4_MATRIX_H

#include <Rcpp.h>
#include <Matrix.h>

namespace Matrix {

    class Matrix {
    public:
	Matrix(Rcpp::S4);
	
	Rcpp::IntegerVector Dim;
	Rcpp::List Dimnames;
    };

    class dMatrix : public Matrix {
    public:
	dMatrix(Rcpp::S4);
	
	Rcpp::NumericVector x;
    };

    class ddenseMatrix : public dMatrix {
    public:
	ddenseMatrix(Rcpp::S4);
    };

    class dgeMatrix : public ddenseMatrix {
    public:
	dgeMatrix(Rcpp::S4);

	Rcpp::List factors;
    };

    class triangularMatrix :  public Matrix {
    public:
	triangularMatrix(Rcpp::S4);

	Rcpp::CharacterVector uplo, diag;
    };
    
	
    class dtrMatrix : public ddenseMatrix, triangularMatrix {
    public:
	dtrMatrix(Rcpp::S4);
    };

    class Cholesky : public dtrMatrix {
    public:
	Cholesky(Rcpp::S4);
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

    class dgCMatrix {
    public:
	dgCMatrix(Rcpp::S4);
	~dgCMatrix(){delete sp;}
	void update(CHM_SP);
	
	Rcpp::IntegerVector i, p, Dim;
	Rcpp::List Dimnames, factors;
	Rcpp::NumericVector x;
	CHM_SP sp;
    };
}

#endif
