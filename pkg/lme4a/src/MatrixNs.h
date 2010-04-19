// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-

// Can't name the file Matrix.h because of the Matrix.h in Matrix/include

#ifndef LME4_MATRIX_H
#define LME4_MATRIX_H

#include <Rcpp.h>
#include <Matrix.h>

extern cholmod_common c;

namespace MatrixNs {

    class Matrix {
    public:
	Rcpp::Dimension Dim;
	Rcpp::List Dimnames;

	Matrix(Rcpp::S4 &xp) :
	    Dim(SEXP(xp.slot("Dim"))),
	    Dimnames(SEXP(xp.slot("Dimnames"))) { }

	int nrow(){return Dim[0];}
	int ncol(){return Dim[1];}
	
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
    public:
	Rcpp::List factors;
	compMatrix(Rcpp::S4 &xp) : factors(SEXP(xp.slot("factors"))) { }
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
	
    class dgeMatrix : public ddenseMatrix, public generalMatrix {
    public:
	dgeMatrix(Rcpp::S4 xp) : ddenseMatrix(xp), generalMatrix(xp) {}

	void dgemv(Trans,double,const Rcpp::NumericVector&,
		   double,Rcpp::NumericVector&);
	void dgemv(char Tr, double alpha, const Rcpp::NumericVector &X,
		   double beta, Rcpp::NumericVector &Y) {
	    dgemv(Trans(Tr), alpha, X, beta, Y);
	}


	void dgemm(Trans,Trans,double,const dgeMatrix&,double,dgeMatrix&);
	void dgemm(char TrA, char TrB,
		   double alpha, const dgeMatrix &B,
		   double beta, dgeMatrix &C) {
	    dgemm(Trans(TrA), Trans(TrB), alpha, B, beta, C);
	}
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

	void update(const dpoMatrix&);

	void update(Trans,double,const dgeMatrix&,double,const dsyMatrix&);
	void update(char Tr, double alpha, const dgeMatrix& A,
		    double beta, const dsyMatrix& C) {
	    update(Trans(Tr), alpha, A, beta, C);
	}

	void dpotrs(Rcpp::NumericVector&);

	double logDet2();
    };

    class chmDn : public cholmod_dense {
    private:
	void init(double*, int, int);
    public:
	chmDn(double *xx, int nr, int nc) : cholmod_dense() {
	    this->init(xx, nr, nc);
	}
	chmDn(std::vector<double> v) : cholmod_dense() {
	    this->init(&v[0], v.size(), 1);
	}
	chmDn(Rcpp::NumericVector v) : cholmod_dense() {
	    this->init(v.begin(), v.size(), 1);
	}
	chmDn(Rcpp::NumericMatrix m) : cholmod_dense() {
	    this->init(m.begin(), m.nrow(), m.ncol());
	}
    	chmDn(ddenseMatrix m) : cholmod_dense() {
	    this->init(m.x.begin(), m.nrow(), m.ncol());
	}

	int nr() const { return nrow; }
	int nc() const { return ncol; }
	double* begin() {return (double*)x;} // template this
	double* end() {return begin() + nrow * ncol;}
    };

    class chmSp : public cholmod_sparse { 
    public:
	chmSp(Rcpp::S4);

	void update(cholmod_sparse&);
	CHM_SP transpose(int values = 1) {
	    return ::M_cholmod_transpose(this, values, &c);
	}
	int dmult(char tr, double alpha, double beta,
		  chmDn &src, chmDn &dest) {
	    return ::M_cholmod_sdmult(this, Trans(tr).TR == 'T', &alpha,
				      &beta, &src, &dest, &c);
	}
	CHM_SP aat() {
	    return ::M_cholmod_aat(this, (int*)NULL, 0/*fsize*/, 1/*mode*/, &c);
	}
    };

    class chmFa : public cholmod_factor {
    public:
	chmFa(Rcpp::S4);

	void update(cholmod_sparse &A, double Imult = 0.) {
	    double ImVec[] = {Imult, 0};
	    ::M_cholmod_factorize_p(&A, ImVec, (int*)NULL, 0, this, &c);
	}
	CHM_DN solve(int sys, CHM_DN b) {
	    return ::M_cholmod_solve(sys, this, b, &c);
	}
	CHM_SP spsolve(int sys, CHM_SP b) {
	    return ::M_cholmod_spsolve(sys, this, b, &c);
	}
    };
}

#endif
