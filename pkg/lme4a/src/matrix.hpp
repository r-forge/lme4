#include "lme4utils.h"		// for CHM_DN, CHM_SP, c

// Classes of double precision matrices

/** abstract class of double precision matrices in CHOLMOD form */
class dMatrix {
public:
    virtual int ncol() = 0;
    virtual int nrow() = 0;
};

class CHM_r : public dMatrix {
public:    
    virtual void drmult(int transpose, double alpha, double beta,
			CHM_DN X, CHM_DN Y) = 0;
};

/** concrete class of double precision dense matrices */
class CHM_rd : public CHM_r {
public:
    CHM_rd(SEXP x);
    ~CHM_rd(){delete A;}
    virtual void drmult(int transpose, double alpha, double beta,
			CHM_DN X, CHM_DN Y);
    virtual int ncol() {return A->ncol;}
    virtual int nrow() {return A->nrow;}
    CHM_DN A;
};

/** concrete class of double precision general matrices */
class CHM_rs : public CHM_r {
public:
    CHM_rs(SEXP x);
    ~CHM_rs(){delete A;}
    virtual void drmult(int transpose, double alpha, double beta,
			CHM_DN X, CHM_DN Y);
    virtual int ncol() {return A->ncol;}
    virtual int nrow() {return A->nrow;}
    CHM_SP A;
};

/** abstract class of double precision Cholesky decompositions */
class Cholesky_r : public dMatrix {
public:
    virtual int update(CHM_r *A) = 0;
};

class Cholesky_rd : public Cholesky_r {
public:
    Cholesky_rd(SEXP x);
    virtual int update(CHM_r *A);
    virtual int ncol() {return n;}
    virtual int nrow() {return n;}
    const char* uplo;
    int n;
    double *X;
};

class Cholesky_rs : public Cholesky_r {
public:
    Cholesky_rs(SEXP x);
    ~Cholesky_rs(){delete F;}
    virtual int update(CHM_r *A);
    virtual int ncol() {return (int)F->n;}
    virtual int nrow() {return (int)F->n;}
    CHM_FR F;
};

class dpoMatrix : public dMatrix {
public:
    dpoMatrix(SEXP x);
    virtual int ncol() {return n;}
    virtual int nrow() {return n;}
    const char* uplo;
    int n;
    double *X;
    SEXP factors;
};

