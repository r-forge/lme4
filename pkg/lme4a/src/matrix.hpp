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
    virtual void freeA() = 0;
    virtual void drmult(int transpose, double alpha, double beta,
			CHM_DN X, CHM_DN Y) = 0;
    virtual CHM_r* solveCHM_FR(CHM_FR L, int sys) = 0;
    virtual void copy_contents(CHM_r *src) = 0;
// FIXME: Change the mernew structure to bypass the special case of
// diagonal Lambda and eliminate the need to pass the double*
    virtual CHM_r* crossprod_SP(CHM_SP) = 0;
};

/** concrete class of double precision dense matrices */
class CHM_rd : public CHM_r {
public:
    CHM_rd(SEXP x);
    CHM_rd(CHM_DN x){A = x;}
    ~CHM_rd(){delete A;}
    virtual void freeA(){M_cholmod_free_dense(&A, &c);}
    virtual void drmult(int transpose, double alpha, double beta,
			CHM_DN X, CHM_DN Y);
    virtual int ncol() {return A->ncol;}
    virtual int nrow() {return A->nrow;}
    virtual CHM_r* solveCHM_FR(CHM_FR L, int sys);
    virtual void copy_contents(CHM_r *src);
    virtual CHM_r* crossprod_SP(CHM_SP);
	
    CHM_DN A;
};

/** concrete class of double precision sparse matrices */
class CHM_rs : public CHM_r {
public:
    CHM_rs(SEXP x);
    CHM_rs(CHM_SP x){A = x;}
    ~CHM_rs(){delete A;}
    virtual void freeA(){M_cholmod_free_sparse(&A, &c);}
    virtual void drmult(int transpose, double alpha, double beta,
			CHM_DN X, CHM_DN Y);
    virtual int ncol() {return A->ncol;}
    virtual int nrow() {return A->nrow;}
    virtual CHM_r* solveCHM_FR(CHM_FR, int);
    virtual void copy_contents(CHM_r*);
    virtual CHM_r* crossprod_SP(CHM_SP);

    CHM_SP A;
};

/** abstract class of double precision Cholesky decompositions */
class Cholesky_r : public dMatrix {
public:
    virtual int update(CHM_r*) = 0;
    virtual void downdate(CHM_r*, double, dMatrix*, double) = 0;
    virtual CHM_DN solveA(CHM_DN) = 0;
};

class Cholesky_rd : public Cholesky_r {
public:
    Cholesky_rd(SEXP);
    virtual int update(CHM_r*);
    virtual void downdate(CHM_r*, double, dMatrix*, double);
    virtual int ncol() {return n;}
    virtual int nrow() {return n;}
    virtual CHM_DN solveA(CHM_DN);
    const char* uplo;
    int n;

    double *X;
};

class Cholesky_rs : public Cholesky_r {
public:
    Cholesky_rs(SEXP);
    ~Cholesky_rs(){delete F;}
    virtual int update(CHM_r*);
    virtual void downdate(CHM_r*, double, dMatrix*, double);
    virtual CHM_DN solveA(CHM_DN);
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

