#include "lme4utils.h"		// for CHM_DN, CHM_SP, c

// Classes of double precision matrices

/** abstract class of double precision matrices in CHOLMOD form */
class CHM_r {
public:
    virtual void drmult(int transpose, double alpha, double beta,
			CHM_DN X, CHM_DN Y) = 0;
    virtual int ncol() = 0;
    virtual int nrow() = 0;
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
class Cholesky_r {
public:
    virtual void update(CHM_r *A) = 0;
    virtual int ncol() = 0;
    virtual int nrow() = 0;
};

class Cholesky_rd : public Cholesky_r {
public:
    Cholesky_rd(SEXP x);
    ~Cholesky_rd(){delete F;}
    virtual void update(CHM_r *A);
    virtual int ncol() {return F->ncol;}
    virtual int nrow() {return F->nrow;}
    CHM_DN F;
    int upper;			/**< Boolean indicator of upper factor */
};

class Cholesky_rs : public Cholesky_r {
public:
    Cholesky_rs(SEXP x);
    ~Cholesky_rs(){delete F;}
    virtual void update(CHM_r *A);
    // virtual int ncol() {return F->ncol;}
    // virtual int nrow() {return F->nrow;}
    virtual int ncol() {return 0;}
    virtual int nrow() {return 0;}
    CHM_FR F;
};

