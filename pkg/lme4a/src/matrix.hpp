// Classes of double precision matrices

/** abstract class of double precision matrices */
class dMatrix {
public:
    double *contents;
    virtual int nrow() = 0;
    virtual int ncol() = 0;
    virtual void rdprod(double *ans, const double* rhs,
		int nrhs, double alpha, double beta,
		int transa, int transb) = 0;
    static const int i1 = 1, i0 = 0;
    static const double d1, d0, dm1;
};

const double dMatrix::d1 = 1;
const double dMatrix::d0 = 0;
const double dMatrix::dm1 = -1;

/** concrete class of double precision diagonal square matrices */
class ddiMatrix : public dMatrix {
public:
    ddiMatrix(double *x, int nr, int nc);
    ddiMatrix(double *x, int *dims);
    virtual int nrow() {return nro;}
    virtual int ncol() {return nco;}
    virtual void rdprod(double *ans, const double* rhs,
		int nrhs, double alpha, double beta,
		int transa);
//    void lprod(double *ans, const double* lhs);
protected:
    int nro;
    int nco;
};

/** concrete class of double precision general matrices */
class dgeMatrix : public dMatrix {
public:
    dgeMatrix(double *x, int nr, int nc);
    dgeMatrix(double *x, int *dims);
    virtual int nrow() {return nro;}
    virtual int ncol() {return nco;}
    virtual void rdprod(double *ans, const double* rhs,
		int nrhs, double alpha, double beta,
		int transa);
protected:
    int nro;
    int nco;
};

/** concrete class of double precision, sparse, compressed-column matrices */
class dgCMatrix : public dMatrix {
    dgCMatrix(CHM_SP A);
    CHM_SP mm;
    virtual int nrow() {return (int)mm->nrow;}
    virtual int ncol(){return (int)mm->ncol;}
    virtual void rdprod(double *ans, const double* rhs,
		int nrhs, double alpha, double beta,
		int transa);
};
