#include "GQ.h"
#include <limits>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <iostream>

namespace GaussQuad {
    inline double r8_sign(const double& x) {return (x < 0) ? -1. : 1.;}

    GQ::GQ(int n, Rule r, double a, double b, double alpha, double beta)
	: d_n(n), d_rule(r), d_a(a), d_b(b), d_alpha(alpha), d_beta(beta),
	  d_t(n), d_wts(n), d_diag(n), d_sub(n) {
	std::vector<int> mlt(d_n), ndx(d_n);
				// Compute the Gauss quadrature formula for default a and b.
				//  Get the Jacobi matrix and zero-th moment.
	class_matrix();
				//  Compute the knots and weights.
	sgqf();
				// Prepare to scale the quadrature formula to other weight
				// function with valid A and B.
	std::fill(mlt.begin(), mlt.end(), 1);
	for (int i = 0; i < d_n; i++) ndx[i] = i + 1;
	scqf(mlt, ndx);
    }

/** 
 * parchk checks parameters alpha and beta for classical weight functions. 
 * 
 * @param m the order of the highest moment to be calculated
 */
    void GQ::parchk(int m) {
            //  Check alpha for Gegenbauer, Jacobi, Laguerre, Hermite, Exponential.
	if (Gegenbauer <= d_rule && d_alpha <= -1.0)
	    throw std::runtime_error("parchk:  3 <= kind and alpha <= -1.");
				//  check beta for Jacobi.
	if (d_rule == Jacobi && d_beta <= -1.0)
	    throw std::runtime_error("parchk: Jacobi rule and beta <= -1.0");
				//  check alpha and beta for rational.
	if (d_rule == Rational) {
	    double tmp = d_alpha + d_beta + double(m) + 1.0;
	    if (0.0 <= tmp || tmp <= d_beta)
		throw std::runtime_error("parchk: Rational rule condition on alpha and beta fails");
	}
    }

/** 
 * Compute the Jacobi matrix for a quadrature rule.
 * 
 * This routine computes the diagonal and sub-diagonal
 * elements of the order d_n tridiagonal symmetric Jacobi matrix
 * associated with the polynomials orthogonal with respect to
 * the weight function specified by d_rule.
 *
 * For weight functions other than Rational, d_n elements are defined in d_sub even
 * though only d_n-1 are needed.  For the Rational weight function, d_sub[d_n-1] is
 * set to zero.
 *
 * Sets the zero-th moment of the weight function, d_zemu.
 */
    void GQ::class_matrix() {
	const double pi = 3.14159265358979323846264338327950;
	double a2b2 = 0., ab = 0., abi = 0, abp2 = 0., apone = d_alpha + 1.;

	parchk(2 * d_n - 1);
	d_zemu = 0.;

	if (500.0 *  std::numeric_limits<double>::epsilon() < std::abs(std::pow(tgamma(0.5), 2) - pi))
	    throw std::runtime_error("Gamma function does not match machine parameters.");
	switch(d_rule) {
	case Legendre:
//	    ab = 0.0;
//	    d_zemu = 2.0 / (ab + 1.0);
	    d_zemu = 2.;
	    std::fill(d_diag.begin(), d_diag.end(), double());
	    for (int i = 1; i <= d_n; i++) {
		double abi = i;// + ab * (i % 2);
		double abj = 2 * i;// + ab;
		d_sub[i-1] = std::sqrt(abi * abi / (abj * abj - 1.0));
	    }
	    break;
	case Chebyshev:
	    d_zemu = pi;
	    std::fill(d_diag.begin(), d_diag.end(), double());
	    std::fill(d_sub.begin(), d_sub.end(), 0.5);
	    d_sub[0] = std::sqrt(0.5);
	    break;
	case Gegenbauer:
	    ab = d_alpha * 2.0;
	    d_zemu = std::pow(2.0, ab + 1.0) * std::pow(tgamma(d_alpha + 1.0), 2)
		/ tgamma (ab + 2.0);
	    std::fill(d_diag.begin(), d_diag.end(), double());
	    d_sub[0] = std::sqrt(1.0 / (2.0 * d_alpha + 3.0));
	    for(int i = 2; i <= d_n; i++) {
		d_sub[i-1] = std::sqrt(i * (i + ab) / (4.0 * std::pow(i + d_alpha, 2) - 1.0));
	    }
	    break;
	case Jacobi:
	    ab = d_alpha + d_beta;
	    abp2 = 2.0 + ab;
	    d_zemu = std::pow(2.0, ab + 1.0) * tgamma (d_alpha + 1.0) 
		* tgamma (d_beta + 1.0) / tgamma (abp2);
	    d_diag[0] = (d_beta - d_alpha) / abp2;
	    d_sub[0] = std::sqrt(4.0 * (1.0 + d_alpha) * (1.0 + d_beta) 
			   / ((abi + 1.0) * abp2 * abp2));
	    a2b2 = d_beta * d_beta - d_alpha * d_alpha;
	    for(int i = 2; i <= d_n; i++) {
		double abi = 2.0 * i + ab;
		d_diag[i-1] = a2b2 / ((abi - 2.0) * abi);
		abi = abi * abi;
		d_sub[i-1] = std::sqrt(4.0 * i * (i + d_alpha) * (i + d_beta) * (i + ab) 
				       / ((abi - 1.0) * abi));
	    }
	    break;
	case Laguerre:
	    d_zemu = tgamma(d_alpha + 1.0);
	    for(int i = 1; i <= d_n; i++) {
		d_diag[i-1] = 2.0 * i - 1.0 + d_alpha;
		d_sub[i-1] = std::sqrt(i * (i + d_alpha));
	    }
	    break;
	case Hermite:
	    d_zemu = tgamma((d_alpha + 1.0) / 2.0);
	    std::fill(d_diag.begin(), d_diag.end(), double());
	    for(int i = 1; i <= d_n; i++) {
		d_sub[i-1] = std::sqrt((i + d_alpha * (i % 2)) / 2.);
	    }
	    break;
	case Exponential:
	    ab = d_alpha;
	    d_zemu = 2.0 / (ab + 1.0);
	    std::fill(d_diag.begin(), d_diag.end(), double());
	    for(int i = 1; i <= d_n; i++) {
		abi = i + ab * (i % 2);
		double abj = 2 * i + ab;
		d_sub[i-1] = std::sqrt(abi * abi / (abj * abj - 1.0));
	    }
	    break;
	case Rational:
	    ab = d_alpha + d_beta;
	    d_zemu = tgamma(d_alpha + 1.0) * tgamma(-(ab + 1.0)) / tgamma(-d_beta);
	    d_diag[0] = -apone / (ab + 2.0);
	    d_sub[0] = -d_diag[0] * (d_beta + 1.0) / (ab + 2.0) / (ab + 3.0);
	    for(int i = 2; i <= d_n; i++) {
		double abti = ab + 2.0 * i;
		d_diag[i-1] = ab * apone + 2.0 * (ab + i) * (i - 1);
		d_diag[i-1] = - d_diag[i-1] / abti / (abti - 2.0);
	    }

	    for(int i = 2; i <= d_n - 1; i++) {
		double abti = ab + 2.0 * i;
		d_sub[i-1] = i * (d_alpha + i) / (abti - 1.0) * (d_beta + i) 
		    / (abti * abti) * (ab + i) / (abti + 1.0);
	    }
	    d_sub[d_n - 1] = 0.0;
	    std::transform(d_sub.begin(), d_sub.end(), d_sub.begin(), ::sqrt);
	    break;
	case Type2:
	    // This case is not handled in John Burhardt's C++ code
	    break;
	}
    }

/** 
 * Diagonalize a symmetric tridiagonal matrix.
 * 
 * This routine is a slightly modified version of the EISPACK routine to 
 * perform the implicit QL algorithm on a symmetric tridiagonal matrix. 
 *
 * The authors thank the authors of EISPACK for permission to use this
 * routine. 
 *
 * Reference:
 *
 * Sylvan Elhay, Jaroslav Kautsky,
 * Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
 * Interpolatory Quadrature,
 * ACM Transactions on Mathematical Software,
 * Volume 13, Number 4, December 1987, pages 399-415.
 *
 * Roger Martin, James Wilkinson,
 * The Implicit QL Algorithm,
 * Numerische Mathematik,
 * Volume 12, Number 5, December 1968, pages 377-383.
 *
 * It has been modified to produce the product Q' * Z, where Z is an input 
 * vector and Q is the orthogonal matrix diagonalizing the input matrix.  
 * The changes consist (essentially) of applying the orthogonal transformations
 * directly to Z as they are generated.
 *
 * @param z the value of Q' * Z where Q is the matrix that
 *          diagonalizes the input symmetric tridiagonal matrix 
 */
    void GQ::imtqlx(std::vector<double>& z) {
	double b;
	double c;
	double f;
	double g;
	int i;
	int itn = 30;
	int m = 1;
	int mml;
	double p;
	double prec = std::numeric_limits<double>::epsilon();
	double r;
	double s;
	std::vector<double> bi(d_sub);

	if (d_n == 1) return;

	bi[d_n - 1] = 0.0;

	for (int l = 1; l <= d_n; l++) {
	    int j = 0;
	    for ( ; ;) {
		for (m = l; m <= d_n; m++) {
		    if (m == d_n) break; 

		    if (std::abs(bi[m-1]) <= prec * (std::abs(d_t[m-1]) + std::abs(d_t[m]))) break;
		}
		p = d_t[l-1];
		if (m == l) break;
		if (itn <= j) throw std::runtime_error("imtqlx: Iteration limit exceeded");
		j++;
		g = (d_t[l] - p) / (2.0 * bi[l-1]);
		r = std::sqrt(g * g + 1.0);
		g = d_t[m-1] - p + bi[l-1] / (g + std::abs(r) * r8_sign(g));
		s = 1.0;
		c = 1.0;
		p = 0.0;
		mml = m - l;
		
		for(int ii = 1; ii <= mml; ii++) {
		    i = m - ii;
		    f = s * bi[i-1];
		    b = c * bi[i-1];
		    
		    if (std::abs(g) <= std::abs(f)) {
			c = g / f;
			r = std::sqrt(c * c + 1.0);
			bi[i] = f * r;
			s = 1.0 / r;
			c = c * s;
		    } else {
			s = f / g;
			r = std::sqrt(s * s + 1.0);
			bi[i] = g * r;
			c = 1.0 / r;
			s = s * c;
		    }
		    g = d_t[i] - p;
		    r = (d_t[i-1] - g) * s + 2.0 * c * b;
		    p = s * r;
		    d_t[i] = g + p;
		    g = c * r - b;
		    f = z[i];
		    z[i] = s * z[i-1] + c * f;
		    z[i-1] = c * z[i-1] - s * f;
		}
		d_t[l-1] = d_t[l-1] - p;
		bi[l-1] = g;
		bi[m-1] = 0.0;
	    }
	}
				//  Sorting.
	for(int ii = 2; ii <= m; ii++) {
	    int i = ii - 1;
	    int k = i;
	    double p = d_t[i-1];
	    
	    for (int j = ii; j <= d_n; j++) {
		if (d_t[j-1] < p) {
		    k = j;
		    p = d_t[j-1];
		}
	    }
	    
	    if (k != i) {
		d_t[k-1] = d_t[i-1];
		d_t[i-1] = p;
		p = z[i-1];
		z[i-1] = z[k-1];
		z[k-1] = p;
	    }
	}
    }


/** 
 * Scale a quadrature formula to a nonstandard interval.
 * 
 * @param mlt the multiplicity of the knots
 * @param ndx 
 */
    void GQ::scqf (const std::vector<int>& mlt, const std::vector<int>& ndx) {
	double al = d_alpha;
	double be = d_beta;
	double shft = (d_a + d_b)/2.;
	double slp  = (d_b - d_a)/2.;

	switch (d_rule) {
	case Legendre:
	case Chebyshev:
	case Gegenbauer:
	case Jacobi:
	case Exponential:
	case Type2:
	    if (std::abs(d_b - d_a) <= std::numeric_limits<double>::epsilon())
		throw std::runtime_error("scqf: |B - A| too small");
	    break;		// unnecessary but it can't hurt
	case Laguerre:
	case Hermite:
	    if (d_b <= 0.0) throw std::runtime_error("scqf: b <= 0");
	    break;
	case Rational:
	    break;
	}

	switch (d_rule) {
	case Legendre:
	    al = 0.0;
	    be = 0.0;
	    break;
	case Chebyshev:
	    al = -0.5;
	    be = -0.5;
	    break;
	case Gegenbauer:
	    be = d_alpha;
	    break;
	case Jacobi:
	    break;
	case Laguerre:
	    shft = d_a;
	    slp = 1.0 / d_b;
	    be = 0.0;
	    break;
	case Hermite:
	    shft = d_a;
	    slp = 1.0 / std::sqrt(d_b);
	    be = 0.0;
	    break;
	case Exponential:
	    be = 0.0;
	    break;
	case Rational:
	    if (d_a + d_b <= 0.0) throw std::runtime_error("scqf: A + B <= 0.");
	    shft = d_a;
	    slp = d_a + d_b;
	    break;
	case Type2:
	    al = 0.5;
	    be = 0.5;
	    break;
	}

	double p = std::pow(slp, al + be + 1.0);

	for (int k = 0; k < d_n; k++) {
	    d_t[k] = shft + slp * d_t[k];
	    int l = std::abs(ndx[k]);
	    if (l != 0) {
		double tmp = p;
		for(int i = l - 1; i <= l - 1 + mlt[k] - 1; i++) {
		    d_wts[i] *=  tmp;
		    tmp = tmp * slp;
		}
	    }
	}
    }

/** 
 * Compute knots and weights of a Gauss Quadrature formula.
 * 
 * @param zemu the zero-th moment of the weight function
 */
    void GQ::sgqf() {
	if (d_zemu <= 0.0) throw std::runtime_error("sgqf: zemu <= 0.");
				//  Set up vectors for imtqlx.
	std::copy(d_diag.begin(), d_diag.end(), d_t.begin());
	std::fill(d_wts.begin(), d_wts.end(), 0.);
	d_wts[0] = std::sqrt(d_zemu);
				//  Diagonalize the Jacobi matrix.
	imtqlx(d_wts);
	for(int i = 0; i < d_n; i++) d_wts[i] *= d_wts[i];
    }
}

