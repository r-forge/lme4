// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
#ifndef LME4_GHQ_H
#define LME4_GHQ_H

#include <vector>

namespace GaussQuad {
    enum Rule{Legendre=1, Chebyshev, Gegenbauer, Jacobi, Laguerre,
	      Hermite, Exponential, Rational, Type2};
    class GQ {
	int                 d_n;
	Rule                d_rule;
	double              d_a, d_b, d_alpha, d_beta, d_zemu;
	std::vector<double> d_t, d_wts, d_diag, d_sub;
    public:
	GQ(int, Rule, double, double, double, double);

	const std::vector<double>& diag() const {return d_diag;}
	const std::vector<double>&  sub() const {return d_sub;}
	const std::vector<double>&    t() const {return d_t;}
	const std::vector<double>&  wts() const {return d_wts;}	

	double                        a() const {return d_a;}
	double                        b() const {return d_b;}
	double                    alpha() const {return d_alpha;}
	double                     beta() const {return d_beta;}
	double                     zemu() const {return d_zemu;}

	int                           n() const {return d_n;}
	int                        rule() const {return int(d_rule);}

	void               class_matrix();
	void                     imtqlx(std::vector<double>&);
	void                       scqf(const std::vector<int>&, const std::vector<int>&);
	void                       sgqf();
	void                     parchk(int);
    };
}

#endif
