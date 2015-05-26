/*
 * tests.cpp
 *
 *  Created on: May 5, 2015
 *      Author: dmarce1
 */

#include "../RKDGFMM.hpp"

extern "C" {
void riemann_(double* rhol, double* pl, double* ul, double* rhor, double* pr, double* ur, double* gamma, double* t,
		double* x, double* rho, double* u, double* p);
}

std::vector<real> sod_shock_tube_analytic(real _x, real y, real z, real _t) {
	double rho, egas, s;
	double rho_l = sod_rho_l;
	double rho_r = sod_rho_r;
	double p_l = sod_p_l;
	double p_r = sod_p_r;
	double u_l = real(0);
	double u_r = real(0);
	double gamma = fgamma;
	double u, p;
	double t = _t;
	double x = _x;
	riemann_(&rho_l, &p_l, &u_l, &rho_r, &p_r, &u_r, &gamma, &t, &x, &rho, &u, &p);
	s = rho * u;
	egas = p / (gamma - real(1)) + s * u / real(2);
	std::vector < real > U(NF);
	U[rho_i] = rho;
	U[s_i + XDIM] = s;
	U[s_i + YDIM] = U[s_i + ZDIM] = real(0);
	U[egas_i] = egas;
	U[tau_i] = std::pow(egas, real(1) / gamma);
	return U;

}
