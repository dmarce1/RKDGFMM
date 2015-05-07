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

void sod_shock_tube_analytic(real x, real t, real* rho, real* s, real* egas) {
	double rho_l = real(1);
	double rho_r = real(1) / real(8);
	double p_l = real(1);
	double p_r = real(1) / real(10);
	double u_l = real(0);
	double u_r = real(0);
	double gamma = real(7) / real(5);
	double u, p;
	riemann_(&rho_l, &p_l, &u_l, &rho_r, &p_r, &u_r, &gamma, &t, &x, rho, &u, &p);
	*s = *rho * u;
	*egas = p / (gamma - real(1)) + *s * u / real(2);

}
