/*
 * initial.cpp
 *
 *  Created on: Apr 24, 2015
 *      Author: dmarce1
 */

#include "initial.hpp"

std::vector<real> sod_shock_tube(real x, real y, real z) {
	std::vector<real> u(NF, real(0));
	if (x < real(0)) {
		u[rho_i] = real(1);
		u[egas_i] = real(2.5);
	} else {
		u[rho_i] = real(.125);
		u[egas_i] = real(0.1);
	}
	u[tau_i] = std::pow(real(u[egas_i]), (real(1) / real(fgamma)));
	return u;
}

std::vector<real> blast_wave(real x, real y, real z) {
	std::vector<real> u(NF, real(0));
	u[rho_i] = real(1);
	if (std::sqrt(x * x + y * y + z * z) < 4.0e-2) {
		u[egas_i] = real(1);
	} else {
		u[egas_i] = real(2.5e-5);
	}
	u[tau_i] = std::pow(real(u[egas_i]), (real(1) / real(fgamma)));
	return u;
}
