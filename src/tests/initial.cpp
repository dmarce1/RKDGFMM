/*
 * initial.cpp
 *
 *  Created on: Apr 24, 2015
 *      Author: dmarce1
 */

#include "initial.hpp"
#include "lane_emden.hpp"
#include <atomic>

const real x0 = 0.0;
const real y0_ = -0.0;
const real z0 = 0.0;
constexpr
real rmax = 3.7;
constexpr
real dr = rmax / 128.0;
const real alpha = real(1) / real(6);

std::vector<real> star(real x, real y, real z) {
	x -= x0;
	y -= y0_;
	z -= z0;
	real theta;
	const real n = real(1) / (fgamma - real(1));
	const real rho_min = 1.0e-4;
	std::vector < real > u(NF, real(0));
	const real r = std::sqrt(x * x + y * y + z * z) / alpha;
	const real theta_min = std::pow(rho_min, real(1) / n);
	const auto c0 = real(4) * real(M_PI) * alpha * alpha / (n + real(1));
	if (r <= rmax) {
		theta = lane_emden(r, dr);
		theta = std::max(theta, theta_min);
	} else {
		theta = theta_min;
	}
	u[rho_i] = std::pow(theta, n);
	u[egas_i] = std::pow(theta, fgamma * n) * c0 / (fgamma - real(1));
	u[tau_i] = std::pow(u[egas_i], (real(1) / real(fgamma)));
	return u;
}

void star_force(real x, real y, real z, real& fx, real& fy, real& fz) {
	x -= x0;
	y -= y0_;
	z -= z0;
	constexpr
	real mtot = 3.410854709920873e+01;
	real m;
	const real r = std::sqrt(x * x + y * y + z * z) / alpha;
	if (r < rmax) {
		lane_emden(r, dr, &m);
	} else {
		m = mtot;
	}
	fx = -m * x / (r * r * r);
	fy = -m * y / (r * r * r);
	fz = -m * z / (r * r * r);
}

std::vector<real> sod_shock_tube(real x, real y, real z) {
	std::vector < real > u(NF, real(0));
	if (x < real(0)) {
		u[rho_i] = sod_rho_l;
		u[egas_i] = sod_p_l / (fgamma - real(1));
	} else {
		u[rho_i] = sod_rho_r;
		u[egas_i] = sod_p_r / (fgamma - real(1));
	}
	u[tau_i] = std::pow(real(u[egas_i]), (real(1) / real(fgamma)));
	return u;
}

std::vector<real> blast_wave(real x, real y, real z) {
	std::vector < real > u(NF, real(0));
	u[rho_i] = real(1);
	if (std::sqrt(x * x + y * y + z * z) < 2.5e-2) {
		u[egas_i] = real(1);
	} else {
		u[egas_i] = real(2.5e-5);
	}
	u[tau_i] = std::pow(real(u[egas_i]), (real(1) / real(fgamma)));
	return u;
}
