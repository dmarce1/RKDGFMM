/*
 * initial.hpp
 *
 *  Created on: Apr 24, 2015
 *      Author: dmarce1
 */

#ifndef INITIAL_HPP_
#define INITIAL_HPP_

#include "../RKDGFMM.hpp"

constexpr real sod_rho_l = real(1);
constexpr real sod_rho_r = real(1) / real(8);
constexpr real sod_p_l = real(1);
constexpr real sod_p_r = real(1);

std::vector<real> sod_shock_tube(real x, real y, real z);
std::vector<real> blast_wave(real x, real y, real z);
std::vector<real> star(real x, real y, real z);
void star_force(real x, real y, real z, real& fx, real& fy, real& fz);

#endif /* INITIAL_HPP_ */
