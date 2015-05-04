/*
 * legendre.hpp
 *
 *  Created on: Apr 11, 2015
 *      Author: dmarce1
 */

#ifndef LEGENDRE_HPP_
#define LEGENDRE_HPP_

#include "RKDGFMM.hpp"

#include <array>
#include <cmath>
#include <vector>

std::vector<real> LegendreP(real x, integer N);
std::vector<real> AssociatedLegendreP(real x, integer N);
std::vector<real> SphericalHarmonic(real x, real y ,real z, integer N);
real LegendreP(integer n, real x);
real dLegendreP_dx(integer n, real x);
std::vector<real> dLegendreP_dx(real x, integer N);
std::vector<real> LegendreP_roots(integer nn);
std::vector<real> dLegendreP_dx_roots(integer nn);
real LegendreP_norm(integer);
integer factorial(integer);

#endif /* LEGENDRE_HPP_ */
