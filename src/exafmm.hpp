/*
 * exafmm.hpp
 *
 *  Created on: May 3, 2015
 *      Author: dmarce1
 */

#ifndef EXAFMM_HPP_
#define EXAFMM_HPP_

#include "RKDGFMM.hpp"
#include <complex>
#include <valarray>
#include <array>

using complex = std::complex<real>;

struct exafmm {
	std::valarray<complex> M2M(const std::valarray<complex>& Mj, std::array<real, NDIM> dist);
	std::valarray<complex> M2L(const std::valarray<complex>& Mj, std::array<real, NDIM> dist);
	std::valarray<complex> L2L(const std::valarray<complex>& Lj, std::array<real, NDIM> dist);
	std::valarray<complex> P2M(std::array<real, NDIM> dist);
	std::valarray<real> L2P(const std::valarray<complex>& Li, std::array<real, NDIM> dist);
	exafmm();
};

#endif /* EXAFMM_HPP_ */
