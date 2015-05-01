/*
 * quadrature.hpp
 *
 *  Created on: Apr 25, 2015
 *      Author: dmarce1
 */

#ifndef QUADRATURE_HPP_
#define QUADRATURE_HPP_

#include "RKDGFMM.hpp"

constexpr integer PMAX = 4;

real quad_points[PMAX][PMAX + 1] = { { -real(1), +real(1) }, { -real(1), real(0), +real(1) }, { -real(1), -real(1)
		/ std::sqrt(real(5)), +real(1) / std::sqrt(real(5)), +real(1) }, { -real(1), -std::sqrt(real(3) / real(7)),
		real(0), +std::sqrt(real(3) / real(7)), +real(1) } };

real quad_weights[PMAX][PMAX + 1] = { { real(1), real(1) }, { real(1) / real(3), real(4) / real(3), real(1) / real(3) },
		{ real(1) / real(6), real(5) / real(6), real(5) / real(6), real(1) / real(6) }, { real(1) / real(10), real(49)
				/ real(90), real(32) / real(45), real(49) / real(90), real(1) / real(10) } };

#endif /* QUADRATURE_HPP_ */
