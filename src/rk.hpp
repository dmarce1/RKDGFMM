/*
 * rk.hpp
 *
 *  Created on: Mar 5, 2015
 *      Author: dmarce1
 */

#ifndef RK_HPP_
#define RK_HPP_

#include "RKDGFMM.hpp"

constexpr integer NRK_MAX = 4;

constexpr real alpha_rk[NRK_MAX][NRK_MAX][NRK_MAX] = { { { real(1) } }, { { real(1) }, { real(1) / real(2), real(1)
		/ real(2) } }, { { real(1) }, { real(3) / real(4), real(1) / real(4) }, { real(1) / real(3), real(0), real(2)
		/ real(3) } }, { { real(1) }, { real(1) / real(2), real(1) / real(2) }, { real(1) / real(9), real(2) / real(9),
		real(2) / real(3) }, { real(0), real(1) / real(3), real(1) / real(3), real(1) / real(3) } } };

constexpr real beta_rk[NRK_MAX][NRK_MAX][NRK_MAX] = { { { real(1) } }, { { real(1) }, { real(0), real(1) / real(2) } },
		{ { real(1) }, { real(0), real(1) / real(4) }, { real(0), real(0), real(2) / real(3) } }, {
				{ real(1) / real(2) }, { -real(1) / real(4), real(1) / real(2) }, { -real(1) / real(9), -real(1)
						/ real(3), real(1) }, { real(0), real(1) / real(6), real(0), real(1) / real(6) } } };

constexpr real cfl[NRK_MAX] = { real(1), real(1), real(1), real(2) / real(3) };

#endif /* RK_HPP_ */
