/*
 * RKDGFMM.hpp
 *
 *  Created on: Apr 11, 2015
 *      Author: dmarce1
 */

#ifndef RKDGFMM_HPP_
#define RKDGFMM_HPP_

#include <hpx/hpx_init.hpp>
#include <hpx/include/components.hpp>
#include <hpx/lcos/wait_all.hpp>
#include <hpx/lcos/when_all.hpp>
#include <hpx/lcos/when_any.hpp>
#include <hpx/hpx_fwd.hpp>
#include <hpx/lcos/local/mutex.hpp>
#include <hpx/runtime/actions/component_action.hpp>
#include <hpx/runtime/components/server/managed_component_base.hpp>



#include <array>
#include <valarray>

//#define ANG_MOM

using real = double;
using integer = long long int;

#define NOGRAVITY

const real rho_floor = 0.0;

#define ALTENERGY

static constexpr integer P = 2;
static constexpr integer NRK = P;
static constexpr integer INX = 16;

static constexpr integer LMAX = P+2;
static constexpr integer L2 = LMAX*LMAX;
static constexpr integer NF = 6;
static constexpr integer BW = 2;
static constexpr integer NX = INX + 2 * BW;
static const real fgamma = real(5) / real(3);
static constexpr real dual_energy_switch1 = 1.0e-1;
static constexpr real dual_energy_switch2 = 1.0e-3;
static constexpr integer DNMAX = 3;

static constexpr integer N3 = NX * NX * NX;
static constexpr integer P3 = (P + 2) * (P + 1) * P / 6;
static constexpr integer G3 = P * P * P;
static constexpr integer G2 = P * P;

constexpr integer NDIM = 3;
constexpr integer NFACE = 2 * NDIM;
constexpr integer NVERTEX = 8;

static constexpr integer rho_i = 0;
static constexpr integer egas_i = 1;
static constexpr integer tau_i = 2;
static constexpr integer s_i = 3;

//static constexpr integer rho_i = 0;
static constexpr integer p_i = 1;
static constexpr integer ei_i = 2;
static constexpr integer ek_i = 3;
static constexpr integer h_i = 4;
static constexpr integer c_i = 5;
static constexpr integer v_i = 6;

#include "tests/initial.hpp"

enum dimension {
	XDIM = 0, YDIM = 1, ZDIM = 2
};

enum face {
	XM = 0, XP = 1, YM = 2, YP = 3, ZM = 4, ZP = 5
};

using xpoint = std::array<real,NDIM>;
using simd_vector = std::valarray<real>;

inline bool float_eq(real a, real b) {
	return std::abs(a - b) < 1.0e-10;
}

inline bool xpoint_eq(const xpoint& a, const xpoint& b) {
	bool rc = true;
	for (integer d = 0; d != NDIM; ++d) {
		if (!float_eq(a[d], b[d])) {
			rc = false;
			break;
		}
	}
	return rc;
}



constexpr integer chunk_size = 1000;

#endif /* RKDGFMM_HPP_ */


