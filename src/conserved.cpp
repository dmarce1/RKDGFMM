/*
 * conserved.cpp
 *
 *  Created on: Apr 28, 2015
 *      Author: dmarce1
 */

#include "conserved.hpp"
#include "primitive.hpp"
#include <array>
#include <assert.h>

primitive_vars conserved_vars::to_primitive() const {
	const integer sz = rho().size();
	primitive_vars v(sz);
	v.rho() = rho();
	v.ek() = real(0);
	for (integer dim = 0; dim != NDIM; ++dim) {
		v.v(dimension(dim)) = s(dimension(dim)) / rho();
		v.ek() += v.v(dimension(dim)) * s(dimension(dim)) / real(2);
	}
	for (integer i = 0; i != sz; ++i) {
		if ((egas()[i] - v.ek()[i]) < dual_energy_switch2 * egas()[i]) {
			v.ei()[i] = std::pow(tau()[i], fgamma);
		} else {
			v.ei()[i] = std::max(egas()[i] - v.ek()[i], real(0));
		}
	}
	v.p() = (fgamma - real(1)) * v.ei();
	v.h() = (egas() + v.p()) / rho();
	v.c() = std::sqrt(fgamma * v.p() / v.rho());
	return v;
}

std::vector<simd_vector> HLLC_flux(const conserved_vars& UL, const conserved_vars& UR, const simd_vector& phi,
		dimension dim) {
	std::vector<simd_vector> F(NF, simd_vector (G2));
	integer shear1, shear2;
	if (dim == XDIM) {
		shear1 = s_i + YDIM;
		shear2 = s_i + ZDIM;
	} else if (dim == YDIM) {
		shear1 = s_i + XDIM;
		shear2 = s_i + ZDIM;
	} else {
		shear1 = s_i + XDIM;
		shear2 = s_i + YDIM;
	}
	auto qs = [](real p_star, real ps) {
		real gp1over2g = (fgamma + real(1)) / (real(2) * fgamma);
		const real H = p_star / ps;
		real r;
		if( H < real(1)) {
			r= real(1);
		} else {
			r= std::sqrt(real(1) + gp1over2g * (H- real(1)));
		}
		return r;
	};

	constexpr
	real hf = real(1) / real(2);
	real _1over8 = real(1) / real(8);
	const auto VR = UR.to_primitive();
	const auto VL = UL.to_primitive();
	const simd_vector& uR = VR.v(dim);
	const simd_vector& uL = VL.v(dim);
	const simd_vector& pR = VR.p();
	const simd_vector& pL = VL.p();
	const simd_vector& aR = VR.c();
	const simd_vector& aL = VL.c();
	const simd_vector& rhoR = VR.rho();
	const simd_vector& rhoL = VL.rho();
	simd_vector u_star = hf * (uL + uR) - (real(2) * (pR - pL)) / ((rhoR + rhoL) * (aR + aL));
	simd_vector p_star = hf * (pL + pR) - (_1over8 * (uR - uL) * (rhoR + rhoL) * (aR + aL));
	simd_vector& sM = u_star;
	simd_vector qR(G2), qL(G2);
	for (integer g = 0; g != G2; ++g) {
		qR[g] = qs(p_star[g], pR[g]);
		qL[g] = qs(p_star[g], pL[g]);
	}
	simd_vector sL = uL - aL * qL;
	simd_vector sR = uR + aR * qR;
	for (integer g = 0; g != G2; ++g) {
		sL[g] = std::min(real(0), sL[g]);
		sR[g] = std::max(real(0), sR[g]);
		sM[g] = std::min(sR[g], sM[g]);
		sM[g] = std::max(sL[g], sM[g]);
	}

	const auto fR = UR.flux(VR, dim);
	const auto fL = UL.flux(VL, dim);
	std::vector<simd_vector> Q(NF, simd_vector (G2));
	std::vector<simd_vector> R(NF, simd_vector (G2));
	for (integer f = 0; f != NF; ++f) {
		Q[f] = sL * UL[f] - fL[f];
		R[f] = sR * UR[f] - fR[f];
	}

	std::array < simd_vector, NF > U_starR;
	std::array < simd_vector, NF > U_starL;

	p_star = ((sM * Q[rho_i] - Q[s_i + dim]) + (sM * R[rho_i] - R[s_i + dim])) / real(2);
	p_star /= real(2);
	p_star = (p_star + std::abs(p_star));
	U_starL[rho_i] = Q[rho_i] / (sL - sM);
	U_starR[rho_i] = R[rho_i] / (sR - sM);
	U_starL[s_i + dim] = U_starL[rho_i] * u_star;
	U_starR[s_i + dim] = U_starR[rho_i] * u_star;
	U_starL[shear1] = UL[shear1];
	U_starR[shear1] = UR[shear1];
	U_starL[shear2] = UL[shear2];
	U_starR[shear2] = UR[shear2];
	U_starL[egas_i] = p_star / (fgamma - real(1))
			+ hf
					* (std::pow(U_starL[s_i + XDIM], real(2)) + std::pow(U_starL[s_i + YDIM], real(2))
							+ std::pow(U_starL[s_i + ZDIM], real(2))) / U_starL[rho_i];
	U_starR[egas_i] = p_star / (fgamma - real(1))
			+ hf
					* (std::pow(U_starR[s_i + XDIM], real(2)) + std::pow(U_starR[s_i + YDIM], real(2))
							+ std::pow(U_starR[s_i + ZDIM], real(2))) / U_starR[rho_i];
	U_starL[tau_i] = U_starR[tau_i] = std::pow(p_star / (fgamma - real(1)), real(1) / fgamma);

	for (integer f = 0; f != NF; ++f) {
		for (integer g = 0; g != G2; ++g) {
			assert(sR[g] >= sM[g]);
			assert(sM[g] >= sL[g]);
			if (sL[g] >= real(0) && sM[g] >= real(0)) {
				F[f][g] = fL[f][g];
			} else if (sL[g] < real(0) && sR[g] > real(0)) {
				if (sM[g] >= real(0)) {
					F[f][g] = fL[f][g] + sL[g] * (U_starL[f][g] - UL[f][g]);
				} else if (sM[g] <= real(0)) {
					F[f][g] = fR[f][g] + sR[g] * (U_starR[f][g] - UR[f][g]);
				}
			} else if (sM[g] <= real(0) && sR[g] <= real(0)) {
				F[f][g] = fR[f][g];
			} else {
				assert(false);
			}
		}
	}
	const integer f = egas_i;
	for (integer g = 0; g != G2; ++g) {
		assert(sR[g] >= sM[g]);
		assert(sM[g] >= sL[g]);
		if (sL[g] >= real(0) && sM[g] >= real(0)) {
			F[f][g] += UL[s_i + dim][g] * phi[g];
		} else if (sL[g] < real(0) && sR[g] > real(0)) {
			if (sM[g] >= real(0)) {
				F[f][g] += phi[g] * (UL[s_i + dim][g] + sL[g] * (U_starL[rho_i][g] - UL[rho_i][g]));
			} else if (sM[g] <= real(0)) {
				F[f][g] += phi[g] * (UR[s_i + dim][g] + sR[g] * (U_starR[rho_i][g] - UR[rho_i][g]));
			}
		} else if (sM[g] <= real(0) && sR[g] <= real(0)) {
			F[f][g] += UR[s_i + dim][g] * phi[g];
		} else {
			assert(false);
		}
	}
	return F;
}

std::vector<simd_vector> conserved_vars::flux(const primitive_vars& v, dimension dim) const {
	std::vector < simd_vector > f(nfields);
	f[rho_i] = rho() * v.v(dim);
	f[egas_i] = rho() * v.h() * v.v(dim);
	for (integer d = 0; d != NDIM; ++d) {
		const dimension this_dim = dimension(d);
		f[s_i + this_dim] = s(this_dim) * v.v(dim);
	}
	f[s_i + dim] += v.p();
	f[tau_i] = tau() * v.v(dim);
	return f;
}

