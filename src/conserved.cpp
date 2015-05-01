/*
 * conserved.cpp
 *
 *  Created on: Apr 28, 2015
 *      Author: dmarce1
 */

#include "conserved.hpp"
#include "primitive.hpp"
#include <array>

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

std::vector<simd_vector> conserved_vars::flux(const primitive_vars& v, dimension dim) const {
	std::vector<simd_vector> f(nfields);
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

