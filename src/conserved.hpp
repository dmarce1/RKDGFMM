/*
 * conserved.hpp
 *
 *  Created on: Apr 28, 2015
 *      Author: dmarce1
 */

#ifndef CONSERVED_HPP_
#define CONSERVED_HPP_

#include "RKDGFMM.hpp"

class primitive_vars;

class conserved_vars {
private:
	static constexpr integer nfields = NF;
	std::vector<simd_vector> vars;
public:
	conserved_vars(integer sz) :
			vars(nfields, simd_vector(sz)) {
	}
	const simd_vector& operator[](integer i) const {
		return vars[i];
	}
	simd_vector& operator[](integer i) {
		return vars[i];
	}
	const simd_vector& rho() const {
		return vars[rho_i];
	}
	simd_vector& rho() {
		return vars[rho_i];
	}
	const simd_vector& s(dimension dim) const {
		return vars[s_i + dim];
	}
	simd_vector& s(dimension dim) {
		return vars[s_i + dim];
	}
	const simd_vector& egas() const {
		return vars[egas_i];
	}
	simd_vector& egas() {
		return vars[egas_i];
	}
	const simd_vector& tau() const {
		return vars[tau_i];
	}
	simd_vector& tau() {
		return vars[tau_i];
	}
	std::vector<simd_vector> flux(const primitive_vars&, dimension dim) const;
	primitive_vars to_primitive() const;
	real max_signal_speed() const;
};

std::vector<simd_vector> HLLC_flux(const conserved_vars& UL, const conserved_vars& UR, const simd_vector&, dimension dim);

#endif /* CONSERVED_HPP_ */
