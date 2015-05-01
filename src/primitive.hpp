/*
 * primitive.hpp
 *
 *  Created on: Apr 28, 2015
 *      Author: dmarce1
 */

#ifndef PRIMITIVE_HPP_
#define PRIMITIVE_HPP_

#include "RKDGFMM.hpp"

class conserved_vars;

class primitive_vars {
private:
	static constexpr integer nfields = 9;
	std::vector<simd_vector> vars;
public:
	primitive_vars(integer sz) :
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
	const simd_vector& v(dimension dim) const {
		return vars[v_i + dim];
	}
	simd_vector& v(dimension dim) {
		return vars[v_i + dim];
	}
	const simd_vector& p() const {
		return vars[p_i];
	}
	simd_vector& p() {
		return vars[p_i];
	}
	const simd_vector& h() const {
		return vars[h_i];
	}
	simd_vector& h() {
		return vars[h_i];
	}
	const simd_vector& c() const {
		return vars[c_i];
	}
	simd_vector& c() {
		return vars[c_i];
	}
	const simd_vector& ei() const {
		return vars[ei_i];
	}
	simd_vector& ei() {
		return vars[ei_i];
	}
	const simd_vector& ek() const {
		return vars[ek_i];
	}
	simd_vector& ek() {
		return vars[ek_i];
	}
	conserved_vars to_conserved() const;
	std::vector<real> roe_avg() const;
};

#endif /* PRIMITIVE_HPP_ */
