/*
 * fourier_legendre.hpp
 *
 *  Created on: Apr 23, 2015
 *      Author: dmarce1
 */

#ifndef FOURIER_LEGENDRE_HPP_
#define FOURIER_LEGENDRE_HPP_

#include "RKDGFMM.hpp"

#include "legendre.hpp"

class fourier_legendre {
private:
	const integer N;
	const integer N3;
	integer N_Lobatto;
	std::vector<real> quadrature_point;
	std::vector<real> quadrature_weight;
	std::vector<real> lobatto_quadrature_point;
	std::vector<real> lobatto_quadrature_weight;
	std::vector<std::array<real,NDIM>> quadrature_point_3d;
	std::vector<simd_vector> transform_coefficient;
	std::vector<simd_vector> inverse_transform_coefficient;
	std::vector<std::vector<simd_vector>> lobatto_inverse_transform_coefficient;
	std::vector<std::vector<simd_vector>> volume_transform_coefficient;
	std::vector<std::vector<simd_vector>> volume_inverse_transform_coefficient;
	std::vector<std::vector<simd_vector>> surface_transform_coefficient;
	std::vector<std::vector<simd_vector>> surface_inverse_transform_coefficient;
public:
	integer pindex(integer l, integer m, integer n) const;
	fourier_legendre(integer);
	simd_vector transform(const simd_vector&) const;
	simd_vector inverse_transform(const simd_vector&) const;
	integer lobatto_point_count() const;
	real lobatto_edge_weight() const;
	simd_vector lobatto_inverse_transform(const simd_vector&, dimension) const;
	simd_vector volume_transform(dimension, const simd_vector&) const;
	simd_vector volume_inverse_transform(dimension, const simd_vector&) const;
	simd_vector surface_transform(face, const simd_vector&) const;
	simd_vector surface_inverse_transform(face, const simd_vector&) const;
	const std::vector<std::array<real,NDIM>>& quadrature_points() const;
	const std::vector<real>& quadrature_weights() const;
};



#endif /* FOURIER_LEGENDRE_HPP_ */
