/*
 * fourier_legendre.hpp
 *
 *  Created on: Apr 23, 2015
 *      Author: dmarce1
 */

#ifndef FOURIER_LEGENDRE_HPP_
#define FOURIER_LEGENDRE_HPP_

#include "RKDGFMM.hpp"

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/valarray.hpp>
#include <boost/serialization/vector.hpp>
#include <fstream>

#include "legendre.hpp"

class fourier_legendre {
private:
	static integer PHI;
	static integer PHI3;
	static integer N_Lobatto;
	static std::vector<real> qpt;
	static std::vector<real> hires_qwt;
	static std::vector<real> hires_qpt;
	static std::vector<real> qwt;
	static std::vector<real> lobatto_qpt;
	static std::vector<real> lobatto_qwt;
	static std::vector<std::array<real, NDIM>> qpt_3d;
	static std::vector<simd_vector> transform_coefficient;
	static std::vector<simd_vector> vertex_inverse_transform_coefficient;
	static std::vector<simd_vector> inverse_transform_coefficient;
	static std::vector<std::vector<simd_vector>> dinverse_transform_coefficient_dx;
	static std::vector<std::vector<simd_vector>> lobatto_inverse_transform_coefficient;
	static std::vector<std::vector<simd_vector>> volume_transform_coefficient;
	static std::vector<std::vector<simd_vector>> volume_inverse_transform_coefficient;
	static std::vector<std::vector<simd_vector>> surface_transform_coefficient;
	static std::vector<std::vector<simd_vector>> surface_inverse_transform_coefficient;
	static std::vector<std::vector<simd_vector>> restrict_coefficients;
	static std::vector<std::vector<simd_vector>> prolong_coefficients;
	static std::vector<std::vector<simd_vector>> P2_phi;
	static std::vector<std::vector<simd_vector>> P2_gx;
	static std::vector<std::vector<simd_vector>> P2_gy;
	static std::vector<std::vector<simd_vector>> P2_gz;
	static std::vector<simd_vector> rho_2M;
	static std::vector<simd_vector> L2_phi;
	static std::vector<simd_vector> L2_gx;
	static std::vector<simd_vector> L2_gy;
	static std::vector<simd_vector> L2_gz;
	static std::vector<std::vector<simd_vector>> M2M;
	static std::vector<std::vector<simd_vector>> M2L;
	static std::vector<std::vector<simd_vector>> L2L;
	static simd_vector rpow_l;
	static simd_vector rpow_p;

	static void write(const char* filename);
	static void read(const char* filename);
	static void allocate();
public:
	static integer pindex(integer l, integer m, integer n);
	fourier_legendre();
	static simd_vector transform(const simd_vector&);
	static simd_vector inverse_transform(const simd_vector&);
	static simd_vector vertex_inverse_transform(const simd_vector&);
	static integer lobatto_point_count();
	static real lobatto_edge_weight();
	static simd_vector lobatto_inverse_transform(const simd_vector&, dimension);
	static simd_vector volume_transform(dimension, const simd_vector&);
	static simd_vector dinverse_transform_dx(dimension, const simd_vector&);
	static simd_vector volume_inverse_transform(dimension, const simd_vector&);
	static simd_vector surface_transform(face, const simd_vector&);
	static simd_vector surface_inverse_transform(face, const simd_vector&);
	static const std::vector<std::array<real, NDIM>>& quadrature_points();
	static const std::vector<real>& quadrature_weights();
	static simd_vector prolong(const simd_vector&, integer);
	static simd_vector _restrict(const simd_vector&, integer);
	static simd_vector p2_phi_transform(integer, integer, integer, const simd_vector&, real dx);
	static simd_vector p2_gx_transform(integer, integer, integer, const simd_vector&, real dx);
	static simd_vector p2_gy_transform(integer, integer, integer, const simd_vector&, real dx);
	static simd_vector p2_gz_transform(integer, integer, integer, const simd_vector&, real dx);
	static simd_vector m2l_transform(integer, integer, integer, const simd_vector&, real dx);
	static simd_vector m2m_transform(integer, const simd_vector&, real dx);
	static simd_vector l2l_transform(integer, const simd_vector&, real dx);
	static simd_vector l2p_transform(const simd_vector&, real dx);
	static simd_vector dl2p_transform_dx(dimension, const simd_vector&, real dx);
	static simd_vector p2m_transform(const simd_vector&, real dx);
};

#endif /* FOURIER_LEGENDRE_HPP_ */
