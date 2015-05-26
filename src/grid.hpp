/*
 * grid.hpp
 *
 *  Created on: Apr 23, 2015
 *      Author: dmarce1
 */

#include "RKDGFMM.hpp"

#include "fourier_legendre.hpp"
#include "conserved.hpp"
#include "primitive.hpp"

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/valarray.hpp>
#include <boost/serialization/vector.hpp>

class grid {
private:
	std::vector<std::vector<conserved_vars>> U_p;
	real t, grid_amax;
	real output_dt;
	integer output_cnt;
	integer dx_i, dy_i, dz_i;
	std::array<integer, NDIM> d_i;
	const real h;
	const real hinv;
	fourier_legendre Fourier;
	integer nlevel;
	std::vector<std::vector<simd_vector>> rho_l;
	std::vector<simd_vector> phi_p_analytic;
	std::vector<simd_vector> gx_h_analytic;
	std::vector<simd_vector> gy_h_analytic;
	std::vector<simd_vector> gz_h_analytic;

	std::vector<std::vector<simd_vector>> phi_h;
	std::array<std::vector<std::vector<simd_vector>>, NDIM> g_h;
	std::vector<std::vector<simd_vector>>& gx_h;
	std::vector<std::vector<simd_vector>>& gy_h;
	std::vector<std::vector<simd_vector>>& gz_h;

	std::vector<std::vector<simd_vector>> phi_l;
	std::vector<conserved_vars> U_h_analytic;
	std::vector<std::vector<std::vector<simd_vector>>>dU_dt_p;
	std::vector<std::vector<real>> dO_dt;
	std::vector<std::vector<real>> O;
	std::vector<std::vector<std::vector<simd_vector>>> F_p;
	std::vector<std::vector<simd_vector>> S_p;
	std::vector<bool> is_interior;
	std::vector<bool> is_on_edge;
	std::vector<real> cell_x;
	std::vector<real> cell_y;
	std::vector<real> cell_z;

	std::vector<simd_vector> con_to_prim(const std::vector<simd_vector>&) const;
	std::vector<simd_vector> prim_to_con(const std::vector<simd_vector>&) const;

	void apply_limiter(conserved_vars& U0_p, const conserved_vars& UR_p, const conserved_vars& UL_p, const simd_vector fg_p, dimension dim) const;
	std::vector<simd_vector> prim_to_roe_average(const std::vector<simd_vector>&) const;

public:
	grid();
	real step();
	void initialize(const char*);
	void initialize(std::function<std::vector<real>(real, real, real)>&&);
	void compute_analytic(std::function<std::vector<real>(real, real, real, real)>&&, real t);
	void enforce_positivity(integer rk);
	void compute_flux(integer rk);
	real compute_cfl_condition();
	void compute_du(integer rk);
	void compute_next_u(integer rk, real dt);
	integer project(std::vector<conserved_vars>&, const std::vector<simd_vector>&, const std::vector<simd_vector>&,const std::vector<simd_vector>&, const std::vector<simd_vector>& );
	void enforce_positivity(std::vector<conserved_vars>&, const std::vector<simd_vector>& );
	integer project(integer rk);
	void output( const char*) const;
	void enforce_boundaries(std::vector<conserved_vars>&);
	void enforce_boundaries(integer rk);
	void diagnostics(real t);
	void compute_multipoles(integer rk);
	void compute_interactions(integer rk);
	void expand_phi(integer rk);
	void compute_force(integer rk);
	void filter_fmm(integer rk);
	void compute_fmm(integer rk) {
		compute_multipoles(rk);
		compute_interactions(rk);
		expand_phi(rk);
		compute_force(rk);
	//	filter_fmm(rk);
	}
	template<class Arc>
	void serialize( Arc& arc, const unsigned)  {
		arc & t;
		arc & grid_amax;
		arc & U_p[0];
		arc & phi_h;
		arc & gx_h;
		arc & gy_h;
		arc & gz_h;
		arc & O[0];
	}
	void save( const char* filename) const {
		std::ofstream ofs(filename);
		boost::archive::binary_oarchive arc(ofs);
		arc << *this;
		ofs.close();
	}
	void load( const char* filename) {
		std::ifstream ifs(filename);
		boost::archive::binary_iarchive arc(ifs);
		arc >> *this;
		ifs.close();
	}
};
