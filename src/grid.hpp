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

class grid {
private:

	integer dx_i, dy_i, dz_i;
	std::array<integer, NDIM> d_i;
	const real dx;
	fourier_legendre Fourier;
	integer nlevel;
	std::vector<std::vector<simd_vector>> rho_l;
	std::vector<simd_vector> phi_p;
	std::vector<std::vector<simd_vector>> phi_l;
	std::vector<std::vector<conserved_vars>> U_p;
	std::vector<std::vector<std::vector<simd_vector>>>dU_dt_p;
	std::vector<std::vector<std::vector<simd_vector>>> F_p;
	std::vector<std::vector<simd_vector>> S_p;
	std::vector<bool> is_interior;
	std::vector<bool> is_on_edge;
	std::vector<real> cell_x;
	std::vector<real> cell_y;
	std::vector<real> cell_z;

	std::vector<simd_vector> con_to_prim(const std::vector<simd_vector>&) const;
	std::vector<simd_vector> prim_to_con(const std::vector<simd_vector>&) const;

	void apply_limiter(conserved_vars& U0_p, const conserved_vars& UR_p, const conserved_vars& UL_p, dimension dim) const;
	std::vector<simd_vector> prim_to_roe_average(const std::vector<simd_vector>&) const;


public:
	grid();
	void initialize(std::function<std::vector<real>(real, real, real)>&&);
	real enforce_positivity(integer rk);
	void compute_flux(integer rk);
	void compute_du(integer rk);
	void compute_next_u(integer rk, real dt);
	void project(integer rk);
	void output( const char*) const;
	void enforce_boundaries(integer rk);
	void diagnostics(integer rk );
	void compute_multipoles(integer rk);
	void compute_interactions(integer rk);
	void expand_phi(integer rk);
	void compute_fmm(integer rk) {
		compute_multipoles(rk);
		compute_interactions(rk);
		expand_phi(rk);
	}
};
