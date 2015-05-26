/*O.
 * grid.cpp
 *
 *  Created on: Apr 23, 2015
 *      Author: dmarce1
 */

#include "grid.hpp"
#include "rk.hpp"
#include <cassert>
#include <list>
#include <limits>
#include <silo.h>
#include <set>
#include <utility>
#include <mpi.h>


real minmod(real a, real b) {
	real c;
	if( a*b <= real(0)) {
		c =real(0);
	} else if( std::abs(a) < std::abs(b)) {
		c= a;
	} else {
		c= b;
	}
	//real c1 = (std::copysign(real(1), a) + std::copysign(real(1), b)) * std::min(std::abs(a), std::abs(b)) / real(2);
	//if( c1 != c ) {
	//	printf( "%e %e %e %e\n", a, b, c, c1);
//	}
	return c;
}

simd_vector minmod(const simd_vector& a, const simd_vector& b) {
	const integer sz = a.size();
	simd_vector c(sz);
	for (integer i = 0; i != sz; ++i) {
		c[i] = (std::copysign(real(1), a[i]) + std::copysign(real(1), b[i])) * std::min(std::abs(a[i]), std::abs(b[i]))
				/ real(2);
	}
	return c;
}

real grid::step() {
#ifdef NDEBUG
		real tstart = MPI_Wtime();
#endif
		diagnostics(t);
		real dt = ((real(2)/real(NX)) / grid_amax) * cfl[NRK - 1];
		printf("t = %e dt = %e ", double(t + dt), double(dt));
		for (integer rk = 0; rk != NRK; ++rk) {
			const integer rk2 = rk == NRK - 1 ? 0 : rk + 1;
			compute_flux(rk);
			compute_du(rk);
			compute_next_u(rk, dt);
			enforce_boundaries(rk2);
			phi_h[rk2] = phi_h[rk];
			gx_h[rk2] = gx_h[rk];
			gy_h[rk2] = gy_h[rk];
			gz_h[rk2] = gz_h[rk];
			enforce_positivity(rk2);
			project(rk2);
//			enforce_positivity(rk2);
			compute_fmm(rk2);
			if( rk2 == 0) {
				grid_amax = compute_cfl_condition();
			}
		}
		t += dt;
		if (output_cnt <= integer(t / output_dt)) {
			char* tmp_ptr;
			if (!asprintf(&tmp_ptr, "X.%i.silo", int(output_cnt))) {
				abort();
			}
			printf("*");
			output(tmp_ptr);
			++output_cnt;
			free(tmp_ptr);
		}
#ifdef NDEBUG
		printf( " wt = %e ", MPI_Wtime() - tstart);
#endif
		printf( "\n");
	return t;
}

grid::grid() :
		h(real(1) / real(INX)), hinv(real(INX)), gx_h(g_h[XDIM]), gy_h(g_h[YDIM]), gz_h(g_h[ZDIM]) {

	O.resize(NRK+1,std::vector<real>(NF,real(0)));
	dO_dt.resize(NRK+1,std::vector<real>(NF,real(0)));
	U_p.resize(NRK + 1, std::vector < conserved_vars > (N3, conserved_vars(P3)));
	U_h_analytic.resize(N3, conserved_vars(G3));
	F_p.resize(NFACE,
			std::vector < std::vector < simd_vector >> (N3, std::vector < simd_vector > (NF, simd_vector(P3))));
	S_p.resize(N3, std::vector < simd_vector > (NF, simd_vector(P3)));
	dU_dt_p.resize(NRK,
			std::vector < std::vector < simd_vector >> (N3, std::vector < simd_vector > (NF, simd_vector(P3))));
	nlevel = 0;
	for (integer inx = INX; inx > 1; inx /= 2) {
		++nlevel;
	}
	phi_h.resize(NRK+1,std::vector<simd_vector>(NX * NX * NX, simd_vector(real(0), G3)));
	gx_h.resize(NRK+1,std::vector<simd_vector>(NX * NX * NX, simd_vector(real(0), G3)));
	gy_h.resize(NRK+1,std::vector<simd_vector>(NX * NX * NX, simd_vector(real(0), G3)));
	gz_h.resize(NRK+1,std::vector<simd_vector>(NX * NX * NX, simd_vector(real(0), G3)));
	phi_l.resize(nlevel);
	phi_p_analytic.resize(NX * NX * NX, simd_vector(real(0), P3));
	gx_h_analytic.resize(NX * NX * NX, simd_vector(real(0), G3));
	gy_h_analytic.resize(NX * NX * NX, simd_vector(real(0), G3));
	gz_h_analytic.resize(NX * NX * NX, simd_vector(real(0), G3));
	rho_l.resize(nlevel);
	nlevel = 0;
	for (integer inx = INX; inx > 1; inx /= 2) {
		const integer nx = inx + 2 * BW;
		rho_l[nlevel].resize(nx * nx * nx, simd_vector(L2));
		phi_l[nlevel].resize(nx * nx * nx, simd_vector(L2));
		++nlevel;
	}

	d_i[XDIM] = dx_i = NX * NX;
	d_i[YDIM] = dy_i = NX;
	d_i[ZDIM] = dz_i = 1;
	is_interior.resize(N3, false);
	is_on_edge.resize(N3, false);
	cell_x.resize(N3);
	cell_y.resize(N3);
	cell_z.resize(N3);
	for (integer i = 0; i != NX; ++i) {
		for (integer j = 0; j != NX; ++j) {
			is_on_edge[i * NX * NX + j * NX + (0)] = true;
			is_on_edge[i * NX * NX + j * NX + (NX - 1)] = true;
			is_on_edge[i * NX * NX + j + (0) * NX] = true;
			is_on_edge[i * NX * NX + j + (NX - 1) * NX] = true;
			is_on_edge[i * NX + j + (0) * NX * NX] = true;
			is_on_edge[i * NX + j + (NX - 1) * NX * NX] = true;
		}
	}
	for (integer i = BW; i != NX - BW; ++i) {
		for (integer j = BW; j != NX - BW; ++j) {
			for (integer k = BW; k != NX - BW; ++k) {
				const integer index = i * dx_i + j * dy_i + k * dz_i;
				is_interior[index] = true;
			}
		}
	}
	for (integer i = 0; i != NX; ++i) {
		for (integer j = 0; j != NX; ++j) {
			for (integer k = 0; k != NX; ++k) {
				const integer index = i * dx_i + j * dy_i + k * dz_i;
				cell_x[index] = real(2 * i - 2 * BW - INX + 1) * h;
				cell_y[index] = real(2 * j - 2 * BW - INX + 1) * h;
				cell_z[index] = real(2 * k - 2 * BW - INX + 1) * h;
			}
		}
	}
}


void grid::enforce_positivity(integer rk) {
	return enforce_positivity(U_p[rk], phi_h[rk]);
}

void grid::enforce_positivity(std::vector<conserved_vars>& this_u_p, const std::vector<simd_vector>& this_phi_h ) {
	const integer G3L = Fourier.lobatto_point_count();
	std::list<hpx::future<std::pair<real,real>>> futs;

	for( integer i_chunk = 0; i_chunk < N3; i_chunk += chunk_size ) {
		futs.push_back(hpx::async([&](integer this_i_chunk) {
			real eps_rho = std::numeric_limits<real>::max();
			conserved_vars u_h(G3);
			for (integer iii = this_i_chunk; iii != N3 && iii != this_i_chunk + chunk_size; ++iii) {
				auto& u_p = this_u_p[iii];
				for (integer f = 0; f != NF; ++f) {
					u_h[f] = Fourier.inverse_transform(u_p[f]);
				}
				for (integer iii = this_i_chunk; iii != N3 && iii != this_i_chunk + chunk_size; ++iii) {
					const simd_vector rho = u_h.rho();
					eps_rho = std::min(eps_rho, u_p.rho()[0]);
				}
			}
			return std::make_pair(eps_rho,eps_rho);
		}, i_chunk));
	}

	real eps_rho = std::numeric_limits<real>::max();
	for( auto i = futs.begin(); i != futs.end(); ++i) {
		auto mns = i->get();
		eps_rho = std::min(eps_rho, mns.first);
	}
	eps_rho *= real(1.00);



	std::list<hpx::future<void>> futs2;

	for( integer i_chunk = 0; i_chunk < N3; i_chunk += chunk_size ) {
		futs2.push_back(hpx::async([&](integer this_i_chunk) {
			std::vector<conserved_vars> u_h(NDIM, conserved_vars(G3L));
			for (integer iii = this_i_chunk; iii != N3 && iii != this_i_chunk + chunk_size; ++iii) {
				auto& u_p = this_u_p[iii];
					for (integer dim = 0; dim != NDIM; ++dim) {
					for (integer f = 0; f != NF; ++f) {
						u_h[dim][f] = Fourier.lobatto_inverse_transform(u_p[f], dimension(dim));
					}
				}
				real rho_min = std::numeric_limits<real>::max();
				for (integer dim = 0; dim != NDIM; ++dim) {
					rho_min = std::min(rho_min, u_h[dim].rho().min());
				}
				const real this_rho = u_p.rho()[0];
				const real theta = this_rho != rho_min ? (this_rho - eps_rho) / (this_rho - rho_min) : real(1);
				if (theta < real(1)) {
					const simd_vector rho_before_p = u_p.rho();
					for (integer f = 0; f != NF; ++f) {
						for (integer ppp = 1; ppp < P3; ++ppp) {
							u_p[f][ppp] *= theta;
						}
					}
					const simd_vector rho_after_p = u_p.rho();
					const simd_vector drho_h = Fourier.inverse_transform(rho_after_p - rho_before_p);
					u_p.egas() -= Fourier.transform(drho_h * this_phi_h[iii]);

				}
			}
		}, i_chunk));
	}
	hpx::wait_all(futs2.begin(), futs2.end());


	futs.clear();
	futs2.clear();
	for( integer i_chunk = 0; i_chunk < N3; i_chunk += chunk_size ) {
		futs.push_back(hpx::async([&](integer this_i_chunk) {
			conserved_vars u_h(G3);
			real eps_tau = std::numeric_limits<real>::max();
			for (integer iii = this_i_chunk; iii != N3 && iii != this_i_chunk + chunk_size; ++iii) {
				auto& u_p = this_u_p[iii];
				for (integer f = 0; f != NF; ++f) {
					u_h[f] = Fourier.inverse_transform(u_p[f]);
				}
				const simd_vector egas = u_h.egas();
				const simd_vector rho = u_h.rho();
				const simd_vector u = simd_vector(u_h.s(XDIM) / rho);
				const simd_vector v = simd_vector(u_h.s(YDIM) / rho);
				const simd_vector w = simd_vector(u_h.s(ZDIM) / rho);
				const simd_vector ei = simd_vector(egas - rho * (u * u + v * v + w * w) / real(2));
				for (integer ggg = 0; ggg != G3; ++ggg) {
					if ((ei[ggg] > dual_energy_switch1) * egas[ggg] && (ei[ggg] > real(0))) {
						u_h.tau()[ggg] = std::pow(ei[ggg], real(1) / fgamma);
					}
				}
				u_p.tau() = Fourier.transform(u_h.tau());
				eps_tau = std::min(eps_tau, u_p.tau()[0]);
			}
			return std::make_pair(eps_tau,eps_tau);
		}, i_chunk));
	}

	real eps_tau = std::numeric_limits<real>::max();
	for( auto i = futs.begin(); i != futs.end(); ++i) {
		auto mns = i->get();
		eps_tau = std::min(eps_tau, mns.second);
	}
	eps_tau *= real(1.00);



	for( integer i_chunk = 0; i_chunk < N3; i_chunk += chunk_size ) {
		futs2.push_back(hpx::async([&](integer this_i_chunk) {
			std::vector<conserved_vars> u_h(NDIM, conserved_vars(G3L));
			for (integer iii = this_i_chunk; iii != N3 && iii != this_i_chunk + chunk_size; ++iii) {
				auto& u_p = this_u_p[iii];
				for (integer dim = 0; dim != NDIM; ++dim) {
					for (integer f = 0; f != NF; ++f) {
						u_h[dim][f] = Fourier.lobatto_inverse_transform(u_p[f], dimension(dim));
					}
				}
				real tau_min = std::numeric_limits<real>::max();
				for (integer dim = 0; dim != NDIM; ++dim) {
					tau_min = std::min(tau_min, u_h[dim].tau().min());
				}
				const real this_tau = u_p.tau()[0];
				const real theta = this_tau != tau_min ? (this_tau - eps_tau) / (this_tau - tau_min) : real(1);
				if (theta < real(1)) {
					const simd_vector rho_before_p = u_p.rho();
					for (integer f = 0; f != NF; ++f) {
						for (integer ppp = 1; ppp < P3; ++ppp) {
							u_p[f][ppp] *= theta;
						}
						for (integer dim = 0; dim != NDIM; ++dim) {
							u_h[dim][f] = Fourier.lobatto_inverse_transform(u_p[f], dimension(dim));
						}
					}
					const simd_vector rho_after_p = u_p.rho();
					const simd_vector drho_h = Fourier.inverse_transform(rho_after_p - rho_before_p);
					u_p.egas() -= Fourier.transform(drho_h * this_phi_h[iii]);

				}
				assert(u_h[XDIM].tau().min() > real(0));
				assert(u_h[YDIM].tau().min() > real(0));
				assert(u_h[ZDIM].tau().min() > real(0));
			}
		}, i_chunk));
	}
	hpx::wait_all(futs2.begin(), futs2.end());
}


real grid::compute_cfl_condition() {
	std::vector<conserved_vars>& this_u_p = U_p[0];
	const real cfl = Fourier.lobatto_edge_weight() / real(4);
	const integer G3L = Fourier.lobatto_point_count();
	std::list<hpx::future<real>> futs2;

	real amax = real(0);
//	for( integer i_chunk = 0; i_chunk < N3; i_chunk += chunk_size ) {
	//	futs2.push_back(hpx::async([&](integer this_i_chunk) {
			std::vector<conserved_vars> u_h(NDIM, conserved_vars(G3L));
			for (integer iii = 0; iii != N3; ++iii) {
				auto& u_p = this_u_p[iii];
				for (integer dim = 0; dim != NDIM; ++dim) {
					for (integer f = 0; f != NF; ++f) {
						u_h[dim][f] = Fourier.lobatto_inverse_transform(u_p[f], dimension(dim));
					}
				}
				for (integer dim = 0; dim != NDIM; ++dim) {
					const auto prims = u_h[dim].to_primitive();
					const real a_x = (std::abs(prims.v(XDIM)) + prims.c()).max();
					const real a_y = (std::abs(prims.v(YDIM)) + prims.c()).max();
					const real a_z = (std::abs(prims.v(ZDIM)) + prims.c()).max();
					amax = std::max(amax, a_x + a_y + a_z);
				}
			}
		//	return amax;
	//	}, i_chunk));
//	}
//	real amax = real(0);
//	for( auto i = futs2.begin(); i != futs2.end(); ++i) {
//		amax = std::max(i->get(), amax);
//	}

	return amax / cfl;
}


void grid::enforce_boundaries(integer rk) {
	enforce_boundaries(U_p[rk]);
}

void grid::enforce_boundaries(std::vector<conserved_vars>& u) {

	for (integer j = BW; j != NX - BW; ++j) {
		for (integer k = BW; k != NX - BW; ++k) {
			const integer iiia = BW * dx_i + j * dy_i + k * dz_i;
			const integer im1a = iiia - dx_i;
			const integer im2a = iiia - 2 * dx_i;
			const integer iiib = (NX - BW - 1) * dx_i + j * dy_i + k * dz_i;
			const integer ip1b = iiib + dx_i;
			const integer ip2b = iiib + 2 * dx_i;
			for (integer f = 0; f != NF; ++f) {
				for (integer l = 0; l != P; ++l) {
					for (integer m = 0; m != P - l; ++m) {
						for (integer n = 0; n != P - l - m; ++n) {
							const integer pppx = Fourier.pindex(l, m, n);
								u[im1a][f][pppx] = u[iiia][f][pppx];
								u[im2a][f][pppx] = u[iiia][f][pppx];
								u[ip1b][f][pppx] = u[iiib][f][pppx];
								u[ip2b][f][pppx] = u[iiib][f][pppx];
						}
					}
				}
			}
			if (u[iiia][s_i + XDIM][0] > real(0)) {
				u[im1a][s_i + XDIM][0] = real(0);
				u[im2a][s_i + XDIM][0] = real(0);
			}
			if (u[iiib][s_i + XDIM][0] < real(0)) {
				u[ip1b][s_i + XDIM][0] = real(0);
				u[ip2b][s_i + XDIM][0] = real(0);
			}
		}
	}

	for (integer j = 0; j != NX; ++j) {
		for (integer k = BW; k != NX - BW; ++k) {
			const integer jjja = j * dx_i + BW * dy_i + k * dz_i;
			const integer jm1a = jjja - dy_i;
			const integer jm2a = jjja - 2 * dy_i;
			const integer jjjb = j * dx_i + (NX - BW - 1) * dy_i + k * dz_i;
			const integer jp1b = jjjb + dy_i;
			const integer jp2b = jjjb + 2 * dy_i;
			for (integer f = 0; f != NF; ++f) {
				for (integer l = 0; l != P; ++l) {
					for (integer m = 0; m != P - l; ++m) {
						for (integer n = 0; n != P - l - m; ++n) {
							const integer pppy = Fourier.pindex(m, l, n);
								u[jm1a][f][pppy] = u[jjja][f][pppy];
								u[jm2a][f][pppy] = u[jjja][f][pppy];
								u[jp1b][f][pppy] = u[jjjb][f][pppy];
								u[jp2b][f][pppy] = u[jjjb][f][pppy];
						}
					}
				}
			}
			if (u[jjja][s_i + YDIM][0] > real(0)) {
				u[jm1a][s_i + YDIM][0] = real(0);
				u[jm2a][s_i + YDIM][0] = real(0);
			}
			if (u[jjjb][s_i + YDIM][0] < real(0)) {
				u[jp1b][s_i + YDIM][0] = real(0);
				u[jp2b][s_i + YDIM][0] = real(0);
			}
		}
	}
	for (integer j = 0; j != NX; ++j) {
		for (integer k = 0; k != NX; ++k) {
			const integer kkka = k * dx_i + j * dy_i + BW * dz_i;
			const integer km1a = kkka - dz_i;
			const integer km2a = kkka - 2 * dz_i;
			const integer kkkb = k * dx_i + j * dy_i + (NX - BW - 1) * dz_i;
			const integer kp1b = kkkb + dz_i;
			const integer kp2b = kkkb + 2 * dz_i;
			for (integer f = 0; f != NF; ++f) {
				for (integer l = 0; l != P; ++l) {
					for (integer m = 0; m != P - l; ++m) {
						for (integer n = 0; n != P - l - m; ++n) {
							const integer pppz = Fourier.pindex(m, n, l);
								u[km1a][f][pppz] = u[kkka][f][pppz];
								u[km2a][f][pppz] = u[kkka][f][pppz];
								u[kp1b][f][pppz] = u[kkkb][f][pppz];
								u[kp2b][f][pppz] = u[kkkb][f][pppz];
						}
					}
				}
			}
			if (u[kkka][s_i + ZDIM][0] > real(0)) {
				u[km1a][s_i + ZDIM][0] = real(0);
				u[km2a][s_i + ZDIM][0] = real(0);
			}
			if (u[kkkb][s_i + ZDIM][0] < real(0)) {
				u[kp1b][s_i + ZDIM][0] = real(0);
				u[kp2b][s_i + ZDIM][0] = real(0);
			}
		}
	}
}

integer grid::project(integer rk ) {
	integer cnt;
	project(U_p[rk], phi_h[rk], gx_h[rk], gy_h[rk], gz_h[rk]);
	return cnt;
}

integer grid::project(std::vector<conserved_vars>& u_p, const std::vector<simd_vector>& this_phi_h,
		const std::vector<simd_vector>& this_gx_h,
		const std::vector<simd_vector>& this_gy_h,
		const std::vector<simd_vector>& this_gz_h ) {
	integer cnt = 0;


		std::list<hpx::future<void>> futs;
	std::vector<conserved_vars> u0_p(N3, conserved_vars(P3));

	for( integer i_chunk = 0; i_chunk < N3; i_chunk += chunk_size ) {
		futs.push_back(hpx::async([&](integer this_i_chunk) {
			for (integer iii = this_i_chunk; iii != N3 && iii != this_i_chunk + chunk_size; ++iii) {
				u0_p[iii] = u_p[iii];
			}
		}, i_chunk));
	}
	hpx::wait_all(futs.begin(),futs.end());
	futs.clear();


	for( integer i_chunk = 0; i_chunk < N3; i_chunk += chunk_size ) {
		futs.push_back(hpx::async([&](integer this_i_chunk) {
			for (integer iii = this_i_chunk; iii != N3 && iii != this_i_chunk + chunk_size; ++iii) {
				auto ux_p = u0_p[iii];
				auto uy_p = u0_p[iii];
				auto uz_p = u0_p[iii];
				if (!is_on_edge[iii]) {
					simd_vector rho_h = Fourier.inverse_transform(u0_p[iii].rho());
					simd_vector fgx = Fourier.transform(rho_h*this_gx_h[iii]);
					simd_vector fgy = Fourier.transform(rho_h*this_gy_h[iii]);
					simd_vector fgz = Fourier.transform(rho_h*this_gz_h[iii]);
					apply_limiter(ux_p, u0_p[iii + dx_i], u0_p[iii - dx_i], fgx, XDIM);
					apply_limiter(uy_p, u0_p[iii + dy_i], u0_p[iii - dy_i], fgy, YDIM);
					apply_limiter(uz_p, u0_p[iii + dz_i], u0_p[iii - dz_i], fgz, ZDIM);
				}
				for (integer f = 0; f != NF; ++f) {
					u_p[iii][f] = minmod(minmod(ux_p[f], uy_p[f]), uz_p[f]);
					for( integer ppp = 0; ppp != P3; ++ppp) {
						if( std::abs(u_p[iii][f][ppp] - u0_p[iii][f][ppp]) > 1.0e-6) {
							++cnt;
						}
					}
				}
				const simd_vector rho_before_p = u0_p[iii].rho();
				const simd_vector rho_after_p = u_p[iii].rho();
				const simd_vector drho_h = Fourier.inverse_transform(rho_after_p - rho_before_p);
				u_p[iii].egas() -= Fourier.transform(drho_h * this_phi_h[iii]);
			}
		}, i_chunk));
	}
	hpx::wait_all(futs.begin(),futs.end());
	return cnt;
}

void grid::compute_flux(integer rk) {
	constexpr
	real hf = real(1) / real(2);
	const auto& u_p = U_p[rk];
	for (integer i = 0; i != N3; ++i) {
		for (integer f = 0; f != NF; ++f) {
			S_p[i][f] = real(0);
		}
	}
	std::list<hpx::future<void>> futs;
	for (integer i_chunk = NX * NX; i_chunk < N3; i_chunk += chunk_size) {
		futs.push_back(hpx::async([=]() {
			conserved_vars Uv_h(G3);
			conserved_vars Usp_h(G2);
			conserved_vars Usm_h(G2);
			for( integer i = i_chunk; i != i_chunk+chunk_size && i < N3; ++i) {
				for (integer dim = 0; dim != NDIM; ++dim) {
					face fcp = face(2 * dim + 1);
					face fcm = face(2 * dim);
					if (is_interior[i] || is_interior[i - d_i[dim]]) {
						if (is_interior[i] || is_interior[i - d_i[dim]]) {
							for (integer f = 0; f != NF; ++f) {
								Usm_h[f] = Fourier.surface_inverse_transform(fcp, u_p[i - d_i[dim]][f]);
								Usp_h[f] = Fourier.surface_inverse_transform(fcm, u_p[i][f]);
							}
							const simd_vector phi_m_p = Fourier.transform(phi_h[rk][i - d_i[dim]]);
							const simd_vector phi_p_p = Fourier.transform(phi_h[rk][i]);
							const simd_vector phi_sm_h = Fourier.surface_inverse_transform(fcp, phi_m_p);
							const simd_vector phi_sp_h = Fourier.surface_inverse_transform(fcm, phi_p_p);
							simd_vector phi_f = (phi_sm_h + phi_sp_h) * hf;
							auto surface_flux1 = KT_flux(Usm_h, Usp_h, phi_f, dimension(dim));
							for (integer f = 0; f != NF; ++f) {
								simd_vector pflux1 = Fourier.surface_transform(fcp, surface_flux1[f]) * hinv;
								simd_vector mflux1 = Fourier.surface_transform(fcm, surface_flux1[f]) * hinv;
								F_p[fcp][i - d_i[dim]][f] = pflux1;
								F_p[fcm][i][f] = mflux1;
							}
						} else {
							for (integer f = 0; f != NF; ++f) {
								F_p[fcp][i - d_i[dim]][f] = real(0);
								F_p[fcm][i][f] = real(0);
							}
						}
					}
					if (is_interior[i]) {
						for (integer f = 0; f != NF; ++f) {
							Uv_h[f] = Fourier.volume_inverse_transform(dimension(dim), u_p[i][f]);
						}
						const auto Vv_h = Uv_h.to_primitive();
						auto volume_flux = Uv_h.flux(Vv_h, dimension(dim));
						for (integer f = 0; f != NF; ++f) {
							S_p[i][f] += Fourier.volume_transform(dimension(dim), volume_flux[f]) * hinv;
						}
						simd_vector this_rho_p = u_p[i].rho();
						simd_vector rho_h = Fourier.inverse_transform(this_rho_p);
						S_p[i][s_i + dim] += Fourier.transform(rho_h * g_h[dimension(dim)][rk][i]);
						const simd_vector phi_vflux = phi_h[rk][i] * Uv_h[s_i + dim];
						S_p[i][egas_i] += Fourier.volume_transform(dimension(dim), phi_vflux) * hinv;
					}
				}
			}
		}));
	}
	hpx::wait_all(futs.begin(), futs.end());
}

void grid::compute_du(integer rk) {
	const real dV = real(8) * h * h * h;
	auto& du_dt_p = dU_dt_p[rk];
	for (integer f = 0; f != NF; ++f) {
		dO_dt[rk][f] = real(0);
	}
	std::list<hpx::future<void>> futs;
	for( integer i_chunk = 0; i_chunk < N3; i_chunk += chunk_size ) {
		futs.push_back(hpx::async([&](integer this_i_chunk) {
			for (integer iii = this_i_chunk; iii != N3 && iii != this_i_chunk + chunk_size; ++iii) {
				if (is_interior[iii]) {
					for (integer f = 0; f != NF; ++f) {
						const simd_vector src = S_p[iii][f];
						std::array<simd_vector*,NDIM> fp, fm;
						fp[0] = &F_p[XP][iii][f];
						fp[1] = &F_p[YP][iii][f];
						fp[2] = &F_p[ZP][iii][f];
						fm[0] = &F_p[XM][iii][f];
						fm[1] = &F_p[YM][iii][f];
						fm[2] = &F_p[ZM][iii][f];
						const simd_vector flux = *(fp[0]) + *(fp[1]) + *(fp[2]) - *(fm[0]) - *(fm[1]) - *(fm[2]);
						du_dt_p[iii][f] = src - flux;
					}
				}
				const simd_vector ep = Fourier.transform((Fourier.inverse_transform(du_dt_p[iii][rho_i]) * phi_h[rk][iii]));
				du_dt_p[iii][egas_i] -= ep;
			}
		}, i_chunk));
	}
	hpx::wait_all(futs.begin(),futs.end());
	for (integer iii = 0; iii != N3; ++iii) {
		if (is_interior[iii]) {
			for (integer f = 0; f != NF; ++f) {
				std::array<simd_vector*,NDIM> fp, fm;
				fp[0] = &F_p[XP][iii][f];
				fp[1] = &F_p[YP][iii][f];
				fp[2] = &F_p[ZP][iii][f];
				fm[0] = &F_p[XM][iii][f];
				fm[1] = &F_p[YM][iii][f];
				fm[2] = &F_p[ZM][iii][f];
				for( integer d = 0; d != NDIM; ++d) {
					if( !is_interior[iii-d_i[d]] ) {
						dO_dt[rk][f] -= (*fm[d])[0]*dV;
					}
					if( !is_interior[iii+d_i[d]] ) {
						dO_dt[rk][f] += (*fp[d])[0]*dV;
					}
				}
			}
		}
	}
}

void grid::compute_next_u(integer rk, real dt) {
	const auto& alpha = alpha_rk[NRK - 1];
	const auto& beta = beta_rk[NRK - 1];
	std::list<hpx::future<void>> futs;

	for( integer i_chunk = 0; i_chunk < N3; i_chunk += chunk_size ) {
		futs.push_back(hpx::async([&](integer this_i_chunk) {
			for (integer iii = this_i_chunk; iii != N3 && iii != this_i_chunk + chunk_size; ++iii) {
				for (integer f = 0; f != NF; ++f) {
					U_p[rk + 1][iii][f] = real(0);
					for (integer k = 0; k < rk + 1; ++k) {
						U_p[rk + 1][iii][f] += alpha[rk][k]*U_p[k][iii][f] + beta[rk][k] * dU_dt_p[k][iii][f] * dt;
					}
				}
				if (rk + 1 == NRK) {
					for (integer f = 0; f != NF; ++f) {
						U_p[0][iii][f] = U_p[rk + 1][iii][f];
					}
				}
			}
		}, i_chunk));
	}
	hpx::wait_all(futs.begin(),futs.end());
	for (integer f = 0; f != NF; ++f) {
		O[rk + 1][f] = real(0);
		for (integer k = 0; k < rk + 1; ++k) {
			O[rk + 1][f] += alpha[rk][k] * O[k][f] + beta[rk][k] * dO_dt[k][f] * dt	;
		}
	}
	if (rk + 1 == NRK) {
		for (integer f = 0; f != NF; ++f) {
			O[0][f] = O[rk + 1][f];
		}
	}
}


void grid::apply_limiter(conserved_vars& U0_p, const conserved_vars& UR_p, const conserved_vars& UL_p,
		const simd_vector fg_p, dimension dim) const {
	const dimension i1 = dim;
	const dimension i2 = ((dim == XDIM ? YDIM : XDIM));
	const dimension i3 = ((dim == ZDIM ? YDIM : ZDIM));
	const real rho = U0_p.rho()[0];
	const real u = U0_p.s(i1)[0] / rho;
	const real v = U0_p.s(i2)[0] / rho;
	const real w = U0_p.s(i3)[0] / rho;
	real ei = (U0_p.egas()[0] - rho * (u * u + v * v + w * w) / real(2));
	assert(U0_p.tau()[0] > real(0));
	assert(U0_p.rho()[0] > real(0));
	if (ei < dual_energy_switch2 * U0_p.egas()[0]) {
		ei = std::pow(U0_p.tau()[0], fgamma);
	}
	const real p = (fgamma - real(1)) * ei;
	const real c = std::sqrt(fgamma * p / rho);
	const real h0 = (p + U0_p.egas()[0]) / rho;

	auto char_decomp =
			[&](simd_vector& dC, const real& drho, const real& ds1, const real& ds2, const real& ds3, const real& degas, const real& dtau, const real& fg) {
				const auto du = (ds1 - u * drho) / rho;
				const auto dv = (ds2 - v * drho) / rho;
				const auto dw = (ds3 - w * drho) / rho;
				auto dp = (fgamma - real(1)) * (degas - rho * (u * du + v * dv + w * dw) - drho * (u * u + v * v + w * w) / real(2));
	//			dp -= fg*h;
				dC[0] = drho - dp / (c * c);
				dC[1] = dv;
				dC[2] = dw;
				dC[3] = du + dp / (rho * c);
				dC[4] = du - dp / (rho * c);
				dC[5] = dtau;
			};

	auto char_recomp =
			[&](simd_vector& dC, real& drho, real& ds1, real& ds2, real& ds3, real& degas, real& dtau, real fg) {
		//		dC[0] -= fg*h / (c*c);
		//		dC[3] += fg*h/(rho*c);
		//		dC[4] -= fg*h/(rho*c);
				drho = dC[0] + rho / (real(2) * c) * (dC[3] - dC[4]);
				ds1 = dC[0] * u + rho / (real(2) * c) * (dC[3] * (u + c) - dC[4] * (u - c));
				ds2 = dC[0] * v + rho * dC[1] + rho * v / (real(2) * c) * (dC[3]- dC[4]);
				ds3 = dC[0] * w + rho * dC[2] + rho * w / (real(2) * c) * (dC[3]- dC[4]);
				degas = dC[0] * (u * u + v * v + w * w) / real(2) +
				rho * (v * dC[1] + w * dC[2]) +
				rho / (real(2) * c) * (dC[3] * (h0 + u * c) - dC[4] * (h0 - u * c));
				dtau = dC[5];
			};

	simd_vector C0(NF), CR(NF), CL(NF);
	simd_vector CP(NF), CM(NF);
	integer l3[NDIM];
	integer l3p1[NDIM];
	integer & l = l3[XDIM];
	integer & m = l3[YDIM];
	integer & n = l3[ZDIM];
	integer & lp1 = l3p1[XDIM];
	integer & mp1 = l3p1[YDIM];
	integer & np1 = l3p1[ZDIM];
	for (l3[i2] = 0; l3[i2] < P - 1; ++l3[i2]) {
		l3p1[i2] = l3[i2];
		for (l3[i3] = 0; l3[i3] < P - 1 - l3[i2]; ++l3[i3]) {
			l3p1[i3] = l3[i3];
			for (l3[i1] = P - 2 - l3[i2] - l3[i3]; l3[i1] >= 0; --l3[i1]) {
				l3p1[i1] = l3[i1] + 1;
				const integer p0 = Fourier.pindex(l, m, n);
				const integer p1 = Fourier.pindex(lp1, mp1, np1);
				const real norm = real(2 * l3[i1] + 1);
				char_decomp(C0, U0_p.rho()[p1], U0_p.s(i1)[p1], U0_p.s(i2)[p1], U0_p.s(i3)[p1], U0_p.egas()[p1],
						U0_p.tau()[p1], fg_p[p0]);
				char_decomp(CR, //
						(UR_p.rho()[p0] - U0_p.rho()[p0]), //
						(UR_p.s(i1)[p0] - U0_p.s(i1)[p0]), //
						(UR_p.s(i2)[p0] - U0_p.s(i2)[p0]), //
						(UR_p.s(i3)[p0] - U0_p.s(i3)[p0]), //
						(UR_p.egas()[p0] - U0_p.egas()[p0]), (UR_p.tau()[p0] - U0_p.tau()[p0]), fg_p[p0]);
				char_decomp(
						CP, //
						(UR_p.rho()[p0] - U0_p.rho()[p0]) - norm * UR_p.rho()[p1], //
						(UR_p.s(i1)[p0] - U0_p.s(i1)[p0]) - norm * UR_p.s(i1)[p1], //
						(UR_p.s(i2)[p0] - U0_p.s(i2)[p0]) - norm * UR_p.s(i2)[p1], //
						(UR_p.s(i3)[p0] - U0_p.s(i3)[p0]) - norm * UR_p.s(i3)[p1], //
						(UR_p.egas()[p0] - U0_p.egas()[p0]) - norm * UR_p.egas()[p1],
						(UR_p.tau()[p0] - U0_p.tau()[p0]) - norm * UR_p.tau()[p1], fg_p[p0]);
				char_decomp(CL, //
						(U0_p.rho()[p0] - UL_p.rho()[p0]), //
						(U0_p.s(i1)[p0] - UL_p.s(i1)[p0]), //
						(U0_p.s(i2)[p0] - UL_p.s(i2)[p0]), //
						(U0_p.s(i3)[p0] - UL_p.s(i3)[p0]), //
						(U0_p.egas()[p0] - UL_p.egas()[p0]), (U0_p.tau()[p0] - UL_p.tau()[p0]), fg_p[p0]);
				char_decomp(
						CM, //
						(U0_p.rho()[p0] - UL_p.rho()[p0]) - norm * UL_p.rho()[p1], //
						(U0_p.s(i1)[p0] - UL_p.s(i1)[p0]) - norm * UL_p.s(i1)[p1], //
						(U0_p.s(i2)[p0] - UL_p.s(i2)[p0]) - norm * UL_p.s(i2)[p1], //
						(U0_p.s(i3)[p0] - UL_p.s(i3)[p0]) - norm * UL_p.s(i3)[p1], //
						(U0_p.egas()[p0] - UL_p.egas()[p0]) - norm * UL_p.egas()[p1],
						(U0_p.tau()[p0] - UL_p.tau()[p0]) - norm * UL_p.tau()[p1], fg_p[p0]);
				for (integer f = 0; f != NF; ++f) {
					const real clim1 = minmod(CL[f], CR[f]) / norm;
					const real clim2 = minmod(CP[f], CM[f]) / norm;
					C0[f] = std::max(C0[f], std::min(clim1, clim2));
					C0[f] = std::min(C0[f], std::max(clim1, clim2));
				}
				char_recomp(C0, U0_p.rho()[p1], U0_p.s(i1)[p1], U0_p.s(i2)[p1], U0_p.s(i3)[p1], U0_p.egas()[p1],
						U0_p.tau()[p1], fg_p[p0]);
			}
		}
	}

}


void grid::initialize(std::function<std::vector<real>(real, real, real)>&& func) {
	printf( "Initializing\n");
	const auto quad_point = Fourier.quadrature_points();
	std::vector<simd_vector> U(NF, simd_vector(G3));
	simd_vector gx_h(G3);
	simd_vector gy_h(G3);
	simd_vector gz_h(G3);
	for (integer i = 0; i != N3; ++i) {
		fflush(stdout);
		for (integer g = 0; g != G3; ++g) {
			const real x = cell_x[i] + quad_point[g][XDIM] * h;
			const real y = cell_y[i] + quad_point[g][YDIM] * h;
			const real z = cell_z[i] + quad_point[g][ZDIM] * h;
			const auto this_u = func(x, y, z);
			for (integer f = 0; f != NF; ++f) {
				U[f][g] = this_u[f];
			}
		//	star_force(x,y,z, gx_h_analytic[i][g], gy_h_analytic[i][g], gz_h_analytic[i][g]);
		}
		for (integer f = 0; f != NF; ++f) {
			U_p[0][i][f] = Fourier.transform(U[f]);
		}
		if( i % (NX*NX) == 0 ) {
			printf( ".");
		}
	}
	printf( "\nDone\n");
	enforce_boundaries(0);
	compute_fmm(0);
	enforce_positivity(0);
	project(0);
	t = real(0);
	enforce_boundaries(0);
	enforce_positivity(0);
	grid_amax = compute_cfl_condition();
	output_dt = 1.0e-2;
	output_cnt = 0;
	compute_fmm(0);
	output("X.0.silo");
	++output_cnt;
	printf("t = 0.0 *\n");
}


void grid::initialize(const char* fname) {
	printf( "Loading checkpoint\n");
	load(fname);
	output_dt = 1.0e-2;
	output_cnt = integer(t / output_dt) + 1;


}

void grid::compute_analytic(std::function<std::vector<real>(real, real, real, real)>&& func, real t) {
	const auto quad_point = Fourier.quadrature_points();
	for (integer i = 0; i != N3; ++i) {
		fflush(stdout);
		for (integer g = 0; g != G3; ++g) {
			const real x = cell_x[i] + quad_point[g][XDIM] * h;
			const real y = cell_y[i] + quad_point[g][YDIM] * h;
			const real z = cell_z[i] + quad_point[g][ZDIM] * h;
			const auto this_u = func(x, y, z, t);
			for (integer f = 0; f != NF; ++f) {
				U_h_analytic[i][f][g] = this_u[f];
			}
		}
	}
}

struct node_point {
	xpoint pt;
	integer index;
	bool operator==(const node_point& other) const {
		return xpoint_eq(other.pt, pt);
	}
	bool operator<(const node_point& other) const {
		bool rc = false;
		for (integer d = 0; d != NDIM; ++d) {
			if (!float_eq(pt[d], other.pt[d])) {
				rc = (pt[d] < other.pt[d]);
				break;
			}
		}
		return rc;
	}
};

void grid::output(const char* filename) const {
	if( system( "cp hello.chk.dat goodbye.chk.dat\n") == 0 ) {

	}
	save("hello.chk.dat");

	constexpr
	integer vertex_order[8] = { 0, 1, 3, 2, 4, 5, 7, 6 };
	std::set<node_point> node_list;

	const integer rk = 0;
	std::list < integer > zone_list;
	const auto& quad_weights = Fourier.quadrature_weights();
	std::vector < real > quad_points(P + 1);
	quad_points[0] = -h;
	for (integer p = 0; p != P; ++p) {
		quad_points[p + 1] = quad_points[p] + quad_weights[p] * h;
	}
	constexpr integer this_bw = BW;
	for (integer i = this_bw; i != NX - this_bw; ++i) {
		for (integer j = this_bw; j != NX - this_bw; ++j) {
			for (integer k = this_bw; k != NX - this_bw; ++k) {
				for (integer gx = 0; gx != P; ++gx) {
					for (integer gy = 0; gy != P; ++gy) {
						for (integer gz = 0; gz != P; ++gz) {
							const integer iii = NX * NX * i + NX * j + k;
							for (integer ci = 0; ci != NVERTEX; ++ci) {
								const integer vi = vertex_order[ci];
								const integer xi = (vi >> 0) & 1;
								const integer yi = (vi >> 1) & 1;
								const integer zi = (vi >> 2) & 1;
								node_point this_x;
								this_x.pt[0] = cell_x[iii] + quad_points[gx + xi];
								this_x.pt[1] = cell_y[iii] + quad_points[gy + yi];
								this_x.pt[2] = cell_z[iii] + quad_points[gz + zi];
								auto iter = node_list.find(this_x);
								integer index;
								if (iter != std::end(node_list)) {
									index = iter->index;
								} else {
									index = node_list.size();
									this_x.index = index;
									node_list.insert(this_x);
								}
								zone_list.push_back(index);
							}
						}
					}
				}
			}
		}
	}

	const int nzones = zone_list.size() / NVERTEX;
	std::vector<int> zone_nodes(nzones * NVERTEX);
	integer index = 0;
	for (auto iter = std::begin(zone_list); iter != std::end(zone_list); ++iter) {
		zone_nodes[index] = *iter;
		++index;
	}

	const int nnodes = node_list.size();
	std::vector<double> x_coord(nnodes);
	std::vector<double> y_coord(nnodes);
	std::vector<double> z_coord(nnodes);
	std::array<double*, NDIM> node_coords = { x_coord.data(), y_coord.data(), z_coord.data() };
	for (auto iter = std::begin(node_list); iter != std::end(node_list); ++iter) {
		x_coord[iter->index] = iter->pt[0];
		y_coord[iter->index] = iter->pt[1];
		z_coord[iter->index] = iter->pt[2];
	}

	constexpr
	int nshapes = 1;
	int shapesize[1] = { NVERTEX };
	int shapetype[1] = { DB_ZONETYPE_HEX };
	int shapecnt[1] = { nzones };
	const char* coord_names[NDIM] = { "x", "y", "z" };

	DBfile *db = DBCreateReal(filename, DB_CLOBBER, DB_LOCAL, "Euler Mesh", DB_PDB);
	DBPutZonelist2(db, "zones", nzones, int(NDIM), zone_nodes.data(), nzones * NVERTEX, 0, 0, 0, shapetype, shapesize,
			shapecnt, nshapes, nullptr);
	DBPutUcdmesh(db, "mesh", int(NDIM), coord_names, node_coords.data(), nnodes, nzones, "zones", nullptr, DB_DOUBLE,
			nullptr);

	const char* field_names[] = { "rho", "egas", "tau", "sx", "sy", "sz", "phi", "gx", "gy", "gz", "rho_analytic",
			"egas_analytic", "tau_analytic", "sx_analytic", "sy_analytic", "sz_analytic", "phi_analytic", "gx_analytic",
			"gy_analytic", "gz_analytic" };
	constexpr
	integer I3 = std::pow(NX-2*this_bw,3);
	std::vector<double> rho(I3 * G3), sx(I3 * G3), sy(I3 * G3), sz(I3 * G3), egas(I3 * G3), tau(I3 * G3), phi(I3 * G3),
			gx(I3 * G3), gy(I3 * G3), gz(I3 * G3), phi_analytic(I3 * G3), gx_analytic(I3 * G3), gy_analytic(I3 * G3),
			gz_analytic(I3 * G3), rho_analytic(I3 * G3), sx_analytic(I3 * G3), sy_analytic(I3 * G3), sz_analytic(
					I3 * G3), egas_analytic(I3 * G3), tau_analytic(I3 * G3);
	std::array<double*, 2 * (NF + 1 + NDIM)> u_data = { rho.data(), sx.data(), sy.data(), sz.data(), egas.data(),
			tau.data(), phi.data(), gx.data(), gy.data(), gz.data(), rho_analytic.data(), egas_analytic.data(),
			tau_analytic.data(), sx_analytic.data(), sy_analytic.data(), sz_analytic.data(), phi_analytic.data(),
			gx_analytic.data(), gy_analytic.data(), gz_analytic.data() };
	index = 0;
	for (integer i = this_bw; i != NX - this_bw; ++i) {
		for (integer j = this_bw; j != NX - this_bw; ++j) {
			for (integer k = this_bw; k != NX - this_bw; ++k) {
				const integer iii = NX * NX * i + NX * j + k;
				for (integer f = 0; f != NF; ++f) {
					simd_vector U_h = Fourier.inverse_transform(U_p[rk][iii][f]);
					for (integer ggg = 0; ggg != G3; ++ggg) {
						u_data[f][index + ggg] = U_h[ggg];
						u_data[f + (NF + 1 + NDIM)][index + ggg] = U_h_analytic[iii][f][ggg];
					}
				}
				simd_vector phi_h_analytic = Fourier.inverse_transform(phi_p_analytic[iii]);
				for (integer ggg = 0; ggg != G3; ++ggg) {
					u_data[NF][index + ggg] = phi_h[0][iii][ggg];
					u_data[NF + 1 + XDIM][index + ggg] = gx_h[0][iii][ggg];
					u_data[NF + 1 + YDIM][index + ggg] = gy_h[0][iii][ggg];
					u_data[NF + 1 + ZDIM][index + ggg] = gz_h[0][iii][ggg];
					u_data[2 * NF + 1 + NDIM][index + ggg] = phi_h_analytic[ggg];
					u_data[2 * NF + 2 + NDIM + XDIM][index + ggg] = gx_h_analytic[iii][ggg];
					u_data[2 * NF + 2 + NDIM + YDIM][index + ggg] = gy_h_analytic[iii][ggg];
					u_data[2 * NF + 2 + NDIM + ZDIM][index + ggg] = gz_h_analytic[iii][ggg];
				}
				index += G3;
			}
		}
	}
	for (int f = 0; f != 2 * (NF + 1 + NDIM); ++f) {
		DBPutUcdvar1(db, field_names[f], "mesh", u_data[f], nzones, nullptr, 0, DB_DOUBLE, DB_ZONECENT, nullptr);
	}

	DBClose(db);
}

void grid::compute_multipoles(integer rk) {
	integer lev = 0;
	for (integer inx = INX; inx > 1; inx >>= 1) {
		const integer nxp = inx + 2 * BW;
		const integer nxc = (2 * inx) + 2 * BW;
		std::fill(std::begin(rho_l[lev]), std::end(rho_l[lev]), simd_vector(real(0), L2));
		std::list<hpx::future<void>> futs;
		for (integer ip = BW; ip != nxp - BW; ++ip) {
			futs.push_back(hpx::async([=](){
				for (integer jp = BW; jp != nxp - BW; ++jp) {
					for (integer kp = BW; kp != nxp - BW; ++kp) {
						if( lev != 0 ) {
							for (integer ci = 0; ci != NVERTEX; ++ci) {
								const integer ic = (2 * ip - BW) + ((ci >> 0) & 1);
								const integer jc = (2 * jp - BW) + ((ci >> 1) & 1);
								const integer kc = (2 * kp - BW) + ((ci >> 2) & 1);
								const integer iiic = nxc * nxc * ic + nxc * jc + kc;
								const integer iiip = nxp * nxp * ip + nxp * jp + kp;
								const simd_vector this_rho = rho_l[lev - 1][iiic];
								const real this_dx = real(1 << (lev - 1)) * h;
								rho_l[lev][iiip] += Fourier.m2m_transform(ci, this_rho, this_dx );
							}
						} else {
							const integer iii = nxp * nxp * ip + nxp * jp + kp;
							simd_vector this_rho = U_p[rk][iii].rho();
							this_rho[0] -= rho_floor;
							if (is_interior[iii]) {
								rho_l[lev][iii] = Fourier.p2m_transform(this_rho, h);
							} else {
								rho_l[lev][iii] = real(0);
							}
						}
					}
				}
			}));
		}
		hpx::wait_all(futs.begin(), futs.end());
		++lev;
	}
	lev = 1;
}

void grid::compute_interactions(integer rk) {
	integer lev = nlevel - 1;
	std::list<hpx::future<void>> futs;
	for (integer inx = 2; inx <= INX; inx <<= 1) {
		const integer nx = inx + 2 * BW;
		for (integer i0 = 1; i0 != nx-1 ; ++i0) {
			futs.push_back(hpx::async([=]() {
				for (integer j0 = 1; j0 != nx-1; ++j0) {
					for (integer k0 = 1; k0 != nx-1; ++k0) {
						const integer iii0 = i0 * nx * nx + j0 * nx + k0;
						phi_l[lev][iii0] = real(0);
						const integer imin = std::max(integer(0),2 * ((i0 / 2) - 1));
						const integer imax = std::min(integer(nx-1),2 * ((i0 / 2) + 1) + 1);
						const integer jmin = std::max(integer(0),2 * ((j0 / 2) - 1));
						const integer jmax = std::min(integer(nx-1),2 * ((j0 / 2) + 1) + 1);
						const integer kmin = std::max(integer(0),2 * ((k0 / 2) - 1));
						const integer kmax = std::min(integer(nx-1),2 * ((k0 / 2) + 1) + 1);
						for (integer i1 = imin; i1 <= imax; ++i1) {
							for (integer j1 = jmin; j1 <= jmax; ++j1) {
								for (integer k1 = kmin; k1 <= kmax; ++k1) {
									const integer iii1 = i1 * nx * nx + j1 * nx + k1;
									integer max_dist = std::max(std::abs(i0 - i1), std::abs(j0 - j1));
									max_dist = std::max(std::abs(k0 - k1), max_dist);
									if (max_dist > 1  && lev != 0 ) {
										phi_l[lev][iii0] += Fourier.m2l_transform(i0 - i1, j0 - j1, k0 - k1,
												rho_l[lev][iii1], real(1 << lev) * h);
									}
								}
							}
						}
					}
				}
			}));
		}
		--lev;
	}
	hpx::wait_all(futs.begin(),futs.end());
}

void grid::expand_phi(integer rk) {
	integer lev = nlevel - 1;
	for (integer inx = 4; inx <= INX; inx <<= 1) {
		--lev;
		const integer nxp = (inx / 2) + 2 * BW;
		const integer nxc = inx + 2 * BW;
		std::list<hpx::future<void>> futs;
		for (integer ip = 1; ip != nxp - 1; ++ip) {
			futs.push_back(hpx::async([=]() {
				for (integer jp = 1; jp != nxp - 1; ++jp) {
					for (integer kp = 1; kp != nxp - 1; ++kp) {
						const integer iiip = nxp * nxp * ip + nxp * jp + kp;
						for (integer ci = 0; ci != NVERTEX; ++ci) {
							const integer ic = (2 * ip - BW) + ((ci >> 0) & 1);
							const integer jc = (2 * jp - BW) + ((ci >> 1) & 1);
							const integer kc = (2 * kp - BW) + ((ci >> 2) & 1);
							const integer iiic = nxc * nxc * ic + nxc * jc + kc;
							phi_l[lev][iiic] += Fourier.l2l_transform(ci, phi_l[lev + 1][iiip], real(1 << lev) * h);
						}
					}
				}
			}));
		}
		hpx::wait_all(futs.begin(),futs.end());
	}
}

void grid::compute_force(integer rk) {

	std::list<hpx::future<void>> futs;
	for (integer i1 = 1; i1 != NX-1 ; ++i1) {
		futs.push_back(hpx::async([=]() {
			for (integer j1 = 1; j1 != NX-1 ; ++j1) {
				for (integer k1 = 1; k1 != NX-1; ++k1) {
					const integer iii1 = i1 * dx_i + j1 * dy_i + k1 * dz_i;
					phi_h[rk][iii1] = real(0);
					gx_h[rk][iii1] = real(0);
					gy_h[rk][iii1] = real(0);
					gz_h[rk][iii1] = real(0);
					const integer imin = std::max(integer(0), 2*(i1/2 - 1));
					const integer jmin = std::max(integer(0), 2*(j1/2 - 1));
					const integer kmin = std::max(integer(0), 2*(k1/2 - 1));
					const integer imax = std::min(integer(NX-1), 2*(i1/2 + 1)+1);
					const integer jmax = std::min(integer(NX-1), 2*(j1/2 + 1)+1);
					const integer kmax = std::min(integer(NX-1), 2*(k1/2 + 1)+1);
					for (integer i2 = imin; i2 <= imax; ++i2) {
						for (integer j2 = jmin; j2 <= jmax; ++j2) {
							for (integer k2 = kmin; k2 <= kmax; ++k2) {
								const integer iii2 = i2 * dx_i + j2 * dy_i + k2 * dz_i;
								if (is_interior[iii2]) {
									const integer di = i1 - i2;
									const integer dj = j1 - j2;
									const integer dk = k1 - k2;
									phi_h[rk][iii1] -= Fourier.p2_phi_transform(di, dj, dk, U_p[rk][iii2].rho(), h);
									gx_h[rk][iii1] -= Fourier.p2_gx_transform(di, dj, dk, U_p[rk][iii2].rho(), h) * hinv;
									gy_h[rk][iii1] -= Fourier.p2_gy_transform(di, dj, dk, U_p[rk][iii2].rho(), h) * hinv;
									gz_h[rk][iii1] -= Fourier.p2_gz_transform(di, dj, dk, U_p[rk][iii2].rho(), h) * hinv;
								}
							}
						}
					}
					phi_h[rk][iii1] -= Fourier.l2p_transform(phi_l[0][iii1], h);
					gx_h[rk][iii1] -= Fourier.dl2p_transform_dx(XDIM, phi_l[0][iii1], h) * hinv;
					gy_h[rk][iii1] -= Fourier.dl2p_transform_dx(YDIM, phi_l[0][iii1], h) * hinv;
					gz_h[rk][iii1] -= Fourier.dl2p_transform_dx(ZDIM, phi_l[0][iii1], h) * hinv;
				}
			}
		}));
	}
	hpx::wait_all(futs.begin(), futs.end());
}

void grid::diagnostics(real t) {
	std::vector<real> sum_in(NF,real(0));
	const real dV = real(8) * h * h * h;
//	real gx_err = real(0);
//	real gx_norm = real(0);
	for (integer iii = 0; iii != N3; ++iii) {
		if (is_interior[iii]) {
			const simd_vector rho_h = Fourier.inverse_transform((U_p[0][iii].rho()));
			const auto rho_phi = Fourier.transform(phi_h[0][iii]*rho_h)[0];
			for( integer f = 0; f != NF; ++f) {
				sum_in[f] += U_p[0][iii][f][0]*dV;
			}
			sum_in[egas_i] += rho_phi/real(2)*dV;
	//		gx_err += std::pow(gx_h[0][iii] - gx_h_analytic[iii],real(2)).sum();
//			gx_norm += std::pow(gx_h_analytic[iii],real(2)).sum();
		}
	}
	FILE* fp = fopen("diag.dat", "at");
	fprintf(fp, "%e ", double(t));
	for( integer f = 0; f != NF; ++f) {
		fprintf(fp, "%.9e %.9e ", double(sum_in[f]), double(sum_in[f]+O[0][f]));
	}
	fprintf(fp, "\n");
	fclose(fp);
//	printf( "%e\n", gx_err/gx_norm);
//	sleep(5);
//	abort();
}

