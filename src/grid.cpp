/*
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

real minmod(real a, real b) {
	return (std::copysign(real(1), a) + std::copysign(real(1), b)) * std::min(std::abs(a), std::abs(b)) / real(2);
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

grid::grid() :
		dx(real(2) / real(INX)), Fourier(P) {
	U_p.resize(NRK + 1, std::vector<conserved_vars>(N3, conserved_vars(P3)));
	F_p.resize(NFACE, std::vector<std::vector<simd_vector>>(N3, std::vector<simd_vector>(NF, simd_vector(P3))));
	S_p.resize(N3, std::vector<simd_vector>(NF, simd_vector(P3)));
	dU_dt_p.resize(NRK, std::vector<std::vector<simd_vector>>(N3, std::vector<simd_vector>(NF, simd_vector(P3))));
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
				cell_x[index] = (real(i - BW - INX / 2) + real(1) / real(2)) * dx;
				cell_y[index] = (real(j - BW - INX / 2) + real(1) / real(2)) * dx;
				cell_z[index] = (real(k - BW - INX / 2) + real(1) / real(2)) * dx;
			}
		}
	}
}

real grid::enforce_positivity(integer rk) {
	real amax = real(0);
	real eps_rho = std::numeric_limits<real>::max();
	real eps_tau = std::numeric_limits<real>::max();
	const real cfl = Fourier.lobatto_edge_weight() / real(4);
	const integer G3L = Fourier.lobatto_point_count();

	for (integer iii = 0; iii != N3; ++iii) {
		if (!is_on_edge[iii]) {
			auto& u_p = U_p[rk][iii];
			conserved_vars u_h(G3);
			for (integer f = 0; f != NF; ++f) {
				u_h[f] = Fourier.inverse_transform(u_p[f]);
			}
			const auto egas = u_h.egas();
			const auto rho = u_h.rho();
			const auto u = simd_vector(u_h.s(XDIM) / rho);
			const auto v = simd_vector(u_h.s(YDIM) / rho);
			const auto w = simd_vector(u_h.s(ZDIM) / rho);
			const auto ei = simd_vector(egas - rho * (u * u + v * v + w * w) / real(2));
			for (integer ggg = 0; ggg != G3; ++ggg) {
				if ((ei[ggg] > dual_energy_switch1) * egas[ggg] && (ei[ggg] > real(0))) {
					u_h.tau()[ggg] = std::pow(ei[ggg], real(1) / fgamma);
				}
			}

			u_p.tau() = Fourier.transform(u_h.tau());
			eps_rho = std::min(eps_rho, u_p.rho()[0]);
			eps_tau = std::min(eps_tau, u_p.tau()[0]);
		}
	}
	eps_rho /= real(2);
	eps_tau /= real(2);
	for (integer iii = 0; iii != N3; ++iii) {
		if (!is_on_edge[iii]) {
			std::vector<conserved_vars> u_h(NDIM, conserved_vars(G3L));
			auto& u_p = U_p[rk][iii];
			for (integer dim = 0; dim != NDIM; ++dim) {
				for (integer f = 0; f != NF; ++f) {
					u_h[dim][f] = Fourier.lobatto_inverse_transform(u_p[f], dimension(dim));
				}
			}
			real rho_min = std::numeric_limits<real>::max();
			real tau_min = std::numeric_limits<real>::max();
			for (integer dim = 0; dim != NDIM; ++dim) {
				rho_min = std::min(rho_min, u_h[dim].rho().min());
				tau_min = std::min(tau_min, u_h[dim].tau().min());
			}
			const real this_rho = u_p.rho()[0];
			const real this_tau = u_p.tau()[0];
			const real theta1 = this_rho != rho_min ? (this_rho - eps_rho) / (this_rho - rho_min) : real(1);
			const real theta2 = this_tau != tau_min ? (this_tau - eps_tau) / (this_tau - tau_min) : real(1);
			const real theta = std::min(theta1, theta2);
			if (theta < real(1)) {
				for (integer f = 0; f != NF; ++f) {
					for (integer ppp = 1; ppp < P3; ++ppp) {
						u_p[f][ppp] *= theta;
					}
					for (integer dim = 0; dim != NDIM; ++dim) {
						u_h[dim][f] = Fourier.lobatto_inverse_transform(u_p[f], dimension(dim));
					}
				}
			}
			assert(u_h[XDIM].tau().min() > real(0));
			assert(u_h[YDIM].tau().min() > real(0));
			assert(u_h[ZDIM].tau().min() > real(0));
			for (integer dim = 0; dim != NDIM; ++dim) {
				const auto prims = u_h[dim].to_primitive();
				const real a_x = (std::abs(prims.v(XDIM)) + prims.c()).max();
				const real a_y = (std::abs(prims.v(YDIM)) + prims.c()).max();
				const real a_z = (std::abs(prims.v(ZDIM)) + prims.c()).max();
				amax = std::max(amax, a_x + a_y + a_z);
			}
		}
	}
	return amax / cfl;
}

void grid::enforce_boundaries(integer rk) {
	auto& u = U_p[rk];

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
							if (l == 0) {
								u[im1a][f][pppx] = u[iiia][f][pppx];
								u[im2a][f][pppx] = u[iiia][f][pppx];
								u[ip1b][f][pppx] = u[iiib][f][pppx];
								u[ip2b][f][pppx] = u[iiib][f][pppx];
							} else {
								u[im1a][f][pppx] = real(0);
								u[im2a][f][pppx] = real(0);
								u[ip1b][f][pppx] = real(0);
								u[ip2b][f][pppx] = real(0);
							}
						}
					}
				}
			}
			if (u[iiia][s_i + XDIM][0] < real(0)) {
				u[im1a][s_i + XDIM][0] = real(0);
				u[im2a][s_i + XDIM][0] = real(0);
			}
			if (u[iiib][s_i + XDIM][0] > real(0)) {
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
							if (l == 0) {
								u[jm1a][f][pppy] = u[jjja][f][pppy];
								u[jm2a][f][pppy] = u[jjja][f][pppy];
								u[jp1b][f][pppy] = u[jjjb][f][pppy];
								u[jp2b][f][pppy] = u[jjjb][f][pppy];
							} else {
								u[jm1a][f][pppy] = real(0);
								u[jm2a][f][pppy] = real(0);
								u[jp1b][f][pppy] = real(0);
								u[jp2b][f][pppy] = real(0);
							}
						}
					}
				}
			}
			if (u[jjja][s_i + YDIM][0] < real(0)) {
				u[jm1a][s_i + YDIM][0] = real(0);
				u[jm2a][s_i + YDIM][0] = real(0);
			}
			if (u[jjjb][s_i + YDIM][0] > real(0)) {
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
							if (l == 0) {
								u[km1a][f][pppz] = u[kkka][f][pppz];
								u[km2a][f][pppz] = u[kkka][f][pppz];
								u[kp1b][f][pppz] = u[kkkb][f][pppz];
								u[kp2b][f][pppz] = u[kkkb][f][pppz];
							} else {
								u[km1a][f][pppz] = real(0);
								u[km2a][f][pppz] = real(0);
								u[kp1b][f][pppz] = real(0);
								u[kp2b][f][pppz] = real(0);
							}
						}
					}
				}
			}
			if (u[kkka][s_i + ZDIM][0] < real(0)) {
				u[km1a][s_i + ZDIM][0] = real(0);
				u[km2a][s_i + ZDIM][0] = real(0);
			}
			if (u[kkkb][s_i + ZDIM][0] > real(0)) {
				u[kp1b][s_i + ZDIM][0] = real(0);
				u[kp2b][s_i + ZDIM][0] = real(0);
			}
		}
	}
}

void grid::project(integer rk) {
	auto& u_p = U_p[rk];
	auto ux_p = U_p[rk];
	auto uy_p = U_p[rk];
	auto uz_p = U_p[rk];
	for (integer iii = 0; iii != N3; ++iii) {
		if (!is_on_edge[iii]) {
			apply_limiter(ux_p[iii], u_p[iii + dx_i], u_p[iii - dx_i], XDIM);
			apply_limiter(uy_p[iii], u_p[iii + dy_i], u_p[iii - dy_i], YDIM);
			apply_limiter(uz_p[iii], u_p[iii + dz_i], u_p[iii - dz_i], ZDIM);
		}
	}
	for (integer iii = 0; iii != N3; ++iii) {
		if (!is_on_edge[iii]) {
			for (integer f = 0; f != NF; ++f) {
				u_p[iii][f] = minmod(ux_p[iii][f], minmod(uy_p[iii][f], uz_p[iii][f]));
			}
		}
	}
}

void grid::compute_flux(integer rk) {
	constexpr real hf = real(1) / real(2);
	const auto& u_p = U_p[rk];
	for (integer i = 0; i != N3; ++i) {
		for (integer f = 0; f != NF; ++f) {
			S_p[i][f] = real(0);
		}
	}
	for (integer dim = 0; dim != NDIM; ++dim) {
		face fcp = face(2 * dim + 1);
		face fcm = face(2 * dim);
		conserved_vars Uv_h(G3);
		conserved_vars Usp_h(G2);
		conserved_vars Usm_h(G2);
		for (integer i = NX * NX; i != N3; ++i) {
			if (is_interior[i] || is_interior[i - d_i[dim]]) {
				for (integer f = 0; f != NF; ++f) {
					Usm_h[f] = Fourier.surface_inverse_transform(fcp, u_p[i - d_i[dim]][f]);
					Usp_h[f] = Fourier.surface_inverse_transform(fcm, u_p[i][f]);
				}
				const auto Vsp_h = Usp_h.to_primitive();
				const auto Vsm_h = Usm_h.to_primitive();
				const real ap = simd_vector(Vsp_h.c() + std::abs(Vsp_h.v(dimension(dim)))).max();
				const real am = simd_vector(Vsm_h.c() + std::abs(Vsm_h.v(dimension(dim)))).max();
				const real a = std::max(ap, am);
				auto fp = Usp_h.flux(Vsp_h, dimension(dim));
				auto fm = Usm_h.flux(Vsm_h, dimension(dim));
				std::vector<simd_vector> surface_flux(NF, simd_vector(G2));
				for (integer f = 0; f != NF; ++f) {
					surface_flux[f] = ((fp[f] + fm[f]) - a * (Usp_h[f] - Usm_h[f])) * hf;
				}
				for (integer f = 0; f != NF; ++f) {
					F_p[fcp][i - d_i[dim]][f] = Fourier.surface_transform(fcp, surface_flux[f]);
					F_p[fcm][i][f] = Fourier.surface_transform(fcm, surface_flux[f]);
				}
			}
		}
		for (integer i = 0; i != N3; ++i) {
			if (is_interior[i]) {
				for (integer f = 0; f != NF; ++f) {
					Uv_h[f] = Fourier.volume_inverse_transform(dimension(dim), u_p[i][f]);
				}
				const auto Vv_h = Uv_h.to_primitive();
				auto volume_flux = Uv_h.flux(Vv_h, dimension(dim));
				for (integer f = 0; f != NF; ++f) {
					S_p[i][f] += Fourier.volume_transform(dimension(dim), volume_flux[f]);
				}
			}
		}
	}
}

void grid::compute_du(integer rk) {
	const real dxinv = real(1) / real(dx);
	auto& du_dt_p = dU_dt_p[rk];
	for (integer i = 0; i != N3; ++i) {
		if (is_interior[i]) {
			for (integer f = 0; f != NF; ++f) {
				du_dt_p[i][f] = S_p[i][f] * real(2) * dxinv;
				du_dt_p[i][f] -= (F_p[XP][i][f] - F_p[XM][i][f]) * real(2) * dxinv;
				du_dt_p[i][f] -= (F_p[YP][i][f] - F_p[YM][i][f]) * real(2) * dxinv;
				du_dt_p[i][f] -= (F_p[ZP][i][f] - F_p[ZM][i][f]) * real(2) * dxinv;
			}
		}
	}
}

void grid::compute_next_u(integer rk, real dt) {
	const auto& alpha = alpha_rk[NRK - 1];
	const auto& beta = beta_rk[NRK - 1];
	for (integer i = 0; i != N3; ++i) {
		if (is_interior[i]) {
			for (integer f = 0; f != NF; ++f) {
				U_p[rk + 1][i][f] = real(0);
				for (integer k = 0; k < rk + 1; ++k) {
					U_p[rk + 1][i][f] += alpha[rk][k] * U_p[k][i][f] + beta[rk][k] * dU_dt_p[k][i][f] * dt;
				}
			}
		}
	}
	if (rk + 1 == NRK) {
		for (integer i = 0; i != N3; ++i) {
			if (is_interior[i]) {
				for (integer f = 0; f != NF; ++f) {
					U_p[0][i][f] = U_p[rk + 1][i][f];
				}
			}
		}
	}
}

void grid::apply_limiter(conserved_vars& U0_p, const conserved_vars& UR_p, const conserved_vars& UL_p,
		dimension dim) const {
	const dimension i1 = dim;
	const dimension i2 = ((dim == XDIM ? YDIM : XDIM));
	const dimension i3 = ((dim == ZDIM ? YDIM : ZDIM));
	const real rho = U0_p.rho()[0];
	const real u = U0_p.s(i1)[0] / rho;
	const real v = U0_p.s(i2)[0] / rho;
	const real w = U0_p.s(i3)[0] / rho;
	real ei = (U0_p.egas()[0] - rho * (u * u + v * v + w * w) / real(2));
	ei = std::max(ei, real(0));
	if (ei < dual_energy_switch2 * U0_p.egas()[0]) {
		ei = std::pow(U0_p.tau()[0], fgamma);
	}
	const real p = (fgamma - real(1)) * ei;
	const real c = std::sqrt(fgamma * p / rho);
	const real h = (p + U0_p.egas()[0]) / rho;

	auto char_decomp =
			[&](simd_vector& dC, const real& drho, const real& ds1, const real& ds2, const real& ds3, const real& degas, const real& dtau) {
				const auto du = (ds1 - u * drho) / rho;
				const auto dv = (ds2 - v * drho) / rho;
				const auto dw = (ds3 - w * drho) / rho;
				const auto dp = (fgamma - real(1)) * (degas - rho * (u * du + v * dv + w * dw) - drho * (u * u + v * v + w * w) / real(2));
				dC[0] = drho - dp / (c * c);
				dC[1] = dv;
				dC[2] = dw;
				dC[3] = du + dp / (rho * c);
				dC[4] = du - dp / (rho * c);
				dC[5] = dtau;
			};

	auto char_recomp =
			[&](const simd_vector& dC, real& drho, real& ds1, real& ds2, real& ds3, real& degas, real& dtau) {
				drho = dC[0] + rho / (real(2) * c) * (dC[3] - dC[4]);
				ds1 = dC[0] * u + rho / (real(2) * c) * (dC[3] * (u + c) - dC[4] * (u - c));
				ds2 = dC[0] * v + rho * dC[1] + rho * v / (real(2) * c) * (dC[3]- dC[4]);
				ds3 = dC[0] * w + rho * dC[2] + rho * w / (real(2) * c) * (dC[3]- dC[4]);
				degas = dC[0] * (u * u + v * v + w * w) / real(2) +
				rho * (v * dC[1] + w * dC[2]) +
				rho / (real(2) * c) * (dC[3] * (h + u * c) - dC[4] * (h - u * c));
				dtau = dC[5];
			};

	simd_vector C0(NF), CR(NF), CL(NF);
	simd_vector CP(NF), CM(NF);
	integer l3[NDIM];
	integer l3p1[NDIM];
	integer& l = l3[XDIM];
	integer& m = l3[YDIM];
	integer& n = l3[ZDIM];
	integer& lp1 = l3p1[XDIM];
	integer& mp1 = l3p1[YDIM];
	integer& np1 = l3p1[ZDIM];
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
						U0_p.tau()[p1]);
				char_decomp(CR, //
						(UR_p.rho()[p0] - U0_p.rho()[p0]), //
						(UR_p.s(i1)[p0] - U0_p.s(i1)[p0]), //
						(UR_p.s(i2)[p0] - U0_p.s(i2)[p0]), //
						(UR_p.s(i3)[p0] - U0_p.s(i3)[p0]), //
						(UR_p.egas()[p0] - U0_p.egas()[p0]), (UR_p.tau()[p0] - U0_p.tau()[p0]));
				char_decomp(
						CP, //
						(UR_p.rho()[p0] - U0_p.rho()[p0]) - norm * UR_p.rho()[p1], //
						(UR_p.s(i1)[p0] - U0_p.s(i1)[p0]) - norm * UR_p.s(i1)[p1], //
						(UR_p.s(i2)[p0] - U0_p.s(i2)[p0]) - norm * UR_p.s(i2)[p1], //
						(UR_p.s(i3)[p0] - U0_p.s(i3)[p0]) - norm * UR_p.s(i3)[p1], //
						(UR_p.egas()[p0] - U0_p.egas()[p0]) - norm * UR_p.egas()[p1],
						(UR_p.tau()[p0] - U0_p.tau()[p0]) - norm * UR_p.tau()[p1]);
				char_decomp(CL, //
						(U0_p.rho()[p0] - UL_p.rho()[p0]), //
						(U0_p.s(i1)[p0] - UL_p.s(i1)[p0]), //
						(U0_p.s(i2)[p0] - UL_p.s(i2)[p0]), //
						(U0_p.s(i3)[p0] - UL_p.s(i3)[p0]), //
						(U0_p.egas()[p0] - UL_p.egas()[p0]), (U0_p.tau()[p0] - UL_p.tau()[p0]));
				char_decomp(
						CM, //
						(U0_p.rho()[p0] - UL_p.rho()[p0]) - norm * UL_p.rho()[p1], //
						(U0_p.s(i1)[p0] - UL_p.s(i1)[p0]) - norm * UL_p.s(i1)[p1], //
						(U0_p.s(i2)[p0] - UL_p.s(i2)[p0]) - norm * UL_p.s(i2)[p1], //
						(U0_p.s(i3)[p0] - UL_p.s(i3)[p0]) - norm * UL_p.s(i3)[p1], //
						(U0_p.egas()[p0] - UL_p.egas()[p0]) - norm * UL_p.egas()[p1],
						(U0_p.tau()[p0] - UL_p.tau()[p0]) - norm * UL_p.tau()[p1]);
				for (integer f = 0; f != NF; ++f) {
					const real clim1 = minmod(CL[f], CR[f]) / norm;
					const real clim2 = minmod(CP[f], CM[f]) / norm;
					C0[f] = std::max(C0[f], std::min(clim1, clim2));
					C0[f] = std::min(C0[f], std::max(clim1, clim2));
				}
				char_recomp(C0, U0_p.rho()[p1], U0_p.s(i1)[p1], U0_p.s(i2)[p1], U0_p.s(i3)[p1], U0_p.egas()[p1],
						U0_p.tau()[p1]);

			}
		}
	}

}

void grid::initialize(std::function<std::vector<real>(real, real, real)>&& f) {
	const auto quad_point = Fourier.quadrature_points();
	std::vector<simd_vector> U(NF, simd_vector(G3));
	for (integer i = 0; i != N3; ++i) {
		for (integer g = 0; g != G3; ++g) {
			const real x = cell_x[i] + quad_point[g][XDIM] * dx / real(2);
			const real y = cell_y[i] + quad_point[g][YDIM] * dx / real(2);
			const real z = cell_z[i] + quad_point[g][ZDIM] * dx / real(2);
			const auto this_u = f(x, y, z);
			for (integer f = 0; f != NF; ++f) {
				U[f][g] = this_u[f];
			}
		}
		for (integer f = 0; f != NF; ++f) {
			U_p[0][i][f] = Fourier.transform(U[f]);
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
	constexpr integer vertex_order[8] = { 0, 1, 3, 2, 4, 5, 7, 6 };
	std::set<node_point> node_list;

	const integer rk = 0;
	std::list<integer> zone_list;
	const auto& quad_weights = Fourier.quadrature_weights();
	std::vector<real> quad_points(P + 1);
	quad_points[0] = -dx / real(2);
	for (integer p = 0; p != P; ++p) {
		quad_points[p + 1] = quad_points[p] + quad_weights[p] * dx / real(2);
	}
	for (integer i = BW; i != NX - BW; ++i) {
		for (integer j = BW; j != NX - BW; ++j) {
			for (integer k = BW; k != NX - BW; ++k) {
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

	constexpr int nshapes = 1;
	int shapesize[1] = { NVERTEX };
	int shapetype[1] = { DB_ZONETYPE_HEX };
	int shapecnt[1] = { nzones };
	const char* coord_names[NDIM] = { "x", "y", "z" };

	DBfile *db = DBCreateReal(filename, DB_CLOBBER, DB_LOCAL, "Euler Mesh", DB_PDB);
	DBPutZonelist2(db, "zones", nzones, int(NDIM), zone_nodes.data(), nzones * NVERTEX, 0, 0, 0, shapetype, shapesize,
			shapecnt, nshapes, nullptr);
	DBPutUcdmesh(db, "mesh", int(NDIM), coord_names, node_coords.data(), nnodes, nzones, "zones", nullptr, DB_DOUBLE,
			nullptr);

	const char* field_names[] = { "rho", "egas", "tau", "sx", "sy", "sz" };
	constexpr integer I3 = INX * INX * INX;
	std::vector<double> rho(I3 * G3), sx(I3 * G3), sy(I3 * G3), sz(I3 * G3), egas(I3 * G3), tau(I3 * G3);
	std::array<double*, NF> u_data = { rho.data(), sx.data(), sy.data(), sz.data(), egas.data(), tau.data() };
	index = 0;
	for (integer i = BW; i != NX - BW; ++i) {
		for (integer j = BW; j != NX - BW; ++j) {
			for (integer k = BW; k != NX - BW; ++k) {
				const integer iii = NX * NX * i + NX * j + k;
				for (integer f = 0; f != NF; ++f) {
					auto U_h = Fourier.inverse_transform(U_p[rk][iii][f]);
					for (integer ggg = 0; ggg != G3; ++ggg) {
						u_data[f][index + ggg] = U_h[ggg];
					}
				}
				index += G3;
			}
		}
	}
	for (int f = 0; f != NF; ++f) {
		DBPutUcdvar1(db, field_names[f], "mesh", u_data[f], nzones, nullptr, 0, DB_DOUBLE, DB_ZONECENT, nullptr);
	}

	DBClose(db);
}



void grid::diagnostics(integer rk) {
	const auto& u_p = U_p[rk];
	printf("\n");
	for (integer f = 0; f != NF; ++f) {
		real sum = real(0);
		for (integer i = BW; i != NX - BW; ++i) {
			for (integer j = BW; j != NX / 2; ++j) {
				for (integer k = BW; k != NX - BW; ++k) {
					const real iii1 = i * NX * NX + j * NX + k;
					const real iii2 = i * NX * NX + (NX - 1 - j) * NX + k;
					for (integer l = 0; l != P; ++l) {
						for (integer m = 0; m != P - l; ++m) {
							for (integer n = 0; n != P - l - m; ++n) {
								integer ppp = Fourier.pindex(l, m, n);
								real v1 = std::abs(u_p[iii1][f][ppp]);
								real v2 = std::abs(u_p[iii2][f][ppp]);
								sum += std::pow(v1 - v2, 2);
							}
						}
					}
				}
			}
		}
		printf("- %16.8e ", std::sqrt(sum / real(INX * INX * INX)));
	}
	printf("\n");
}
