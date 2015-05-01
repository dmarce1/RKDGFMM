/*
 * fourier_legendre.cpp
 *
 *  Created on: Apr 23, 2015
 *      Author: dmarce1
 */

#include "fourier_legendre.hpp"
#include "legendre.hpp"

#include <stdio.h>

integer fourier_legendre::pindex(integer l, integer m, integer n) const {
	const integer lmn = l + m + n;
	const integer mn = n + m;
	return ((lmn + 2) * (lmn + 1) * lmn) / 6 + ((mn + 1) * mn) / 2 + n;

}

const std::vector<real>& fourier_legendre::quadrature_weights() const {
	return quadrature_weight;
}

static integer gindex(integer x, integer y, integer z, integer N) {
	return x * N * N + y * N + z;
}

static integer gindex(integer y, integer z, integer N) {
	return y * N + z;
}

const std::vector<std::array<real, NDIM>>& fourier_legendre::quadrature_points() const {
	return quadrature_point_3d;
}

fourier_legendre::fourier_legendre(integer n) :
		N(n), N3((n + 2) * (n + 1) * n / 6) {
	N_Lobatto = std::max(integer((n + 4) / real(2)), integer(2));
	quadrature_point = LegendreP_roots(N);
	quadrature_weight.resize(N);
	auto tmp = dLegendreP_dx_roots(N_Lobatto - 1);
	lobatto_quadrature_point.resize(N_Lobatto);
	lobatto_quadrature_point[0] = -real(1);
	for (integer qi = 0; qi < N_Lobatto - 2; ++qi) {
		lobatto_quadrature_point[qi + 1] = tmp[qi];
	}
	lobatto_quadrature_point[N_Lobatto - 1] = +real(1);
	lobatto_quadrature_weight.resize(N_Lobatto);
	transform_coefficient.resize(N3);
	surface_transform_coefficient.resize(NFACE, std::vector<simd_vector>(N3));
	volume_transform_coefficient.resize(NDIM, std::vector<simd_vector>(N3));
	for (integer i = 0; i != N3; ++i) {
		transform_coefficient[i] = simd_vector(real(0), N * N * N);
		for (integer fc = 0; fc != NFACE; ++fc) {
			surface_transform_coefficient[fc][i] = simd_vector(real(0), N * N);
		}
		for (integer dim = 0; dim != NDIM; ++dim) {
			volume_transform_coefficient[dim][i] = simd_vector(real(0), N * N * N);
		}
	}
	inverse_transform_coefficient.resize(N * N * N);
	lobatto_inverse_transform_coefficient.resize(NDIM, std::vector<simd_vector>(N * N * N_Lobatto));
	surface_inverse_transform_coefficient.resize(NFACE, std::vector<simd_vector>(N * N));
	volume_inverse_transform_coefficient.resize(NDIM, std::vector<simd_vector>(N * N * N));
	for (integer i = 0; i != N * N * N_Lobatto; ++i) {
		for (integer dim = 0; dim != NDIM; ++dim) {
			lobatto_inverse_transform_coefficient[dim][i] = simd_vector(real(0), N3);
		}
	}
	for (integer i = 0; i != N * N * N; ++i) {
		inverse_transform_coefficient[i] = simd_vector(real(0), N3);
		for (integer fc = 0; fc != NFACE; ++fc) {
			for (integer dim = 0; dim != NDIM; ++dim) {
				volume_inverse_transform_coefficient[dim][i] = simd_vector(real(0), N3);
			}
		}
	}
	for (integer i = 0; i != N * N; ++i) {
		for (integer fc = 0; fc != NFACE; ++fc) {
			surface_inverse_transform_coefficient[fc][i] = simd_vector(real(0), N3);
		}
	}
	for (integer n = 0; n != N; ++n) {
		const real x = quadrature_point[n];
		const real dP_dx = dLegendreP_dx(N, x);
		quadrature_weight[n] = real(2) / ((real(1) - x * x) * dP_dx * dP_dx);
	}
	for (integer n = 0; n != N_Lobatto; ++n) {
		const real x = lobatto_quadrature_point[n];
		const real P = LegendreP(N_Lobatto - 1, x);
		lobatto_quadrature_weight[n] = real(2) / (real(N_Lobatto * (N_Lobatto - 1)) * P * P);
	}

	printf("Gauss-Legendre Quadrature Points\n");
	for (integer i = 0; i != N; ++i) {
		printf("%24.16e %24.16e\n", quadrature_point[i], quadrature_weight[i]);
	}

	printf("Lobatto-Legendre Quadrature Points\n");
	for (integer i = 0; i != N_Lobatto; ++i) {
		printf("%24.16e %24.16e\n", lobatto_quadrature_point[i], lobatto_quadrature_weight[i]);
	}

	for (integer gx = 0; gx != N; ++gx) {
		for (integer gy = 0; gy != N; ++gy) {
			for (integer gz = 0; gz != N; ++gz) {
				const integer i = gindex(gx, gy, gz, N);
				const auto px = LegendreP(quadrature_point[gx], N);
				const auto py = LegendreP(quadrature_point[gy], N);
				const auto pz = LegendreP(quadrature_point[gz], N);
				const auto dpx_dx = dLegendreP_dx(quadrature_point[gx], N);
				const auto dpy_dy = dLegendreP_dx(quadrature_point[gy], N);
				const auto dpz_dz = dLegendreP_dx(quadrature_point[gz], N);
				for (integer l = 0; l != N; ++l) {
					for (integer m = 0; m != N - l; ++m) {
						for (integer n = 0; n != N - l - m; ++n) {
							const auto wx = quadrature_weight[gx] * px[l];
							const auto wy = quadrature_weight[gy] * py[m];
							const auto wz = quadrature_weight[gz] * pz[n];
							const auto dwx_dx = quadrature_weight[gx] * dpx_dx[l];
							const auto dwy_dy = quadrature_weight[gy] * dpy_dy[m];
							const auto dwz_dz = quadrature_weight[gz] * dpz_dz[n];
							const auto nx = LegendreP_norm(l);
							const auto ny = LegendreP_norm(m);
							const auto nz = LegendreP_norm(n);
							const integer p = pindex(l, m, n);
							inverse_transform_coefficient[i][p] += px[l] * py[m] * pz[n];
							transform_coefficient[p][i] += (wx * wy * wz) / (nx * ny * nz);
							for (integer dim = 0; dim != NDIM; ++dim) {
								volume_inverse_transform_coefficient[dim][i][p] += px[l] * py[m] * pz[n];
							}
							volume_transform_coefficient[XDIM][p][i] += (dwx_dx * wy * wz) / (nx * ny * nz);
							volume_transform_coefficient[YDIM][p][i] += (wx * dwy_dy * wz) / (nx * ny * nz);
							volume_transform_coefficient[ZDIM][p][i] += (wx * wy * dwz_dz) / (nx * ny * nz);
						}
					}
				}
			}
		}
	}
	const integer NL = N_Lobatto;
	const integer NG = N;
	for (integer g1 = 0; g1 != NL; ++g1) {
		for (integer g2 = 0; g2 != NG; ++g2) {
			for (integer g3 = 0; g3 != NG; ++g3) {
				const integer i1 = (g1 * NG + g2) * NG + g3;
				const integer i2 = (g2 * NL + g1) * NG + g3;
				const integer i3 = (g2 * NG + g3) * NL + g1;
				const auto P1 = LegendreP(lobatto_quadrature_point[g1], N);
				const auto P2 = LegendreP(quadrature_point[g2], N);
				const auto P3 = LegendreP(quadrature_point[g3], N);
				for (integer l1 = 0; l1 != N; ++l1) {
					for (integer l2 = 0; l2 != N - l1; ++l2) {
						for (integer l3 = 0; l3 != N - l1 - l2; ++l3) {
							const integer p1 = pindex(l1, l2, l3);
							const integer p2 = pindex(l2, l1, l3);
							const integer p3 = pindex(l2, l3, l1);
							lobatto_inverse_transform_coefficient[XDIM][i1][p1] += P1[l1] * P2[l2] * P3[l3];
							lobatto_inverse_transform_coefficient[YDIM][i2][p2] += P1[l1] * P2[l2] * P3[l3];
							lobatto_inverse_transform_coefficient[ZDIM][i3][p3] += P1[l1] * P2[l2] * P3[l3];
						}
					}
				}
			}
		}
	}

	for (integer g1 = 0; g1 != N; ++g1) {
		for (integer g2 = 0; g2 != N; ++g2) {
			const integer i = gindex(g1, g2, N);
			const auto p1 = LegendreP(quadrature_point[g1], N);
			const auto p2 = LegendreP(quadrature_point[g2], N);
			for (integer l1 = 0; l1 != N; ++l1) {
				for (integer l2 = 0; l2 != N - l1; ++l2) {
					for (integer n = 0; n != N - l1 - l2; ++n) {
						const auto w1 = quadrature_weight[g1] * p1[l1];
						const auto w2 = quadrature_weight[g2] * p2[l2];
						const auto n1 = LegendreP_norm(l1);
						const auto n2 = LegendreP_norm(l2);
						const real nsgn = n % 2 == 0 ? +real(1) : -real(1);
						const real nnorm = LegendreP_norm(n);
						const integer pi = pindex(n, l1, l2);
						const integer pj = pindex(l1, n, l2);
						const integer pk = pindex(l1, l2, n);
						surface_inverse_transform_coefficient[XP][i][pi] += p1[l1] * p2[l2];
						surface_inverse_transform_coefficient[XM][i][pi] += p1[l1] * p2[l2] * nsgn;
						surface_inverse_transform_coefficient[YP][i][pj] += p1[l1] * p2[l2];
						surface_inverse_transform_coefficient[YM][i][pj] += p1[l1] * p2[l2] * nsgn;
						surface_inverse_transform_coefficient[ZP][i][pk] += p1[l1] * p2[l2];
						surface_inverse_transform_coefficient[ZM][i][pk] += p1[l1] * p2[l2] * nsgn;
						surface_transform_coefficient[XP][pi][i] += (w1 * w2) / (n1 * n2 * nnorm);
						surface_transform_coefficient[XM][pi][i] += (w1 * w2) / (n1 * n2 * nnorm) * nsgn;
						surface_transform_coefficient[YP][pj][i] += (w1 * w2) / (n1 * n2 * nnorm);
						surface_transform_coefficient[YM][pj][i] += (w1 * w2) / (n1 * n2 * nnorm) * nsgn;
						surface_transform_coefficient[ZP][pk][i] += (w1 * w2) / (n1 * n2 * nnorm);
						surface_transform_coefficient[ZM][pk][i] += (w1 * w2) / (n1 * n2 * nnorm) * nsgn;
					}
				}
			}
		}
	}
	/*for(integer dim = 0; dim != NFACE; ++dim ){
	for (integer gx = 0; gx != N; ++gx) {
		for (integer gy = 0; gy != N; ++gy) {
				printf( "\n");
				for( integer p = 0; p != 4; ++p)
					printf( "%16.7e ", surface_transform_coefficient[dim][p][gindex(gx,gy,N)]);
				printf( "\n");
		}
	}
	printf( "\n");
	}*/
	quadrature_point_3d.resize(N * N * N);
	for (integer gx = 0; gx != N; ++gx) {
		for (integer gy = 0; gy != N; ++gy) {
			for (integer gz = 0; gz != N; ++gz) {
				const integer i = gindex(gx, gy, gz, N);
				quadrature_point_3d[i][0] = quadrature_point[gx];
				quadrature_point_3d[i][1] = quadrature_point[gy];
				quadrature_point_3d[i][2] = quadrature_point[gz];
			}
		}
	}
}

integer fourier_legendre::lobatto_point_count() const {
	return N * N * N_Lobatto;
}

real fourier_legendre::lobatto_edge_weight() const {
	return lobatto_quadrature_weight[0];
}

simd_vector fourier_legendre::transform(const simd_vector& U_h) const {
	simd_vector U_p(N3);
	for (integer p = 0; p != N3; ++p) {
		U_p[p] = (transform_coefficient[p] * U_h).sum();
	}
	return U_p;
}

simd_vector fourier_legendre::inverse_transform(const simd_vector& U_p) const {
	simd_vector U_h(N * N * N);
	for (integer i = 0; i != N * N * N; ++i) {
		U_h[i] = (inverse_transform_coefficient[i] * U_p).sum();
	}
	return U_h;
}

simd_vector fourier_legendre::lobatto_inverse_transform(const simd_vector& U_p, dimension dim) const {
	simd_vector U_h(N * N * N_Lobatto);
	for (integer i = 0; i != N * N * N_Lobatto; ++i) {
		U_h[i] = (lobatto_inverse_transform_coefficient[dim][i] * U_p).sum();
	}
	return U_h;
}

simd_vector fourier_legendre::volume_transform(dimension dim, const simd_vector& U_h) const {
	simd_vector U_p(N3);
	for (integer p = 0; p != N3; ++p) {
		U_p[p] = (volume_transform_coefficient[dim][p] * U_h).sum();
	}
	return U_p;
}

simd_vector fourier_legendre::volume_inverse_transform(dimension dim, const simd_vector& U_p) const {
	simd_vector U_h(N * N * N);
	for (integer i = 0; i != N * N * N; ++i) {
		U_h[i] = (volume_inverse_transform_coefficient[dim][i] * U_p).sum();
	}
	return U_h;
}

simd_vector fourier_legendre::surface_transform(face fc, const simd_vector& U_h) const {
	simd_vector U_p(N3);
	for (integer p = 0; p != N3; ++p) {
		U_p[p] = (surface_transform_coefficient[fc][p] * U_h).sum();
	}
	return U_p;
}

simd_vector fourier_legendre::surface_inverse_transform(face fc, const simd_vector& U_p) const {
	simd_vector U_h(N * N);
	for (integer i = 0; i != N * N; ++i) {
		U_h[i] = (surface_inverse_transform_coefficient[fc][i] * U_p).sum();
	}
	return U_h;
}

