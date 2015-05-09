/*
 * fourier_legendre.cpp
 *
 *  Created on: Apr 23, 2015
 *      Author: dmarce1
 */

#include "fourier_legendre.hpp"
#include "legendre.hpp"
#include <stdio.h>
#include <atomic>
#include "exafmm.hpp"

integer fourier_legendre::PHI = P+1;
integer fourier_legendre::PHI3 = (PHI + 2) * (PHI + 1) * PHI / 6;
integer fourier_legendre::N_Lobatto;
std::vector<std::vector<simd_vector>> fourier_legendre::P2_phi;
std::vector<std::vector<simd_vector>> fourier_legendre::P2_gx;
std::vector<std::vector<simd_vector>> fourier_legendre::P2_gy;
std::vector<std::vector<simd_vector>> fourier_legendre::P2_gz;
std::vector<real> fourier_legendre::qpt;
std::vector<real> fourier_legendre::hires_qwt;
std::vector<real> fourier_legendre::hires_qpt;
std::vector<real> fourier_legendre::qwt;
std::vector<real> fourier_legendre::lobatto_qpt;
std::vector<real> fourier_legendre::lobatto_qwt;
std::vector<std::array<real, NDIM>> fourier_legendre::qpt_3d;
std::vector<simd_vector> fourier_legendre::transform_coefficient;
std::vector<simd_vector> fourier_legendre::inverse_transform_coefficient;
std::vector<std::vector<simd_vector>> fourier_legendre::dinverse_transform_coefficient_dx;
std::vector<std::vector<simd_vector>> fourier_legendre::lobatto_inverse_transform_coefficient;
std::vector<std::vector<simd_vector>> fourier_legendre::volume_transform_coefficient;
std::vector<std::vector<simd_vector>> fourier_legendre::volume_inverse_transform_coefficient;
std::vector<std::vector<simd_vector>> fourier_legendre::surface_transform_coefficient;
std::vector<std::vector<simd_vector>> fourier_legendre::surface_inverse_transform_coefficient;
std::vector<std::vector<simd_vector>> fourier_legendre::restrict_coefficients;
std::vector<std::vector<simd_vector>> fourier_legendre::prolong_coefficients;
std::vector<simd_vector> fourier_legendre::rho_2M;
std::vector<simd_vector> fourier_legendre::L2_phi;
std::vector<simd_vector> fourier_legendre::L2_gx;
std::vector<simd_vector> fourier_legendre::L2_gy;
std::vector<simd_vector> fourier_legendre::L2_gz;
std::vector<std::vector<simd_vector>> fourier_legendre::M2M;
std::vector<std::vector<simd_vector>> fourier_legendre::M2L;
std::vector<std::vector<simd_vector>> fourier_legendre::L2L;
simd_vector fourier_legendre::rpow_l;
simd_vector fourier_legendre::rpow_p;

integer fourier_legendre::pindex(integer l, integer m, integer n) {
	const integer lmn = l + m + n;
	const integer mn = n + m;
	return ((lmn + 2) * (lmn + 1) * lmn) / 6 + ((mn + 1) * mn) / 2 + n;

}

const std::vector<real>& fourier_legendre::quadrature_weights() {
	return qwt;
}

static integer gindex(integer x, integer y, integer z, integer N) {
	return x * N * N + y * N + z;
}

static integer gindex(integer y, integer z, integer N) {
	return y * N + z;
}

const std::vector<std::array<real, NDIM>>& fourier_legendre::quadrature_points() {
	return qpt_3d;
}

void fourier_legendre::allocate() {
	const integer ML_size = (2 * DNMAX + 1) * (2 * DNMAX + 1) * (2 * DNMAX + 1);
	rho_2M.resize(L2, simd_vector(real(0), P3));
	L2_phi.resize(G3, simd_vector(real(0), L2));
	L2_gx.resize(G3, simd_vector(real(0), L2));
	L2_gy.resize(G3, simd_vector(real(0), L2));
	L2_gz.resize(G3, simd_vector(real(0), L2));
	M2M.resize(NVERTEX, std::vector < simd_vector > (L2, simd_vector(real(0), L2)));
	L2L.resize(NVERTEX, std::vector < simd_vector > (L2, simd_vector(real(0), L2)));
	M2L.resize(ML_size, std::vector < simd_vector > (L2, simd_vector(real(0), L2)));
	P2_phi.resize((2 * DNMAX + 1) * (2 * DNMAX + 1) * (2 * DNMAX + 1),
			std::vector < simd_vector > (G3, simd_vector(real(0), P3)));
	P2_gx.resize((2 * DNMAX + 1) * (2 * DNMAX + 1) * (2 * DNMAX + 1),
			std::vector < simd_vector > (G3, simd_vector(real(0), P3)));
	P2_gy.resize((2 * DNMAX + 1) * (2 * DNMAX + 1) * (2 * DNMAX + 1),
			std::vector < simd_vector > (G3, simd_vector(real(0), P3)));
	P2_gz.resize((2 * DNMAX + 1) * (2 * DNMAX + 1) * (2 * DNMAX + 1),
			std::vector < simd_vector > (G3, simd_vector(real(0), P3)));
	qpt_3d.resize(P * P * P);
	qpt.resize(P);
	qwt.resize(P);
	hires_qwt.resize(PHI);
	hires_qpt.resize(PHI);
	restrict_coefficients.resize(P3, std::vector < simd_vector > (NVERTEX, simd_vector(P3)));
	prolong_coefficients.resize(P3, std::vector < simd_vector > (NVERTEX, simd_vector(P3)));
	lobatto_qpt.resize(N_Lobatto);
	lobatto_qwt.resize(N_Lobatto);
	transform_coefficient.resize(P3);
	surface_transform_coefficient.resize(NFACE, std::vector < simd_vector > (P3));
	volume_transform_coefficient.resize(NDIM, std::vector < simd_vector > (P3));
	for (integer i = 0; i != P3; ++i) {
		transform_coefficient[i] = simd_vector(real(0), P * P * P);
		for (integer fc = 0; fc != NFACE; ++fc) {
			surface_transform_coefficient[fc][i] = simd_vector(real(0), P * P);
		}
		for (integer dim = 0; dim != NDIM; ++dim) {
			volume_transform_coefficient[dim][i] = simd_vector(real(0), P * P * P);
		}
	}
	inverse_transform_coefficient.resize(P * P * P);
	lobatto_inverse_transform_coefficient.resize(NDIM, std::vector < simd_vector > (P * P * N_Lobatto));
	surface_inverse_transform_coefficient.resize(NFACE, std::vector < simd_vector > (P * P));
	volume_inverse_transform_coefficient.resize(NDIM, std::vector < simd_vector > (P * P * P));
	dinverse_transform_coefficient_dx.resize(NDIM, std::vector < simd_vector > (P * P * P));
	for (integer i = 0; i != P * P * N_Lobatto; ++i) {
		for (integer dim = 0; dim != NDIM; ++dim) {
			lobatto_inverse_transform_coefficient[dim][i] = simd_vector(real(0), P3);
		}
	}
	for (integer i = 0; i != P * P * P; ++i) {
		inverse_transform_coefficient[i] = simd_vector(real(0), P3);
		for (integer dim = 0; dim != NDIM; ++dim) {
			volume_inverse_transform_coefficient[dim][i] = simd_vector(real(0), P3);
			dinverse_transform_coefficient_dx[dim][i] = simd_vector(real(0), P3);
		}
	}
	for (integer i = 0; i != P * P; ++i) {
		for (integer fc = 0; fc != NFACE; ++fc) {
			surface_inverse_transform_coefficient[fc][i] = simd_vector(real(0), P3);
		}
	}

}

fourier_legendre::fourier_legendre() {

	static std::atomic<int> initialization_begun(0);
	static std::atomic<int> initialization_complete(0);

	if (initialization_begun++ != 0) {
		while (initialization_complete == 0) {
		}
		return;
	}

	N_Lobatto = std::max(integer((P + 4) / real(2)), integer(2));
	allocate();

	rpow_l.resize(L2);
	rpow_p.resize(P3);
	for (integer l = 0; l != LMAX; ++l) {
		for (integer m = -l; m <= l; ++m) {
			rpow_l[l * l + l + m] = real(l);
		}
	}
	for (integer l = 0; l != P; ++l) {
		for (integer m = 0; m != P - l; ++m) {
			for (integer n = 0; n != P - l - m; ++n) {
				rpow_p[pindex(l, m, n)] = real(l + m + n);
			}
		}
	}

	char* filename;
	if (!asprintf(&filename, "basis.%i.%i.dat", int(P), int(LMAX))) {
		abort();
	}
	FILE* test = fopen(filename, "rb");
	if (test != NULL) {
		printf("Found transform coefficient file\n");
		fclose(test);
		read(filename);
		free(filename);
	} else {

		qpt = LegendreP_roots(P);
		for (integer gx = 0; gx != P; ++gx) {
			for (integer gy = 0; gy != P; ++gy) {
				for (integer gz = 0; gz != P; ++gz) {
					const integer i = gindex(gx, gy, gz, P);
					qpt_3d[i][0] = qpt[gx];
					qpt_3d[i][1] = qpt[gy];
					qpt_3d[i][2] = qpt[gz];
				}
			}
		}
		hires_qpt = LegendreP_roots(PHI);
		auto tmp = dLegendreP_dx_roots(N_Lobatto - 1);
		lobatto_qpt[0] = -real(1);
		for (integer qi = 0; qi < N_Lobatto - 2; ++qi) {
			lobatto_qpt[qi + 1] = tmp[qi];
		}
		lobatto_qpt[N_Lobatto - 1] = +real(1);
		for (integer n = 0; n != P; ++n) {
			const real x = qpt[n];
			const real dP_dx = dLegendreP_dx(P, x);
			qwt[n] = real(2) / ((real(1) - x * x) * dP_dx * dP_dx);
		}
		for (integer n = 0; n != PHI; ++n) {
			const real x = hires_qpt[n];
			const real dP_dx = dLegendreP_dx(PHI, x);
			hires_qwt[n] = real(2) / ((real(1) - x * x) * dP_dx * dP_dx);
		}
		for (integer n = 0; n != N_Lobatto; ++n) {
			const real x = lobatto_qpt[n];
			const real Pn = LegendreP(N_Lobatto - 1, x);
			lobatto_qwt[n] = real(2) / (real(N_Lobatto * (N_Lobatto - 1)) * Pn * Pn);
		}

		printf("Gauss-Legendre Quadrature Points\n");
		for (integer i = 0; i != P; ++i) {
			printf("%24.16e %24.16e\n", qpt[i], qwt[i]);
		}

		printf("Lobatto-Legendre Quadrature Points\n");
		for (integer i = 0; i != N_Lobatto; ++i) {
			printf("%24.16e %24.16e\n", lobatto_qpt[i], lobatto_qwt[i]);
		}

		printf("Computing Legendre transform and inverse transform coefficients...\n");
		for (integer gx = 0; gx != P; ++gx) {
			for (integer gy = 0; gy != P; ++gy) {
				for (integer gz = 0; gz != P; ++gz) {
					const integer i = gindex(gx, gy, gz, P);
					const auto px = LegendreP(qpt[gx], P);
					const auto py = LegendreP(qpt[gy], P);
					const auto pz = LegendreP(qpt[gz], P);
					const auto dpx_dx = dLegendreP_dx(qpt[gx], P);
					const auto dpy_dy = dLegendreP_dx(qpt[gy], P);
					const auto dpz_dz = dLegendreP_dx(qpt[gz], P);
					for (integer l = 0; l != P; ++l) {
						for (integer m = 0; m != P - l; ++m) {
							for (integer n = 0; n != P - l - m; ++n) {
								const auto wx = qwt[gx] * px[l];
								const auto wy = qwt[gy] * py[m];
								const auto wz = qwt[gz] * pz[n];
								const auto dwx_dx = qwt[gx] * dpx_dx[l];
								const auto dwy_dy = qwt[gy] * dpy_dy[m];
								const auto dwz_dz = qwt[gz] * dpz_dz[n];
								const auto nx = LegendreP_norm(l);
								const auto ny = LegendreP_norm(m);
								const auto nz = LegendreP_norm(n);
								const integer p = pindex(l, m, n);
								inverse_transform_coefficient[i][p] += px[l] * py[m] * pz[n];
								transform_coefficient[p][i] += (wx * wy * wz) / (nx * ny * nz);
								for (integer dim = 0; dim != NDIM; ++dim) {
									volume_inverse_transform_coefficient[dim][i][p] += px[l] * py[m] * pz[n];
								}
								dinverse_transform_coefficient_dx[XDIM][i][p] += dpx_dx[l] * py[m] * pz[n];
								dinverse_transform_coefficient_dx[YDIM][i][p] += px[l] * dpy_dy[m] * pz[n];
								dinverse_transform_coefficient_dx[ZDIM][i][p] += px[l] * py[m] * dpz_dz[n];
								volume_transform_coefficient[XDIM][p][i] += (dwx_dx * wy * wz) / (nx * ny * nz);
								volume_transform_coefficient[YDIM][p][i] += (wx * dwy_dy * wz) / (nx * ny * nz);
								volume_transform_coefficient[ZDIM][p][i] += (wx * wy * dwz_dz) / (nx * ny * nz);
							}
						}
					}
				}
			}
		}

		for (integer g1 = 0; g1 != P; ++g1) {
			for (integer g2 = 0; g2 != P; ++g2) {
				const integer i = gindex(g1, g2, P);
				const auto p1 = LegendreP(qpt[g1], P);
				const auto p2 = LegendreP(qpt[g2], P);
				for (integer l1 = 0; l1 != P; ++l1) {
					for (integer l2 = 0; l2 != P - l1; ++l2) {
						for (integer n = 0; n != P - l1 - l2; ++n) {
							const auto w1 = qwt[g1] * p1[l1];
							const auto w2 = qwt[g2] * p2[l2];
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

		const integer NL = N_Lobatto;
		const integer NG = P;
		printf("Computing Lobatto inverse transforms...\n");
		for (integer g1 = 0; g1 != NL; ++g1) {
			for (integer g2 = 0; g2 != NG; ++g2) {
				for (integer g3 = 0; g3 != NG; ++g3) {
					const integer i1 = (g1 * NG + g2) * NG + g3;
					const integer i2 = (g2 * NL + g1) * NG + g3;
					const integer i3 = (g2 * NG + g3) * NL + g1;
					const auto P1 = LegendreP(lobatto_qpt[g1], P);
					const auto P2 = LegendreP(qpt[g2], P);
					const auto P3 = LegendreP(qpt[g3], P);
					for (integer l1 = 0; l1 != P; ++l1) {
						for (integer l2 = 0; l2 != P - l1; ++l2) {
							for (integer l3 = 0; l3 != P - l1 - l2; ++l3) {
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

		printf("Computing legendre restrict and prolong coefficients...\n");
		for (integer ci = 0; ci != NVERTEX; ci++) {
			const integer xi = (ci >> 0) & 1;
			const integer yi = (ci >> 1) & 1;
			const integer zi = (ci >> 2) & 1;
			for (integer gx = 0; gx != PHI; ++gx) {
				for (integer gy = 0; gy != PHI; ++gy) {
					for (integer gz = 0; gz != PHI; ++gz) {
						const real xc = hires_qpt[gx];
						const real yc = hires_qpt[gy];
						const real zc = hires_qpt[gz];
						const real xp = xc / real(2) + real(xi) - real(1) / real(2);
						const real yp = yc / real(2) + real(yi) - real(1) / real(2);
						const real zp = zc / real(2) + real(zi) - real(1) / real(2);
						const auto pcx = LegendreP(xc, P);
						const auto pcy = LegendreP(yc, P);
						const auto pcz = LegendreP(zc, P);
						const auto ppx = LegendreP(xp, P);
						const auto ppy = LegendreP(yp, P);
						const auto ppz = LegendreP(zp, P);
						for (integer lc = 0; lc != P; ++lc) {
							for (integer mc = 0; mc != P - lc; ++mc) {
								for (integer nc = 0; nc != P - lc - mc; ++nc) {
									for (integer lp = 0; lp != P; ++lp) {
										for (integer mp = 0; mp != P - lp; ++mp) {
											for (integer np = 0; np != P - lp - mp; ++np) {
												const real wpx = hires_qwt[gx] / (real(NVERTEX) * LegendreP_norm(lp));
												const real wpy = hires_qwt[gy] / (real(NVERTEX) * LegendreP_norm(mp));
												const real wpz = hires_qwt[gz] / (real(NVERTEX) * LegendreP_norm(np));
												const real wcx = hires_qwt[gx] / LegendreP_norm(lc);
												const real wcy = hires_qwt[gy] / LegendreP_norm(mc);
												const real wcz = hires_qwt[gz] / LegendreP_norm(nc);
												const integer pp = pindex(lp, mp, np);
												const integer pc = pindex(lc, mc, nc);
												restrict_coefficients[pp][ci][pc] += pcx[lc] * pcy[mc] * pcz[nc]
														* ppx[lp] * ppy[mp] * ppz[np] * wpx * wpy * wpz;
												prolong_coefficients[pc][ci][pp] += pcx[lc] * pcy[mc] * pcz[nc]
														* ppx[lp] * ppy[mp] * ppz[np] * wcx * wcy * wcz;
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}

		exafmm fmm;

		printf("Computing Legendre-to-Spherical transform coefficients...\n");
		for (integer ppp = 0; ppp != P3; ++ppp) {
			std::valarray<real> rho_p(real(0), P3);
			rho_p[ppp] = real(1);
			auto rho_h = inverse_transform(rho_p);
			integer g = 0;
			for (integer g1 = 0; g1 != P; ++g1) {
				for (integer g2 = 0; g2 != P; ++g2) {
					for (integer g3 = 0; g3 != P; ++g3) {
						const real wt = qwt[g1] * qwt[g2] * qwt[g3];
						auto M = fmm.P2M( { qpt_3d[g][0], qpt_3d[g][1], qpt_3d[g][2] });
						for (integer l = 0; l != LMAX; ++l) {
							for (integer m = 0; m <= l; ++m) {
								rho_2M[l * l + l + m][ppp] += M[l * (l + 1) / 2 + m].real() * rho_h[g] * wt;
								rho_2M[l * l + l - m][ppp] += M[l * (l + 1) / 2 + m].imag() * rho_h[g] * wt;
							}
						}
						++g;
					}
				}
			}
		}
		printf("Computing Spherical-to-Legendre transform coefficients...\n");
		for (integer l = 0; l != LMAX; ++l) {
			for (integer m = -l; m <= +l; ++m) {
				std::valarray<complex> L(complex(0, 0), L2);
				if (m >= 0) {
					L[l * (l + 1) / 2 + m].real(1);
				} else {
					L[l * (l + 1) / 2 - m].imag(1);
				}
				std::valarray<real> phi_h(real(0), G3);
				std::valarray<real> gx_h(real(0), G3);
				std::valarray<real> gy_h(real(0), G3);
				std::valarray<real> gz_h(real(0), G3);
				for (integer g = 0; g != G3; ++g) {
					auto this_L = fmm.L2L(L, qpt_3d[g]);
					phi_h[g] = this_L[0].real();
					gx_h[g] = +this_L[2].real() * std::sqrt(2);
					gy_h[g] = -this_L[2].imag() * std::sqrt(2);
					gz_h[g] = -this_L[1].real();
				}
				const integer i = l * (l + 1) + m;
				const auto phi_p = transform(phi_h);
				for (integer g = 0; g != G3; ++g) {
					L2_phi[g][i] += phi_h[g];
					L2_gx[g][i] += gx_h[g];
					L2_gy[g][i] += gy_h[g];
					L2_gz[g][i] += gz_h[g];
				}
			}
		}

		printf("Computing spherical restrict and prolong coefficients...\n");
		for (integer ci = 0; ci != NVERTEX; ++ci) {
			std::array < real, NDIM> dist_m = {real(2 * ((ci >> 0) & 1) - 1), real(2 * ((ci >> 1) & 1) - 1), real(
						2 * ((ci >> 2) & 1) - 1)};
			std::array<real, NDIM> dist_l = {real(2 * ((ci >> 0) & 1) - 1), real(2 * ((ci >> 1) & 1) - 1), real(
						2 * ((ci >> 2) & 1) - 1)};
			for (integer d = 0; d != NDIM; ++d) {
				dist_m[d] *= -real(1);
				//			dist_l[d] /= +real(2);
			}
			for (integer l = 0; l != LMAX; ++l) {
				for (integer m = -l; m <= +l; ++m) {
					std::valarray<complex> A(complex(0, 0), LMAX * (LMAX + 1) / 2);
					const integer i_complex = l * (l + 1) / 2 + std::abs(m);
					const integer i_real = l * (l + 1) + m;
					if (m >= 0) {
						A[i_complex].real(1);
					} else {
						A[i_complex].imag(1);
					}
					const auto m2m = fmm.M2M(A, dist_m);
					const auto l2l = fmm.L2L(A, dist_l);
					for (integer j = 0; j != LMAX; ++j) {
						for (integer k = 0; k <= j; ++k) {
							const integer jk = j * (j + 1) / 2 + k;
							const integer jkp = j * j + j + k;
							const integer jkm = j * j + j - k;
							M2M[ci][jkp][i_real] += m2m[jk].real();
							L2L[ci][jkp][i_real] += l2l[jk].real();
							if (jkp != jkm) {
								M2M[ci][jkm][i_real] += m2m[jk].imag();
								L2L[ci][jkm][i_real] += l2l[jk].imag();
							}
						}
					}

				}
			}
		}

		printf("Computing spherical multipole interactions...\n");
		std::array < real, NDIM > dist;
		for (integer di = -DNMAX; di <= DNMAX; ++di) {
			dist[0] = real(2 * di);
			for (integer dj = -DNMAX; dj <= DNMAX; ++dj) {
				printf("	%i %i %i-%i\n", int(di), int(dj), int(-DNMAX), int(DNMAX));
				dist[1] = real(2 * dj);
				for (integer dk = -DNMAX; dk <= DNMAX; ++dk) {
					if (std::abs(di) + std::abs(dj) + std::abs(dk) != 0) {
						dist[2] = real(2 * dk);
						const integer index = ((di + DNMAX) * (2 * DNMAX + 1) + (dj + DNMAX)) * (2 * DNMAX + 1)
								+ (dk + DNMAX);
						for (integer l = 0; l != LMAX; ++l) {
							for (integer m = -l; m <= +l; ++m) {
								std::valarray<complex> A(complex(0, 0), LMAX * (LMAX + 1) / 2);
								const integer i_complex = l * (l + 1) / 2 + std::abs(m);
								const integer i_real = l * (l + 1) + m;
								if (m >= 0) {
									A[i_complex].real(1);
								} else {
									A[i_complex].imag(1);
								}
								const auto m2l = fmm.M2L(A, dist);
								for (integer j = 0; j != LMAX; ++j) {
									for (integer k = 0; k <= j; ++k) {
										const integer jk = j * (j + 1) / 2 + k;
										const integer jkp = j * j + j + k;
										const integer jkm = j * j + j - k;
										M2L[index][jkp][i_real] += m2l[jk].real();
										if (jkp != jkm) {
											M2L[index][jkm][i_real] += m2l[jk].imag();
										}
									}
								}
							}
						}
					}
				}
			}
		}

		printf("Computing gravity interaction coefficients...\n");
		for (integer di = -DNMAX; di <= DNMAX; ++di) {
			for (integer dj = -DNMAX; dj <= DNMAX; ++dj) {
				printf("	%i %i %i-%i\n", int(di), int(dj), int(-DNMAX), int(DNMAX));
				for (integer dk = -DNMAX; dk <= DNMAX; ++dk) {
					const integer iii = ((di + DNMAX) * (2 * DNMAX + 1) + (dj + DNMAX)) * (2 * DNMAX + 1)
							+ (dk + DNMAX);
					for (integer g1x = 0; g1x != P; ++g1x) {
						for (integer g1y = 0; g1y != P; ++g1y) {
							for (integer g1z = 0; g1z != P; ++g1z) {
								const integer g = gindex(g1x, g1y, g1z, P);
								for (integer g2x = 0; g2x != P; ++g2x) {
									for (integer g2y = 0; g2y != P; ++g2y) {
										for (integer g2z = 0; g2z != P; ++g2z) {
											const auto P2x = LegendreP(qpt[g2x], P);
											const auto P2y = LegendreP(qpt[g2y], P);
											const auto P2z = LegendreP(qpt[g2z], P);
											const real w2x = qwt[g2x];
											const real w2y = qwt[g2y];
											const real w2z = qwt[g2z];
											const real dx = real(2) * real(di) + qpt[g1x] - qpt[g2x];
											const real dy = real(2) * real(dj) + qpt[g1y] - qpt[g2y];
											const real dz = real(2) * real(dk) + qpt[g1z] - qpt[g2z];
											const real w2 = (w2x * w2y * w2z);
											const real r = std::sqrt(dx * dx + dy * dy + dz * dz);
											if (!(di == 0 && dj == 0 && dk == 0 && g1x == g2x && g1y == g2y
													&& g1z == g2z)) {
												const real kernel_phi = real(1) / r;
												const real kernel_gx = real(dx) / (r * r * r);
												const real kernel_gy = real(dy) / (r * r * r);
												const real kernel_gz = real(dz) / (r * r * r);
												for (integer l2 = 0; l2 != P; ++l2) {
													for (integer m2 = 0; m2 != P - l2; ++m2) {
														for (integer n2 = 0; n2 != P - l2 - m2; ++n2) {
															const integer p = pindex(l2, m2, n2);
															const real rho = P2x[l2] * P2y[m2] * P2z[n2];
															P2_phi[iii][g][p] += rho * w2 * kernel_phi;
															P2_gx[iii][g][p] += rho * w2 * kernel_gx;
															P2_gy[iii][g][p] += rho * w2 * kernel_gy;
															P2_gz[iii][g][p] += rho * w2 * kernel_gz;
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}

		write(filename);
		free(filename);
	}
	++initialization_complete;
}

simd_vector fourier_legendre::prolong(const simd_vector& up_p, integer ci) {
	simd_vector up_c(P3);
	for (integer ppp = 0; ppp != P3; ++ppp) {
		up_c[ppp] = (prolong_coefficients[ppp][ci] * up_p).sum();
	}
	return up_c;
}

simd_vector fourier_legendre::_restrict(const simd_vector& up_c, integer ci) {
	simd_vector up_p(P3);
	for (integer ppp = 0; ppp != P3; ++ppp) {
		up_p[ppp] = (prolong_coefficients[ppp][ci] * up_c).sum();
	}
	return up_p;
}

integer fourier_legendre::lobatto_point_count() {
	return P * P * N_Lobatto;
}

real fourier_legendre::lobatto_edge_weight() {
	return lobatto_qwt[0];
}

simd_vector fourier_legendre::transform(const simd_vector& U_h) {
	simd_vector U_p(P3);
	for (integer p = 0; p != P3; ++p) {
		U_p[p] = (transform_coefficient[p] * U_h).sum();
	}
	return U_p;
}

simd_vector fourier_legendre::inverse_transform(const simd_vector& U_p) {
	simd_vector U_h(P * P * P);
	for (integer i = 0; i != P * P * P; ++i) {
		U_h[i] = (inverse_transform_coefficient[i] * U_p).sum();
	}
	return U_h;
}

simd_vector fourier_legendre::lobatto_inverse_transform(const simd_vector& U_p, dimension dim) {
	simd_vector U_h(P * P * N_Lobatto);
	for (integer i = 0; i != P * P * N_Lobatto; ++i) {
		U_h[i] = (lobatto_inverse_transform_coefficient[dim][i] * U_p).sum();
	}
	return U_h;
}

simd_vector fourier_legendre::volume_transform(dimension dim, const simd_vector& U_h) {
	simd_vector U_p(P3);
	for (integer p = 0; p != P3; ++p) {
		U_p[p] = (volume_transform_coefficient[dim][p] * U_h).sum();
	}
	return U_p;
}

simd_vector fourier_legendre::dinverse_transform_dx(dimension dim, const simd_vector& U_p) {
	simd_vector U_h(P * P * P);
	for (integer i = 0; i != P * P * P; ++i) {
		U_h[i] = (dinverse_transform_coefficient_dx[dim][i] * U_p).sum();
	}
	return U_h;
}

simd_vector fourier_legendre::volume_inverse_transform(dimension dim, const simd_vector& U_p) {
	simd_vector U_h(P * P * P);
	for (integer i = 0; i != P * P * P; ++i) {
		U_h[i] = (volume_inverse_transform_coefficient[dim][i] * U_p).sum();
	}
	return U_h;
}

simd_vector fourier_legendre::surface_transform(face fc, const simd_vector& U_h) {
	simd_vector U_p(P3);
	for (integer p = 0; p != P3; ++p) {
		U_p[p] = (surface_transform_coefficient[fc][p] * U_h).sum();
	}
	return U_p;
}

simd_vector fourier_legendre::surface_inverse_transform(face fc, const simd_vector& U_p) {
	simd_vector U_h(P * P);
	for (integer i = 0; i != P * P; ++i) {
		U_h[i] = (surface_inverse_transform_coefficient[fc][i] * U_p).sum();
	}
	return U_h;
}

simd_vector fourier_legendre::p2_phi_transform(integer i, integer j, integer k, const simd_vector& rho_p, real dx) {
	simd_vector phi_h(real(0), G3);
	const real dx3 = dx * dx * dx;
	const integer index = (i + DNMAX) * (2 * DNMAX + 1) * (2 * DNMAX + 1) + (j + DNMAX) * (2 * DNMAX + 1) + (k + DNMAX);
	for (integer g = 0; g != G3; ++g) {
		phi_h[g] += ((P2_phi[index][g] * (dx3 / dx)) * rho_p).sum();
	}
	return phi_h;
}

simd_vector fourier_legendre::p2_gx_transform(integer i, integer j, integer k, const simd_vector& rho_p, real dx) {
	simd_vector phi_h(real(0), G3);
	const real dx3 = dx * dx * dx;
	const integer index = (i + DNMAX) * (2 * DNMAX + 1) * (2 * DNMAX + 1) + (j + DNMAX) * (2 * DNMAX + 1) + (k + DNMAX);
	for (integer g = 0; g != G3; ++g) {
		phi_h[g] += ((P2_gx[index][g] * (dx3 / dx)) * rho_p).sum();
	}
	return phi_h;
}

simd_vector fourier_legendre::p2_gy_transform(integer i, integer j, integer k, const simd_vector& rho_p, real dx) {
	simd_vector phi_h(real(0), G3);
	const real dx3 = dx * dx * dx;
	const integer index = (i + DNMAX) * (2 * DNMAX + 1) * (2 * DNMAX + 1) + (j + DNMAX) * (2 * DNMAX + 1) + (k + DNMAX);
	for (integer g = 0; g != G3; ++g) {
		phi_h[g] += ((P2_gy[index][g] * (dx3 / dx)) * rho_p).sum();
	}
	return phi_h;
}

simd_vector fourier_legendre::p2_gz_transform(integer i, integer j, integer k, const simd_vector& rho_p, real dx) {
	simd_vector phi_h(real(0), G3);
	const real dx3 = dx * dx * dx;
	const integer index = (i + DNMAX) * (2 * DNMAX + 1) * (2 * DNMAX + 1) + (j + DNMAX) * (2 * DNMAX + 1) + (k + DNMAX);
	for (integer g = 0; g != G3; ++g) {
		phi_h[g] += ((P2_gz[index][g] * (dx3 / dx)) * rho_p).sum();
	}
	return phi_h;
}

simd_vector fourier_legendre::l2p_transform(const simd_vector& l_in, real dx) {
	simd_vector phi_h(G3);
	simd_vector dx_np = std::pow(dx, +rpow_l);
	for (integer g = 0; g != G3; ++g) {
		phi_h[g] = (L2_phi[g] * (l_in * dx_np)).sum();
	}
	return phi_h;
}

simd_vector fourier_legendre::dl2p_transform_dx(dimension d, const simd_vector& l_in, real dx) {
	simd_vector g_out(G3);
	simd_vector dx_np = std::pow(dx, +rpow_l);
	switch (d) {
	case XDIM:
		for (integer g = 0; g != G3; ++g) {
			g_out[g] = (L2_gx[g] * (l_in * dx_np)).sum();
		}
		break;
	case YDIM:
		for (integer g = 0; g != G3; ++g) {
			g_out[g] = (L2_gy[g] * (l_in * dx_np)).sum();
		}
		break;
	case ZDIM:
		for (integer g = 0; g != G3; ++g) {
			g_out[g] = (L2_gz[g] * (l_in * dx_np)).sum();
		}
		break;
	}
	return g_out;
}

simd_vector fourier_legendre::p2m_transform(const simd_vector& rho_in, real dx) {
	simd_vector m_out(L2);
	const real dx3 = dx * dx * dx;
	simd_vector dx_np = std::pow(dx, +rpow_l);
	for (integer l = 0; l != L2; ++l) {
		m_out[l] = (rho_2M[l] * dx_np[l] * (rho_in * dx3)).sum();
	}
	return m_out;
}

simd_vector fourier_legendre::m2l_transform(integer i, integer j, integer k, const simd_vector& p_in, real dx) {
	simd_vector p_out(real(0), L2);
	const integer index = (i + DNMAX) * (2 * DNMAX + 1) * (2 * DNMAX + 1) + (j + DNMAX) * (2 * DNMAX + 1) + (k + DNMAX);
	simd_vector dx_nm = std::pow(dx, -rpow_l);
	simd_vector p2 = dx_nm * p_in;
	for (integer l = 0; l != L2; ++l) {
		p_out[l] += ((M2L[index][l] * dx_nm[l] / dx) * p2).sum();
	}
	return p_out;
}

simd_vector fourier_legendre::m2m_transform(integer ci, const simd_vector& m_in, real dx) {
	simd_vector m_out(real(0), L2);
	simd_vector dx_np = std::pow(dx, +rpow_l);
	simd_vector dx_nm = std::pow(dx, -rpow_l);
	simd_vector m2 = dx_nm * m_in;
	for (integer l = 0; l != L2; ++l) {
		m_out[l] += ((M2M[ci][l] * dx_np[l]) * m2).sum();
	}
	return m_out;
}

simd_vector fourier_legendre::l2l_transform(integer ci, const simd_vector& l_in, real dx) {
	simd_vector l_out(real(0), L2);
//	dx = dx / real(2);
	simd_vector dx_np = std::pow(dx, +rpow_l);
	simd_vector dx_nm = std::pow(dx, -rpow_l);
	simd_vector l2 = dx_np * l_in;
	for (integer l = 0; l != L2; ++l) {
		l_out[l] += ((L2L[ci][l] * dx_nm[l]) * l2).sum();
	}
	return l_out;
}

void fourier_legendre::write(const char* filename) {
	std::ofstream ofs(filename);
	boost::archive::binary_oarchive arc(ofs);
	arc << qpt;
	arc << hires_qwt;
	arc << hires_qpt;
	arc << qwt;
	arc << lobatto_qpt;
	arc << lobatto_qwt;
	arc << qpt_3d;
	arc << transform_coefficient;
	arc << inverse_transform_coefficient;
	arc << lobatto_inverse_transform_coefficient;
	arc << volume_transform_coefficient;
	arc << volume_inverse_transform_coefficient;
	arc << dinverse_transform_coefficient_dx;
	arc << surface_transform_coefficient;
	arc << surface_inverse_transform_coefficient;
	arc << restrict_coefficients;
	arc << prolong_coefficients;
	arc << P2_phi;
	arc << P2_gx;
	arc << P2_gy;
	arc << P2_gz;
	arc << rho_2M;
	arc << L2_phi;
	arc << L2_gx;
	arc << L2_gy;
	arc << L2_gz;
	arc << M2M;
	arc << M2L;
	arc << L2L;
	ofs.close();
}

void fourier_legendre::read(const char* filename) {
	std::ifstream ifs(filename);
	boost::archive::binary_iarchive arc(ifs);
	arc >> qpt;
	arc >> hires_qwt;
	arc >> hires_qpt;
	arc >> qwt;
	arc >> lobatto_qpt;
	arc >> lobatto_qwt;
	arc >> qpt_3d;
	arc >> transform_coefficient;
	arc >> inverse_transform_coefficient;
	arc >> lobatto_inverse_transform_coefficient;
	arc >> volume_transform_coefficient;
	arc >> volume_inverse_transform_coefficient;
	arc >> dinverse_transform_coefficient_dx;
	arc >> surface_transform_coefficient;
	arc >> surface_inverse_transform_coefficient;
	arc >> restrict_coefficients;
	arc >> prolong_coefficients;
	arc >> P2_phi;
	arc >> P2_gx;
	arc >> P2_gy;
	arc >> P2_gz;
	arc >> rho_2M;
	arc >> L2_phi;
	arc >> L2_gx;
	arc >> L2_gy;
	arc >> L2_gz;
	arc >> M2M;
	arc >> M2L;
	arc >> L2L;
	ifs.close();
}

