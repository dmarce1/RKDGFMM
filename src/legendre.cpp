/*
 * legendre.cpp
 *
 *  Created on: Apr 23, 2015
 *      Author: dmarce1
 */

#include "legendre.hpp"

std::vector<real> LegendreP(real x, integer N) {
	std::vector<real> P(N);
	P[0] = real(1);
	if (N > 1) {
		P[1] = x;
	}
	for (integer n = 1; n < N - 1; ++n) {
		P[n + 1] = (real(2 * n + 1) * x * P[n] - real(n) * P[n - 1]) / real(n + 1);
	}
	return P;
}

integer factorial(integer n) {
	integer rc;
	if (n != 0) {
		rc = integer(n) * factorial(n - 1);
	} else {
		rc = integer(1);
	}
	return rc;
}
;

std::vector<real> SphericalHarmonic(real x, real y, real z, integer N) {
	std::vector<real> Ynm(N * (N + 1) / 2);
	const real R = std::sqrt(x * x + y * y);
	const real r = std::sqrt(R * R + z * z);
	const real cos_theta = z / r;
	const real cos_phi = x / R;
	const real sin_phi = y / R;

	const auto Pnm = AssociatedLegendreP(cos_theta, N);
	for (integer l = 0; l != N; ++l) {
		for (integer m = 0; m != l; ++m) {
			const integer p_even = l * (l + 1) + m;
			const integer p_odd = l * (l + 1) - m;
			const integer lm = l * (l + 1) / 2 + m;
			const real norm = std::sqrt(real(factorial(l - m)) / real(factorial(l + m)));
			Ynm[p_even] = norm * Pnm[lm] * cos_phi;
			Ynm[p_odd] = norm * Pnm[lm] * sin_phi;
		}
	}
	return Ynm;
}

std::vector<real> AssociatedLegendreP(real x, integer N) {
	std::vector<real> Pnm(N * (N + 1) / 2);
	auto Pn = LegendreP(x, N);
	for (integer l = 0; l != N; ++l) {
		const integer p = l * (l + 1) / 2;
		Pnm[p] = Pn[l];
	}
	for (integer l = 1; l < N; ++l) {
		for (integer m = 0; m != l; ++m) {
			const integer p0 = l * (l + 1) / 2 + m;
			if (x != real(0)) {
				const integer pm = l * (l - 1) / 2 + m;
				const integer pp = l * (l + 1) / 2 + (m + 1);
				Pnm[pp] = (real(l - m) * x * Pnm[p0] - real(l + m) * Pnm[pm]) / std::sqrt(real(1) - x * x);
			} else {
				Pnm[p0] = real(0);
			}
		}
	}
	return Pnm;
}

real LegendreP(integer n, real x) {
	real Pnm1, Pn, Pnp1;
	Pnm1 = real(0);
	Pn = real(1);
	for (integer i = 0; i != n; ++i) {
		Pnp1 = (real(2 * i + 1) * x * Pn - real(i) * Pnm1) / real(i + 1);
		Pnm1 = Pn;
		Pn = Pnp1;
	}
	return Pn;
}

real LegendreP_norm(integer n) {
	return real(2) / real(2 * n + 1);
}

std::vector<real> dLegendreP_dx(real x, integer N) {
	std::vector<real> dP(N);
	const real xxm1 = (x * x - real(1));
	if (xxm1 != real(0)) {
		std::vector<real> P(N);
		P = LegendreP(x, N);
		for (integer n = 1; n != N; ++n) {
			dP[n] = (x * P[n] - P[n - 1]) * real(n) / xxm1;
		}
		dP[0] = real(0);
	} else {
		for (integer n = 0; n != N; ++n) {
			dP[n] = real(n * (n + 1) * std::pow(x, real(n + 1)) / 2);
		}
	}
	return dP;
}

real dLegendreP_dx(integer n, real x) {
	const real xxm1 = (x * x - real(1));
	real dPdx = real(0);
	if (xxm1 != real(0)) {
		real Pnm1, Pn, Pnp1;
		Pnm1 = real(0);
		Pn = real(1);
		for (integer i = 0; i != n; ++i) {
			Pnp1 = (real(2 * i + 1) * x * Pn - real(i) * Pnm1) / real(i + 1);
			Pnm1 = Pn;
			Pn = Pnp1;
		}
		dPdx = (x * Pn - Pnm1) * real(n) / xxm1;
	} else {
		dPdx = real(n * (n + 1) * std::pow(x, real(n + 1)) / 2);
	}
	return dPdx;
}

std::vector<real> LegendreP_roots(integer nn) {
	real all_roots[nn + 1][nn];
	all_roots[1][0] = real(0);
	for (integer n = 1; n != nn; ++n) {
		for (integer m = 0; m <= n; ++m) {
			real x;
			real this_min = (m == 0) ? -real(1) : all_roots[n][m - 1];
			real this_max = (m == n) ? +real(1) : all_roots[n][m];
			do {
				x = (this_min + this_max) / real(2);
				const real xval = LegendreP(n + 1, x);
				if (xval != real(0)) {
					const real sgn = xval * LegendreP(n + 1, this_min);
					if (sgn < real(0)) {
						this_max = x;
					} else {
						this_min = x;
					}
				} else {
					break;
				}
			} while (this_max - this_min > real(1.0e-14));
			all_roots[n + 1][m] = x;
		}

	}
	std::vector<real> rc(nn);
	for (integer n = 0; n != nn; ++n) {
		rc[n] = all_roots[nn][n];
	}
	return rc;
}

std::vector<real> dLegendreP_dx_roots(integer nn) {
	real all_roots[nn + 1][nn];
	all_roots[2][0] = real(0);
	for (integer n = 2; n < nn; ++n) {
		for (integer m = 0; m <= n - 1; ++m) {
			real x;
			real this_min = (m == 0) ? -real(1) : all_roots[n][m - 1];
			real this_max = (m == n - 1) ? +real(1) : all_roots[n][m];
			do {
				x = (this_min + this_max) / real(2);
				const real xval = dLegendreP_dx(n + 1, x);
				if (xval != real(0)) {
					const real sgn = xval * dLegendreP_dx(n + 1, this_min);
					if (sgn < real(0)) {
						this_max = x;
					} else {
						this_min = x;
					}
				} else {
					break;
				}
			} while (this_max - this_min > real(1.0e-14));
			all_roots[n + 1][m] = x;
		}

	}
	std::vector<real> rc(nn);
	for (integer n = 0; n != nn; ++n) {
		rc[n] = all_roots[nn][n];
	}
	return rc;
}
