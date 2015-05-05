/*
 Copyright (C) 2011 by Rio Yokota, Simon Layton, Lorena Barba

 LMAXermission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.

 THE SOFTWARE IS LMAXROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXLMAXRESS OR
 IMLMAXLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A LMAXARTICULAR LMAXURLMAXOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COLMAXYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 THE SOFTWARE.
 */
#include "exafmm.hpp"

#define ODDEVEN(n) ((((n) & 1) == 1) ? -1 : 1)

real norm(const std::array<real, NDIM>& v) {
	return v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
}
//! Get r,theta,phi from x,y,z
void cart2sph(real& r, real& theta, real& phi, std::array<real, NDIM> dist) {
	r = std::sqrt(norm(dist));                           // r = sqrt(x^2 + y^2 + z^2)
	if (r == real(0)) {                                             // If r == 0
		theta = 0;                                                //  theta can be anything so we set it to 0
	} else {                                                    // If r != 0
		theta = std::acos(dist[2] / r);                                //  theta = acos(z / r)
	}                                                           // End if for r == 0
	if (std::abs(dist[0]) + std::abs(dist[1]) == real(0)) {                 // If |x| < eps & |y| < eps
		phi = 0;                                                  //  phi can be anything so we set it to 0
	} else {                                  // If x > 0
		phi = std::atan2(dist[1], dist[0]);                            //  phi = atan(y / x)
	}
}

//! Spherical to cartesian coordinates
void sph2cart(real r, real theta, real phi, std::array<real, NDIM>& spherical, std::array<real, NDIM> &cartesian) {
	cartesian[0] = sin(theta) * cos(phi) * spherical[0]         // x component (not x itself)
	+ cos(theta) * cos(phi) / r * spherical[1] - sin(phi) / r / sin(theta) * spherical[2];
	cartesian[1] = sin(theta) * sin(phi) * spherical[0]         // y component (not y itself)
	+ cos(theta) * sin(phi) / r * spherical[1] + cos(phi) / r / sin(theta) * spherical[2];
	cartesian[2] = cos(theta) * spherical[0]                    // z component (not z itself)
	- sin(theta) / r * spherical[1];
}

real *factorial;                                              //!< Factorial
real *prefactor;                                              //!< \f$ \sqrt{ \frac{(n - |m|)!}{(n + |m|)!} } \f$
real *Anm;                                                    //!< \f$ (-1)^n / \sqrt{ \frac{(n + m)!}{(n - m)!} } \f$
complex *Cnm;                                                 //!< M2L translation matrix \f$ C_{jn}^{km} \f$

exafmm::exafmm() {
	const complex I(0., 1.);                                     // Imaginary unit
	factorial = new real[LMAX];                                  // Factorial
	prefactor = new real[LMAX * LMAX];                                // sqrt( (n - |m|)! / (n + |m|)! )
	Anm = new real[LMAX * LMAX];                                // (-1)^n / sqrt( (n + m)! / (n - m)! )
	Cnm = new complex[LMAX * LMAX * LMAX * LMAX];                          // M2L translation matrix Cjknm

	factorial[0] = 1;                                           // Initialize factorial
	for (int n = 1; n != LMAX; ++n) {                                 // Loop to LMAX
		factorial[n] = factorial[n - 1] * n;                        //  n!
	}                                                           // End loop to LMAX

	for (int n = 0; n != LMAX; ++n) {                                 // Loop over n in Anm
		for (int m = -n; m <= n; ++m) {                              //  Loop over m in Anm
			int nm = n * n + n + m;                                       //   Index of Anm
			int nabsm = std::abs(m);                                     //   |m|
			real fnmm = real(1);                                        //   Initialize (n - m)!
			for (int i = 1; i <= n - m; ++i)
				fnmm *= i;                  //   (n - m)!
			real fnpm = real(1);                                        //   Initialize (n + m)!
			for (int i = 1; i <= n + m; ++i)
				fnpm *= i;                  //   (n + m)!
			real fnma = 1.0;                                        //   Initialize (n - |m|)!
			for (int i = 1; i <= n - nabsm; ++i)
				fnma *= i;              //   (n - |m|)!
			real fnpa = 1.0;                                        //   Initialize (n + |m|)!
			for (int i = 1; i <= n + nabsm; ++i)
				fnpa *= i;              //   (n + |m|)!
			prefactor[nm] = std::sqrt(fnma / fnpa);                   //   sqrt( (n - |m|)! / (n + |m|)! )
			Anm[nm] = ODDEVEN(n) / std::sqrt(fnmm * fnpm);              //   (-1)^n / sqrt( (n + m)! / (n - m)! )
		}                                                         //  End loop over m in Anm
	}                                                           // End loop over n in Anm

	for (int j = 0, jk = 0, jknm = 0; j != LMAX; ++j) {                   // Loop over j in Cjknm
		for (int k = -j; k <= j; ++k, ++jk) {                         //  Loop over k in Cjknm
			for (int n = 0, nm = 0; n != LMAX; ++n) {                       //   Loop over n in Cjknm
				for (int m = -n; m <= n; ++m, ++nm, ++jknm) {            //    Loop over m in Cjknm
					if (j + n < LMAX) {
						const int jnkm = (j + n) * (j + n) + j + n + m - k;               //     Index C_{j+n}^{m-k}
						Cnm[jknm] = std::pow(I, real(std::abs(k - m) - std::abs(k) - std::abs(m)))          //     Cjknm
						* real(ODDEVEN(j) * Anm[nm] * Anm[jk] / Anm[jnkm]);
					}
				}                                                     //    End loop over m in Cjknm
			}                                                       //   End loop over n in Cjknm
		}                                                         //  End loop over in k in Cjknm
	}                                                           // End loop over in j in Cjknm
}

//! Evaluate solid harmonics \f$ r^n Y_{n}^{m} \f$
void evalMultipole(real rho, real alpha, real beta, complex *Ynm, complex *YnmTheta) {
	const complex I(0., 1.);                                     // Imaginary unit
	real x = std::cos(alpha);                                   // x = cos(alpha)
	real y = std::sin(alpha);                                   // y = sin(alpha)
	real fact = 1;                                              // Initialize 2 * m + 1
	real pn = 1;                                                // Initialize Legendre polynomial LMAXn
	real rhom = 1;                                              // Initialize rho^m
	for (int m = 0; m != LMAX; ++m) {                                 // Loop over m in Ynm
		complex eim = std::exp(I * real(m * beta));               //  exp(i * m * beta)
		real p = pn;                                              //  Associated Legendre polynomial LMAXnm
		int npn = m * m + 2 * m;                                  //  Index of Ynm for m > 0
		int nmn = m * m;                                          //  Index of Ynm for m < 0
		Ynm[npn] = rhom * p * prefactor[npn] * eim;               //  rho^m * Ynm for m > 0
		Ynm[nmn] = std::conj(Ynm[npn]);                           //  Use conjugate relation for m < 0
		real p1 = p;                                              //  LMAXnm-1
		p = x * (2 * m + 1) * p1;                                 //  LMAXnm using recurrence relation
		if (y != real(0)) {
			YnmTheta[npn] = rhom * (p - (m + 1) * x * p1) / y * prefactor[npn] * eim;   // theta derivative of r^n * Ynm
		}
		rhom *= rho;                                              //  rho^m
		real rhon = rhom;                                         //  rho^n
		for (int n = m + 1; n != LMAX; ++n) {                             //  Loop over n in Ynm
			int npm = n * n + n + m;                                //   Index of Ynm for m > 0
			int nmm = n * n + n - m;                                //   Index of Ynm for m < 0
			Ynm[npm] = rhon * p * prefactor[npm] * eim;             //   rho^n * Ynm
			Ynm[nmm] = std::conj(Ynm[npm]);                         //   Use conjugate relation for m < 0
			real p2 = p1;                                           //   LMAXnm-2
			p1 = p;                                                 //   LMAXnm-1
			p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);           //   LMAXnm using recurrence relation
			if (y != real(0)) {
				YnmTheta[npm] = rhon * ((n - m + 1) * p - (n + 1) * x * p1) / y * prefactor[npm] * eim; // theta derivative
			}
			rhon *= rho;                                            //   Update rho^n
		}                                                         //  End loop over n in Ynm
		pn = -pn * fact * y;                                      //  LMAXn
		fact += 2;                                                //  2 * m + 1
	}
}

//! Evaluate singular harmonics \f$ r^{-n-1} Y_n^m \f$
void evalLocal(real rho, real alpha, real beta, complex *Ynm, complex *YnmTheta) {
	const complex I(0., 1.);                                     // Imaginary unit
	real x = std::cos(alpha);                                   // x = cos(alpha)
	real y = std::sin(alpha);                                   // y = sin(alpha)
	real fact = 1;                                              // Initialize 2 * m + 1
	real pn = 1;                                                // Initialize Legendre polynomial LMAXn
	real rhom = 1.0 / rho;                                      // Initialize rho^(-m-1)
	for (int m = 0; m != LMAX; ++m) {                                 // Loop over m in Ynm
		complex eim = std::exp(I * real(m * beta));               //  exp(i * m * beta)
		real p = pn;                                              //  Associated Legendre polynomial LMAXnm
		int npn = m * m + 2 * m;                                  //  Index of Ynm for m > 0
		int nmn = m * m;                                          //  Index of Ynm for m < 0
		Ynm[npn] = rhom * p * prefactor[npn] * eim;               //  rho^(-m-1) * Ynm for m > 0
		Ynm[nmn] = std::conj(Ynm[npn]);                           //  Use conjugate relation for m < 0
		real p1 = p;                                              //  LMAXnm-1
		p = x * (2 * m + 1) * p1;                                 //  LMAXnm using recurrence relation
		YnmTheta[npn] = rhom * (p - (m + 1) * x * p1) / y * prefactor[npn] * eim;       // theta derivative of r^n * Ynm
		rhom /= rho;                                              //  rho^(-m-1)
		real rhon = rhom;                                         //  rho^(-n-1)
		for (int n = m + 1; n != LMAX; ++n) {                             //  Loop over n in Ynm
			int npm = n * n + n + m;                                //   Index of Ynm for m > 0
			int nmm = n * n + n - m;                                //   Index of Ynm for m < 0
			Ynm[npm] = rhon * p * prefactor[npm] * eim;             //   rho^n * Ynm for m > 0
			Ynm[nmm] = std::conj(Ynm[npm]);                         //   Use conjugate relation for m < 0
			real p2 = p1;                                           //   LMAXnm-2
			p1 = p;                                                 //   LMAXnm-1
			p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);               //   LMAXnm using recurrence relation
			YnmTheta[npm] = rhon * ((n - m + 1) * p - (n + 1) * x * p1) / y * prefactor[npm] * eim;  // theta derivative
			rhon /= rho;                                            //   rho^(-n-1)
		}                                                         //  End loop over n in Ynm
		pn = -pn * fact * y;                                      //  LMAXn
		fact += 2;                                                //  2 * m + 1
	}                                                           // End loop over m in Ynm
}

std::valarray<complex> exafmm::M2M(const std::valarray<complex>& Mj, std::array<real, NDIM> dist) {
	std::valarray<complex> Mi(complex(real(0), real(0)), LMAX * (LMAX + 1) / 2);
	const complex I(0., 1.);
	complex Ynm[LMAX * LMAX], YnmTheta[LMAX * LMAX];
	real rho, alpha, beta;
	cart2sph(rho, alpha, beta, dist);
	evalMultipole(rho, alpha, -beta, Ynm, YnmTheta);
	for (int j = 0; j != LMAX; ++j) {
		for (int k = 0; k <= j; ++k) {
			int jk = j * j + j + k;
			int jks = j * (j + 1) / 2 + k;
			complex M = 0;
			for (int n = 0; n <= j; ++n) {
				for (int m = -n; m <= std::min(k - 1, n); ++m) {
					if (j - n >= k - m) {
						int jnkm = (j - n) * (j - n) + j - n + k - m;
						int jnkms = (j - n) * (j - n + 1) / 2 + k - m;
						int nm = n * n + n + m;
						M += Mj[jnkms] * std::pow(I, real(m - std::abs(m))) * Ynm[nm]
								* real(ODDEVEN(n) * Anm[nm] * Anm[jnkm] / Anm[jk]);
					}
				}
				for (int m = k; m <= n; ++m) {
					if (j - n >= m - k) {
						int jnkm = (j - n) * (j - n) + j - n + k - m;
						int jnkms = (j - n) * (j - n + 1) / 2 - k + m;
						int nm = n * n + n + m;
						M += std::conj(Mj[jnkms]) * Ynm[nm] * real(ODDEVEN(k + n + m) * Anm[nm] * Anm[jnkm] / Anm[jk]);
					}
				}
			}
			Mi[jks] += M;
		}
	}
	return Mi;
}

std::valarray<complex> exafmm::M2L(const std::valarray<complex>& Mj, std::array<real, NDIM> dist) {
	std::valarray<complex> Li(complex(0, 0), LMAX * LMAX);
	complex Ynm[LMAX * LMAX], YnmTheta[LMAX * LMAX];
	real rho, alpha, beta;
	cart2sph(rho, alpha, beta, dist);
	evalLocal(rho, alpha, beta, Ynm, YnmTheta);
	for (int j = 0; j != LMAX; ++j) {
		for (int k = 0; k <= j; ++k) {
			int jk = j * j + j + k;
			int jks = j * (j + 1) / 2 + k;
			complex L = 0;
			for (int n = 0; n != LMAX - j; ++n) {
				for (int m = -n; m < 0; ++m) {
					int nm = n * n + n + m;
					int nms = n * (n + 1) / 2 - m;
					int jknm = jk * LMAX * LMAX + nm;
					int jnkm = (j + n) * (j + n) + j + n + m - k;
					L += std::conj(Mj[nms]) * Cnm[jknm] * Ynm[jnkm];
				}
				for (int m = 0; m <= n; ++m) {
					int nm = n * n + n + m;
					int nms = n * (n + 1) / 2 + m;
					int jknm = jk * LMAX * LMAX + nm;
					int jnkm = (j + n) * (j + n) + j + n + m - k;
					L += Mj[nms] * Cnm[jknm] * Ynm[jnkm];
				}
			}
			Li[jks] += L;
		}
	}
	return Li;
}

std::valarray<complex> exafmm::L2L(const std::valarray<complex>& Lj, std::array<real, NDIM> dist) {
	std::valarray<complex> Li(complex(0, 0), LMAX * (LMAX + 1) / 2);
	const complex I(0., 1.);
	complex Ynm[LMAX * LMAX], YnmTheta[LMAX * LMAX];
	real rho, alpha, beta;
	cart2sph(rho, alpha, beta, dist);
	evalMultipole(rho, alpha, beta, Ynm, YnmTheta);
	for (int j = 0; j != LMAX; ++j) {
		for (int k = 0; k <= j; ++k) {
			int jk = j * j + j + k;
			int jks = j * (j + 1) / 2 + k;
			complex L = 0;
			for (int n = j; n != LMAX; ++n) {
				for (int m = j + k - n; m < 0; ++m) {
					int jnkm = (n - j) * (n - j) + n - j + m - k;
					int nm = n * n + n - m;
					int nms = n * (n + 1) / 2 - m;
					L += std::conj(Lj[nms]) * Ynm[jnkm] * real(ODDEVEN(k) * Anm[jnkm] * Anm[jk] / Anm[nm]);
				}
				for (int m = 0; m <= n; ++m) {
					if (n - j >= std::abs(m - k)) {
						int jnkm = (n - j) * (n - j) + n - j + m - k;
						int nm = n * n + n + m;
						int nms = n * (n + 1) / 2 + m;
						L += Lj[nms] * std::pow(I, real(m - k - std::abs(m - k))) * Ynm[jnkm] * Anm[jnkm] * Anm[jk]
								/ Anm[nm];
					}
				}
			}
			Li[jks] += L;
		}
	}
	return Li;
}

std::valarray<complex> exafmm::P2M(std::array<real, NDIM> dist) {
//	printf( "%e %e %e\n", dist[0], dist[1], dist[2]);
	std::valarray<complex> Mj(complex(0, 0), LMAX * (LMAX + 1) / 2);
	complex Ynm[LMAX * LMAX], YnmTheta[LMAX * LMAX];
	real rho, alpha, beta;
	cart2sph(rho, alpha, beta, dist);
	evalMultipole(rho, alpha, -beta, Ynm, YnmTheta);
	for (int n = 0; n != LMAX; ++n) {
		for (int m = 0; m <= n; ++m) {
			int nm = n * n + n + m;
			int nms = n * (n + 1) / 2 + m;
			Mj[nms] += Ynm[nm];
		}
	}
//	printf( "		%e %e %e %e\n",Ynm[0].real(), Ynm[1].real(), Ynm[2].real(), Ynm[3].real());
//	printf( "		%e %e %e %e\n",Ynm[0].imag(), Ynm[1].imag(), Ynm[2].imag(), Ynm[3].imag());
	return Mj;
}

real exafmm::L2P(const std::valarray<complex>& Li, std::array<real, NDIM> dist) {
	real trg = real(0);
	const complex I(0., 1.);                                       // Imaginary unit
	complex Ynm[LMAX * LMAX], YnmTheta[LMAX * LMAX];
	real r, theta, phi;
	cart2sph(r, theta, phi, dist);
	evalMultipole(r, theta, phi, Ynm, YnmTheta);
	for (int n = 0; n != LMAX; ++n) {
		int nm = n * n + n;
		int nms = n * (n + 1) / 2;
		trg += std::real(Li[nms] * Ynm[nm]);
		for (int m = 1; m <= n; ++m) {
			nm = n * n + n + m;
			nms = n * (n + 1) / 2 + m;
			trg += 2 * std::real(Li[nms] * Ynm[nm]);
		}
	}
	return trg;
}

