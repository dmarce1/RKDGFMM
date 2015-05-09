/*
 * main.cpp
 *
 *  Created on: Apr 11, 2015
 *      Author: dmarce1
 */

#include "grid.hpp"
#include "rk.hpp"
#include "tests/initial.hpp"
#include "tests/tests.hpp"
#include <fenv.h>

#define TMAX 10.0

int main(void) {
#ifndef NDEBUG
	feenableexcept(FE_DIVBYZERO);
	feenableexcept(FE_INVALID);
	feenableexcept(FE_OVERFLOW);
#endif
	grid g;
	g.initialize(star);
	g.enforce_boundaries(0);
	g.project(0);
	real t = real(0);
	real dt;
	real dx = real(2) / real(NX);
	g.enforce_boundaries(0);
	real amax = g.enforce_positivity(0);
	real output_dt = 2.0e-2;
	integer output_cnt = 0;
	g.compute_fmm(0);
	g.compute_analytic(sod_shock_tube_analytic, real(0));
	g.output("X.0.silo");
	++output_cnt;

	printf("t = 0.0 *\n");
	while (t < TMAX) {
		g.diagnostics(t);
		dt = dx / amax * cfl[NRK - 1];
		const real tremain = output_cnt * output_dt - t;
		dt = tremain / real(integer(tremain / dt) + 1);
		printf("t = %e dt = %e ", double(t + dt), double(dt));
		for (integer rk = 0; rk != NRK; ++rk) {
			const integer rk2 = rk == NRK - 1 ? 0 : rk + 1;
			g.compute_flux(rk);
			g.compute_du(rk);
			g.compute_next_u(rk, dt);
			g.enforce_boundaries(rk2);
			g.project(rk2);
			amax = g.enforce_positivity(rk2);
			g.compute_fmm(rk2);
		}
		t += dt;
		if (output_cnt <= integer(t / output_dt + 1.0e-6)) {
			char* tmp_ptr;
			if (!asprintf(&tmp_ptr, "X.%i.silo", int(output_cnt))) {
				abort();
			}
			printf("*");
			g.compute_analytic(sod_shock_tube_analytic, t);
			g.output(tmp_ptr);
			++output_cnt;
			free(tmp_ptr);
		}
		printf("\n");
	}
	return 0;
}

