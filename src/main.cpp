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
#include <boost/chrono.hpp>
#include <mpi.h>

#define TMAX 100.0

int hpx_main(int argc, char* argv[]) {
#ifndef NDEBUG
	feenableexcept(FE_DIVBYZERO);
	feenableexcept(FE_INVALID);
	feenableexcept(FE_OVERFLOW);
#endif
	grid g;
	if( argc > 1 ) {
		g.initialize(argv[1]);
	} else {
		g.initialize(star);
	}
	while (g.step() < TMAX) {
	}
	return 0;
}

