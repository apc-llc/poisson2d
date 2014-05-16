/*
 * poisson2d - solver for 2D Poisson problem
 *             with Dirichlet or Neumann boundary conditions
 *
 * Copyright (C) 2011 Dmitry Mikushin
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "wrapper.h"
#include "fft.h"
#include "shutter.h"

#include <assert.h>
#include <malloc.h>
#include <math.h>
#include <omp.h>
#include <stdlib.h>

// Defines internal structure for fft solver.
struct poisson2d_fft_solver_t
{
	unsigned int m, n;
	real hx, hy;
	
	// Arrays for problem right hand side and solution.
	// Internally also used to keep transformed data and
	// shutter result.
	real *rhs, *solution;
	
	// Arrays for transformed boundary conditions.
	real *cby, *cey;

	// Shutter coefficients.
	real *alpha, *beta;
	
	fft_plan *plan_main, *plan_bc, *plan_ec;
};

#define MAX(a, b) ((a) > (b) ? (a) : (b))

// Initialize 2D Poisson equation FFT solver for the
// specified problem size and data arrays.
poisson2d_fft_solver poisson2d_fft_solver_init(
	unsigned int m, unsigned int n, real hx, real hy,
	real* bx, real* ex, real* by, real* ey,
	real* rhs, real* solution)
{
	fft_init_threads();
	fft_plan_with_nthreads(omp_get_num_threads());

	// Create and populate solver configuration structure.
	struct poisson2d_fft_solver_t* solver =
		(struct poisson2d_fft_solver_t*)malloc(
			sizeof(struct poisson2d_fft_solver_t));
	solver->m = m; solver->n = n; solver->hx = hx; solver->hy = hy;
	solver->rhs = rhs; solver->solution = solution;

	// Create main transform pass plan and optionally
	// benchmark it to let FFT select algorithm with
	// optimal performance.
	solver->plan_main = fft_create_multi(m, n,
		rhs, solution, m, m, FFTW_RODFT00, FFT_MEASURE);
	if (!solver->plan_main)
	{
		breeze2d_set_error(BREEZE2D_FFT_PLAN_CREATION_FAILED);
		return NULL;
	}

	// Create arrays for shutter coefficients.
	solver->alpha = (real*)fft_malloc(sizeof(real) * n);
	solver->beta = (real*)fft_malloc(sizeof(complex) * n);
	
	// Allocate arrays to hold transformed boundary conditions.
	solver->cby = (real*)fft_malloc(m * sizeof(real));
	solver->cey = (real*)fft_malloc(m * sizeof(real));

	// Create plans to transform boundary conditions.
	solver->plan_bc = fft_create(m, by, solver->cby,
		FFTW_RODFT00, FFT_MEASURE);
	if (!solver->plan_bc)
	{
		breeze2d_set_error(BREEZE2D_FFT_PLAN_CREATION_FAILED);
		return NULL;
	}
	solver->plan_ec = fft_create(m, ey, solver->cey,
		FFTW_RODFT00, FFT_WISDOM_ONLY | FFT_MEASURE);
	if (!solver->plan_ec)
	{
		breeze2d_set_error(BREEZE2D_FFT_PLAN_CREATION_FAILED);
		return NULL;
	}
	
	return (poisson2d_fft_solver)solver;
}

// Release resources used by the specified fft solver instance.
void poisson2d_fft_solver_dispose(poisson2d_fft_solver desc)
{
	struct poisson2d_fft_solver_t* solver =
		(struct poisson2d_fft_solver_t*)desc;
	
	fft_free(solver->alpha); fft_free(solver->beta);
	fft_free(solver->cby); fft_free(solver->cey);
	
	fft_dispose(solver->plan_main);
	fft_dispose(solver->plan_bc);
	fft_dispose(solver->plan_ec);
	
	free(solver);
}

// Solve 2D Poisson equation with the given right hand
// side using 1D fast Fourier transform by X and shutter by Y.
// Place result to the specified output array.
// Note m and n are the numbers of INNER grid points,
// i.e. including boundaries the total number is + 2.
void poisson2d_fft_solve(poisson2d_fft_solver desc)
{
	struct poisson2d_fft_solver_t* solver =
		(struct poisson2d_fft_solver_t*)desc;
	
	int n = solver->n, m = solver->m;
	real hx = solver->hx, hy = solver->hy;

	// Compute coefficients for the right hand side.
	fft_forward(solver->plan_main);

	// Compute coefficients for boundary conditions.
	fft_forward(solver->plan_bc);
	fft_forward(solver->plan_ec);

	// Solve m/2 3-diagonal systems of n equations
	// using shutter method.	
	poisson2d_shutter_r(m, n, hx, hy,
		solver->solution, solver->rhs,
		solver->alpha, solver->beta,
		solver->cby, solver->cey);

	// Compute result using inverse transform
	// on 3-diagonal systems solutions.
	fft_inverse(solver->plan_main);
}
