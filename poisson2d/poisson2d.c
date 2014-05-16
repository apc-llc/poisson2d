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

#include <breeze2d.h>
#include <malloc.h>

#include "fft/fft.h"

// Defines internal structure for solver.
struct breeze2d_poisson_solver_t
{
	int mode;
	void* desc; // nested solver descriptor
};

// Initialize 2D Poisson equation solver for the
// specified problem size and data arrays.
// Note m and n are the numbers of INNER grid points,
// i.e. including boundaries the total number is + 2.
breeze2d_poisson_solver breeze2d_poisson_solver_init(int mode,
	unsigned int m, unsigned int n, real hx, real hy,
	real* bx, real* ex, real* by, real* ey,
	real* rhs, real* solution)
{
	struct breeze2d_poisson_solver_t* solver =
		(struct breeze2d_poisson_solver_t*)malloc(
			sizeof(struct breeze2d_poisson_solver_t));
	solver->mode = mode;

	switch (solver->mode)
	{
	case BREEZE2D_POISSON_SOLVER_FFT :
		solver->desc = poisson2d_fft_solver_init(
			m, n, hx, hy, bx, ex, by, ey, rhs, solution);
		break;
	default :
		free(solver);
		breeze2d_set_error(BREEZE2D_UNDEFINED_POISSION_SOLVER_METHOD);
		return NULL;
	}
	
	return (breeze2d_poisson_solver)solver;
}

// Release resources used by the specified fft solver instance.
void breeze2d_poisson_solver_dispose(breeze2d_poisson_solver desc)
{
	struct breeze2d_poisson_solver_t* solver =
		(struct breeze2d_poisson_solver_t*)desc;

	switch (solver->mode)
	{
	case BREEZE2D_POISSON_SOLVER_FFT :
		poisson2d_fft_solver_dispose((poisson2d_fft_solver)solver->desc);
		break;
	default :
		breeze2d_set_error(BREEZE2D_UNDEFINED_POISSION_SOLVER_METHOD);
	}

	free(solver);
}

// Solve 2D Poisson equation with the given right hand side.
// Place result to the output array specified in solver
// configuration.
void breeze2d_poisson_solve(breeze2d_poisson_solver desc)
{
	struct breeze2d_poisson_solver_t* solver =
		(struct breeze2d_poisson_solver_t*)desc;

	switch (solver->mode)
	{
	case BREEZE2D_POISSON_SOLVER_FFT :
		poisson2d_fft_solve((poisson2d_fft_solver)solver->desc);
		break;
	default :
		breeze2d_set_error(BREEZE2D_UNDEFINED_POISSION_SOLVER_METHOD);
	}
}

// Count the number of floating-point operations
// used by the specified solver configuration.
double breeze2d_poisson_flops(breeze2d_poisson_solver desc)
{
	struct breeze2d_poisson_solver_t* solver =
		(struct breeze2d_poisson_solver_t*)desc;

	switch (solver->mode)
	{
	case BREEZE2D_POISSON_SOLVER_FFT :
		return poisson2d_fft_flops((poisson2d_fft_solver)solver->desc);
	default :
		breeze2d_set_error(BREEZE2D_UNDEFINED_POISSION_SOLVER_METHOD);
		return 0;
	}
}

