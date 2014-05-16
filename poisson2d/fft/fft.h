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

#ifndef FFT_H
#define FFT_H

#include <breeze2d.h>

/**
 * The 2D Poisson euqation FFT solver descriptor.
 */
typedef void* poisson2d_fft_solver;

/**
 * Initialize 2D Poisson equation FFT solver for the
 * specified problem size and data arrays.
 * Note m and n are the numbers of INNER grid points,
 * i.e. including boundaries the total number is + 2.
 * @param m - The problem X grid dimension, excluding boundaries
 * @param n - The problem Y grid dimension, excluding boundaries
 * @param hx - The problem X grid step
 * @param hy - The problem Y grid step
 * @param bx - The X left side boundary m x 1 array
 * @param ex - The X right side boundary m x 1 array
 * @param by - The Y lower side boundary n x 1 array
 * @param ey - The Y upper side boundary n x 1 array
 * @param rhs - The right hand side m x n array
 * @param solution - The problem solution m x n array
 * @return The solver configuration.
 */
poisson2d_fft_solver poisson2d_fft_solver_init(
	unsigned int m, unsigned int n, real hx, real hy,
	real* bx, real* ex, real* by, real* ey,
	real* rhs, real* solution);

/**
 * Release resources used by the specified fft solver instance.
 * @param desc - The solver configuration
 */
void poisson2d_fft_solver_dispose(poisson2d_fft_solver desc);

/**
 * Solve 2D Poisson equation with the given right hand
 * side using 1D fast Fourier transform by X and shutter by Y.
 * Place result to the output array specified in solver
 * configuration.
 * @param desc - The solver configuration
 */
void poisson2d_fft_solve(poisson2d_fft_solver desc);

#endif // FFT_H

