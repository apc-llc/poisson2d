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

#ifndef BREEZE2D_POISSON_H
#define BREEZE2D_POISSON_H

#ifndef BREEZE2D_H
#error Please always include <breeze2d.h>, and never include other BREEZE2D headers
#endif

/**
 * Defines identifier for Poisson solver
 * based on FFT.
 */
#define BREEZE2D_POISSON_SOLVER_FFT	0

/**
 * Defines identifier for Poisson solver
 * based on finite differences.
 */
#define BREEZE2D_POISSON_SOLVER_FDIFFS	1

/**
 * The 2D Poisson euqation solver descriptor.
 */
typedef void* breeze2d_poisson_solver;

/**
 * Initialize 2D Poisson equation solver for the
 * specified problem size and data arrays.
 * Note m and n are the numbers of INNER grid points,
 * i.e. including boundaries the total number is + 2.
 * @param mode - The Poisson equation solver to use
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
 */
breeze2d_poisson_solver breeze2d_poisson_solver_init(int mode,
	unsigned int m, unsigned int n, real hx, real hy,
	real* bx, real* ex, real* by, real* ey,
	real* rhs, real* solution);

/**
 * Release resources used by the specified solver instance.
 * @param desc - The solver configuration
 */
void breeze2d_poisson_solver_dispose(breeze2d_poisson_solver desc);

/**
 * Solve 2D Poisson equation with the given right hand side.
 * Place result to the output array specified in solver
 * configuration.
 * @param desc - The solver configuration
 */
void breeze2d_poisson_solve(breeze2d_poisson_solver desc);

/**
 * Count the number of floating-point operations
 * used by the specified solver configuration.
 * @param desc - The solver configuration
 * @return The number of FLOPs. 
 */
double breeze2d_poisson_flops(breeze2d_poisson_solver desc);

#endif // BREEZE2D_POISSON_H

