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

#include "shutter.h"

#include <math.h>

// Solve m/2 3-diagonal systems of n equations
// using shutter method in complex space.
extern "C" void poisson2d_shutter_c(
	int m, int n, real hx, real hy,
	complex* rhs, complex* solution,
	real* alpha, complex* beta, complex* bc, complex* ec)
{
	real r = hy / hx;
	real a = 1.0, c = 1.0;
	real invm = 0.5 / m;
	m /= 2;
	for (int p = 0; p < m; p++)
	{
		real val = r * sin(M_PI * (p + 1) * 2 * invm);
		real b = 2.0 + 4.0 * val * val;
		
		// ground b.c.
		{
			alpha[0] = 0.0;
			beta[0] = bc[p];
		}
		
		for (int k = 1; k < n; k++)
		{
			real val = 1.0 / (b - c * alpha[k - 1]);

			alpha[k] = a * val;
			beta[k] = (c * beta[k - 1] - 
				hy * hy * rhs[p + (k - 1) * m]) * val;
		}

		// top b.c.
		{
			solution[p + (n - 1) * m] = invm * ec[p];
		}
		
		for (int k = n - 1; k >= 1; k--)
		{
			solution[p + (k - 1) * m] =
				alpha[k] * solution[p + k * m] + invm * beta[k];
		}
	}
}

