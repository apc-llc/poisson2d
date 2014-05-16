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

#ifndef SHUTTER_H
#define SHUTTER_H

#ifdef __cplusplus

#include <complex>

using namespace std;

#include "breeze2d.h"

#define complex complex<real>

extern "C"
{
#else

#include "breeze2d.h"

typedef real complex[2];

#endif // __cplusplus

void poisson2d_shutter_r(
	int m, int n, real hx, real hy,
	real* rhs, real* solution,
	real* alpha, real* beta, real* bc, real* ec);

void poisson2d_shutter_c(
	int m, int n, real hx, real hy,
	complex* rhs, complex* solution,
	real* alpha, complex* beta, complex* bc, complex* ec);

#ifdef __cplusplus
}
#endif // __cplusplus

#endif // SHUTTER_H

