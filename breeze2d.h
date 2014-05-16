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

#ifndef BREEZE2D_H
#define BREEZE2D_H

#ifdef HAVE_SINGLE
#define real float
#endif

#ifdef HAVE_DOUBLE
#define real double
#endif

#ifdef HAVE_INTEGER
#define real long
#endif

#if !defined(HAVE_SINGLE) && !defined(HAVE_DOUBLE)
#error "The single or double floating-precision must be set"
#endif

#if defined(HAVE_SINGLE) && defined(HAVE_DOUBLE)
#error "The single or double floating-precision must be set"
#endif

#if !defined(HAVE_FFTW) && !defined(HAVE_FFTW_MKL)
#error "The FFTW or MKL library must be turned on as FFT backend"
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif

#ifdef __cplusplus
extern "C"
{
#endif // __cplusplus

#include <breeze2d_poisson.h>
#include <breeze2d_interop.h>
#include <breeze2d_status.h>
#include <breeze2d_timing.h>

#ifdef __cplusplus
}
#endif // __cplusplus

#endif // BREEZE2D_H

