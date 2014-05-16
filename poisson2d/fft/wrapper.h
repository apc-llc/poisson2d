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

#ifndef FFT_WRAPPER_H
#define FFT_WRAPPER_H

#if (defined(HAVE_FFTW) || defined(HAVE_FFTW_MKL))
#include <fftw3.h>
#endif
#ifdef HAVE_FFTW_MKL
#include <fftw3_mkl.h>
#endif

#include <breeze2d.h>

#if defined(HAVE_FFTW) || defined(HAVE_FFTW_MKL)
typedef fftw_r2r_kind fft_kind;
#endif

#ifdef HAVE_FFTW
#define FFT_ESTIMATE	FFTW_ESTIMATE
#define FFT_MEASURE	FFTW_MEASURE
#define FFT_EXHAUSTIVE	FFTW_EXHAUSTIVE
#define FFT_WISDOM_ONLY	FFTW_WISDOM_ONLY
#else
#define FFT_ESTIMATE	0
#define FFT_MEASURE	0
#define FFT_EXHAUSTIVE	0
#define FFT_WISDOM_ONLY	0
#endif

// Defines extended fft plan structure,
// incorporating settings specific to different
// fft libraries.
typedef struct
{
	real *in, *out;
	fft_kind kind;
	int n, nplans, howmany;
	int idist, odist;
#if defined(HAVE_FFTW) || defined(HAVE_FFTW_MKL)
#ifdef HAVE_SINGLE
	fftwf_plan *forward, *inverse;
#endif
#ifdef HAVE_DOUBLE
	fftw_plan *forward, *inverse;
#endif
#endif
}
fft_plan;

// Create fft processing plan.
fft_plan* fft_create(
	int n, real* in, real* out,
	fft_kind kind, unsigned flags);

// Create batched fft processing plan.
fft_plan* fft_create_multi(int n, int howmany,
	real* in, real* out, int idist, int odist,
	fft_kind kind, unsigned flags);

// Execute fft plan forward transform.
void fft_forward(fft_plan* plan);

// Execute fft plan inverse transform.
void fft_inverse(fft_plan* plan);

// Destroy the fft processing plan.
void fft_dispose(fft_plan* plan);

// Initialize threading support for further
// fft processing.
void fft_init_threads();

// Set the number of threads
// to use for further fft processing plans.
void fft_plan_with_nthreads(int nthreads);

// Malloc aligned data array of the specified size.
void* fft_malloc(size_t size);

// Release data array previously allocated with fft_malloc.
void fft_free(void* desc);

#endif // FFT_WRAPPER_H

