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

#include <assert.h>
#include <malloc.h>
#include <math.h>
#include <string.h>

#ifdef HAVE_FFTW
// Filename to store FFTW wisdom.
#define FFTW_WISDOM_FILENAME ".wisdom"
static FILE* wisdom_file = NULL;
#endif

#ifdef HAVE_SINGLE
#define FFTW(call) fftwf_##call
#else
#define FFTW(call) fftw_##call
#endif 

// Create fft processing plan.
fft_plan* fft_create(
	int n, real* in, real* out,
	fft_kind kind, unsigned flags)
{
	// TODO: link existing plans **.
	fft_plan* plan = (fft_plan*)malloc(sizeof(fft_plan) +
		sizeof(plan->forward) + sizeof(plan->inverse));
#if defined(HAVE_FFTW) || defined(HAVE_FFTW_MKL)
	plan->forward = (FFTW(plan)*)(plan + 1);
	plan->inverse = plan->forward + 1;
#endif

	plan->in = in; plan->out = out;
#if defined(HAVE_FFTW) || defined(HAVE_FFTW_MKL)
#ifdef HAVE_FFTW
	if (!wisdom_file)
	{
		wisdom_file = fopen(FFTW_WISDOM_FILENAME, "r");
		if (wisdom_file)
		{
			FFTW(import_wisdom_from_file(wisdom_file));
			fclose(wisdom_file);
		}
	}
#endif
	plan->forward[0] = FFTW(plan_r2r_1d(n, in, out, FFTW_RODFT00, flags));
	if (!plan->forward)
	{
		free(plan);
		return NULL;
	}
	plan->inverse[0] = plan->forward[0];
#ifdef HAVE_FFTW
	if (!wisdom_file)
	{
		wisdom_file = fopen(FFTW_WISDOM_FILENAME, "w");
		FFTW(export_wisdom_to_file(wisdom_file));
		fclose(wisdom_file);
	}
#endif
#endif
	plan->n = n;
	plan->howmany = 1;
	plan->idist = 0;
	plan->odist = 0;
	plan->nplans = 1;
	plan->kind = kind;

	return plan;
}

// Create fft processing plan.
fft_plan* fft_create_multi(int n, int howmany,
	real* in, real* out, int idist, int odist,
	fft_kind kind, unsigned flags)
{
#ifdef HAVE_FFTW_MKL
	int nplans = howmany;
#else
	int nplans = 1;
#endif
	fft_plan* plan = (fft_plan*)malloc(sizeof(fft_plan) +
		nplans * (sizeof(plan->forward) + sizeof(plan->inverse)));
#if defined(HAVE_FFTW) || defined(HAVE_FFTW_MKL)
	plan->forward = (FFTW(plan)*)(plan + 1);
	plan->inverse = plan->forward + nplans;
#endif

	plan->in = in; plan->out = out;
#ifdef HAVE_FFTW
	int* nmany = (int*)malloc(sizeof(int) * howmany);
	fft_kind* kindmany = (FFTW(r2r_kind)*)malloc(sizeof(FFTW(r2r_kind)) * howmany);
	for (int i = 0; i < howmany; i++) { kindmany[i] = kind; nmany[i] = n; }
	plan->in = in; plan->out = out;
	plan->forward[0] = FFTW(plan_many_r2r(1, nmany, howmany,
		in, NULL, 1, idist, out, NULL, 1, odist,
		kindmany, flags));
	free(kindmany);
	free(nmany);
	if (!plan->forward[0])
	{
		free(plan);
		return NULL;
	}
	plan->inverse[0] = plan->forward[0];
#endif
#ifdef HAVE_FFTW_MKL
	for (int i = 0; i < howmany; i++) 
	{
		plan->forward[i] = FFTW(plan_r2r_1d(n,
			in + i * idist, out + i * odist, 
			kind, FFT_EXHAUSTIVE));
		plan->inverse[i] = plan->forward[i];
		if (!plan->forward[i])
		{
			for (int k = 0; k < i; k++)
				FFTW(destroy_plan(plan->forward[k]));
			free(plan);
			return NULL;
		}
	}
#endif
	plan->n = n;
	plan->howmany = howmany;
	plan->idist = idist;
	plan->odist = odist;
	plan->nplans = nplans;
	plan->kind = kind;

	return plan;
}

// Execute fft forward transform.
void fft_forward(fft_plan* plan)
{
#if defined(HAVE_FFTW) || defined(HAVE_FFTW_MKL)
#ifdef HAVE_FFTW_MKL
	#pragma omp parallel for
#endif
	for (int i = 0; i < plan->nplans; i++)
		FFTW(execute(plan->forward[i]));
#endif
}

// Execute fft forward transform.
void fft_inverse(fft_plan* plan)
{
#if defined(HAVE_FFTW) || defined(HAVE_FFTW_MKL)
#ifdef HAVE_FFTW_MKL
	#pragma omp parallel for
#endif
	for (int i = 0; i < plan->nplans; i++)
		FFTW(execute(plan->inverse[i]));
#endif
}

// Destroy the fft processing plan.
void fft_dispose(fft_plan* plan)
{
#if defined(HAVE_FFTW) || defined(HAVE_FFTW_MKL)
	for (int i = 0; i < plan->nplans; i++)
		FFTW(destroy_plan(plan->forward[i]));
#endif
	free(plan);
}

// Initialize threading support for further
// fft processing.
void fft_init_threads()
{
#if defined(HAVE_FFTW) || defined(HAVE_FFTW_MKL)
#if defined(HAVE_FFTW_THREADS)
	FFTW(init_threads());
#endif
#endif
}

// Set the number of devices (cpu threads or gpus)
// to use for further fft processing plans.
void fft_plan_with_nthreads(int nthreads)
{
#if defined(HAVE_FFTW) || defined(HAVE_FFTW_MKL)
#if defined(HAVE_FFTW_THREADS)
	FFTW(plan_with_nthreads(nthreads));
#endif
#endif
}

// Malloc aligned data array of the specified size.
void* fft_malloc(size_t size)
{
#if defined(HAVE_FFTW) || defined(HAVE_FFTW_MKL)
	return FFTW(malloc(size));
#endif
}

// Release data array previously allocated with fft_malloc.
void fft_free(void* desc)
{
#if defined(HAVE_FFTW) || defined(HAVE_FFTW_MKL)
	FFTW(free(desc));
#endif
}

// Count the number of floating-point operations
// used by the specified plan.
double fft_flops(fft_plan* plan)
{
	double base = plan->howmany * 
		plan->n * log(plan->n) / log(2.0);
#if defined(HAVE_FFTW) || defined(HAVE_FFTW_MKL)
	// TODO: update FLOPS for RODFT
	return 2.5 * base;
#endif
}

