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

#include <assert.h>
#include <malloc.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define ABS(x) (((x) > 0) ? (x) : -(x))

#define x0 0 /*0.5 * M_PI*/
#define y0 0

// Known exact solution (to check residual).
void solution(int m, int n, real hx, real hy, real* phi)
{
	for (int j = 0; j < n; j++)
		for (int i = 0; i < m; i++)
		{
			real x = x0 + hx * (i + 1);
			real y = y0 + hy * (j + 1);
			phi[j * m + i] = sin(x) * cos(y);
		}
}

// Set right hand side f(x,y).
void init_f(int m, int n, real hx, real hy, real* f)
{
	for (int j = 0; j < n; j++)
		for (int i = 0; i < m; i++)
		{
			real x = x0 + hx * (i + 1);
			real y = y0 + hy * (j + 1);
			f[j * m + i] = -2.0 * sin(x) * cos(y);
		}
}

// Set boundary conditions g.
void init_g(int m, int n, real hx, real hy,
	real* gbx, real* gex, real* gby, real* gey)
{
	real xN = x0 + hx * m;
	real yN = y0 + hy * n;

	for (int j = 0; j < n; j++)
	{
		real y = y0 + hy * (j + 1);
		gbx[j] = sin(x0) * cos(y);
		gex[j] = sin(xN) * cos(y);
	}

	for (int i = 0; i < m; i++)
	{
		real x = x0 + hx * (i + 1);
		gby[i] = sin(x) * cos(y0);
		gey[i] = sin(x) * cos(yN);
	}
}

#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define MAX(a,b) (((a) > (b)) ? (a) : (b))

void minmaxsum_diff(int n, real* A0, real* A1,
	real* min, real* max, real* sum)
{
	*min = *max = *sum = A0[0] - A1[0];
	if (n <= 32)
	{
		for (int i = 1; i < n; i++)
			printf("%f\t%f\t%f\n", A0[i], A1[i], A0[i] / A1[i]);
	}
	for (int i = 1; i < n; i++)
	{
		real val = A0[i] - A1[i];
		*min = MIN(*min, val);
		*max = MAX(*max, val);
		*sum += val;
	}
}

int main(int argc, char* argv[])
{
	printf("Solve 2D Poisson equation\n");
	printf("with Dirichlet b.c. by X,\n");
	printf("and Neumann b.c. by Y:\n");
	printf("Lx = f in D, phi = g on dD.\n\n");
	printf("Method: 1d fft + shutter\n\n");

#define USAGE() \
	{ \
		printf("Usage: %s <m> <n>, where\n", argv[0]); \
		printf("m, n - problem dimensions\n"); \
		printf("Note m and n denote the number of INNER grid points,\n"); \
		printf("i.e. including boundaries the total number is (m + 2) x (n + 2)\n"); \
		return 0; \
	}
	
	if (argc != 3) USAGE();

#ifdef HAVE_COACCEL
	char* ndevices = getenv("COACCEL_NUM_DEVICES");
	if (ndevices)
		printf("Solver ndevices = %d %s\n\n", atoi(ndevices),
			"(COACCEL)");
	else
	{
		printf("Solver ndevices = 1 %s\n\n",
			"(COACCEL)");
	}
#else
	char* nthreads = getenv("OMP_NUM_THREADS");
	if (nthreads)
		printf("Solver ncores = %d\n\n", atoi(nthreads));
	else
		printf("Solver ncores = 1\n\n");
#endif
	
	int m = atoi(argv[1]), n = atoi(argv[2]);
	if ((m <= 0) || (n <= 0)) USAGE();
	
	real hx = 2.0 * M_PI / (m + 1);
	real hy = 2.0 * M_PI / (n + 1);

	struct timespec start, finish;
	breeze2d_get_time(&start);

	real* phi1 = (real*)malloc(m * n * sizeof(real));
	real* phi2 = (real*)malloc(m * n * sizeof(real));
	real* f = (real*)malloc(m * n * sizeof(real));

	real* gbx = (real*)malloc(n * sizeof(real));
	real* gex = (real*)malloc(n * sizeof(real));
	real* gby = (real*)malloc(m * sizeof(real));
	real* gey = (real*)malloc(m * sizeof(real));

	breeze2d_poisson_solver solver = breeze2d_poisson_solver_init(
		BREEZE2D_POISSON_SOLVER_FFT,
		m, n, hx, hy, gbx, gex, gby, gey,
		(real*)f, phi1);

	init_f(m, n, hx, hy, f);
	init_g(m, n, hx, hy, gbx, gex, gby, gey);

	breeze2d_get_time(&finish);
	
	printf("Init time = %f\n", breeze2d_get_time_diff(
		start, finish));

	breeze2d_get_time(&start);

	breeze2d_poisson_solve(solver);

	breeze2d_get_time(&finish);
	
	double solver_time =
		breeze2d_get_time_diff(start, finish);
	printf("Solver time = %f\n", solver_time);
	
	printf("Solver gflops = %f\n",
		breeze2d_poisson_flops(solver) / solver_time * 1e-9);

	breeze2d_get_time(&start);

	solution(m, n, hx, hy, phi2);
	real min, max, sum;
	minmaxsum_diff(m * n, phi1, phi2, &min, &max, &sum);
	printf("residual min = %f, max = %f, sum = %f\n",
		min, max, sum);

	breeze2d_get_time(&finish);
	
	printf("Check time = %f\n", breeze2d_get_time_diff(
		start, finish));

	breeze2d_get_time(&start);

	FILE* fp1 = fopen("poisson2d_fft_phi1.bin", "w");
	fclose(fp1);
	breeze2d_dump2db("poisson2d_fft_phi1.bin", m, n,
		phi1, m, n, sizeof(real));
	breeze2d_create_grads_ctl(m, n, "poisson2d_fft_phi1");
	breeze2d_create_grads_gs(m, n, 1, "poisson2d_fft_phi1");
	breeze2d_create_grads_pl("poisson2d_fft_phi1");

	FILE* fp2 = fopen("poisson2d_fft_phi2.bin", "w");
	fclose(fp2);
	breeze2d_dump2db("poisson2d_fft_phi2.bin", m, n,
		phi2, m, n, sizeof(real));
	breeze2d_create_grads_ctl(m, n, "poisson2d_fft_phi2");
	breeze2d_create_grads_gs(m, n, 1, "poisson2d_fft_phi2");
	breeze2d_create_grads_pl("poisson2d_fft_phi2");

	breeze2d_get_time(&finish);
	
	printf("Output time = %f\n", breeze2d_get_time_diff(
		start, finish));

	breeze2d_get_time(&start);

	breeze2d_poisson_solver_dispose(solver);

	free(phi1); free(phi2); free(f);
	free(gbx); free(gex); free(gby); free(gey);

	breeze2d_get_time(&finish);

	printf("Deinit time = %f\n", breeze2d_get_time_diff(
		start, finish));
	
	return 0;
}

