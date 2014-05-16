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
#include <math.h>
#include <stdio.h>
#include <string.h>

#define ABS(x) (((x) > 0) ? (x) : -(x))

// Note N is the number of INNER grid points,
// i.e. including boundaries the total number is N + 2
#define x0 -2.0 * M_PI
#define y0 -2.0 * M_PI
#define N 32
#define h (ABS(x0) * 2.0 / (N + 2))

real phi1[N * N];
real phi2[N * N];
real f[N][N];

real gbx[N], gex[N];
real gby[N], gey[N];

// Known exact solution (to check residual).
void solution(real phi[N][N])
{
	for (int j = 0; j < N; j++)
		for (int i = 0; i < N; i++)
		{
			real x = x0 + h * (i + 1);
			real y = y0 + h * (j + 1);
			phi[j][i] = sin(x) * sin(y);
		}
}

// Set right hand side f(x,y).
void init_f()
{
	for (int j = 0; j < N; j++)
		for (int i = 0; i < N; i++)
		{
			real x = x0 + h * (i + 1);
			real y = y0 + h * (j + 1);
			f[j][i] = 2.0 * sin(x) * sin(y);
		}
}

// Set boundary conditions g.
void init_g()
{
	real xN = x0 + h * (N + 1);
	real yN = y0 + h * (N + 1);

	for (int j = 0; j < N; j++)
	{
		real y = y0 + h * (j + 1);
		gbx[j] = sin(x0) * sin(y);
		gex[j] = sin(xN) * sin(y);
	}

	for (int i = 0; i < N; i++)
	{
		real x = x0 + h * (i + 1);
		gby[i] = sin(x) * sin(y0);
		gey[i] = sin(x) * sin(yN);
	}
}

#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define MAX(a,b) (((a) > (b)) ? (a) : (b))

void minmaxsum_diff(int n, real* A0, real* A1,
	real* min, real* max, real* sum)
{
	*min = *max = *sum = A0[0] - A1[0];
	for (int i = 1; i < n; i++)
	{
		printf("%f\t%f\n", A0[i], A1[i]);
	
		real val = A0[i] - A1[i];
		*min = MIN(*min, val);
		*max = MAX(*max, val);
		*sum += val;
	}
}

int main(int argc, char* argv[])
{
	printf("Solve 2D Poisson equation\n");
	printf("with Dirichlet boundary conditions:\n");
	printf("Lx = f in D, phi = g on dD.\n\n");
	printf("Method: finite differences\n\n");

	breeze2d_poisson_solver solver = breeze2d_poisson_solver_init(
		BREEZE2D_POISSON_SOLVER_FDIFFS,
		N, N, h, h, gbx, gex, gby, gey,
		(real*)f, phi1);

	init_f();
	init_g();
	
	breeze2d_poisson_solve(solver);
	breeze2d_poisson_solver_dispose(solver);

	solution(phi2);
	real min, max, sum;
	minmaxsum_diff(N * N, phi1, phi2, &min, &max, &sum);
	printf("residual min = %f, max = %f, sum = %f\n",
		min, max, sum);

	FILE* fp1 = fopen("poisson2d_fdiffs_phi1.bin", "w");
	fclose(fp1);
	breeze2d_dump2db("poisson2d_fdiffs_phi1.bin", N, N,
		phi1, N, N, sizeof(real));
	breeze2d_dump2db("poisson2d_fdiffs_phi2.bin", N, N,
		phi2, N, N, sizeof(real));

	FILE* fp2 = fopen("poisson2d_fdiffs_phi2.bin", "w");
	fclose(fp2);
	breeze2d_create_grads_ctl(N, N, "poisson2d_fdiffs_phi1");
	breeze2d_create_grads_gs(N, N, 1, "poisson2d_fdiffs_phi1");
	breeze2d_create_grads_pl("poisson2d_fdiffs_phi1");

	breeze2d_create_grads_ctl(N, N, "poisson2d_fdiffs_phi2");
	breeze2d_create_grads_gs(N, N, 1, "poisson2d_fdiffs_phi2");
	breeze2d_create_grads_pl("poisson2d_fdiffs_phi2");
	
	return 0;
}

