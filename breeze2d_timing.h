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

#ifndef BREEZE2D_TIMING_H
#define BREEZE2D_TIMING_H

#ifndef BREEZE2D_H
#error Please always include <breeze2d.h>, and never include other BREEZE2D headers
#endif

#include <time.h>

#ifdef CLOCK_REALTIME
#define CLOCKID CLOCK_REALTIME
//#define CLOCKID CLOCK_MONOTONIC
//#define CLOCKID CLOCK_PROCESS_CPUTIME_ID
//#define CLOCKID CLOCK_THREAD_CPUTIME_ID
#else
#define CLOCKID 0
#define CLOCK_GETTIME_NOT_IMPLEMENTED
struct timespec
{
	long int tv_sec;            /* Seconds.  */
	long int tv_nsec;           /* Nanoseconds.  */
};
#endif

// Get the built-in timer value.
void breeze2d_get_time(struct timespec* t);

// Get the built-in timer measured values difference.
double breeze2d_get_time_diff(
	struct timespec t1, struct timespec t2);

// Print the built-in timer measured values difference.
void breeze2d_print_time_diff(
	struct timespec t1, struct timespec t2);

#endif // BREEZE2D_TIMING_H

