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

#ifndef BREEZE2D_STATUS_H
#define BREEZE2D_STATUS_H

#include <assert.h>
#include <stdio.h>

#define breeze2d_set_error(status) { printf("Error %d\n", status); assert(0); }

#define BREEZE2D_UNDEFINED_POISSION_SOLVER_METHOD	10
#define BREEZE2D_FFT_PLAN_CREATION_FAILED		11

#endif // BREEZE2D_STATUS_H

