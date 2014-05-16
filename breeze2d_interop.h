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

#ifndef BREEZE2D_INTEROP_H
#define BREEZE2D_INTEROP_H

#ifndef BREEZE2D_H
#error Please always include <breeze2d.h>, and never include other BREEZE2D headers
#endif

/**
 * Read configuration file and populate initial structures.
 * @param filename - The configuration filename
 */
void breeze2d_opencfg(const char* filename);

/**
 * Read configuration parameter with the specified name.
 * @param value - The value pointer (filled on exit)
 * @param name - The configuration parameter name
 */
void breeze2d_readpar(void* value, char* name);

/**
 * Read (src_nx, src_ny) 2D array from text file,
 * and place it in the center of (dst_nx, dst_ny) array.
 * @param dst - The destination data array
 * @param dst_nx - The destination data array X dimension
 * @param dst_ny - The destination data array Y dimension
 * @param src_nx - The source data file X dimension
 * @param src_ny - The source data file Y dimension
 * @param filename - The source data filename
 * @param src_skip - The number of data frames to skip
 * from the beginning of the file
 * @param format - The source data value text format
 * @param size - The source data value size
 */
void breeze2d_read2dt(
	void* dst, int dst_nx, int dst_ny,
	char* filename, int src_nx, int src_ny,
	int src_skip, char* format, int size);

// Append the specified data array to the file in binary mode.
void breeze2d_dump2db(
	char* filename, int dst_nx, int dst_ny,
	void* src, int src_nx, int src_ny,
	int size);

// Create GrADS visualization ctl file.
void breeze2d_create_grads_ctl(int nx, int ny, const char* name);

// Create GrADS visualization gs script.
void breeze2d_create_grads_gs(int nx, int ny, int nt, const char* name);

// Create GrADS visualization perl script.
void breeze2d_create_grads_pl(const char* name);

#endif // BREEZE2D_INTEROP_H

