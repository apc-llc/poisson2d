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

#include <assert.h>
#include <breeze2d.h>
#include <stdio.h>

// Append the specified data array to the file in binary mode.
void breeze2d_dump2db(
	char* filename, int dst_nx, int dst_ny,
	void* src, int src_nx, int src_ny,
	int size)
{
	assert(dst_nx >= 0); assert(dst_ny >= 0);
	assert(src_nx >= dst_nx); assert(src_ny >= dst_ny);
	
	// Open file in binary mode.
	FILE* fp = fopen(filename, "a");
	assert(fp);
	
	// Write data to file.
	int offset_x = (src_nx - dst_nx) / 2;
	int offset_y = (src_ny - dst_ny) / 2;
	for (int j = 0; j < dst_ny; j++)
	{
		int count = fwrite((char*)src + 
			((offset_y + j) * src_nx + offset_x) * size,
			size, dst_nx, fp);
		assert(count == dst_nx);
	}
	
	fclose(fp);
}

