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
#include <malloc.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h> // mkdir

// Create GrADS visualization ctl file.
void breeze2d_create_grads_ctl(int nx, int ny, const char* name)
{
	assert(nx > 0); assert(ny > 0);
	assert(name);

	char* filename = (char*)malloc(strlen(name) + 5);
	sprintf(filename, "%s.ctl", name);
	FILE* fp = fopen(filename, "w");

	fprintf(fp, "DSET %s.bin\n", name);
	fprintf(fp, "XDEF %d LINEAR 0 1\n", nx);
	fprintf(fp, "YDEF %d LINEAR 0 1\n", ny);
	fprintf(fp, "ZDEF 1 LINEAR 0 1\n");
	
	fprintf(fp, "EDEF 1 NAMES c\n");
	fprintf(fp, "TDEF 4001 LINEAR 18:00Z04jul2000 1hr\n");
	fprintf(fp, "UNDEF 1E+300\n");
	
	fprintf(fp, "VARS 1\n");
	fprintf(fp, "c 1 -1 source\n");
	fprintf(fp, "ENDVARS\n");
	
	fclose(fp);
	
	free(filename);
}

// Create GrADS visualization gs script.
void breeze2d_create_grads_gs(int nx, int ny, int nt, const char* name)
{
	assert(nx > 0); assert(ny > 0);
	assert(nt > 0);	assert(name);

	char* filename = (char*)malloc(strlen(name) + 5);
	sprintf(filename, "%s.gs", name);
	FILE* fp = fopen(filename, "w");

	fprintf(fp, "'open '\"%s.ctl\"\n", name);
	for (int it = 0; it < nt; it++)
	{
		fprintf(fp, "'set parea '0.5' '10.5' '1.5' '7.5\n");
		fprintf(fp, "'set xaxis '0' '%d' '%d\n", nx, nx / 16);
		fprintf(fp, "'set yaxis '0' '%d' '%d\n", ny, ny / 16);
		
		fprintf(fp, "'set mproj 'off\n");
		fprintf(fp, "'set mpdraw 'off\n");
		fprintf(fp, "'set gxout 'shaded\n");

		fprintf(fp, "'set clevs -'%f' -'%f' -'%f' -'%f' -'%f'",
			1.0, 0.8, 0.6, 0.4, 0.2);
		fprintf(fp, "'%f' '%f' '%f' '%f' '%f' '%f\n",
			0.0, 0.2, 0.4, 0.6, 0.8, 1.0);

		fprintf(fp, "'set t '%d\n", it + 1);
		fprintf(fp, "'display 'c\n");
		fprintf(fp, "'draw title '\"Source in %s flow - step %04d\"\n",
			name, it);
		fprintf(fp, "'run 'cbarm.gs\n");

		fprintf(fp, "'printim '\"%s_%04d.png\"' 'png' 'x800' 'y1024\n",
			name, it);
		fprintf(fp, "'clear 'graphics\n");
	}

	fclose(fp);

	free(filename);
}

// Create GrADS visualization perl script.
void breeze2d_create_grads_pl(const char* name)
{
	assert(name);

	char* filename = (char*)malloc(strlen(name) + 5);
	sprintf(filename, "%s.pl", name);
	FILE* fp = fopen(filename, "w");

	fprintf(fp, "#!/usr/bin/perl -w\n");
	fprintf(fp, "system(\"grads -lbxc \\\"run %s.gs\\\"\");\n", name);

	fclose(fp);

	free(filename);
}

