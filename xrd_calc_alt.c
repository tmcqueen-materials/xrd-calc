/*  This file is part of xrd-calc
 *
 *  xrd-calc is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  xrd-calc is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with xrd-calc.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Copyright (c) 2010 T. M. McQueen.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "defines.h"

// Calculate F values for the same hkls as those specified in the second
// input file.
int main(int argc, char *argv[]) {
	unitcell *cell;
	hklF *F;
	int cr;
	BOOL writeF = BFALSE;
	FILE *ipf;

	if (argc < 4 || argc > 5) {
		fprintf(stderr, "Usage: xrd_calc_alt input.ins input.hkl output.hkl <f-type>\n");
		fprintf(stderr, "      <f-type> is optional, but if present, 0 = write F2 and sigF2, 1 = write magF and phF.\n");
		fprintf(stderr, "\n");
		exit(1);
	}
	if (argc >= 5) { if (atoi(argv[4])) writeF = BTRUE; else writeF = BFALSE; }

        cell = calloc(1, sizeof(unitcell));
        if (!cell) { fprintf(stderr, "Out of memory!\n"); exit(2); }

	ipf = fopen(argv[1], "r");
	if (!ipf) {
		fprintf(stderr, "Error opening %s for input. Does it exist?\n", argv[1]);
		free(cell);
		exit(3);
	}

	fprintf(stdout, "Reading %s...\n", argv[1]);
	if (!read_ins_cell(cell, ipf)) { fprintf(stderr, "Error during read!\n"); free(cell); fclose(ipf); exit(4); }
	fclose(ipf);

	print_cell(cell);
	fprintf(stdout, "Done. INS file successfully read. Reading reflections list...\n");

        ipf = fopen(argv[2], "r");
        if (!ipf) {
                fprintf(stderr, "Error opening %s for input. Does it exist?\n", argv[2]);
		free(cell->sym); free(cell->a); free(cell);
                exit(5);
        }

        F = calloc(1, sizeof(hklF));
        if (!F) { fprintf(stderr, "Out of memory!\n"); free(cell->sym); free(cell->a); free(cell); exit(6); }
        if (!read_hklF2(F, ipf)) { fprintf(stderr, "Error during read!\n"); free(F); free(cell->sym); free(cell->a); free(cell); fclose(ipf); exit(7); }
        fclose(ipf);

	fprintf(stdout, "Done. Generating F values...\n");
        if (!hklF_fill(F, cell)) { fprintf(stderr, "Error during F calulations!\n"); free(F->refs); free(F); free(cell->sym); free(cell->a); free(cell); exit(6); }
	fprintf(stdout, "Done. Opening %s for output...\n", argv[3]);
	ipf = fopen(argv[3], "w");
	if (!ipf) { fprintf(stderr, "Error opening output file!\n"); free(F->refs); free(F); free(cell->sym); free(cell->a); free(cell); exit(6); }
	if (writeF) {
		fprintf(stdout, "Writing Reflections magF/phF...\n");
		print_hklF(F, ipf);
	} else {
		fprintf(stdout, "Writing Reflections...\n");
        	print_hklF2(F, ipf);
	}
	fclose(ipf);
	fprintf(stdout, "Success.\n");

	free(F->refs); free(F); free(cell->sym); free(cell->a); free(cell);
	return 0;
}
