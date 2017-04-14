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

// Determine an overall scale factor from an INS file and a given hklF2 data file
int main(int argc, char *argv[]) {
	unitcell *cell;
	hklF *F, *F2;
	int cr;
	FILE *ipf;
	FP osf, GooF;

	if (argc < 3 || argc > 3) {
		fprintf(stderr, "Usage: find_osf input.ins input.hkl\n");
		fprintf(stderr, "\n");
		exit(1);
	}

        cell = calloc(1, sizeof(unitcell));
        if (!cell) { fprintf(stderr, "Out of memory!\n"); exit(2); }

	ipf = fopen(argv[1], "r");
	if (!ipf) {
		free(cell);
		fprintf(stderr, "Error opening %s for input. Does it exist?\n", argv[1]);
		exit(3);
	}

	fprintf(stdout, "Reading %s...\n", argv[1]);
	if (!read_ins_cell(cell, ipf)) { fprintf(stderr, "Error during read!\n"); free(cell); fclose(ipf); exit(4); }
	fclose(ipf);

	print_cell(cell);
	fprintf(stdout, "Done. INS file successfully read. Reading reflection file...\n");

        F = calloc(1, sizeof(hklF));
        if (!F) { fprintf(stderr, "Out of memory!\n"); free(cell->sym); free(cell->a); free(cell); exit(5); }

	ipf = fopen(argv[2], "r");
	if (!ipf) {
		free(F); free(cell->sym); free(cell->a); free(cell);
		fprintf(stderr, "Error opening %s for input. Does it exist?\n", argv[2]);
		exit(6);
	}

	if (!read_hklF2(F, ipf)) { fprintf(stderr, "Error during read!\n"); free(F); free(cell->sym); free(cell->a); free(cell); fclose(ipf); exit(7); }
	fclose(ipf);

        F2 = calloc(1, sizeof(hklF));
        if (!F2) { fprintf(stderr, "Out of memory!\n"); free(F->refs); free(F); free(cell->sym); free(cell->a); free(cell); exit(8); }

        F2->nrefs = F->nrefs;
        F2->refs = calloc(F2->nrefs, sizeof(hklF_uno));
        if (!(F2->refs)) { fprintf(stderr, "Out of memory!\n"); free(cell->sym); free(cell->a); free(cell); free(F); free(F->refs); free(F2); exit(9); }
	cr = 0;
        for (cr = 0; cr < F->nrefs; cr++) {
		F2->refs[cr].hkl[0] = F->refs[cr].hkl[0];
		F2->refs[cr].hkl[1] = F->refs[cr].hkl[1];
		F2->refs[cr].hkl[2] = F->refs[cr].hkl[2];
	}

	fprintf(stdout, "Done. Generating F values from INS Model...\n");
        if (!hklF_fill(F2, cell)) { fprintf(stderr, "Error during F calulations!\n"); free(F2->refs); free(F2); free(F->refs); free(F); free(cell->sym); free(cell->a); free(cell); exit(10); }
	fprintf(stdout, "Done. Finding OSF...\n");

	osf = 2.0;
	if (!find_osf(&osf, &GooF, F, F2, 0.1, 0.0, 1e-5)) { fprintf(stderr, "Error during finding OSF!\n"); free(F2->refs); free(F2); free(F->refs); free(F); free(cell->sym); free(cell->a); free(cell); exit(11); }
	fprintf(stdout, "Done. Results:\n");
	fprintf(stdout, "  OSF = %.5E\n", osf);
	fprintf(stdout, " GooF = %.5E\n", GooF);

	free(F2->refs); free(F2); free(F->refs); free(F); free(cell->sym); free(cell->a); free(cell);
	return 0;
}
