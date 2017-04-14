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

#define MIN_ZN_FP   -10.0
#define MAX_ZN_FP    -8.0
#define STEP_ZN_FP    0.1
#define MIN_ZN_FPP    0.0
#define MAX_ZN_FPP    5.0
#define STEP_ZN_FPP   0.5
#define MIN_CU_FP    -4.0
#define MAX_CU_FP     0.0
#define STEP_CU_FP    0.2
#define MIN_CU_FPP    0.0
#define MAX_CU_FPP    5.0
#define STEP_CU_FPP   0.5

// Bespoke routine to find OSF and GooF values for a range of f' and f'' values for Zn and Cu sites in HBS
// For a given input model and hkl file
int main(int argc, char *argv[]) {
	unitcell *cell;
	hklF *F, *F2;
	int cr, i;
	FILE *ipf;
	FP osf = 2.0, GooF, ZnFP, ZnFPP, CuFP, CuFPP;

	if (argc < 3 || argc > 3) {
		fprintf(stderr, "Usage: canvas_hbs_f input.ins input.hkl\n");
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

        for (cr = 0; cr < F->nrefs; cr++) {
		F2->refs[cr].hkl[0] = F->refs[cr].hkl[0];
		F2->refs[cr].hkl[1] = F->refs[cr].hkl[1];
		F2->refs[cr].hkl[2] = F->refs[cr].hkl[2];
	}
	fprintf(stdout, "Done. Starting Canvas...\n");

	for (ZnFP = MIN_ZN_FP; ZnFP <= MAX_ZN_FP; ZnFP += STEP_ZN_FP)
	  for (ZnFPP = MIN_ZN_FPP; ZnFPP <= MAX_ZN_FPP; ZnFPP += STEP_ZN_FPP)
	    for (CuFP = MIN_CU_FP; CuFP <= MAX_CU_FP; CuFP += STEP_CU_FP)
	      for (CuFPP = MIN_CU_FPP; CuFPP <= MAX_CU_FPP; CuFPP += STEP_CU_FPP) {
		// Set new values on appropriate atoms
		for (i = 0; i < cell->natom; i++) {
			if (strncmp(cell->a[i].label, "ZN", 2) == 0) {
				cell->a[i].f[11] = ZnFP;
				cell->a[i].f[12] = ZnFPP;
			} else if (strncmp(cell->a[i].label, "CU", 2) == 0) {
				cell->a[i].f[11] = CuFP;
				cell->a[i].f[12] = CuFPP;
			}
		}
		// Generate F values
        	if (!hklF_fill(F2, cell)) { fprintf(stderr, "Error during F calulations!\n"); free(F2->refs); free(F2); free(F->refs); free(F); free(cell->sym); free(cell->a); free(cell); exit(10); }
		// Find OSF
		if (!find_osf(&osf, &GooF, F, F2, 0.1, 0.0, 1e-8)) { fprintf(stderr, "Error during finding OSF!\n"); free(F2->refs); free(F2); free(F->refs); free(F); free(cell->sym); free(cell->a); free(cell); exit(11); }
		// Output Line
		fprintf(stdout, "%8.4f,%8.4f,%8.4f,%8.4f,%.6E,%.6E\n", ZnFP, ZnFPP, CuFP, CuFPP, osf, GooF);
	      }

	free(F2->refs); free(F2); free(F->refs); free(F); free(cell->sym); free(cell->a); free(cell);
	return 0;
}
