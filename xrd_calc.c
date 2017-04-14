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

// Calculate F values for a given range of h, k, and l and a given INS file
int main(int argc, char *argv[]) {
	unitcell *cell;
	hklF *F; int h, k, l;
	int cr;
	int minh, maxh, mink, maxk, minl, maxl;
	BOOL writeF = BFALSE;
	FILE *ipf;

	if (argc < 9 || argc > 10) {
		fprintf(stderr, "Usage: xrd_calc input.ins output.hkl min-h max-h min-k max-k min-l max-l <f-type>\n");
		fprintf(stderr, "      <f-type> is optional, but if present, 0 = write F2 and sigF2, 1 = write magF and phF.\n");
		fprintf(stderr, "\n");
		exit(1);
	}

	minh = atoi(argv[3]); maxh = atoi(argv[4]); if (maxh < minh) { h = maxh; maxh = minh; minh = h; }
	mink = atoi(argv[5]); maxk = atoi(argv[6]); if (maxk < mink) { k = maxk; maxk = mink; mink = k; }
	minl = atoi(argv[7]); maxl = atoi(argv[8]); if (maxl < minl) { l = maxl; maxl = minl; minl = l; }
	if (fabs(((double)minh-(double)maxh)*((double)mink-(double)maxk)*((double)minl-(double)maxl)) > (double)MAX_REFS) {
		fprintf(stderr, "Usage: xrd_calc input.ins output.hkl min-h max-h min-k max-k min-l max-l\n");
		fprintf(stderr, "UNUSUALLY LARGE number of reflections predicted (> %i). Were the h,k,l limits correctly parsed?:\n", MAX_REFS);
		fprintf(stderr, "  minh = %i, maxh = %i\n", minh, maxh);
                fprintf(stderr, "  mink = %i, maxk = %i\n", mink, maxk);
                fprintf(stderr, "  minl = %i, maxl = %i\n", minl, maxl);
		fprintf(stderr, "Ended.\n");
                exit(1);
	}
	if (argc >= 10) { if (atoi(argv[9])) writeF = BTRUE; else writeF = BFALSE; }

        cell = calloc(1, sizeof(unitcell));
        if (!cell) { fprintf(stderr, "Out of memory!\n"); exit(2); }

	ipf = fopen(argv[1], "r");
	if (!ipf) {
		fprintf(stderr, "Error opening %s for input. Does it exist?\n", argv[1]);
		exit(3);
	}

	fprintf(stdout, "Reading %s...\n", argv[1]);
	if (!read_ins_cell(cell, ipf)) { fprintf(stderr, "Error during read!\n"); free(cell); fclose(ipf); exit(4); }
	fclose(ipf);

	print_cell(cell);
	fprintf(stdout, "Done. INS file successfully read. Setting up reflection list...\n");

        F = calloc(1, sizeof(hklF));
        if (!F) { fprintf(stderr, "Out of memory!\n"); free(cell->sym); free(cell->a); free(cell); exit(5); }
        F->nrefs = (maxh-minh+1)*(maxk-mink+1)*(maxl-minl+1);
        F->refs = calloc(F->nrefs, sizeof(hklF_uno));
        if (!(F->refs)) { fprintf(stderr, "Out of memory!\n"); free(cell->sym); free(cell->a); free(cell); free(F); exit(5); }
	cr = 0;
        for (h = minh; h <= maxh; h++)
          for (k = mink; k <= maxk; k++)
            for (l = minl; l <= maxl; l++) {
                F->refs[cr].hkl[0] = h;
                F->refs[cr].hkl[1] = k;
                F->refs[cr].hkl[2] = l;
		if (h != 0 || k != 0 || l != 0) cr++;
            }
	F->nrefs = cr;
	fprintf(stdout, "Done. Generating F values...\n");
        if (!hklF_fill(F, cell)) { fprintf(stderr, "Error during F calulations!\n"); free(F->refs); free(F); free(cell->sym); free(cell->a); free(cell); exit(6); }
	fprintf(stdout, "Done. Opening %s for output...\n", argv[2]);
	ipf = fopen(argv[2], "w");
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
