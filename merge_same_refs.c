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

// Merge all identical reflections and output a new hkl file
int main(int argc, char *argv[]) {
	hklF *F;
	int cr;
	BOOL pairUp = BFALSE;
	FILE *ipf;

	if (argc < 3 || argc > 4) {
		fprintf(stderr, "Usage: merge_same_refs input.hkl output.hkl <friedel-pairs>\n");
		fprintf(stderr, "  <friedel-pairs> is optional, but if set to a non-zero value,\n");
		fprintf(stderr, "    puts F2 and sig-F2 for both members of a Friedel pair on the same\n");
		fprintf(stderr, "    output line.\n");
		fprintf(stderr, "\n");
		exit(1);
	}
	if (argc >= 4) { if (atoi(argv[3])) pairUp = BTRUE; else pairUp = BFALSE; }

        F = calloc(1, sizeof(hklF));
        if (!F) { fprintf(stderr, "Out of memory!\n"); exit(2); }
 
	ipf = fopen(argv[1], "r");
	if (!ipf) {
		fprintf(stderr, "Error opening %s for input. Does it exist?\n", argv[1]);
		exit(3);
	}

	fprintf(stdout, "Reading in reflections...\n");

        if (!read_hklF2(F, ipf)) { fprintf(stderr, "Error during read!\n"); free(F); fclose(ipf); exit(4); }
        fclose(ipf);

	fprintf(stdout, "%6i reflections read. Merging |F|^2 values...\n", F->nrefs);
	if (!hklF_merge_identical(F)) { fprintf(stderr, "Error during Merge!\n"); free(F->refs); free(F); exit(5); }
	fprintf(stdout, "Done. %6i reflections after merging. Opening %s for output...\n", F->nrefs, argv[2]);
	ipf = fopen(argv[2], "w");
	if (!ipf) { fprintf(stderr, "Error opening output file!\n"); free(F->refs); free(F); exit(6); }
	if (!pairUp) {
		fprintf(stdout, "Writing Reflections...\n");
		print_hklF2(F, ipf);
	} else {
		fprintf(stdout, "Writing Reflections With Friedel Pairs on Same Line...\n");
        	print_hklF2_paired(F, ipf);
	}
	fclose(ipf);
	fprintf(stdout, "Success.\n");

	free(F->refs); free(F);
	return 0;
}
