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
#include <time.h>
#include <math.h>
#include "defines.h"

#define DAMPING    0.01	// Initial Gamma value
#define DERIV_STEP 1e-4
#define OSF_TOL    1e-8

// Bespoke routine to find minimal values of f' and f'' for the A and B
// sites of HBS. Requires the HBS model in INS form, and the relevant
// hkl one.
int main(int argc, char *argv[]) {
	unitcell *cell;
	hklF *F, *F2;
	int cr, i, s, bests, iter = 0;
	FILE *ipf;
	FP tmp;
	FP osf, GooF, GooFs[8], Gamma[4], LastGooF, ZnFP, ZnFPP, CuFP, CuFPP;

	if (argc < 3 || argc > 3) {
		fprintf(stderr, "Usage: search_hbs_f input.ins input.hkl\n");
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
	fprintf(stdout, "Done. Starting Search...\n");

	// Before going into iterative loop, extract the starting values of ZnFP, ZnFPP
	// CuFP, and CuFPP from the atom list
	for (i = 0; i < cell->natom; i++) {
		if (strncmp(cell->a[i].label, "ZN", 2) == 0) {
			ZnFP = cell->a[i].f[11];
			ZnFPP = cell->a[i].f[12];
		} else if (strncmp(cell->a[i].label, "CU", 2) == 0) {
			CuFP = cell->a[i].f[11];
			CuFPP = cell->a[i].f[12];
		}
	}
	LastGooF = 1e200;
	for (i = 0; i < 4; i++) Gamma[i] = DAMPING;

	// Iterative loop. Each iteration, we:
	// 1. Set f' and f'' on both sites to the current values, and compute
	//    the GooF and osf.
	// 2. Then we vary f' and f'' each site individually by +/- DERIV_STEP
	//    and compute new GooFs (and osfs).
	// 3. If there is no improvement in GooF along any direction, we are done.
	//    Otherwise, take a step of Gamma*(GooF-(GooF-better))
	//    along each improvement direction.
	// 4. Goto 1.
	// FUTURE: While there is no problem with this algorithm*, it can be
	// replaced by a most robust (against false minima) and faster algorithm
	// * It is gradient descent!
	while (1 == 1) {
		iter++;
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
		// Find new 'good' OSF and GooF
		if (!find_osf(&osf, &GooF, F, F2, 0.1, 0.0, OSF_TOL)) { fprintf(stderr, "Error during finding OSF!\n"); free(F2->refs); free(F2); free(F->refs); free(F); free(cell->sym); free(cell->a); free(cell); exit(11); }
                fprintf(stdout, "%4i, %8.4f, %8.4f, %8.4f, %8.4f, %12.8f, %.5E, %.5E\n", iter, ZnFP, ZnFPP, CuFP, CuFPP, GooF, GooF-LastGooF, osf);
		fflush(stdout);
		// Test all steps
		for (s = 0; s < 8; s++) {
			switch (s) {
				case 0:
					ZnFP -= DERIV_STEP;
					break;
				case 1:
					ZnFP += 2.0*DERIV_STEP;
					break;
				case 2:
					ZnFP -= DERIV_STEP;
					ZnFPP -= DERIV_STEP;
					break;
				case 3:
					ZnFPP += 2.0*DERIV_STEP;
					break;
				case 4:
					ZnFPP -= DERIV_STEP;
					CuFP -= DERIV_STEP;
					break;
				case 5:
					CuFP += 2.0*DERIV_STEP;
					break;
				case 6:
					CuFP -= DERIV_STEP;
					CuFPP -= DERIV_STEP;
					break;
				case 7:
					CuFPP += 2.0*DERIV_STEP;
					break;
			}
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
        	        if (!hklF_fill(F2, cell)) { fprintf(stderr, "Error during F calulations!\n"); free(F2->refs); free(F2); free(F->refs); free(F); free(cell->sym); free(cell->a); free(cell); exit(11); }
	                // Find new GooF
        	        if (!find_osf(&tmp, GooFs+s, F, F2, 0.1, 0.0, OSF_TOL)) { fprintf(stderr, "Error during finding OSF!\n"); free(F2->refs); free(F2); free(F->refs); free(F); free(cell->sym); free(cell->a); free(cell); exit(12); }
		}
		// Reset CuFPP
		CuFPP -= DERIV_STEP;
		// Calculate differences and find direction of greatest descent
		bests = -1; tmp = 0.0;
		for (s = 0; s < 8; s++)
			if ((GooF-GooFs[s])/fabs(GooF)-OSF_TOL > tmp) {
				bests = s; tmp = GooFs[s];
			}
		if (bests == -1 && fabs(LastGooF-GooF) < OSF_TOL) break; // DONE! -- no direction improves GooF
		// Take a step with the gradient along each direction
		for (s = 0; s < 8; s += 2) {
			// Update Gammas. Increase if GooF decreased, Decrease otherwise
			if (LastGooF-GooF > 0.0)
				Gamma[s>>1] *= 1.1;
			else
				Gamma[s>>1] /= 2.0;
			tmp = Gamma[s>>1]*(GooFs[s]-GooFs[s+1])/(2.0*DERIV_STEP);
			switch (s) {
				case 0: ZnFP  += tmp; break;
				case 2: ZnFPP += tmp; break;
				case 4: CuFP  += tmp; break;
				case 6: CuFPP += tmp; break;
			}
		}
		LastGooF = GooF;
	}
	fprintf(stdout, "Done. Best Values Obtained:\n");
	fprintf(stdout, "  GooF = %.5E\n", GooF);
	fprintf(stdout, "   OSF = %.5E\n", osf);
	fprintf(stdout, " A  f' = %8.4f f'' = %8.4f\n", ZnFP, ZnFPP);
	fprintf(stdout, " B  f' = %8.4f f'' = %8.4f\n", CuFP, CuFPP);

	free(F2->refs); free(F2); free(F->refs); free(F); free(cell->sym); free(cell->a); free(cell);
	return 0;
}
