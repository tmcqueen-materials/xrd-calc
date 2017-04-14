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
#include <string.h>
#include "defines.h"

void print_matrix(mtx33 m) {
	if (!m) return;
	fprintf(stdout, " [ %11.5e %11.5e %11.5e ]\n", m[0][0], m[0][1], m[0][2]);
        fprintf(stdout, " [ %11.5e %11.5e %11.5e ]\n", m[1][0], m[1][1], m[1][2]);
        fprintf(stdout, " [ %11.5e %11.5e %11.5e ]\n", m[2][0], m[2][1], m[2][2]);
	return;
}

void print_cell(unitcell *cell) {
	if (!cell) return;
	fprintf(stdout, "Real space cell:\n");
	fprintf(stdout, "  a = %11.6f A, alpha = %9.5f deg\n", cell->abc[0], cell->albega[0]*fromRad);
        fprintf(stdout, "  b = %11.6f A,  beta = %9.5f deg\n", cell->abc[1], cell->albega[1]*fromRad);
        fprintf(stdout, "  c = %11.6f A, gamma = %9.5f deg\n", cell->abc[2], cell->albega[2]*fromRad);
	fprintf(stdout, " Volume = %11.5f A^3\n\n", cell->vol);
        fprintf(stdout, "Reciprocal space cell:\n");
        fprintf(stdout, "  a = %11.8f A, alpha = %9.5f deg\n", cell->invabc[0], cell->invalbega[0]*fromRad);
        fprintf(stdout, "  b = %11.8f A,  beta = %9.5f deg\n", cell->invabc[1], cell->invalbega[1]*fromRad);
        fprintf(stdout, "  c = %11.8f A, gamma = %9.5f deg\n", cell->invabc[2], cell->invalbega[2]*fromRad);
        fprintf(stdout, " Volume = %11.7f A^-3\n", cell->invvol);
	fprintf(stdout, " G* = \n");
	print_matrix(cell->invG);
	fprintf(stdout, " A* = \n");
	print_matrix(cell->invA);
	fprintf(stdout, "Num Symmetry Elements %5i\n", cell->nsym);
	fprintf(stdout, "Num Atoms = %5i\n", cell->natom);
	return;
}

// hklF file has magF and phaseF (deg)
// files are terminated by an entry of 0 0 0, so we skip those until the end
void print_hklF(hklF *F, FILE *file) {
	char lastline[4096] = { '\0' };
	int i;
	hklF_uno *cF;
	FP magF, phF;

	if (!F || !file) return;
	for (i = 0; i < F->nrefs; i++) {
		cF = &(F->refs[i]);
		magF = CMMAG(cF->F);
		phF = CMPH(cF->F);
		if (fabs(magF) < 1e-2) { magF = 0.0; phF = 0.0; } // flush values that are zero in format to 0
		if ((int)(cF->hkl[0]) != 0 || (int)(cF->hkl[1]) != 0 || (int)(cF->hkl[2]) != 0) {
		  if (magF < 99999.99)
		    fprintf(file, "%4i%4i%4i%8.2f%8.2f\n", (int)(cF->hkl[0]), (int)(cF->hkl[1]), (int)(cF->hkl[2]), magF, phF*fromRad);
		  else if (magF < 9999999.)
                    fprintf(file, "%4i%4i%4i%8.0f%8.2f\n", (int)(cF->hkl[0]), (int)(cF->hkl[1]), (int)(cF->hkl[2]), magF, phF*fromRad);
		  else
		    fprintf(stdout, "Could not write reflection (%4i%4i%4i) to hkl file. magF = %.8E > 9999999!\n", (int)(cF->hkl[0]), (int)(cF->hkl[1]), (int)(cF->hkl[2]), magF);
		} else { // a 0 0 0 entry, so add to lastline
		    snprintf(lastline+strlen(lastline),4095-strlen(lastline), "%4i%4i%4i%8.0f%8.0f\n", 0, 0, 0, magF, phF*fromRad);
		}
	}
	// Terminating line
	if (strlen(lastline) > 0)
		fprintf(file, "%s\n", lastline);
	else
		fprintf(file, "%4i%4i%4i%8.2f%8.2f\n\n", 0, 0, 0, 0.0, 0.0);
	return;
}

// hklF2 file is usual modern hkl file with |F|^2 and sigma-|F|^2
// files are terminated by an entry of 0 0 0, so we skip those until the end
void print_hklF2(hklF *F, FILE *file) {
        char lastline[4096] = { '\0' };
        int i;
        hklF_uno *cF;

        if (!F || !file || !F->refs || F->nrefs < 1 || F->nrefs > MAX_REFS) return;
        for (i = 0; i < F->nrefs; i++) {
                cF = &(F->refs[i]);
		if ((int)(cF->hkl[0]) != 0 || (int)(cF->hkl[1]) != 0 || (int)(cF->hkl[2]) != 0) {
		  if (cF->F2 > 9999999. || cF->sigF2 > 9999999.)
		    fprintf(stdout, "Could not write reflection (%4i%4i%4i) to hkl file. F2 = %.8E or sigF2 = %.8E > 9999999!\n", (int)(cF->hkl[0]), (int)(cF->hkl[1]), (int)(cF->hkl[2]), cF->F2, cF->sigF2);
		  else if (cF->F2 < -999999. || cF->sigF2 < -999999.)
		    fprintf(stdout, "Could not write reflection (%4i%4i%4i) to hkl file. F2 = %.8E or sigF2 = %.8E < -999999!\n", (int)(cF->hkl[0]), (int)(cF->hkl[1]), (int)(cF->hkl[2]), cF->F2, cF->sigF2);
		  else if (cF->F2 < -9999.99 || cF->sigF2 < -9999.99 || cF->F2 > 99999.99 || cF->sigF2 > 99999.99)
                    fprintf(file, "%4i%4i%4i%8.0f%8.0f\n", (int)(cF->hkl[0]), (int)(cF->hkl[1]), (int)(cF->hkl[2]), cF->F2, cF->sigF2);
		  else
                    fprintf(file, "%4i%4i%4i%8.2f%8.2f\n", (int)(cF->hkl[0]), (int)(cF->hkl[1]), (int)(cF->hkl[2]), cF->F2, cF->sigF2);
		} else { // a 0 0 0 entry, so add to lastline
		    snprintf(lastline+strlen(lastline),4095-strlen(lastline), "%4i%4i%4i%8.0f%8.0f\n", 0, 0, 0, cF->F2, cF->sigF2);
		}
        }
        // Terminating line
        if (strlen(lastline) > 0)
                fprintf(file, "%s\n", lastline);
        else
                fprintf(file, "%4i%4i%4i%8.2f%8.2f\n\n", 0, 0, 0, 0.0, 0.0);
        return;
}

// Just like above, but put Friedel pairs on the same line
void print_hklF2_paired(hklF *F, FILE *file) {
        char lastline[4096] = { '\0' };
        int i, j;
        hklF_uno *cF, *cF2;
	BOOL *used;

        if (!F || !file || !F->refs || F->nrefs < 1 || F->nrefs > MAX_REFS) return;

	used = calloc(F->nrefs, sizeof(BOOL));

        for (i = 0; i < F->nrefs; i++) {
	    if (!used[i]) {
                cF = &(F->refs[i]);
		cF2 = NULL;
		for (j = i+1; j < F->nrefs; j++)
			if ((int)(cF->hkl[0]) == -1*(int)(F->refs[j].hkl[0]) &&
			    (int)(cF->hkl[1]) == -1*(int)(F->refs[j].hkl[1]) &&
			    (int)(cF->hkl[2]) == -1*(int)(F->refs[j].hkl[2])) break;
		if (j < F->nrefs) { cF2 = &(F->refs[j]); used[j] = BTRUE; }

		// Do cF
                if ((int)(cF->hkl[0]) != 0 || (int)(cF->hkl[1]) != 0 || (int)(cF->hkl[2]) != 0) {
                  if (cF->F2 > 9999999. || cF->sigF2 > 9999999.)
                    fprintf(stdout, "Could not write reflection (%4i%4i%4i) to hkl file. F2 = %.8E or sigF2 = %.8E > 9999999!\n", (int)(cF->hkl[0]), (int)(cF->hkl[1]), (int)(cF->hkl[2]), cF->F2, cF->sigF2);
                  else if (cF->F2 < -999999. || cF->sigF2 < -999999.)
                    fprintf(stdout, "Could not write reflection (%4i%4i%4i) to hkl file. F2 = %.8E or sigF2 = %.8E < -999999!\n", (int)(cF->hkl[0]), (int)(cF->hkl[1]), (int)(cF->hkl[2]), cF->F2, cF->sigF2);
                  else if (cF->F2 < -9999.99 || cF->sigF2 < -9999.99 || cF->F2 > 99999.99 || cF->sigF2 > 99999.99)
                    fprintf(file, "%4i%4i%4i%8.0f%8.0f", (int)(cF->hkl[0]), (int)(cF->hkl[1]), (int)(cF->hkl[2]), cF->F2, cF->sigF2);
                  else
                    fprintf(file, "%4i%4i%4i%8.2f%8.2f", (int)(cF->hkl[0]), (int)(cF->hkl[1]), (int)(cF->hkl[2]), cF->F2, cF->sigF2);
                } else { // a 0 0 0 entry, so add to lastline
                    snprintf(lastline+strlen(lastline),4095-strlen(lastline), "%4i%4i%4i%8.0f%8.0f\n", 0, 0, 0, cF->F2, cF->sigF2);
                }

		// then CF2
		if (cF2) {
                 if ((int)(cF2->hkl[0]) != 0 || (int)(cF2->hkl[1]) != 0 || (int)(cF2->hkl[2]) != 0) {
                  if (cF2->F2 > 9999999. || cF2->sigF2 > 9999999.)
                    fprintf(stdout, "Could not write reflection (%4i%4i%4i) to hkl file. F2 = %.8E or sigF2 = %.8E > 9999999!\n", (int)(cF2->hkl[0]), (int)(cF2->hkl[1]), (int)(cF2->hkl[2]), cF2->F2, cF2->sigF2);
                  else if (cF2->F2 < -999999. || cF2->sigF2 < -999999.)
                    fprintf(stdout, "Could not write reflection (%4i%4i%4i) to hkl file. F2 = %.8E or sigF2 = %.8E < -999999!\n", (int)(cF2->hkl[0]), (int)(cF2->hkl[1]), (int)(cF2->hkl[2]), cF2->F2, cF2->sigF2);
                  else if (cF2->F2 < -9999.99 || cF2->sigF2 < -9999.99 || cF2->F2 > 99999.99 || cF2->sigF2 > 99999.99)
                    fprintf(file, "%8.0f%8.0f", cF2->F2, cF2->sigF2);
                  else
                    fprintf(file, "%8.2f%8.2f", cF2->F2, cF2->sigF2);
		  fprintf(file, "\n");
                 } else { // a 0 0 0 entry, so add to lastline
                    snprintf(lastline+strlen(lastline),4095-strlen(lastline), "%4i%4i%4i%8.0f%8.0f\n", 0, 0, 0, cF2->F2, cF2->sigF2);
                 }
		} else fprintf(file, "\n");
	    }
        }

        // Terminating line
        if (strlen(lastline) > 0)
                fprintf(file, "%s\n", lastline);
        else
                fprintf(file, "%4i%4i%4i%8.2f%8.2f\n\n", 0, 0, 0, 0.0, 0.0);

	free(used);

        return;
}
