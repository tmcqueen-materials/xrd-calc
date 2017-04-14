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
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "defines.h"

// Internal helper to read a line and automatically convert to uppercase
// It also handles the fortran convention of an = in column 79 meaning
// a line continue.
int freadline(char *result, int maxresult, FILE *file) {
	int c, i = 0, j = 0;
	char *lre;
	if (!result || !file || maxresult < 1 || feof(file)) return EOF;

read_a_line:
	c = getc(file);
	while (c != '\n' && c != '\r' && c != EOF && i < maxresult-1) {
		result[i] = c; i++;
		if (result[i] >= 97 && result[i] <= 122) result[i] -= 32; // a-z to upper case
		c = getc(file);
	}
	while (c == '\n' || c == '\r') { c = getc(file); } fseek(file, -1, SEEK_CUR); // gracefully deal with CR/LF, NL, CR line endings
	//OLD strict one if (i >= 79 && result[78] == '=') { result[78] = ' '; i = 79; goto read_a_line; } // read another line
        // be like new SHELX and less paranoid about finding an = sign at the right column, instead only requiring it at the end of the line
	result[i] = '\0'; 
	if (i == maxresult-1) goto finished_reading_line;
	lre = strrchr(result, '=');
	if (lre != NULL) {
		for (j = 0; lre[j] != '\0'; j++)
			if (lre[j] != ' ' && lre[j] != '=') goto finished_reading_line;
		i = lre-result; result[i] = ' '; i++; // replace '=' with ' ' and set position at next character
		goto read_a_line;
	}
 finished_reading_line:
	result[i] = '\0';
	return i;
}

// Internal helper to parse a portion of the SYMM line
// Right now this only handles forms like -X or X-Y or -X+Y
// (no coefficients)
BOOL parse_symm_el(FP *rv, char *str, int i, int ilast) {
        if (i == ilast || (i == ilast+1 && str[ilast] == '+')) {
                rv[0] = 1.0;
        } else if (i == ilast+1 && str[ilast] == '-') {
                rv[0] = -1.0;
        } else {
        	// TODO!
		return BFALSE;
        }
	return BTRUE;
}
BOOL parse_symm(symm *rv, char *sym, int idx) {
	char rstr[4096] = { '\0' };
	int i, ilast;
	BOOL sx = BFALSE, sy = BFALSE, sz = BFALSE;
	if (!rv || !sym || idx < 0 || idx > 2) return BFALSE;
	rv->r[idx][0] = 0.0; rv->r[idx][1] = 0.0; rv->r[idx][2] = 0.0;
	for (i = 0, ilast = 0; i < strlen(sym); i++)
		if (sym[i] != ' ') { rstr[ilast] = sym[i]; ilast++; }
	for (i = 0, ilast = 0; i < strlen(rstr); i++) {
		switch (rstr[i]) {
			case 'X':
				if (sx) return BFALSE;
				sx = BTRUE;
				if (!parse_symm_el(&(rv->r[idx][0]), rstr, i, ilast)) return BFALSE;
				ilast = i+1;
				break;
 			case 'Y':
                                if (sy) return BFALSE;
                                sy = BTRUE;
                                if (!parse_symm_el(&(rv->r[idx][1]), rstr, i, ilast)) return BFALSE;
				ilast = i+1;
 				break;
			case 'Z':
                                if (sz) return BFALSE;
                                sz = BTRUE;
                                if (!parse_symm_el(&(rv->r[idx][2]), rstr, i, ilast)) return BFALSE;
				ilast = i+1;
				break;
			default:
				break;
		}
	}
	return BTRUE;
}

// Generate all symmetry-related atoms with the given input parameters, and store in a
BOOL atoms_fill(int *ngen, atom *a, int maxgen, double x, double y, double z, double sof, double U11, double U22, double U33, double U23, double U13, double U12, char *label, symm *sym, int nsym, scatf f, unitcell *cell) {
	int cursym, i;
	atom fa;
	mtx33 sot;
	if (!ngen || !a || !f || maxgen < 1 || nsym < 1 || !sym) return BFALSE;
	*ngen = 0;
	if (sof < 0.0 || sof > 1.0) return BFALSE;

	// Then generate symmetry equivalents
	for (cursym = 0; cursym < nsym && *ngen < maxgen; cursym++) {
        	// Fill atom
        	a[*ngen].occ = sof*(double)(nsym);
	        a[*ngen].xyz[0] = x; a[*ngen].xyz[1] = y; a[*ngen].xyz[2] = z;
        	a[*ngen].uij[0][0] = U11; a[*ngen].uij[0][1] = U12; a[*ngen].uij[0][2] = U13;
        	a[*ngen].uij[1][0] = U12; a[*ngen].uij[1][1] = U22; a[*ngen].uij[1][2] = U23;
	        a[*ngen].uij[2][0] = U13; a[*ngen].uij[2][1] = U23; a[*ngen].uij[2][2] = U33;
        	strncpy(a[*ngen].label, label, MAX_LABEL_LEN-1);
        	if (!SCATFCOPY(a[*ngen].f, f)) return BFALSE;
		// Apply symmetry operator
		if (!MVpV2(a[*ngen].xyz, sym[cursym].r, a[*ngen].xyz, sym[cursym].t)) return BFALSE;
		// Check for uniqueness
		for (i = 0; i < *ngen; i++) if (VEQ_FRAC(a[*ngen].xyz, a[i].xyz)) goto next_sym;
		// Apply rotation part of transform to Uij
		if (!Mt(sot, sym[cursym].r)) return BFALSE;
		if (!M3M2M(a[*ngen].uij, sym[cursym].r, a[*ngen].uij, sot)) return BFALSE;
		*ngen += 1;
next_sym:
		cursym = cursym;
	}

	// Deal with special case of this being the last loaded atom
	if (cursym < nsym) {
            for (;cursym < nsym; cursym++) {
                // Fill atom
                fa.xyz[0] = x; fa.xyz[1] = y; fa.xyz[2] = z;
                // Apply symmetry operator
                if (!MVpV2(fa.xyz, sym[cursym].r, fa.xyz, sym[cursym].t)) return BFALSE;
                // Check for uniqueness
                for (i = 0; i < *ngen; i++) if (VEQ_FRAC(fa.xyz, a[i].xyz)) goto next_sym_2;
		fprintf(stdout, "Generated more atoms than were listed in the UNIT instruction.\n"); return BFALSE;
next_sym_2:
		cursym = cursym;
            }
	}

	if (*ngen < 1) return BFALSE;

	// Adjust occupancies on all created atoms
	for (i = 0; i < *ngen; i++) a[i].occ /= (FP)(*ngen);
        if (a[0].occ < 0.0 || fabs(a[0].occ-2.0*COORD_TOL) > 1.0) fprintf(stdout, "Warning: Atom %s has unphysical occupancy of %.2f (versus special position)!\n", label, a[0].occ); 

	return BTRUE;
}

BOOL read_ins_cell(unitcell *cell, FILE *ipf) {
	char curline[4096], scratch[4096];
	char *ptr;
	int LATT = 0;
	int nt, i, j;
	int natomgen = 0;
	BOOL centro = BTRUE;
	BOOL sCELL = BFALSE, sUNIT = BFALSE;
	float tmp[12];
	float lastU[6] = { 0.0 };
	float fv[MAX_FREE_VARS] = { 0.0 };
	int nFV = 0;
	symm sym[MAX_SYMM];
	int nsym = 0;
	scatf f[MAX_ATOM_TYPES];
	char fl[MAX_ATOM_TYPES][MAX_LABEL_LEN] = { { 0 } };
	int nf = 0;

	if (!cell || !ipf || cell->natom != 0 || cell->a != NULL) return BFALSE;
	while (!sUNIT && freadline(curline, 4096, ipf) != EOF) {

		if (strncmp(curline, "CELL", 4) == 0) {
			if (sCELL) { fprintf(stdout, "Invalid INS file. Has multiple CELL lines!\n"); return BFALSE; }
			sCELL = BTRUE;
			// wavelength, a, b, c, alpha, beta, gamma
			if (sscanf(curline, "CELL %f %f %f %f %f %f %f", tmp, tmp+1, tmp+2, tmp+3, tmp+4, tmp+5, tmp+6) != 7) {
				fprintf(stdout, "Malformed CELL line in INS file.\n"); return BFALSE; }
			cell->abc[0] = tmp[1]; cell->abc[1] = tmp[2]; cell->abc[2] = tmp[3];
			if (cell->abc[0] < 0.0 || cell->abc[1] < 0.0 || cell->abc[2] < 0.0) {
				fprintf(stdout, "Malformed CELL line in INS file.\n"); return BFALSE; }
			cell->albega[0] = tmp[4]*toRad; cell->albega[1] = tmp[5]*toRad; cell->albega[2] = tmp[6]*toRad;
                        if (cell->albega[0] < 0.0 || cell->albega[1] < 0.0 || cell->albega[2] < 0.0) {
                                fprintf(stdout, "Malformed CELL line in INS file.\n"); return BFALSE; }
		}

		if (strncmp(curline, "LATT", 4) == 0) {
			if (LATT != 0) { fprintf(stdout, "Invalid INS file. Has multiple LATT lines!\n"); return BFALSE; }
			if (sscanf(curline, "LATT %i", &LATT) != 1) {
				fprintf(stdout, "Malformed LATT line in INS file.\n"); return BFALSE; }
			if (LATT < 0) { LATT = -1*LATT; centro = BFALSE; }
			if (LATT < 1 || LATT > 7) {
                                fprintf(stdout, "Malformed LATT line in INS file.\n"); return BFALSE; }
		}

		if (strncmp(curline, "SYMM", 4) == 0) {
			ptr = curline+4;
			sym[nsym].r[0][0] = 0.0; sym[nsym].r[0][1] = 0.0; sym[nsym].r[0][2] = 0.0;
                        sym[nsym].r[1][0] = 0.0; sym[nsym].r[1][1] = 0.0; sym[nsym].r[1][2] = 0.0;
                        sym[nsym].r[2][0] = 0.0; sym[nsym].r[2][1] = 0.0; sym[nsym].r[2][2] = 0.0;
			sym[nsym].t[0] = 0.0; sym[nsym].t[1] = 0.0; sym[nsym].t[2] = 0.0;
			if (!strchr(ptr, ',')) {
                                fprintf(stdout, "Malformed SYMM line in INS file.\n"); return BFALSE; }
			else { strchr(ptr, ',')[0] = '\0'; }
			if (!parse_symm(&(sym[nsym]), ptr, 0)) {
                                fprintf(stdout, "Malformed SYMM line in INS file.\n"); return BFALSE; }
                        ptr = ptr+strlen(ptr)+1;
                        if (!strchr(ptr, ',')) {
                                fprintf(stdout, "Malformed SYMM line in INS file.\n"); return BFALSE; }
                        else { strchr(ptr, ',')[0] = '\0'; }
                        if (!parse_symm(&(sym[nsym]), ptr, 1)) {
                                fprintf(stdout, "Malformed SYMM line in INS file.\n"); return BFALSE; }
  			ptr = ptr+strlen(ptr)+1;
                        if (!parse_symm(&(sym[nsym]), ptr, 2)) {
                                fprintf(stdout, "Malformed SYMM line in INS file.\n"); return BFALSE; }
			nsym++;
		}

		if (strncmp(curline, "SFAC", 4) == 0) {
			if (sUNIT) { fprintf(stdout, "Unexpected SFAC after UNIT in INS file.\n"); return BFALSE; }
			// Try long format
			if (sscanf(curline, "SFAC %s %f", scratch, tmp) == 2) {
				// Long
				if (sscanf(curline, "SFAC %s %f %f %f %f %f %f %f %f %f %f %f %f", fl[nf], tmp, tmp+1, tmp+2, tmp+3,
				    tmp+4, tmp+5, tmp+6, tmp+7, tmp+8, tmp+9, tmp+10, tmp+11) != 12) {
					fprintf(stdout, "Malformed SFAC line in INS file. a1...b4,c,f',f'', and mu are all REQUIRED.\n"); return BFALSE; }
				f[nf][0] = tmp[0]; f[nf][1] = tmp[1];
				f[nf][2] = tmp[2]; f[nf][3] = tmp[3];
				f[nf][4] = tmp[4]; f[nf][5] = tmp[5];
				f[nf][6] = tmp[6]; f[nf][7] = tmp[7];
				f[nf][8] = 0.0;    f[nf][9] = 0.0;
				f[nf][10] = tmp[8];
				f[nf][11] = tmp[9];
				f[nf][12] = tmp[10];
				nf++;
			} else {
				// Short
				nt = sscanf(curline, "SFAC %s %s %s %s %s %s %s %s %s %s %s %s", scratch,
					scratch+16,scratch+32,scratch+48,scratch+64,scratch+80,scratch+96,
					scratch+112,scratch+128,scratch+144,scratch+160,scratch+172);
				for (i = 0; i < nt; i++) {
					strncpy(fl[nf], scratch+16*i, MAX_LABEL_LEN-1);
					if (!scatf_fill(f[nf], fl[nf])) {
						fprintf(stdout, "Could not find scattering factor for %s in internal database! Use the long form of SFAC!\n", scratch+16*i);
						return BFALSE;
					}
					nf++;
				}
			} 
		}

		if (strncmp(curline, "DISP", 4) == 0) {
                        if (sUNIT) { fprintf(stdout, "Unexpected DISP after UNIT in INS file.\n"); return BFALSE; }
			if (sscanf(curline, "DISP %s %f %f %f", scratch, tmp, tmp+1, tmp+2) != 4) {
				fprintf(stdout, "Malformed DISP instruction in INS file!\n"); return BFALSE; }
			nt = 0;
			for (i = 0; i < nf && !nt; i++) {
				if (strcmp(fl[i], scratch) == 0) {
					f[i][11] = tmp[0]; f[i][12] = tmp[1]; nt = 1;
				}
			}
			if (!nt) {
                                fprintf(stdout, "Could not find element %s from DISP instruction in INS file!\n", scratch); return BFALSE; }
		}

                if (strncmp(curline, "UNIT", 4) == 0) {
			sUNIT = BTRUE;
			if (nf < 1) {
				fprintf(stdout, "Unexpected UNIT instruction (it must be after SFAC and DISP lines)!\n"); return BFALSE; }
			ptr = curline+4;
			nt = nf;
			while (nt > 0) {
				if (!sscanf(ptr, "%i", &i)) {
					fprintf(stdout, "Malformed UNIT instruction!\n"); return BFALSE; }
				cell->natom += i;
				nt--; ptr = strchr(ptr+1, ' ');
                                if (!ptr && nt > 0) {
                                        fprintf(stdout, "Malformed UNIT instruction!\n"); return BFALSE; }
			}
                }
	}

	if (!sCELL || nf == 0 || !sUNIT || LATT == 0 || cell->natom == 0) return BFALSE;

	if (!cell_fill(cell)) return BFALSE;

//        fprintf(stdout, "Global INS Information:\n");
//        print_cell(cell);
//        fprintf(stdout, "LATT = %i, nsym_read = %i\n", LATT, nsym);

	// Generate all symmetry operations
	// add unity
	sym[nsym].r[0][0] = 1.0; sym[nsym].r[0][1] = 0.0; sym[nsym].r[0][2] = 0.0;
        sym[nsym].r[1][0] = 0.0; sym[nsym].r[1][1] = 1.0; sym[nsym].r[1][2] = 0.0;
        sym[nsym].r[2][0] = 0.0; sym[nsym].r[2][1] = 0.0; sym[nsym].r[2][2] = 1.0;
	sym[nsym].t[0] = 0.0; sym[nsym].t[1] = 0.0; sym[nsym].t[2] = 0.0;
	nsym++;
	nt = nsym;
	// I-centered
	for (i = 0; i < nt && LATT == 2; i++) {
		if (!SYMMCOPY(&(sym[nsym]), &(sym[i]))) return BFALSE;
		sym[nsym].t[0] += 0.5; sym[nsym].t[1] += 0.5; sym[nsym].t[2] += 0.5;
		nsym++;
	}
	// R-centered (trigonal cell)
        for (i = 0; i < nt && LATT == 3; i++) {
                if (!SYMMCOPY(&(sym[nsym]), &(sym[i]))) return BFALSE; 
                sym[nsym].t[0] += 2.0/3.0; sym[nsym].t[1] += 1.0/3.0; sym[nsym].t[2] += 1.0/3.0;
                nsym++;
                if (!SYMMCOPY(&(sym[nsym]), &(sym[i]))) return BFALSE;
                sym[nsym].t[0] += 1.0/3.0; sym[nsym].t[1] += 2.0/3.0; sym[nsym].t[2] += 2.0/3.0;
                nsym++;
        }
 	// F-centered
        for (i = 0; i < nt && LATT == 4; i++) {
                if (!SYMMCOPY(&(sym[nsym]), &(sym[i]))) return BFALSE; 
                sym[nsym].t[0] += 0.5; sym[nsym].t[1] += 0.5; sym[nsym].t[2] += 0.0;
                nsym++;
                if (!SYMMCOPY(&(sym[nsym]), &(sym[i]))) return BFALSE;
                sym[nsym].t[0] += 0.5; sym[nsym].t[1] += 0.0; sym[nsym].t[2] += 0.5;
                nsym++;
                if (!SYMMCOPY(&(sym[nsym]), &(sym[i]))) return BFALSE;
                sym[nsym].t[0] += 0.0; sym[nsym].t[1] += 0.5; sym[nsym].t[2] += 0.5;
                nsym++;
        }
	// A-centered
        for (i = 0; i < nt && LATT == 5; i++) {
                if (!SYMMCOPY(&(sym[nsym]), &(sym[i]))) return BFALSE;
                sym[nsym].t[0] += 0.0; sym[nsym].t[1] += 0.5; sym[nsym].t[2] += 0.5;
                nsym++;
        }
	// B-centered
        for (i = 0; i < nt && LATT == 6; i++) {
                if (!SYMMCOPY(&(sym[nsym]), &(sym[i]))) return BFALSE;
                sym[nsym].t[0] += 0.5; sym[nsym].t[1] += 0.0; sym[nsym].t[2] += 0.5;
                nsym++;
        }
	// C-centered
        for (i = 0; i < nt && LATT == 7; i++) {
                if (!SYMMCOPY(&(sym[nsym]), &(sym[i]))) return BFALSE;
                sym[nsym].t[0] += 0.5; sym[nsym].t[1] += 0.5; sym[nsym].t[2] += 0.0;
                nsym++;
        }
	// Inversion to ALL of the above
	nt = nsym;
        for (i = 0; i < nt && centro; i++) {
                if (!SYMMCOPY(&(sym[nsym]), &(sym[i]))) return BFALSE;
		sym[nsym].r[0][0] *= -1.0; sym[nsym].r[0][1] *= -1.0; sym[nsym].r[0][2] *= -1.0;
                sym[nsym].r[1][0] *= -1.0; sym[nsym].r[1][1] *= -1.0; sym[nsym].r[1][2] *= -1.0;
                sym[nsym].r[2][0] *= -1.0; sym[nsym].r[2][1] *= -1.0; sym[nsym].r[2][2] *= -1.0;
                sym[nsym].t[0] *= -1.0; sym[nsym].t[1] *= -1.0; sym[nsym].t[2] *= -1.0;
                nsym++;
        }
 
//	fprintf(stdout, "Total symmetry operations = %i\n", nsym);

	cell->a = calloc(cell->natom, sizeof(atom));
	if (!cell->a) return BFALSE;
	cell->sym = calloc(nsym, sizeof(symm));
	if (!cell->sym) { free(cell->a); return BFALSE; }
	for (i = 0; i < nsym; i++)
		if (!SYMMCOPY(cell->sym+i,sym+i)) { free(cell->sym); free(cell->a); return BFALSE; }
	cell->nsym = nsym;

	// Continue processing to get atoms
	while (freadline(curline, 4096, ipf) != EOF && strncmp(curline, "HKLF", 4) != 0) {
		if (strncmp(curline, "RESI", 4) == 0) {
			fprintf(stdout, "Structures defined in multiple parts (RESI) not supported!\n"); return BFALSE; }
		if (strncmp(curline, "FVAR", 4) == 0) {
			if (nFV) {
				fprintf(stdout, "Unexpected second FVAR line!\n"); return BFALSE; }
			ptr = curline+4;
			while (nFV < MAX_FREE_VARS && ptr != NULL) {
				if (sscanf(ptr, "%f", fv+nFV) == 1)
					nFV++;
				ptr = strchr(ptr, '.');
				if (ptr) ptr = strchr(ptr, ' ');
			}
			if (nFV >= MAX_FREE_VARS && ptr != NULL)
				fprintf(stdout, "Warning: Maximum number of free variables exceeded.\n");
		}
		if (strncmp(curline, "REM", 3) != 0) { // skip if a comment, otherwise process
			if (sscanf(curline, "%s %i %f %f %f %f %f %f %f %f %f %f", scratch, &i, 
			    tmp, tmp+1, tmp+2, tmp+3, tmp+4, tmp+5, tmp+6, tmp+7, tmp+8, tmp+9) == 12) {
				// a long atom definition line
				if (i < 1 || i > nf) {
					fprintf(stdout, "Malformed atom line %s. Invalid scattering number. IGNORING.\n", scratch);
				} else {
					for (j = 0; j < 10; j++) {
						if (fabs(tmp[j]) > 10.0) {
							nt = (int)(fabs(tmp[j])/10.0);
							tmp[10] = fabs(tmp[j])-10.0*(double)nt;
							if (nt-2 < nFV && nt-2 >= 0) {
								if (tmp[j] < 0.0)
									tmp[j] = tmp[10]*(1.0-fv[nt-2]);
								else
									tmp[j] = tmp[10]*fv[nt-2];
							} else tmp[j] = tmp[10];
							if (j > 3) lastU[j-4] = tmp[j];
						} else if (j > 3) {
							if (tmp[j] < -0.49999 && tmp[j] > -5.00001)
								tmp[j] = fabs(tmp[j])*lastU[j-4];
							else
								lastU[j-4] = tmp[j];
						}
					}
					nt = 0;
					if (!atoms_fill(&nt, cell->a+natomgen, cell->natom-natomgen, tmp[0], tmp[1], tmp[2], tmp[3], tmp[4], tmp[5], tmp[6], tmp[7], tmp[8], tmp[9], scratch, sym, nsym, f[i-1], cell)) 
						fprintf(stdout, "Unable to generate symmetry equivalents of %s. IGNORING.\n", scratch);
					else natomgen += nt;
				}
			} else if (sscanf(curline, "%s %i %f %f %f %f %f", scratch, &i, tmp, tmp+1, tmp+2, tmp+3, tmp+4) == 7) {
				// a short atom definition line
                                if (i < 1 || i > nf) {
                                        fprintf(stdout, "Malformed atom line %s. Invalid scattering number. IGNORING.\n", scratch);
                                } else {
                                        for (j = 0; j < 5; j++) {
                                                if (fabs(tmp[j]) > 10.0) {
                                                        nt = (int)(fabs(tmp[j])/10.0);
                                                        tmp[10] = fabs(tmp[j])-10.0*(double)nt;
                                                        if (nt-2 < nFV && nt-2 >= 0) {
                                                                if (tmp[j] < 0.0)
                                                                        tmp[j] = tmp[10]*(1.0-fv[nt-2]);
                                                                else
                                                                        tmp[j] = tmp[10]*fv[nt-2];
							} else tmp[j] = tmp[10];
                                                        if (j > 3) lastU[j-4] = tmp[j];
                                                } else if (j > 3) {
                                                        if (tmp[j] < -0.49999 && tmp[j] > -5.00001)
                                                                tmp[j] = fabs(tmp[j])*lastU[j-4];
                                                        else
                                                                lastU[j-4] = tmp[j];
                                                }
                                        }
                                        nt = 0;
					// compute U in non-cartesian basis
					tmp[5] = tmp[4]; tmp[6] = tmp[4];
					tmp[7] = tmp[4]*cos(cell->invalbega[0]);
					tmp[8] = tmp[4]*cos(cell->invalbega[1]);
					tmp[9] = tmp[4]*cos(cell->invalbega[2]);
                                        if (!atoms_fill(&nt, cell->a+natomgen, cell->natom-natomgen, tmp[0], tmp[1], tmp[2], tmp[3], tmp[4], tmp[5], tmp[6], tmp[7], tmp[8], tmp[9], scratch, sym, nsym, f[i-1], cell))
                                                fprintf(stdout, "Unable to generate symmetry equivalents of %s. IGNORING.\n", scratch);
					else natomgen += nt;
                                }
			} else if (sscanf(curline, "%s %i %f %f %f", scratch, &i, tmp, tmp+1, tmp+2) == 5) {
				// catch possible atom specification lines and warn
				fprintf(stdout, "Warning: %s IGNORED. You must have sof and Uiso or Uaniso.\n", scratch);
			}
		}
	}

	// Warn if atom numbers don't math
	if (cell->natom != natomgen) {
		fprintf(stdout, "Warning: Total atoms specified on UNIT line %i but total atoms generated %i.\n\n", cell->natom, natomgen);
		cell->natom = natomgen;
	}

	return BTRUE;
}


// Read a 'normal' HKL file (h, k, l, F2, sig-F2)
// Note we CANNOT just use sscanf to process lines directly, as it treats whitespace
// differently than fortran, which this is compatible with.
#define STEP_REFS 100000
BOOL read_hklF2(hklF *rv, FILE *ipf) {
	char curline[4096] = { '\0' }, tmpchr, totlen;
	int nrefalloc = 0, nr = 0;
	BOOL end = BFALSE;
	int hkl[3];
	float tmp[2];
	hklF_uno *tmpPtr;
	if (!rv || !ipf || rv->nrefs != 0 || rv->refs != NULL) return BFALSE;

	// We use realloc judiciously to make our life less painful here, but allocate
	// a sane number of reflections beyond what is normally encountered to make it
	// less likely to reach a memory wall.
	rv->nrefs = 0;
	rv->refs = calloc(STEP_REFS, sizeof(hklF_uno));
	if (!rv->refs) return BFALSE;
	nrefalloc = STEP_REFS;

	// keep reading until we reach a line with all zeros (0 0 0 0.00 0.00),
	// or until we reach the end-of-file
	while (freadline(curline, 4096, ipf) != EOF && !end) {
		nr = 0;
		totlen = strlen(curline);
		tmpchr = curline[4]; curline[4] = 0;
		if (sscanf(curline, "%4i", hkl) && totlen > 0) {
		  curline[4] = tmpchr; tmpchr = curline[8]; curline[8] = 0;
		  if (sscanf(curline+4, "%4i", hkl+1) && totlen > 4) {
		    curline[8] = tmpchr; tmpchr = curline[12]; curline[12] = 0;
		    if (sscanf(curline+8, "%4i", hkl+2) && totlen > 8) {
			curline[12] = tmpchr; nr = 3;
			tmpchr = curline[20];
			if (sscanf(curline+12, "%8f", tmp) && totlen > 12) nr++;
			curline[20] = tmpchr;
			tmpchr = curline[28];
			if (sscanf(curline+20, "%8f", tmp+1) && nr == 4 && totlen > 20) nr++;
			curline[28] = tmpchr;
		    }
		  }
		}
		if (nr >= 3 && hkl[0] == 0 && hkl[1] == 0 && hkl[2] == 0) { 
			end = BTRUE; // but we can still process it as an entry if tmp or tmp+1 is non-zero
		} 
		if (nr >= 4) {
			if (rv->nrefs >= nrefalloc) { // realloc to get more memory
				tmpPtr = realloc(rv->refs, sizeof(hklF_uno)*(nrefalloc+STEP_REFS));
				if (!tmpPtr) { fprintf(stdout, "Out of memory on realloc in read_hklF2.\n"); free(rv->refs); rv->nrefs = 0; return BFALSE; }
				rv->refs = tmpPtr; nrefalloc += STEP_REFS;
			}
			if (rv->nrefs >= nrefalloc) { fprintf(stdout, "Seriour error in read_hklF2. Realloc did not return enough memory.\n"); free(rv->refs); rv->nrefs = 0; return BFALSE; }
			rv->refs[rv->nrefs].hkl[0] = hkl[0];
                        rv->refs[rv->nrefs].hkl[1] = hkl[1];
                        rv->refs[rv->nrefs].hkl[2] = hkl[2];
			rv->refs[rv->nrefs].F2 = tmp[0];
			if (nr >= 5)
				rv->refs[rv->nrefs].sigF2 = tmp[1];
			else
				rv->refs[rv->nrefs].sigF2 = sqrt(fabs(tmp[0]));
			// Can't actually know F (only F2), so just set F = sqrt(|F|^2)
			if (tmp[0] > 0.0)
				rv->refs[rv->nrefs].F[0] = sqrt(tmp[0]);
			else
				rv->refs[rv->nrefs].F[0] = 0.0;
			rv->refs[rv->nrefs].F[1] = 0.0;
			if (!end || (end && fabs(tmp[0]) > 0.0))
				rv->nrefs += 1;
		}
	}

	// Do final realloc to save memory
	if (rv->nrefs < nrefalloc) {
		tmpPtr = realloc(rv->refs, sizeof(hklF_uno)*(rv->nrefs));
		if (!tmpPtr) { fprintf(stdout, "Out of memory on realloc in read_hklF2.\n"); free(rv->refs); rv->nrefs = 0; return BFALSE; }
		rv->refs = tmpPtr;
	}

	return BTRUE;
}
