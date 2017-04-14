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

char *pass = "Pass";
char *fail = "Fail";

char *test_answer(BOOL tv) {
	if (tv) return pass;
	return fail;
}

// Test 1: Complex functions (multiplication, mag, phase)
BOOL test01() {
	CMPLX a, b, c;
	a[0] = 0.707; a[1] = -1.324;
	b[0] = 1.123; b[1] = 0.123;
	if (!cmplx_mul(c,a,b)) return BFALSE;
	fprintf(stdout, "Test 01 Result: %.14f + %.14f*i\n", c[0], c[1]);
	if (fabs(c[0]-0.956813)>EXACT_TOL) return BFALSE;
	if (fabs(c[1]+1.399891)>EXACT_TOL) return BFALSE;
	b[0] = -0.707; b[1] = -1.324;
	if (fabs(CMMAG(a)-CMMAG(b))>EXACT_TOL) return BFALSE;
	if (fabs(540.0-CMPH(a)-CMPH(b))>TRIG_TOL) return BFALSE;
	a[0] = 1.123; a[1] = 0.123;
	b[0] = -1.123; b[1] = 0.123;
        if (fabs(CMMAG(a)-CMMAG(b))>EXACT_TOL) return BFALSE;
        if (fabs(180.0-CMPH(a)-CMPH(b))>TRIG_TOL) return BFALSE;
	return BTRUE;
}

// Test 2: Vector-Vector and Vector-Matrix Routines
BOOL test02() {
	vtr3 a, b;
	mtx33 m;
	FP rv;
	a[0] = 1.65; a[1] = -6.57; a[2] = 9.88;
	b[0] = -100.0; b[1] = 10.0; b[2] = 0.01;
	m[0][0] = 1.0; m[0][1] = 0.1; m[0][2] = 3.0;
	m[1][0] = 0.2; m[1][1] = 5.0; m[1][2] = -6.0;
	m[2][0] = -8.0; m[2][1] = 0.0; m[2][2] = 3.5;
	rv = DOT(a,b);
	fprintf(stdout, "Test 02 Result: DOT=%.6e", rv);
	if (fabs(rv+230.6012) > MATRIX_TOL) return BFALSE;
	rv = VDOTMV(a,m,a);
	fprintf(stdout, ", VDOTMV=%.7e", rv);
	if (fabs(rv-864.90485) > MATRIX_TOL) return BFALSE;
	if (!MtVM(a, b, m)) return BFALSE;
	fprintf(stdout, ", MtVM=[%.6e %.6e %.6e]\n", a[0], a[1], a[2]);
	if (fabs(a[0]-2789.64) > MATRIX_TOL || fabs(a[1]-190.192) > MATRIX_TOL || fabs(a[2]+1794.1175) > MATRIX_TOL) return BFALSE;
	return BTRUE;
}

// Test 3: Matrix-Matrix Routines
BOOL test03() {
	mtx33 m3, m2, m, rv;
        m[0][0] = 1.0; m[0][1] = 0.1; m[0][2] = 3.0;
        m[1][0] = 0.2; m[1][1] = 5.0; m[1][2] = -6.0;
        m[2][0] = -8.0; m[2][1] = 0.0; m[2][2] = 3.5;
        m2[0][0] = 10.0; m2[0][1] = 0.01; m2[0][2] = 30.0;
        m2[1][0] = 0.02; m2[1][1] = 50.0; m2[1][2] = -60.0;
        m2[2][0] = -80.0; m2[2][1] = 0.0; m2[2][2] = 35.0;
        m3[0][0] = 1.0; m3[0][1] = 1.0; m3[0][2] = 3.0;
        m3[1][0] = -1.0; m3[1][1] = 2.0; m3[1][2] = -2.0;
        m3[2][0] = -2.0; m3[2][1] = 2.0; m3[2][2] = 3.0;
	fprintf(stdout, "Test 03 Result:\n");
	if (!M2M(rv, m2, m)) return BFALSE;
	fprintf(stdout, " M2M=\n");
	print_matrix(rv);
	if (fabs(rv[0][0]+229.998) > MATRIX_TOL || fabs(rv[0][1]-1.05) > MATRIX_TOL || fabs(rv[0][2]-134.94) > MATRIX_TOL) return BFALSE;
        if (fabs(rv[1][0]-490.02) > MATRIX_TOL || fabs(rv[1][1]-250.002) > MATRIX_TOL || fabs(rv[1][2]+509.94) > MATRIX_TOL) return BFALSE;
        if (fabs(rv[2][0]+360.0) > MATRIX_TOL || fabs(rv[2][1]+8.0) > MATRIX_TOL || fabs(rv[2][2]+117.5) > MATRIX_TOL) return BFALSE;

	if (!M3M2M(rv, m2, m3, m)) return BFALSE;
	fprintf(stdout, " M3M2M=\n");
	print_matrix(rv);
        if (fabs(rv[0][0]+995.846) > MATRIX_TOL || fabs(rv[0][1]-345.099) > MATRIX_TOL || fabs(rv[0][2]+150.22) > MATRIX_TOL) return BFALSE;
        if (fabs(rv[1][0]-2305.544) > MATRIX_TOL || fabs(rv[1][1]+92.898) > MATRIX_TOL || fabs(rv[1][2]+649.85) > MATRIX_TOL) return BFALSE;
        if (fabs(rv[2][0]-928.0) > MATRIX_TOL || fabs(rv[2][1]+65.0) > MATRIX_TOL || fabs(rv[2][2]+862.5) > MATRIX_TOL) return BFALSE;

        if (!Minv(m3, m2)) return BFALSE;
        fprintf(stdout, " Minv=\n");
        print_matrix(m3);
        if (!M2M(rv, m3, m2)) return BFALSE;
        fprintf(stdout, " Minv*M=\n");
        print_matrix(rv);
        if (fabs(rv[0][0]-1.0) > MATRIX_TOL || fabs(rv[0][1]) > MATRIX_TOL || fabs(rv[0][2]) > MATRIX_TOL) return BFALSE;
        if (fabs(rv[1][1]-1.0) > MATRIX_TOL || fabs(rv[1][0]) > MATRIX_TOL || fabs(rv[1][2]) > MATRIX_TOL) return BFALSE;
        if (fabs(rv[2][2]-1.0) > MATRIX_TOL || fabs(rv[2][1]) > MATRIX_TOL || fabs(rv[2][0]) > MATRIX_TOL) return BFALSE; 

	return BTRUE;
}

// Test 4: Unit cell fill works. One (odd) triclinic cell is sufficient
// Cell is for Re2O5Cl2, and 'correct' values from output of gsas and sexie
BOOL test04() {
        unitcell *cell;
	mtx33 m;
        cell = calloc(1, sizeof(unitcell));
        if (!cell) return BFALSE;
        cell->abc[0] = 5.78; cell->abc[1] = 6.02; cell->abc[2] = 8.08;
        cell->albega[0] = 114.8*toRad; cell->albega[1] = 96.2*toRad; cell->albega[2] = 95.0*toRad;
        if (!cell_fill(cell)) { free(cell); return BFALSE; }
        fprintf(stdout, "Test 04 Result:\n");
        print_cell(cell);
        if (fabs(cell->invabc[0]-0.1759336) > LEN_TOL) { free(cell); return BFALSE; }
        if (fabs(cell->invabc[1]-0.1849920) > LEN_TOL) { free(cell); return BFALSE; }
        if (fabs(cell->invabc[2]-0.1381115) > LEN_TOL) { free(cell); return BFALSE; }
        if (fabs(cell->invalbega[0]-64.33962*toRad) > ANG_TOL) { free(cell); return BFALSE; }
        if (fabs(cell->invalbega[1]-80.80175*toRad) > ANG_TOL) { free(cell); return BFALSE; }
        if (fabs(cell->invalbega[2]-81.56014*toRad) > ANG_TOL) { free(cell); return BFALSE; }
        if (fabs(cell->vol-250.97969) > VOL_TOL) { free(cell); return BFALSE; }
        if (fabs(cell->invvol-0.00398438619) > INVVOL_TOL) { free(cell); return BFALSE; }
        if (fabs(cell->invG[0][0]-3.0953e-2) > 1e-6) { free(cell); return BFALSE; }
        if (fabs(cell->invG[1][0]-4.777e-3 ) > 1e-6) { free(cell); return BFALSE; }
        if (fabs(cell->invG[1][1]-3.4223e-2) > 1e-6) { free(cell); return BFALSE; }
        if (fabs(cell->invG[2][0]-3.884e-3 ) > 1e-6) { free(cell); return BFALSE; }
        if (fabs(cell->invG[2][1]-1.1064e-2) > 1e-6) { free(cell); return BFALSE; }
        if (fabs(cell->invG[2][2]-1.9075e-2) > 1e-6) { free(cell); return BFALSE; }
        if (fabs(cell->invG[1][0]-cell->invG[0][1]) > MATRIX_TOL) { free(cell); return BFALSE; }
        if (fabs(cell->invG[2][0]-cell->invG[0][2]) > MATRIX_TOL) { free(cell); return BFALSE; }
        if (fabs(cell->invG[2][1]-cell->invG[1][2]) > MATRIX_TOL) { free(cell); return BFALSE; }
	if (!Mt(m, cell->invA)) { free(cell); return BFALSE; }
	if (!M2M(m, m, cell->invA)) { free(cell); return BFALSE; }
	if (fabs(m[0][0]-cell->invG[0][0]) > MATRIX_TOL) { free(cell); return BFALSE; }
        if (fabs(m[0][1]-cell->invG[0][1]) > MATRIX_TOL) { free(cell); return BFALSE; }
        if (fabs(m[0][2]-cell->invG[0][2]) > MATRIX_TOL) { free(cell); return BFALSE; }
        if (fabs(m[1][0]-cell->invG[1][0]) > MATRIX_TOL) { free(cell); return BFALSE; }
        if (fabs(m[1][1]-cell->invG[1][1]) > MATRIX_TOL) { free(cell); return BFALSE; }
        if (fabs(m[1][2]-cell->invG[1][2]) > MATRIX_TOL) { free(cell); return BFALSE; }
        if (fabs(m[2][0]-cell->invG[2][0]) > MATRIX_TOL) { free(cell); return BFALSE; }
        if (fabs(m[2][1]-cell->invG[2][1]) > MATRIX_TOL) { free(cell); return BFALSE; }
        if (fabs(m[2][2]-cell->invG[2][2]) > MATRIX_TOL) { free(cell); return BFALSE; }

        free(cell);
        return BTRUE;
}

// Test 5: sin(Th)/lambda calculation. Also compare to values computed using the G* matrix
BOOL test05() {
	unitcell *cell;
	vtr3 hkl;
	FP rv, rv2;

	cell = calloc(1, sizeof(unitcell));
	if (!cell) return BFALSE;
        cell->abc[0] = 5.78; cell->abc[1] = 6.02; cell->abc[2] = 8.08;
        cell->albega[0] = 114.8*toRad; cell->albega[1] = 96.2*toRad; cell->albega[2] = 95.0*toRad;
        if (!cell_fill(cell)) { free(cell); return BFALSE; }

	// 0 0 1 is 0.06906
	hkl[0] = 0; hkl[1] = 0; hkl[2] = 1;
	if (!sinThInvLam(&rv, hkl, cell)) { free(cell); return BFALSE; }
	fprintf(stdout, "Test 05 Result:\n");
	fprintf(stdout, " (001)   = %.5f", rv);
	if (fabs(rv-0.06906) > 1e-5) { free(cell); return BFALSE; }
	rv2 = sqrt(VDOTMV(hkl, cell->invG, hkl)/4.0);
	fprintf(stdout, " (%.5f)\n", rv2);
	if (fabs(rv-rv2) > MATRIX_TOL) { free(cell); return BFALSE; }
	// 0 1 0 is 0.09250
        hkl[0] = 0; hkl[1] = 1; hkl[2] = 0;
        if (!sinThInvLam(&rv, hkl, cell)) { free(cell); return BFALSE; } 
        fprintf(stdout, " (010)   = %.5f", rv);
        if (fabs(rv-0.09250) > 1e-5) { free(cell); return BFALSE; }
        rv2 = sqrt(VDOTMV(hkl, cell->invG, hkl)/4.0);
        fprintf(stdout, " (%.5f)\n", rv2);
        if (fabs(rv-rv2) > MATRIX_TOL) { free(cell); return BFALSE; }
        // 1 0 0 is 0.08797
        hkl[0] = 1; hkl[1] = 0; hkl[2] = 0;
        if (!sinThInvLam(&rv, hkl, cell)) { free(cell); return BFALSE; }
        fprintf(stdout, " (100)   = %.5f", rv);
        if (fabs(rv-0.08797) > 1e-5) { free(cell); return BFALSE; }
        rv2 = sqrt(VDOTMV(hkl, cell->invG, hkl)/4.0);
        fprintf(stdout, " (%.5f)\n", rv2);
        if (fabs(rv-rv2) > MATRIX_TOL) { free(cell); return BFALSE; }
 	// -2 3 -5 is 0.38637
        hkl[0] = -2; hkl[1] = 3; hkl[2] = -5;
        if (!sinThInvLam(&rv, hkl, cell)) { free(cell); return BFALSE; }
        fprintf(stdout, " (-23-5) = %.5f", rv);
        if (fabs(rv-0.38637) > 1e-5) { free(cell); return BFALSE; }
        rv2 = sqrt(VDOTMV(hkl, cell->invG, hkl)/4.0);
        fprintf(stdout, " (%.5f)\n", rv2);
        if (fabs(rv-rv2) > MATRIX_TOL) { free(cell); return BFALSE; }
 	// 18 15 20 is 3.05667
        hkl[0] = 18; hkl[1] = 15; hkl[2] = 20;
        if (!sinThInvLam(&rv, hkl, cell)) { free(cell); return BFALSE; }
        fprintf(stdout, " (010)   = %.5f", rv);
        if (fabs(rv-3.05667) > 1e-5) { free(cell); return BFALSE; }
        rv2 = sqrt(VDOTMV(hkl, cell->invG, hkl)/4.0);
        fprintf(stdout, " (%.5f)\n", rv2);
        if (fabs(rv-rv2) > MATRIX_TOL) { free(cell); return BFALSE; }
 
	free(cell);
	return BTRUE;
}

// Test 6: Test FHatom and fcoeff by computing F values for simple cubic
// polonium and comparing to that obtained from powder prediction software
// First, Diamond 3 (no U, f', or f''), then SHELX (with U, f', f'')
BOOL test06() {
	unitcell *cell;
	vtr3 hkl;
	CMPLX rv;

        cell = calloc(1, sizeof(unitcell));
	if (!cell) return BFALSE;
	cell->a = calloc(1, sizeof(atom));
        if (!cell->a) { free(cell); return BFALSE; }
        cell->abc[0] = 3.359; cell->abc[1] = 3.359; cell->abc[2] = 3.359;
        cell->albega[0] = 90.0*toRad; cell->albega[1] = 90.0*toRad; cell->albega[2] = 90.0*toRad;
        if (!cell_fill(cell)) { free(cell->a); free(cell); return BFALSE; }
	fprintf(stdout, "Test 06 Result:\n");

	cell->a->occ = 1.0;
	cell->a->xyz[0] = 0.0; cell->a->xyz[1] = 0.0; cell->a->xyz[2] = 0.0;
	// Neglect U's initially
	cell->a->uij[0][0] = 0.00; cell->a->uij[0][1] = 0.00; cell->a->uij[0][2] = 0.00;
	cell->a->uij[1][0] = 0.00; cell->a->uij[1][1] = 0.00; cell->a->uij[1][2] = 0.00;
	cell->a->uij[2][0] = 0.00; cell->a->uij[2][1] = 0.00; cell->a->uij[2][2] = 0.00;
	cell->natom = 1;

	// Use International Tables Vol. C four exponential versions for conformance
	// And neglect f' and f'' as powder prediction software Diamond 3 ignores those terms
	cell->a->f[0] = 34.6726; cell->a->f[1] =  0.700999; cell->a->f[2] = 15.4733; cell->a->f[3] =  3.55078;
	cell->a->f[4] = 13.1138; cell->a->f[5] =  9.55642; cell->a->f[6] =  7.02588; cell->a->f[7] = 47.0045;
	cell->a->f[8] =  0.0; cell->a->f[9] = 0.0; cell->a->f[10] = 13.6770;
	cell->a->f[11] = 0.0; cell->a->f[12] = 0.0;

	print_cell(cell);

	// 1 0 0
	hkl[0] = 1; hkl[1] = 0; hkl[2] = 0;
	if (!FHatom(rv, hkl, 0, cell)) { free(cell->a); free(cell); return BFALSE; }
	fprintf(stdout, "F(100)       = %14.5f %9.5f\n", CMMAG(rv), CMPH(rv)*fromRad);
	if (fabs(CMMAG(rv)-75.21) > 1e-2 || fabs(CMPH(rv)) > TRIG_TOL) { free(cell->a); free(cell); return BFALSE; }

        // 0 0 -2
        hkl[0] = 0; hkl[1] = 0; hkl[2] = -2;
        if (!FHatom(rv, hkl, 0, cell)) { free(cell->a); free(cell); return BFALSE; }
        fprintf(stdout, "F(00-2)      = %14.5f %9.5f\n", CMMAG(rv), CMPH(rv)*fromRad);
        if (fabs(CMMAG(rv)-63.29) > 1e-2 || fabs(CMPH(rv)) > TRIG_TOL) { free(cell->a); free(cell); return BFALSE; }

        // 1 1 1
        hkl[0] = 1; hkl[1] = 1; hkl[2] = 1;
        if (!FHatom(rv, hkl, 0, cell)) { free(cell->a); free(cell); return BFALSE; }
        fprintf(stdout, "F(111)       = %14.5f %9.5f\n", CMMAG(rv), CMPH(rv)*fromRad);
        if (fabs(CMMAG(rv)-66.25) > 1e-2 || fabs(CMPH(rv)) > TRIG_TOL) { free(cell->a); free(cell); return BFALSE; }

        // 2 3 5
        hkl[0] = 2; hkl[1] = 3; hkl[2] = 5;
        if (!FHatom(rv, hkl, 0, cell)) { free(cell->a); free(cell); return BFALSE; }
        fprintf(stdout, "F(235)       = %14.5f %9.5f\n", CMMAG(rv), CMPH(rv)*fromRad);
        if (fabs(CMMAG(rv)-33.68) > 1e-2 || fabs(CMPH(rv)) > TRIG_TOL) { free(cell->a); free(cell); return BFALSE; }

	fprintf(stdout, "With U, f', and f'':\n");
        // Now use a U
        cell->a->uij[0][0] = 0.01; cell->a->uij[0][1] = 0.00; cell->a->uij[0][2] = 0.00;
        cell->a->uij[1][0] = 0.00; cell->a->uij[1][1] = 0.01; cell->a->uij[1][2] = 0.00;
        cell->a->uij[2][0] = 0.00; cell->a->uij[2][1] = 0.00; cell->a->uij[2][2] = 0.01;
	// and f' and f'' from Int. Tab. Vol. C for Mo-Ka
        cell->a->f[11] = -5.121; cell->a->f[12] = 11.0496;

	// Compare to values from SHELX for a different set of reflections
        // 1 0 0
        hkl[0] = 1; hkl[1] = 0; hkl[2] = 0;
        if (!FHatom(rv, hkl, 0, cell)) { free(cell->a); free(cell); return BFALSE; }
        fprintf(stdout, "F(100)       = %14.5f %9.5f\n", CMMAG(rv), CMPH(rv)*fromRad);
        if (fabs(CMMAG(rv)-69.72) > 1e-2 || fabs(CMPH(rv)-8.96*toRad) > 1e-2*toRad) { free(cell->a); free(cell); return BFALSE; }

        // -1 0 0
        hkl[0] = -1; hkl[1] = 0; hkl[2] = 0;
        if (!FHatom(rv, hkl, 0, cell)) { free(cell->a); free(cell); return BFALSE; }
        fprintf(stdout, "F(-100)      = %14.5f %9.5f\n", CMMAG(rv), CMPH(rv)*fromRad);
        if (fabs(CMMAG(rv)-69.72) > 1e-2 || fabs(CMPH(rv)-8.96*toRad) > 1e-2*toRad) { free(cell->a); free(cell); return BFALSE; }

        // 2 3 2
        hkl[0] = 2; hkl[1] = 3; hkl[2] = 2;
        if (!FHatom(rv, hkl, 0, cell)) { free(cell->a); free(cell); return BFALSE; }
        fprintf(stdout, "F(232)       = %14.5f %9.5f\n", CMMAG(rv), CMPH(rv)*fromRad);
        if (fabs(CMMAG(rv)-30.54) > 1e-2 || fabs(CMPH(rv)-15.59*toRad) > 1e-2*toRad) { free(cell->a); free(cell); return BFALSE; }

        // -2 -3 -2
        hkl[0] = -2; hkl[1] = -3; hkl[2] = -2;
        if (!FHatom(rv, hkl, 0, cell)) { free(cell->a); free(cell); return BFALSE; }
        fprintf(stdout, "F(-2-3-2)    = %14.5f %9.5f\n", CMMAG(rv), CMPH(rv)*fromRad);
        if (fabs(CMMAG(rv)-30.54) > 1e-2 || fabs(CMPH(rv)-15.59*toRad) > 1e-2*toRad) { free(cell->a); free(cell); return BFALSE; }

        // 5 4 3
        hkl[0] = 5; hkl[1] = 4; hkl[2] = 3;
        if (!FHatom(rv, hkl, 0, cell)) { free(cell->a); free(cell); return BFALSE; }
        fprintf(stdout, "F(543)       = %14.5f %9.5f\n", CMMAG(rv), CMPH(rv)*fromRad);
        if (fabs(CMMAG(rv)-11.32) > 1e-2 || fabs(CMPH(rv)-24.01*toRad) > 1e-2*toRad) { free(cell->a); free(cell); return BFALSE; }

	fprintf(stdout, "\n");
	free(cell->a); free(cell);
	return BTRUE;
}

// Test 7: test reading cell data from an INS-format file and writing to an hklF2 file
BOOL test07() {
	unitcell *cell;
	hklF *F; int h, k, l;
	FILE *ipf;

	ipf = fopen("test07.ins", "r");
	if (!ipf) return BFALSE;
	cell = calloc(1, sizeof(unitcell));
	if (!cell) { fclose(ipf); return BFALSE; }

	fprintf(stdout, "Test 07 Result:\n");
	if (!read_ins_cell(cell, ipf)) { free(cell); fclose(ipf); return BFALSE; }
	fclose(ipf);

	print_cell(cell);

	if (cell->nsym != 36) { free(cell->sym); free(cell->a); free(cell); return BFALSE; }
	if (cell->natom != 54) { free(cell->sym); free(cell->a); free(cell); return BFALSE; }
	if (fabs(cell->abc[0]-6.819400)>COORD_TOL) { free(cell->sym); free(cell->a); free(cell); return BFALSE; }
	if (fabs(cell->abc[1]-6.819400)>COORD_TOL) { free(cell->sym); free(cell->a); free(cell); return BFALSE; }
	if (fabs(cell->abc[2]-14.024400)>COORD_TOL) { free(cell->sym); free(cell->a); free(cell); return BFALSE; }
	if (fabs(cell->albega[0]-90.00000*toRad)>TRIG_TOL) { free(cell->sym); free(cell->a); free(cell); return BFALSE; }
        if (fabs(cell->albega[1]-90.00000*toRad)>TRIG_TOL) { free(cell->sym); free(cell->a); free(cell); return BFALSE; }
        if (fabs(cell->albega[2]-120.00000*toRad)>TRIG_TOL) { free(cell->sym); free(cell->a); free(cell); return BFALSE; }

	fprintf(stdout, "INS file successfully read. Writing to HKLF2...\n");

	// Write HKLF file for double-checking by tester with a SHELX refinement
        F = calloc(1, sizeof(hklF));
        if (!F) { free(cell->sym); free(cell->a); free(cell); return BFALSE; }
        F->nrefs = 19*19*19;
        F->refs = calloc(F->nrefs, sizeof(hklF_uno));
        if (!(F->refs)) { free(cell->sym); free(cell->a); free(cell); free(F); return BFALSE; }
        for (h = -9; h <= 9; h++)
          for (k = -9; k <= 9; k++)
            for (l = -9; l <= 9; l++) {
                F->refs[19*19*(h+9)+19*(k+9)+(l+9)].hkl[0] = h;
                F->refs[19*19*(h+9)+19*(k+9)+(l+9)].hkl[1] = k;
                F->refs[19*19*(h+9)+19*(k+9)+(l+9)].hkl[2] = l;
            }

        if (!hklF_fill(F, cell)) { free(F->refs); free(F); free(cell->sym); free(cell->a); free(cell); return BFALSE; }
	ipf = fopen("test07.hkl", "w");
        print_hklF2(F, ipf);
	fclose(ipf);
	fprintf(stdout, "Success.\n");

	free(F->refs); free(F); free(cell->sym); free(cell->a); free(cell);
	return BTRUE;
}

// Test 08: Reading an hklF2 file (reads the one generated by test07)
BOOL test08() {
	FILE *ipf;
	hklF *F;
	int i;

	ipf = fopen("test07.hkl", "r");
	if (!ipf) return BFALSE;
	F = calloc(1, sizeof(hklF));
	if (!F) { fclose(ipf); return BFALSE; }

	fprintf(stdout, "Test 08 Result:\n");
	if (!read_hklF2(F, ipf)) { fclose(ipf); free(F); return BFALSE; }
	fclose(ipf);

	fprintf(stdout, "%i reflections read.\n", F->nrefs);
	if (F->nrefs != 19*19*19) { free(F->refs); free(F); return BFALSE; }

	for (i = 0; i < F->nrefs; i++) {
		if ((int)(F->refs[i].hkl[0]) == -9 && (int)(F->refs[i].hkl[1]) == -9 && (int)(F->refs[i].hkl[2]) == -3) {
			if (fabs(F->refs[i].F2-687.38) > 1e-2) { free(F->refs); free(F); return BFALSE; }
			if (fabs(F->refs[i].sigF2-26.22) > 1e-2) { free(F->refs); free(F); return BFALSE; }
			if (fabs(F->refs[i].F[0]-sqrt(687.38)) > 1e-2) { free(F->refs); free(F); return BFALSE; }
			if (fabs(F->refs[i].F[1]) > EXACT_TOL) { free(F->refs); free(F); return BFALSE; }
		}
                if ((int)(F->refs[i].hkl[0]) == -3 && (int)(F->refs[i].hkl[1]) == 2 && (int)(F->refs[i].hkl[2]) == 7) {
                        if (fabs(F->refs[i].F2-30525.69) > 1e-2) { free(F->refs); free(F); return BFALSE; }
                        if (fabs(F->refs[i].sigF2-174.72) > 1e-2) { free(F->refs); free(F); return BFALSE; }
                        if (fabs(F->refs[i].F[0]-sqrt(30525.69)) > 1e-2) { free(F->refs); free(F); return BFALSE; }
                        if (fabs(F->refs[i].F[1]) > EXACT_TOL) { free(F->refs); free(F); return BFALSE; }
                }
                if ((int)(F->refs[i].hkl[0]) == 0 && (int)(F->refs[i].hkl[1]) == 5 && (int)(F->refs[i].hkl[2]) == 2) {
                        if (fabs(F->refs[i].F2) > EXACT_TOL) { free(F->refs); free(F); return BFALSE; }
                        if (fabs(F->refs[i].sigF2) > EXACT_TOL) { free(F->refs); free(F); return BFALSE; }
                        if (fabs(F->refs[i].F[0]) > EXACT_TOL) { free(F->refs); free(F); return BFALSE; }
                        if (fabs(F->refs[i].F[1]) > EXACT_TOL) { free(F->refs); free(F); return BFALSE; }
                }
                if ((int)(F->refs[i].hkl[0]) == 4 && (int)(F->refs[i].hkl[1]) == -2 && (int)(F->refs[i].hkl[2]) == -6) {
                        if (fabs(F->refs[i].F2-46291.05) > 1e-2) { free(F->refs); free(F); return BFALSE; }
                        if (fabs(F->refs[i].sigF2-215.15) > 1e-2) { free(F->refs); free(F); return BFALSE; }
                        if (fabs(F->refs[i].F[0]-sqrt(46291.05)) > 1e-2) { free(F->refs); free(F); return BFALSE; }
                        if (fabs(F->refs[i].F[1]) > EXACT_TOL) { free(F->refs); free(F); return BFALSE; }
                }
                if ((int)(F->refs[i].hkl[0]) == 0 && (int)(F->refs[i].hkl[1]) == 0 && (int)(F->refs[i].hkl[2]) == 0) {
                        if (fabs(F->refs[i].F2-382141.00) > 1e-2) { free(F->refs); free(F); return BFALSE; }
                        if (fabs(F->refs[i].sigF2-618.00) > 1e-2) { free(F->refs); free(F); return BFALSE; }
                        if (fabs(F->refs[i].F[0]-sqrt(382141.00)) > 1e-2) { free(F->refs); free(F); return BFALSE; }
                        if (fabs(F->refs[i].F[1]) > EXACT_TOL) { free(F->refs); free(F); return BFALSE; }
                }
	}

	fprintf(stdout, "\n");
	free(F->refs); free(F);
	return BTRUE;
}

// Routines to test all the functions we write
int main() {
	fprintf(stderr, "Test 01: %s\n", test_answer(test01()));
	fprintf(stderr, "Test 02: %s\n", test_answer(test02()));
	fprintf(stderr, "Test 03: %s\n", test_answer(test03()));
	fprintf(stderr, "Test 04: %s\n", test_answer(test04()));
	fprintf(stderr, "Test 05: %s\n", test_answer(test05()));
	fprintf(stderr, "Test 06: %s\n", test_answer(test06()));
	fprintf(stderr, "Test 07: %s\n", test_answer(test07()));
	fprintf(stderr, "Test 08: %s\n", test_answer(test08()));
	return 0;
}
