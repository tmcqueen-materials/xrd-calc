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
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "defines.h"

BOOL VEQ_FRAC(vtr3 v, vtr3 v2) {
	vtr3 tmp;
	int i;
	if (!v || !v2) return BFALSE;
	tmp[0] = fabs(v[0] - v2[0]); i = (int)(tmp[0]+COORD_TOL); tmp[0] -= (double)i;
	tmp[1] = fabs(v[1] - v2[1]); i = (int)(tmp[1]+COORD_TOL); tmp[1] -= (double)i;
	tmp[2] = fabs(v[2] - v2[2]); i = (int)(tmp[2]+COORD_TOL); tmp[2] -= (double)i;
	if (fabs(tmp[0]) > COORD_TOL) return BFALSE;
        if (fabs(tmp[1]) > COORD_TOL) return BFALSE;
        if (fabs(tmp[2]) > COORD_TOL) return BFALSE;
	return BTRUE;
}

BOOL ATOMCOPY(atom *rv, atom *a) {
	if (!rv || !a) return BFALSE;
	if (!VCOPY_FZ(rv->xyz, a->xyz)) return BFALSE;
	if (!MCOPY_FZ(rv->uij, rv->uij)) return BFALSE;
	if (!SCATFCOPY(rv->f, a->f)) return BFALSE;
	strncpy(rv->label, a->label, MAX_LABEL_LEN-1);
	rv->occ = a->occ;
	return BTRUE;
}

BOOL SCATFCOPY(scatf rv, scatf f) {
	int i;
	if (!rv || !f) return BFALSE;
	for (i = 0; i < 13; i++) rv[i] = f[i];
	return BTRUE;
}

BOOL SYMMCOPY(symm *rv, symm *sym) {
	if (!rv || !sym) return BFALSE;
	if (!MCOPY_FZ(rv->r, sym->r)) return BFALSE;
	if (!VCOPY_FZ(rv->t, sym->t)) return BFALSE;
	return BTRUE;
}

BOOL HKLFUNOCOPY(hklF_uno *rv, hklF_uno *u) {
	if (!rv || !u) return BFALSE;
	rv->hkl[0] = u->hkl[0];
        rv->hkl[1] = u->hkl[1];
        rv->hkl[2] = u->hkl[2];
	rv->F[0] = u->F[0];
	rv->F[1] = u->F[1];
	rv->F2 = u->F2;
	rv->sigF2 = u->sigF2;
	return BTRUE;
}

// copy vector to vector, zeroing if below exact_tol
BOOL VCOPY_FZ(vtr3 rv, vtr3 v) {
	if (!rv || !v) return BFALSE;
	rv[0] = FZ(v[0]);
	rv[1] = FZ(v[1]);
	rv[2] = FZ(v[2]);
	return BTRUE;
}

// copy matrix to matrix, zeroing if below exact_tol
BOOL MCOPY_FZ(mtx33 rv, mtx33 m) {
	if (!rv || !m) return BFALSE;
	rv[0][0] = FZ(m[0][0]); rv[0][1] = FZ(m[0][1]); rv[0][2] = FZ(m[0][2]);
        rv[1][0] = FZ(m[1][0]); rv[1][1] = FZ(m[1][1]); rv[1][2] = FZ(m[1][2]);
        rv[2][0] = FZ(m[2][0]); rv[2][1] = FZ(m[2][1]); rv[2][2] = FZ(m[2][2]);
	return BTRUE;
}

// compute M times V
BOOL MV(vtr3 rv, mtx33 m, vtr3 v) {
	vtr3 tmp;
	if (!rv || !m || !v) return BFALSE;
	tmp[0] = m[0][0]*v[0] + m[0][1]*v[1] + m[0][2]*v[2];
	tmp[1] = m[1][0]*v[0] + m[1][1]*v[1] + m[1][2]*v[2];
	tmp[2] = m[2][0]*v[0] + m[2][1]*v[1] + m[2][2]*v[2];
	VCOPY_FZ(rv, tmp);
	return BTRUE;
}

// compute M times V plus V2
BOOL MVpV2(vtr3 rv, mtx33 m, vtr3 v, vtr3 v2) {
	vtr3 tmp;
	if (!rv || !m || !v || !v2) return BFALSE;
	if (!MV(tmp, m, v)) return BFALSE;
	rv[0] = tmp[0] + v2[0]; rv[1] = tmp[1] + v2[1]; rv[2] = tmp[2] + v2[2];
	return BTRUE;
}

// compute M-transpose
BOOL Mt(mtx33 rv, mtx33 m) {
	FP tmp;
	if (!rv || !m) return BFALSE;
	rv[0][0] = m[0][0]; rv[1][1] = m[1][1]; rv[2][2] = m[2][2];
	tmp = m[0][1]; rv[0][1] = m[1][0]; rv[1][0] = tmp;
	tmp = m[0][2]; rv[0][2] = m[2][0]; rv[2][0] = tmp;
	tmp = m[1][2]; rv[1][2] = m[2][1]; rv[2][1] = tmp;
	return BTRUE;
}

// compute M-inverse
BOOL Minv(mtx33 rv, mtx33 m) {
	mtx33 tmp;
	FP tmp2;
	if (!rv || !m) return BFALSE;
	tmp[0][0] = m[2][2]*m[1][1]-m[2][1]*m[1][2];
        tmp[1][0] = m[1][0]*m[2][2]-m[2][0]*m[1][2];
        tmp[2][0] = m[1][0]*m[2][1]-m[2][0]*m[1][1];
        tmp2 = FZ(m[0][0]*tmp[0][0]-m[0][1]*tmp[1][0]+m[0][2]*tmp[2][0]);
        if (tmp2 <= 0.0) return BFALSE;
        tmp[0][0] = (m[1][1]*m[2][2]-m[2][1]*m[1][2])/tmp2;
        tmp[0][1] = (m[0][2]*m[2][1]-m[2][2]*m[0][1])/tmp2;
        tmp[0][2] = (m[0][1]*m[1][2]-m[1][1]*m[0][2])/tmp2;
        tmp[1][0] = (m[1][2]*m[2][0]-m[2][2]*m[1][0])/tmp2;
        tmp[1][1] = (m[0][0]*m[2][2]-m[2][0]*m[0][2])/tmp2;
        tmp[1][2] = (m[0][2]*m[1][0]-m[1][2]*m[0][0])/tmp2;
        tmp[2][0] = (m[1][0]*m[2][1]-m[2][0]*m[1][1])/tmp2;
        tmp[2][1] = (m[0][1]*m[2][0]-m[2][1]*m[0][0])/tmp2;
        tmp[2][2] = (m[0][0]*m[1][1]-m[1][0]*m[0][1])/tmp2;
        if (!MCOPY_FZ(rv,tmp)) return BFALSE;
        return BTRUE;
}

// compute M-transpose times (V times M)
BOOL MtVM(vtr3 rv, vtr3 v, mtx33 m) {
	vtr3 tmp;
	mtx33 tmp2;
	if (!rv || !v || !m) return BFALSE;
	tmp[0] = (v[0]*m[0][0]+v[1]*m[1][0]+v[2]*m[2][0]);
        tmp[1] = (v[0]*m[0][1]+v[1]*m[1][1]+v[2]*m[2][1]);
        tmp[2] = (v[0]*m[0][2]+v[1]*m[1][2]+v[2]*m[2][2]);
	if (!Mt(tmp2, m)) return BFALSE;
	rv[0] = tmp2[0][0]*tmp[0]+tmp2[0][1]*tmp[1]+tmp2[0][2]*tmp[2];
	rv[1] = tmp2[1][0]*tmp[0]+tmp2[1][1]*tmp[1]+tmp2[1][2]*tmp[2];
	rv[2] = tmp2[2][0]*tmp[0]+tmp2[2][1]*tmp[1]+tmp2[2][2]*tmp[2];
	return BTRUE;
}

// compute M2 times M
BOOL M2M(mtx33 rv, mtx33 m2, mtx33 m) {
	mtx33 tmp;
	if (!rv || !m2 || !m) return BFALSE;
	// M2 times M
        tmp[0][0] = m2[0][0]*m[0][0]+m2[0][1]*m[1][0]+m2[0][2]*m[2][0];
        tmp[0][1] = m2[0][0]*m[0][1]+m2[0][1]*m[1][1]+m2[0][2]*m[2][1];
        tmp[0][2] = m2[0][0]*m[0][2]+m2[0][1]*m[1][2]+m2[0][2]*m[2][2];
        tmp[1][0] = m2[1][0]*m[0][0]+m2[1][1]*m[1][0]+m2[1][2]*m[2][0];
        tmp[1][1] = m2[1][0]*m[0][1]+m2[1][1]*m[1][1]+m2[1][2]*m[2][1];
        tmp[1][2] = m2[1][0]*m[0][2]+m2[1][1]*m[1][2]+m2[1][2]*m[2][2];
        tmp[2][0] = m2[2][0]*m[0][0]+m2[2][1]*m[1][0]+m2[2][2]*m[2][0];
        tmp[2][1] = m2[2][0]*m[0][1]+m2[2][1]*m[1][1]+m2[2][2]*m[2][1];
        tmp[2][2] = m2[2][0]*m[0][2]+m2[2][1]*m[1][2]+m2[2][2]*m[2][2];
 	if (!MCOPY_FZ(rv, tmp)) return BFALSE;
	return BTRUE;
}

// compute M3 times (M2 times M)
BOOL M3M2M(mtx33 rv, mtx33 m3, mtx33 m2, mtx33 m) {
	mtx33 tmp;
	if (!rv || !m3 || !m2 || !m) return BFALSE;
	// M2 times M
	if (!M2M(tmp, m2, m)) return BFALSE;
	// M3 times above
	if (!M2M(rv, m3, tmp)) return BFALSE;
	return BTRUE;
}

// f0 + f' + i*f''
// f0 = c + sum(0to4, ai*exp(-bi*sin(Th)^2/lam^2))
BOOL fcoeff(CMPLX rv, scatf f, FP sin2ThInvLam2) {
	int i;
	if (!f || !rv) return BFALSE;
	rv[0] = f[11]; rv[1] = f[12];
	rv[0] += f[10];
	for (i = 0; i < 5; i++) 
		rv[0] += f[2*i]*exp(-1.0*f[2*i+1]*sin2ThInvLam2);
	return BTRUE;
}

// FH = (f0+f'+i*f'')*exp(2*pi*i*H*X - 2*pi^2*H*U*H)
BOOL FHatom(CMPLX rv, vtr3 hkl, int natom, unitcell *cell) {
	FP tmp, sin2ThInvLam2;
	CMPLX f, ept;
	mtx33 invuij;
	atom *a;
        int i,j;
	if (!cell || !rv || !hkl || natom < 0 || natom >= cell->natom || !cell->a) return BFALSE;
	a = &(cell->a[natom]);
	if (!sinThInvLam(&tmp, hkl, cell)) return BFALSE;
	sin2ThInvLam2 = tmp*tmp;
	if (!fcoeff(f, a->f, sin2ThInvLam2)) return BFALSE;
	tmp = p2i*DOT(hkl,a->xyz);
	ept[0] = cos(tmp);
	ept[1] = sin(tmp);

        // Eq. 2.1.24 from http://ww1.iucr.org/comm/cnom/adp/finrepone/finrepone.html
        tmp = 0.0;
	for (i=0;i<3;i++) for (j=0;j<3;j++) tmp += hkl[i]*cell->invabc[i]*a->uij[i][j]*cell->invabc[j]*hkl[j];
        tmp = exp(-1.0*p2i2*tmp);

	ept[0] *= tmp;
	ept[1] *= tmp;
	cmplx_mul(rv, f, ept);
	return BTRUE;
}

BOOL cmplx_mul(CMPLX rv, CMPLX a, CMPLX b) {
	CMPLX tmp;
	if (!rv || !a || !b) return BFALSE;
	tmp[0] = CMR(a,b);
	tmp[1] = CMI(a,b);
	rv[0] = FZ(tmp[0]);
	rv[1] = FZ(tmp[1]);
	return BTRUE;
}

// Use formula from international tables, vol. C
BOOL sinThInvLam(FP *rv, vtr3 hkl, unitcell *cell) {
	FP tmp = 0.0;
	if (!rv || !hkl || !cell) return BFALSE;
	tmp = hkl[0]*hkl[0]*cell->invabc[0]*cell->invabc[0]; // h^2 * (a*)^2
	tmp += hkl[1]*hkl[1]*cell->invabc[1]*cell->invabc[1]; // + k^2 * (b*)^2
	tmp += hkl[2]*hkl[2]*cell->invabc[2]*cell->invabc[2]; // + l^2 * (c*)^2
	tmp += 2.0*hkl[0]*hkl[1]*cell->invabc[0]*cell->invabc[1]*cos(cell->invalbega[2]); // + 2hk(a*)(b*)cos(gamma*)
        tmp += 2.0*hkl[1]*hkl[2]*cell->invabc[1]*cell->invabc[2]*cos(cell->invalbega[0]); // + 2kl(b*)(c*)cos(alpha*)
        tmp += 2.0*hkl[2]*hkl[0]*cell->invabc[2]*cell->invabc[0]*cos(cell->invalbega[1]); // + 2lh(c*)(a*)cos(beta*)
        rv[0] = sqrt(tmp)/2.0;
	return BTRUE;
}

// Formulas from international tables, vol. C
BOOL cell_fill(unitcell *cell) {
	FP tmp;
	if (!cell) return BFALSE;

	// Volume first
	tmp = 1.0-cos(cell->albega[0])*cos(cell->albega[0]);
	tmp -= cos(cell->albega[1])*cos(cell->albega[1]);
	tmp -= cos(cell->albega[2])*cos(cell->albega[2]);
	tmp += 2.0*cos(cell->albega[0])*cos(cell->albega[1])*cos(cell->albega[2]);
	cell->vol = cell->abc[0]*cell->abc[1]*cell->abc[2]*sqrt(tmp);
	cell->invvol = 1.0/cell->vol;

	// a*, b*, c*
	cell->invabc[0] = cell->invvol*cell->abc[1]*cell->abc[2]*sin(cell->albega[0]);
	cell->invabc[1] = cell->invvol*cell->abc[0]*cell->abc[2]*sin(cell->albega[1]);
	cell->invabc[2] = cell->invvol*cell->abc[0]*cell->abc[1]*sin(cell->albega[2]);

	// alpha*, beta*, gamma* are slightly harder, due to arcsin being defined over -pi/2,pi/2
	// and arccos over 0,pi, but angles conventionally [0,2pi).
	// So first calculate from inverse cosine, then check and possibly correct
	// with inverse sine.
        tmp = cos(cell->albega[1])*cos(cell->albega[2])-cos(cell->albega[0]);
        cell->invalbega[0] = acos(tmp/sin(cell->albega[1])/sin(cell->albega[2]));
	tmp = cell->abc[0]*cell->abc[1]*cell->abc[2]*sin(cell->albega[1])*sin(cell->albega[2]);
	if (tmp+TRIG_TOL < 0.0) cell->invalbega[0] = p2i - cell->invalbega[0];
        tmp = cos(cell->albega[0])*cos(cell->albega[2])-cos(cell->albega[1]);
        cell->invalbega[1] = acos(tmp/sin(cell->albega[0])/sin(cell->albega[2]));
        tmp = cell->abc[0]*cell->abc[1]*cell->abc[2]*sin(cell->albega[0])*sin(cell->albega[2]);
        if (tmp+TRIG_TOL < 0.0) cell->invalbega[1] = p2i - cell->invalbega[1];
        tmp = cos(cell->albega[0])*cos(cell->albega[1])-cos(cell->albega[2]);
        cell->invalbega[2] = acos(tmp/sin(cell->albega[0])/sin(cell->albega[1]));
        tmp = cell->abc[0]*cell->abc[1]*cell->abc[2]*sin(cell->albega[0])*sin(cell->albega[1]);
        if (tmp+TRIG_TOL < 0.0) cell->invalbega[2] = p2i - cell->invalbega[2];

	// Calculate the inverse cell volume from the reciprocal cell
	tmp = 1.0-cos(cell->invalbega[0])*cos(cell->invalbega[0]);
	tmp -= cos(cell->invalbega[1])*cos(cell->invalbega[1]);
	tmp -= cos(cell->invalbega[2])*cos(cell->invalbega[2]);
	tmp += 2.0*cos(cell->invalbega[0])*cos(cell->invalbega[1])*cos(cell->invalbega[2]);
	cell->invvol = cell->invabc[0]*cell->invabc[1]*cell->invabc[2]*sqrt(tmp);
	// Sanity Check
	if (fabs(cell->invvol-1.0/cell->vol) > TRIG_TOL) {
		fprintf(stderr, "cell_fill: Oops! Inverse volume of real space cell and calculated volume of inverse cell do not match!\n");
		print_cell(cell);
		return BFALSE;
	}

	// Calculate G* metric tensor
	cell->invG[0][0] = FZ(cell->invabc[0]*cell->invabc[0]);
	cell->invG[0][1] = FZ(cell->invabc[0]*cell->invabc[1]*cos(cell->invalbega[2]));
	cell->invG[0][2] = FZ(cell->invabc[0]*cell->invabc[2]*cos(cell->invalbega[1]));
	cell->invG[1][0] = cell->invG[0][1];
	cell->invG[1][1] = FZ(cell->invabc[1]*cell->invabc[1]);
	cell->invG[1][2] = FZ(cell->invabc[1]*cell->invabc[2]*cos(cell->invalbega[0]));
	cell->invG[2][0] = cell->invG[0][2];
	cell->invG[2][1] = cell->invG[1][2];
	cell->invG[2][2] = FZ(cell->invabc[2]*cell->invabc[2]);

	// And A* metric tensor
	cell->invA[0][0] = FZ(cell->invabc[0]);
	cell->invA[0][1] = FZ(cell->invabc[1]*cos(cell->invalbega[2]));
	cell->invA[0][2] = FZ(cell->invabc[2]*cos(cell->invalbega[1]));
	cell->invA[1][0] = 0.0;
	cell->invA[1][1] = FZ(sqrt(cell->invabc[1]*cell->invabc[1]-cell->invA[0][1]*cell->invA[0][1]));
	cell->invA[1][2] = FZ((cell->invabc[1]*cell->invabc[2]*cos(cell->invalbega[0])-cell->invA[0][1]*cell->invA[0][2])/cell->invA[1][1]);
	cell->invA[2][0] = 0.0;
	cell->invA[2][1] = 0.0;
	cell->invA[2][2] = FZ(sqrt(cell->invabc[2]*cell->invabc[2]-cell->invA[0][2]*cell->invA[0][2]-cell->invA[1][2]*cell->invA[1][2]));

	// We do not deal with the atom or symmetry parts

	return BTRUE;
}

// Usual definition, from Int. Tables Vol. C
BOOL hklF_fill(hklF *rv, unitcell *cell) {
	int i, j;
	CMPLX tmp;
	if (!rv || !rv->refs || !cell || !cell->a || cell->natom < 1 || cell->natom > MAX_ATOMS) return BFALSE;
	if (rv->nrefs < 1 || rv->nrefs > MAX_REFS) return BFALSE;

	for (i = 0; i < rv->nrefs; i++) {
		rv->refs[i].F[0] = 0.0;
		rv->refs[i].F[1] = 0.0;
		for (j = 0; j < cell->natom; j++) {
			if (!FHatom(tmp, rv->refs[i].hkl, j, cell)) return BFALSE;
			rv->refs[i].F[0] += tmp[0]; rv->refs[i].F[1] += tmp[1];
		}
		rv->refs[i].F2 = CMMAG(rv->refs[i].F)*CMMAG(rv->refs[i].F);
		rv->refs[i].sigF2 = CMMAG(rv->refs[i].F); // = sqrt(I)
	}

	return BTRUE;
}

// For two sets of reflection values determine an overall scale
// factor, OSF, such that sum(w*(|F1|^2-OSF*(|F2|^2))^2)^(1/2)
// with w = 1/(sig|F1|^2 + (aP)^2 + bP) and
// P = (2*OSF*|F2|^2 + Max(|F1|^2,0)) / 3,
// is minimized.
// This is just an unreduced Goodness of Fit, and
// can be used to calculate GooF and wR2.
// This assumes that the number of reflections in F1 and F2 are
// the same, and that they are in the same order!
BOOL find_osf(FP *osf, FP *GooF, hklF *F1, hklF *F2, FP a, FP b, FP tol) {
	int cref;
	FP lastGooF = 1e200, curGooF, curOSF, stepOSF;
	FP P, w, tmp;

	if (!osf || !GooF || !F1 || !F2 || fabs(a) < TRIG_TOL || b < 0.0 || F1->nrefs != F2->nrefs || tol < EXACT_TOL) return BFALSE;

	// Go through the reflections once to make sure that all
	// hkl's are the same, and determine an initial overall
	// scale factor guess based on sum(|F2|^2)/sum(|F1|^2)
	P = 0.0; w = 0.0;
	for (cref = 0; cref < F1->nrefs; cref++) {
                if ((int)(F1->refs[cref].hkl[0])-(int)(F2->refs[cref].hkl[0]) != 0 ||
                    (int)(F1->refs[cref].hkl[1])-(int)(F2->refs[cref].hkl[1]) != 0 ||
                    (int)(F1->refs[cref].hkl[2])-(int)(F2->refs[cref].hkl[2]) != 0)
                        return BFALSE; // reflections not the same
		P += F1->refs[cref].F2;
		w += F2->refs[cref].F2;
	}
	if (fabs(P) > 0.0) curOSF = w/P;
	else curOSF = tol;
	stepOSF = 0.05*curOSF;
	if (fabs(stepOSF) <= 2.0*tol) stepOSF = 4.0*tol;

	while (fabs(stepOSF/curOSF) >= tol) {
		curGooF = 0.0;
		// Iterate over all reflections and calculate new GooF
		for (cref = 0; cref < F1->nrefs; cref++) {
                        P = 2.0*curOSF*F2->refs[cref].F2; if (F1->refs[cref].F2 > 0.0) P += F1->refs[cref].F2;  
			P /= 3.0;
                        w = a*a*P*P+b*P + fabs(F1->refs[cref].sigF2);
                        if (w <= 0.0) w = 0.0;
                        else w = 1.0/w;
                        tmp = F1->refs[cref].F2 - curOSF*F2->refs[cref].F2;
                        curGooF += w*tmp*tmp;
		}
		if (curGooF < 0.0) return BFALSE; // should never happen
		curGooF = sqrt(curGooF);
		if (curGooF < lastGooF) { curOSF += stepOSF; stepOSF *= 1.05; }
		else { stepOSF /= -2.0; curOSF += stepOSF; }
		while (curOSF <= 0.0) { stepOSF = fabs(stepOSF); curOSF += stepOSF; }
		lastGooF = curGooF;
	}
	*osf = curOSF;
	*GooF = curGooF;
	return BTRUE;
}

// Merge all reflections with identical h,k,l indices into a single
// value. Works in place, so if it returns with an error, the
// state of F is undefined.
BOOL hklF_merge_identical(hklF *rv) {
	int i, j, nr;
	BOOL *used;
        hklF_uno *tmp;

	if (!rv || !rv->refs || rv->nrefs < 1 || rv->nrefs > MAX_REFS) return BFALSE;

	used = calloc(rv->nrefs, sizeof(BOOL));
	if (!used) return BFALSE;

	// Loop over all reflections
	for (i = 0; i < rv->nrefs; i++) {
		if (!used[i]) {
			// New reflection
			// We DO NOT mark these as used. All unused reflections
			// at the end are the ones to keep!
			nr = 1;
			rv->refs[i].sigF2 = rv->refs[i].sigF2*rv->refs[i].sigF2;
			// Find all identical ones and add
			for (j = i+1; j < rv->nrefs; j++)
				if ((int)(rv->refs[i].hkl[0]) == (int)(rv->refs[j].hkl[0]) &&
				    (int)(rv->refs[i].hkl[1]) == (int)(rv->refs[j].hkl[1]) &&
				    (int)(rv->refs[i].hkl[2]) == (int)(rv->refs[j].hkl[2])) {
					used[j] = BTRUE;
					rv->refs[i].F2 += rv->refs[j].F2;
					rv->refs[i].sigF2 += (rv->refs[j].sigF2*rv->refs[j].sigF2);
					nr++;
				}
			// Calculate final numbers, adding error in quadrature
			rv->refs[i].F2 /= (FP)nr;
			rv->refs[i].sigF2 = sqrt(rv->refs[i].sigF2)/(FP)nr;
			if (rv->refs[i].F2 > 0.0) rv->refs[i].F[0] = sqrt(rv->refs[i].F2);
			else rv->refs[i].F[0] = 0.0;
			rv->refs[i].F[1] = 0.0;
		}
	}

	// Condense and shrink
	nr = 0;
	for (i = 0; i < rv->nrefs; i++) {
		if (!used[i]) {
			if (!HKLFUNOCOPY(&(rv->refs[nr]), &(rv->refs[i]))) { free(used); return BFALSE; }
			nr++;
		}
	}
	rv->nrefs = nr;

	free(used);

	// Reallocate the array. On failure, ignore, as having a too-big allocation is not fatal
        tmp = NULL;
	tmp = realloc(rv->refs, sizeof(hklF_uno)*(rv->nrefs));
	if (tmp != NULL) { rv->refs = tmp; }

	// Done
	return BTRUE;
}
