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
#ifndef __DEFINES_H_
#define __DEFINES_H_

// The following is a hard-coded limit. The code could be changed
// so this was no longer a fixed limit.
#define MAX_FREE_VARS 32
// These three ARE hard-coded limits, but all exceed the maximum 
// number of possibilities, and thus should never have to change
#define MAX_LABEL_LEN 16
#define MAX_ATOM_TYPES 128
#define MAX_SYMM  512
// Neither of these is a hard limit, and are used simply for 
// sanity checking. They should be way beyond what is normally seen.
// Can set to max_int to disable the sanity checks
#define MAX_ATOMS 1000000
#define MAX_REFS  10000000

#define COORD_TOL 1e-5
#define VOL_TOL 1e-5
#define INVVOL_TOL 1e-9
#define LEN_TOL 1e-7
#define ANG_TOL 1e-5
#define FLUSH_TOL 1e-12
#define EXACT_TOL 1e-14
#define MATRIX_TOL FLUSH_TOL
#define TRIG_TOL 1e-7
#define BTRUE 1
#define BFALSE 0
#define BOOL int
#define FP double
#define p2i (2.0*3.14159265358979323846)
#define p2i2 (p2i*3.14159265358979323846)
#define toRad (p2i/360.0)
#define fromRad (1.0/toRad)

// Complex number a + b*i
typedef FP CMPLX[2]; 
// 3 vector
typedef FP vtr3[3];
// 3x3 matrix row, then column
typedef FP mtx33[3][3];
// 5-gaussian scattering factor a0,b0,a1,b1,...,a4,b4,c,f',f'' (see fcoeff in compute.c)
typedef FP scatf[13];

// This is a single symmetry element
typedef struct {
	mtx33 r;
	vtr3 t;
} symm;

// This is the structure used to compute the various
// contributions to the structure factor
typedef struct {
	FP occ; // occupancy as a fraction of the site (NOT general position)
	vtr3 xyz; // x,y, and z in fractions of the unit cell
	mtx33 uij; // the 3x3 symmetric Uij matrix 
	scatf f; // scattering factor for this atom (see above)
	char label[MAX_LABEL_LEN]; // label for this atom
} atom;

// This is all the information for a single reflection
typedef struct {
	vtr3 hkl;
	CMPLX F;
	FP F2;
	FP sigF2;
} hklF_uno;

// This is the information for a set of reflections
typedef struct {
	int nrefs;
	hklF_uno *refs;
} hklF;

// This is the information about the unit cell
typedef struct {
	vtr3 abc;    // a, b, c (angstroms)
	vtr3 albega; // alpha, beta, gamma (radians)
	vtr3 invabc; // a*, b*, c* (inverse angstroms)
	vtr3 invalbega; // alpha*, beta*, gamma* (radians)
	FP vol; // volume (angstroms^3)
	FP invvol; // inverse volume (angstroms^-3)
	mtx33 invG; // G* metric tensor (= (A*)'(A*))
	mtx33 invA; // A* metric tensor
	// symmetry specifics
	int nsym;  // = multiplicity of the general site
	symm *sym;
	// atom specifics
	int natom;
	atom *a;
} unitcell;

// vector dot vector
#define DOT(a,b) (((a)[0]*(b)[0])+((a)[1]*(b)[1])+((a)[2]*(b)[2]))
// vector dot (matrix times vector)
#define VDOTMV(a,m,b) (((a)[0]*((m)[0][0]*(b)[0]+(m)[0][1]*(b)[1]+(m)[0][2]*(b)[2]))+((a)[1]*((m)[1][0]*(b)[0]+(m)[1][1]*(b)[1]+(m)[1][2]*(b)[2]))+((a)[2]*((m)[2][0]*(b)[0]+(m)[2][1]*(b)[1]+(m)[2][2]*(b)[2])))
// complex multiply, real part
#define CMR(a,b) ((a)[0]*(b)[0]-(a)[1]*(b)[1])
// imaginary part
#define CMI(a,b) ((a)[0]*(b)[1]+(a)[1]*(b)[0])
// magnitude
#define CMMAG(a) (sqrt((a)[0]*(a)[0] + (a)[1]*(a)[1]))
// phase
#define CMPH(a) (fabs((a)[1]) < TRIG_TOL) ? ( ((a)[0] < 0.0) ? (180.0*toRad) : (0.0) ) : ( (fabs((a)[0]) < TRIG_TOL) ? ( ((a)[1] < 0.0) ? (270.0*toRad) : (90.0*toRad) ) : ( ((a)[1] < 0.0) ? (p2i-acos((a)[0]/(CMMAG(a)))) : (acos((a)[0]/(CMMAG(a)))) ) )
// Flush to zero if < FLUSH_TOL
#define FZ(a) ((fabs(a) < FLUSH_TOL) ? (0.0) : (a))

//
// FROM compute.c
//

// Compare two vectors assuming +/- integral values on any coordinate
// are equivalent (so are these the same point in a unit cell)
BOOL VEQ_FRAC(vtr3 v, vtr3 v2);

// Copy an atom
BOOL ATOMCOPY(atom *rv, atom *a);

// Copy a scattering factor
BOOL SCATFCOPY(scatf rv, scatf f);

// Copy a symmetry element
BOOL SYMMCOPY(symm *rv, symm *sym);

// Copy a hklF_uno
BOOL HKLFUNOCOPY(hklF_uno *rv, hklF_uno *u);

// copy v to rv, flushing small values to zero
BOOL VCOPY_FZ(vtr3 rv, vtr3 v);

// copy m to rv, flushing small values to zero
BOOL MCOPY_FZ(mtx33 rv, mtx33 m);

// compute M times V
BOOL MV(vtr3 rv, mtx33 m, vtr3 v);

// compute M times V plus V2
BOOL MVpV2(vtr3 rv, mtx33 m, vtr3 v, vtr3 v2);

// compute M-transpose
BOOL Mt(mtx33 rv, mtx33 m);

// compute M-transpose times (V times M)  -- small values flushed to zero
BOOL MtVM(vtr3 rv, vtr3 v, mtx33 m);

// compute M2 times M  -- small values flushed to zero
BOOL M2M(mtx33 rv, mtx33 m2, mtx33 m);

// compute M3 times (M2 times M)  -- small values flushed to zero
BOOL M3M2M(mtx33 rv, mtx33 m3, mtx33 m2, mtx33 m);

// complex multiply  -- small values flushed to zero
BOOL cmplx_mul(CMPLX rv, CMPLX a, CMPLX b);

// compute scattering f coefficient from atom scattering factor plus (sin(Th)/lam)^2
BOOL fcoeff(CMPLX rv, scatf f, FP sin2ThInvLam2);

// compute contribution to F from one atom
BOOL FHatom(CMPLX rv, vtr3 hkl, int natom, unitcell *cell);

// compute sin(Th)/lam for a given hkl and unit cell
BOOL sinThInvLam(FP *rv, vtr3 hkl, unitcell *cell);

// compute the volumes and reciprocal space vectors for a given real space unit cell
BOOL cell_fill(unitcell *cell);

// compute all the F factors over the given range
BOOL hklF_fill(hklF *rv, unitcell *cell);

// For two sets of reflection values, find a number, OSF, to
// multiply all the reflections in F2 by to minimize their difference
// to F1. This assumes that the number of reflections in F1 and F2 are
// the same, and that they are in the same order!
BOOL find_osf(FP *osf, FP *GooF, hklF *F1, hklF *F2, FP a, FP b, FP tol);

// Merge identical (i.e. same h,k,l) reflections, carrying error properly
BOOL hklF_merge_identical(hklF *rv);

//
// From print.c
//
void print_cell(unitcell *cell);

void print_matrix(mtx33 m);

void print_hklF(hklF *F, FILE *file);

void print_hklF2(hklF *F, FILE *file);

void print_hklF2_paired(hklF *F, FILE *file);

//
// From read.c
//
BOOL read_ins_cell(unitcell *cell, FILE *ipf);

BOOL read_hklF2(hklF *rv, FILE *ipf);

//
// From scatdb.c
//
BOOL scatf_fill(scatf rv, char *atom);

#endif
