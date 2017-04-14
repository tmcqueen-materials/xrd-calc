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
#include <math.h>
#include <string.h>
#include "defines.h"

BOOL scatf_fill(scatf rv, char *atom) {
	if (!rv || !atom) return BFALSE;
	// All two-letter elements must come before one letter ones!
	if (strncmp(atom, "CL", 2) == 0) {
		rv[0] =  1.446071; rv[1] =  0.052357;
		rv[2] =  6.870609; rv[3] =  1.193165;
		rv[4] =  6.151801; rv[5] = 18.343416;
		rv[6] =  1.750347; rv[7] = 46.398394;
		rv[8] =  0.634168; rv[9] =  0.401005;
		rv[10] = 0.146773;
                rv[11] = 0.0; rv[12] = 0.0; // f' and f'' to 0
                return BTRUE;
	} else if (strncmp(atom, "CU", 2) == 0) {
		rv[0] = 14.014192; rv[1] =  3.738280;
		rv[2] =  4.784577; rv[3] =  0.003744;
		rv[4] =  5.056806; rv[5] = 13.034982;
		rv[6] =  1.457971; rv[7] = 72.554793;
		rv[8] =  6.932996; rv[9] =  0.265666;
		rv[10] = -3.254477;
                rv[11] = 0.0; rv[12] = 0.0; // f' and f'' to 0
                return BTRUE;
	} else if (strncmp(atom, "ZN", 2) == 0) {
		rv[0] = 14.741002; rv[1] =  3.388232;
		rv[2] =  6.907748; rv[3] =  0.243315;
		rv[4] =  4.642337; rv[5] = 11.903689;
		rv[6] =  2.191766; rv[7] = 63.312130;
		rv[8] = 38.424042; rv[9] =  0.000397;
		rv[10] = -36.915828;
                rv[11] = 0.0; rv[12] = 0.0; // f' and f'' to 0
                return BTRUE;
	} else if (strncmp(atom, "PO", 2) == 0) {
		rv[0] = 16.289164; rv[1] =  0.098121;
		rv[2] = 32.807170; rv[3] =  0.966265;
		rv[4] = 21.095164; rv[5] =  6.046622;
		rv[6] =  2.505901; rv[7] = 76.598071;
		rv[8] =  7.254589; rv[9] = 28.096128;
		rv[10] = 4.046556;
		rv[11] = 0.0; rv[12] = 0.0; // f' and f'' to 0
		return BTRUE;
	} else if (strncmp(atom, "H", 1) == 0) {
		rv[0] =  0.489918; rv[1] = 10.5109;
		rv[2] =  0.322912; rv[3] = 26.1257;
		rv[4] =  0.140191; rv[5] = 3.14236;
		rv[6] =  0.040810; rv[7] = 57.7997;
		rv[8] =  0.0;      rv[9] = 0.0;
		rv[10] = 0.003038;
		rv[11] = 0.0; rv[12] = 0.0;
		return BTRUE;
	} else if (strncmp(atom, "O", 1) == 0) {
		rv[0] =  2.960427; rv[1] = 14.182259;
		rv[2] =  2.508818; rv[3] =  5.936858;
		rv[4] =  0.637853; rv[5] =  0.112726;
		rv[6] =  0.722838; rv[7] = 34.958481;
		rv[8] =  1.142756; rv[9] =  0.390240;
		rv[10] = 0.027014;
                rv[11] = 0.0; rv[12] = 0.0; // f' and f'' to 0
                return BTRUE;
	}
	return BFALSE;
}
