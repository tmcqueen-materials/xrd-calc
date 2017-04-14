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
		rv[0] =  11.4604; rv[1] =  0.01040;
		rv[2] =  7.19640; rv[3] =  1.16620;
		rv[4] =  6.25560; rv[5] =  18.5194;
		rv[6] =  1.64550; rv[7] =  47.7784;
		rv[8] =  0.00000; rv[9] =  0.00000;
		rv[10] = -9.5574;
                rv[11] = 0.0; rv[12] = 0.0; // f' and f'' to 0
                return BTRUE;
	} else if (strncmp(atom, "CU", 2) == 0) {
		rv[0] =  13.3380; rv[1] =  3.58280;
		rv[2] =  7.16760; rv[3] =  0.24700;
		rv[4] =  5.61580; rv[5] =  11.3966;
		rv[6] =  1.67350; rv[7] =  64.8126;
		rv[8] =  0.00000; rv[9] = 0.00000;
		rv[10] = 1.19100;
                rv[11] = 0.0; rv[12] = 0.0; // f' and f'' to 0
                return BTRUE;
	} else if (strncmp(atom, "ZN", 2) == 0) {
		rv[0] =  14.0743; rv[1] =  3.26550;
		rv[2] =  7.03180; rv[3] =  0.23330;
		rv[4] =  5.16520; rv[5] =  10.3163;
		rv[6] =  2.41000; rv[7] =  58.7097;
		rv[8] =  0.00000; rv[9] =  0.00000;
		rv[10] = 1.30410;
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
/*RHF version
		rv[0] =  0.489918; rv[1] =  20.6593;
		rv[2] =  0.262003; rv[3] =  7.74039;
		rv[4] =  0.196767; rv[5] =  49.5519;
		rv[6] =  0.049879; rv[7] =  2.20159;
		rv[8] =  0.000000; rv[9] =  0.0;
		rv[10] = 0.001305;
*/
		//SDS version
                rv[0] =  0.493002; rv[1] =  10.5109;
                rv[2] =  0.322912; rv[3] =  26.1257;
                rv[4] =  0.140191; rv[5] =  3.14236;
                rv[6] =  0.040810; rv[7] =  57.7997;
                rv[8] =  0.000000; rv[9] =  0.0;
                rv[10] = 0.003038;
		rv[11] = 0.0; rv[12] = 0.0;
		return BTRUE;
	} else if (strncmp(atom, "O", 1) == 0) {
		rv[0] =  3.04850; rv[1] =  13.2771;
		rv[2] =  2.28680; rv[3] =  5.70110;
		rv[4] =  1.54630; rv[5] =  0.32390;
		rv[6] =  0.86700; rv[7] =  32.9089;
		rv[8] =  0.00000; rv[9] =  0.00000;
		rv[10] = 0.25080;
                rv[11] = 0.0; rv[12] = 0.0; // f' and f'' to 0
                return BTRUE;
	}
	return BFALSE;
}
