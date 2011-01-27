/**********************************************************************************
 *
 * Copyright (c) 2009, 2010 Alin Murarasu
 *
 * This program is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program. If not, see <http://www.gnu.org/licenses/>.
 *
 * For any other enquiries send an email to Alin Murarasu, murarasu@in.tum.de.
 *
 * When publishing work that is based on this program please cite:
 * A. Murarasu, J. Weidendorfer, G. Buse, D. Butnaru, and D. Pflueger:
 * "Compact Data Structure and Scalable Algorithms for the Sparse Grid Technique"
 * PPoPP, Feb. 2011
 *
 *********************************************************************************/

#include <iostream>

#include "Helper.h"

#include "Converter.h"

using namespace fsg;

/* returns the number of grid points of a 0-boundary sparse grid, d-dimensional, level of refinement n */
int Helper::zerob_size(int d, int n)
{
	int j;
	int s0b = 0;

	if (d == 0)
		return 1;

	for (j = 0; j < n; j++)
		s0b += (1 << j) * combi(d - 1 + j, j);

	return s0b;
}

int Helper::combi(int n, int k)
{
	int i, c = 1;

	for (i = k + 1; i <= n; i++) {
		c *= i;
		c /= i - k;
	}

	return c;
}

int Helper::generate_grid_points(sparse_grid_t sg, float* gp, int crt_d, int n, Function* f)
{
	int i, j, count = 0;
	int levels[sg.d], indices[sg.d];
	int val;

	if (crt_d == -1) {
		Converter::coord2li(gp, levels, indices, sg.d);

		val = Converter::gp2idx(levels, indices, sg.d, sg.l);
		// fill vector
		sg.sg1d[val] = f->getValue(gp, sg.d);
		
		return 1;
	} else {
		gp[crt_d] = 0.0f;
		count += generate_grid_points(sg, gp, crt_d - 1, n, f);
		gp[crt_d] = 1.0f;
		count += generate_grid_points(sg, gp, crt_d - 1, n, f);

		/* for each level */
		for (i = 0; i < n; i++) {
			/* for each point on the level in dimension crt_d */
			for (j = 0; j < (1 << i); j++) {
				gp[crt_d] = (1.0f / (1 << i)) * (j + 0.5f);
				count += generate_grid_points(sg, gp, crt_d - 1, n - i, f);
			}
		}
	}

	return count;
}
