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


#include "DataStructure.h"
#include "Function.h"

#ifndef HELPER_H_
#define HELPER_H_

namespace fsg {
	/**
	 * @class Helper
	 *
	 * @brief Collection of helper functions used by the sparse grid methods
	 *
	 *
	 * @author Alin Murarasu
	 *
	 */
	class Helper {
	public:
		/**
		 * @param n Number of elements in a set
		 * @param k Number of combinations
		 * @return Result of computation
		 */
		static int combi(int n, int k);
		/**
		 * Returns the number of grid points of a 0-boundary sparse grid, d-dimensional, level of refinement n
		 * @param d Number of dimensions
		 * @param n Level of refinement
		 * @return Result of computation
		 */
		static int zerob_size(int d, int n);
		/**
		 * Recursive function that generates the points on the sparse grid
		 * @param sg The sparse grid structure in which the result will be stored
		 * @param gp Vector of size d to store generated coordinates in a recursive call
		 * @param crt_d The dimension to start with (d-1)
		 * @param n Level of refinement
		 * @param f The function giving the value for the sparse grid points
		 * @return Number of points generated
		 */
		static int generate_grid_points(Sparse_grid_t sg, float* gp, int crt_d, int n, Function* f);
	};
}

#endif /* HELPER_H_ */
