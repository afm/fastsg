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

#ifndef SGFUNCTIONS_H_
#define SGFUNCTIONS_H_

namespace fsg
{
	/**
	* @class SparseGrid
	*
	* @brief Sparse grid functions
	*
	* Note: for all methods available it's up to the user to make sure he is in the [0,1]^d domain
	*
	* @author Alin Murarasu
	*
	*/
	class SparseGrid
	{
		public:
			/**
			 * Class constructor
			 * @param d Number of dimensions
			 * @param l Level of refinement
			 * @param f Function that gives the values for the sparse grid
			 */
			SparseGrid(int l, Function* f);
			virtual ~SparseGrid();

			/**
			 * @param sg The sparse grid structure
			 * @param coords The point in which we want to compute the value
			 * Evaluates (or interpolates) the sparse grid at point coords inside the [0, 1]^d domain
			 * @return The result of the evaluation
			 */
			float evaluate(float *coords);
			/**
			 * @param sg The sparse grid structure
			 * @param levels The l vector
			 * @param indices The i vector
			 * Evaluates (or interpolates) the sparse grid at point (l,i) inside the [0, 1]^d domain
			 * @return The result of the evaluation
			 */
			float evaluate(int *levels, int *indices);
			/**
			 * @param sg The sparse grid structure
			 * Computes the hierarchical coefficients for a d-dimesional, level n, non-0 boundary sparse grid.
		 	 * Initially, sg1d contains function values.
			 * @return Returns 0 if successful
			 */
			int hierarchize();
			/**
			 * @param sg The sparse grid structure
			 * @param levels The l vector of the child
			 * @param indices The i vector of the child
			 * @param plevels The l vector of the left parent
			 * @param pindices The i vector of the left parent
			 * @param cd The dimension for which we are computing the values
			 * returns the (l, i) of the left parent in dimension cd
			 * @return Returns 0 if successful
			 */
			int get_lparent(int *levels, int *indices, int *plevels, int *pindices, int cd);
			/**
			 * @param sg The sparse grid structure
			 * @param levels The l vector of the child
			 * @param indices The i vector of the child
			 * @param plevels The l vector of the right parent
			 * @param pindices The i vector of the right parent
			 * @param cd The dimension for which we are computing the values
			 * returns the (l, i) of the right parent in dimension cd
			 * @return Returns 0 if successful
			 */
			int get_rparent(int *levels, int *indices, int *plevels, int *pindices, int cd);
			/**
			 * @param sg The sparse grid structure
			 * @param coords The coords vector of the child
			 * @param pcoords The coords vector of the left parent
			 * @param cd The dimension for which we are computing the values
			 * returns the (coords) of the left parent in dimension cd
			 * @return Returns 0 if successful
			 */
			int get_lparent(float *coords, float *pcoords, int cd);
			/**
			 * @param sg The sparse grid structure
			 * @param coords The coords vector of the child
			 * @param pcoords The coords vector of the right parent
			 * @param cd The dimension for which we are computing the values
			 * returns the (coords) of the right parent in dimension cd
			 * @return Returns 0 if successful
			 */
			int get_rparent(float *coords, float *pcoords, int cd);
			/**
			 * @param sg The sparse grid structure
			 * @param crt_levels
			 * @param crt_indices
			 * @param next_levels
			 * @param next_indices
			 * Returns the (l, i) pair corresponding to the beginning of the next sparse grid
			 * @return Returns 0 if successful
			 */
			int next(int *crt_levels, int *crt_indices, int *next_levels, int *next_indices);
			int getNumOfGridPoints() const;
			/**
			 * @param d The number of dimensions
			 * @param n The level of refinement
			 * The size of a non-zero boundary, d-dimensional, n-refined sparse grid
			 * @return The size
			 */

			int size();
			
			static int size(int d, int n);
			
			int getD();
			
			int getL();
			
		private:
			int numOfGridPoints;
			float *sg1d;
			int d, l;
	};
}

#endif /* SGFUNCTIONS_H_ */
