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
			 * @param l Level of refinement
			 * @param f Function to be represented using the sparse grid technique
			 */
			SparseGrid(int l, Function* f);

			/**
			 * Class destructor
			 */
			virtual ~SparseGrid();

			/**
			 * @param coords The point at which we evaluate (interpolate) the sparse grid
			 * Evaluates (or interpolates) the sparse grid at point coords inside the [0, 1]^d domain
			 * @return The result of the evaluation
			 */
			float evaluate(float *coords);

			/**
			 * @param coords The set of points at which we evaluate (interpolate) the sparse grid
			 * @param n The size of the set
			 * @param vals The results of the evaluation
			 * Evaluates (or interpolates) the sparse grid at points stored in coords inside the [0, 1]^d domain
			 * @return Returns 0 if successfull
			 */
			int evaluate(float *coords, int n, float *vals);

			/**
			 * Computes the hierarchical coefficients for a d-dimesional, level n, non-0 boundary sparse grid.
		 	 * Initially, the sparse grid contains function values at required grid's coordinates.
			 * @return Returns 0 if successful
			 */
			int hierarchize();

			/**
			 * @param levels The l vector of the child
			 * @param indices The i vector of the child
			 * @param plevels The l vector of the left parent in dimension cd
			 * @param pindices The i vector of the left parent in dimension cd
			 * @param cd The dimension for which we are computing the neighbors
			 * returns the (l, i) representation of the left parent in dimension cd
			 * @return Returns 0 if successful
			 */
			int getLeftParent(int *levels, int *indices, int *plevels, int *pindices, int cd);

			/**
			 * @param levels The l vector of the child
			 * @param indices The i vector of the child
			 * @param plevels The l vector of the right parent in dimension cd
			 * @param pindices The i vector of the right parent in dimension cd
			 * @param cd The dimension for which we are computing the neighbors 
			 * returns the (l, i) representation of the right parent in dimension cd
			 * @return Returns 0 if successful
			 */			 
			int getRightParent(int *levels, int *indices, int *plevels, int *pindices, int cd);

			/**
			 * @param coords The coords vector of the child
			 * @param pcoords The coords vector of the left parent
			 * @param cd The dimension for which we are computing the neighbors
			 * returns the coords of the left parent in dimension cd
			 * @return Returns 0 if successful
			 */
			int getLeftParent(float *coords, float *pcoords, int cd);

			/**
			 * @param coords The coords vector of the child
			 * @param pcoords The coords vector of the right parent
			 * @param cd The dimension for which we are computing the neighbors
			 * returns the coords of the right parent in dimension cd
			 * @return Returns 0 if successful
			 */
			int getRightParent(float *coords, float *pcoords, int cd);

			/**
			 * @param crt_levels
			 * @param crt_indices
			 * @param next_levels
			 * @param next_indices
			 * Returns the (l, i) pair corresponding to the beginning of the next sparse grid
			 * @return Returns 0 if successful
			 */
			int next(int *crt_levels, int *crt_indices, int *next_levels, int *next_indices);

			/**
			 * The number of grid points composing the sparse grid
			 * @return The size of the sparse grid
			 */
			int size() const;
			
			/**
			 * @param d The number of dimensions
			 * @param n The level of refinement
			 * The size of a non-0 boundary, d-dimensional, level n sparse grid
			 * @return The size
			 */
			static int size(int d, int n);

			/**
			 * The number of dimensions of the sparse grid
			 * @return The dimensionality of the sparse grid
			 */			
			int getD();

			/**
			 * The refinement level of the sparse grid
			 * @return The refinement level of the sparse grid
			 */			
			int getL();
			
		private:
			int numOfGridPoints;
			float *sg1d;
			int d, l;
	};
}

#endif /* SGFUNCTIONS_H_ */
