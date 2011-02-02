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

#ifndef COORDINATES_H_
#define COORDINATES_H_

namespace fsg
{
	/**
	* @class Converter
	*
	* @brief Useful methods for working with various representations of sparse grid points
	*
	*
	* @author Alin Murarasu
	*
	*/
	class Converter
	{
		public:
			/**
			 * @param levels The l component (of size d)
			 * @param indices The i component (of size d)
			 * @param d Number of dimensions
			 * @param n Level of refinement
			 * @return Index corresponding to the (levels, indices) pair
			 */
			static int gp2idx(int *levels, int *indices, int d, int n);

			/**
			 * @param index The index from the sparse grid that is converted to (l, i) representation
			 * @param levels The computed l component (of size d)
			 * @param indices The computed i component (of size d)
			 * @param d Number of dimensions
			 * @param n Level of refinement
			 * @return If successful, returns 0
			 */
			static int idx2gp(int index, int *levels, int *indices, int d, int n);

			/**
			 * @param coords Vector of coords to be converted into (l, i)
			 * @param levels The computed l component (of size d)
			 * @param indices The computed i component (of size d)
			 * @param d Number of dimensions
			 * @return If successful, returns 0
			 */
			static int coord2li(float *coords, int *levels, int *indices, int d);

			/**
			 * @param levels The l component (of size d)
			 * @param indices The i component (of size d)
			 * @param coords Vector of computed coords
			 * @param d Number of dimensions
			 * @return If successful, returns 0
			 */
			static int li2coord(int *levels, int *indices, float *coords, int d);

			/**
			 * @param coords The coordinates of the sparse grid point
			 * @param d Number of dimensions
			 * @param n Level of refinement
			 * @return Index corresponding to coords
			 */
			static int gp2idx(float *coords, int d, int n);

			/**
			 * @param index The index from the sparse grid that is converted to coordinates
			 * @param coords The computed coordinates of size d
			 * @param d Number of dimensions
			 * @param n Level of refinement
			 * @return If successful, returns 0
			 */
			static int idx2gp(int index, float *coords, int d, int n);

			/**
			 * @param x A real number in the interval [a, b]
			 * @param a The left margin of the interval
			 * @param b The right margin of the interval
			 * @return Returns the refinement level corresponding to x
			 */
			static int computeLevel(float x, float a, float b);

			/**
			 * @param levels The l component (of size d)
			 * @param indices The i component (of size d)
			 * @param d Number of dimensions
			 * @return Index corresponding to the (levels, indices) pair in a zero boundary sparse grid
			 */
			static int zb_gp2idx(int *levels, int *indices, int d);

			/**
			 * @param coords The coordinates of the sparse grid point
			 * @param d Number of dimensions
			 * @return Index corresponding to coords in a zero boundary sparse grid
			 */
			static int zb_gp2idx(float *coords, int d);

			/**
			 * @param index The index from the zero boundary sparse grid that is converted to (l, i)
			 * @param levels The computed l component (of size d)
			 * @param indices The computed i component (of size d)
			 * @param d Number of dimensions
			 * @return If successful, returns 0
			 */
			static int zb_idx2gp(int index, int *levels, int *indices, int d);

			/**
			 * @param index The index from the zero boundary sparse grid that is converted to coordinates
			 * @param coords The computed coordinates of size d
			 * @param d Number of dimensions
			 * @return If successful, returns 0
			 */
			static int zb_idx2gp(int index, float *coords, int d);
		};
}

#endif /* COORDINATES_H_ */
