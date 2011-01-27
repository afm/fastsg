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
 
/**
 * Description: 
 *	Dynamic programming for computing the number of vectors for which 
 *  their component sum is l and a restriction (maximum) for each individual 
 *  component of the vector is given.
 *	This example contains three solutions for computing the number of
 *  solutions. The count_valid_levels is based on brute force. The 
 *  dyn_count_valid_levels function exploits the recursive property:
 *  num_solutions(d, l) = sum(num_solutions(d - 1, i)), i = 0..l. The last
 *  is an improvement reducing the number of for loops.
 *
 * The routines can be used to build a bijection for dimensionally adaptive 
 * sparse grids.
 */
 
#include <stdio.h>

int limits[10] = {1, 2, 1, 2, 3, 1, 4, 5, 6, 1};
int levels[10];

#define MIN(x, y)	((x < y)? x: y)

/****************************************************************************
 *  brute force function
 ****************************************************************************/
int count_valid_levels(int *limits, int *levels, int cd, int cl, int d, int l)
{
	int i, count = 0;
	
	if (cd == 0) {
		if (cl > limits[0])
			return 0;
		levels[0] = cl;

		return 1;
	}

	for (i = 0; i <= MIN(cl, limits[cd]); i++) {
		levels[cd] = i;
		count += count_valid_levels(limits, levels, cd - 1, cl - i, d, l);
	}
		
	return count;
}

/****************************************************************************
 *  dynamic programming solution
 ****************************************************************************/
int dyn_count_valid_levels(int *limits, int *levels, int d, int l)
{
	int count, i, j, k;
	int a[d][l + 1];
	
	// compute a matrix a[][] that has the meaning:
	// a[i][j] represents the number of possibilities 
	// of obtaining sum j using the first i elements 
	// of levels given the restrictions in limits

	for (i = 0; i < d; i++) {
		for (j = 0; j <= l; j++) {
			a[i][j] = 0;
		}
	}

	// base case
	for (j = 0; j <= MIN(l, limits[0]); j++)
		a[0][j] = 1;
	for (i = 1; i < d; i++)
		a[i][0] = 1;
	
	// loop over dimensions
	for (i = 1; i < d; i++) {
		// loop over sums
		for (j = 1; j <= l; j++) {
			for (k = 0; k <= MIN(j, limits[i]); k++)
				a[i][j] += a[i - 1][j - k];
		}
	}

	count = a[d - 1][l];

	return count;
}

/****************************************************************************
 *  optimized dynamic programming solution
 ****************************************************************************/
int opt_count_valid_levels(int *limits, int *levels, int d, int l)
{
	int count, i, j, k;
	int a[d][l + 1];
	
	// compute a matrix a[][] that has the meaning:
	// a[i][j] represents the number of possibilities 
	// of obtaining sum j using the first i elements 
	// of levels given the restrictions in limits

	for (i = 0; i < d; i++) {
		for (j = 0; j <= l; j++) {
			a[i][j] = 0;
		}
	}

	// base case
	for (j = 0; j <= MIN(l, limits[0]); j++)
		a[0][j] = 1;
	for (i = 1; i < d; i++)
		a[i][0] = 1;

	// loop over dimensions
	for (i = 1; i < d; i++) {
		// loop over sums
		for (j = 1; j <= MIN(l, limits[i]); j++) {
			a[i][j] = a[i][j - 1] + a[i - 1][j];	
		}
		
		for (j = MIN(l, limits[i]) + 1; j <= l; j++) {
			a[i][j] = a[i][j - 1] - a[i - 1][j - 1 - limits[i]] + a[i - 1][j];
		}
	}

	count = a[d - 1][l];

	return count;
}

int main()
{
	printf("Recursive: num. of valid levels ................ %d\n", count_valid_levels(limits, levels, 9, 5, 10, 5));
	
	printf("Dynamic programming: num. of valid levels ...... %d\n", dyn_count_valid_levels(limits, levels, 10, 5));

	printf("Optimized: num. of valid levels ................ %d\n", opt_count_valid_levels(limits, levels, 10, 5));

	return 0;
}

