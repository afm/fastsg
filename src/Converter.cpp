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

#include "Converter.h"
#include "Helper.h"

using namespace fsg;

/* zero boundary gp2idx */
int Converter::zb_gp2idx(int *levels, int *indices, int d)
{
	int index1, index2, index3, i, sum;

	sum = 0;
	index2 = 0;

	index1 = indices[0];
		for (i = 1; i < d; i++)
			index1 = (index1 << levels[i]) + indices[i];

	for (i = 0; i < d - 1; i++) {
		sum += levels[i];
		if (sum > 0)
			index2 += Helper::combi(i + sum, sum - 1);
	}
	sum += levels[i];
	index2 <<= sum;
	index3 = 0;
	for (i = 0; i < sum; i++) {
		index3 += (1 << i) * Helper::combi(d - 1 + i, i);
	}

	return index1 + index2 + index3;
}

/* zero boundary idx2gp */
int Converter::zb_idx2gp(int index, int *levels, int *indices, int d)
{
	int i, j, f, isum, sum, level, dindex, rest;

	f = 1;
	isum = 0;
	i = 0;
	while (index >= isum + Helper::combi(d - 1 + i, i) * f) {
		isum += Helper::combi(d - 1 + i, i) * f;
		f *= 2;
		i++;
	}

	sum = i;
	index -= isum;
	rest = index % (1 << i);
	index /= (1 << i);

	for (i = d - 2; i >= 0; i--) {
		isum = 0;
		j = 0;
		while (index >= isum + Helper::combi(i + j, j)) {
			isum += Helper::combi(i + j, j);
			j++;
		}
		level = sum - j;
		sum = j;
		f = (1 << level);
		dindex = rest % f;
		rest /= f;
		levels[i + 1] = level;
		indices[i + 1] = dindex;
		index -= isum;
	}

	level = sum;
	f = (1 << level);
	dindex = rest % f;
	rest /= f;
	levels[0] = level;
	indices[0] = dindex;

	return 0;
}

/* zero boundary gp2idx */
int Converter::zb_gp2idx(float *coords, int d) {
	float index1;
	int index2, index3, i, sum, level;

	sum = 0;
	index1 = index2 = 0;
	for (i = 0; i < d - 1; i++) {
		level = computeLevel(coords[i], 0, 1);
		index1 = (index1 + (coords[i] - 0.0) / (1.0 - 0.0)) * (1 << level) - 0.5 / (1.0 - 0.0);
		sum += level;
		index2 += Helper::combi(i + sum, sum - 1);
	}

	level = computeLevel(coords[i], 0, 1);
	index1 = (index1 + (coords[i] - 0.0) / (1.0 - 0.0)) * (1 << level) - 0.5 / (1.0 - 0.0);
	sum += level;
	index2 <<= sum;

	index3 = 0;
	for (i = 0; i < sum; i++) {
		index3 += (1 << i) * Helper::combi(d - 1 + i, i);
	}

	return (int) (index1 + 0.5) + index2 + index3;
}

/* zero boundary idx2gp */
int Converter::zb_idx2gp(int index, float *coords, int d) {
	int i, j, f, isum, sum, level, dindex, rest;

	f = 1;
	isum = 0;
	i = 0;
	while (index >= isum + Helper::combi(d - 1 + i, i) * f) {
		isum += Helper::combi(d - 1 + i, i) * f;
		f *= 2;
		i++;
	}

	sum = i;
	index -= isum;
	rest = index % (1 << i);
	index /= (1 << i);

	for (i = d - 2; i >= 0; i--) {
		isum = 0;
		j = 0;
		while (index >= isum + Helper::combi(i + j, j)) {
			isum += Helper::combi(i + j, j);
			j++;
		}
		level = sum - j;
		sum = j;
		f = (1 << level);
		dindex = rest % f;
		rest /= f;
		coords[i + 1] = 0 + (1 - 0) / (2.0 * f) + dindex * (1 - 0) / (1.0 * f);
		index -= isum;
	}

	level = sum;
	f = (1 << level);
	dindex = rest % f;
	rest /= f;
	coords[0] = 0 + (1 - 0) / (2.0 * f) + dindex * (1 - 0) / (1.0 * f);

	return 0;
}

/* non-zero gp2idx, wrapper around zb_gp2idx */
int Converter::gp2idx(int *levels, int *indices, int d, int n)
{
	int index1, index2, index3;
	int pd, n01;
	int plevels[d], pindices[d];
	int i;

	/* select points on the boundary */
	pd = 0;
	for (i = 0; i < d; i++)
		if (levels[i] != -1) {
			plevels[pd] = levels[i];
			pindices[pd++] = indices[i];
		}
	if (pd) {
		index1 = Converter::zb_gp2idx(plevels, pindices, pd);
	}
	else
		index1 = 0;

	/* select the right 0-boundary sparse grid (its beginning) */
	index2 = 0;
	n01 = d - pd;
	for (i = 0; i < d; i++) {
		if (levels[i] != -1) {
			index2 += (1 << n01) * Helper::combi(d - i - 1, n01 - 1);
		}
		else {
			n01--;

			if (indices[i] == 1)
				index2 += (1 << n01) * Helper::combi(d - i - 1, n01);
			}
	}
	index2 *= Helper::zerob_size(pd, n);

	/* count the number of grid points preceeding the group of 0-boundary sparse grids that contains the grid point */
	index3 = 0;

	for (i = 0; i < d - pd; i++) {
		index3 += (1 << i) * Helper::combi(d, i) * Helper::zerob_size(d-i, n);
	}
	return index1 + index2 + index3;
}

/* for a given index, returns equivalent (levels, indices) representation */
int Converter::idx2gp(int index, int *levels, int *indices, int d, int n)
{
	int i = 0, j;
	int n01, pd;
	int index1, index2;
	int plevels[d], pindices[d];

	/* after the while, i contains the number of -1 components in levels */
	while (index >= (1 << i) * Helper::combi(d, i) * Helper::zerob_size(d - i, n)) {

		index -= (1 << i) * Helper::combi(d, i) * Helper::zerob_size(d - i, n);

		i++;
	}

	n01 = i;
	for (i = 0; i < d; i++)
		levels[i] = 0;

	/* pd is the dimensionality of the projection that contains the grid point given by index */
	pd = d - n01;

	/* index1 is the index inside the sparse grid */
	index1 = index % Helper::zerob_size(pd, n);
	/* index2 is the index of the sparse grid from the beginning of its group */
	index2 = index / Helper::zerob_size(pd, n);

	/* convert index1 to (l, i) representation for the projection */

	Converter::zb_idx2gp(index1, plevels, pindices, pd);

	/* find the positions in levels of the -1 components and more... */
	j = 0;
	for (i = 0; i < d; i++) {
		if (index2 >= (1 << n01) * Helper::combi(d - i - 1, n01 - 1)) {
			levels[i] = plevels[j];
			indices[i] = pindices[j++];
			index2 -= (1 << n01) * Helper::combi(d - i - 1, n01 - 1);
		} else {
			levels[i] = -1;
			n01--;
			if (index2 >= (1 << n01) * Helper::combi(d - i - 1, n01)) {
				indices[i] = 1;
				index2 -= (1 << n01) * Helper::combi(d - i - 1, n01);
			} else {
				indices[i] = 0;
			}
		}
	}

	return 0;
}

/* converts coords (floats) to (l, i) representation (integers) */
int Converter::coord2li(float *coords, int *levels, int *indices, int d)
{
	int i;
	float cc;

	for (i = 0; i < d; i++) {
		if (coords[i] == 0.0f) {
			levels[i] = -1;
			indices[i] = 0;
		} else if (coords[i] == 1.0f) {
			levels[i] = -1;
			indices[i] = 1;
		} else {
			levels[i] = -1;
			cc = coords[i];
			while (cc != (int) cc) {
				cc *= 2;
				levels[i]++;
			}
			indices[i] = (int) ((cc - 1.0f) / 2.0f);
		}
	}

	return 0;
}

/* converts (l, i) to coords (floats) */
int Converter::li2coord(int *levels, int *indices, float *coords, int d)
{
	int i;

	for (i = 0; i < d; i++) {
		if (levels[i] == -1) {
			if (indices[i] == 0)
				coords[i] = 0.0f;
			else
				coords[i] = 1.0f;
		} else {
			coords[i] = (1.0f / (1 << levels[i])) * (indices[i] + 0.5f);
		}
	}

	return 0;
}

/* computes the refinement level of x located between in the interval [a, b] */
int Converter::computeLevel(float x, float a, float b)
{
	int i;
	float r = (x - a) / (b - a), d = 0.5f;

	i = 0;
	while (r != d) {
		if (r > d)
			r -= d;
		d /= 2.0f;
		i++;
	}

	return i;
}

/* returns the 1d index of the grid point coords */
int Converter::gp2idx(float *coords, int d, int n)
{
	int levels[d], indices[d];
	
	coord2li(coords, levels, indices, d);
	
	return gp2idx(levels, indices, d, n);
}

/* returns the coords of the grid point with index */
int Converter::idx2gp(int index, float *coords, int d, int n)
{
	int levels[d], indices[d];
	
	idx2gp(index, levels, indices, d, n);
	li2coord(levels, indices, coords, d);
}
