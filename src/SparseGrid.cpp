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

#include "SparseGrid.h"
#include "DataStructure.h"
#include "Converter.h"
#include "Helper.h"

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

using namespace fsg;

SparseGrid::SparseGrid(int l, Function* f)
{
	float gp[f->getD()];
	int count;
	int i;

	d = f->getD();

	try {
		if (d < 0 || l < 0)
			throw 1;
		this->d = d;
		this->l = l;
		numOfGridPoints = size(d, l);
		
		sg1d = (float*) malloc(numOfGridPoints * sizeof(float));

		for (i = 0; i < numOfGridPoints; i++) {
			Converter::idx2gp(i, gp, d, l);
			sg1d[i] = f->getValue(gp);
		}
	} catch (int e) {
		std::cout
				<< "Exception: number of dimensions and refinement level must be positive!"
				<< std::endl;
	}
}

SparseGrid::~SparseGrid()
{
	free(sg1d);
}

/* evaluates (or interpolates) the sparse grid at point coords inside the [0, 1]^d domain */
float SparseGrid::evaluate(float *coords)
{
	int k, i, index1, index2, t0, pd, kk;
	float left, prod, val = 0, div, m, prod0;
	int indices[d], plevels[d], levels[d];
	float pcoords[d];
	float *sg1d = this->sg1d;

	try {
		for (i = 0; i < d; i++)
			if (coords[i] > 1 || coords[i] < 0)
				throw 1;
		val = 0.0f;
		index1 = 0;

		/* loop over groups of sparse grids of the same dimensionality
		 pd = projection dimensionality */
		for (pd = d; pd >= 0; pd--) {
			/* loop over sparse grids of the same dimensionality */
			for (kk = 0; kk < (1 << (d - pd)) * Helper::combi(d, d - pd); kk++) {
				/* convert index pointing to the current sparse grid to (l, i) */
				Converter::idx2gp(index1, levels, indices, d, l);

				/* move index to next sparse grid in the group */
				index1 += Helper::zerob_size(pd, l);

				/* prod0 is the same for all the regular grids composing the current sparse grid */
				prod0 = 1.0f;
				i = 0;
				for (k = 0; k < d; k++) {
					if (levels[k] == -1) {
						if (indices[k] == 0)
							prod0 *= (1 - coords[k]);
						else
							prod0 *= coords[k];
					} else {
						pcoords[i++] = coords[k];
					}
				}

				/* no need to proceed if the sparse grids are 0-dimensional */
				if (pd == 0) {
					val += prod0 * sg1d[0];
					sg1d++;
					continue;
				}

				/* initialize plevels with 0
				 plevels = projection levels */
				memset(plevels, 0, pd * sizeof(int));

				/* start evaluation of 0-boundary sparse grids */
				for (i = 0; i < l; i++) {
					plevels[0] = 0;
					plevels[pd - 1] = i;
					do {
						/* initilize production with initial product! */
						prod = prod0;
						index2 = 0;
						/* multiply pd 1-dimensional hat functions */
						for (k = 0; k < pd; k++) {
							div = (1.0f - 0.0f) / (1 << plevels[k]);
							index2 = index2 * (1 << plevels[k])
									+ (int) ((pcoords[k] - 0.0f) / div);
							left = (int) ((pcoords[k] - 0.0f) / div) * div;
							m = (2.0f * (pcoords[k] - left) - div) / div;
							prod *= 1.0f + m * ((m < 0.0f) - !(m < 0.0f));
						}

						/* multiply with corresponding hierarchical coefficient */
						prod *= sg1d[index2];
						/* add contribution to the interpolation result */
						val += prod;

						/* move to the next regular (full) grid of the current sparse grid of dimensionality pd */
						sg1d += 1 << i;

						/* if the end of the group of regular grids is reached, stop */
						if (plevels[0] == i)
							break;

						/* otherwise, use iterator to generate the next valid levels */
						k = 1;
						while (plevels[k] == 0)
							k++;
						plevels[k]--;
						t0 = plevels[0];
						plevels[0] = 0;
						plevels[k - 1] = t0 + 1;
					} while (1);
				}
				/* end evaluation of 0-boundary sparse grids */
			}
		}
		
	} catch (int i) {
		std::cout << "The coordinates are not in [0,1]^d domain" << std::endl;
	}
	return val;
}

/* evaluates (or interpolates) the sparse grid at points stored in coords inside the [0, 1]^d domain */
int SparseGrid::evaluate(float *coords, int n, float *vals)
{
	int k, i, j, index1, index2, t0, pd, kk;
	float left, prod, div, m, prod0s[n];
	int indices[d], plevels[d], levels[d];
	float pcoords[n][d];
	float *sg1d = this->sg1d;
	float (*nxcoords)[d] = (float (*)[d]) coords;

	for (j = 0; j < n; j++)
		vals[j] = 0;

	try {
		for (j = 0; j < n; j++)
			for (i = 0; i < d; i++)
				if (nxcoords[j][i] > 1 || nxcoords[j][i] < 0)
					throw 1;
		index1 = 0;

		/* loop over groups of sparse grids of the same dimensionality
		 pd = projection dimensionality */
		for (pd = d; pd >= 0; pd--) {
			/* loop over sparse grids of the same dimensionality */
			for (kk = 0; kk < (1 << (d - pd)) * Helper::combi(d, d - pd); kk++) {
				/* convert index pointing to the current sparse grid to (l, i) */
				Converter::idx2gp(index1, levels, indices, d, l);

				/* move index to next sparse grid in the group */
				index1 += Helper::zerob_size(pd, l);

				for (j = 0; j < n; j++) {
					/* for a given point, prod0 is the same for all the regular grids composing the current sparse grid */
					prod0s[j] = 1.0f;
					i = 0;
					for (k = 0; k < d; k++) {
						if (levels[k] == -1) {
							if (indices[k] == 0)
								prod0s[j] *= (1 - nxcoords[j][k]);
							else
								prod0s[j] *= nxcoords[j][k];
						} else {
							pcoords[j][i++] = nxcoords[j][k];
						}
					}
				}

				/* no need to proceed if the sparse grids are 0-dimensional */
				if (pd == 0) {
					for (j = 0; j < n; j++)
						vals[j] += prod0s[j] * sg1d[0];
					sg1d++;
					continue;
				}

				/* initialize plevels with 0
				 plevels = projection levels */
				memset(plevels, 0, pd * sizeof(int));

				/* start evaluation of 0-boundary sparse grids */
				for (i = 0; i < l; i++) {
					plevels[0] = 0;
					plevels[pd - 1] = i;
					do {
						/* for each evaluation point */
						for (j = 0; j < n; j++) {
							/* initilize production with initial product! */
							prod = prod0s[j];
							index2 = 0;
							/* multiply pd 1-dimensional hat functions */
							for (k = 0; k < pd; k++) {
								div = (1.0f - 0.0f) / (1 << plevels[k]);
								index2 = index2 * (1 << plevels[k])
										+ (int) ((pcoords[j][k] - 0.0f) / div);
								left = (int) ((pcoords[j][k] - 0.0f) / div) * div;
								m = (2.0f * (pcoords[j][k] - left) - div) / div;
								prod *= 1.0f + m * ((m < 0.0f) - !(m < 0.0f));
							}

							/* multiply with corresponding hierarchical coefficient */
							prod *= sg1d[index2];
							/* add contribution to the interpolation result */
							vals[j] += prod;
						}

						/* move to the next regular (full) grid of the current sparse grid of dimensionality pd */
						sg1d += 1 << i;

						/* if the end of the group of regular grids is reached, stop */
						if (plevels[0] == i)
							break;

						/* otherwise, use iterator to generate the next valid levels */
						k = 1;
						while (plevels[k] == 0)
							k++;
						plevels[k]--;
						t0 = plevels[0];
						plevels[0] = 0;
						plevels[k - 1] = t0 + 1;
					} while (1);
				}
				/* end evaluation of 0-boundary sparse grids */
			}
		}
		
	} catch (int i) {
		std::cout << "The coordinates are not in [0,1]^d domain" << std::endl;
		
		return -1;
	}
	
	return 0;
}

float SparseGrid::evaluate(int *levels, int *indices)
{
	float coords[d];

	Converter::li2coord(levels, indices, coords, d);

	return evaluate(coords);
}

/* 
 * computes the hierarchical coefficients for a d-dimesional, level n, non-0 boundary sparse grid
 * initially, sg1d contains function values 
 */
int SparseGrid::hierarchize()
{
	int i, j;
	float val1, val2;
	int levels[d], indices[d];
	int plevels[d], pindices[d];
	float *sg1d = this->sg1d;

	/* loop over dimensions */
	for (i = 0; i < d; i++)
		/* loop over grid points */
		for (j = size(d, l) - 1; j >= 0; j--) {
			/* convert index to (l, i) */
			Converter::idx2gp(j, levels, indices, d, l);

			/* retrieve left parent's value from sparse grid */
			if (get_lparent(levels, indices, plevels, pindices, i) != -1)
				val1 = sg1d[Converter::gp2idx(plevels, pindices, d, l)];
			else
				val1 = 0;

			/* retrieve right parent's value from sparse grid */
			if (get_rparent(levels, indices, plevels, pindices, i) != -1)
				val2 = sg1d[Converter::gp2idx(plevels, pindices, d, l)];
			else
				val2 = 0;
			/* update current hierarchical coefficient (at position j) */
			sg1d[j] = sg1d[j] - (val1 + val2) / 2.0f;
		}

	return 0;
}

/* returns the (l, i) of the left parent in dimension cd */
int SparseGrid::get_lparent(int *levels, int *indices, int *plevels, int *pindices, int cd)
{
	int i;
	float pc;

	/* if the point is located on the border, it has no parent */
	if (levels[cd] == -1)
		return -1;

	/* the (l, i) of the parent are the same as the point's (l, i) excepting the cd-th component */
	for (i = 0; i < d; i++) {
		plevels[i] = levels[i];
		pindices[i] = indices[i];
	}

	/* if point's index is 0, the left parent is  (l = -1, i = 0) */
	if (indices[cd] == 0) {
		plevels[cd] = -1;
		pindices[cd] = 0;

		return 0;
	}

	/* compute the coordinate of the left parent in dimension cd */
	pc = (1.0f / (1 << levels[cd])) * indices[cd];

	/* convert coordinate to (l, i) only for dimension cd */
	Converter::coord2li(&pc, plevels + cd, pindices + cd, 1);

	return 0;
}

/* returns the (l, i) of the right parent in dimension cd */
int SparseGrid::get_rparent(int *levels, int *indices, int *plevels, int *pindices, int cd)
{
	int i;
	float pc;

	/* if the point is located on the border then it has no parent */
	if (levels[cd] == -1)
		return -1;

	/* the (l, i) of the parent are the same as the point's (l, i) excepting the cd-th component */
	for (i = 0; i < d; i++) {
		plevels[i] = levels[i];
		pindices[i] = indices[i];
	}

	/* if the point preceds the right boder in dimension cd, its right parent is on the border */
	if (indices[cd] == (1 << levels[cd]) - 1) {
		plevels[cd] = -1;
		pindices[cd] = 1;

		return 0;
	}

	/* compute the coordinate of the right parent in dimension cd */
	pc = (1.0f / (1 << levels[cd])) * (indices[cd] + 1);

	/* convert coordinate to (l, i) only for dimension cd */
	Converter::coord2li(&pc, plevels + cd, pindices + cd, 1);

	return 0;
}

int SparseGrid::get_lparent(float *coords, float *pcoords, int cd)
{
	int levels[d], indices[d], plevels[d], pindices[d];
	Converter::coord2li(coords, levels, indices, d);
	get_lparent(levels, indices, plevels, pindices, cd);
	Converter::li2coord(plevels, pindices, pcoords, d);

	return 0;
}

int SparseGrid::get_rparent(float *coords, float *pcoords, int cd)
{
	int levels[d], indices[d], plevels[d], pindices[d];
	Converter::coord2li(coords, levels, indices, d);
	get_rparent(levels, indices, plevels, pindices, cd);
	Converter::li2coord(plevels, pindices, pcoords, d);

	return 0;
}

int SparseGrid::next(int *crt_levels, int *crt_indices, int *next_levels, int *next_indices)
{
	int index = Converter::gp2idx(crt_levels, crt_indices, d, l);
	int i, pd = 0;

	for (i = 0; i < d; i++)
		if (crt_levels[i] != -1)
			pd++;

	index += Helper::zerob_size(pd, l);

	Converter::idx2gp(index, next_levels, next_indices, d, l);

	return 0;
}

int SparseGrid::getNumOfGridPoints() const
{
	return numOfGridPoints;
}

/* the size of a non-zero boundary, d-dimensional, n-refined sparse grid */
int SparseGrid::size(int d, int n)
{
	int i, s = 0;

	/* for each dimension; 0-dimensional sparse grids are valid! */
	for (i = 0; i <= d; i++) {
		s += (1 << i) * Helper::combi(d, i) * Helper::zerob_size(d - i, n);
	}

	return s;
}

/* the size of the sparse grid */
int SparseGrid::size()
{
	return numOfGridPoints;
}

int SparseGrid::getD()
{
	return d;
}

int SparseGrid::getL()
{
	return l;
}

