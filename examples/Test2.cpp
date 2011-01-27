/**********************************************************************************
 *
 * Copyright (c) 2010, 2011 Alin Murarasu, Aurora Mirea
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

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <iostream>
#include <vector>
#include <set>

#include "SparseGrid.h"
#include "Converter.h"
#include "Helper.h"

using namespace std;

using namespace fsg;

class SampleFct : public Function
{
public:
	float getValue (float *coords, int d)
	{
		int i;
		float prod = 1;

		for (i = 0; i < d; i++)
			prod *= coords[i] * (2 - coords[i]);

		return prod;
	}
};

/*
 * test if n0gp2idx generates all consecutive indices from 0 to nrGridPoints-1
 */
int testgp2idx () {
	int b = 0;

	SampleFct fct;
	// specify d = number of dimensions and l = refinement level
	int d = 5, l = 4;

	// create a SparseGrid object
	SparseGrid sgf = SparseGrid(d, l, &fct);

	for (int i = 0; i < sgf.getNumOfGridPoints(); i++) {
		if (!Helper::visited[i]) {
			printf("Visited empty on position %i\n", i);
			b = 1;
		}
	}
	if (b) {
		printf("******************** Test gp2idx failed! ************************\n");
		return 1;
	}
	else
	{
		printf("************** Test gp2idx passed! ****************************\n");
		return 0;
	}
	return 0;
}


int dim;
typedef std::pair<int*,int*> Pair;
/*
 * Comparator for (l,i) pairs
 */
struct CompareVectors {
  bool operator ()(const Pair& p1, const Pair& p2) const {
	  int* l1 = (int*) p1.first;
	  int* l2 = (int*) p2.first;
	  int* i1 = (int*) p1.second;
	  int* i2 = (int*) p2.second;
	  int i = 0;
	  for (i = 0; i < dim; i++) {
		if (l1[i] != l2[i] || i1[i] != i2[i]) {

			return (l1[i] != l2[i] || i1[i] != i2[i]);
		}
	  }
	  return 1==1;
  }
};

/*
 * use a set to test if n0idx2gp generates unique pairs (l,i)
 */
int testidx2gp () {
	int i;
	std::set<Pair, CompareVectors> mapli;
	int *lev, *ind;
	std::set<Pair, CompareVectors>::iterator it;
	// create an object which represents the function you want to use
	SampleFct fct;
	// specify d = number of dimensions and l = refinement level
	int d = 5, l = 4;
	// create a SparseGrid object
	SparseGrid sgf = SparseGrid(d, l, &fct);

	int nrGridPoints = sgf.getNumOfGridPoints();
	dim = d;
	for (i = 0; i < nrGridPoints; i++) {
		lev = (int*)malloc(d*sizeof(int));
		ind = (int*)malloc(d*sizeof(int));
		if (!lev || !ind) {
			printf("Allocation error!");
			return 1;
		}
		Converter::idx2gp(i, lev, ind, sgf.sg.d, sgf.sg.l);
		mapli.insert( std::make_pair(lev, ind));
	}

	if (mapli.size() != nrGridPoints) {
		printf("******** Test idx2gp failed! Size is %i, expected size is %i ************\n", mapli.size(), nrGridPoints);
		return 1;
	}
	else {
		printf("*********** Test idx2gp passed! Size of set is %i ***********************\n", mapli.size());
		return 0;
	}
	return 0;
}

/*
 * test if n0idx2gp(n0gp2idx(point_on_grid)) = index
 */
int testBijection() {
	int b = 0, i;
	// create an object which represents the function you want to use
	SampleFct fct;
	// specify d = number of dimensions and l = refinement level
	int d = 5, l = 4;
	// create a SparseGrid object
	SparseGrid sgf = SparseGrid(d, l, &fct);
	int nrGridPoints = sgf.getNumOfGridPoints();
	int lev[d], ind[d];

	for (i = 0; i < nrGridPoints; i++) {
		Converter::idx2gp(i, lev, ind, sgf.sg.d, sgf.sg.l);
		if (i != Converter::gp2idx(lev, ind, sgf.sg.d, sgf.sg.l)) {
			b = 1;
			break;
		}
	}
	if (!b) {
		printf ("************************** Bijection test passed! *****************************\n");
	}
	else {
		printf("************************** Bijection test failed! **************************\n);");
		return 1;
	}
	return 0;
}

/*
 * test hierarchization and evaluation return correct results
 */
int testSparseGridOps(int d, int l) {
	int b = 0, i;
	// create an object which represents the function you want to use
	SampleFct fct;
	int di, li;

	for (di = 1; di <= d; di++) {
		for (li = 1; li <= 4; li++) {
			// create a SparseGrid object
			SparseGrid sgf = SparseGrid(di, li, &fct);
			int nrGridPoints = sgf.getNumOfGridPoints();
			int lev[di], ind[di];
			float coords[di];
			
			sgf.hierarchize();

			for (i = 0; i < nrGridPoints; i++) {
				Converter::idx2gp(i, lev, ind, sgf.sg.d, sgf.sg.l);
				Converter::li2coord(lev, ind, coords, sgf.sg.d);
				if (sgf.evaluate(coords) != fct.getValue(coords, sgf.sg.d)) {
					b = 1;
					goto stop;
				}
			}
		}
	}

	stop:
	if (!b) {
		cout << "************************** Operation test passed! *****************************" << endl;
	}
	else {
		cout << "************************** Operation test failed! **************************" << endl;
		return 1;
	}

	return 0;
}

int main()
{
	try {
		if (testgp2idx()) throw 1;

		if (testidx2gp()) throw 2;

		if (testBijection()) throw 3;
		
		if (testSparseGridOps(5, 4)) throw 4;
	}
	catch (int e) {
		std::cout<<"Test number "<<e<<" failed"<<std::endl;
	}
    return 0;
}
