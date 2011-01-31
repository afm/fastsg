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
	private:
		int d;

	public:
		SampleFct(int d)
		{
			this->d = d;
		}
	
		int getD()
		{
			return d;
		}
	
		float getValue(float *coords)
		{
			int i;
			float prod = 1;

			for (i = 0; i < d; i++)
				prod *= coords[i] * (2 - coords[i]);

			return prod;
		}
};

std::vector<int> visited;
int generate_points(int const_d, int const_l, float* gp, int crt_d,  int n, int numGridPoints)
{
	int i, j, count = 0;
	int levels[const_d], indices[const_d];
	int val;
	if (crt_d == -1) {
		Converter::coord2li(gp, levels, indices, const_d);

		val = Converter::gp2idx(levels, indices, const_d, const_l);

		if (val < numGridPoints) {
			visited[val] = 1;
		}
		else
			std::cout<<"Index out of range [" << val << "]" << std::endl;

		return 1;
	} else {
		gp[crt_d] = 0.0f;
		count += generate_points(const_d, const_l, gp, crt_d - 1, n, numGridPoints);
		gp[crt_d] = 1.0f;
		count += generate_points(const_d, const_l, gp, crt_d - 1, n, numGridPoints);

		/* for each level */
		for (i = 0; i < n; i++) {
			/* for each point on the level in dimension crt_d */
			for (j = 0; j < (1 << i); j++) {
				gp[crt_d] = (1.0f / (1 << i)) * (j + 0.5f);
				count += generate_points(const_d, const_l, gp, crt_d - 1, n - i, numGridPoints);
			}
		}
	}

	return count;
}

/*
 * test if n0gp2idx generates all consecutive indices from 0 to nrGridPoints-1
 */
int testgp2idx (int d, int l) {
	int b = 0;
	int i, numGridPoints = 0;
	float gp[d];
	/* for each dimension; 0-dimensional sparse grids are valid! */
	for (i = 0; i <= d; i++) {
		numGridPoints += (1 << i) * Helper::combi(d, i) * Helper::zerob_size(d - i, l);
	}
	
	visited.resize(numGridPoints);
	generate_points(d, l, gp, d - 1, l, numGridPoints);
	
	for (i = 0; i < numGridPoints; i++) {
		if (!visited[i]) {
			b = 1;
			goto stop;
		}
	}
	
	stop:
	if (b) {
		cout << "gp2idx test .............................. [failed]" << endl;
		return 1;
	} else {
		cout << "gp2idx test .............................. [passed]" << endl;
		return 0;
	}
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
	  return 1;
  }
};

/*
 * use a set to test if n0idx2gp generates unique pairs (l,i)
 */
int testidx2gp(int d, int l) {
	int i;
	std::set<Pair, CompareVectors> mapli;
	int *lev, *ind;
	std::set<Pair, CompareVectors>::iterator it;
	// create an object which represents the function you want to use
	SampleFct fct(d);
	// create a SparseGrid object
	SparseGrid sgf = SparseGrid(l, &fct);

	int nrGridPoints = sgf.getNumOfGridPoints();
	dim = d;
	for (i = 0; i < nrGridPoints; i++) {
		lev = (int*)malloc(d*sizeof(int));
		ind = (int*)malloc(d*sizeof(int));
		if (!lev || !ind) {
			cout << "Allocation error" << endl;
			return 1;
		}
		Converter::idx2gp(i, lev, ind, sgf.getD(), sgf.getL());
		mapli.insert( std::make_pair(lev, ind));
	}

	stop:
	if (mapli.size() != nrGridPoints) {
		cout << "idx2gp test .............................. [failed]" << endl;
		cout << "Error: size is " << mapli.size() << " , expected size is " << nrGridPoints; 
		return 1;
	} else {
		cout << "idx2gp test .............................. [passed]" << endl;
		return 0;
	}
}

/*
 * test if n0idx2gp(n0gp2idx(point_on_grid)) = index
 */
int testBijection(int d, int l) {
	int b = 0, i;
	// create an object which represents the function you want to use
	SampleFct fct(d);
	// create a SparseGrid object
	SparseGrid sgf = SparseGrid(l, &fct);
	int nrGridPoints = sgf.getNumOfGridPoints();
	int lev[d], ind[d];

	for (i = 0; i < nrGridPoints; i++) {
		Converter::idx2gp(i, lev, ind, sgf.getD(), sgf.getL());
		if (i != Converter::gp2idx(lev, ind, sgf.getD(), sgf.getL())) {
			b = 1;
			break;
		}
	}

	stop:
	if (!b) {
		cout << "Bijection test ........................... [passed]" << endl;
		return 0;
	} else {
		cout << "Bijection test ........................... [failed]" << endl;
		return 1;
	}
}

/*
 * test hierarchization and evaluation return correct results
 */
int testSparseGridOps(int d, int l) {
	int b = 0, i;
	// create an object which represents the function you want to use
	SampleFct fct(d);

	// create a SparseGrid object
	SparseGrid sgf = SparseGrid(l, &fct);
	int nrGridPoints = sgf.getNumOfGridPoints();
	int lev[d], ind[d];
	float coords[d];
	
	sgf.hierarchize();

	for (i = 0; i < nrGridPoints; i++) {
		Converter::idx2gp(i, coords, sgf.getD(), sgf.getL());
		if (sgf.evaluate(coords) != fct.getValue(coords)) {
			b = 1;
			goto stop;
		}
	}

	stop:
	if (!b) {
		cout << "Operation test ........................... [passed]" << endl;
		return 0;
	} else {
		cout << "Operation test ........................... [failed]" << endl;
		return 1;
	}
}

int main()
{
	int maxDim = 5, maxL = 5;
	
	try {
		for (int d = 1; d <= maxDim; d++)
			for (int l = 1; l <= maxL; l++) {
				cout << "Testing d = " << d << ", l = " << l << endl;
				cout << "---------------------------------------------------" << endl;
				
				if (testgp2idx(d, l)) throw 1;
				if (testidx2gp(d, l)) throw 2;
				if (testBijection(d, l)) throw 3;
				if (testSparseGridOps(d, l)) throw 4;
		
				cout << endl;
			}
		
		cout << "Success";
	}
	catch (int e) {
		std::cout << "Test number "<< e <<" failed" << endl;
	}
    
    return 0;
}

