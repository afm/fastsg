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

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include <iostream>

#include "SparseGrid.h"
#include "Converter.h"
#include "Helper.h"
#include "Function.h"


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
				prod *= coords[i] * (3 - coords[i]);

			return prod;
		}
};

int testDemoFunc()
{
	// specify d = number of dimensions and l = refinement level
	int d = 5, l = 4;
	// create an object which represents the function you want to use
	SampleFct fct(d);
	
	// create a SparseGrid object
	SparseGrid sgf = SparseGrid(l, &fct);

	/********** demo of methods contained in SparseGrid class *****************/
	printf ("\n######### Demo of SparseGrid functions #################\n");
	// compute the hierarchical coefficients
	assert(sgf.hierarchize() == 0);

	int lc[sgf.getD()], ic[sgf.getD()];
	float val_ev, coords[d];

	Converter::idx2gp(1, lc, ic, sgf.getD(), sgf.getL());
	Converter::li2coord(lc, ic, coords, sgf.getD());

	// evaluate a point specified by coords
	val_ev = sgf.evaluate(coords);
	printf("Value from evaluation is %f; expecting: %f\n", val_ev, fct.getValue(coords));

	int ln[sgf.getD()],in[sgf.getD()];
	Converter::idx2gp(2, lc, ic, sgf.getD(), sgf.getL());

	// go to the next sparse grid
	sgf.next(lc, ic, ln, in);

	// print result
	printf("Next level starts at:\n");
	printf("l = ");
	for (int i = 0; i < d; i++) {
		printf ("%i ", ln[i]);
	}
	printf("\ni = ");
	for (int i = 0; i < d; i++) {
		printf ("%i ", in[i]);
	}

	return 0;
}

int demoConverter()
{
	int index = 2;
	int d = 3 , l = 3;
	int lev[d], idx[d];
	float coords[d];

	printf ("\n######### Demo of conversion functions #################\n");
	// convert a index from the linearized structure to (l,i) pair

	Converter::idx2gp(index, lev, idx, d, l);
	printf("Index %i corresponds to ", index);
	for (int i = 0; i < d; i++) {
		printf ("(%i %i)", lev[i], idx[i]);
	}
	printf ("\n and in coordinates is ");
	Converter::li2coord(lev, idx, coords, d);
	for (int i = 0; i < d; i++)
		printf("%f ", coords[i]);

	index = Converter::gp2idx(lev, idx, d, l);
	printf("\nWe convert (l, i) back into index (should be initial index) ");
	printf("\n Result of conversion is %i\n", index);


	return  0;
}

int main()
{
	try {
		if (testDemoFunc()) throw 1;

		if (demoConverter()) throw 2;
	}
	catch (int e) {
		std::cout<<"Test number "<<e<<" failed"<<std::endl;
	}
    return 0;
}

