#include <iostream>
#include "FermiBasis.hpp"

#define cerr std::cerr
#define cout std::cout
#define endl std::endl 

FermiBasis::FermiBasis(int ss, int ps) {

	particles = ps;
	sites = ss;

	//Check the system is physical
	if (particles > sites || particles < 1 || sites <1) {
		cout<< "Error: FermiBasis: more fermions than sites" << endl << endl;
		throw "Error: more sites than particles";
	}

	basisSize = getBasisSize(sites, particles);

	basis = intMatrix::Zero(basisSize, sites);
	populateBasisStates();
}

int FermiBasis::getBasisSize(int ss, int ps) {

	//Calculate [ss Choose ps] in a stable way:
	if (ps > ss)
		return 0;
	if (ps * 2 > ss)
		ps = ss - ps;
	if (ps == 0)
		return 1;

	int result = ss;
	for (int i = 2; i <= ps; ++i) {
		result *= (ss - i + 1);
		result /= i;
	}
	return result;
}

void FermiBasis::populateBasisStates(void) {

	int p[particles]; //Contains the positions of the particles

	//Initialise particles on the left of the lattice
	for (int i = 0; i < particles; ++i) {
		p[i] = i;
	}

	//Build all basis states, of which there are basisSize
	for (int j = 0; j < basisSize; ++j) {

		for (int i = 0; i < particles; ++i) {
			//Put particles in basis state j
			basis(j, p[i]) = 1; 
		}

		//To determine if the particle can be iterated no further
		int breaker = 0; 

		/* To identify the next particle to moved, up to the leftmost
		 * particle, with shift=0 being the rightmost particle.*/
		int shift = 0;   

		while (breaker == 0 && shift < particles) {
			if (p[particles - shift - 1] < sites - shift - 1) {
				p[particles - shift - 1] = p[particles - shift - 1] + 1;
				for (int i = 0; i < shift; i++) {
					p[particles - shift + i] = p[particles - shift + i - 1] + 1;
				}
				breaker = 1; //particle can move no further
			}
			/* take the particle to the left and repeat until it is next to the
			 * last particle */
			shift += 1; 

		}

	}

}

//Print the integers corresponding to the state stateIndex
void FermiBasis::printBasisState(int stateIndex) {

	cout<< basis.row(stateIndex) << endl;

}

