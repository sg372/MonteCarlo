#include <iostream>
#include "FermiBasis.hpp"

#define cerr std::cerr
#define cout std::cout
#define endl std::endl 

FermiBasis::FermiBasis(unsigned ss, unsigned ps) {

	particles = ps;
	sites = ss;

	if (particles > sites) {
		cout<< "Error: FermiBasis: more fermions than sites" << endl << endl;
		throw "Error: more sites than particles";
	}

	basisSize = getBasisSize(sites, particles);

	basis = UnsignedMatrix::Zero(basisSize, sites);
	populateBasisStates();
}

unsigned FermiBasis::getBasisSize(unsigned ss, unsigned ps) {
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

	unsigned p[particles];

	for (unsigned i = 0; i < particles; ++i) {
		p[i] = i;
	}

	for (int j = 0; j < basisSize; ++j) {
		for (int i = 0; i < particles; ++i) {
			basis(j, p[i]) = 1;
		}

		int breaker = 0;
		int shift = 0;
		while (breaker == 0 && shift < particles) {
			if (p[particles - shift - 1] < sites - shift - 1) {
				p[particles - shift - 1] = p[particles - shift - 1] + 1;
				for (int i = 0; i < shift; i++) {
					p[particles - shift + i] = p[particles - shift + i - 1] + 1;
				}
				breaker = 1;
			}
			shift += 1;

		}

	}

}

void FermiBasis::printBasisState(unsigned stateIndex) {

	cout<< basis.row(stateIndex) << endl;

}

