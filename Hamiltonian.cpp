#include <iostream>
#include <random>
#include <ctime>
#include "Hamiltonian.hpp"

#define cerr std::cerr
#define cout std::cout
#define endl std::endl 

Hamiltonian::Hamiltonian(int ss, int ps, double J, double dis) :
		FermiBasis(ss, ps) {

	hoppingIntegral = J;
	disorderStrength = dis;

	onSiteEnergies = DoubleVector::Zero(sites);
	JMatrix = DoubleMatrix::Zero(sites, sites);
	HMatrix = DoubleMatrix::Zero(basisSize, basisSize);

	makeOnSiteEnergies();
	makeJMatrix();

	makeHamiltonian();

}


/* Member function to print to stdout internal matrices for
inspection and return model parameters including the generated
disorder realisation across the lattice */
void Hamiltonian::printHMatrix() {
	//Check that matrix is small (and so readily inspected)
	if (basisSize < 40) {
		cout<< HMatrix << endl;
	}
	else {
		cout << "Matrix too large for inspection!" << endl;
	}
}

void Hamiltonian::printJMatrix() {
	cout<< JMatrix << endl;
}

double Hamiltonian::getOnSiteEnergy(int site) {
	return onSiteEnergies(site);
}

double Hamiltonian::getHoppingIntegral(){
    return hoppingIntegral;
}

double Hamiltonian::getDisorderStrength(){
    return disorderStrength;
}

DoubleVector Hamiltonian::getOnSiteEnergies() {
	return onSiteEnergies;
}

/* Make matrix of sites x sites which describes topology 
 * of the lattice with periodic boundary conditions and 
 * has magnitude = hopping integral */
void Hamiltonian::makeJMatrix() {
	for (int i = 0; i < sites - 1; ++i) {
		JMatrix(i, i + 1) = hoppingIntegral;
		JMatrix(i + 1, i) = hoppingIntegral;
	}
	JMatrix(0, sites - 1) = hoppingIntegral;
	JMatrix(sites - 1, 0) = hoppingIntegral;

}

/* Populate DoubleVector of random on-site energies from 
 * Gaussian with std dev = disorderStrength */
void Hamiltonian::makeOnSiteEnergies() {

	std::default_random_engine generator(std::time(0));
	std::normal_distribution<double> distribution(0.0, disorderStrength);

	for (int i = 0; i < sites; ++i) {
		onSiteEnergies(i) = distribution(generator);
	}
}

/* Populate elements of basisSize x basisSize Hamiltonian
 * matrix using the specified hopping matrix elements and
 * the generated on-site energies */
void Hamiltonian::makeHamiltonian() {

	for (int i = 0; i < basisSize; ++i) {
		for (int j = 0; j < basisSize; ++j) {
			int s = 0;
			int breaker = 0;
			int match[3] = { 0 };

            //Check for two differences between states so that a hop can occur
			while (breaker < 3 && s < sites) { 
				if (basis(i, s) != basis(j, s)) {
					match[breaker] = s;
					breaker += 1;
				}
				s += 1;
			}

            //Populate elements if a flip can occur (two differences in states)
			if (s == sites && breaker == 2) {

                //Check that we are moving a particle from a hole
				if (basis(i, match[0]) == 0 && basis(i, match[1]) == 1
						&& basis(j, match[0]) == 1 && basis(j, match[1]) == 0) {

                    //Count positions to preserve Fermi commutation relations
					int sum1 = 0;
					for (int m = 0; m < match[1]; ++m) {
						sum1 = sum1 + basis(i, m);
					}

					int sum2 = 0;
					for (int m = 0; m < match[0]; ++m) {
						sum2 = sum2 + basis(j, m);
					}

					int sgn = 2 * (sum1 + sum2) % 2 - 1;

					HMatrix(i, j) = HMatrix(i, j)
							+ JMatrix(match[0], match[1]) * sgn;

                //Repeat for moving a particle TO a hole
				} else if (basis(i, match[0]) == 1 && basis(i, match[1]) == 0
						&& basis(j, match[0]) == 0 && basis(j, match[1]) == 1) {

					int sum1 = 0;
					for (int m = 0; m < match[0]; ++m) {
						sum1 = sum1 + basis(i, m);
					}

					int sum2 = 0;
					for (int m = 0; m < match[1]; ++m) {
						sum2 = sum2 + basis(j, m);
					}

					int sgn = 2 * (sum1 + sum2) % 2 - 1;

					HMatrix(i, j) = HMatrix(i, j)
							+ JMatrix(match[0], match[1]) * sgn;
				}
			}

		}
	}

    // Add diagonal disorder energies 
	for (int i = 0; i < basisSize; ++i) {
		for (int k = 0; k < sites; ++k) {
			if (basis(i, k) == 1) {
				HMatrix(i, i) += onSiteEnergies(k);
			}
		}
	}


}

