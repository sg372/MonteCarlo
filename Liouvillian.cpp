#include <iostream>
#include <random>
#include <cmath>
#include "Liouvillian.hpp"

#define cerr std::cerr
#define cout std::cout
#define endl std::endl


Liouvillian::Liouvillian(double s_val, double cs, int ss, int ps,
		double J, double dis) :
		Hamiltonian(ss, ps, J, dis) {

	s = s_val;
	couplingStrength = cs;

	couplingParameters = DoubleVector::Zero(basisSize);

	W = DoubleMatrix::Zero(basisSize, basisSize);

	makeCouplingParameters();

    /* Find eigenvalues and eigenvectors by creating a
     * SelfAdjointEigenSolver object from the eigen 
     * package */
	eigenSystem = SelfAdjointEigenSolver<DoubleMatrix>(HMatrix);

	makeLiouvillian();

}

/* Not of use currently.  Will be used with a Lanczos
 * diagonaliser for s-ensemble calculations
 */
void Liouvillian::resetSValue(double new_s) {

	this->s = new_s;
	makeLiouvillian();
}

/* Creates a Liouvillian matrix for a given set of random 
 * couplings and the eigenstates of the Hamiltonian
 * Note that the eigenvalue dependence arises due to the 
 * Ohmic nature of the coupling.
 */
void Liouvillian::makeLiouvillian() {

	DoubleVector diag(basisSize);

	for (int eig_i = 0; eig_i < basisSize; ++eig_i) {
		for (int eig_j = 0; eig_j <= eig_i; ++eig_j) {
			for (int basis_k = 0; basis_k < basisSize; ++basis_k) {

				W(eig_i, eig_j) += pow(
						eigenSystem.eigenvectors()(basis_k, eig_i)
								* eigenSystem.eigenvectors()(basis_k, eig_j)
								* couplingParameters(basis_k), 2)
						* couplingStrength
						* (eigenSystem.eigenvalues()(eig_i)
								- eigenSystem.eigenvalues()(eig_j));

				W(eig_j, eig_i) = W(eig_i, eig_j);
			}
		}
	}

	//Infer diagonal elements from probability conservation
	for (int eig_i = 0; eig_i < basisSize; ++eig_i) {
		diag(eig_i) = 0.0 - W.row(eig_i).sum();
	}

	//Apply s bias
	W = W * exp(-s);

	//Insert diagonals back
	for (int eig_i = 0; eig_i < basisSize; ++eig_i) {
		W(eig_i, eig_i) = diag(eig_i);
	}

}


/* Create a vector of random couplings to be used to construct the
 * Liouvillian matrix.
 */
void Liouvillian::makeCouplingParameters() {

	std::default_random_engine generator;
	std::normal_distribution<double> distribution(0.0, 1.0);

	//Assign coupling parameters
	for (int eig_i = 0; eig_i < basisSize; ++eig_i) {
		couplingParameters(eig_i) = distribution(generator);
	}
}

