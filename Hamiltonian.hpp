/* Hamiltonian class contructs a Hamiltonian matrix for a
 * quantum system of fermionic particles on a one dimensional
 * lattice with PBC for a given tunnelling strength and a given 
 * Gaussian disorder with standard deviation to be set. I.e.:
 *
 * $H = \sum_i e_i c^\dag_i c_i + J c^\dag_i c_{i+1} + h.c.$
 *
 * Inherits from FermiBasis, which is used to construct the matrix
 * elements of the Hamiltonian and may later be used to refer to 
 * particular Fock basis states if required via public methods in
 * FermiBasis.
 * 
 * The constructor takes:
 * int ss : number of sites
 * int ps : number of particles
 * double J : hopping integral
 * double dis : disorder strength (std dev of a Gaussian).
 *
 * The matrix data member
 * DoubleMatrix HMatrix 
 * is available to other classes, for example if studying unitary
 * or open system dynamics of the model.
 * 
 * Public members allow inspection of the site hopping matrix,
 * small Hamiltonian matrices and specified parameters.
 */ 

#ifndef HAMILTONIAN_HPP_
#define HAMILTONIAN_HPP_

#include "defs.hpp"
#include "FermiBasis.hpp"

class Hamiltonian : public FermiBasis
{

public:

Hamiltonian(int, int, double, double);



double getOnSiteEnergy(int);

DoubleVector getOnSiteEnergies();

double getHoppingIntegral();

double getDisorderStrength();

void printHMatrix();

void printJMatrix();


private:

DoubleMatrix JMatrix;

DoubleVector onSiteEnergies;

double hoppingIntegral;

double disorderStrength;

void makeJMatrix();

void makeOnSiteEnergies();

void makeHamiltonian();

protected:

DoubleMatrix HMatrix;


};

#endif  /* HAMILTONIAN_HPP_ */
