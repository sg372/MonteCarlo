#ifndef HAMILTONIAN_HPP_
#define HAMILTONIAN_HPP_

#include "defs.hpp"
#include "FermiBasis.hpp"

class Hamiltonian : public FermiBasis
{

public:

Hamiltonian(unsigned, unsigned, double, double, double);

double hoppingIntegral;
double interactionStrength;
double disorderStrength;

double getOnSiteEnergy(unsigned);

void printHMatrix();

void printJMatrix();

DoubleMatrix HMatrix;

private:

DoubleMatrix JMatrix;
DoubleVector onSiteEnergies;


void makeJMatrix();

void makeOnSiteEnergies();

void makeHamiltonian();

protected:



};

#endif  /* HAMILTONIAN_HPP_ */
