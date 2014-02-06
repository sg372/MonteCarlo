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
