#ifndef FERMIBASIS_HPP_
#define FERMIBASIS_HPP_

#include "defs.hpp"

class FermiBasis
{


public:

//parameters defining the physical system
int particles;
int sites;
int basisSize;

//Matrix of basis states
intMatrix basis;


//Constructor
FermiBasis(int, int);


//Print out a basis state
void printBasisState(int);

private:

//Function to get basis size for system
int getBasisSize(int, int);

//Create basis states in basis matrix
void populateBasisStates(void);

};

#endif
