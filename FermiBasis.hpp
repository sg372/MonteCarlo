#ifndef FERMIBASIS_HPP_
#define FERMIBASIS_HPP_

#include "defs.hpp"

class FermiBasis
{


public:

//parameters defining the physical system
unsigned particles;
unsigned sites;
unsigned basisSize;

//Matrix of basis states
UnsignedMatrix basis;


//Constructor
FermiBasis(unsigned, unsigned);


//Print out a basis state
void printBasisState(unsigned);

private:

//Function to get basis size for system
unsigned getBasisSize(unsigned, unsigned);

//Create basis states in basis matrix
void populateBasisStates(void);

};

#endif
