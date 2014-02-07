/* Constructs the basis for ps fermions on ss sites.  The int
 * matrix of states is basis(basisSize, sites), such that a
 * particular state is represented by a row of integers. 
 * 
 * Particles are represented by 1, with holes 0.  E.g. for a 
 * six-site lattice with 3 particles in the first three sites, 
 * the state is represented by (1,1,1,0,0,0).
 * The states are generated iteratively with the rightmost 
 * moveable particle moved to the right for consecutive states.
 * 
 * The constuctor is passed:
 * int ps : the number of particles 
 * int ss : the number of sites (>= ps)
 * and throws a char * if the user asks for an unphysical
 * system or an empty system.
 *
 * The basis size (exponentially large in system size) is 
 * calculated and held in basisSize.
 * 
 * There exists a method to print out a basis state for
 * inspection and members holding the number of sites and
 * particles.
 */


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
