/* Class Liouvillian
 * 
 * Implementation of a Liouvillian master operator for a Pauli master equation
 * for the transitions between eigenstates of a Hamiltonian object.
 * 
 * Couplings to the bath are random in the basis of the Hamiltonian and the 
 * bath is in an infinite-temperature state (i.e. P(i->j) = P(j->i), such
 * that the steady state occupation probability is uniform.
 * 
 * (s parameter is included as Markus may wish to diagonalise and find long 
 * lived steady states, but I haven't implemented this yet.)
 * 
 * Constructor input parameters are:
 * double s_val :  to be implemented later for "s-biased master equations.
 * double cs : coupling strength between Hamiltonian system and bath.
 * int ss : number of sites in quantum sys (passed to Hamiltonian constructor).
 * int ps : number of fermions < ss on the lattice (passed to Hamiltonian).
 * double J : hopping integral (to Hamiltonian).
 * double dis : disorder strength (to Hamiltonian).
 * 
 * DoubleMatrix W is the Liouvillian matrix operator, currently stored as a
 * an eigen v3 matrix (see defs.hpp).
 * 
 * 
 * 
 */

#ifndef LIOUVILLIAN_HPP_
#define LIOUVILLIAN_HPP_

#include "defs.hpp"
#include "Hamiltonian.hpp"



class Liouvillian : public Hamiltonian
{

public:

Liouvillian(double, double, int, int, double, double);

DoubleMatrix W;

/* Methods for finding and changing "s" in an
 * s-biased master equation.  Currently for
 * future use only when exact diagonalisation of
 * an s-generalised Liouvillian is implemented 
 * and a full redefinition of W would be costly */
void resetSValue(double);

double getSValue();

private:

double s;

SelfAdjointEigenSolver<DoubleMatrix> eigenSystem;

void makeLiouvillian();

void makeCouplingParameters();

double couplingStrength;

DoubleVector couplingParameters;

};



#endif /* LIOUVILLIAN_HPP_ */
