#ifndef LIOUVILLIAN_HPP_
#define LIOUVILLIAN_HPP_


#include "defs.hpp"
#include "Hamiltonian.hpp"



class Liouvillian : public Hamiltonian
{

public:

Liouvillian(double, double, unsigned, unsigned, double, double, double);


double couplingStrength;

DoubleVector couplingParameters;
DoubleMatrix W;

SelfAdjointEigenSolver<DoubleMatrix> eigenSystem;

void resetSValue(double);

double getSValue();

private:

double s;

void makeLiouvillian();
void makeCouplingParameters();

};



#endif /* LIOUVILLIAN_HPP_ */
