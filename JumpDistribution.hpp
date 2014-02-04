#ifndef JUMPDISTRIBUTION_HPP_

#define JUMPDISTRIBUTION_HPP_

#include <vector>
#include "defs.hpp"
#include "ContinuousTimeMonteCarlo.hpp"



class JumpDistribution{


public:

JumpDistribution( double, unsigned, double, double, ContinuousTimeMonteCarlo *);


double interval;
int bins;
double range[2];
double binWidth;

ContinuousTimeMonteCarlo * MonteCarloObjectPtr;

std::vector<double> midPointValues;
std::vector<unsigned> frequency;




private:

void makeDistribution();







};










#endif /* JUMPDISTRIBUTION_HPP_ */
