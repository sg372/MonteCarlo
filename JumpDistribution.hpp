/* Class to produce a histogram for the number of quantum jumps
 * in a time interval, using data from a continuous-time Monte
 * Carlo simulation.  
 * 
 * Objects are constructed with
 * double intv : the size of the time interval to sample.
 * double bs : the number of bins for the histogram.
 * double r1, r2 : the range specified for the histogram.
 * double 
 */

#ifndef JUMPDISTRIBUTION_HPP_

#define JUMPDISTRIBUTION_HPP_

#include <vector>
#include "defs.hpp"
#include "ContinuousTimeMonteCarlo.hpp"



class JumpDistribution{


public:

JumpDistribution( double, int, double, double, ContinuousTimeMonteCarlo *);


double interval;
int bins;
double range[2];
double binWidth;

ContinuousTimeMonteCarlo * MonteCarloObjectPtr;

std::vector<double> midPointValues;
std::vector<int> frequency;




private:

void makeDistribution();


};



#endif /* JUMPDISTRIBUTION_HPP_ */
