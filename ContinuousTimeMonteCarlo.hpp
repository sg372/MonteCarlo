#ifndef CONTINUOUSTIMEMONTECARLO_HPP_
#define CONTINUOUSTIMEMONTECARLO_HPP_

#include "defs.hpp"
#include "Liouvillian.hpp"
#include <vector>

class ContinuousTimeMonteCarlo
{

public:

ContinuousTimeMonteCarlo(unsigned, unsigned, Liouvillian *);

ContinuousTimeMonteCarlo(unsigned, double, Liouvillian *);

ContinuousTimeMonteCarlo(std::string);

unsigned initialState;
unsigned timeSteps;
double maxTime;

unsigned basisSize;


std::vector<double> jumpTimes;
std::vector<unsigned> jumpStates;


void writeToFile(std::string);

private:

DoubleMatrix *W;

void timeStep();



double time;
unsigned state;

void readFromFile(std::string);

};


#endif /* CONTINUOUSTIMEMONTECARLO_HPP_ */
