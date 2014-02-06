#ifndef CONTINUOUSTIMEMONTECARLO_HPP_
#define CONTINUOUSTIMEMONTECARLO_HPP_

#include "defs.hpp"
#include "Liouvillian.hpp"
#include <vector>

class ContinuousTimeMonteCarlo
{

public:

ContinuousTimeMonteCarlo(int, int, Liouvillian *);

ContinuousTimeMonteCarlo(int, double, Liouvillian *);

ContinuousTimeMonteCarlo(std::string);

int initialState;
int timeSteps;
double maxTime;

int basisSize;


std::vector<double> jumpTimes;
std::vector<int> jumpStates;


void writeToFile(std::string);

private:

DoubleMatrix *W;

void timeStep();



double time;
int state;

void readFromFile(std::string);

};


#endif /* CONTINUOUSTIMEMONTECARLO_HPP_ */
