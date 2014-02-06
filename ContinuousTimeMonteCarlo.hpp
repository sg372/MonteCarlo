/* Performs continous time (BKL) Monte Carlo for a given
 * Liouvillian passed by a pointer.
 * 
 * Two constructors exist allowing the passing of
 * int lS : the initial Liouville state
 * <<<EITHER:
 * int tS : the number of time steps through continuous time
 * OR:
 * double mT : the maximum time, beyond which the routine stops
 * >>>
 * Liouvillian * : pointer to Liouvillian object to be simulated
 *
 * A third constructor is provided creat the object from data
 * in from a file with
 * std::string fname : a string with the filename
 * (provided in c++11 standard)
 *
 * The data corresponding to the times of jumps and the states
 * jumped to are held in members of type std::vector<int> for 
 * the jumpStates and std::vector<double>
 * 
 */


#ifndef CONTINUOUSTIMEMONTECARLO_HPP_
#define CONTINUOUSTIMEMONTECARLO_HPP_

#include "defs.hpp"
#include "Liouvillian.hpp"
#include <vector>
#include <string>

class ContinuousTimeMonteCarlo
{

public:

//Constructors:
ContinuousTimeMonteCarlo(int, int, Liouvillian *);

ContinuousTimeMonteCarlo(int, double, Liouvillian *);

ContinuousTimeMonteCarlo(std::string);


//Parameters to specify the trajectory
int initialState;
int timeSteps;
double maxTime;

//Size of space passed from Liouvillian
int basisSize;

//Data from simulation
std::vector<double> jumpTimes;
std::vector<int> jumpStates;

//To write data out to file
void writeToFile(std::string);

private:

//Ptr to the Liouvillian
DoubleMatrix *W;

//Step in continous time
void timeStep();

//Current state of the system and time
double time;
int state;

//For constructor to read from file
void readFromFile(std::string);

};


#endif /* CONTINUOUSTIMEMONTECARLO_HPP_ */
