#include <vector>
#include <iostream>

#include "JumpDistribution.hpp"

#define cout std::cout
#define endl std::endl


JumpDistribution::JumpDistribution(double intv, int bs,
    double r1, double r2,  ContinuousTimeMonteCarlo * MC)
{
 
interval = intv;
bins = bs;
range[0] = r1;
range[1] = r2;
MonteCarloObjectPtr = MC;

/* Check that the values input externally produce a histogram
 * and throw exception if not */
if (range[0] > range[1] || bins < 1 || interval < 0){
    cout << "Histogram bins not well defined" << endl << endl;
    throw "Histgram bins not well defined!";
}

binWidth = (range[1]-range[0])/bins;

midPointValues.resize(bins);
frequency.resize(bins);

for (int i=0; i<bins; ++i){
   midPointValues[i] = range[0] + (i+0.5)*binWidth;
}

makeDistribution();


}



void JumpDistribution::makeDistribution(){

double time = interval;
int jump = 0;

/* Loop over data until either the interval moves past
 * the last time, or the last data point is reached */
while (time < MonteCarloObjectPtr->jumpTimes.back() &&
    jump < static_cast<int>( MonteCarloObjectPtr->jumpTimes.size() ) ){

    int jumps = 0;

    //Sum jumps in time interval
    while (MonteCarloObjectPtr->jumpTimes[jump] < time){
            jumps += 1;
            jump += 1;
    }

    //Move to next time interval
    time += interval;

    //Locate the bin...
    int bin = static_cast<int> (jumps - range[0])/binWidth;

    //...and add one to the samples for that bin
    if (bin < bins && bin >=0){
        frequency[bin] += 1;
    }
}


}





















