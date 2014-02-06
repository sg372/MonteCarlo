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

while (time < MonteCarloObjectPtr->jumpTimes.back() &&
    jump < MonteCarloObjectPtr->jumpTimes.size()){

    int jumps = 0;

    while (MonteCarloObjectPtr->jumpTimes[jump] < time){
            jumps += 1;
            jump += 1;
    }

    time += interval;

    int bin = static_cast<int> (jumps - range[0])/binWidth;

    if (bin < bins && bin >=0){
        frequency[bin] += 1;
    }
}


}





















