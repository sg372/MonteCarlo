
#include <iostream>
#include <fstream>
#include <vector>

#include "defs.hpp"
#include "FermiBasis.hpp"
#include "Hamiltonian.hpp"
#include "Liouvillian.hpp"
#include "ContinuousTimeMonteCarlo.hpp"
#include "JumpDistribution.hpp"

//Parameters to specify the quantum model
#define SITES           15
#define PARTICLES       2
#define TUNNELLING      1.0
#define DISORDER_STR_1  1.0
#define DISORDER_STR_2  20.0

//Parameters to specify the master equation when coupling the model to a bath
#define COUPLING_STR    1.0

//Parameters to specify the length of the trajectory
#define TIME_STEPS      10000000
#define INIT_Q_STATE    1

//Parameters to specify the histogramming of jump frequencies
#define NO_OF_BINS      40
#define SAMPLING_TIME   50.0
#define LOWER_BOUND     0.0
#define UPPER_BOUND     150.0

//Filenames for output of histograms
#define FNAME_1         "histogram_weak_disorder.dat"
#define FNAME_2         "histogram_strong_disorder.dat"



int main(){


//Construct a weakly-disordered system and evolve it and histogram results
Liouvillian * Liou = new Liouvillian(COUPLING_STR, SITES, PARTICLES,
        TUNNELLING, DISORDER_STR_1);

ContinuousTimeMonteCarlo * MC = new ContinuousTimeMonteCarlo(INIT_Q_STATE,
        TIME_STEPS,Liou);

JumpDistribution * Dist = new JumpDistribution(SAMPLING_TIME, NO_OF_BINS, 
        LOWER_BOUND, UPPER_BOUND, MC);


//Write the histogram to a file and to stdout
std::cout << "Histogram data for the case of weak disorder:" << std::endl;
std::cout << "JUMPS   FREQUENCY" << std::endl;

std::ofstream histogramFile (FNAME_1);

std::vector<int>::iterator itFreq = Dist->frequency.begin();
std::vector<double>::iterator itMpv = Dist->midPointValues.begin();

for ( ; itFreq < Dist->frequency.end(), itMpv< Dist->midPointValues.end();
         ++itFreq, ++itMpv){

    histogramFile << *itMpv << " " << *itFreq << std::endl;
    std::cout  << *itMpv << " " << *itFreq << std::endl;

}

histogramFile.close();

delete Dist;
delete MC;
delete Liou;


//Construct a strongly-disordered system and evolve it and histogram results
Liou = new Liouvillian(COUPLING_STR, SITES, PARTICLES,
        TUNNELLING, DISORDER_STR_2);

MC = new ContinuousTimeMonteCarlo(INIT_Q_STATE,
        TIME_STEPS,Liou);

Dist = new JumpDistribution(SAMPLING_TIME, NO_OF_BINS, 
        LOWER_BOUND, UPPER_BOUND, MC);

//Write the histogram to a file and to stdout
std::cout << "Histogram data for the case of strong disorder:" << std::endl;
std::cout << "JUMPS   FREQUENCY" << std::endl;

histogramFile.open(FNAME_2);

itFreq = Dist->frequency.begin();
itMpv = Dist->midPointValues.begin();

for ( ; itFreq < Dist->frequency.end(), itMpv< Dist->midPointValues.end();
         ++itFreq, ++itMpv){

    histogramFile << *itMpv << " " << *itFreq << std::endl;
    std::cout  << *itMpv << " " << *itFreq << std::endl;

}

histogramFile.close();

delete Dist;
delete MC;
delete Liou;


return 0;

}
