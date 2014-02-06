
#include <iostream>
#include <fstream>
#include <vector>
#include "defs.hpp"
#include "FermiBasis.hpp"
#include "Hamiltonian.hpp"
#include "Liouvillian.hpp"
#include "ContinuousTimeMonteCarlo.hpp"
#include "JumpDistribution.hpp"

#define cout std::cout
#define endl std::endl


int main(){



Liouvillian * Liou = new Liouvillian(0.0, 1.0, 15, 2, 1.0, 1.0);

/*
cout << "the Liouvillian is:" << endl;
cout << Liou->W << endl << endl;
*/

ContinuousTimeMonteCarlo * MC = new ContinuousTimeMonteCarlo(1,100000000,Liou);


JumpDistribution * Dist2 = new JumpDistribution(50.0, 40, 0, 150, MC);

std::vector<double>::iterator itV = Dist2->midPointValues.begin();
std::vector<int>::iterator itF = Dist2->frequency.begin();
for ( ; itF < Dist2->frequency.end(), itV<Dist2->midPointValues.end() ; 
     ++itF, ++itV){
    cout << *itV << "  " << *itF << endl;
}


std::ofstream histogramFile ("histogram5.dat");

std::vector<int>::iterator itFreq = Dist2->frequency.begin();
std::vector<double>::iterator itMpv = Dist2->midPointValues.begin();
for ( ; itFreq < Dist2->frequency.end(), itMpv< Dist2->midPointValues.end();
         ++itFreq, ++itMpv){
    histogramFile << *itMpv << " " << *itFreq << endl;
}

histogramFile.close();


delete Dist2;
delete MC;
delete Liou;





return 0;

}
