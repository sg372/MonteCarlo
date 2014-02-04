
#include <iostream>
#include "FermiBasis.hpp"
#include "Hamiltonian.hpp"
#include "Liouvillian.hpp"
#include "ContinuousTimeMonteCarlo.hpp"
#include "JumpDistribution.hpp"

#define cout std::cout
#define endl std::endl


int main(){



Liouvillian * Liou = new Liouvillian(0.0, 1.0, 8, 4, 1.0, 1.0, 0.0);

/*
cout << "the Liouvillian is:" << endl;
cout << Liou->W << endl << endl;
*/

ContinuousTimeMonteCarlo * MC = new ContinuousTimeMonteCarlo(1,1000.0,Liou);

/*

cout << MC->state << endl;


for (unsigned i=0; i<MC->jumpTimes.size(); i++){
    cout << MC->jumpTimes[i] <<"  "<<MC->jumpStates[i] << endl;
}


std::vector<double>::iterator itT=MC->jumpTimes.begin();
std::vector<unsigned>::iterator itS=MC->jumpStates.begin();
for ( ; itT < MC->jumpTimes.end(), itS < MC->jumpStates.end(); ++itT, ++itS ){
    cout << *itT <<"  " << *itS <<endl;
}
*/

JumpDistribution * Dist = new JumpDistribution(100.0, 20, 150, 350, MC);

std::vector<double>::iterator itV=Dist->midPointValues.begin();
std::vector<unsigned>::iterator itF=Dist->frequency.begin();
for ( ; itF < Dist->frequency.end(), itV<Dist->midPointValues.end() ; 
    ++itF, ++itV){
        cout << *itV << "  " << *itF << endl;
}



MC->writeToFile("outputTest.dat");

cout << "basisSize  " << MC->basisSize << endl;




ContinuousTimeMonteCarlo * MC2 = new ContinuousTimeMonteCarlo ("outputTest.dat");


JumpDistribution * Dist2 = new JumpDistribution(100.0, 20, 150, 350, MC2);

itV=Dist2->midPointValues.begin();
itF=Dist2->frequency.begin();
for ( ; itF < Dist2->frequency.end(), itV<Dist2->midPointValues.end() ; 
    ++itF, ++itV){
        cout << *itV << "  " << *itF << endl;
}



MC->writeToFile("outputTest2.dat");

cout << "basisSize  " << MC2->basisSize << endl;




delete Dist2;
delete Dist;
delete MC;
delete Liou;




return 0;

}
