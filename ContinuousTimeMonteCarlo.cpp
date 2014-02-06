#include <random>
#include <iostream>
#include <vector>
#include <iterator>
#include <fstream>
#include <string>
#include <ctime>

#include "ContinuousTimeMonteCarlo.hpp"

#define cout std::cout
#define endl std::endl

//Set smallest jump rate below which state is disconnected
#define EPS 1e-10


//Constructor to perform a fixed number of time steps
ContinuousTimeMonteCarlo::ContinuousTimeMonteCarlo(int lS,
    int tS, Liouvillian * stochasticMatrix)
{

    initialState = lS;
    timeSteps = tS;
    state = initialState;
    time = 0.0;
    
    this->basisSize = stochasticMatrix->basisSize;
    this->W = &stochasticMatrix->W;
    
    //Perform required number of time steps
    for (int j=0; j<timeSteps; ++j){
        try{
            timeStep();
        }
        catch (int error) {
            cout << "Error: state " << error << " is disconnected never jump";
            cout << " and will never jump" << endl;
            break;
        }
    }

}


//Constructor to iterate until a desired time has been reached
ContinuousTimeMonteCarlo::ContinuousTimeMonteCarlo(int lS,
    double mT, Liouvillian * stochasticMatrix)
{

    initialState = lS;
    maxTime = mT;
    state = initialState;
    time = 0.0;
    
    this->basisSize = stochasticMatrix->basisSize;
    this->W = &stochasticMatrix->W;
    
    //Step until required time has been reached
    while (time < maxTime){
        try{
            timeStep();
        }
        catch (int error) {
            cout << "Error: state " << error << " is disconnected never jump";
            cout << " and will never jump" << endl;
            break;
        }
    }

}



//Constructor to read in object data from a file
ContinuousTimeMonteCarlo::ContinuousTimeMonteCarlo(std::string fname)
{
    readFromFile(fname);

    initialState = *jumpStates.begin();
    maxTime = *jumpTimes.end();
    timeSteps = jumpTimes.size();

}




/* Perform a single time step and add the time and new state
 * to the stored values jumpTimes and jumpStates.
 * This routine will throw an int int corresponding to an
 * infinitely-lived state, should one be prepared initially */
void ContinuousTimeMonteCarlo::timeStep()
{
//static std::mt19937 generator(std::time(0));
static std::default_random_engine generator(std::time(0));
static std::uniform_real_distribution<double> flatDistribution(0.0,1.0);

std::vector<double> cumul(basisSize);	//table of cumulative rates 

    double cu=0.0;      //for cumulative probability weight
    
    for (int i=0; i<basisSize; ++i){
	    if(i != state){
	        cu += (*W)(i,state);
	        cumul[i] = cu;
        }
	    else if(i == state){	//Exclude current state by not incrementing
	        cumul[i] = cu;	//cumulative rate for this state
	    }
    }

//Select random number in [0, max(cu))
double u = flatDistribution(generator) * cu;	

//Select new state from table cumul (binary search)
state = std::distance(cumul.begin(),
                          std::lower_bound(cumul.begin(), cumul.end(), u));
    
double random_number = 1.0-flatDistribution(generator);  //distributed in (0.0,1.0]


//Check state is not isolated and a jump will happen, else throw
if ( -(*W)(state,state) < EPS ){
    throw state;
}

//Increment time
time += log(random_number)/ (*W)(state,state);

//Record time and corresponding new state;
jumpTimes.push_back(time);
jumpStates.push_back(state);										 
					
}



void ContinuousTimeMonteCarlo::writeToFile(std::string fname){

std::ofstream file (fname);

if (file.is_open() ){
    
    file << "#ContinuousTimeMonteCarlo" << endl;
    file << "#" << basisSize << endl;
    file << "#" << jumpTimes.size() << endl;


    std::vector<double>::iterator itT = jumpTimes.begin();
    std::vector<int>::iterator itS = jumpStates.begin();

    for ( ; itT < jumpTimes.end(), itS < jumpStates.end(); ++itT, ++itS){

        file << *itT << " "<< *itS << endl;

    }
}
else {
    cout << "File open unsuccessful. ";
    cout << " Monte Carlo data not written to file." << endl;
}

file.close(); 

}




void ContinuousTimeMonteCarlo::readFromFile(std::string fname){

char dummyhash;
int sizeOfVectors;
std::string checkstring;
std::ifstream file (fname);


if (file.is_open() ){

    file >> checkstring;

    if ( checkstring.compare("#ContinuousTimeMonteCarlo")==0 ){
        
        file >> dummyhash >> basisSize;
        file >> dummyhash >> sizeOfVectors;

        jumpTimes.resize(sizeOfVectors);
        jumpStates.resize(sizeOfVectors);
        
        std::vector<double>::iterator itT = jumpTimes.begin();
        std::vector<int>::iterator itS = jumpStates.begin();

        for ( ; itT < jumpTimes.end(), itS < jumpStates.end(); ++itT, ++itS){

            file >> *itT >> *itS;

        }

    }
    else cout << "File does not contain ContinousTimeMonteCarlo data" << endl;

}
else cout << "File open unsuccessful" << endl;




}







