



#include <random>





void BKLstep(state, time)
{
std::default_random_engine generator;
std::flat_distribution<double> flatDistribution;

std::vector<double> cumul(basisSize);	//table of cumulative rates
integer :: i,j,k 
real(dp) :: u  //for random number for new state selection
real(dp) :: cu //for cumulative probability weight

    cu=0.0
    for (i=0; i<basisSize; ++i){
	if(i != state){
	    cu += W(i,state)
	    cumul[i] = cu
        }
	else if(i == state){	//Exclude current state by not incrementing
	    cumul[i] = cu	//cumulative rate for this state
	}
    }
    //Select random number in [0, max(cu))
    u = flatDistribution(generator) * cu;	

    //Select new state from table cumul (binary search)
    state = std::lower_bound(cumul.front, cumul.back, u)

    //Add on time increment
    time += log(flatDistribution(generator))/W(state,state)
											 
					
}





unsigned binarySearch(std::vector<double> values, double valToFind)
{
    // Finds the lower bound in at most log(last - first) + 1 comparisons
    double i = ;
    return i; // found
}





