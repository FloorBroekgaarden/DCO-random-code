//
//  uniformDistribution.cpp


#include "uniformDistribution.h"

double uniformDistribution(const gsl_rng *r, double min, double max){
    
    // Generate numbers uniformly between min and max (tested in unifromDistribution project)
    
    double range = max - min;
    
    return min + gsl_rng_uniform(r)*range;
    
}
