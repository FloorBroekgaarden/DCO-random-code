//
//  metallicityDistribution.cpp

#include "metallicityDistribution.h"

double metallicityDistribution(programOptions options, const gsl_rng *r){
    
    if(options.fixedMetallicity){
        
        return options.metallicity;
        
    }
    else{
        
        //std::cout << "Error in setting Metallicity of stars. Setting to solar by default." << std::endl;
        
        return Zsol; // Use solar metallicity by default.
        
    }
    
}