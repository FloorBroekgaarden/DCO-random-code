//
//  spinDistribution.cpp


#include "spinDistribution.h"

double spinDistribution(programOptions const &options, const gsl_rng *r){
    
    if(boost::iequals(options.spinDistribution, "ZERO")){
        
        return 0;
        
    }
    else if(boost::iequals(options.spinDistribution, "FLAT")){
        
        return uniformDistribution(r, options.spinDistributionMin, options.spinDistributionMax);
        
    }
    else if(boost::iequals(options.spinDistribution, "FIXED")){
        
        double spinValue = 0.7;         // Could have this set by the user, only use if they pick fixed
        
        if(spinValue < 0.0 or spinValue > 1.0){
            spinValue = 0.0;
            std::cerr << "ERROR: Spin not between 0 and 1. Setting spin = 0." << std::endl;
        }
        
        return spinValue;
    }
    else{
        // Error checking
        std::cerr << "ERROR : invalid spin distribution - using default 0" << std::endl;
        return 0;
    }
}

void assignMisalignments(double &theta1, double &theta2, double iPrime, programOptions const &options, const gsl_rng *r){
    
    // Assign misalignments to S1 and S2 based on assumptions for spin study
    
    // Set spin misalignment angles - will need some model switch here depending on what you are running - either as given by code for both primary and secondary, secondary given by code and primary = 0 or secondary given by code and primary isotropic (uniform in cos(theta))
    
    if(boost::iequals(options.spinAssumption, "useSame")){
        theta1 = iPrime;
        theta2 = iPrime;
    }
    else if(boost::iequals(options.spinAssumption, "bothAligned")){
        theta1 = 0;
        theta2 = 0;
    }
    else if(boost::iequals(options.spinAssumption, "secondaryMisaligned")){
        theta1 = 0;
        theta2 = iPrime;
    }
    else if(boost::iequals(options.spinAssumption, "bothIsotropic")){
        double costheta1 = gsl_rng_uniform(r)*2 - 1.0;              // Initial misalignment of star 1 uniform in cos(theta)
        double costheta2 = gsl_rng_uniform(r)*2 - 1.0;              // Initial misalignment of star 2 uniform in cos(theta)
        
        theta1 = acos(costheta1);
        theta2 = acos(costheta2);
    }
    else if(boost::iequals(options.spinAssumption, "gerosaInspired")){
        // Inspired by Gerosa et al 2013 who assume that after the 1st SN, the secondary is realigned but the primary is not. After the second SN, the secondary is therefore slightly misaligned whilst the primary is misaligned as a function of both kicks. To model this, we choose the primary misalignment isotropically, whilst the secondary is given from my code.
        double costheta1 = gsl_rng_uniform(r)*2 - 1.0;              // Initial misalignment of star 1 uniform in cos(theta)
        
        theta1 = acos(costheta1);
        theta2 = iPrime;
    }
    else{
        std::cerr << "Error: Something has gone wrong setting misalignment angles" << std::endl;
        theta1  = 0;                                            // Initial misalignment of star 1
        theta2  = 0;                                            // Initial misalignment of star 2
    }
    
}
