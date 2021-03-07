//
//  kickVelocityDistribution.cpp
//
//  Implementation of functions in kickAndOrbitalPlane.cpp
//
//
//

#include "kickVelocityDistribution.h"

double kickVelocityDistribution(programOptions options, const gsl_rng *r, double kickVelocityDistributionSigma, double COCoreMass, unsigned long m_randomSeed, double star_kick_magnitude_random_number, double Mejecta, double Mremnant){
    /*
    Draw a kick velocity from the chosen distribution

    Parameters
    ------------
    options : programOptions
        User specified programOptions
    r : gsl_rng
        Random number generator
    kickVelocityDistributionSigma : double
        Velocity dispersion in km s^-1
    COCoreMass : double
        Carbon Oxygen core mass of exploding star in solar masses
    m_randomSeed : unsigned long
        Random seed
    star_kick_magnitude_random_number : double
        random number used to draw kick
    Mejecta : double
        Change in mass of exploding star (= ejecta mass) in solar masses
    Mremnant : double
        Mass of remnant in solar masses

    Returns
    --------
    vK : double
        Kick velocity in km s^-1
    */
    
	bool	debbugging = false;
	double	kickVelocity = NEVER_SET;
	double 	kickScalingFactor = options.kickScalingFactor;

    if (boost::iequals(options.kickVelocityDistribution, "Maxwellian") or boost::iequals(options.kickVelocityDistribution, "Maxwell")) {
        return kickVelocityDistributionMaxwell(r, kickVelocityDistributionSigma, star_kick_magnitude_random_number)/kickScalingFactor;
    }
    else if(boost::iequals(options.kickVelocityDistribution, "Flat")){
        return kickVelocityDistributionFlat(r, options.kickVelocityDistributionMaximum, star_kick_magnitude_random_number)/kickScalingFactor;
    }
    else if(boost::iequals(options.kickVelocityDistribution, "Zero")){
        return kickVelocityDistribution0()/kickScalingFactor; 
    }
    else if(boost::iequals(options.kickVelocityDistribution, "Fixed")){
        return kickVelocityDistributionSigma/kickScalingFactor; 
    }
    else if(boost::iequals(options.kickVelocityDistribution, "BrayEldridge")){
        return kickVelocityBrayEldridge(Mejecta, Mremnant, ALPHA_BRAY_ELDRIDGE, BETA_BRAY_ELDRIDGE)/kickScalingFactor; // constants defined in constants.h
    }
    else if(boost::iequals(options.kickVelocityDistribution, "muller2016")){
		kickVelocity = MullerRemnantKick(COCoreMass, m_randomSeed)/kickScalingFactor; //defined in SNMuller2016.cpp
		if(debbugging){
			std::cout << "kickVelocityDistribution: muller2016" << std::endl;
			std::cout << "Vkick:\t" << kickVelocity << std::endl;
		}
        return kickVelocity;
    }
    else if(boost::iequals(options.kickVelocityDistribution, "muller2016Maxwellian")){

        double mullerVKick = MullerRemnantKick(COCoreMass, m_randomSeed); // Defined in SNMuller2016.cpp, returns a 3D kick
		double mullerSigma = mullerVKick/sqrt(3.0);
		double kickVelocity = kickVelocityDistributionMaxwell(r, mullerSigma, star_kick_magnitude_random_number)/kickScalingFactor; // Should be fed a 1D sigma

		if(debbugging){
			std::cout << "kickVelocityDistribution: muller2016" << std::endl;
			std::cout << "Vkick, sigma, kickScalingFactor, VkickDraw:\t" << mullerVKick << ", " << mullerSigma << " ," << kickScalingFactor << " , "<< kickVelocity << std::endl;
		}
        return kickVelocity;
    }
    else{

        std::cerr << "ERROR: You did not specify a valid kick distribution. Using default." << std::endl;
        return kickVelocityDistributionMaxwell(r, kickVelocityDistributionSigma, star_kick_magnitude_random_number)/kickScalingFactor;

    }

}
