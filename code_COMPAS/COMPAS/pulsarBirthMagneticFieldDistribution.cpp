#include "pulsarBirthMagneticFieldDistribution.h"

double pulsarBirthMagneticFieldDistribution(programOptions options, const gsl_rng *r){
	/*
	Return log10 of the magnetic field (in G) for a pulsar at birth

	Parameters
	-----------
	options : programOptions
		User specified program options

	Returns
	--------
	log10B : double
		log10 of the birth magnetic field in G

	*/

	double log10B = 0;

	if(boost::iequals(options.pulsarBirthMagneticFieldDistributionString, "ZERO")){
		log10B = 0;
    }
    else if(boost::iequals(options.pulsarBirthMagneticFieldDistributionString, "FIXED")){
    	// Set to a fixed constant value
    	std::cout << "Fixed value not yet implemented. Implement as options.pulsarBirthMagneticFieldFixedValue" << std::endl;
    	log10B = 0;
    }
    else if(boost::iequals(options.pulsarBirthMagneticFieldDistributionString, "FLATINLOG")){
    	// Flat in the log distribution from Oslowski et al 2011 https://arxiv.org/abs/0903.3538 (log10B0min = , log10B0max = )

    	double u = gsl_rng_uniform (r);
    	log10B = options.pulsarBirthMagneticFieldDistributionMin + (u * (options.pulsarBirthMagneticFieldDistributionMax - options.pulsarBirthMagneticFieldDistributionMin));
    }
    else if(boost::iequals(options.pulsarBirthMagneticFieldDistributionString, "UNIFORM")){
        // Flat distribution used in Kiel et al 2008 https://arxiv.org/abs/0805.0059 (log10B0min = 11, log10B0max = 13.5 see section 3.4 and Table 1.)

        double u = gsl_rng_uniform (r);

        double B = pow(10, options.pulsarBirthMagneticFieldDistributionMin) + (u * (pow(10,options.pulsarBirthMagneticFieldDistributionMax) - pow(10,options.pulsarBirthMagneticFieldDistributionMin)));
        
        log10B = log10(B);
    }
    else if(boost::iequals(options.pulsarBirthMagneticFieldDistributionString, "LOGNORMAL")){
    	// log normal distribution from Faucher-Giguere and Kaspi 2006 https://arxiv.org/abs/astro-ph/0512585
        
    	// Values hard-coded for now, can make them options if necessary
    	double pulsarBirthMagneticFieldDistributionFaucherGiguereKaspi2006Mean = 12.65;
    	double pulsarBirthMagneticFieldDistributionFaucherGiguereKaspi2006Std = 0.55;

    	log10B = gsl_ran_gaussian (r, pulsarBirthMagneticFieldDistributionFaucherGiguereKaspi2006Std) + pulsarBirthMagneticFieldDistributionFaucherGiguereKaspi2006Mean;
        
    }
    else{
    	std::cerr << "Invalid pulsarBirthMagneticFieldDistribution provided " << std::endl;
    }

    return log10B;

}
