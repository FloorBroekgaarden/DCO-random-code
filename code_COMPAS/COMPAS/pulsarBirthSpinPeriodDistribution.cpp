#include "pulsarBirthSpinPeriodDistribution.h"

double pulsarBirthSpinPeriodDistribution(programOptions options, const gsl_rng *r){
	/*
	Return the spin period of a pulsar at birth in ms
    
	Parameters
	-----------
	options : programOptions
		User specified program options

	Returns
	--------
	Pspin : double
		Birth spin period of pulsar in ms
    
	*/

	double Pspin = 0;

	//std::cout << "spin distribution = " << options.pulsarBirthSpinPeriodDistributionString << std::endl;


	if(boost::iequals(options.pulsarBirthSpinPeriodDistributionString, "ZERO")){
		Pspin = 0;
    }
    else if(boost::iequals(options.pulsarBirthSpinPeriodDistributionString, "FIXED")){
    	// Set to a fixed constant value as used in default model in
        // Oslowski et al 2011 https://arxiv.org/abs/0903.3538
    //	std::cout << "Fixed value not yet implemented. Implement as options.pulsarBirthSpinPeriodFixedValue" << std::endl;
    	Pspin = 0;
    }
    else if(boost::iequals(options.pulsarBirthSpinPeriodDistributionString, "UNIFORM")){
        // Use uniform distribution between minimum and maximum value as in
        // Oslowski et al 2011 https://arxiv.org/abs/0903.3538 (default Pmin = and Pmax = )
        // and also Kiel et al 2008 https://arxiv.org/abs/0805.0059 (default Pmin = 10 ms and Pmax 100 ms, section 3.4)
        
        double u = gsl_rng_uniform (r);
        
        Pspin = options.pulsarBirthSpinPeriodDistributionMin + (u * (options.pulsarBirthSpinPeriodDistributionMax - options.pulsarBirthSpinPeriodDistributionMin));

	//std::cout << "in uniform" << std::endl;
	//std::cout << "u = " << u << std::endl;
	//std::cout << "Pspin = " << Pspin << std::endl;
    }
    else if(boost::iequals(options.pulsarBirthMagneticFieldDistributionString, "NORMAL")){
    	// Normal distribution from Faucher-Giguere and Kaspi 2006 https://arxiv.org/abs/astro-ph/0512585

    	// Values hard-coded for now, can make them options if necessary
    	double pulsarBirthSpinPeriodDistributionFaucherGiguereKaspi2006Mean = 300.0;
    	double pulsarBirthSpinPeriodDistributionFaucherGiguereKaspi2006Std = 150.0;

        Pspin = -1.0;

        while(Pspin < 0.0){

            Pspin = gsl_ran_gaussian(r, pulsarBirthSpinPeriodDistributionFaucherGiguereKaspi2006Std) + pulsarBirthSpinPeriodDistributionFaucherGiguereKaspi2006Mean;

        }        

    }
    else{
    	std::cerr << "Invalid pulsarBirthSpinPeriodDistribution provided " << std::endl;
    }

    return Pspin;
}
