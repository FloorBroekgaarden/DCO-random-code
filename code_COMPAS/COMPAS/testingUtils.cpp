//
//  testingUtils.cpp
//  
//  Created by Jim Barrett on 12/01/2017.
//
//	utilities to help with unit tests
//

#include "testingUtils.h"

#include "programOptions.h"
#include "initialMassDistribution.h"
#include "star.h"
#include "BinaryStar.h"

#include <boost/math/constants/constants.hpp>

programOptions* getDefaultCompasOptions() {
	//Sets the COMPAS options to their default options as defined in programOptions.cpp

    static bool firstTime = true;

	static programOptions options;

    if (firstTime)
    {
        options.setToFiducialValues();
        firstTime = false;
	}

	return &options;

}

gsl_rng* getGslRng() {
	//Sets up and returns the gsl_rng object for random numbers

    static bool firstTime = true;    

	static gsl_rng* r;
	
	if (firstTime)
	{
		gsl_rng_env_setup();
	
		gsl_rng_default_seed = time(0);
	
		r = gsl_rng_alloc (gsl_rng_default);

        firstTime = false;
	}
	
	return r;

}
// Need to update constructor to work with Floor's AIS  (Floor - 26/05/2018) // FLOOR 
// Star generateRandomStar() {
// 	//Generates a random star

// 	programOptions * options = getDefaultCompasOptions();
	
// 	options->randomSeed = makeRandomSeed();
	
// 	gsl_rng * r = getGslRng();

// 	std::vector<double> mu_M1,  cov_M1; // For Adaptive Importance Sampling

// 	int RandomGaussianDraw = 99999; // For Adaptive importance Sampling

// 	bool DrawingFromAISdistributions = false;

// 	double mass = initialMassDistribution(*options, r, DrawingFromAISdistributions, RandomGaussianDraw, aisvariables);
	
// 	double metallicity = options->metallicity;
	
// 	return Star(mass,metallicity, *options, r);
	
//}


// BinaryStar generateRandomBinary( ) { // Floor 2 
// 	//Generates a random binary
// 	std::vector<double> mu_M1,  mu_loga,  mu_q,  cov_M1,  cov_loga,  cov_q;
// 	programOptions * options = getDefaultCompasOptions();
	
// 	options->randomSeed = makeRandomSeed();
// //    options->randomSeed = 1229114076;
	
// 	gsl_rng * r = getGslRng();

// 	return BinaryStar(r, *options, mu_M1, mu_loga, mu_q, cov_M1, cov_loga, cov_q );

// }

double getPi() {
    
    static bool firstTime = true;
    static double pi;
    
    if (firstTime) pi = boost::math::constants::pi<double>();
    
    return pi;
    
}

unsigned long makeRandomSeed()
{
    //makes a random seed

    gsl_rng* r = getGslRng(); 
    
    double u = gsl_rng_uniform(r);
    
    long seed = std::floor(u * 9999999999);
    
    return seed;
    
}
