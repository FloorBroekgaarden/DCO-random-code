
//  eDistribution.cpp


#include "eDistribution.h"

double eDistribution(programOptions const &options, const gsl_rng *r, unsigned long randomSeed){

    if(boost::iequals(options.eccentricityDistribution, "ZERO")){
        // All systems have are initially circular i.e. have zero eccentricity 
        return 0;
    }
    else if(boost::iequals(options.eccentricityDistribution, "FIXED")){
        // All systems have same initial eccentricity
        std::cout << "Fixed eccentricity not yet implemented" << std::endl;
        return 0;
    }
    else if(boost::iequals(options.eccentricityDistribution, "FLAT")){
        double ePower = 0;
        double eMax   = options.eccentricityDistributionMax;
        double eMin   = options.eccentricityDistributionMin;

        return inverse_sampling_from_power_law(r, ePower, eMax, eMin);
    }
    else if(boost::iequals(options.eccentricityDistribution, "THERMALISED") or boost::iequals(options.eccentricityDistribution, "THERMAL")){
        // Thermal eccentricity distribution p(e) = 2e
        double ePower = 1;
        double eMax   = options.eccentricityDistributionMax;
        double eMin   = options.eccentricityDistributionMin;

        return inverse_sampling_from_power_law(r, ePower, eMax, eMin);
    }
    else if(boost::iequals(options.eccentricityDistribution, "GELLER+2013")){
        // M35 eccentricity distribution from Geller, Hurley and Mathieu 2013
        // Gaussian with mean 0.38 and sigma 0.23
        // http://iopscience.iop.org/article/10.1088/0004-6256/145/1/8/pdf
        // Sampling function taken from binpop.f in NBODY6
        
        double e = -1.0;

        while(e < 0.0 or e > 1.0){

            double x1 = gsl_rng_uniform(r);              // Draw a random number between 0 and 1
            double x2 = gsl_rng_uniform(r);              // Draw a random number between 0 and 1

            e = 0.23 * sqrt(-2.0 * log(x1)) * cos(2.0 * pi * x2) + 0.38; //0.23 * SQRT(-2.D0*LOG(x1))*COS(TWOPI*x2)+0.38

        }

        return e;
    }
    else if(boost::iequals(options.eccentricityDistribution, "DuquennoyMayor1991")){
        // Eccentricity distribution from Duquennoy & Mayor (1991)
        // http://adsabs.harvard.edu/abs/1991A%26A...248..485D
        // Sampling function taken from binpop.f in NBODY6

        double e = -1.0;

        while(e < 0.0 or e > 1.0){

            double x1 = gsl_rng_uniform(r);              // Draw a random number between 0 and 1
            double x2 = gsl_rng_uniform(r);              // Draw a random number between 0 and 1

            e = 0.15 * sqrt(-2.0 * log(x1)) * cos(2.0 * pi * x2) + 0.3; //0.15*SQRT(-2.D0*LOG(x1))*COS(TWOPI*x2)+0.3

        }

        return e;
    }
    else if(boost::iequals(options.eccentricityDistribution, "SANA2012")){
        // Sana et al 2012 (http://science.sciencemag.org/content/sci/337/6093/444.full.pdf) distribution of eccentricities.
        // Taken from table S3 in http://science.sciencemag.org/content/sci/suppl/2012/07/25/337.6093.444.DC1/1223344.Sana.SM.pdf
        // See also de Mink and Belczynski 2015 http://arxiv.org/pdf/1506.03573v2.pdf
        double ePower = -0.42;
        double uncertainty = 0.17;
        double eMax   = options.eccentricityDistributionMax;
        double eMin   = options.eccentricityDistributionMin;

        return inverse_sampling_from_power_law(r, ePower, eMax, eMin);
    }
    else if(boost::iequals(options.eccentricityDistribution, "IMPORTANCE")){
  
        std::cout << "Not yet implemented. Do importance sampling here" << std::endl;
        
        return 0.0;
    }
    else{
        // User has specified an invalid distribution - use default
        std::cerr << "Invalid e distribution - using default" << std::endl;
        return 0;
    }
}
