//
//  qDistribution.cpp
//
//  This function returns a mass ratio value q

#include "qDistribution.h"

double qDistribution(programOptions const &options, const gsl_rng *r,  AISvariables const &aisvariables){ 

    // draw from priors if  DrawingFromAISdistributions = false
    if(!aisvariables.DrawingFromAISdistributions){

        if(boost::iequals(options.massRatioDistribution, "FLAT")){
            double qDistributionPower = 0;
            double qDistributionMin   = options.massRatioDistributionMin;
            double qDistributionMax   = options.massRatioDistributionMax;
            return inverse_sampling_from_power_law(r, qDistributionPower, qDistributionMax, qDistributionMin);
        }
        else if(boost::iequals(options.eccentricityDistribution, "DuquennoyMayor1991")){
            // Mass ratio distribution from Duquennoy & Mayor (1991)
            // http://adsabs.harvard.edu/abs/1991A%26A...248..485D

            double q = -1.0;

            while(q < 0.0 or q > 1.0){

                double x1 = gsl_rng_uniform(r);              // Draw a random number between 0 and 1
                double x2 = gsl_rng_uniform(r);              // Draw a random number between 0 and 1

                q = 0.42 * sqrt(-2.0 * log(x1)) * cos(2.0 * pi * x2) + 0.23;

            }

            return q;
        }
        else if(boost::iequals(options.massRatioDistribution, "SANA2012")){
            // Sana et al 2012 (http://science.sciencemag.org/content/sci/337/6093/444.full.pdf) distribution of eccentricities.
            // Taken from table S3 in http://science.sciencemag.org/content/sci/suppl/2012/07/25/337.6093.444.DC1/1223344.Sana.SM.pdf
            // See also de Mink and Belczynski 2015 http://arxiv.org/pdf/1506.03573v2.pdf
            double qDistributionPower = -0.1;
            double uncertainty        = 0.6;
            double qDistributionMin   = options.massRatioDistributionMin;           // de Mink and Belczynski use 0.1
            double qDistributionMax   = options.massRatioDistributionMax;           // de Mink and Belczynski use 1.0
            return inverse_sampling_from_power_law(r, qDistributionPower, qDistributionMax, qDistributionMin);
        }
        else{
            // User has entered an invalid - reset to default
            
            std::cout << "Invalid q distribution entered - using default." << std::endl;
            
            double qDistributionPower = 0;
            double qDistributionMin   = 0;
            double qDistributionMax   = 1;
            return inverse_sampling_from_power_law(r, qDistributionPower, qDistributionMax, qDistributionMin);
        }
    }

    // draw from AIS distributions if  DrawingFromAISdistributions = true
    else{   
        // Mass ratio distribution from Adaptive Importance Sampling v1 from Broekgaarden et al. (in prep 2018)
        // mu_q and cov_q consist of arrays with resp the gaussian means and covariances of the gaussians of a which where determined in step 2 of AIS 
        
        // Now draw a random q from the random Gaussian chosen with RandomGaussianDraw
        double qMean =  aisvariables.mu_q[aisvariables.RandomGaussianDraw];            // The mean of the RandomGaussianDraw-th Gaussian 
        double qStd  =  aisvariables.cov_q[aisvariables.RandomGaussianDraw];           // The cov of the RandomGaussianDraw-th Gaussin
        double q = gsl_ran_gaussian (r, qStd) + qMean;  // draw random nr from Gaussian 
        return   q; 
    }

}

////////////////////////////////////////////////////////////////////////////////////
//                                  NOTES
////////////////////////////////////////////////////////////////////////////////////


