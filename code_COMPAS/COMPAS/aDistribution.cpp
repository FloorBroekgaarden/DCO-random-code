#include "aDistribution.h"

double aDistribution(programOptions const &options, const gsl_rng *r, double Mass1, double Mass2, AISvariables const &aisvariables){
    /*
     Draw semi-major axis from the distribution specified by the user
     
     Parameters
     -----------
     options : programOptions
        Contains parameters set by the user specifying their physical assumptions
     r : gsl_rng
        Random number generator
     
     Returns
     -------
     a : double
        Semi-major axis in AU

     */

    double power = 0.0;
    double semiMajorAxis = 0.0;

    // draw from priors if  DrawingFromAISdistributions = false
    if(!aisvariables.DrawingFromAISdistributions){

        if(boost::iequals(options.semiMajorAxisDistribution, "FLATINLOG")){         
            // Declare variables
            power = -1;
            return inverse_sampling_from_power_law(r, power, options.semiMajorAxisDistributionMax, options.semiMajorAxisDistributionMin);
            
        }
        else if(boost::iequals(options.semiMajorAxisDistribution, "DuquennoyMayor1991")){
            // Period distribution from Duquennoy & Mayor (1991)
            // http://adsabs.harvard.edu/abs/1991A%26A...248..485D
            // See also the period distribution (Figure 1) of M35 in Geller+ 2013 https://arxiv.org/abs/1210.1575
            // See also the period distribution (Figure 13) of local solar type binaries from Raghavan et al 2010 https://arxiv.org/abs/1007.0414
            // They have log-normal distribution with a mean of 5.03 and a standard deviation of 2.28, with a minimum period of around 0.1 days
            // Sampling function taken from binpop.f in NBODY6
            
            double semiMajorAxis = NEVER_SET;
            double logPeriodInDays = NEVER_SET;
            double periodInDays = NEVER_SET;
            double x1 = NEVER_SET;
            double x2 = NEVER_SET;

            // Make sure that the drawn semi-major axis is in the range specified by the user
            while(semiMajorAxis < options.semiMajorAxisDistributionMin or semiMajorAxis > options.semiMajorAxisDistributionMax){

                x1 = gsl_rng_uniform(r);              // Draw a random number between 0 and 1
                x2 = gsl_rng_uniform(r);              // Draw a random number between 0 and 1

                logPeriodInDays = 2.3 * sqrt(-2.0 * log(x1)) * cos(2.0 * pi * x2) + 4.8; // exp1 = 2.3*SQRT(-2.D0*LOG(x1))*COS(TWOPI*x2)+4.8

                periodInDays = pow(10.0, logPeriodInDays);

                // Convert period in days to semi-major axis in AU for the rest of the code -- to do this, you need to also pass this function the masses of the binary and then calculate this last, should be fine
                semiMajorAxis = convertPeriodInDaysToSemiMajorAxisInAU(Mass1, Mass2, periodInDays);

            }
        
            return semiMajorAxis;

        }
        else if(boost::iequals(options.semiMajorAxisDistribution, "CUSTOM")){
            return inverse_sampling_from_power_law(r, options.semiMajorAxisDistributionPower, options.semiMajorAxisDistributionMax, options.semiMajorAxisDistributionMin);
        }
        else if(boost::iequals(options.semiMajorAxisDistribution, "SANA2012")){
            // Sana et al 2012 (http://science.sciencemag.org/content/sci/337/6093/444.full.pdf) distribution of semi-major axes. Sana et al fit for the orbital period, which we sample in here, before returning the semi major axis
            // Taken from table S3 in http://science.sciencemag.org/content/sci/suppl/2012/07/25/337.6093.444.DC1/1223344.Sana.SM.pdf
            // See also de Mink and Belczynski 2015 http://arxiv.org/pdf/1506.03573v2.pdf
            power = -0.55;                  // Called \Pi in Sana et al 2012
            double uncertainty = 0.2;       // uncertainty from Sana et al 2012
            
            double logPeriodMin = 0.0;      // Variable for smallest initial log period
            double logPeriodMax = 0.0;      // Variable for largest initial log period
            
            if(options.periodDistributionMin > 1.0){
                logPeriodMin = log(options.periodDistributionMin);
            }
            else{
                std::cout << "Period distribution will blow up (1.0/0) for periods <= 1.0 day" << std::endl;
                logPeriodMin = 0.0;
            }
            
            if(options.periodDistributionMax > 1.0){
                logPeriodMax = log(options.periodDistributionMax);
            }
            else{
                std::cout << "Period distribution will blow up (1.0/0) for periods <= 1.0 day, logPeriod < 0.0" << std::endl;
                logPeriodMax = 0.0;
            }
            
            // Draw a period in days from their distribution
            double drawLogOfPeriodInDays = inverse_sampling_from_power_law(r, power, logPeriodMax, logPeriodMin);
            
            double periodInDays = exp(drawLogOfPeriodInDays);
            
            // Convert period in days to semi-major axis in AU for the rest of the code -- to do this, you need to also pass this function the masses of the binary and then calculate this last, should be fine
            semiMajorAxis = convertPeriodInDaysToSemiMajorAxisInAU(Mass1, Mass2, periodInDays);
        
            return semiMajorAxis;
        }
        else{
            // User has specified an invalid distribution
            std::cout << "Invalid choice of semi major axis distribution -- using the default instead." << std::endl;
            power = -1;
            double max   = 100;
            double min   = 0.5;
            return inverse_sampling_from_power_law(r, power, max, min);
        }
    }

    // draw from AIS distributions if  DrawingFromAISdistributions = true
    else{ 
        // Mass ratio distribution from Adaptive Importance Sampling v1 from Broekgaarden et al. (in prep 2018)
        // Function Returns a random semiMajorAxis drawn from one of the random gaussians defined bu vectors mu_loga & cov_loga
        // Notice-> the mu and cov are in log10(a) space so range is e.g. (-1,3) instead of (0.1, 1000). 
       
        // Now draw a random semiMajorAxis from the random Gaussian chosen with RandomGaussianDraw
        double log10semiMajorAxisMean =  aisvariables.mu_loga[aisvariables.RandomGaussianDraw];    // The mean of the RandomGaussianDraw-th Gaussian 
        double log10semiMajorAxisStd  =  aisvariables.cov_loga[aisvariables.RandomGaussianDraw];   // The cov of the RandomGaussianDraw-th Gaussin (this should not be sigma**2)
        double log10semiMajorAxis = gsl_ran_gaussian (r, log10semiMajorAxisStd) + log10semiMajorAxisMean; // draw random nr from Gaussian 

        // convert random nr in log10(a) space back to the initial space defined in pythonSubmit.py
        semiMajorAxis = pow(10, log10semiMajorAxis);
        return  semiMajorAxis;


    }

}


double convertPeriodInDaysToSemiMajorAxisInAU(double M1, double M2, double period){
    /*
     Converts a period in days to a semi-major axis in AU
     
     Parameters
     -----------
     M1 : double
        Primary mass in Msol
     M2 : double
        Secondary mass in Msol
     period : double
        Orbital period in days
     
     Returns
     -------
     a : double
        Semi-major axis in AU

     */
    double a_cubed_SI_top = G * ((M1 * Msol) + (M2 * Msol)) * period * period * day * day;
    double a_cubed_SI_bottom = 4.0 * pi * pi;
    double a_cubed_SI = a_cubed_SI_top / a_cubed_SI_bottom;
    double a_SI = pow(a_cubed_SI, (1.0/3.0));
    double a = a_SI / AU;
    return a;
}

double convertSemiMajorAxisInAUToPeriodInDays(double M1, double M2, double semiMajorAxis){
    /*
     Converts a semi-major axis in AU to orbital period in day
     
     Conjugate function of convertPeriodInDaysToSemiMajorAxisInAU
     
     Parameters
     -----------
     M1 : double
        Primary mass in Msol
     M2 : double
        Secondary mass in Msol
     semiMajorAxis : double
        Semi major axis in AU
     
     Returns
     -------
     period : double
        Orbital period in days

     */
    double period_squared_SI_top = 4.0 * pi * pi * semiMajorAxis * semiMajorAxis * semiMajorAxis * AU * AU * AU;
    double period_squared_SI_bottom = G * ((M1*Msol) + (M2*Msol));
    double period_squared_SI = period_squared_SI_top / period_squared_SI_bottom;
    double period_SI = sqrt(period_squared_SI);
    double period = period_SI / day; // Convert from seconds to days
    return period;
}



// NOTES
