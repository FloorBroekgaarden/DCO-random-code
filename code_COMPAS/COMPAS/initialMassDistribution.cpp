//
//  initialMassDistribution.cpp


#include "initialMassDistribution.h"

double initialMassDistribution(programOptions options, const gsl_rng *r, AISvariables aisvariables){ 
    
    double thisMass;

    // draw from priors if  DrawingFromAISdistributions = false
    if(!aisvariables.DrawingFromAISdistributions){

        if(boost::iequals(options.initialMassFunction, "Salpeter")){

            // Standard Salpeter IMF

            double SalpeterPower = -2.35;
            
            // Draw mass - pass options plus anything you have here to inverse_sampling_from_power_law.

            thisMass = inverse_sampling_from_power_law(r, SalpeterPower, options.initialMassFunctionMax, options.initialMassFunctionMin);

        }
        else if(boost::iequals(options.initialMassFunction, "POWERLAW")){

            thisMass = inverse_sampling_from_power_law(r, options.initialMassFunctionPower, options.initialMassFunctionMax, options.initialMassFunctionMin); 

        }
        else if(boost::iequals(options.initialMassFunction, "UNIFORM")){

            // Convienience function for POWERLAW with slope of 0

            double u = gsl_rng_uniform(r);              // Draw a random number between 0 and 1

            thisMass = u * (options.initialMassFunctionMax - options.initialMassFunctionMin) + options.initialMassFunctionMin;
        }
        else if(boost::iequals(options.initialMassFunction, "KROUPA")){
            
            // I have moved the constants to the constants.h file

            // C is the normalisation constant for the integral of the IMF over the desired range. 
            double C = 0.0;

            // u will hold a uniform random number U(0,1)
            double u = 0.0;

            // Find our where the user specificed their minimum and maximum masses to generate
            if(options.initialMassFunctionMin <= kroupaBreak1 and options.initialMassFunctionMax <= kroupaBreak1){

                // Draw mass using inverse sampling
                thisMass = inverse_sampling_from_power_law(r, kroupaPower1, options.initialMassFunctionMax, options.initialMassFunctionMin);

            }
            else if(options.initialMassFunctionMin > kroupaBreak1 and options.initialMassFunctionMin <= kroupaBreak2 and options.initialMassFunctionMax > kroupaBreak1 and options.initialMassFunctionMax <= kroupaBreak2){

                // Draw mass using inverse sampling
                thisMass = inverse_sampling_from_power_law(r, kroupaPower2, options.initialMassFunctionMax, options.initialMassFunctionMin);

            }
            else if(options.initialMassFunctionMin > kroupaBreak2 and options.initialMassFunctionMax > kroupaBreak2){

                // Draw mass using inverse sampling
                thisMass = inverse_sampling_from_power_law(r, kroupaPower3, options.initialMassFunctionMax, options.initialMassFunctionMin);

            }
            else if(options.initialMassFunctionMin <= kroupaBreak1 and options.initialMassFunctionMax > kroupaBreak1 and options.initialMassFunctionMax <= kroupaBreak2){

                //std::cout << "12" << std::endl;

                u = gsl_rng_uniform(r);              // Draw a random number between 0 and 1

                double term1 = oneOverKroupaPower1Plus1 * (pow(kroupaBreak1, kroupaPowerPlus1_1) - pow(options.initialMassFunctionMin, kroupaPowerPlus1_1));
                double term2 = oneOverKroupaPower2Plus1 * pow(kroupaBreak1, (kroupaPower1 - kroupaPower2)) * (pow(options.initialMassFunctionMax, kroupaPowerPlus1_2) - pow(kroupaBreak1, kroupaPowerPlus1_2));

                double oneOverC1 = term1 + term2;

                double C1 = 1.0 / oneOverC1;
                double C2 = C1 * pow(kroupaBreak1, (kroupaPower1 - kroupaPower2));

                // Calculate value of CDF at breaks
                double CDFatKroupaBreak1 = cdfKroupa(kroupaBreak1, options);

                double A = oneOverKroupaPower1Plus1 * C1 * (pow(kroupaBreak1, kroupaPowerPlus1_1) - pow(options.initialMassFunctionMin, kroupaPowerPlus1_1));

                if(u < CDFatKroupaBreak1){

                    thisMass = pow(u * (kroupaPowerPlus1_1/C1) + pow(options.initialMassFunctionMin, kroupaPowerPlus1_1), oneOverKroupaPower1Plus1);
                }
                else{

                    thisMass = pow((u - A) * (kroupaPowerPlus1_2/C2) + pow(kroupaBreak1, kroupaPowerPlus1_2), oneOverKroupaPower2Plus1);
                }

            }
            else if(options.initialMassFunctionMin <= kroupaBreak1 and options.initialMassFunctionMax > kroupaBreak2){

                //std::cout << "123" << std::endl;

                u = gsl_rng_uniform(r);              // Draw a random number between 0 and 1

                double term1 = oneOverKroupaPower1Plus1 * (pow(kroupaBreak1, kroupaPowerPlus1_1) - pow(options.initialMassFunctionMin, kroupaPowerPlus1_1));
                double term2 = oneOverKroupaPower2Plus1 * pow(kroupaBreak1, (kroupaPower1 - kroupaPower2)) * (pow(kroupaBreak2, kroupaPowerPlus1_2) - pow(kroupaBreak1, kroupaPowerPlus1_2));
                double term3 = oneOverKroupaPower3Plus1 * pow(kroupaBreak1, (kroupaPower1-kroupaPower2)) * pow(kroupaBreak2, kroupaPower2-kroupaPower3) * (pow(options.initialMassFunctionMax, kroupaPowerPlus1_3) - pow(kroupaBreak2, kroupaPowerPlus1_3));

                double oneOverC1 = term1 + term2 + term3;

                double C1 = 1.0 / oneOverC1;
                double C2 = C1 * pow(kroupaBreak1, (kroupaPower1 - kroupaPower2));
                double C3 = C2 * pow(kroupaBreak2, (kroupaPower2 - kroupaPower3));

                // Calculate value of CDF at breaks
                double CDFatKroupaBreak1 = cdfKroupa(kroupaBreak1, options);
                double CDFatKroupaBreak2 = cdfKroupa(kroupaBreak2, options);

                double A = oneOverKroupaPower1Plus1 * C1 * (pow(kroupaBreak1, kroupaPowerPlus1_1) - pow(options.initialMassFunctionMin, kroupaPowerPlus1_1));

                double B = oneOverKroupaPower2Plus1 * C2 * (pow(kroupaBreak2, kroupaPowerPlus1_2) - pow(kroupaBreak1, kroupaPowerPlus1_2));

                //std::cout << "u = " << u << std::endl;
                //std::cout << "CDF at kroupaBreak1 (" << kroupaBreak1 << ") = " << CDFatKroupaBreak1 << std::endl;
                //std::cout << "CDF at kroupaBreak2 (" << kroupaBreak2 << ") = " << CDFatKroupaBreak2 << std::endl;

                if(u < CDFatKroupaBreak1){
                    thisMass = pow(u * (kroupaPowerPlus1_1/C1) + pow(options.initialMassFunctionMin, kroupaPowerPlus1_1), oneOverKroupaPower1Plus1);
                }
                else if(u >= CDFatKroupaBreak1 and u < CDFatKroupaBreak2){
                    thisMass = pow((u - A) * (kroupaPowerPlus1_2/C2) + pow(kroupaBreak1, kroupaPowerPlus1_2), oneOverKroupaPower2Plus1);
                }
                else if(u >= CDFatKroupaBreak2){
                    thisMass = pow((u - A - B) * (kroupaPowerPlus1_3/C3) + pow(kroupaBreak2, kroupaPowerPlus1_3), oneOverKroupaPower3Plus1);
                }
                else{
                    std::cout << "Error in inverse sampling from Kroupa IMF" << std::endl;
                }

            }
            else if(options.initialMassFunctionMin > kroupaBreak1 and options.initialMassFunctionMin <= kroupaBreak2 and options.initialMassFunctionMax > kroupaBreak2){

                //std::cout << "23" << std::endl;

                u = gsl_rng_uniform(r);              // Draw a random number between 0 and 1

                double term1 = oneOverKroupaPower2Plus1 * (pow(kroupaBreak2, kroupaPowerPlus1_2) - pow(options.initialMassFunctionMin, kroupaPowerPlus1_2));
                double term2 = oneOverKroupaPower3Plus1 * pow(kroupaBreak2, (kroupaPower2 - kroupaPower3)) * (pow(options.initialMassFunctionMax, kroupaPowerPlus1_3) - pow(kroupaBreak2, kroupaPowerPlus1_3));

                double oneOverC2 = term1 + term2;

                double C2 = 1.0 / oneOverC2;
                double C3 = C2 * pow(kroupaBreak2, (kroupaPower2 - kroupaPower3));

                // Calculate value of CDF at break
                double CDFatKroupaBreak2 = cdfKroupa(kroupaBreak2, options);

                double B = oneOverKroupaPower2Plus1 * C2 * (pow(kroupaBreak2, kroupaPowerPlus1_2) - pow(options.initialMassFunctionMin, kroupaPowerPlus1_2));

                //std::cout << "u, cdfatkroupabreak2 = " << u << " " << CDFatKroupaBreak2 << std::endl;

                if(u < CDFatKroupaBreak2){

                    //u = 0;

                    thisMass = pow(u * (kroupaPowerPlus1_2/C2) + pow(options.initialMassFunctionMin, kroupaPowerPlus1_2), oneOverKroupaPower2Plus1);

                    //std::cout << "mass should be min = " << thisMass << std::endl;

                }
                else{

                    thisMass = pow((u - B) * (kroupaPowerPlus1_3/C3) + pow(kroupaBreak2, kroupaPowerPlus1_3), oneOverKroupaPower3Plus1);

                }

            }
            else{

                std::cerr << "initialMassDistribution.cpp - shouldn't get here" << std::endl;

                thisMass = 0.0;

            }

        }
        else{
            
            // USER HAS SPECIFIED AN INVALID IMF - USE DEFAULT
            std::cout << "Invalid choice of IMF -- using the default instead." << std::endl;
            
            // Reset to defaults.
            double kroupaPower = -2.35;             //
            double kroupaMin   = 0.5;               //
            double kroupaMax   = 100;               //
            
            // Generate mass using this
            thisMass = inverse_sampling_from_power_law(r, kroupaPower, kroupaMax, kroupaMin);
            
        }
    }
    
    // draw from AIS distributions if  DrawingFromAISdistributions = true
    else{
        // Mass ratio distribution from Adaptive Importance Sampling v1 from Broekgaarden et al. (in prep 2018)
        // Function Returns a random Primary Mass drawn from one of the random gaussians defined bu vectors mu_M1 & cov_M1

        // Now draw a random Mass from the random Gaussian chosen with RandomGaussianDraw
        double thisMassMean =  aisvariables.mu_M1[aisvariables.RandomGaussianDraw];    // The mean of the RandomGaussianDraw-th Gaussian 
        double thisMassStd  =  aisvariables.cov_M1[aisvariables.RandomGaussianDraw];   // The cov of the RandomGaussianDraw-th Gaussin
        thisMass = gsl_ran_gaussian (r, thisMassStd) + thisMassMean; // draw random nr from Gaussian 
    }
    
    return thisMass;
}

double cdfKroupa(double x, programOptions options){
    /*

    Calculate the value of the CDF of the Kroupa (2001) IMF at x

    Parameters
    -----------
    x : double
        Mass value (in solar masses) to calculate the CDF at
    options : programOptions
        User specified program options

    Returns
    --------
    CDF : double
        value of the CDF
    */
    double term1 = 0.0;
    double term2 = 0.0;
    double term3 = 0.0;

    double oneOverC1 = 0.0;

    double C1 = 0.0;
    double C2 = 0.0;
    double C3 = 0.0;

    double CDF = 0;

    if(options.initialMassFunctionMin <= kroupaBreak1 and options.initialMassFunctionMax > kroupaBreak1 and options.initialMassFunctionMax <= kroupaBreak2){

        term1 = oneOverKroupaPower1Plus1 * (pow(kroupaBreak1, kroupaPowerPlus1_1) - pow(options.initialMassFunctionMin, kroupaPowerPlus1_1));
        term2 = oneOverKroupaPower2Plus1 * pow(kroupaBreak1, (kroupaPower1 - kroupaPower2)) * (pow(options.initialMassFunctionMax, kroupaPowerPlus1_2) - pow(kroupaBreak1, kroupaPowerPlus1_2));

        oneOverC1 = term1 + term2;

        C1 = 1.0 / oneOverC1;
        C2 = C1 * pow(kroupaBreak1, (kroupaPower1 - kroupaPower2));

        if(x >= options.initialMassFunctionMin and x < kroupaBreak1){

            CDF = oneOverKroupaPower1Plus1 * C1 * (pow(x, kroupaPowerPlus1_1) - pow(options.initialMassFunctionMin, kroupaPowerPlus1_1));

        }
        else if(x >= kroupaBreak1 and x < kroupaBreak2){

            CDF = oneOverKroupaPower1Plus1 * C1 * (pow(kroupaBreak1, kroupaPowerPlus1_1) - pow(options.initialMassFunctionMin, kroupaPowerPlus1_1)) + 

                  oneOverKroupaPower2Plus1 * C2 * (pow(x, kroupaPowerPlus1_2) - pow(kroupaBreak1, kroupaPowerPlus1_2));

        }
        else{
            std::cout << "Error in kroupa CDF" << std::endl;
        }

    }
    else if(options.initialMassFunctionMin <= kroupaBreak1 and options.initialMassFunctionMax > kroupaBreak2){

        term1 = oneOverKroupaPower1Plus1 * (pow(kroupaBreak1, kroupaPowerPlus1_1) - pow(options.initialMassFunctionMin, kroupaPowerPlus1_1));
        term2 = oneOverKroupaPower2Plus1 * pow(kroupaBreak1, (kroupaPower1 - kroupaPower2)) * (pow(kroupaBreak2, kroupaPowerPlus1_2) - pow(kroupaBreak1, kroupaPowerPlus1_2));
        term3 = oneOverKroupaPower3Plus1 * pow(kroupaBreak1, (kroupaPower1-kroupaPower2)) * pow(kroupaBreak2, kroupaPower2-kroupaPower3) * (pow(options.initialMassFunctionMax, kroupaPowerPlus1_3) - pow(kroupaBreak2, kroupaPowerPlus1_3));

        oneOverC1 = term1 + term2 + term3;

        C1 = 1.0 / oneOverC1;
        C2 = C1 * pow(kroupaBreak1, (kroupaPower1 - kroupaPower2));
        C3 = C2 * pow(kroupaBreak2, (kroupaPower2 - kroupaPower3));

        if(x >= options.initialMassFunctionMin and x < kroupaBreak1){

            CDF = oneOverKroupaPower1Plus1 * C1 * (pow(x, kroupaPowerPlus1_1) - pow(options.initialMassFunctionMin, kroupaPowerPlus1_1));

        }
        else if(x >= kroupaBreak1 and x < kroupaBreak2){

            CDF = oneOverKroupaPower1Plus1 * C1 * (pow(kroupaBreak1, kroupaPowerPlus1_1) - pow(options.initialMassFunctionMin, kroupaPowerPlus1_1)) + 

                  oneOverKroupaPower2Plus1 * C2 * (pow(x, kroupaPowerPlus1_2) - pow(kroupaBreak1, kroupaPowerPlus1_2));

        }
        else if(x >= kroupaBreak2 and x < options.initialMassFunctionMax){

            CDF =   oneOverKroupaPower1Plus1 * C1 * (pow(kroupaBreak1, kroupaPowerPlus1_1) - pow(options.initialMassFunctionMin, kroupaPowerPlus1_1)) + 

                    oneOverKroupaPower2Plus1 * C2 * (pow(kroupaBreak2, kroupaPowerPlus1_2) - pow(kroupaBreak1, kroupaPowerPlus1_2)) + 

                    oneOverKroupaPower3Plus1 * C3 * (pow(x, kroupaPowerPlus1_3) - pow(kroupaBreak2, kroupaPowerPlus1_3));

        }
        else{
            std::cout << "Error in kroupa CDF" << std::endl;
        }

    }
    else if(options.initialMassFunctionMin > kroupaBreak1 and options.initialMassFunctionMin <= kroupaBreak2 and options.initialMassFunctionMax > kroupaBreak2){

        //std::cout << "cdf 23" << std::endl;

        term1 = oneOverKroupaPower2Plus1 * (pow(kroupaBreak2, kroupaPowerPlus1_2) - pow(options.initialMassFunctionMin, kroupaPowerPlus1_2));
        term2 = oneOverKroupaPower3Plus1 * pow(kroupaBreak2, (kroupaPower2 - kroupaPower3)) * (pow(options.initialMassFunctionMax, kroupaPowerPlus1_3) - pow(kroupaBreak2, kroupaPowerPlus1_3));

        double oneOverC2 = term1 + term2;

        C2 = 1.0 / oneOverC2;
        C3 = C2 * pow(kroupaBreak2, (kroupaPower2 - kroupaPower3));

        if(x >= options.initialMassFunctionMin and x < kroupaBreak2){

            CDF = oneOverKroupaPower2Plus1 * C2 * (pow(x, kroupaPowerPlus1_2) - pow(options.initialMassFunctionMin, kroupaPowerPlus1_2));

        }
        else if(x >= kroupaBreak2 and x < options.initialMassFunctionMax){

            CDF = oneOverKroupaPower2Plus1 * C2 * (pow(kroupaBreak2, kroupaPowerPlus1_2) - pow(options.initialMassFunctionMin, kroupaPowerPlus1_2)) + 

                  oneOverKroupaPower3Plus1 * C3 * (pow(x, kroupaPowerPlus1_3) - pow(kroupaBreak2, kroupaPowerPlus1_3))  ;

        }
        else{
            std::cout << "Error in kroupa CDF" << std::endl;
        }

    }

    return CDF;
}
