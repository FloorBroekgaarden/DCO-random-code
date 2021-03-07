//
//  rotationalVelocityDistribution.cpp


#include "rotationalVelocityDistribution.h"

// Implementation goes here

double hurley2000RotationalVelocityDistribution(double MZAMS){
    /*

    Fit from Hurley et al 2000 based on Lang 1992

    Parameters
    ------------
    MZAMS : double
        Zero age main sequence mass in Msol

    Returns
    --------
    ve : double
        Equatorial rotational velocity in km s^-1

    */
    return (330.0 * pow(MZAMS, 3.3))/(15.0 + pow(MZAMS, 3.45));
}

double rotationalVelocityDistribution(double MZAMS, programOptions options, const gsl_rng *r){
    /*
    Return the inital rotational velocity (in km s^{-1} ) of a star with ZAMS mass MZAMS

    Parameters
    ----------
    MZAMS : double
        Zero age main sequence mass in solar masses
    options : programOptions
        Contains user defined program options
    r : gsl_rng
        Random number generator
    
    Returns
    --------
    ve : double
        Initial equatorial rotational velocity in km s^{-1}

    */
    
    if(boost::iequals(options.rotationalVelocityString, "ZERO")){
        return 0;
    }
    else if(boost::iequals(options.rotationalVelocityString, "HURLEY")){
        // Fit from Hurley et al 2000 based on Lang 1992
        return hurley2000RotationalVelocityDistribution(MZAMS);
    }
    else if(boost::iequals(options.rotationalVelocityString, "VLTFLAMES")){
        // Rotational velocity based on VLT-FLAMES survey.
        // For O-stars use results of Ramirez-Agudelo et al (2013) https://arxiv.org/abs/1309.2929 (single stars) 
        // and Ramirez-Agudelo et al. (2015) https://arxiv.org/abs/1507.02286 (spectroscopic binaries)
        // For B-stars use results of Dufton et al. (2013) https://arxiv.org/abs/1212.2424
        // For lower mass stars, I don't know what updated results there are so default back to 
        // Hurley et al 2000 distribution for now
        
        bool debugging = false;

        double u = gsl_rng_uniform(r);              // Draw a random number between 0 and 1

        double ve = 0;                       // Equatorial rotational velocity in km s^-1

        if(MZAMS >= 16.0){

            ve = ramirezAgudelo2013OStarRotationalVelocityAnalyticCDFInverseSampling(u, 0.0, 800.0);

            if(debugging){
                std::cout << "O type ve = " << ve << std::endl;
            }

        }
        else if (MZAMS >= 2.0 and MZAMS < 16.0){

            ve = inverse_sampling_from_tabulated_cdf(u, dufton2013BStarRotationalVelocityCDFTable_velocity, dufton2013BStarRotationalVelocityCDFTable_cdf, nrows_dufton2013BStarRotationalVelocityCDFTable, dufton2013BStarRotationalVelocityCDFTable_maximum_velocity, dufton2013BStarRotationalVelocityCDFTable_minimum_velocity, dufton2013BStarRotationalVelocityCDFTable_maximum_cdf);

            if(debugging){
                std::cout << "B type ve = " << ve << std::endl;
            }

        }
        else{

            // Don't know what better to use for low mass stars so for now default back to Hurley et al 2000 distribution

            ve = hurley2000RotationalVelocityDistribution(MZAMS);

        }

        return ve;
    }
    else{
        // User has specified an invalid distribution - use default
        std::cerr << "Invalid rotational velocity distribution - using 0 by default" << std::endl;
        return 0;
    }
}

double rotationalAngularVelocityDistribution(double MZAMS, double RZAMS, programOptions options, const gsl_rng *r){
    /*
    Return the inital angular velocity (in yr^{-1}) of a star with ZAMS mass and radius MZAMS and RZAMS respectively
    
    Parameters
    ------------
    MZAMS : double
        Zero age main sequence mass in Msol
    RZAMS : double
        Zero age main sequence radius in Rsol
    options : programOptions
        User defined program options
    r : gsl_rng
        Random number generator

    Returns
    ---------
    omega : double
        initial angular velocity (in yr^{-1}) of star
     */

    double vrot = rotationalVelocityDistribution(MZAMS, options, r);

    return 45.35/(vrot / RZAMS);
}