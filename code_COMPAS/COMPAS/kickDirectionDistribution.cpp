//
//  kickDirectionDistribution.cpp

#include "kickDirectionDistribution.h"

void kickDirectionDistribution(const programOptions &options, const gsl_rng *r, double &theta, double &phi, double &kickDirectionPower){
    /*
    
    Draws the angular components of the supernova kick theta and phi according to user specified options.

    Parameters
    -----------
    options : programOptions
        User specified options
    r : gsl_rng
        Used for the random number generation
    theta : double
        Polar angle for kick (pointer)
    phi : double
        Azimuthal angle for kick (pointer)
    kickDirectionPower : double
        Power law for the POWER angle distribution

    Returns
    --------
    updates values of theta and phi
    */

    bool debugging = false;
    //debugging = true;

    // Declare some variables for kick angles -- don't really need to be declared every time around a loop
    double n_degrees    = 1.0;                        // Can be set by user - size of cone/wedge
    double delta        = n_degrees*degree;         // Small angle () in radians - could be set by user in options
    double u            = gsl_rng_uniform(r);       // Generate a random number -- make sure you only use u once per cycle - should not really bother defining this here

    if(boost::iequals(options.kickDirection, "ISOTROPIC")){
        // Draw theta and phi isotropically
        theta = acos(1.0 - 2.0*gsl_rng_uniform(r)) - (pi/2.0);          // Theta, angle out of the plane
        phi   = gsl_rng_uniform(r)*2.0*pi;                              // Phi, angle in the plane
    }
    else if(boost::iequals(options.kickDirection, "POWERLAW")){
        // Draw phi uniform in [0,2pi], theta according to a powerlaw (power law power = 0 = isotropic, +infinity = kick along pole, -infinity = kick in plane)

        // double magnitude_of_cos_theta = inverse_sampling_from_power_law(r, fabs(options.kickDirectionPower), 1.0, 1E-6);    // Choose magnitude of power law distribution -- if using a negative power law that blows up at 0, need a lower cutoff (currently set at 1E-6), check it doesn't affect things too much
        double magnitude_of_cos_theta = inverse_sampling_from_power_law(r, fabs(kickDirectionPower), 1.0, 1E-6);

        // if(options.kickDirectionPower < 0.0){
        if(kickDirectionPower < 0.0){

            magnitude_of_cos_theta = 1.0 - magnitude_of_cos_theta;
        }

        double actual_cos_theta = magnitude_of_cos_theta;

        // By taking the magnitude of cos theta we lost the information about whether it was up or down, so we put that back in randomly here
        if(gsl_rng_uniform(r) < 0.5){
            actual_cos_theta = -magnitude_of_cos_theta;
        }

        // Check costheta > -1.0 and < 1.0
        if(actual_cos_theta < -1.0){
            actual_cos_theta = -1.0;
        }

        if(actual_cos_theta > 1.0){
            actual_cos_theta = 1.0;
        }

        theta = acos(actual_cos_theta);                                          // Set the kick angle out of the plane theta

        phi   = gsl_rng_uniform(r)*2.0*pi;                    // Allow to randomly take an angle 0 - 2pi in the plane
    }
    else if(boost::iequals(options.kickDirection, "INPLANE")){
        // Force the kick to be in the plane theta = 0
        theta = 0.0;                                          // Force the kick to be in the plane - theta = 0
        phi   = gsl_rng_uniform(r)*2.0*pi;                    // Allow to still take an angle 0 - 2pi
    }
    else if(boost::iequals(options.kickDirection, "PERPENDICULAR")){
        // Force kick to be along spin axis

        if( u >= 0.0 && u < 0.5){
            // UP
            theta = 0.5*pi;                                       // Theta = pi/2 -- up
            phi   = gsl_rng_uniform(r)*2.0*pi;                    // Allow to still take an angle 0 - 2pi
        }
        else if(u >= 0.5 && u < 1.0){
            // DOWN
            theta = -0.5*pi;                                      // Theta = -pi/2 -- down
            phi   = gsl_rng_uniform(r)*2.0*pi;                    // Allow to still take an angle 0 - 2pi
        }
        else{
            // ERROR checking in case I did something wrong
            std::cerr << "ERROR: KICK DIRECTION : u = " << u << std::endl;
            theta = 0.0;
            phi   = 0.0;
        }

    }
    else if(boost::iequals(options.kickDirection, "POLES")){
        // Direct the kick in a small cone around the poles
        phi   = gsl_rng_uniform(r)*2.0*pi;                    // Allow to still take an angle 0 - 2pi

        // CAREFUL - NOT SO STRAIGHTFORWARD DUE TO THE CUT-OFF AT -pi/2 / pi/2 - USE ABS(delta) to avoid problems
        if( u >= 0.0 && u < 0.5){
            // UP

            // Draw number (changes value of u) from a gaussian of width delta
            u = gsl_ran_gaussian(r, delta); // CPLB: Change name!

            // Take abs value
            u = fabs(u);

            theta = 0.5*pi - u;                         // slightly less than or equal to pi/2
        }
        else if(u >= 0.5 && u < 1.0){
            // DOWN

            // Draw number from a gaussian of width delta
            u = gsl_ran_gaussian(r, delta); // CPLB: Change name!

            // Take abs value
            u = fabs(u);

            theta = -0.5*pi + u;                        // slightly more than or equal to -pi/2
        }
        else{
            // ERROR checking in case I did something wrong
            std::cerr << "ERROR: KICK DIRECTION : u = " << u << std::endl;
        }
    }
    else if(boost::iequals(options.kickDirection, "WEDGE")){
        // Direct kick into a wedge around the horizon (theta = 0)
        theta = gsl_ran_gaussian(r, delta);               // Gaussian around 0 with a deviation delta
        phi   = gsl_rng_uniform(r)*2.0*pi;                  // Phi, angle in the plane - uniform between 0 and 2pi like normal
    }
    else{
        std::cerr << "Invalid kick direction choice." << std::endl;

        std::cout << "kick direction: " << options.kickDirection << std::endl;

        if(boost::iequals(options.kickDirection, "ISOTROPIC")){
            std::cout << "TRUE" << std::endl;
        }
        else{
            std::cout << "FALSE" << std::endl;
        }

        // Draw theta and phi isotropically
        theta = acos(1.0 - 2.0*u) - (0.5*pi);                   // Theta, angle out of the plane
        phi   = gsl_rng_uniform(r)*2.0*pi;                    // Phi, angle in the plane
    }
}