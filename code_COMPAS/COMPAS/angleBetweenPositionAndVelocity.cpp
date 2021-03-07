//
//  angleBetweenPositionAndVelocity.cpp


#include "angleBetweenPositionAndVelocity.h"

// Implementation

double angleBetweenPositionAndVelocity(double a, double e, double M, double radius, unsigned long randomSeed){
    // Calculate the angle between the position vector and the velocity vector - call this angle beta
    // See for example Equation A2 in Hurley+ 2002:
    // https://arxiv.org/pdf/astro-ph/0201220.pdf
    // Note for a circular orbit (e = 0) this angle is always pi/2 = 90 degrees.
    // This is Equation (13) in the post-SN orbital characteristics 2 document
    // a = semi-major axis (AU)
    // e = eccentricity
    // M = total mass in Msol
    // radius = current separation (AU)

    // Declare variables
    double beta             = 0;

    // In a circular orbit the angle beta is always 90 degrees = pi/2
    if(e == 0){
        beta = pi/2.0;
    }
    else{

        double sinsquaredbeta   = 0;
        double sinbeta          = 0;

        sinsquaredbeta  = (a*a*(1.0-e*e))/(2.0*radius*a - radius*radius);
        sinbeta         = sqrt(sinsquaredbeta);
        beta            = asin(sinbeta);

        // Check that this never returns NAN due to numerical errors/rounding - this function should now be written in a way that doesn't produce NANs due to numerical rounding
        if(beta != beta){
            std::cerr.precision(std::numeric_limits<double>::digits10 + 2);
            std::cerr << randomSeed << TAB << "ERROR: beta is NAN even in nice units." << std::endl;
            std::cerr << "sin^2(beta) = " << sinsquaredbeta << std::endl;
            std::cerr << "sin(beta) = " << sinbeta << std::endl;
            std::cerr << "a:\t" << a << std::endl;
            std::cerr << "e:\t" << e << std::endl;
            std::cerr << "M:\t" << M << std::endl;
            std::cerr << "radius:\t" << radius << std::endl;
        }

    }

    return beta;

}



///////////////////////////////////////////////////////////////////////////////////
//                              NOTES
///////////////////////////////////////////////////////////////////////////////////

//double angleBetweenPositionAndVelocity(double a, double e, double M, double radius, double v){
//    // Calculate the angle between the position vector and the velocity vector - call this angle beta
//    // Note for a circular orbit (e = 0) this angle is always pi/2 = 90 degrees.
//
//    // All should be in SI?
//    // a = semi-major axis
//    // e = eccentricity
//    // M = total mass
//    // radius = current separation (equal to a for e = 0, not in general otherwise)
//    // v = current orbital velocity
//
//    double sinsquaredbeta   = 0;
//    double sinbeta          = 0;
//    double beta             = 0;
//
//    sinsquaredbeta  = (G*M*a*(1.0 - e*e))/(radius*radius*v*v); // for e = 0, should equal 1
//    sinbeta         = sqrt(sinsquaredbeta);
//    beta            = asin(sinbeta);
//
//    // So basically this is probably all down to using SI units.
//    // Can you do this in nice units?
//
//    double sinsquaredbetalt = 0;
//
//    sinsquaredbetalt = 0;
//
//    // DEBUGGING/ error checks
//    // check that long way and short way give same result
//    // if they do, simplify function
//
//    double answer = 0;
//
//    answer =  asin(sqrt((G*M*a*(1.0 - e*e))/(radius*radius*v*v)));
//
//    if(answer != beta){
//        std::cout << "ERROR: different values :\t" << answer << "\t" << beta << std::endl;
//    }
//
//    // Sometimes beta is nan?? why ??? :( CPLB: Round of error making sinbeta > 1? SPS: I think I have fixed this
//    // can test by using fact that equalitites between nans are always false
//    if(beta != beta){
//        std::cout.precision(std::numeric_limits<double>::digits10 + 2);
//        std::cout << "ERROR: beta is NAN." << std::endl;
//        std::cout << "sin^2(beta) = " << sinsquaredbeta << std::endl;
//        std::cout << "sin(beta) = " << sinbeta << std::endl;
//        std::cout << "a:\t" << a << std::endl;
//        std::cout << "e:\t" << e << std::endl;
//        std::cout << "M:\t" << M << std::endl;
//        std::cout << "radius:\t" << radius << std::endl;
//        std::cout << "vrel:\t" << v << std::endl;
//    }
//
//    return beta;
//
//}
