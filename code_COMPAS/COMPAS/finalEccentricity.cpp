//
//  finalEccentricity.cpp

#include "finalEccentricity.h"

double finalEccentricity(double uk, double theta, double phi, double beta, double M, double MPrime, double r, double a){
    /*
     This is Equation (31) in my post-SN orbital characteristics 2 notes - simplifies to Equation (2.8) in Brandt & Podsiadlowski 1995 (http://arxiv.org/abs/astro-ph/9412023) when e = 0, beta = pi/2
     This is also given by Equation A.12 in Hurley et al 2002 (http://arxiv.org/pdf/astro-ph/0201220v1.pdf)
     
     Parameters
     -----------
     uk : double 
        Dimensionless kick velocity vk/vrel
     theta : double
        Angle
     phi : double
        Angle
     beta : double
        Angle the preSN velocity vector makes to the radial position vector
     M : double
        Total system mass before the supernova
     Mprime : double
        Total system mass after the supernova
     r : double
        Current orbital separation in AU (equal to a for eccentricity of 0) at the time of the supernova (sampled randomly for an eccentric orbit)
     a : double
        Semi major axis in AU
     
     Returns
     --------
     ePrime : double
        Orbital eccentricity after a supernova

     */
    
    double MOverMPrime              = M/MPrime;
    double twoOverrMinusOneOvera    = 2.0/r - 1.0/a;
    double quadraticTerm            = 1.0 + 2.0*uk*cos(theta)*cos(phi) + uk*uk;
    double firstSquareBrackets      = uk*uk*sin(theta)*sin(theta) + (uk*cos(theta)*sin(phi)*cos(beta) - sin(beta)*(uk*cos(theta)*cos(phi) + 1))*(uk*cos(theta)*sin(phi)*cos(beta) - sin(beta)*(uk*cos(theta)*cos(phi) + 1.0));
    double secondSquareBrackets     = 2.0/r - MOverMPrime*twoOverrMinusOneOvera*quadraticTerm;
    double oneMinuseSquared         = r*r*MOverMPrime*twoOverrMinusOneOvera*firstSquareBrackets*secondSquareBrackets;
    double eSquared                 = 1.0 - oneMinuseSquared;
    double final_e                  = sqrt(eSquared);
    
    return final_e;
    
}
