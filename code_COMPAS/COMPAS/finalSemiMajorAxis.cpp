//
//  finalSemiMajorAxis.cpp


#include "finalSemiMajorAxis.h"

double finalSemiMajorAxis(double uk, double M, double MPrime, double r, double a, double theta, double phi){
    /*
     This function calculates the post-supernova semi major axis
     
     This is Equation (22) in the post-SN orbital characteristics 2 document
     
     Parameters
     ----------
     uk : double
        Dimensionless kick velocity in units of the preSN orbital velocity
     M : double
        Total system mass before the supernova
     MPrime : double
        Total system mass after the supernova
     r : double
        r is the instantaneous separation at the moment of the SN (accounts for eccentricity) in AU
     a : double
        Semi major axis of the orbit before the supernova
     theta : double
        theta is the kick direction angle out of the plane
     phi : double
        phi is the kick direction angle in the plane
     
     Returns
     -------
     aPrime : double
        Semi major axis of the orbit after the supernova

     */
    
    double MOverMPrime   = M/MPrime;
    double firstBrackets = (2.0/r - 1.0/a);
    double quadraticTerm = 1.0 + 2.0*uk*cos(theta)*cos(phi) + uk*uk;
    double bottom        = 2.0/r - MOverMPrime*firstBrackets*quadraticTerm;
    double final_a       = 1.0/bottom;
    
    return final_a;
    
}
