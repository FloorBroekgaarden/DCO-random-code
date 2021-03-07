//
//  rocheLobeRadius.cpp

#include "rocheLobeRadius.hpp"

double rocheLobeRadius(double Mass_i, double Mass_j){
    /*
     Calculate the roche lobe radius using the chosen approximation
     
     By default uses the fit given in Equation 2 of Eggleton 1983 (Could in principle pass it options and use other fits like Sepinsky+)
     
     Parameters
     -----------
     Mass_i : double
        Primary mass
     Mass_j : double
        Secondary mass
     
     Returns
     --------
     R_L / a : double
     Radius of Roche lobe in units of the semi-major axis a

     */
    double q = Mass_i/Mass_j;
    double top = 0.49;
    double bottom = 0.6 + pow(q, -2.0/3.0) * log(1.0 + pow(q, 1.0/3.0));
    return top/bottom;
}

double rocheLobeResponse(double ma, double md, double fa){
    // Calculation of the adiabatic index acording to the derivation in Woods, Ivanova and van der Sluys (2012), eq. 9-13
    double  q = md/ma;
    double  q1 = pow(q, 1/3);
    double  num1 = q1*((1.2*q1)+(1/(1+q1)));
    double  den1 = 3*((0.6*pow(q,2/3))+log(1+q1));
    double  RLMassRatio = (2/3)-(num1/den1);        // Eq. (12)
    double  RMassDonor = 1+(fa*md/ma);              // Eq. (13)
    double  num2 = (2*md*md)-(2*ma*ma)-(md*ma*(1-fa));
    double  den2 = ma*(md+ma);
    double  aMassRatio = num2/den2;                 // Eq. (11)
    
    return aMassRatio+(RLMassRatio*RMassDonor);     // Eq. (9)
}

