//
//  solveKeplersEquation.cpp

#include "solveKeplersEquation.h"

void solveKeplersEquation(double meanAnomaly, double eccentricity, double & eccentricAnomaly, double & trueAnomaly){
    /*
    Solve Kepler's Equation using root finding techniques. Here we use Newton-Raphson.

    Needs to be passed the mean anomaly and eccentricity and pass by reference eccentric anomaly and true anomaly which are initially empty

    For a definition of all the anomalies see here:

    https://en.wikipedia.org/wiki/Mean_anomaly
    https://en.wikipedia.org/wiki/True_anomaly
    https://en.wikipedia.org/wiki/Eccentric_anomaly
    
    Parameters
    ----------
    meanAnomaly : double
        Mean anomaly
    eccentricity : double
        Eccentricity
    eccentricAnomaly : double
        Eccentric anomaly
    trueAnomaly : double
        True anomaly

    Returns
    --------
    */
    
    double Etmp = meanAnomaly;        // Temporary variable for the eccentric anomaly E - inital guess at E is M - correct for eccentricity = 0
    
    // Repeat the approximation until Etmp is within the specified error of the true value
    while(fabs(keplersEquation(eccentricity, Etmp, meanAnomaly)) >= epsilonNR){
        Etmp = Etmp - keplersEquation(eccentricity, Etmp, meanAnomaly)/keplersEquationPrimed(eccentricity, Etmp, meanAnomaly);
    }
    
    // Set final value of eccentric anomaly E
    eccentricAnomaly = Etmp;
    
    // Convert Eccentric anomaly into true anomaly
    eccentricToTrue(eccentricity, trueAnomaly, eccentricAnomaly);

}

double keplersEquation(double eccentricity, double eccentricAnomaly, double meanAnomaly){
    // Return f(E) = 0
    // This is Equation (92) in my "A simple toy model" document
    
    double e = eccentricity;
    double E = eccentricAnomaly;
    double M = meanAnomaly;
    return E - e*sin(E) - M;
}

double keplersEquationPrimed(double eccentricity, double eccentricAnomaly, double meanAnomaly){
    // Return the derivative of f(E), f'(E)
    // This is Equation (94) in my "A simple toy model" document
    
    double e = eccentricity;
    double E = eccentricAnomaly;
    return 1.0 - e*cos(E);
}

void eccentricToTrue(double eccentricity, double & trueAnomaly, double eccentricAnomaly){
    // Convert Eccentric anomaly into true anomaly.
    // This is Equation (96) in my "A simple toy model" document
    
    double e  = eccentricity;
    double E  = eccentricAnomaly;
    double nu = trueAnomaly;
    
    if(E >= M_PI && E <= 2.0*M_PI){
        nu = 2.0*atan((sqrt((1.0 + e)/(1.0 - e))) * tan(0.5*E)) + 2.0*M_PI;
    }
    else if(E >= 0.0 && E < M_PI){
        nu = 2.0*atan((sqrt((1.0 + e)/(1.0 - e))) * tan(0.5*E));
    }
    else{
        nu = 2.0*atan((sqrt((1.0 + e)/(1.0 - e))) * tan(0.5*E));
        
        std::cerr << "Something went wrong converting E to nu. E was " << E << std::endl;
    }
    
    trueAnomaly = nu;
    
}
