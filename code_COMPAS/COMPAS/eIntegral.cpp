//
//  eIntegral.cpp


#include "eIntegral.h"

void eIntegral(gsl_interp_accel *&my_accel_ptr, gsl_spline *&my_spline_ptr){
    // The *& means we are passing the pointers by reference so that the object they point to can be modified.
    // CPLB: Why not just pass the things? (Am I being stupid?)
    
    // Declare things
    double eMin         = 0;                                        // Minimum eccentricity -- circular
    double eMax         = 1;                                        // Maximum eccentricity -- straight line
    int nPoints         = 999;                                      // Number of points to calculate the integral at
    int nPointsActual   = nPoints + 1;                              // Actual number of points including e = 1
    int nPointsIntegral = 1E5;                                      // Number of points to use in Reimann sum in integral
    double de           = (eMax - eMin) / nPoints;                  // Eccentricity step
    // May want to have nPointsActual and/or nPointsIntegral set by user.
    // Although since you are creating arrays which depend on this, they would have to be dynamically allocated which I'm not sure I remember how to do. -- see DynamicArrays file for revision on that
    
    // Declare arrays
    double e[nPointsActual];                                        // Declare e array
    double Integral[nPointsActual];                                 // Declare integral array
    double tOverTcArray[nPointsActual];                             // Declare array to store t/tc

    // Initialise arrays
    for(int i = 0; i < nPointsActual; i++){
        e[i]                    = 0;                                // Initialise e array
        Integral[i]             = 0;                                // Initialise integral array
        tOverTcArray[i]         = 0;                                // Intialise array for t/tc using numerics
    }
    
    // Inform user
    std::cout << "Calculating e integral numerically at " << nPointsActual << " values of e." << std::endl;

    // Calculate the integral on a relatively fine grid using numerical integration
    for(int i = 0; i < nPointsActual; i++){
        
        // Calculate e
        e[i] = i*de;
        
        // Compute integral at that e
        Integral[i] = e_integral(0, e[i], nPointsIntegral);
        
        // Calculate the ratio of t/tc
        tOverTcArray[i] = tovertc(e[i], Integral[i]);
        
    }
    
    // Allocate pointers:
    // Check if the pointer is Null (initially it should be):
    if(my_accel_ptr == NULL){
        my_accel_ptr  = gsl_interp_accel_alloc();                                       // Allocate 'interpolation accelerator'
    }
    if(my_spline_ptr == NULL){
        my_spline_ptr = gsl_spline_alloc(gsl_interp_cspline, nPointsActual);            // Allocate 'spline object'
    }
    
    // Now check that they were correctly allocated using:
    if(my_accel_ptr == NULL){
        std::cerr << "ERROR: ACCEL ptr is NULL" << std::endl;                           // Tell user the pointer is null
    }
    if(my_spline_ptr == NULL){
        std::cerr << "ERROR: SPLINE ptr is NULL" << std::endl;                          // Tell user the pointer is null
        
    }
    
    // Initialise the spline for the results of t/tc:                                   // Should you do this here or in main?
    gsl_spline_init(my_spline_ptr, e, tOverTcArray, nPointsActual);                     // It requires e, tovertcarray so i hope i
                                                                                        // can do it here or i have wasted all this time
}

// Implementation of other functions
// CPLB: Formulae checked!
double e_integrand(double e){
    // Calculate the integrand of the e integral at the required point
    // This is Equation (5.14) in Peters 1964 (http://journals.aps.org/pr/pdf/10.1103/PhysRev.136.B1224)
    
    // ERROR PREVIOUSLY - NOW CHANGED
    // CHECK THIS -- SHOULD THIS LAST TERM BE -1.5 NOT 1.5 ??
    // TRY CHECKING AGAINST THE PETERS ODEINT INTEGRATOR
    // double x = pow(e, (29.0/19.0))*pow((1.0 + (121.0/304.0)*e*e), (1181.0/2299.0))*pow((1.0 - e*e), 1.5);
    double x = pow(e, (29.0/19.0))*pow((1.0 + (121.0/304.0)*e*e), (1181.0/2299.0))*pow((1.0 - e*e), -1.5);
    return x;
}

// CPLB: We should be able to do better than a simple Euler method!
// SPS: But do we need to? If so then I can fix this.
double e_integral(double lowerLimit=0.01, double upperLimit=0.99, int nPointsIntegral=1E5){
    // Integral in the coalescence expression in terms of eccentricity e. Must be integrated numerically, performed here using the Riemann Sum method.
    double sum      = 0.0;
    double e        = 0.0;
    double deltae   = (upperLimit - lowerLimit)/nPointsIntegral;
    for(int i = 0; i<nPointsIntegral; i++){
        e    = i*deltae;
        sum += e_integrand(e)*deltae;
    }
    return sum;
}

double tovertc(double e0, double integral){
    // Ratio of time to coalescence over the time for a circular orbit
    // This is derived from Equation (5.10), (5.11) and (5.14) in Peters 1964 (http://journals.aps.org/pr/pdf/10.1103/PhysRev.136.B1224)
    // This is also Equation (113) in my "A simple toy model" notes
    
    if(e0 == 0.0){
        return 1.0;
    }
    else if(e0 > 0.0 and e0 <= 1.0){ // CPLB: Only for e0 < 1    // SPS: fixed
        double squareBrackets = (1.0 + (121.0/304.0)*e0*e0);
        return (48.0/19.0)*pow((1.0 - e0*e0), 4.0)*pow(squareBrackets, (-3480.0/2299.0))*pow(e0, (-48.0/19.0))*integral;
    }
    else{
        std::cerr << "Error : e outside range [0, 1] e = " << e0 << std::endl;
        return 0;
    }
    
}

double polynomialFitToToverTc(double e0){
    /*
     Calculate a polynomial fit to the ratio of the inspiral time for a general eccentric orbit to that of a circular one
     
     Parameters
     -----------
     e0 : double
        Initial eccentricity
     
     Returns
     --------
     tovertc : double
        Ratio of inspiral time for a general eccentric orbit to that of a circular one

     */
    //constants are in constants.h
    double fitToToverTc = POLYNOMIAL_TOVERTC_FIT_C0 + (POLYNOMIAL_TOVERTC_FIT_C1 * e0) + (POLYNOMIAL_TOVERTC_FIT_C2 * e0 * e0) + (POLYNOMIAL_TOVERTC_FIT_C3 * e0 * e0 * e0) + (POLYNOMIAL_TOVERTC_FIT_C4 * e0 * e0 * e0 * e0) + (POLYNOMIAL_TOVERTC_FIT_C5 * e0 * e0 * e0 * e0 * e0) + (POLYNOMIAL_TOVERTC_FIT_C6 * e0 * e0 * e0 * e0 * e0 * e0) + (POLYNOMIAL_TOVERTC_FIT_C7 * e0 * e0 * e0 * e0 * e0 * e0 * e0) + (POLYNOMIAL_TOVERTC_FIT_C8 * e0 * e0 * e0 * e0 * e0 * e0 * e0 * e0) + (POLYNOMIAL_TOVERTC_FIT_C9 * e0 * e0 * e0 * e0 * e0 * e0 * e0 * e0 * e0);
    
    // Fix edge cases
    if(fitToToverTc > 1.0){
        std::cerr << "t/tc > 1.0" << std::endl;
        return 1.0;
    }
    else if(fitToToverTc < 0.0){
        std::cerr << "t/tc < 0.0" << std::endl;
        return 0.0;
    }
    else{
        return fitToToverTc;
    }
}


///////////////////////////////////////////////////////
//          UNUSED FUNCTIONS
///////////////////////////////////////////////////////

//double calculateC0(double a0=AU, double e0=0){
//    // Calculate the constant c0 based on initial conditions a0 and e0.
//    return a0*(1-e0*e0)*pow((1 + (121.0/304.0)*e0*e0), (-870.0/2299.0))*pow(e0, (-12.0/19.0));
//}
//

//double tovertcinterpolated(double e0, gsl_interp_accel *my_accel_ptr, gsl_spline *my_spline_ptr){
//    // Ratio of time to coalescence over the time for a circular orbit
//    double squareBrackets = (1 + (121.0/304.0)*e0*e0);
//    double integral = gsl_spline_eval(my_spline_ptr, e0, my_accel_ptr);
//    return (48.0/19.0)*pow((1 - e0*e0), 4.0)*pow(squareBrackets, (-3480.0/2299.0))*pow(e0, (-48.0/19.0))*integral;
//}

//double timeToCoalescence(double a0=AU, double e0=0.0, double m1=Msol, double m2=Msol){
//    // Calculate the time to coalescence of the eccentric orbit due to emission in gravitational radiation
//    double c0   = calculateC0(a0, e0);
//    double beta = calculateBeta(m1, m2);
//    double t    = 0;
//    if(e0 == 0){
//        t = pow(a0, 4)/(4.0*beta);
//    }
//    else{
//        t = (12.0/19.0)*c0*c0*c0*c0*e_integral(0, e0)/beta; // and integral
//    }
//    return t;
//}
//
//double timeToCoalescenceApprox1(double a0=AU, double e0=0.0, double m1=Msol, double m2=Msol){
//    // Calculate the time to coalescence using the approximation given
//    // Supposedly valid for LARGE e ~ 1
//    double c0       = calculateC0(a0, e0);              // Get the constant c0
//    double beta     = calculateBeta(m1, m2);            // Get the constant beta
//    double tovertc  = pow((1-e0*e0), (7.0/2.0));        // Ratio of time to circular time using approximation.
//    double t        = 0;                                // Initialise t
//    t               = c0*c0*c0*c0/(4*beta)*tovertc;     // We want the actual time to coalesce
//    return t;
//}
//
//double timeToCoalescenceApprox2(double a0=AU, double e0=0.0, double m1=Msol, double m2=Msol){
//    // Calculate the time to coalescence using the approximation given
//    // Supposedly valid for LARGE e ~ 1
//    double c0       = calculateC0(a0, e0);              // Get the constant c0
//    double beta     = calculateBeta(m1, m2);            // Get the constant beta
//    double tovertc  = (768.0/425.0)*pow((1-e0*e0), (7.0/2.0));                                // Ratio of time to circular time using approximation.
//    double t        = 0;                                // Initialise t
//    t               = c0*c0*c0*c0/(4*beta)*tovertc;     // We want the actual time to coalesce
//    return t;
//}



///////////////////////////////////////////////////////
//                  NOTES
///////////////////////////////////////////////////////

// here we have already declared and initialised pointers
// now want to prepare the pointers (which is the whole point of this function)
// Interpolate arrays:
// Probably don't need to output to user here
// output data to user
//ofstream outs("/Users/simons/Documents/C++/linearInterpolation/Results/NumericalResults.txt");
//
// Header for output
//outs << "i\te\tintegral\tt/tc" << endl;
//
// Save results to file -- not needed here
//outs << i << TAB << e[i] << TAB << Integral[i] << TAB << tOverTcArray[i] << NL;


//
//    my_accel_ptr  = gsl_interp_accel_alloc();                                       // Allocate 'interpolation accelerator'
//    my_spline_ptr = gsl_spline_alloc(gsl_interp_cspline, nPointsActual);            // Allocate 'spline object'

// Initialise the spline for the results of the Integral:
//gsl_spline_init(my_spline_ptr, e, Integral, nPointsActual);

//double tOverTcInterpArray[nPointsActual];                       // Declare array to store t/tc
//tOverTcInterpArray[i]   = 0;                                // Intialise array for t/tc using interpolation
