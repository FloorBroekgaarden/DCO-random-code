//
//  timeToCoalescence.cpp


#include "timeToCoalescence.h"

double timeToCoalescenceUsingInterpolation(double a0, double e0, double m1, double m2){
    /*
     Calculate the time to coalescence for a binary with arbitrary eccentricity using interpolation
     
     This is Equation 5.14 in Peters 1964 http://journals.aps.org/pr/pdf/10.1103/PhysRev.136.B1224

     Parameters
     -----------
     a0 : double
        Initial semi-major axis in SI units
     e0 : double
        Initial eccentricity
     m1 : double 
        Primary mass in SI units
     m2 : double
        Secondary mass in SI units
     
     Returns
     --------
     t_coalesce : double
        Time to coalescence in SI units (s)

     */
    
    // Declare some variables
    double tc       = 0.0;                                          // Time for a circular binary to coalesce
    double beta     = 0.0;                                          // Beta constant (function of masses)
    double tovertc  = 0.0;                                          // Result of interpolation of t/tc function
	bool	debugging = false;
	
	if(debugging){
		std::cout << "a0 [AU]\t" << a0/AU << std::endl;
		std::cout << "e0\t" << e0 << std::endl;
		std::cout << "m1 [Msol]\t" << m1/Msol << std::endl;
		std::cout << "m2 [Msol]\t" << m2/Msol << std::endl;
	}
    // Calculate circular time to coalescence
    beta = calculateBeta(m1, m2);                                   // Masses should be in SI
    tc = a0*a0*a0*a0/(4.0*beta);                                    // Calculate time for a circular binary to merge
    
    // calculate t/tc using the interpolated function
    //tovertc = polynomialFitToToverTc(e0);                           // Ratio of inspiral time to that of circular system
    if((e0==0.0)or(e0==0)){
		if(debugging){
			std::cout << "tc\t" << tc <<std::endl;
		}
		return	tc;
	}
				
    double c0=a0*(1.0-e0*e0)*pow(e0,-12.0/19.0)*pow(1.0+(121.0*e0*e0/304.0), -870.0/2299.0);
		
		
    if(e0<0.01){ 
		if(debugging){
			std::cout << "e0<0.01" << std::endl;
			std::cout << "beta\t" << beta << std::endl;
			std::cout << "tc\t" << tc << std::endl;
			std::cout << "a0\t" << a0 << std::endl;
			std::cout << "c0\t" << c0 << std::endl;
			std::cout << "t [s]\t" << c0*c0*c0*c0*pow(e0,48.0/19.0)/(4.0*beta) << std::endl;
		}

		return c0*c0*c0*c0*pow(e0,48.0/19.0)/(4.0*beta);
    }

    if(e0>0.99){
		// Approximation of eq. 5.14 of Peters 1964, for high eccentricities
		return (768.0/425.0)*tc*pow((1.0-(e0*e0)),3.5);
    }	

    double sum=0;
    double de=e0/10000;
    
	for(double e=0; e<e0; e=e+de){
		sum=sum+de*pow(e,29.0/19.0)*pow((1.0+(121.0/304.0)*e*e),1181.0/2299.0)/pow((1-e*e),1.5);
    }
	
	if(debugging){
		std::cout << "e0>=0.01" << std::endl;
		std::cout << "beta\t" << beta << std::endl;
		std::cout << "tc\t" << tc << std::endl;
		std::cout << "t [s]\t" << 12.0/19.0*c0*c0*c0*c0/beta*sum << std::endl;
	}
		
    return 12.0/19.0*c0*c0*c0*c0/beta*sum;

    //return tovertc * tc;                                            // t /tc * tc = t, the actual time to coalescence, right?
}

double calculateBeta(double m1, double m2){
    /*
     // Calculate the constant beta based on masses as defined in Equation 5.9 in Peters 1964 http://journals.aps.org/pr/pdf/10.1103/PhysRev.136.B1224
     
     Parameters
     -----------
     m1 : float
        Primary mass in SI units
     m2 : float
        Secondary mass in SI units
     
     Returns
     --------
     beta : float
        Mass parameter beta
 
     */
    double 	M = m1 + m2;
	bool	debugging = false;

	if(debugging)
		std::cout << "beta\t" << (64.0/5.0)*G*G*G*m1*m2*M*pow(c,-5.0) << std::endl;

    return (64.0/5.0)*G*G*G*m1*m2*M*pow(c,-5.0);
}


//////////////////////////////////////////////////////
//                  NOTES/OLD STUFF
//////////////////////////////////////////////////////

//double timeToCoalescenceUsingInterpolation(double a0, double e0, double m1, double m2, gsl_interp_accel *my_accel_ptr, gsl_spline *my_spline_ptr){
//    // Calculate the time to coalescene for a binary with arbitrary eccentricity, using interpolation
//    // Result is in SI (s)
//
//    // Declare some variables
//    double tc       = 0.0;                                            // Time for a circular binary to coalesce
//    double beta     = 0.0;                                            // Beta constant (function of masses)
//    double tovertc  = 0.0;                                            // Result of interpolation of t/tc function
//
//    // Calculate circular time to coalescence
//    beta = calculateBeta(m1, m2);                                   // Masses should be in SI
//    tc = a0*a0*a0*a0/(4.0*beta);                                    // Calculate time for a circular binary to merge
//
//    // calculate t/tc using the interpolated function
//    tovertc = gsl_spline_eval(my_spline_ptr, e0, my_accel_ptr);     // Interpolated t/tc
//    if(tovertc < 0){ tovertc = 0;}                                  // KLUDGE fix to fix negative values
//
//    return tovertc * tc;                                            // t /tc * tc = t, the actual time to coalescence, right?
//
//}
