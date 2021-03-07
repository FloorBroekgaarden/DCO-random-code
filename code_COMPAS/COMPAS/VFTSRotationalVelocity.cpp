#include "VFTSRotationalVelocity.h"

double ramirezAgudelo2013OStarRotationalVelocityAnalyticCDF(double ve){
	/*
	Calculate the analytic cumulative distribution function (CDF) for the 
	equatorial rotational velocity of single O stars from Equation 1-4 of
	Ramirez-Agudelo et al 2013 https://arxiv.org/abs/1309.2929

	Modelled as a mixture of a gamma component and a normal component

	Parameters
	-----------
	ve : double
		Rotational velocity (in km s^-1) at which to calculate cdf
	
	Returns
	--------
	
	*/

    double alpha = 4.82;
    double beta = 1.0/25.0;
    double mu = 205.0; //km s^-1
    double sigma = 190.0; //km s^-1
    double Igamma = 0.43;

    // Use the boost gamma function from
    // https://www.boost.org/doc/libs/1_53_0/libs/math/doc/sf_and_dist/html/math_toolkit/dist/dist_ref/dists/inverse_gamma_dist.html

    boost::math::inverse_gamma_distribution<> gammaComponent(alpha, beta); // (shape, scale) = (alpha, beta)

    // Use the boost normal distribution from
    // https://www.boost.org/doc/libs/1_47_0/libs/math/doc/sf_and_dist/html/math_toolkit/dist/dist_ref/dists/normal_dist.html

    boost::math::normal_distribution<> normalComponent(mu, sigma);

	return Igamma * boost::math::cdf(gammaComponent, ve) + (1.0 - Igamma) * boost::math::cdf(normalComponent, ve);

}

struct ramirezAgudelo2013OStarRotationalVelocityParams {
    double u; // Value of CDF, draw in U(0,1)
};

typedef struct ramirezAgudelo2013OStarRotationalVelocityParams RotationalVelocityParams;

double ramirezAgudelo2013OStarRotationalVelocityAnalyticCDFRootFind(double ve, void * p){
    /*
    Calculate the inverse of ramirezAgudelo2013OStarRotationalVelocityAnalyticCDF 

    Parameters
    -----------
    x : double
        Value of the kick vk which we want to find
    p : RotationalVelocityParams 
        Structure containing y, the CDF draw U(0,1)
    
    Returns
    --------
    z : double
        Should be zero when x = vk, the value of the kick to draw
    */
    
    RotationalVelocityParams * params = (RotationalVelocityParams *)p;
    
    double y = (params->u);             // Value of CDF, should be drawn as U(0,1)
    
    return ramirezAgudelo2013OStarRotationalVelocityAnalyticCDF(ve) - y;
}

double ramirezAgudelo2013OStarRotationalVelocityAnalyticCDFInverseSampling(double u, double xmin, double xmax){
	/*
	Use inverse sampling and root finding to draw a rotational velocity from the CDF
	
	Parameters
	-----------
	u : double
		Random number in (0,1)
	xmin : double
		Minimum value to search for root in
	xmax : double
		Maximum value to search for root in
	
	Returns
	--------
	ve : double
		Rotational velocity in km s^-1

	*/
	bool debugging = false;
    
    int status;
    int iter = 0, max_iter = 100;
    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;
    double result = 0;

	gsl_function F;

    double maximum_inverse = ramirezAgudelo2013OStarRotationalVelocityAnalyticCDF(xmax);
    double minimum_inverse = ramirezAgudelo2013OStarRotationalVelocityAnalyticCDF(xmin);

    if(debugging){

    	std::cout << "u, cdf(xmin), cdf(xmax) = " << u << " " << minimum_inverse << " " << maximum_inverse << std::endl;

    }

    while(u > maximum_inverse){

        xmax *= 2.0;

        maximum_inverse = ramirezAgudelo2013OStarRotationalVelocityAnalyticCDF(xmax);

        if(debugging){
            std::cout << "random number past end point, increasing end point" << std::endl;
            std::cout << "new maximum inverse is " << maximum_inverse << std::endl;
        }

    }

    if(u > maximum_inverse){

        u = maximum_inverse;

        if(debugging){
            std::cout << "random number past end point, resetting to maximum" << std::endl;
        }
    
    }

    if(u < minimum_inverse){

    	result = xmin;

    }
    else{

    	RotationalVelocityParams params = {u}; // u

	    F.function = &ramirezAgudelo2013OStarRotationalVelocityAnalyticCDFRootFind;
	    F.params = &params;

	    // gsl_root_fsolver_brent
	    // gsl_root_fsolver_bisection
	    T = gsl_root_fsolver_brent;
	    s = gsl_root_fsolver_alloc (T);

	    gsl_root_fsolver_set (s, &F, xmin, xmax);

	    status = GSL_CONTINUE;

    	while(status == GSL_CONTINUE and iter < max_iter){

        	//std::cout << "in while loop" << std::endl;

        	iter++;
        	status = gsl_root_fsolver_iterate (s);
        	result = gsl_root_fsolver_root (s);
        	xmin = gsl_root_fsolver_x_lower (s);
        	xmax = gsl_root_fsolver_x_upper (s);
        	status = gsl_root_test_interval (xmin, xmax, 0, 0.001);

        	// if(debugging and status == GSL_SUCCESS){
        	// 	std::cout << "Converged:\n" << std::endl;
        	// 	std::cout << iter << " " << xmin << " " << xmax << " " << r << " " << xmax - xmin << std::endl;
        	// }

        }

        // De-allocate memory for root solver
    	gsl_root_fsolver_free (s);

    } 

    if(debugging){
        std::cout << "x = " << result << std::endl;
    }

    return result;
}

// SIMON: I tried to add a generic inverse_sampling_from_cdf_using_root_finding to inversionSampling
// that could be used for this and drawing kicks from maxwellCDF but I
// couldn't work out how to pass it a generic function with arbitrary structure of parameters.
// I'll try again another time.