//
//  inversionSampling.cpp

#include "inversionSampling.h"

double inverse_sampling_from_tabulated_cdf(double u, const double x[], const double cdf[], const int nrows, const double xmax, const double xmin, const double cdfmax){
	/*
    Draw samples from a distribution defined by a tabulated cdf

    Uses simple linear interpolation

    Note that this function creates a copy of x and cdf (both of type std::vector<double>) so that 
    it can check if the cdf reaches 1, and renormalise it if it doesn't
	
    Parameters
    -----------
    u : double
        Random number between 0 and 1
	x : const double[]
		Variable for which you want to sample the inverse of cdf(x)
	cdf : const double[]
		Tabulated cdf(x)

	Returns
	---------
	xinterp : double
		Value of x corresponding to cdf(xinterp) = u
	*/

	bool debugging = false; //true;
	
	// Scale y axix to be between 0 and cdfmax
	u *= cdfmax; // uniform in (0, cdfmax)
	
    // Declare and initialise variables
    double xbelow = 0;
    double xabove = 0;
    
    double cdfbelow = 0;
    double cdfabove = 0;

    int ibelow = 0;
    int iabove = 0;

    double gradient = 0.0;
    double xinterp = 0.0;
    
    int i = 0;
    
    while(cdfabove < u){
        
        ibelow = i - 1; // will be -1 for i = 0?
        iabove = i;
        
        if(i > 0){
	        xbelow = x[ibelow];  
	        xabove = x[iabove];

	        cdfbelow = cdf[ibelow];
	        cdfabove = cdf[iabove];
	    }
        
        if(debugging){
		    std::cout << "below : i, x, cdf: " << ibelow << " " << xbelow << " " << cdfbelow << std::endl;
			std::cout << "above : i, x, cdf: " << iabove << " " << xabove << " " << cdfabove << std::endl;
		}

        i = i + 1;
    }

	gradient = (cdfabove - cdfbelow)/(xabove - xbelow);
	xinterp = xbelow + (1.0/gradient * (u - cdfbelow));

    if(debugging){
        std::cout << "u = " << u << std::endl;
        std::cout << "gradient = " << gradient << std::endl;
        std::cout << "xinterp = " << xinterp << std::endl;
    }
    
    return xinterp;

}

double inverse_sampling_from_power_law(const gsl_rng *r, double power, double xmax, double xmin){
    
    // This function draws samples from a power law distribution p(x) ~ x^(n) between xmin and xmax
    // Checked against python code for the same function
    
    double u = gsl_rng_uniform(r);              // Draw a random number between 0 and 1
    
    if(power == -1.0){
        return exp(u*log(xmax/xmin))*xmin;
    }
    else{
        return pow((u*(pow(xmax, power+1.0) - pow(xmin, power+1.0)) + pow(xmin, power+1.0)), 1.0/(power+1.0));
    }
}

double sampleFromBrokenPowerLaw(const gsl_rng *r, double xMax, double xMin, std::vector<double> indices, std::vector<double> breaks)
{
	//draws a sample from a broken power law, given a vector of power law indiced and where the breaks in the power law are.
	//requires both an xmin and xmax so that the sampling area is finite

	//check that we're given correctly sized lists
	if (indices.size() != breaks.size()+1)
	{
		std::cerr << "index lists and breaks list passed to broken power law sampler incorrectly shaped";
		throw;
	}

	//trim down so that we only have indices/breaks above xMin and below xMax, and add xMin and xMax so we have a list 
	//of all the limits
	std::vector<double> limits, trimmedIndices;
	limits.push_back(xMin);
	for (int i=0; i<breaks.size(); ++i)
	{
		if (breaks[i] > xMin and breaks[i] < xMax)
		{
			limits.push_back(breaks[i]);
			trimmedIndices.push_back(indices[i]);
		}

	}
	limits.push_back(xMax);
	//need to do the last index manually
	if (xMin > breaks[breaks.size()-1]) trimmedIndices.push_back(indices[indices.size()-1]);


	//calculate the integral in each region and compute the CDF
	//this could be sped up by only doing this once, but not exactly expensive like this
	double total = 0.;
	std::vector<double> cdf;
	for (int i=0; i<limits.size()-1; i++)
	{
		double a = limits[i];
		double b = limits[i+1];

		double oneMinusAlpha = 1. - trimmedIndices[i];

		double integral = (std::pow(b,oneMinusAlpha) - std::pow(a,oneMinusAlpha))/oneMinusAlpha;

		total += integral; 

		cdf.push_back(total);
	}


	//normalise the cdf
	for (int i=0; i<cdf.size(); i++)
	{
		cdf[i] /= cdf[cdf.size()-1];
	}

	//draw a random number
	double u = gsl_rng_uniform(r);

	//and draw from the relevant power law
	for (int i =0; i< cdf.size(); i++)
	{
		if (u < cdf[i]) return inverse_sampling_from_power_law(r,-1.*trimmedIndices[i],limits[i],limits[i+1]);
	}

	std::cerr<<"something went wrong with the broken power law sampling, shouldn't get here"<<std::endl;
	throw;

}


// SIMON: I tried to add a generic inverse_sampling_from_cdf_using_root_finding to inversionSampling
// that could be used for this and drawing kicks from maxwellCDF but I
// couldn't work out how to pass it a generic function with arbitrary structure of parameters.
// I'll try again another time.
