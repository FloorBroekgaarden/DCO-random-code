#include "NSMassRadiusRelation.h"

// implementation goes here

double momentOfInertia(double mass, double radius){
    /*
    Model independent relation between moment of inertia and mass and 
    radius of a neutron star
    
    Equation 8 in Raithel et al 2016 https://arxiv.org/abs/1603.06594
    
    https://tap.arizona.edu/sites/tap.arizona.edu/files/Raithel_2015_manuscript.pdf
    
    Parameters
    ------------
    M : float
        Mass in solar masses
    R : float
        Radius in km
    
    Returns
    ---------
    I : float
        Moment of inertia in g cm^2
    */
   // std::cout << "Mass, radius = " << mass << " " << radius << std::endl;
    return 0.237 * mass * MsolTog * pow((radius * km_in_cm), 2.0) * (1.0 + (4.2 * (mass/radius)) + 90 * pow((mass/radius), 4.0));
}

double interpolatedMassRadiusRelation(double mass, const double massArray[], const double radiusArray[], const int nRowsArray){
	/*
	Interpolate neutron star radius from tabulate mass-radius relation corresponding to a
	given neutron star equation of state.

	Parameters
	-----------
	mass : double
		Mass at which to calculate interpolated radius in solar masses
	massArray : double
		Array of masses in km for a given equation of state
	radiusArray : double
		Array of radii in km for a given equation of state
	nRowsArray : int 
		Number of elements in massArray and radiusArray for a given equation of state

	Returns
	--------
	radiusInterpolated : double
		Neutron star radius interpolated from tabulated equation of state in km

	*/

	// Allocate memory to the gsl objects
	gsl_interp *interpolation = gsl_interp_alloc (gsl_interp_linear, nRowsArray);
	gsl_interp_accel * accelerator =  gsl_interp_accel_alloc();

	// Initialise interpolation
	gsl_interp_init(interpolation, massArray, radiusArray, nRowsArray);

	// Evaluate interpolant at desired mass
	double interpolatedRadius = gsl_interp_eval(interpolation, massArray, radiusArray, mass, accelerator);

	// Free memory
	gsl_interp_free (interpolation);
	gsl_interp_accel_free (accelerator);

	return interpolatedRadius;
}

double neutronStarRadius(double mass, const programOptions &options){
	/*
	Wrapper function to calculate neutron star radius according to selected equation of state

	pass program options

	Parameters
	----------
	mass : double
		Neutron star mass in Msol
	options : programOptions
		User specified program options

	Returns
	-------
	radius : double
		Neutron star radius in km

	*/
	if(options.neutronStarEquationOfState == NS_EOS_SSE){
		return 10.0;
	}
	else if(options.neutronStarEquationOfState == NS_EOS_ARP3){

		// We don't extrapolate so masses outside table just set to extreme values
		if(mass < mass_radius_relation_ARP3_minimum_mass){
			return mass_radius_relation_ARP3_radius[0];
		}
		else if(mass > mass_radius_relation_ARP3_maximum_mass){
			return mass_radius_relation_ARP3_radius[nrows_ARP3-1];
		}
		else{
			return interpolatedMassRadiusRelation(mass, mass_radius_relation_ARP3_mass, mass_radius_relation_ARP3_radius, nrows_ARP3);
		}
	}
	else{
		std::cerr << "Unrecognised NS EOS -- using default NS radius" << std::endl;
		return 10.0;
	}
}
