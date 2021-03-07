#include "lambdaKruckow.h"

double lambdaKruckow(double Radius, double powerAlpha){
	/*
	Proposed fit for the common envelope lambda parameter from Fig. 1 of Kruckow et al. 2016 (arXiv:1610.04417)
	
	Spectrum fit to the regiou bounded by the upper and lower limits as shown in Kruckow+ 2016.
	Made by Alejandro Vigna-Gomez
	lambda(Radius, Alpha) = 1600*(0.00125^-\[Alpha])*(Radius^\[Alpha])
	-1.0 <= Alpha <= -2.0/3

	Parameters
	-----------
	Radius : double
		Radius in solar radii
	powerAlpha : double
		Power 

	Returns
	--------
	lambda : double
		Common envelope lambda parameter
	*/
	bool    debugging = false;
	//debugging = true;

	// Setting the limits for Alpha
	if(powerAlpha < -1.0){powerAlpha = -1.0;}
	if(powerAlpha > -2.0/3.0){powerAlpha = -2.0/3.0;}

	double  lambdaKruckow = 1600.0*pow(0.00125,-powerAlpha)*pow(Radius, powerAlpha);

    if(debugging){
        std::cout << "Radius, A, powerAlpha, lambda: " << Radius << "\t" << 1600.0*pow(0.00125,-powerAlpha) << "\t" << powerAlpha << "\t"<< lambdaKruckow << std::endl;
    }

    return  lambdaKruckow;
}