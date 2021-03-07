//
//
//
//
//  Created by Coen Neijssel 27/11/2016
//
//
// This file contains presctriptions from Dr. Bernhard Muller
// that describe the remnant mass, type, and kick of given
// the CO-Core mass at the time of the supernova.
// (see Muller et al 2016) (add link)
//
//

#include "SNMuller2016.h"

int MullerRemnantType(double COCoreMass,  unsigned long randomSeed){

    int		stellarType;
	bool	debugging = false;

    if (COCoreMass < Mch){
		stellarType = OXYGEN_NEON_WHITE_DWARF;
		std::cerr << randomSeed << "\tError in MullerRemnantType function. COCoreMass < Mch"
                                   " and still considered as a potential NS or BH." << std::endl;
	}
    if (COCoreMass < 3.6){stellarType =NEUTRON_STAR;}

    else if (COCoreMass < 4.05){stellarType =BLACK_HOLE;}   

    else if (COCoreMass < 4.6){stellarType =NEUTRON_STAR;}

    else if (COCoreMass < 5.7){stellarType =BLACK_HOLE;}

    else if (COCoreMass < 6.0){stellarType =NEUTRON_STAR;}

    else if (COCoreMass >= 6.0){stellarType =BLACK_HOLE;}

    	else{
		std::cerr << randomSeed << "\tError in MullerRemnantType function."
                                   " Shouldn't get here." << std::endl;
	}

	if(debugging){std::cout << "remant type:\t" << stellarType << std::endl;}

    return stellarType;
}//closing MullerRemnantType



void MullerRemnant(double COCoreMass, double *Mass,  unsigned long randomSeed){
	// ALEJANDRO - 14/03/2017 - Updated prescription of Bernhard Muller's models, incorporating lower metallicity models for core masses 1.372 <= m_{C/O} < 1.65
	// Where m_{C/O} = M_{C/O}/M_{\odot}
	// First approach is to model the lower limit using Chnadrasekhar mass, m_{C/O} > Mch = 1.44
    double	remnantMass=1.4; 					// Limit mass for a White Dwarf units Msun.
	double 	neutrinoLossFallbackFactor = 1.0;	// Factor which accounts for mass loss in neutrino winds during a supernovae. Should be made a flag and added to pythonSubmit.py
	bool	debugging = false;

    if (COCoreMass < 1.44){
		if(debugging){std::cout << "Inside MullerRemnantMass. COCoreMass < 1.44, therefore remnantMass = 1.4, limit mass for a White Dwarf." << std::endl;} 
	}
	else if (COCoreMass < 1.49){
		remnantMass = 1.21 - 0.4*(COCoreMass - 1.372);
	}
	else if (COCoreMass < 1.65){
		remnantMass = 1.16;
    }
    else if (COCoreMass < 2.4){
		remnantMass = 1.32 + 0.3*(COCoreMass - 1.65);
    }
    else if (COCoreMass < 3.2){
		remnantMass = 1.42 + 0.7*(COCoreMass - 2.4);
    }
    else if (COCoreMass < 3.6){
		remnantMass = 1.32 + 0.25*(COCoreMass - 3.2);
    }
    else if (COCoreMass < 4.05){
		// Going to be a black-hole
		remnantMass = *Mass * neutrinoLossFallbackFactor;
    }
    else if (COCoreMass < 4.6){
		remnantMass = 1.5;
    }
    else if (COCoreMass < 5.7){
		// Going to be a black-hole
		remnantMass = *Mass * neutrinoLossFallbackFactor;
    }
    else if (COCoreMass < 6.0){
		remnantMass = 1.64 - 0.2*(COCoreMass - 5.7);
    }
    else if (COCoreMass >= 6.0){
		// Going to be a black-hole
		remnantMass = *Mass * neutrinoLossFallbackFactor;
    }
	else{
		std::cerr << randomSeed << "\tError in MullerRemnantMass function. Shouldn't get here." << std::endl;
	}
	
	if(debugging){std::cout << "remant mass:\t" << remnantMass << std::endl;}

    *Mass = remnantMass;
}



double MullerRemnantKick(double COCoreMass,  unsigned long randomSeed){
    //for BH assume complete fallback so no ejecta hence no rocket effect
	// If BH doesnt assumes complete fallback, e.g. neutrino mass loss, 
    //a Blauuw Kick should be calculated
	// ALEJANDRO - 17/03/2017 - It seems from Bernhard's notes he gives us 
    //V_kick = V_kick_3D. Checking with some of the highest kicks he gives, we get V_kick>1000 km s^-1
    
    double	remnantKick = NEVER_SET;	//units km/s
	double	lowerRegimeKick = 70.0;		// Bernhard proposes to use 10 km s^-1 to replicate ECSN. 
                                        //He quotes Bray & Eldridge 2016 on this.
	double	oneDimKick = sqrt(3.0);		// Used to change given 3D kick to 1D kick, 
                                        //where VKick_1D = VKick_3D/sqrt(3)
	double	remnantKick1D = NEVER_SET;
	bool	debugging = false;
	
	
	if (COCoreMass < 1.44){
		if(debugging){std::cout << "Inside MullerRemnantKick. kick = 0.0 km s^-1." << std::endl;}
		remnantKick = 0.0;
	}
	else if (COCoreMass < 1.49){
		remnantKick = lowerRegimeKick + 2000.0*(COCoreMass - 1.372);
	}
    else if (COCoreMass < 1.65){
		remnantKick = 180.0 + 1300.0*(COCoreMass - 1.49);
    }
	else if (COCoreMass < 2.4){
		remnantKick = 250.0 + 350.0*(COCoreMass - 1.65);
    }
    else if (COCoreMass < 3.2){
		remnantKick = 400.0 + 1100.0*(COCoreMass - 2.4);
    }
    else if (COCoreMass < 3.6){
		remnantKick = 160.0 + 240.0*(COCoreMass - 3.2);
    }
    else if (COCoreMass < 4.05){
		remnantKick = 0.0; 
    }
    else if (COCoreMass < 4.6){
		remnantKick = 700.0 + 100.0*(COCoreMass - 4.05);
    }
    else if (COCoreMass < 5.7){
		remnantKick = 0.0;
    }
    else if (COCoreMass < 6.0){
		remnantKick = 550.0 - 600.0*(COCoreMass - 5.7);
    }
    else if (COCoreMass >= 6.0){
		remnantKick = 0.0;
    }
	else{
		std::cerr << randomSeed << "\tError in MullerRemnantKick function."
                                   " Shouldn't get here." << std::endl;		
	}
	
	if(debugging){
        std::cout << "remant kick 3D:\t" << remnantKick << std::endl;
        std::cout << "CO COreMass:\t" << remnantKick << std::endl;
    }
    
    return remnantKick;
}//closing MullerRemnantKick


