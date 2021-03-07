//
//  supernovaFunctions.cpp

#include "supernovaFunctions.h"

void supernovaMain(double *totalMassAtCompactObjectFormation, double *HeCoreMassAtCompactObjectFormation,
                   double *COCoreMassAtCompactObjectFormation,double *coreMassAtCompactObjectFormation,
                   double *Mass0, double *Age, double *envMass, double *coreMass, 
                   double *HeCoreMass, double *COCoreMass, double *coreRadius,
                   double *Mass, double MZAMS, double McBAGB, double *fallbackFraction,
                   int *stellarType, double *Radius, double *Luminosity, 
                   double *Temperature, bool *flagSN, bool *flagECSN, bool *flagExperiencedECSN,
                   bool *flagExperiencedCCSN, bool *flagExperiencedPISN, bool *flagExperiencedPPISN, double *drawnKickVelocity, double *kickVelocity, 
                   bool *error, bool *flagHrichSN, bool *flagHpoorSN, int fryerSNengine, unsigned long randomSeed,
                   double *angularMomentum, double *pulsarSpinPeriod, double *pulsarSpinFrequency, double *pulsarMagneticField, double *pulsarSpinDownRate, double *neutronStarMomentOfInertia,
                   const programOptions &options, const gsl_rng *r){

    //This is the main supernova function called in Star.cpp:
    //  - Decides the type of the supernova and passes the pointer
    //    variables to there to be updated.
    //    Returns error if no supernova type is found for given Mass
    //  - It delegates the parameters which need to be (re)set    
    
    bool debugging = false; //false; //true;

    supernovaSetFormationParameters(totalMassAtCompactObjectFormation, HeCoreMassAtCompactObjectFormation,
                                     COCoreMassAtCompactObjectFormation, coreMassAtCompactObjectFormation,
                                     *Mass, *coreMass, *HeCoreMass, *COCoreMass);
                                    //for extra safety I only pass Masses as value since we still need this.    
    

    //For TPAGB there is an exception, we dont use the Mass parameter to determine
    //supernova type, but the mass at Base Asympotic Giant Branch.
    double M = 0.0;
    if(*stellarType == THERMALLY_PULSING_ASYMPTOTIC_GIANT_BRANCH){ M = McBAGB;}
    else if(*stellarType == OXYGEN_NEON_WHITE_DWARF){M = 5.0;} //force ONeWD to ccSN, 5.0 doesn't change a physical
                                                          //parameter of the star.
    else{ M = *Mass;}//dereference for value

	// ALEJANDRO - 04/05/2018 - Check if the SN is H-rich or H-poor. For now, classify it for all possible SNe and not only CCSN forming NS.
	amountOfHydrogenInSupernova(stellarType, flagHrichSN, flagHpoorSN, randomSeed);

    if(options.usePulsationalPairInstability and (*HeCoreMass >= options.pulsationalPairInstabilityLowerLimit) and 
         (*HeCoreMass <= options.pulsationalPairInstabilityUpperLimit)){

        debugging = false; //true;
        
        if(debugging){
            std::cout << "Pulsational Pair Instability SuperNova" << std::endl;
        }

        supernovaPulsationalPairInstability(COCoreMass, HeCoreMass, envMass, Mass, options, flagSN, flagExperiencedPPISN, stellarType, Luminosity, Radius, Temperature, fallbackFraction);
        
    }
    else if(options.usePairInstabilitySupernovae and (*HeCoreMass >= options.pairInstabilityLowerLimit) and 
         (*HeCoreMass <= options.pairInstabilityUpperLimit)){

        debugging = false; //true;

        if(debugging){std::cout<<"Pair Instability SuperNova"<<std::endl;}
        
        supernovaPairInstability(Luminosity, Radius, drawnKickVelocity,
                                 stellarType, Temperature, kickVelocity, fallbackFraction, flagSN, 
                                 flagExperiencedPISN);

    }
    else if(M < 1.6){

        if(debugging){std::cout<<"Type IIa SuperNova"<<std::endl;}

        supernovaTypeIIa(stellarType, Mass, Radius, Luminosity, Temperature, 
                         drawnKickVelocity, kickVelocity);
    }

    else if(M >= 1.6 and M < 2.25){

        if(debugging){std::cout<<"Electron Capture SuperNova"<<std::endl;}
        supernovaElectronCapture(stellarType, Mass, Radius, Luminosity, drawnKickVelocity,
                                 kickVelocity, flagECSN, flagExperiencedECSN, Age, options);
    }
    
    else if(M >= 2.25){

        if(debugging){
            std::cout << "Core Collapse SuperNova" << std::endl;
            std::cout << "Mass = " << M << std::endl;
            std::cout << "COCoreMass = " << *COCoreMass << std::endl;
            std::cout << "HeCoreMass = " << *HeCoreMass << std::endl; 
        }

        supernovaCoreCollapse(envMass, coreMass, COCoreMass, Mass, MZAMS, fallbackFraction, 
                              stellarType, Radius, Luminosity, Temperature, flagSN, 
                              flagExperiencedCCSN, drawnKickVelocity, kickVelocity, 
                              error, fryerSNengine, randomSeed, options);
    }

    else{
        std::cerr << randomSeed << "\tError: No Supernova type found for star of mass= " << *Mass <<
                                  ",  in function supernovaMain() in supernovaFunctions.cpp" << std::endl;

        *error=true;
    }

    if(debugging){
        std::cout << "Seed---------------" << randomSeed << "\n"
                "in function supernovaMain() choosing supernova type for mass=:" << M << std::endl;
    }

    if(debugging){
        std::cout << "Total mass = " << *Mass << std::endl;
        std::cout << "COCoreMass = " << *COCoreMass << std::endl;
        std::cout << "HeCoreMass = " << *HeCoreMass << std::endl;
        std::cout << "envMass = " << *envMass << std::endl;
    }

    if(debugging){
        std::cout << "Exiting supernovaMain() with remnant properties:" << std::endl;
        std::cout << "encountered error =\t" << std::boolalpha << *error<<std::endl;
        std::cout << "remnant Type =\t" << *stellarType << std::endl;
        std::cout << "remnant Mass =\t" << *Mass << std::endl;
        std::cout << "remnant Radius =\t" << *Radius << std::endl;
        std::cout << "remnant Luminosity =\t" << *Luminosity << std::endl;
        std::cout << "remnant drawn kick =\t" << *drawnKickVelocity << std::endl;
        std::cout << "remnant kick =\t" << *kickVelocity << std::endl;
        std::cout << "flagSN = " << *flagSN << std::endl;
        std::cout << "flagPISN = " << *flagExperiencedPISN << std::endl;
    }
    
    assignPulsarParameters(stellarType, Mass, Radius, angularMomentum, pulsarSpinPeriod, pulsarSpinFrequency, pulsarMagneticField, pulsarSpinDownRate, neutronStarMomentOfInertia, options, r);
    
    if(debugging){
    	std::cout << "Birth magnetic field = " << *pulsarMagneticField << std::endl;
    	std::cout << "Birth spin period = " << *pulsarSpinPeriod << std::endl;
    	std::cout << "Birth spin frequency = " << *pulsarSpinFrequency << std::endl;
    }

    // SIMON : What is this line for? Why 'reset' things to 0?
    supernovaResetParameters(Mass0, Age, envMass, coreMass, HeCoreMass, COCoreMass, coreRadius);

}


void supernovaSetFormationParameters(double *totalMassAtCompactObjectFormation,
                                     double *HeCoreMassAtCompactObjectFormation,
                                     double *COCoreMassAtCompactObjectFormation,
                                     double *coreMassAtCompactObjectFormation,double Mass,
                                     double coreMass, double HeCoreMass, double COCoreMass){
    //Dereference here so I can assign value
    *totalMassAtCompactObjectFormation	= Mass;
    *HeCoreMassAtCompactObjectFormation = HeCoreMass;
    *COCoreMassAtCompactObjectFormation = COCoreMass;
    *coreMassAtCompactObjectFormation 	= coreMass; 
}

void supernovaResetParameters(double *Mass0, double *Age, double *envMass, double *coreMass, 
                              double *HeCoreMass, double *COCoreMass, double *coreRadius){
    //set values to zero
    //Dereference here so I can assign value
    *COCoreMass  = 0.0;
    *HeCoreMass  = 0.0;
    *coreMass    = 0.0;
    *Mass0       = 0.0;
    *envMass     = 0.0;
    *coreRadius  = 0.0;
    *Age         = 0.0;
}


void supernovaPairInstability(double *Luminosity, double *Radius, double *drawnKickVelocity,
                              int *stellarType, double *Temperature, double *kickVelocity, double *fallbackFraction, bool *flagSN, bool *flagExperiencedPISN){
    //Pair Instability Supernova
    //short handwavy story, The core is hot and massive enough that there is a significant 
    //amount of high-energetic gamma rays which can create an electron-positron pair.
    //When a significant amounts of pairs are created the photons taken away 
    //reduce the radiative pressure. The core contracts becomes hotter creating
    //more pairs. This runaway proces will end in a SN that explodes the entire core without
    //leaving a remnant.

    //std::cout << "In supernovaPairInstability" << std::endl;

    //Dereference here so I can assign value
    *Luminosity          = 0.0;
    *Radius              = 0.0;
    *Temperature         = 0.0;
    *drawnKickVelocity   = 0.0;
    *kickVelocity        = 0.0;
    *fallbackFraction    = 0.0;
    *stellarType         = MASSLESS_REMNANT;
    *flagExperiencedPISN = true; //This flag remains true to track the history of a star
    *flagSN              = true; // only remains true this timestep, makes sure to print to file, update orbit etc

    //std::cout << "stellarType remnant (should be 15) = " << *stellarType << std::endl;
    //std::cout << "flagPISN = " << *flagExperiencedPISN << std::endl;
    //std::cout << "flagSN = " << *flagSN << std::endl;
}

void woosley2017fit(double *Mass, double *HeCoreMass, const programOptions &options){
    /*
    This is a fit to the models of Woosley 2017 https://arxiv.org/abs/1608.08939
    */
    double Mbary = *HeCoreMass; // Strip off the hydrogen envelope if any was left (factor of 0.9 applied in blackHoleFormationNeutrinoMassLoss)
    *Mass = blackHoleFormationNeutrinoMassLoss(Mbary, options); // Convert to gravitational mass due to neutrino mass loss
}

void belczynski2016fit(double *Mass, const programOptions &options){
    /*
    This function implements the model from Belczynski et al 2016 https://arxiv.org/abs/1607.03116
    */
    double Mbary = 45.0; // (factor of 0.9 applied in blackHoleFormationNeutrinoMassLoss)
    *Mass = blackHoleFormationNeutrinoMassLoss(Mbary, options); // Convert to gravitational mass due to neutrino mass loss
    //*Mass = 0.9 * 45.0;
}

void marchant2018fit(double *Mass, double *HeCoreMass, const programOptions &options){
    /*
    Fit to the models of Marchant et al 2018 https://arxiv.org/abs/1810.13412
    */

    double ratioOfRemnantToHeCoreMass = -1.63057326e-08 * pow(*HeCoreMass, 7) + 
                                         5.36316755e-06 * pow(*HeCoreMass, 6) + 
                                        -7.52206933e-04 * pow(*HeCoreMass, 5) +
                                         5.83107626e-02 * pow(*HeCoreMass, 4) + 
                                        -2.69801221e+00 * pow(*HeCoreMass, 3) +
                                         7.45060098e+01 * pow(*HeCoreMass, 2) +
                                        -1.13694590e+03 * pow(*HeCoreMass, 1) +
                                         7.39643451e+03;

    if (ratioOfRemnantToHeCoreMass > 1.0){
        ratioOfRemnantToHeCoreMass = 1.0;
    }

    if (ratioOfRemnantToHeCoreMass < 0.0){
        ratioOfRemnantToHeCoreMass = 0.0;
    }
    
    double Mbary = ratioOfRemnantToHeCoreMass * *HeCoreMass;
    
    *Mass = blackHoleFormationNeutrinoMassLoss(Mbary, options); // Convert to gravitational mass due to neutrino mass loss
}

void supernovaPulsationalPairInstability(double *COCoreMass, double *HeCoreMass, double *envMass, double *Mass, const programOptions &options, bool *flagSN, bool *flagExperiencedPPISN, int *stellarType, double *Luminosity, double *Radius, double *Temperature, double *fallbackFraction){

    *flagSN               = true; // Flag so that the orbit gets updated and prints to supernovae.txt
    *flagExperiencedPPISN = true; // This flag remains true to track the history of a star

    if(options.pulsationalPairInstabilityPrescription == COMPASPPI){
        woosley2017fit(Mass, HeCoreMass, options);
    }
    else if(options.pulsationalPairInstabilityPrescription == STARTRACKPPI){
        belczynski2016fit(Mass, options);
    }
    else if(options.pulsationalPairInstabilityPrescription == MARCHANTPPI){
        marchant2018fit(Mass, HeCoreMass, options);
    }
    else{
        *Mass = *HeCoreMass; // Just use helium core mass
    }
    
    if(*Mass <= 0.0){
        *stellarType = MASSLESS_REMNANT;
    }
    else{
        *stellarType = BLACK_HOLE;
    }
    
    *Luminosity  = luminosityBlackHole(*Mass);                   // Luminosity of BH
    *Radius      = SchwarzschildRadiusBlackHole(*Mass);          // Schwarzschild radius (not correct for rotating BH)
    *Temperature = calculateTemperature(*Luminosity, *Radius);   // Temperature of BH
    *fallbackFraction = 1.0;        // Fraction of mass that falls back

}

void supernovaElectronCapture(int *stellarType, double *Mass, double *Radius, double *Luminosity, 
                              double *drawnKickVelocity, double *kickVelocity, bool *flagECSN, 
                              bool *flagExperiencedECSN, double *Age, const programOptions &options){
    //Electron Capture supernova
    //Short hand wavy story. The core ignites ONeMg and collapses through electron
    //capture. This core is less dense than the heavier Fe Core collapses and thus
    //people expect the neutrinos to escape easier resulting in a zero kick
    //TODO a nice reference (Nomoto 1984 for thorough discussion and a summary in Nomoto 1987)
    //Dereference here so I can assign value
    *stellarType       = NEUTRON_STAR;
    *Mass              = Mecs_rem;            //defined in constant.h
    //*Radius            = radiusNeutronStar(); //see remnantStructure.cpp
    *Radius            = neutronStarRadius(*Mass, options) * km_in_rsol; // neutronStarRadius in km
    *Luminosity        = luminosityNeutronStar(*Age, *Mass);
    *flagECSN            = true; // This flag only stays true during this timestep see Binarystar.cpp
    *flagExperiencedECSN = true; //This flag remains true to track the history of a star
}

void supernovaTypeIIa(int *stellarType, double *Mass, double *Radius, double *Luminosity, 
                       double *Temperature, double *drawnKickVelocity, double *kickVelocity){

    // Type IIa supernova -- leaves a massless remnant
    *stellarType       = MASSLESS_REMNANT;
    *Mass              = 0.0;
    *Radius            = 0.0; 
    *Luminosity        = 0.0;
    *Temperature       = 0.0;
    *drawnKickVelocity = 0.0;
    *kickVelocity      = 0.0;
}

void supernovaCoreCollapse(double *envMass, double *coreMass, double *COCoreMass, double *Mass, 
                           double MZAMS, double *fallbackFraction, int *stellarType, 
                           double *Radius, double *Luminosity, double *Temperature, bool *flagSN, 
                           bool *flagExperiencedCCSN, double *drawnKickVelocity, double *kickVelocity, 
                           bool *error, int fryerSNEngine, unsigned long randomSeed, const programOptions &options){
    
    //This function is called by supernovaMain()
    //It decides which prescription is used for the core collapse SN 
    //given by the program options.

    //It passes the pointers to the file containing the prescription
    //which updates the following parameters:
    //         Mass, stellarType, drawnKickVelocity, kickVelocity


    //At the end of this function we set the following parameters which are (so far) 
    // independent of the ccSN prescriptions (but do depend on the parameters above):
    //         Luminosity, Radius, Temperature, flagSN, flagexperiencedCCSN


    //Tip for pointers. *means passing value to function, no star means it is a pointer
    //except when it was already passed as a value (e.g. fryer SNEngine / randomseed
    //Pointer on pointers :p To point out the point of pointers

    bool    debugging = false; //true; //false;
    double	maximumNeutronStarMass = options.maximumNeutronStarMass; // MAXIMUM_NS_MASS;
    std::string remnantMassPrescription = options.remnantMassPrescription;
    if(debugging){std::cout<<"Seed---------------"<<randomSeed<<"\n"
                            "in function supernovaCoreCollapse() choosing"
                            " ccSN prescription type:"<<std::endl;}


    
    if(boost::iequals(remnantMassPrescription, "postitNote")){
        if(debugging){std::cout << "remnant prescription = SIMONSN" << std::endl;}
        *Mass = postItNotBlackHoleMass(MZAMS);
        *fallbackFraction = 0.0;     // Not defined
    }


    else if (boost::iequals(remnantMassPrescription, "hurley2000")){
        if(debugging){std::cout << "remnant prescription = HURLEYSN" << std::endl;}
        *Mass = neutronStarRemnantMass(*coreMass);
        *fallbackFraction = 0.0;     // Not defined
    }


    else if (boost::iequals(remnantMassPrescription, "belczynski2002")){
        if(debugging){std::cout << "remnant prescription = BELCZYNSKISN" << std::endl;}
                              //pointers are mass/fallback
        Belczynski2002Remnant(Mass, fallbackFraction, *coreMass); 
    }


    else if (boost::iequals(remnantMassPrescription, "fryer2012")){

        //debugging = true;

        if(debugging){
            std::cout << "remnant prescription = FRYERSN" << std::endl;
            std::cout << "Mass = " << *Mass << std::endl;
            std::cout << "COCoreMass = " << *COCoreMass << std::endl;
        }
        
        //pointers are mass/fallback
        FryerRemnant(Mass, *COCoreMass, fallbackFraction, fryerSNEngine, randomSeed, options);
        
        debugging = false;

    }
    
    else if ((boost::iequals(remnantMassPrescription, "muller2016"))||
             (boost::iequals(remnantMassPrescription, "muller2016Maxwellian"))){
        if(debugging){std::cout << "remnant prescription = MULLERSN" << std::endl;}
                     //pointer is Mass
        MullerRemnant(*COCoreMass, Mass, randomSeed);
    }

    
    else{
        *Mass = 0.0;
        *fallbackFraction = 0.0;
        std::cerr << randomSeed <<  "\tError: Unrecognised remnant mass prescription in "
                                    "supernovaCoreCollapse() in supernovaFunctions.cpp" << std::endl;
        *error= true;
    }

        // Set the stellar type (either use prescription or MAXIMUM_NS_MSS)
    if (boost::iequals(remnantMassPrescription, "muller2016")){
        *stellarType = MullerRemnantType(*COCoreMass, randomSeed);}
    else if(*Mass > maximumNeutronStarMass){
        *stellarType = BLACK_HOLE;}
    else if(*Mass <= maximumNeutronStarMass){
        *stellarType = NEUTRON_STAR;}

    else{
        std::cerr << randomSeed <<  "\tError: Unrecognised remnant mass type in supernovaRemnantMass function" << std::endl;}

    
    // Set radius, luminosity and temperature of compact remnant
    if (*stellarType == BLACK_HOLE){
        *Luminosity = luminosityBlackHole(*Mass);                  // Luminosity of BH
        *Radius = SchwarzschildRadiusBlackHole(*Mass);             // Schwarzschild radius (not correct for rotating BH)
        *Temperature = calculateTemperature(*Luminosity, *Radius);         // Temperature of BH

    }
    else if (*stellarType == NEUTRON_STAR){        
        *Luminosity = luminosityNeutronStar(0.0, *Mass);           // Luminosity of a NS as it cools
        //*Radius = radiusNeutronStar();                                   // Radius of NS set to 10km
        *Radius = neutronStarRadius(*Mass, options) * km_in_rsol;          // neutronStarRadius in km   
        *Temperature = calculateTemperature(*Luminosity, *Radius);         // Temperature of NS
        
    }
    else{
        std::cerr << randomSeed <<  "\tError: Can't set L,R,T of unknown remnant mass type in"
                                    "supernovaCoreCollapse() in supernovaFunctions.cpp" << std::endl;
        *error=true;
    }
     
    *flagSN = true;  // Flag only remains true this timestep see Binarystar.cpp
    *flagExperiencedCCSN = true; //Flag that remains true for the history of the star
    if(debugging){std::cout<<"exit - supernovaCoreCollapse()\n ------"<<std::endl;}

}


void amountOfHydrogenInSupernova(int *stellarType, bool *flagHrichSN, bool *flagHpoorSN, unsigned long randomSeed){
	// Function to classify if SNe are H-poor or H-rich. 
	// This is a simple classification depending on stellar type,
	// and we currently do it for all SN.
	
	bool	debugging = false;
	
	if(debugging){
		std::cout << randomSeed << "Inside amountOfHydrogenInSupernova function" << std::endl;
		std::cout << randomSeed << "\tStellartype:\t" << *stellarType << std::endl;
		std::cout << randomSeed << "\flagHrichSN:\t" << *flagHrichSN << std::endl;
		std::cout << randomSeed << "\flagHpoorSN:\t" << *flagHpoorSN << std::endl;
	}
	
	if(*stellarType <= THERMALLY_PULSING_ASYMPTOTIC_GIANT_BRANCH){
		*flagHrichSN = true;
	}
	else if(*stellarType >= NAKED_HELIUM_STAR_MS or *stellarType <= OXYGEN_NEON_WHITE_DWARF){
		*flagHpoorSN = true;		
	}
	else{
		std::cerr << randomSeed << "Supernova from a NS, BH or MR progenitor. Shouldn't get here." << std::endl;
	}
	
	if(debugging){
		std::cout << randomSeed << "Flags at the end of amountOfHydrogenInSupernova function" << std::endl;
		std::cout << randomSeed << "\flagHrichSN:\t" << *flagHrichSN << std::endl;
		std::cout << randomSeed << "\flagHpoorSN:\t" << *flagHpoorSN << std::endl;
	}	
	
}

void assignPulsarParameters(int *stellar_type, double *Mass, double *Radius, double *angularMomentum, double *pulsarSpinPeriod, double *pulsarSpinFrequency, double *pulsarMagneticField, double *pulsarSpinDownRate, double *neutronStarMomentOfInertia, const programOptions &options, const gsl_rng *r){
    /*
    Function to assign pulsar parameters

    Parameters
    ------------
    stellar_type : int
        Stellar type as defined by Hurley et al 2000
    Mass : double
        Stellar mass in solar masses
    Radius : double
        Stellar radius in solar radii
    pulsarSpinPeriod : double
        Pulsar spin period in ms
    pulsarMagneticField : double
        Pulsar magnetic field strength in G
    neutronStarMomentOfInertia : double
        Moment of inertia of neutron star in kg m^2
    options : programOptions
        User specified program options
    r : gsl_rng
        Random number generator
    */
    if(*stellar_type == NEUTRON_STAR){
        *pulsarMagneticField = pow(10.0, pulsarBirthMagneticFieldDistribution(options, r)) * gauss_to_tesla ; //Magnetic field in Gauss -> convert to Tesla
        *pulsarSpinPeriod = pulsarBirthSpinPeriodDistribution(options, r); // Spin period in ms
        *pulsarSpinFrequency = 2.0 * pi / (*pulsarSpinPeriod * ms);
	*neutronStarMomentOfInertia = momentOfInertia(*Mass, *Radius * RsolToKm) ; // in CGS g cm^2
        double alpha = 1.0;

	// Note we convert neutronStarMomentOfInertia from CGS to SI here
        *pulsarSpinDownRate = calculateSpinDownRate(*pulsarSpinFrequency, *neutronStarMomentOfInertia * g_to_kg * pow (cm_to_m,2.0), *pulsarMagneticField, *Radius * Rsol, alpha);
 
       // std::cout << "Pulsar spin down rate at birth = " << *pulsarSpinDownRate << std::endl;
       //	std::cout << "NS moment of inertia = " << *neutronStarMomentOfInertia << " g cm^2" << std::endl;
	
        *angularMomentum = 2.0 * pi * (*neutronStarMomentOfInertia) / (*pulsarSpinPeriod * ms) * g_to_kg * pow(cm_to_m,2.0); // in kg m^2 sec^-1

        //std::cout << "angular momentum assigned here: " << *angularMomentum << std::endl;
    }
}
