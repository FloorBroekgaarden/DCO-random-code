#include "lambdaLoveridge.h"

double lambdaLoveridgeMassCut(double metallicity, double Mass, unsigned long randomSeed){
	/*
	Function for determining the low/high mass cut according to metallicity required for calculating 
	common envelope lambda according to Loveridge et al 2011 ()

	Electronic tables, program and further information in: http://astro.ru.nl/~sluys/index.php?title=BE

	This function used to be called whatIsTheMassCut
	JIM BARRETT - 28/11/2016 - Passing through the random seed for error tracking purposes

	Parameters
	-----------
	metallicity : double
		Stars metallicity
	Mass : double
		Stars mass in Msol
	randomSeed : unsigned long
		This stars random seed

	Returns
	--------
	cutOffMass :double
		Cut off mass in Msol
	*/
    bool    debugging = false;
    double  massCutOff = 0.0;

    // SIMON: In general the metallicity won't match any of these, so need to use ranges.
    if(metallicity == metallicityGridLoveridge[0]){
        massCutOff = LMHMcutMassz000010;
    }
    else if(metallicity == metallicityGridLoveridge[1]){
        massCutOff = LMHMcutMassz000100;
    }
    else if(metallicity == metallicityGridLoveridge[2]){
        massCutOff = LMHMcutMassz001000;
    }
    else if(metallicity == metallicityGridLoveridge[3]){
        massCutOff = LMHMcutMassz001500;
    }
    else if(metallicity == metallicityGridLoveridge[4]){
        massCutOff = LMHMcutMassz002000;
    }
    else if(metallicity == metallicityGridLoveridge[5]){
        massCutOff = LMHMcutMassz003000;
    }
    else{
        std::cerr << randomSeed << TAB << "Non valid metallicity type for binding energy fits. Shouldn't get here." << std::endl;
    }

    return  massCutOff;

}

// JIM BARRETT - 28/11/2016 - Passing through the random seed for error tracking purposes
double  calculateLogRbound(double metallicity, double Mass, unsigned long randomSeed){
	/*
	Function required for calculating common envelope lambda according to Loveridge et al 2011 ()

	Function for determining the Radius boundary on low mass stars 

	Electronic tables, program and further information in: http://astro.ru.nl/~sluys/index.php?title=BE

	Parameters
	-----------
	metallicity : double
		Stars metallicity
	Mass : double
		Stars mass in Msol
	randomSeed : unsigned long
		This stars random seed
	
	Returns
	---------
	Rbound : double
		Boundary radius in Rsol
	*/

    bool    debugging = false;
    double  deltaM  = 1E-5;
    double  LogRbound  = 0.0;

    // SIMON: Again, we will need to make this into a range of metallicities
    if(metallicity == metallicityGridLoveridge[0]){
        for(int i=0; i < ncolsaCoefficients; i++){
            LogRbound += aCoefficientsz000010[i]*pow(log10(Mass+deltaM),i);
        }
    }
    else if(metallicity == metallicityGridLoveridge[1]){
        for(int i=0; i < ncolsaCoefficients; i++){
            LogRbound += aCoefficientsz000100[i]*pow(log10(Mass+deltaM),i);
        }
    }
    else if(metallicity == metallicityGridLoveridge[2]){
        for(int i=0; i < ncolsaCoefficients; i++){
            LogRbound += aCoefficientsz001000[i]*pow(log10(Mass+deltaM),i);
        }
    }
    else if(metallicity == metallicityGridLoveridge[3]){
        for(int i=0; i < ncolsaCoefficients; i++){
            LogRbound += aCoefficientsz001500[i]*pow(log10(Mass+deltaM),i);
        }
    }
    else if(metallicity == metallicityGridLoveridge[4]){
        for(int i=0; i < ncolsaCoefficients; i++){
            LogRbound += aCoefficientsz002000[i]*pow(log10(Mass+deltaM),i);
            if(debugging){
                std::cout << "a, M, logM ,LogRbound: " << aCoefficientsz002000[i] << "\t" << Mass+deltaM << "\t" << log10(Mass+deltaM) << "\t" << LogRbound << std::endl;
            }
        }
    }
    else if(metallicity == metallicityGridLoveridge[5]){
        for(int i=0; i < ncolsaCoefficients; i++){
            LogRbound += aCoefficientsz003000[i]*pow(log10(Mass+deltaM),i);
        }
    }
    else{
        std::cerr << randomSeed << TAB << "Non valid metallicity type for calcualting Rbound. Shouldn't get here." << std::endl;
    }

    if(debugging){
        std::cout << "Rbound, LogRbound: " << pow(10.0, LogRbound) << "\t" << LogRbound << std::endl;
    }

    return LogRbound;
}

double  massGroupAndEvolutionaryStage(double metallicity, double MassZAMS, double Mass, double COcoreMass, double Radius, unsigned long randomSeed){
	/*
	Function for determining if the group and evolutionary stage of the star

	Electronic tables, program and further information in: http://astro.ru.nl/~sluys/index.php?title=BE

	JIM BARRETT - 28/11/2016 - Passing through the random seed for error tracking purposes
	
	Parameters
	-----------
	metallicity : double
		Stars metallicity
	MassZAMS : double
		Stars zero age main sequence mass in Msol
	Mass : double
		Stars mass in Msol
	COcoreMass : double
		Stars carbon oxygen core mass in Msol
	Radius : double
		Stars radius in solar radii
	randomSeed : unsigned long
		This stars random seed

	Returns
	---------
	evolutionaryStage : double
		Evolutionary stage
	*/

    bool    debugging = false;
    double  cutMass = lambdaLoveridgeMassCut(metallicity, Mass,randomSeed);
    double  evolutionaryStage;
    double  GB;
    double  LogRbound = -7.3;
    const int RGB = 1;
    const int AGB = 2;

    if(COcoreMass > 0.0){
    	// AGB star if CO core exists, as stated in Loveridge code. To be checked.
        GB = AGB;
    }
    else{
        GB = RGB;
    }

    if(debugging){
        std::cout << "Mass, cutMass: " << Mass << "\t" << cutMass << std::endl;
    }
    if(Mass <= cutMass){

        if(GB == RGB){
            if(debugging){
                std::cout << "RGB loop" << std::endl;
            }
            evolutionaryStage = LMR1;

            LogRbound = calculateLogRbound(metallicity, Mass, randomSeed);

            if(log10(Radius) > LogRbound)
                evolutionaryStage = LMR2;

        }
        else{
            if(debugging){
                std::cout << "LMA loop" << std::endl;
            }
            evolutionaryStage = LMA;
        }

    }
    else{
        evolutionaryStage = HM;
    }

    if(debugging){
        std::cout << "Metal: " << metallicity << std::endl;
        std::cout << "MZAMS: " << MassZAMS << std::endl;
        std::cout << "Mass: " << Mass << std::endl;
        std::cout << "COcoreMass: " << COcoreMass << std::endl;
        std::cout << "Radius: " << Radius << std::endl;
        std::cout << "cutMass: " << cutMass << std::endl;
        std::cout << "GB: " << GB << std::endl;
        std::cout << "LogRadius: " << log10(Radius) << std::endl;
        std::cout << "LogRbound: " << LogRbound<< std::endl;
        std::cout << "EvStage: " << evolutionaryStage<< std::endl;
    }

    return evolutionaryStage;
}


double  calculateLogBindingEnergyLoveridge(double metallicity, double MassZAMS, double Mass, double COcoreMass, double Radius, bool isMassLoss, unsigned long randomSeed){
	/*
	Function for calculating the binding energy of the envelope according to detailed calculations done by Loveridge et al. 2011

	This function computes log[BE/erg] as a function of log[Z], Mzams, M, log[R/Ro] and GB.

	Electronic tables, program and further information in: http://astro.ru.nl/~sluys/index.php?title=BE

	Parameters
	-----------
	metallicity : double
		Stars metallicity
	MassZAMS : double
		Stars zero age main sequence mass in Msol
	Mass : double
		Stars mass in Msol
	COcoreMass : double
		Stars carbon oxygen core mass in Msol
	Radius : double
		Stars radius in solar radii
	isMassLoss : bool
		Whether mass loss is active
	randomSeed : unsigned long
		This stars random seed
	
	Returns
	--------
	logBE : double
		log binding energy in ergs
	*/

    bool    debugging = false;

    // Calculate closest metallicity on the grid
    int metallicityIndex = 0;
    double  deltaR = 1E-5;
    double  testingMetallicity = std::abs(metallicity - metallicityGridLoveridge[metallicityIndex]);
    double  logBE0 = 33.29866;

    for(int i = 1; i < ncolsMetallicities; i++){

        if(std::abs(metallicity - metallicityGridLoveridge[i])<testingMetallicity){
            testingMetallicity = std::abs(metallicity - metallicityGridLoveridge[i]);
            metallicityIndex   = i;
        }

    }

    // Determine the evolutionary stage of the star
    if(debugging){
        std::cout << "Z, MassZAMS, MASS, COcoreMass, Radius: " << metallicityGridLoveridge[metallicityIndex] << "\t" << MassZAMS << "\t" << Mass << "\t" << COcoreMass << "\t" << Radius << std::endl;
    }
    double evolutionaryStage = massGroupAndEvolutionaryStage(metallicityGridLoveridge[metallicityIndex], MassZAMS, Mass, COcoreMass, Radius, randomSeed);

    if(debugging){
        std::cout << "evolStage: " << evolutionaryStage << std::endl;
    }

    // Calcualte the 10-logarithm of the binding energy
    double  dLogBindingEnergy = 0.0;
    double  calculate_logBindingEnergy = 0.0;

    // zXXXXXXCoefficients[nrowszXXXXXX][ncols] = {GB, m, r, alpha}

    int l=0;

    if(debugging){
        std::cout << "j\tm\tr\talpha\tdLogBE\tlogBE\n";
    }

    if(metallicityIndex == 0){
        if(debugging){
            std::cout << "If 0)" << std::endl;
        }
        for(int j=0; j < nrowsz000010; j++){
            if(z000010Coefficients[j][0]==evolutionaryStage){
                dLogBindingEnergy = (z000010Coefficients[j][3])*(pow(log10(Mass),z000010Coefficients[j][1]))*(pow(log10(Radius+deltaR),z000010Coefficients[j][2]));
                calculate_logBindingEnergy += dLogBindingEnergy;
                l++;

                if(debugging){
                    std::cout << j << "\t" << z000010Coefficients[j][1] << "\t" << z000010Coefficients[j][2] << "\t" << z000010Coefficients[j][3] << "\t" << dLogBindingEnergy << "\t" << calculate_logBindingEnergy << std::endl;
                }
            }
        }
    }
    else if(metallicityIndex == 1){
        if(debugging){
            std::cout << "If 1)" << std::endl;
        }
        for(int j=0; j < nrowsz000100; j++){
            if(z000100Coefficients[j][0]==evolutionaryStage){
                dLogBindingEnergy = (z000100Coefficients[j][3])*(pow(log10(Mass),z000100Coefficients[j][1]))*(pow(log10(Radius+deltaR),z000100Coefficients[j][2]));
                calculate_logBindingEnergy += dLogBindingEnergy;
                l++;

                if(debugging){
                    std::cout << j << "\t" << z000100Coefficients[j][1] << "\t" << z000100Coefficients[j][2] << "\t" << z000100Coefficients[j][3] << "\t" << dLogBindingEnergy << "\t" << calculate_logBindingEnergy << std::endl;
                }
            }
        }
    }
    else if(metallicityIndex == 2){
        if(debugging){
            std::cout << "If 2)" << std::endl;
        }
        for(int j=0; j < nrowsz001000; j++){
            if(z001000Coefficients[j][0]==evolutionaryStage){
                dLogBindingEnergy = (z001000Coefficients[j][3])*(pow(log10(Mass),z001000Coefficients[j][1]))*(pow(log10(Radius+deltaR),z001000Coefficients[j][2]));
                calculate_logBindingEnergy += dLogBindingEnergy;
                l++;

                if(debugging){
                    std::cout << j << "\t" << z001000Coefficients[j][1] << "\t" << z001000Coefficients[j][2] << "\t" << z001000Coefficients[j][3] << "\t" << dLogBindingEnergy << "\t" << calculate_logBindingEnergy << std::endl;
                }

            }
        }

    }
    else if(metallicityIndex == 3){
        if(debugging){
            std::cout << "If 3)" << std::endl;
        }
        for(int j=0; j < nrowsz001500; j++){
            if(z001500Coefficients[j][0]==evolutionaryStage){
                dLogBindingEnergy = (z001500Coefficients[j][3])*(pow(log10(Mass),z001500Coefficients[j][1]))*(pow(log10(Radius+deltaR),z001500Coefficients[j][2]));
                calculate_logBindingEnergy += dLogBindingEnergy;
                l++;

                if(debugging){
                    std::cout << j << "\t" << z001500Coefficients[j][1] << "\t" << z001500Coefficients[j][2] << "\t" << z001500Coefficients[j][3] << "\t" << dLogBindingEnergy << "\t" << calculate_logBindingEnergy << std::endl;
                }
            }
        }
    }
    else if(metallicityIndex == 4){
        if(debugging){
            std::cout << "If 4)" << std::endl;
        }
        for(int j=0; j < nrowsz002000; j++){
            if(z002000Coefficients[j][0]==evolutionaryStage){
                dLogBindingEnergy = (z002000Coefficients[j][3])*(pow(log10(Mass),z002000Coefficients[j][1]))*(pow(log10(Radius+deltaR),z002000Coefficients[j][2]));
                calculate_logBindingEnergy += dLogBindingEnergy;
                l++;

                if(debugging){
                    std::cout << j << "\t" << z002000Coefficients[j][1] << "\t" << z002000Coefficients[j][2] << "\t" << z002000Coefficients[j][3] << "\t" << dLogBindingEnergy << "\t" << calculate_logBindingEnergy << std::endl;
                }

            }
        }
    }
    else if(metallicityIndex == 5){
        if(debugging) {
            std::cout << "If 5)" << std::endl;
        }
        for(int j=0; j < nrowsz003000; j++){
            if(z003000Coefficients[j][0]==evolutionaryStage){
                dLogBindingEnergy = (z003000Coefficients[j][3])*(pow(log10(Mass),z003000Coefficients[j][1]))*(pow(log10(Radius+deltaR),z003000Coefficients[j][2]));
                calculate_logBindingEnergy += dLogBindingEnergy;
                l++;

                if(debugging){
                    std::cout << j << "\t" << z003000Coefficients[j][1] << "\t" << z003000Coefficients[j][2] << "\t" << z003000Coefficients[j][3] << "\t" << dLogBindingEnergy << "\t" << calculate_logBindingEnergy << std::endl;
                }
            }
        }

    }
    else{
        std::cerr << randomSeed << TAB << "Wrong metallicity index. Shouldn't get here." << std::endl;
    }


    /*
     ! Compute the 10-logarithm of the binding energy:
     calc_logBE = 0.0_dbl
     do ii=1,ndat(iz,ig)
        dlogBE = alphas(iz,ig,ii) * (logM+tiny(logM))**dble(ms(iz,ig,ii)) * (logR+tiny(logR))**dble(rs(iz,ig,ii))  ! Avoid 0**0
            calc_logBE = calc_logBE + dlogBE
     end do


           */
    // Compute and apply the mass-loss correction factor Lambda
    double  Lambda = 1.0;

    if(isMassLoss){
        Lambda = 1.0 + 0.25*pow((MassZAMS-Mass)/MassZAMS,2.0);
    }


    calculate_logBindingEnergy *= Lambda;
    calculate_logBindingEnergy += logBE0;
//    static const int    nrowsz000010        = 469;
 //   static const double z000010Coefficients[nrowsz000010][ncols];

//                        calc_logBE = lambda * calc_logBE
//
//
//                        ! BE was originally expressed in erg/solar mass to avoid large numbers, so we need to convert to erg here:
//                        calc_logBE = calc_logBE + logBE0

    if(debugging){
        for(int i = 0; i < ncolsMetallicities; i++){
            std::cout << "Z" << i << ": " << metallicityGridLoveridge[i] << std::endl;

        }
        std::cout << "Z: " << metallicity << std::endl;
        std::cout << "Zclose: " << metallicityGridLoveridge[metallicityIndex] << std::endl;
        std::cout << "evolutionaryStage: " << evolutionaryStage << std::endl;
        std::cout << "Lambda: " << Lambda << std::endl;
        std::cout << "l: " << l << std::endl;

    }

//    std::cout <<metallicity << "\t" << metallicityGridLoveridge[metallicityIndex] << "\t" << MassZAMS << "\t\t" << Mass << "\t" << Radius << "\t" << evolutionaryStage << "\t" << calculate_logBindingEnergy << std::endl;
    if(evolutionaryStage==0){
        evolutionaryStage=1;
    }

    if(debugging){
        std::cout <<metallicity << "\t" << metallicityGridLoveridge[metallicityIndex] << "\t" << MassZAMS << "\t\t" << Mass << "\t" << Radius << "\t" << evolutionaryStage << "\t" << calculate_logBindingEnergy << std::endl;
    }
//       std::cout <<metallicity << "\t" << metallicityGridLoveridge[metallicityIndex] << "\t" << MassZAMS << "\t\t" << Mass << "\t" << Radius << "\t" << evolutionaryStage << "\t" << calculate_logBindingEnergy << std::endl;
	return	 calculate_logBindingEnergy;
}


double  calculateLogRecombinationEnergyLoveridge(double metallicity, double Mass, double Radius, unsigned long randomSeed){
    // Function for calculating the recombination energy of the envelope according to detailed calculations done by Loveridge et al. 2011
    // This function computes log[BE_recom/erg] as a function of log[Z], Mzams, M, log[R/Ro] and GB.
    // Electronic tables, program and further information in: http://astro.ru.nl/~sluys/index.php?title=BE

    bool    debugging = false;
    // debugging = true;

    // Calculate closest metallicity on the grid
    int metallicityIndex = 0;
    double  deltaR = 1E-5;
    double  testingMetallicity = std::abs(metallicity - metallicityGridLoveridge[metallicityIndex]);
    double  logBE0 = 33.29866;
    double  evolutionaryStage = Recom;

    for(int i = 1; i < ncolsMetallicities; i++){

        if(std::abs(metallicity - metallicityGridLoveridge[i])<testingMetallicity){
            testingMetallicity = std::abs(metallicity - metallicityGridLoveridge[i]);
            metallicityIndex   = i;
        }

    }

    // Determine the evolutionary stage of the star
    if(debugging){
        std::cout << "Z, MASS, Radius: " << metallicityGridLoveridge[metallicityIndex] << "\t" << Mass << "\t" << Radius << std::endl;
    }


    // Calcualte the 10-logarithm of the binding energy
    double  dLogRecombinationEnergy = 0.0;
    double  calculate_logRecombinationEnergy = 0.0;

    int l=0;

    if(debugging){
        std::cout << "j\tm\tr\talpha\tdLogBE\tlogBE\n";
    }

    if(metallicityIndex == 0){
        if(debugging){
            std::cout << "If 0)" << std::endl;
        }
        for(int j=0; j < nrowsz000010; j++){
            if(z000010Coefficients[j][0]==evolutionaryStage){
                dLogRecombinationEnergy = (z000010Coefficients[j][3])*(pow(log10(Mass),z000010Coefficients[j][1]))*(pow(log10(Radius+deltaR),z000010Coefficients[j][2]));
                calculate_logRecombinationEnergy += dLogRecombinationEnergy;
                l++;

                if(debugging){
                    std::cout << j << "\t" << z000010Coefficients[j][1] << "\t" << z000010Coefficients[j][2] << "\t" << z000010Coefficients[j][3] << "\t" << dLogRecombinationEnergy << "\t" << calculate_logRecombinationEnergy << std::endl;
                }
            }
        }
    }
    else if(metallicityIndex == 1){
        if(debugging){
            std::cout << "If 1)" << std::endl;
        }
        for(int j=0; j < nrowsz000100; j++){
            if(z000100Coefficients[j][0]==evolutionaryStage){
                dLogRecombinationEnergy = (z000100Coefficients[j][3])*(pow(log10(Mass),z000100Coefficients[j][1]))*(pow(log10(Radius+deltaR),z000100Coefficients[j][2]));
calculate_logRecombinationEnergy += dLogRecombinationEnergy;
                l++;

                if(debugging){
                    std::cout << j << "\t" << z000100Coefficients[j][1] << "\t" << z000100Coefficients[j][2] << "\t" << z000100Coefficients[j][3] << "\t" << dLogRecombinationEnergy << "\t" << calculate_logRecombinationEnergy << std::endl;
                }
            }
        }
    }
    else if(metallicityIndex == 2){
        if(debugging){
            std::cout << "If 2)" << std::endl;
        }
        for(int j=0; j < nrowsz001000; j++){
            if(z001000Coefficients[j][0]==evolutionaryStage){
                dLogRecombinationEnergy = (z001000Coefficients[j][3])*(pow(log10(Mass),z001000Coefficients[j][1]))*(pow(log10(Radius+deltaR),z001000Coefficients[j][2]));
calculate_logRecombinationEnergy += dLogRecombinationEnergy;
                l++;

                if(debugging){
                    std::cout << j << "\t" << z001000Coefficients[j][1] << "\t" << z001000Coefficients[j][2] << "\t" << z001000Coefficients[j][3] << "\t" << dLogRecombinationEnergy << "\t" << calculate_logRecombinationEnergy << std::endl;
                }

            }
        }

    }
    else if(metallicityIndex == 3){
        if(debugging){
            std::cout << "If 3)" << std::endl;
        }
        for(int j=0; j < nrowsz001500; j++){
            if(z001500Coefficients[j][0]==evolutionaryStage){
                dLogRecombinationEnergy = (z001500Coefficients[j][3])*(pow(log10(Mass),z001500Coefficients[j][1]))*(pow(log10(Radius+deltaR),z001500Coefficients[j][2]));
calculate_logRecombinationEnergy += dLogRecombinationEnergy;
                l++;

                if(debugging){
                    std::cout << j << "\t" << z001500Coefficients[j][1] << "\t" << z001500Coefficients[j][2] << "\t" << z001500Coefficients[j][3] << "\t" << dLogRecombinationEnergy << "\t" << calculate_logRecombinationEnergy << std::endl;
                }
            }
        }
    }
    else if(metallicityIndex == 4){
        if(debugging){
            std::cout << "If 4)" << std::endl;

        }
        for(int j=0; j < nrowsz002000; j++){
            if(z002000Coefficients[j][0]==evolutionaryStage){
                dLogRecombinationEnergy = (z002000Coefficients[j][3])*(pow(log10(Mass),z002000Coefficients[j][1]))*(pow(log10(Radius+deltaR),z002000Coefficients[j][2]));
calculate_logRecombinationEnergy += dLogRecombinationEnergy;
                l++;

                if(debugging){
                    std::cout << j << "\t" << z002000Coefficients[j][1] << "\t" << z002000Coefficients[j][2] << "\t" << z002000Coefficients[j][3] << "\t" << dLogRecombinationEnergy << "\t" << calculate_logRecombinationEnergy << std::endl;
                }

            }
        }
    }
    else if(metallicityIndex == 5){
        if(debugging) {
            std::cout << "If 5)" << std::endl;
        }
        for(int j=0; j < nrowsz003000; j++){
            if(z003000Coefficients[j][0]==evolutionaryStage){
                dLogRecombinationEnergy = (z003000Coefficients[j][3])*(pow(log10(Mass),z003000Coefficients[j][1]))*(pow(log10(Radius+deltaR),z003000Coefficients[j][2]));
calculate_logRecombinationEnergy += dLogRecombinationEnergy;
                l++;

                if(debugging){
                    std::cout << j << "\t" << z003000Coefficients[j][1] << "\t" << z003000Coefficients[j][2] << "\t" << z003000Coefficients[j][3] << "\t" << dLogRecombinationEnergy << "\t" << calculate_logRecombinationEnergy << std::endl;
                }
            }
        }

    }
    else{
        std::cerr << randomSeed << TAB << "Wrong metallicity index. Shouldn't get here." << std::endl;
    }


    /*
     ! Compute the 10-logarithm of the binding energy:
     calc_logBE = 0.0_dbl
     do ii=1,ndat(iz,ig)
     dlogBE = alphas(iz,ig,ii) * (logM+tiny(logM))**dble(ms(iz,ig,ii)) * (logR+tiny(logR))**dble(rs(iz,ig,ii))  ! Avoid 0**0
     calc_logBE = calc_logBE + dlogBE
     end do


     */


    calculate_logRecombinationEnergy += logBE0;

    if(debugging){
        for(int i = 0; i < ncolsMetallicities; i++){
            std::cout << "Z" << i << ": " << metallicityGridLoveridge[i] << std::endl;
        }
        std::cout << "Z: " << metallicity << std::endl;
        std::cout << "Zclose: " << metallicityGridLoveridge[metallicityIndex] << std::endl;
        std::cout << "l: " << l << std::endl;

    }

//    std::cout <<metallicity << "\t" << metallicityGridLoveridge[metallicityIndex] << "\t\t" << Mass << "\t" << Radius << "\t" << calculate_logRecombinationEnergy << std::endl;
    if(debugging){
        std::cout << "\t\t\t\t\t\t\t\t" << calculate_logRecombinationEnergy << std::endl;
    }

    return calculate_logRecombinationEnergy;
}

double lambdaLoveridgeEnergyFormalism(double metallicity, double MassZAMS, double Mass, double envelopeMass, double COcoreMass, double Radius, bool isMassLoss, unsigned long randomSeed){
	/*
	Function for calculating lambda parameter from the so-called energy formalism of CE (Webbink 1984).

	Binding energy from detailed models (Loveridge et al. 2011) is given in [E]=ergs, so use cgs

	JIM BARRETT - 28/11/2016 - Passing through the random seed for error tracking purposes

	Parameters
	-----------
	metallicity : double
		Stars metallicity
	MassZAMS : double
		Stars zero age main sequence mass in Msol
	Mass : double
		Stars mass in Msol
	envelopeMass : double
		Stars envelope mass in Msol
	COcoreMass : double
		Stars carbon oxygen core mass in Msol
	Radius : double
		Stars radius in solar radii
	isMassLoss : bool
		Whether mass loss is active
	randomSeed : unsigned long
		This stars random seed

	Returns
	--------
	lambdaLoveridge : double
		Common envelope lambda parameter
	*/

    double  lambda = pow(10.0,-20);	// Avoid 0
    bool    debugging = false;

	if(debugging){
		std::cout << "Z: " << metallicity << std::endl;
		std::cout << "MassZAMS: " << MassZAMS << std::endl;
		std::cout << "Mass: " << Mass << std::endl;
		std::cout << "envelopeMass: " << envelopeMass << std::endl;
		std::cout << "COcoremass: " << COcoreMass << std::endl;
		std::cout << "Radius: " << Radius << std::endl;
	}


    double	logBindingEnergy = calculateLogBindingEnergyLoveridge(metallicity, MassZAMS, Mass, COcoreMass, Radius, false, randomSeed);
    double  bindingEnergy = pow(10.0,logBindingEnergy);


    if(bindingEnergy > 0.0){
        lambda = Gcgs*(Mass*MsolTog)*(envelopeMass*MsolTog)/((Radius*RsolToAU*AUTocm)*bindingEnergy);

		if(debugging){
            std::cout << "Binding energy > 0. Calculating binding energy." << std::endl;
        }
    }

    if(debugging){
		std::cout << "logBindingEnergy: " << logBindingEnergy << "ergs" << std::endl;
        std::cout << "bindingEnergy: " << bindingEnergy << "ergs" << std::endl;
        std::cout << "Gcgs: " << Gcgs << std::endl;
        std::cout << "Mass: " << Mass << std::endl;
        std::cout << "Mass [cgs]: " << Mass*MsolTog << std::endl;
        std::cout << "envelopeMass: " << envelopeMass << std::endl;
        std::cout << "envelopeMass [cgs]: " << envelopeMass*MsolTog << std::endl;
        std::cout << "lambda: " << lambda << std::endl;
    }

    return lambda;
}

