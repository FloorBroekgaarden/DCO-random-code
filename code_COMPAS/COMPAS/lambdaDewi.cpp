#include "lambdaDewi.h"

double	lambdaDewi(double Mass, double coreMass, double Radius, double RadiusZAMS, double Luminosity, double stellarType, unsigned long randomSeed){
// 	ALEJANDRO - 17/05/2017 - Not fully tested, nor fully coded, nor fully trusted. Any other \lambda prescription is personally preferred. 
//	Fit from Appendix A in Claeys+2014, based on the method described in Dewi & Tauris 2000
//	arXiv:1401.2895 for Claeys+2014
//	arXiv:0007034 for Dewi and Tauris 2000
//	[Radius] = [RadiusZAMS] = Rsol
//	[Mass] = [envMass] = Msol
// 	[Luminosity] = Lsol

	
	double	lambda1 = NEVER_SET;
	double	lambda2 = NEVER_SET;
	double	lambda3 = NEVER_SET;
	double	lambda4 = NEVER_SET;
	double	lambda5 = NEVER_SET;
	double	lambdaCE = NEVER_SET;
	double	massEnv = NEVER_SET;
	bool	debugging = false;
	
	if((Mass>coreMass)and(coreMass>0.0)){
		massEnv = Mass-coreMass;
	}
	else{
		massEnv = 0.0;
	}
	
	if(debugging){
		std::cout << "Calculatig lambda Dewi" << std::endl;
		std::cout << "Mass [Msol]:\t" << Mass << std::endl;
		std::cout << "massEnv [Msol]:\t" << massEnv << std::endl;
		std::cout << "coreMass [Msol]:\t" << coreMass << std::endl;		
		std::cout << "Radius [Rsol]:\t" << Radius << std::endl;
		std::cout << "RadiusZAMS [Rsol]:\t" << RadiusZAMS << std::endl;
		std::cout << "Luminosity [Lsol]:\t" << Luminosity << std::endl;
		std::cout << "StellarType:\t" << stellarType << std::endl;
	}
	
	// (A.2) Claeys+2014
	lambda2 = 0.42*pow(RadiusZAMS/Radius,0.4);
	
	if((stellarType == HERTZSPRUNG_GAP)or(stellarType == FIRST_GIANT_BRANCH)){
		if(debugging){std::cout << "StellarType HG or GB" << std::endl;}
		double	firstTerm =	3.0/(2.4+pow(Mass,-3.0/2.0));
		double	secondTerm = -0.15*log10(Luminosity);

		// (A.3) Claeys+2014
		lambda1 = std::min(0.80,firstTerm+secondTerm);
	}
	else if((stellarType == CORE_HELIUM_BURNING)or(stellarType == EARLY_ASYMPTOTIC_GIANT_BRANCH)or(stellarType == THERMALLY_PULSING_ASYMPTOTIC_GIANT_BRANCH)){
		if(debugging){std::cout << "StellarType CHeB, EAGB or TPAGB" << std::endl;}
		double	firstTerm = -0.9;
		double 	secondTerm = 0.58+0.75*log10(Mass);
		double	thirdTerm = - 0.08*log10(Luminosity);
		// (A.4) Claeys+2014
		lambda3 = std::min(firstTerm,secondTerm)+thirdTerm;
		
		if((stellarType == CORE_HELIUM_BURNING)or(stellarType == EARLY_ASYMPTOTIC_GIANT_BRANCH)){
			double	fourthTerm = 0.8;
			double	fifthTerm = 1.25 - 0.15*log10(Luminosity);
			double	sixthTerm = std::min(fourthTerm,fifthTerm);
			// (A.5) Top, Claeys+2014
			lambda1 = std::min(sixthTerm,lambda3);
		}
		else{
			double seventhTerm = -3.5 - 0.75*log10(Mass) + log10(Luminosity);
			// (A.5) Bottom, Claeys+2014
			double	tempLambda1 = std::max(seventhTerm,lambda3);
			double	TPAGBcap = 1.0;
			lambda1 = std::max(tempLambda1, TPAGBcap);
		}
	}
	else{
		// Do nothing
	}
	
	
	if((stellarType >= HERTZSPRUNG_GAP)and(stellarType <= THERMALLY_PULSING_ASYMPTOTIC_GIANT_BRANCH)){
		if(massEnv == 0.0){
			// (A.1) Top, Claeys+2014	
			lambdaCE = 2.0*lambda2;			
		}
		else if((massEnv > 0.0)and(massEnv < 1.0)){
			double	eighthTerm = std::pow(massEnv,0.5)*(lambda1 - lambda2);
			// (A.1) Mid, Claeys+2014
			lambdaCE = 2.0*(lambda2+eighthTerm);
		}
		else if(massEnv >= 1.0 ){
			// (A.1) Bottom, Claeys+2014
			lambdaCE = 2.0*lambda1;			
		}
		else{
			std::cerr << randomSeed << "\tInside lambdaDewi. Shouldn't get here." << std::endl;
		}
		
	}
	else if((stellarType == NAKED_HELIUM_STAR_MS)or(stellarType == NAKED_HELIUM_STAR_HETZSPRUNG_GAP)or(stellarType == NAKED_HELIUM_STAR_GIANT_BRANCH)){
		lambdaCE = 0.5;
	}
	else{
		// Do nothing
	}

	if(debugging){
		std::cout << "lambda1:\t" << lambda1 << std::endl;
		std::cout << "lambda2:\t" << lambda2 << std::endl;
		std::cout << "lambda3:\t" << lambda3 << std::endl;
		std::cout << "lambda4:\t" << lambda4 << std::endl;
		std::cout << "lambda5:\t" << lambda5 << std::endl;
		std::cout << "lambdaCE:\t" << lambdaCE << std::endl << std::endl;		
	}
	
	// Missing (A.6),(A.7),(A.8),(A.9),(A.10),(A.11) and (A.12), which have to do with ionization energy

	return	lambdaCE;
}