//
//  massTransfer.cpp

#include "massTransfer.h"


double numericalZRocheLobe(double Md, double Ma, double a, double &fa, double jloss, double MTcoefficient, double Ra, double thermalRateAccretor, double thermalRateDonor, double stellarTypeAccretor, double stellarTypeDonor, programOptions optionscopy){
    // Numerical calculation of the Roche Lobe after mass transfer as in StarTrack. Described in Belczynski et al. 2008. Used for a regular star accretor and non-conservative Mass Transfer.
    // dJ=Beta*((1.0-Fa)*(Ma2-Ma1)/(Ma1+Mb1))*Jorb1;
    // a2=((Ma2+Mb2)*pow(Jorb1+dJ,2.0))/(GGG*Ma2*Ma2*Mb2*Mb2);  /* non-conservative MT assumption */
    // jloss = Beta: "Podsiadlowski et al. 1992 -specific angular momentum of matter [2Pia^2/P]"
    // fa: fraction of mass accreted. fa = 1, MT is conservative. fa = 0, MT is fully non-conservative.

    double  ZRL      = 0.0;
    bool    debugging = false;

    /*
     //Startrack prescription for compact objects. Hardcoded from StarTrack.
     if(Kb>=10) {                            // compact accretor, material lost with J of compact accretor //
     q=Ma/Mb;
     F=(0.4*q+(0.333333*pow(q,0.6666667))/(1.0+pow(q,0.333333)))/(0.6*q+pow(q,0.333333)*log(1.0+pow(q,0.333333)));

     if(Ka<10)                             // non-compact donor //
        beta=Fa;
     else if(dMmta_old<acc)
        beta=1.0;
     else
        beta=min(1.0,fabs(dMmtb_old)/dMmta_old);

     dlnR=(2.0*q*q-(1.0-beta)*q-2.0)/(1.0+q)+(1.0+beta*q)*(0.6666667-F);
     }
     */

    double Mdfinal  = (1.0 - MTcoefficient)*Md;                     // Final mass of donor
    double Mafinal  = Ma + fa*MTcoefficient*Md;                     // Final mass of accretor
    double J        = (Md*Ma)*sqrt(G1*(Md+Ma)*a)/(Md+Ma);
    double dJ       = jloss*((1.0-fa)*(Mdfinal-Md)/(Md+Ma))*J;
    double afinal   = (Mdfinal+Mafinal)*pow(J+dJ,2.0)/(G1*Mdfinal*Mdfinal*Mafinal*Mafinal);
    double RL       = a*rocheLobeRadius(Md, Ma);
    double RLfinal  = afinal*rocheLobeRadius(Mdfinal, Mafinal);

    ZRL =(log(RLfinal)-log(RL))/(log(Mdfinal)-log(Md));

//    debugging = true;
    if(debugging){
        std::cout << "fa: " << fa << std::endl;
        std::cout << "J: " << J << std::endl;
        std::cout << "dJ: " << dJ << std::endl;
        std::cout << "Md: " << Md << std::endl;
        std::cout << "Mdfinal: " << Mdfinal << std::endl;
        std::cout << "Ma: " << Ma << std::endl;
        std::cout << "Mafinal: " << Mafinal << std::endl;
        std::cout << "a: " << a << std::endl;
        std::cout << "afinal: " << afinal << std::endl;
        std::cout << "ZRL: " << ZRL << std::endl;
    }

    return ZRL;
}

double ZRocheLobe(double Md, double Ma, double fa){
    // Response of the donor's RL radius  as calculated in Woods et al. (2012). Parameter Beta defines how conservative mass transfer (MT) is. Beta = 1 for coservative MT and Beta = 1 for completely non-conservative MT. Formula from M. Sluys notes "Binary evolution in a nutshell".
    // 0 <= Beta <= 1

    double  q = Md/Ma;
    double  Beta = fa;

    double  K1 = ((2.0*pow(Md,2.0))-(2.0*pow(Ma,2.0))-(Md*Ma*(1.0-Beta)))/(Ma*(Md+Ma));
    double  K2 = (2.0/3.0)-((pow(q,1.0/3.0)*((1.2*pow(q,1.0/3.0))+(1.0/(1.0+pow(q,1.0/3.0)))))/(3.0*((0.6*pow(q,2.0/3.0))+(log(1+pow(q,1.0/3.0))))));
    double  K3 = 1.0+(Beta*Md/Ma);

    return K1+(K2*K3);
}

double massTransferRateEq(double Zdon, double Zlob, double Zevl, double Mass){
    // Calculation of mass transfer rate, eq. (45), Belczynski et al. 2008
    // tau_GR and tau_MB are very large, as in 'infinity'
    // [Zevl] = [dt^-1] = Myr^-1
    // [Mass] = Msol
    // Zdon and Zlob are dimensionless
    return -(Zevl*Mass)/(Zdon-Zlob);
}

double massTransferTimescaleEq(double Mdoteq, double Mass){
    // Calculation of mass transfer timescale, eq. (46), Belczynski et al. 2008
    return  -Mass/Mdoteq;
}

double massTransferRateThermal(double Mass, double tau){
    // Calculation of mass transfer thermal rate, eq. (48), Belczynski et al. 2008
    return  -Mass/tau;
}

double EddingtonCriticalRate(double Radius, double stellarTypeAccretor, double stellarTypeDonor){
    // Calculation of Eddington critical rate, eq. (30), Belczynski et al. 2008
    // Default case if accretor is a BH, NS or a WD
    // [Radius] = Rsol
    // [MDotEdd] = Msol Myr^-1
    double Racc = Radius;
    double epsilon = 1.0;
    double X = 0.0;
    double MDotEdd = 0.0;
    bool debbuging= false;
    
	/*
	// ALEJANDRO - 27/03/2017 - Belczsynski approximation of the Eddington accretion limit, as in Belczysnki+2008
	if (stellarTypeDonor <= THERMALLY_PULSING_ASYMPTOTIC_GIANT_BRANCH) {       // Donor is H-rich
        X=0.7;
    }
	if(stellarTypeAccretor == BLACK_HOLE){    // Case accretor is a BH
        Racc = 3*Radius;
        epsilon =0.5;
		MDotEdd = 2.08797e-03*Racc/(epsilon*(1.0+X))*MyearToyear; // Msol Myr^-1
    }
	*/
	
	// ALEJANDRO - 27/03/2017 - Sluys approximation of the Eddington accretion limit, as in "Binary evolution in a nusthell" notes
	if(stellarTypeAccretor == BLACK_HOLE){    // Case accretor is a BH
        Racc = 3*Radius*RsolToKm;
	}
	else if((stellarTypeAccretor <= NEUTRON_STAR)and((stellarTypeAccretor >= HELIUM_WHITE_DWARF))){
		Racc = Radius*RsolToKm;
	}
	else{
		std::cerr << "Error in EddingtonCriticalRate. Calculating the Eddington accretion limit for a non-compact accretor." << std::endl;
	}
    
	MDotEdd = (1.5E-8)*(Racc/10.0)*MyearToyear;
	
    if(debbuging){
    std::cout << "Radius : " << Radius << std::endl;
    std::cout << "Racc : " << Racc << std::endl;
    std::cout << "X : " << X << std::endl;
    std::cout << "Epsilon : " << epsilon << std::endl;
    std::cout << "MDotEdd : " << MDotEdd << std::endl;
    }

    return MDotEdd;
}

double calculateMassTransferOrbitDegenerateAccretor(double md, double ma, double fa, double MdDot, double a, double worb, double e, double dt, double J, double jloss, double Ra, double thermalRateAccretor, double thermalRateDonor, double stellarTypeAccretor, double stellarTypeDonor, const programOptions & options, unsigned long m_randomSeed){
    // Calculation of semi-major axis due to angular momentum loss, eq. (32), Belczynski et al. 2008
    // Solution for a degenerate accretor
	// [md]		= Msol 
	// [ma]		= Msol
	// [a]		= AU
	// [MdDot]	= Msol Myr^-1
	// [dt]		= Myr
	// [Ra]		=
	
	
    double 	Jorb = ((ma*md)/(ma+md))*sqrt(G1*(ma+md)*a);
    double 	Rcom = 0.0;
    double 	Jprime = 0.0;
    double 	mdPrime = 0.0;
    double 	maPrime = 0.0;
    double 	MPrime = 0.0;
	double	aprime = NEVER_SET;
	bool	debugging = false;
	// debugging = true;


	if(debugging){
		std::cout << "Jorb:\t" << Jorb << std::endl;
		std::cout << "a [AU]:\t" << a << std::endl;
		std::cout << "Md [Msol]:\t" << md << std::endl;
		std::cout << "Ma [Msol]:\t" << ma << std::endl;
		std::cout << "dt [Myrs]:\t" << dt << std::endl;
		std::cout << "fa:\t" << fa << std::endl;
		std::cout << "jloss:\t" << jloss << std::endl;
		std::cout << "MdDot:\t" << MdDot << std::endl;
		std::cout << "MdDot*dt:\t" << MdDot*dt << std::endl;
		std::cout << jloss*Jorb*(1.0-fa)/(ma+md) << std::endl;
	}
	
	// ALEJANDRO - 29/11/2016 - Next line was commented, as it didn't allow case A MT to solve properly for the orbit, specially in the fast phase. This should be revisited in order to make it efficient.
//    if(options.massTransferPrescription == BELCZYNSKI_MASS_TRANSFER or envelopeType(stellarTypeDonor, options) == RADIATIVE_ENVELOPE){
    if(options.massTransferPrescription == BELCZYNSKI_MASS_TRANSFER){
		if(debugging){
			std::cout << "calculateMassTransferOrbitDegenerateAccretor function. Belczynski MT." <<std::endl;
		}
        Rcom = a*md/(md+ma);         // distance from the compact object to the CM
        Jprime = pow(Rcom, 2.0)*worb*(1-fa)*MdDot*dt + Jorb;
        mdPrime = (MdDot*dt)+md;     // mass of the donor after mass transfer
        maPrime = (-fa*MdDot*dt)+ma; // mass of the accretor after mass transfer
        MPrime  = mdPrime + maPrime;


        return pow(Jprime,2.0)*MPrime/(G1*pow(mdPrime*maPrime,2.0)); // Change in the orbit due to angular momentum loss from MT
    }

    else if(options.massTransferPrescription == DEMINK_MASS_TRANSFER){

        int		niter = 1000;      // Number of itterations to solve for the orbit. Abitrarly fixed to 1000.
        double 	delta_t = dt/niter;
		
		if(debugging){
			std::cout << "calculateMassTransferOrbitDegenerateAccretor function. DeMink MT." <<std::endl;
            std::cout << "Jorb, a: " << Jorb << " " << a << std::endl;
		}

        for(int i=1; i<= niter ; i++){
            // Jprime = pow((a*md/(md+ma)),2.0)*(sqrt(G1*(md+ma)/pow(a,3.0)))*(1.0-fa)*MdDot*delta_t + Jorb;
            Jprime = (jloss*Jorb*(1.0-fa)*MdDot/(ma+md))*delta_t + Jorb;
            aprime = (((-2.0*(MdDot/md))*(1.0 - (fa*(md/ma)) - ((1.0-fa)*(jloss+0.5)*(md/(md+ma)))))*a*delta_t) + a;
			
            md = md+(MdDot*delta_t);
            ma = ma-(MdDot*delta_t*fa);
            jloss = gammaAngularMomentumLoss(md, ma, options);
            massAcceptanceRate(Ra, thermalRateAccretor, thermalRateDonor, stellarTypeAccretor, stellarTypeDonor, options, fa, m_randomSeed);
            Jorb = Jprime;
			a	= aprime;
			
           if(debugging){
                std::cout << "Jprime, md, ma, jloss, fa:\t" << Jprime << " " << md << " " << ma << " " << jloss << " " << fa << std::endl;
                std::cout << "aprime, aprimeVar:\t" << pow(Jprime/(md*ma),2.0)*(md+ma)/G1 << ", " << aprime << std::endl;
            }
        }

        mdPrime = md;
        maPrime = ma;
        MPrime = mdPrime + maPrime;

		aprime = pow(Jprime/(mdPrime*maPrime),2.0)*MPrime/G1;
		if(debugging){
			std::cout << "a'[]:\t" << aprime << std::endl;
		}
        return aprime; // Change in the orbit due to angular momentum loss from MT
    }
    else{
        std::cout << "Invalid or no MT prescription for compact accretor. Missing orbit calculation. Return a.\n";
        return a;
    }

}

double calculateMassTransferOrbitNonDegenerateAccretor(double md, double mdEnv, double ma, double fa, double MdDot, double a, double worb, double e, double dt, double tau, double J, double jloss, double Ra, double thermalRateAccretor, double thermalRateDonor, double stellarTypeAccretor, double stellarTypeDonor, const programOptions & options, unsigned long m_randomSeed){
    // Calculation of semi-major axis due to angular momentum loss, eq. (33), Belczynski et al. 2008
    // Solution for a non-degenerate accretor
    // jloss is the specific angular momentum
    // fa is the fraction of matter accreted: fa = 0 for fully conservative and fa = 1 for fully non-conservative
    double Jorb = ((ma*md)/(ma+md))*sqrt(G1*(ma+md)*a);
    double Jprime = 0.0;
    double mdPrime = 0.0;
    double maPrime = 0.0;
    double MPrime = 0.0;
    double aprime = 0.0;

    bool   debugging = false;
	//debugging = true;
	
    if(options.massTransferPrescription == BELCZYNSKI_MASS_TRANSFER or envelopeType(stellarTypeDonor, options, m_randomSeed) == RADIATIVE_ENVELOPE){

        Jprime = (jloss*Jorb*(1.0-fa)*MdDot/(ma+md))*dt + Jorb;
        mdPrime = (MdDot*dt)+md;     // mass of the donor after mass transfer
        maPrime = (-fa*MdDot*dt)+ma; // mass of the accretor after mass transfer
        MPrime  = mdPrime + maPrime;

        return pow(Jprime/(mdPrime*maPrime),2.0)*MPrime/G1; // Change in the orbit due to angular momentum loss from MT
    }
    else if(options.massTransferPrescription == DEMINK_MASS_TRANSFER){

        int     niter = 1000;      // Number of itterations to solve for the orbit. Abitrarly fixed to 1000.
        double delta_t = dt/niter;

        if(debugging){
            std::cout << "Jorb, a: " << Jorb << " " << a << std::endl;

        }
        for(int i=1; i<= niter ; i++){
            Jprime = (jloss*Jorb*(1.0-fa)*MdDot/(ma+md))*delta_t + Jorb;

            // Checking semi-major axis prime
            aprime = (((-2.0*(MdDot/md))*(1.0 - (fa*(md/ma)) - ((1.0-fa)*(jloss+0.5)*(md/(md+ma)))))*a*delta_t) + a;

            if(debugging){
                std::cout << "delta_t: " << delta_t << std::endl;
                std::cout << "MdDot: " << MdDot << std::endl;
                std::cout << "a*delta_t" << a*delta_t << std::endl;
                std::cout << "factor: " << (-2.0*(MdDot/md)*(1.0 - (fa*(md/ma)) - ((1.0-fa)*(jloss+0.5)*(md/(md+ma)))))*a*delta_t << std::endl;
                std::cout << "aprimeTest: " << aprime << std::endl;
            }
            //
            md = md+(MdDot*delta_t);
            ma = ma-(MdDot*delta_t*fa);
            jloss = gammaAngularMomentumLoss(md, ma, options);
            massAcceptanceRate(Ra, thermalRateAccretor, thermalRateDonor, stellarTypeAccretor, stellarTypeDonor, options, fa, m_randomSeed);
            Jorb = Jprime;
            a = aprime;

            if(debugging){
                std::cout << "Jprime, md, ma, jloss, fa: " << Jprime << " " << md << " " << ma << " " << jloss << " " << fa << std::endl;
                std::cout << "aprime: " << pow(Jprime/(md*ma),2.0)*(md+ma)/G1 << std::endl;
            }

        }

        mdPrime = md;
        maPrime = ma;
        MPrime = mdPrime + maPrime;

        return pow(Jprime/(mdPrime*maPrime),2.0)*MPrime/G1; // Change in the orbit due to angular momentum loss from MT
    }
    else{
        std::cout << "Invalid or no MT prescription for non-compact accretor. Missing orbit calculation. Return a.\n";
        return a;
    }
}


double calculateZetaThermal(Star starcopy, double percentageMassLoss, bool addMass, const programOptions &options, bool & calculateZetaThermalErrorMassFlag, bool & calculateZetaThermalErrorRadiusFlag, const gsl_rng *r){
    /*
     Calculate the radius-mass exponent zeta, assuming the star has had time to recover its thermal equilibrium. This is calculated using a fake mass loss step and recomputing the stars radius.

     Parameters
     ----------
     starcopy : Star
     Star to calculate zeta thermal for
     percentageMassLoss : double
     Percentage of mass to artificially lose from the star while attempting to calcluate zeta_thermal
     addMass : bool
     Whether to add mass or remove it when calculating zeta_thermal

     Returns
     --------
     Zeta_thermal : double
     Radius-Mass exponent Zeta for the thermal timescale

     */
    bool debugging = false;
    if(addMass){
        percentageMassLoss *= -1.0;
    }

    double zero_time = 0.0;             // We want to evolve the star for zero time, but update its radius
    bool useFakeTimestep = true;        // We want to use a fake timestep

    double radiusBeforePreEvolve = starcopy.m_Radius;
    double massBeforePreEvolve = starcopy.m_Mass;

    if(debugging){
        std::cout << "Radius before pre evolve (zeta_thermal) = " << radiusBeforePreEvolve << std::endl;
        std::cout << "Mass before pre evolve (zeta_thermal) = " << massBeforePreEvolve << std::endl;
    }
    // Allow star to respond to previous mass loss changes
    starcopy.evolveOneTimestep(zero_time, useFakeTimestep, options, r);
    starcopy.evolveOneTimestep(zero_time, useFakeTimestep, options, r); // what happens if you do this again? Seems to fix some issues, should come back and look at order of updating radius

    // Calculate properties of the star before fake mass loss
    double radiusBeforeMassLoss = starcopy.m_Radius;
    double massBeforeMassLoss = starcopy.m_MassPrev;                // Due to order of updating radius bug

    // Initialise some variables for the mass and radius after mass loss
    double radiusAfterMassLoss = 0.0;
    double massAfterMassLoss = 0.0;

    // To calculate radius-mass exponent need logarithmns of radii, masses
    double logRadiusBeforeMassLoss = 0.0;
    double logMassBeforeMassLoss = 0.0;
    double logRadiusAfterMassLoss = 0.0;
    double logMassAfterMassLoss = 0.0;

    // Variables to calculate difference in mass, difference in radius
    double deltaR = 0.0;
    double deltaM = 0.0;

    if(debugging){
        std::cout << "Radius before mass loss (zeta_thermal) = " << radiusBeforeMassLoss << std::endl;
        std::cout << "Mass before mass loss (zeta_thermal) = " << massBeforeMassLoss << std::endl;
    }
    // Check that radius, mass are sensible
    if(radiusBeforeMassLoss > 0.0){
        logRadiusBeforeMassLoss = log(radiusBeforeMassLoss);
    }
    else{
		if(calculateZetaThermalErrorRadiusFlag==false){
            std::cerr << starcopy.m_randomSeed << "\tError in calculateZetaThermal, radius before fake mass loss < 0.0" << std::endl;
			calculateZetaThermalErrorRadiusFlag = true;
		}
    }
    if(massBeforeMassLoss > 0.0){
        logMassBeforeMassLoss = log(massBeforeMassLoss);
    }
    else{
		if(calculateZetaThermalErrorMassFlag==false){
			std::cerr << starcopy.m_randomSeed << "\tError in calculateZetaThermal, mass before fake mass loss < 0.0" << std::endl;
			calculateZetaThermalErrorMassFlag = true;
		}
    }

    // Reduce mass of star by percentageMassLoss and recalculate the radius
    starcopy.m_Mass = starcopy.m_Mass * (1.0 - (percentageMassLoss/100.0));

    massAfterMassLoss = starcopy.m_Mass;

    if(debugging){std::cout << "Mass after mass loss (zeta_thermal) = " << massAfterMassLoss << std::endl;}

    // Recalculate radius of star
    starcopy.evolveOneTimestep(zero_time, useFakeTimestep, options, r);

    radiusAfterMassLoss = starcopy.m_Radius;

    if(debugging)
        std::cout << "Radius after mass loss (zeta_thermal) = " << radiusAfterMassLoss << std::endl;

    // If mass and radius are both sensible, take the log of them
    if(radiusAfterMassLoss > 0.0){
        logRadiusAfterMassLoss = log(radiusAfterMassLoss);
    }
    else{
		if(calculateZetaThermalErrorRadiusFlag==false){
            std::cerr << starcopy.m_randomSeed << "\tError in calculateZetaThermal, radius after fake mass loss < 0.0" << std::endl;
			calculateZetaThermalErrorRadiusFlag = true;
		}
    }

    if(massAfterMassLoss > 0.0){
        logMassAfterMassLoss = log(massAfterMassLoss);
    }
    else{
		if(calculateZetaThermalErrorMassFlag==false){
            std::cerr << starcopy.m_randomSeed << "\tError in calculateZetaThermal, mass after fake mass loss < 0.0" << std::endl;
			calculateZetaThermalErrorMassFlag = true;
		}
    }

    double Zeta_Thermal = (logRadiusAfterMassLoss - logRadiusBeforeMassLoss)/(logMassAfterMassLoss - logMassBeforeMassLoss);

    deltaM = massAfterMassLoss - massBeforeMassLoss;
    deltaR = radiusAfterMassLoss - radiusBeforeMassLoss;

    if(debugging){
		std::cout << "Stellartype:\t" << starcopy.m_stellarType << std::endl;
        std::cout << "Final Zeta_thermal from function = " << Zeta_Thermal << std::endl;
        std::cout << "DeltaR from function = " << deltaR << std::endl;
        std::cout << "DeltaM from function = " << deltaM << std::endl;
    }

    // Return zeta
    return Zeta_Thermal;
}

double determineConvergedMassStepZetaThermal(const Star & star_to_calculate_zeta_thermal, const programOptions &options, bool & calculateZetaThermalErrorMassFlag, bool & calculateZetaThermalErrorRadiusFlag, const gsl_rng *r){
    /*
     Calculate the radius response of a star to mass loss on the thermal timescale as characterised by the mass-radius exponent zeta_thermal

     Paramters
     ----------
     star_to_calculate_zeta_thermal : Star
     Star to calculate zeta_thermal for

     Returns
     --------
     zeta_thermal : double
     Mass-radius exponent zeta_thermal = dlnR/dlnM

     */
    int niter = 0;                                                          // Current mass iteration
    int niter_max = 10;                                                     // Maximum number of mass step iterations to make before failing loudly
    bool debugging = false;

    if(star_to_calculate_zeta_thermal.m_stellarType == NEUTRON_STAR or star_to_calculate_zeta_thermal.m_stellarType == BLACK_HOLE){
        // I guess this is not really right for NS but is correct for BHs.
        return 1.0;
    }
    else if(star_to_calculate_zeta_thermal.m_stellarType == FIRST_GIANT_BRANCH or star_to_calculate_zeta_thermal.m_stellarType == CORE_HELIUM_BURNING or star_to_calculate_zeta_thermal.m_stellarType == EARLY_ASYMPTOTIC_GIANT_BRANCH or star_to_calculate_zeta_thermal.m_stellarType == THERMALLY_PULSING_ASYMPTOTIC_GIANT_BRANCH or star_to_calculate_zeta_thermal.m_stellarType == NAKED_HELIUM_STAR_GIANT_BRANCH or star_to_calculate_zeta_thermal.m_stellarType == NAKED_HELIUM_STAR_HETZSPRUNG_GAP){
        // As described in BSE paper (Hurley et al 2002) section 2.6.1, the equation is Equaion 47 of Hurley et al 2000
        // Then modified according to Equation 56 of Hurley et al 2002
        double x = calculateRadiusGiantBranchConstantX(star_to_calculate_zeta_thermal.m_logMetallicityXi, star_to_calculate_zeta_thermal.m_randomSeed);
        double exponent = -x + (2.0 * pow((star_to_calculate_zeta_thermal.m_coreMass/star_to_calculate_zeta_thermal.m_Mass), (5.0)));

        // Debugging
        if(debugging){
        std::cout << "Zeta_thermal_giant : x = " << x << std::endl;
        std::cout << "Zeta_thermal_giant : exponent = " << exponent << std::endl;
        }
        return exponent;
    }
    else{

        double Zeta_Thermal = 0.0;                                              // Stores final value of mass-radius exponent
        double Zeta_Thermal_Previous = 0.0;                                     // Stores value of mass-radius exponent at previous mass step

        //    double initialMassLossPercentage = 1.0;                                 // Initial percentage of mass to try removing from the star to calculate zeta_thermal
        double initialMassLossPercentage = 1E-3;            // try initially smaller value?
        double currentMassLossPercentage = initialMassLossPercentage;           // Current percentage of mass to try removing from the star to calculate zeta_thermal
        //    double currentMassLossPercentage = 0.0;                                 // Check what happens at end of EAGB for 0 mass loss
        // Yeah same thing happens, you remove some (0) mass and it collapses to a black hole. Why? This is what the pre-evolution step was supposed to solve?

        double Zeta_Thermal_tolerance = 1E-3;                                   // Tolerance between mass steps

        //    if(star_to_calculate_zeta_thermal.m_stellarType == EARLY_ASYMPTOTIC_GIANT_BRANCH){
        //        std::cout << "Broken near the end of EAGB?. Constant deltaR and decreasing deltaM" << std::endl;
        //    }

        if(debugging){
            std::cout << "Previous, current stellar type (thermal) = " << star_to_calculate_zeta_thermal.m_stellarTypePrev << " " << star_to_calculate_zeta_thermal.m_stellarType << std::endl;

            if(star_to_calculate_zeta_thermal.m_stellarType == HERTZSPRUNG_GAP){
                std::cout << "HERTZSPRUNG GAP" << std::endl;
            }

            if(star_to_calculate_zeta_thermal.m_stellarType == CORE_HELIUM_BURNING){
                std::cout << "CORE HELIUM BURNING" << std::endl;
            }
        }

        // Seems to get dR/RM > 0 which means shrinks with mass loss which doesn't seem right

        while(niter < niter_max){

            // Calculate zeta using this functon
            Zeta_Thermal = calculateZetaThermal(star_to_calculate_zeta_thermal, currentMassLossPercentage, true, options, calculateZetaThermalErrorMassFlag, calculateZetaThermalErrorRadiusFlag, r);

            if((1.0 - fabs(Zeta_Thermal_Previous/ Zeta_Thermal)) < Zeta_Thermal_tolerance){
                if(debugging){
                std::cout << "Zeta_Thermal converged with:" << std::endl;
                std::cout << "Zeta_Thermal_Previous = " << Zeta_Thermal_Previous << std::endl;
                std::cout << "Zeta_Thermal = " << Zeta_Thermal << std::endl;
                std::cout << "test < tol = " << (1.0 - fabs(Zeta_Thermal_Previous/ Zeta_Thermal)) << " < " << Zeta_Thermal_tolerance << std::endl;
                }
                break;
            }
            else{
                // Reduce mass step and increment loop counter
                if(debugging){
                std::cout << "Zeta_Thermal_Previous = " << Zeta_Thermal_Previous << std::endl;
                std::cout << "Zeta_Thermal = " << Zeta_Thermal << std::endl;
                std::cout << "Mass loss percentage = " << currentMassLossPercentage << std::endl;
                std::cout << "niter = " << niter << std::endl;
                }
                currentMassLossPercentage *= 0.1;
                niter += 1;
                Zeta_Thermal_Previous = Zeta_Thermal;
                if(debugging){
                std::cout << "Zeta_Thermal_Previous = " << Zeta_Thermal_Previous << std::endl;
                std::cout << "Zeta_Thermal = " << Zeta_Thermal << std::endl;
                std::cout << "Mass loss percentage = " << currentMassLossPercentage << std::endl;
                std::cout << "niter = " << niter << std::endl;
                }
            }
        }

        if(niter > niter_max){
            std::cerr << star_to_calculate_zeta_thermal.m_randomSeed << "\tFailed to converge calculating Zeta_Thermal" << std::endl;
            if(debugging){
			// std::cout << "Failed to converge calculating Zeta_Thermal" << std::endl;
            std::cout << "Zeta_Thermal_Previous = " << Zeta_Thermal_Previous << std::endl;
            std::cout << "Zeta_Thermal = " << Zeta_Thermal << std::endl;
            std::cout << "Mass loss percentage = " << currentMassLossPercentage << std::endl;
            std::cout << "niter = " << niter << std::endl;
            }
        }

        return Zeta_Thermal;
    }
}


double calculateZetaNuclear(Star starcopy, double timestep, const programOptions &options, bool & calculateZetaNuclearErrorRadiusFlag, bool & calculateZetaNuclearErrorAgeFlag, const gsl_rng *r){
    /*
     Calculate the radius-time exponent zeta, which is derived by the nuclear evolution of the star according to SSE. This is calculated using a fake time step and recomputing the stars radius.

     Parameters
     ----------
     starcopy : Star
     Copy of the star to calculate Zeta Nuclear for

     Returns
     --------
     Zeta_nuclear : double
     Radius-Time exponent Zeta for the nuclear timescale

     */

    // We now want to use fake timesteps to estimate zeta rather than actually evolve the star
    bool debugging = false;
    bool useFakeTimestep = false;
    int niter = 0;                              // Initialise a counter
    double zero_time = 0.0;                     // Time = 0, not sure if evolution routine will like that

    double Zeta_Nuclear = 0;                    // dlnR/dt | nuclear -- radius-time exponent for star on the nuclear timescale

    double radiusBeforeTimeStep = starcopy.m_Radius;
    double ageBeforeTimeStep    = starcopy.m_Age;

    double radiusAfterTimeStep = 0.0;
    double ageAfterTimeStep = 0.0;

    // To calculate radius-time exponent need logarithmns of radii
    double logRadiusBeforeTimeStep  = 0.0;
    double logRadiusAfterTimeStep   = 0.0;

    double deltaR = 0.0;
    double deltaT = 0.0;

    // Check that radius, timestep are sensible
    if(radiusBeforeTimeStep > 0.0){
        logRadiusBeforeTimeStep = log(radiusBeforeTimeStep);
    }
    else{
		if(calculateZetaNuclearErrorRadiusFlag==false){
            std::cerr << starcopy.m_randomSeed << "\tError in calculateZetaNuclear function, radius before fake mass loss < 0.0" << std::endl;
			calculateZetaNuclearErrorRadiusFlag=true;
		}
    }
    if(ageBeforeTimeStep < 0.0){
		if(calculateZetaNuclearErrorAgeFlag==false){
            std::cerr << starcopy.m_randomSeed << "\tError in calculateZetaNuclear function, age < 0.0" << std::endl;
			calculateZetaNuclearErrorAgeFlag = true;
		}
    }

    double initialRadiusBeforeTimeStep = radiusBeforeTimeStep;
    double initialAgeBeforeTimeStep = ageBeforeTimeStep;

    starcopy.evolveOneTimestep(0.0, useFakeTimestep, options, r);

    // With the current way of doing things, this star has already lost a bunch of mass due to the mass loss caller, need to update radius of star to respond to this, before we evolve it for a short amount of time assuming no mass loss. When Alejandro moves the massLossCaller function, we will need to do things differently, and this will need fixing.
    radiusBeforeTimeStep = starcopy.m_Radius;
    ageBeforeTimeStep = starcopy.m_Age;

    starcopy.evolveOneTimestep(timestep, useFakeTimestep, options, r);

    radiusAfterTimeStep = starcopy.m_Radius;
    ageAfterTimeStep = starcopy.m_Age;

    if(radiusAfterTimeStep > 0.0){
        logRadiusAfterTimeStep = log(radiusAfterTimeStep);
    }
    else{
		if(calculateZetaNuclearErrorRadiusFlag==false){
			std::cerr << starcopy.m_randomSeed << "\tError in calculateZetaNuclear, radius after timestep < 0.0" << std::endl;
			calculateZetaNuclearErrorRadiusFlag=true;
		}
    }

    if(ageAfterTimeStep < 0.0){
		if(calculateZetaNuclearErrorAgeFlag==false){
			std::cerr << starcopy.m_randomSeed << "\tError in calculateZetaNuclear, age after timestep < 0.0" << std::endl;
			calculateZetaNuclearErrorAgeFlag=true;
		}
    }

    deltaR = radiusAfterTimeStep - initialRadiusBeforeTimeStep;
    deltaT = ageAfterTimeStep - ageBeforeTimeStep;

    Zeta_Nuclear = (logRadiusAfterTimeStep - logRadiusBeforeTimeStep)/(ageAfterTimeStep - ageBeforeTimeStep);

    // Debugging
    //std::cout << "Radius, age after timestep while attempting to calculate zeta_nuclear = " << starcopy.m_Radius << " " << starcopy.m_Age << std::endl;
    //std::cout << "DeltaR from function = " << deltaR << std::endl;
    //std::cout << "DeltaT from function = " << deltaT << std::endl;
    //std::cout << "Zeta nuclear         = " << Zeta_Nuclear << std::endl;

    // Return zeta
    return Zeta_Nuclear;

}

double determineConvergedTimestepZetaNuclear(const Star & star_to_calculate_zeta_nuclear, const programOptions &options, bool & calculateZetaNuclearErrorRadiusFlag, bool & calculateZetaNuclearErrorAgeFlag, const gsl_rng *r){
    /*
     Calculate the radius-time exponent zeta, which is derived by the nuclear evolution of the star according to SSE. This is calculated using a fake time step and recomputing the stars radius.

     Iterate evolving the star for smaller timesteps until the value of Zeta nuclear doesn't change by more than some tolerance

     Parameters
     -----------
     star_to_calculate_zeta_nuclear : Star
     Star to calculate zeta_nuclear for

     Returns
     --------
     Zeta_Nuclear : double
     dlnR/dt | nuc

     */
	 
	 bool	debugging = false;
	 
    if(star_to_calculate_zeta_nuclear.m_stellarType == BLACK_HOLE){
        return 0.0;
    }
    else{
        // We now want to use fake timesteps to estimate zeta rather than actually evolve the star
        int niter = 0;                              // Initialise a counter
        int niter_max = 10;                         // Maximum number of times to try

        double Zeta_Nuclear = 0;                    // dlnR/dt | nuclear -- radius-time exponent for star on the nuclear timescale
        double Zeta_Nuclear_Previous = 0.0;         // dlbR/dt | nuclear at previous size timestep

        double timestep_initial = 1E-3;         // Initial timestep to use
        //    double timestep_initial = 1.0;         // Initial timestep to use
        double timestep = timestep_initial;     // Current timestep to use
        double zeta_nuclear_tolerance = 1E-3;   // Fractional difference to tolerate when reducing stepsize calculating zeta_nuclear

        // Check convergence of Zeta nuclear
        // Want to make it so that the function works out what timestep makes sense to use itself

        while(niter < niter_max){

            // Calculate zeta using this functon
            Zeta_Nuclear = calculateZetaNuclear(star_to_calculate_zeta_nuclear, timestep, options, calculateZetaNuclearErrorRadiusFlag, calculateZetaNuclearErrorAgeFlag, r);

            if((1.0 - fabs(Zeta_Nuclear_Previous / Zeta_Nuclear)) < zeta_nuclear_tolerance){
                //            std::cout << "Zeta_nuclear converged with:" << std::endl;
                //            std::cout << "Zeta_Nuclear_Previous = " << Zeta_Nuclear_Previous << std::endl;
                //            std::cout << "Zeta_Nuclear = " << Zeta_Nuclear << std::endl;
                //            std::cout << "test < tol = " << (1.0 - fabs(Zeta_Nuclear_Previous / Zeta_Nuclear)) << " < " << zeta_nuclear_tolerance << std::endl;
                break;
            }
            else{
                // Reduce timestep and increment loop counter
                //            std::cout << "zeta_nuclear_previous = " << Zeta_Nuclear_Previous << std::endl;
                //            std::cout << "timestep = " << timestep << std::endl;
                //            std::cout << "niter = " << niter << std::endl;
                timestep *= 0.1;
                niter += 1;
                Zeta_Nuclear_Previous = Zeta_Nuclear;
                //            std::cout << "zeta_nuclear_previous = " << Zeta_Nuclear_Previous << std::endl;
                //            std::cout << "timestep = " << timestep << std::endl;
                //            std::cout << "niter = " << niter << std::endl;
            }
        }

        if(niter > niter_max){
            std::cerr << star_to_calculate_zeta_nuclear.m_randomSeed << "\tFailed to converge calculating Zeta_Nuclear" << std::endl;
			if(debugging){
				std::cout << "zeta_nuclear_previous = " << Zeta_Nuclear_Previous << std::endl;
				std::cout << "timestep = " << timestep << std::endl;
				std::cout << "niter = " << niter << std::endl;
			}
        }

        // Return zeta
        return Zeta_Nuclear;
    }
}

int massTransferCase(int stellarType){
	// Function to determine if the star which initiated mass transfer will go through case A, B or C.
    if(stellarType == MS_MORE_THAN_07 or stellarType == MS_LESS_THAN_07 or stellarType == NAKED_HELIUM_STAR_MS){
        return CASE_A;
    }
    else if(stellarType == HERTZSPRUNG_GAP or stellarType == NAKED_HELIUM_STAR_HETZSPRUNG_GAP){
        return CASE_B;
    }
    else if((stellarType >= FIRST_GIANT_BRANCH and stellarType <= THERMALLY_PULSING_ASYMPTOTIC_GIANT_BRANCH) or NAKED_HELIUM_STAR_GIANT_BRANCH){
        return CASE_C;
    }
    else{
        return NO_MASS_TRANSFER;
    }
}


int initialMassTransferCase(int stellarType, const programOptions & options){
        // Function to determine if the star that just filled it's RL and initiated mass transfer goes through case A, B or C.
    if(stellarType == MS_MORE_THAN_07 or stellarType == MS_LESS_THAN_07 or stellarType == NAKED_HELIUM_STAR_MS){
        return CASE_A;
    }
    else if(stellarType == HERTZSPRUNG_GAP or stellarType == NAKED_HELIUM_STAR_HETZSPRUNG_GAP){
        return CASE_B;
    }
    else if((stellarType >= FIRST_GIANT_BRANCH and stellarType <= THERMALLY_PULSING_ASYMPTOTIC_GIANT_BRANCH) or (stellarType == NAKED_HELIUM_STAR_GIANT_BRANCH)){
        return CASE_C;
    }
    else{
        return NO_MASS_TRANSFER;
    }
}


//Mass Transfer functions for deMink mass transfer
double gammaAngularMomentumLoss(double md, double ma, const programOptions & options){
    /* Function to determine the value of gamma (as in Pols's notes) or jloss (as in Belczynski et al. 2008, which is the fraction of specific angular momentum with which the non-accreted mass leaves the system.
    */
	bool	debugging = false;

	double	gamma = 1.0;
	
    if(options.massTransferAngularMomentumLossPrescription == JEANS_MASS_TRANSFER_ANGULAR_MOMENTUM_LOSS){
		if(debugging){std::cout << "JEANS_MASS_TRANSFER_ANGULAR_MOMENTUM_LOSS prescription chosen for angular momentum loss prescription." << std::endl;}
		gamma = ma/md;
    }
    else if(options.massTransferAngularMomentumLossPrescription == ISOTROPIC_RE_EMISSION_MASS_TRANSFER_ANGULAR_MOMENTUM_LOSS){
		if(debugging){std::cout<<"ISOTROPIC_RE_EMISSION_MASS_TRANSFER_ANGULAR_MOMENTUM_LOSS prescription chosen for angular momentum loss prescription."<<std::endl;}
        gamma = md/ma;
    }
    else if(options.massTransferAngularMomentumLossPrescription == CIRCUMBINARY_RING_MASS_TRANSFER_ANGULAR_MOMENTUM_LOSS){
        // Based on the assumption that a_ring ~= a*sqrt(2) in Evernote Notebook based on talks with deMink, or as tricky people call it, "private communication"
		if(debugging){std::cout<<"CIRCUMBINARY_RING_MASS_TRANSFER_ANGULAR_MOMENTUM_LOSS chosen prescription for angular momentum loss prescription."<<std::endl;}
        gamma = sqrt(2)*pow(md+ma,2.0)*pow(md*ma,-1.0);
    }
	else if(options.massTransferAngularMomentumLossPrescription == ARBITRARY_MASS_TRANSFER_ANGULAR_MOMENTUM_LOSS){
		if(debugging){std::cout<<"ARBITRARY_MASS_TRANSFER_ANGULAR_MOMENTUM_LOSS prescription chosen for angular momentum loss prescription."<<std::endl;}
        gamma = options.massTransferJloss;	
	}
    else{
		if(debugging){std::cout<<"No prescription chosen for angular momentum loss prescription. Return default value gamma = 1"<<std::endl;}
        gamma = 1.0;
    }

	if(debugging){std::cout << "Gamma = " << gamma << std::endl;}
    return  gamma;
}

double massAcceptanceRate(double Ra, double thermalRateAccretor, double massRateDonor, int stellarTypeAccretor, int stellarTypeDonor, const programOptions & options, double &fa, unsigned long m_randomSeed){
    /*
    This function does two things. It firstly determines the maximum mass acceptance rate of the accretor star during mass transfer. It then uses this, along with the mass donation rate to determine the accretion efficieny parameter fa, which it then updates.

    We determine the maximum acceptance rate of the accretor star during mass transfer, based on the stellar type. We assume that mass transfer is Eddington limited for BHs and NSs. We use the formalism of Nomoto/Claeys for WDs. For non compact objects:

     1) Kelvin-Helmholtz (thermal) timescale if THERMALLY_LIMITED_MASS_TRANSFER
     2) Choose a fraction of the mass rate that will be effectively accreted for FIXED_FRACTION_MASS_TRANSFER (as in StarTrack)
     3) Disk vs impact accretion for CENTRIFUGALLY_LIMITED_MASS_TRANSFER

    Parameters
    -----------
    Ra : double
        Radius of the accretor
    thermalRateAccretor : double
        Thermal rate of accretor (M/tau_thermal)
    massRateDonor : double
        Mass transfer rate of donor
    stellarTypeAccretor : int
        Stellar type of accretor according to Hurley+ 2000
    stellarTypeDonor : int
        Stellar type of donor according to Hurley+ 2000
    options : programOptions
        Object containing user specified program options
    fa : double
        Mass transfer accretion efficieny, updated by this function

    Returns
    --------
    rateAccretor : double
        Maximum accretion rate of accretor. Value not used elsewhere
     */

    // Declare variables
    double  MDotEddington = 0.0;
    bool    debugging = false;

    // If accreting star is a main sequence star, or evolved off of the main sequence but not yet remnant:
    if(stellarTypeAccretor >= MS_LESS_THAN_07 and stellarTypeAccretor <= NAKED_HELIUM_STAR_GIANT_BRANCH){

        if(debugging){
			std::cout << "massAcceptanceRate function: Thermally limited mass transfer" << std::endl;
		}

        // If using thermally limited mass transfer:
        if(options.massTransferAccretionEfficiencyPrescription == THERMALLY_LIMITED_MASS_TRANSFER){

            if(debugging){
				std::cout << "Thermally limited mass transfer" << std::endl;
				std::cout << "thermalRateAccretor: " << thermalRateAccretor << std::endl;
				std::cout << "tmassRateDonor: " << massRateDonor << std::endl;
			}

            // If the accretor cannot accrete at the rate the donor is providing mass:
            if(thermalRateAccretor < massRateDonor){

                if(debugging){
					std::cout << "thermalRateAccretor < massRateDonor" << std:: endl;
                    std::cout << "fa0: " << fa << std::endl;
                    std::cout << "fa: " << (thermalRateAccretor/massRateDonor) << std::endl;
                }

				if(options.massTransferThermallyLimitedVariation == THERMAL_C_FACTOR){

					if(debugging){
						std::cout << "THERMAL_C_FACTOR" << std::endl;
						std::cout << "Factor C: " << options.massTransferCParameter << std::endl; // Default value is 10.0 (Hurley+ 2002)
						std::cout << "Calculated fa: " << options.massTransferCParameter*(thermalRateAccretor/massRateDonor) << std::endl;
					}

					fa = fmin(1.0,options.massTransferCParameter*(thermalRateAccretor/massRateDonor));
				}
				else if(options.massTransferThermallyLimitedVariation == THERMAL_RADIUS_TO_ROCHELOBE){

					if(debugging){
						std::cout << "THERMAL_RADIUS_TO_ROCHELOBE" << std::endl;
                        std::cout << "Factor C: " << options.massTransferCParameter << std::endl; // Default value is 1.0
						std::cout << "Calculated fa: " << thermalRateAccretor/massRateDonor << std::endl;
					}

                    fa = fmin(1.0,options.massTransferCParameter*(thermalRateAccretor/massRateDonor));
				}
                else{

                    if(debugging){
                        std::cerr << m_randomSeed << "\tUnrecognised mass transfer thermallly limited variation" << std::endl;
                    }

                }

                if(debugging){

					std::cout << "fa': " << fa << std::endl;
					std::cout << "calculated, massRateDonor: " << options.massTransferCParameter*thermalRateAccretor << "\t" << massRateDonor << std::endl;

				}

                return fmin(options.massTransferCParameter*thermalRateAccretor, massRateDonor);
            }
            // Accretor can accrete more than the donor is providing
            else{
				if(debugging){
					std::cout << "thermalRateAccretor >= massRateDonor. Return 1.0" << std:: endl;
                }
                fa = 1.0;
                return massRateDonor;
            }

        }
        // If using fixed fraction of mass accreted, as in StarTrack
        else if(options.massTransferAccretionEfficiencyPrescription == FIXED_FRACTION_MASS_TRANSFER){
				if(debugging){
					std::cout << "massAcceptanceRate function: fixed fraction mass transfer" << std::endl;
				}
                return fmin(massRateDonor, fa*massRateDonor);
        }
        // If using centrifugally limited mass transfer
        else if(options.massTransferAccretionEfficiencyPrescription == CENTRIFUGALLY_LIMITED_MASS_TRANSFER){
				if(debugging){
					std::cout << "massAcceptanceRate function: centrifugally limited mass transfer. Not yet implemented." << std::endl;
				}
                return 0;
        }
        // Catch
        else{
                std::cerr << m_randomSeed << "\tError in massAcceptanceRate. Shouldn't get here.\n";
                return 0;
        }

    }
    // If accretor is a white dwarf
    else if(stellarTypeAccretor == HELIUM_WHITE_DWARF or stellarTypeAccretor == CARBON_OXYGEN_WHITE_DWARF or stellarTypeAccretor == OXYGEN_NEON_WHITE_DWARF){

        if(debugging){
            std::cout << "WD accretor rate as Eddington. Nomoto/Claeys to be implemented.\n";
        }
        MDotEddington = EddingtonCriticalRate(Ra, stellarTypeAccretor, stellarTypeDonor);

        if(MDotEddington < massRateDonor){
            fa = fmin(1.0,MDotEddington/massRateDonor);
            return fmin(MDotEddington, massRateDonor);
        }
        else{
            fa = 1.0;
            return massRateDonor;
        }

    }
    // If accretor is a compact remnant NS/BH
    else if(stellarTypeAccretor == NEUTRON_STAR or stellarTypeAccretor == BLACK_HOLE){
        
            //Determine the maximum accretion rate for NS/BH
            //based on Eddington luminosity   Units Msun Myr-1
            MDotEddington = EddingtonCriticalRate(Ra, stellarTypeAccretor, stellarTypeDonor);
            //Multiply Eddington limit by double (set in pythonSubmit)
            // EddingtonAccretionFraction == 1  equals Eddington limit   
            // EddingtonAccretionFraction > 1   is super Eddington accretion
            // EddingtonAccretionFraction 0-1   is fraction of Eddington
            // EddingtonAccretionFraction < 0   is error
            double maxAccretrionRate = MDotEddington * options.eddingtonAccretionFactor; //Msun Myr-1
            if(options.eddingtonAccretionFactor < 0){
                std::cerr << m_randomSeed 
                          << "\tError in value SuperEddingtonAccretion. Shouldn't get here. Return 0.\n";
            }
            if(debugging){
                std::cout << "\nMassAcceptanceRate function\n";
                std::cout << "Ra: \t" << Ra << std::endl;
                std::cout << "MDotEdd: \t" << MDotEddington << std::endl;
                std::cout << "EddingtonAccretionFraction: \t" 
                          << options.eddingtonAccretionFactor << std::endl;
                std::cout << "massRateDonor: \t" << massRateDonor << std::endl;
            }

            

            if(maxAccretrionRate < massRateDonor){ 
                // if the mass rate of donor is larger then the maximum Accretion rate
                // update pointer fa, (by construction should always be smaller than 1
                fa = fmin(1.0,(maxAccretrionRate)/massRateDonor); //fraction 0-1
                // return the massAccretion rate of the accretor
                // by construction should be maxAccretrionRate
                return  fmin(maxAccretrionRate, massRateDonor); //Msun Myr-1
            }
            else{ 
                // if the mass rate of donor is smaller then the maximum Accretion rate
                // update pointer fa, (by construction should always be 1
                 fa = fmin(1.0,(maxAccretrionRate)/massRateDonor); //fraction 0-1
                // return the massAccretion rate of the accretor
                // by construction should be massRateDonor
                return  fmin(maxAccretrionRate, massRateDonor);//Msun Myr-1
            }
        }//Close accretor is a compact remnant NS/BH

    else{
        std::cerr << m_randomSeed << "\tError in massAcceptanceRate. Shouldn't get here. Return 0.\n";
        return 0;

    }

}

double massTransferRateAccretorBelczynski(double Ra, double fa, double MDot, int stellarTypeDonor, int stellarTypeAccretor, unsigned long randomSeed){
		// ALEJANDRO - 03/03/2017 - Function to calculate the mass acceptance rate for Belczynski MT accretors
		
		double 	massTransferRate = NEVER_SET;
		double 	massTransferRateEddington = EddingtonCriticalRate(Ra, stellarTypeAccretor, stellarTypeDonor);
		bool	debugging = false;
		
		if(debugging){
			std::cout << "Part of re-written Belczynski function. Still needs to be tested." << std::endl;
			std::cout << "massTransferRateAccretorBelczynski function:" << std::endl;
			std::cout << "Radius accretor []:\t" << Ra << std::endl;
			std::cout << "fa:\t" << fa << std::endl;
			std::cout << "MDot []:\t" << MDot << std::endl;
			std::cout << "MDotEdd []:\t" << massTransferRateEddington << std::endl;
			std::cout << "Stellar type donor:\t" << stellarTypeDonor << std::endl;
			std::cout << "Stellar type accretor:\t" << stellarTypeAccretor << std::endl;
		}
		
		if(stellarTypeAccretor <  HELIUM_WHITE_DWARF){
			if(debugging){
				std::cout << "Accretor is a non-degerate object" << std::endl;
			}
			massTransferRate = fa*MDot;
		}
		else if((stellarTypeAccretor >= HELIUM_WHITE_DWARF)and(stellarTypeAccretor <= BLACK_HOLE)){
			if(debugging){
				std::cout << "Accretor is a degerate object" << std::endl;
			}
			massTransferRate = fmin(massTransferRateEddington,fa*MDot);
		}		
		else{
			std::cerr << randomSeed << "\tError in massTransferRateAccretorBelczynski. Accretor is a massless remnant." << std::endl;
		}
		
		return massTransferRate;		
}

double adaptiveRocheLobeOverFlow(Star starcopyDonor, double Ma, double semiMajorAxis, const programOptions &options, double &fa, double &jloss, unsigned long m_randomSeed, bool &m_aRocheLobeOverFlowErrorPrinting, bool & adaptiveRocheLobeOverFlowErrorMassFlag, bool & adaptiveRocheLobeOverFlowErrorRadiusFlag, const gsl_rng *r){
    /*
     Calculate the ammount of mass that a star with radiative envelope needs to lose in order to lose in order to just fill its Roche Lobe. Based on the function –rlof_method in binary_c.

     Parameters
     ----------
     starcopy : Star
     Star to calculate zeta thermal for
     percentageMassLoss : double
     Percentage of mass to artificially lose from the star while attempting to calcluate zeta_thermal
     addMass : bool
     Whether to add mass or remove it when calculating zeta_thermal

     Returns
     --------
     deltaM : double
     Amount of mass lost from the donor in order to barely fill its own Roche Lobe
     */
	 	 

    double Md       = starcopyDonor.m_Mass;                 // Initial donor mass
    double a        = semiMajorAxis;                        // Initial semi-major axis
    double J        = (Md*Ma)*sqrt(G1*(Md+Ma)*a)/(Md+Ma);   // Initial orbital angular momentum
    double RL       = a*rocheLobeRadius(Md, Ma);            // Initial roche lobe
    double RLRsol   = RL*AUToRsol;                          // Initial roche lobe in Rsol
    Star   starcopyDonorTemp = starcopyDonor;

    double Mdfinal  = 0.0;                                  // Final
    double Mafinal  = 0.0;
    double dJ       = 0.0;
    double afinal   = 0.0;
    double RLfinal  = 0.0;
    double RLfinalRsol = 0.0;
    double RLfinalRsolTemp = 0.0;

    bool   	debugging = false;
	bool	fastPhasePrintError = false;

    if(debugging){std::cout << "\nAdaptive RLOF\n";}

    double maxPercentageMassLoss = options.maxPercentageAdaptiveMassTransfer;

    if(debugging){std::cout << "maxPercentageMassLoss: " << maxPercentageMassLoss << std::endl;}

    double zero_time = 0.0;                // We want to evolve the star for zero time, but update its radius
    bool   useFakeTimestep = true;         // We want to use a fake timestep
    bool   fixedRL = false;                // Default false, which deals with both conservative and non-conservative MT
    fixedRL = true;
    //debugging= true;

    // Allow star to respond to previous mass loss changes
     starcopyDonor.evolveOneTimestep(zero_time, useFakeTimestep, options, r);

     // Calculate properties of the star before fake mass loss
     double radiusBeforeMassLoss = starcopyDonor.m_RadiusPrev;
     double massBeforeMassLoss = starcopyDonor.m_MassPrev;

    if(debugging){std::cout << "mPrev, rPrev: " << massBeforeMassLoss  << " " << radiusBeforeMassLoss << std::endl;}

     // Initialise some variables for the mass and radius after mass loss
     double radiusAfterMassLoss = 0.0;
     double radiusAfterMassLossTemp = 0.0;
     double massAfterMassLoss = 0.0;
     double massAfterMassLossTemp = 0.0;

     // Variables to calculate difference in mass, difference in radius
     double deltaR = 0.0;
     double deltaM = 0.0;

     // Check that radius, mass are sensible
    if(radiusBeforeMassLoss <= 0.0){
		if(adaptiveRocheLobeOverFlowErrorRadiusFlag==false){
			std::cerr << starcopyDonor.m_randomSeed << "\tError in adaptiveRocheLobeOverFlow function, radius before fake mass loss < 0.0" << std::endl;
			adaptiveRocheLobeOverFlowErrorRadiusFlag = true;
		}
	}
    if(massBeforeMassLoss <= 0.0){
		if(adaptiveRocheLobeOverFlowErrorMassFlag==false){
			std::cerr << starcopyDonor.m_randomSeed << "\tError in adaptiveRocheLobeOverFlow function, mass before fake mass loss < 0.0" << std::endl;
			adaptiveRocheLobeOverFlowErrorMassFlag = true;
		}
	}

    // Error values
    int niter=100;
    double errorThreshold = 0.001;
    double absoluteError = 1.0;
    double absoluteErrorTemp = 1.0;
    double percentageMassLoss=maxPercentageMassLoss/niter;
//    double percentageMassLoss=50.0/niter;

    for(int i=1; i <= niter; i++){

//		std::cout << "i inside aRLOF: " << i << std::endl;
        starcopyDonorTemp = starcopyDonor;

        if(fixedRL == false){
            // Calculate new Roche Lobe Radius according to the mass lost
            Mdfinal  = (1.0 - percentageMassLoss)*Md;
            Mafinal  = Ma + fa*percentageMassLoss*Md;
            dJ       = jloss*((1.0-fa)*(Mdfinal-Md)/(Md+Ma))*J;
            afinal   = (Mdfinal+Mafinal)*pow(J+dJ,2.0)/(G1*Mdfinal*Mdfinal*Mafinal*Mafinal);
            RLfinal  = afinal*rocheLobeRadius(Mdfinal, Mafinal);
            RLfinalRsolTemp = RLfinal*AUToRsol;
        }
        else{
            Mdfinal  = (1.0 - percentageMassLoss)*Md;
            Mafinal  = Ma + fa*percentageMassLoss*Md;
        }
        // Reduce mass of star by percentageMassLoss and recalculate the radius
        starcopyDonorTemp.m_Mass = starcopyDonorTemp.m_Mass * (1.0 - percentageMassLoss);
        starcopyDonorTemp.m_Mass0 = starcopyDonorTemp.m_Mass0 * (1.0 - percentageMassLoss);
        massAfterMassLossTemp = starcopyDonorTemp.m_Mass;

        starcopyDonorTemp.m_Age = timeChangeAfterMassLoss(starcopyDonorTemp.m_an_coefficients, starcopyDonorTemp.m_Mass0, starcopyDonorTemp.m_Mass, starcopyDonorTemp.m_coreMass, starcopyDonorTemp.m_Age, starcopyDonorTemp.m_timescales, starcopyDonorTemp.m_stellarType, starcopyDonorTemp.m_logMetallicityXi, false);

        // Recalculate radius of star
        starcopyDonorTemp.evolveOneTimestep(zero_time, useFakeTimestep, options, r);
        radiusAfterMassLossTemp = starcopyDonorTemp.m_Radius;

        if(radiusAfterMassLoss > radiusBeforeMassLoss){
			if(m_aRocheLobeOverFlowErrorPrinting == false){
				m_aRocheLobeOverFlowErrorPrinting = true;
				std::cerr << m_randomSeed << "\tError in adaptive RL function, slow phase. Star radii increases after removing mass. Stellartype, Mass, Mass':\t" << starcopyDonorTemp.m_stellarType << ", " << Md << ", " << massAfterMassLossTemp << std::endl;
			}
        }

        if(fixedRL == false){
            absoluteErrorTemp = fabs(radiusAfterMassLossTemp-RLfinalRsolTemp);

            if(i==1){
                absoluteError = absoluteErrorTemp;
                RLfinalRsol = RLfinalRsolTemp;
                massAfterMassLoss = massAfterMassLossTemp;
                radiusAfterMassLoss = radiusAfterMassLossTemp;
            }

            if(absoluteErrorTemp <= errorThreshold){
                if(absoluteErrorTemp < absoluteError){
                    absoluteError = absoluteErrorTemp;
                    RLfinalRsol = RLfinalRsolTemp;
                    massAfterMassLoss = massAfterMassLossTemp;
                    radiusAfterMassLoss = radiusAfterMassLossTemp;
                }
            }
        }
        else{
            absoluteErrorTemp = fabs(radiusAfterMassLossTemp-RLRsol);

            if(i==1){
                absoluteError = absoluteErrorTemp;
                RLfinalRsol = RLRsol;
                massAfterMassLoss = massAfterMassLossTemp;
                radiusAfterMassLoss = radiusAfterMassLossTemp;
            }

            if(absoluteErrorTemp <= errorThreshold){
                if(absoluteErrorTemp < absoluteError){
                    absoluteError = absoluteErrorTemp;
                    RLfinalRsol = RLRsol;
                    massAfterMassLoss = massAfterMassLossTemp;
                    radiusAfterMassLoss = radiusAfterMassLossTemp;
                }
            }

        }

        percentageMassLoss+=(maxPercentageMassLoss/niter);

        if(debugging){
            if(fixedRL == false){
            std::cout << "m', r', RL', dif, rET: " << "\t" << massAfterMassLossTemp<< "\t" << radiusAfterMassLossTemp << "\t" <<RLfinalRsolTemp << "\t" << radiusAfterMassLossTemp-RLfinalRsolTemp << "\t" << absoluteErrorTemp << std::endl;
            }
            else{
            std::cout << "m', r', RL, dif, rET: " << "\t" << massAfterMassLossTemp<< "\t" << radiusAfterMassLossTemp << "\t" <<RLRsol << "\t" << radiusAfterMassLossTemp-RLRsol << "\t" << absoluteErrorTemp << std::endl;
            }
        }
    }


    if(debugging){
        std::cout << "m', r', RL', dif, rET: " << "\t" << massAfterMassLoss<< "\t" << radiusAfterMassLoss << "\t" <<RLfinalRsol << "\t" << radiusAfterMassLoss-RLfinalRsol << "\t" << absoluteError << std::endl;
    }


    // If mass and radius are both sensible, take the log of them
	if(radiusAfterMassLoss <= 0.0){
		if(adaptiveRocheLobeOverFlowErrorRadiusFlag==false){
			std::cout << "Error in adaptiveRocheLobeOverFlow, radius after fake mass loss < 0.0" << std::endl;
			adaptiveRocheLobeOverFlowErrorRadiusFlag = true;
		}
	}
    if(massAfterMassLoss <= 0.0){
		if(adaptiveRocheLobeOverFlowErrorMassFlag==false){
			std::cout << "Error in adaptiveRocheLobeOverFlow, mass after fake mass loss < 0.0" << std::endl;
			adaptiveRocheLobeOverFlowErrorMassFlag = true;
		}
	}

     deltaM = massAfterMassLoss - massBeforeMassLoss;
     deltaR = radiusAfterMassLoss - radiusBeforeMassLoss;

    if(debugging){
        std::cout << "a, a':" << a << " " << afinal << std::endl;
        std::cout << "RL, RL':" << RLRsol << " " << RLfinalRsol << std::endl;
        std::cout << "R, R': " << radiusBeforeMassLoss << " " << radiusAfterMassLoss << std::endl;
        std::cout << "M, M': " << massBeforeMassLoss << " " << massAfterMassLoss << std::endl;
        std::cout << "dM: " << options.massTransferAdaptiveAlphaParameter*deltaM << std::endl;
    }

    // Return mass
    return deltaM;
}



double massTransferFastPhaseCaseA(Star starcopyDonor, double Ma, double Ra, int stellarTypeAccretor, double semiMajorAxis, double orbitalVelocity, double eccentricity, const programOptions & options, double &fa, double &jloss, double thermalRateDonor, double thermalRateAccretor, double dt, unsigned long m_randomSeed, bool &m_aFastPhaseRocheLobeOverFlowErrorPrinting, bool &massTransferFastPhaseCaseAErrorMassFlag, bool &massTransferFastPhaseCaseAErrorRadiusFlag, const gsl_rng *r){
    /*
     Calculate the ammount of mass that a star with radiative envelope needs to lose in order to lose in order to just fill its Roche Lobe. Based on the function –rlof_method in binary_c.
	
     Parameters
     ----------
     starcopy : Star
     Star to calculate zeta thermal for
     percentageMassLoss : double
     Percentage of mass to artificially lose from the star while attempting to calcluate zeta_thermal
     addMass : bool
     Whether to add mass or remove it when calculating zeta_thermal

     Returns
     --------
     deltaM : double
     Amount of mass lost from the donor in order to barely fill its own Roche Lobe
	 
	 For fast phase case A MT, we solve the orbit numerically for the thermal timescale of the donor.

     */

    double Md       = starcopyDonor.m_Mass;                 // Initial donor mass
    double a        = semiMajorAxis;                        // Initial semi-major axis
	double w		= orbitalVelocity;						// Initial orbital velocity
	double e		= eccentricity;							// Initial eccentricity
    double J        = (Md*Ma)*sqrt(G1*(Md+Ma)*a)/(Md+Ma);   // Initial orbital angular momentum
    double RL       = a*rocheLobeRadius(Md, Ma);            // Initial roche lobe
    double RLRsol   = RL*AUToRsol;                          // Initial roche lobe in Rsol
	double Fa		= fa;									// Initial/original values of fa
	double Jloss 	= jloss;								// Initial/original values of Jloss
    Star   starcopyDonorTemp = starcopyDonor;

    double Mdfinal  		= 0.0;                                  // Final
    double Mafinal  		= 0.0;
	double dM				= 0.0;
	double MDot				= 0.0;
    double dJ       		= 0.0;
    double afinal   		= 0.0;
    double RLfinal  		= 0.0;
    double RLfinalRsol 		= 0.0;
    double RLfinalRsolTemp 	= 0.0;
	double aChosen			= 0.0;
	double aChosenTemp		= 0.0;

    bool	debugging = false;
    //debugging = true;
	

    if(debugging){
		std::cout << "\nfastPhase RLOF\n";
		std::cout << "a\t[AU]:\t" << a << std::endl;
		std::cout << "a\t[Rsol]:\t" << a*AUToRsol << std::endl;
		std::cout << "Md\t[Msol]:\t" << Md << std::endl;
		std::cout << "Ma\t[Msol]:\t" << Ma << std::endl;
		std::cout << "e :\t" << e << std::endl;
		std::cout << "dt [Myr]:\t" << dt << std::endl;
		std::cout << "J\t[]:\t" << J << std::endl;
		std::cout << "fa:\t" << fa << std::endl;
		std::cout << "jloss:\t" << jloss << std::endl;
		std::cout << "Rd\t[]:\t" << starcopyDonor.m_Radius << std::endl;
		std::cout << "Ra\t[]:\t" << Ra << std::endl;	
		std::cout << "thrateD[] :\t" << thermalRateDonor << std::endl;
		std::cout << "thrateA[]:\t" << thermalRateAccretor << std::endl;
		std::cout << "RL\t[AU]:\t" << RL << std::endl;
		std::cout << "RL\t[Rsol]:\t" << RLRsol << std::endl;
		
//					std::cout << "MdDot:\t" << MdDot << std::endl;
		
		if(stellarTypeAccretor <= NAKED_HELIUM_STAR_GIANT_BRANCH){std::cout << "Fast phase case A aRLOF, NON-degenerate accretor." << std::endl;}
		else if(stellarTypeAccretor >= HELIUM_WHITE_DWARF and stellarTypeAccretor <= BLACK_HOLE){std::cout << "Fast phase case A aRLOF, degenerate accretor." << std::endl;}
		else{std::cerr << m_randomSeed << "\tError in massTransferFastPhaseCaseA function. Massless remnant accreting mass. Return 'a'." << std::endl;}
	}

	// In fast phase case A MT, we allow to strip almost all of the star in order to see if it fits in it's Roche lobe.
    double maxPercentageMassLoss = 0.99;

    if(debugging){std::cout << "maxPercentageMassLoss: " << maxPercentageMassLoss << std::endl;}

    double zero_time = 0.0;                // We want to evolve the star for zero time, but update its radius
    bool   useFakeTimestep = true;         // We want to use a fake timestep
    bool   fixedRL = false;                // Default false, which deals with both conservative and non-conservative MT

    // Allow star to respond to previous mass loss changes
    starcopyDonor.evolveOneTimestep(zero_time, useFakeTimestep, options, r);

    // Calculate properties of the star before fake mass loss
	//    double radiusBeforeMassLoss = starcopyDonor.m_RadiusPrev;
	//    double massBeforeMassLoss = starcopyDonor.m_MassPrev;
    double radiusBeforeMassLoss = starcopyDonor.m_Radius;
    double massBeforeMassLoss = starcopyDonor.m_Mass;
    if(debugging){std::cout << "m, r: " << massBeforeMassLoss  << " " << radiusBeforeMassLoss << std::endl;}

    // Initialise some variables for the mass and radius after mass loss
    double radiusAfterMassLoss = 0.0;
    double radiusAfterMassLossTemp = 0.0;
    double massAfterMassLoss = 0.0;
    double massAfterMassLossTemp = 0.0;

    // Variables to calculate difference in mass, difference in radius
    double deltaR = 0.0;
    double deltaM = 0.0;

    // Check that radius, mass are sensible
    if(radiusBeforeMassLoss <= 0.0){
		if(massTransferFastPhaseCaseAErrorRadiusFlag==false){
			std::cerr << "Error in massTransferFastPhaseCaseA, radius before fake mass loss < 0.0" << std::endl;
			massTransferFastPhaseCaseAErrorRadiusFlag = true;
		}
	}
    if(massBeforeMassLoss <= 0.0){
		if(massTransferFastPhaseCaseAErrorMassFlag==false){
			std::cerr << "Error in massTransferFastPhaseCaseA, mass before fake mass loss < 0.0" << std::endl;
			massTransferFastPhaseCaseAErrorMassFlag = true;
		}
	}

    // Error values
    int niter=100;	// niter = 100
    double absoluteErrorThreshold = 0.001;
    double absoluteError = 1.0;
    double absoluteErrorTemp = 1.0;
    double relativeErrorThreshold = 0.1;
    double relativeError = 1.0;
    double relativeErrorTemp = 1.0;
    double percentageMassLoss=maxPercentageMassLoss/niter;
	

    for(int i=1; i <= niter; i++){

        starcopyDonorTemp = starcopyDonor;
			
			fa=Fa;
			jloss=Jloss;

            // Calculate new Roche Lobe Radius according to the mass lost
			dM		= percentageMassLoss*Md;
			MDot	= -dM/dt;
            Mdfinal  = (1.0 - percentageMassLoss)*Md;
            Mafinal  = Ma + fa*percentageMassLoss*Md;

			if(debugging){
				std::cout << "dM [Msol]:\t\t" << dM << std::endl;
				std::cout << "Md' [Msol]:\t\t" << Mdfinal << std::endl;
				std::cout << "Ma' [Msol]:\t\t" << Mafinal << std::endl;	
				std::cout << "MDot [Msol Myr-1]:\t" << MDot<< std::endl;			
			}

		    if(stellarTypeAccretor <= NAKED_HELIUM_STAR_GIANT_BRANCH){
				// For fast phase, dt = tauThermalDonor 
				afinal = calculateMassTransferOrbitNonDegenerateAccretor(Md, starcopyDonor.m_envMass, Ma, fa, MDot, a, w, e, dt, dt, J, jloss, Ra, thermalRateAccretor, thermalRateDonor, stellarTypeAccretor, starcopyDonor.m_stellarType, options, m_randomSeed);       // Eq (33) Belczynski et al. (2008)
				// double   calculateMassTransferOrbitNonDegenerateAccretor(double md, double mdEnv, double ma, double fa, double MdDot, double a, double worb, double e, double dt, double tau, double JPrev, double jloss, double Ra, double thermalRateAccretor, double thermalRateDonor, double stellarTypeAccretor, double stellarTypeDonor, const programOptions & options);
			}
            else if(stellarTypeAccretor >= HELIUM_WHITE_DWARF and stellarTypeAccretor <= BLACK_HOLE){
				afinal = calculateMassTransferOrbitDegenerateAccretor(Md, Ma, fa, MDot, a, w, e, dt, J, jloss, Ra, thermalRateAccretor, thermalRateDonor, stellarTypeAccretor, starcopyDonor.m_stellarType, options, m_randomSeed);				
//                    afinal = calculateMassTransferOrbitDegenerateAccretor(md, ma, fa, MdDot, a, w, e, dt, J, jloss, Ra, thermalRateAccretor, thermalRateDonor, stellarTypeAccretor, stellarTypeDonor,options);
				//double   calculateMassTransferOrbitDegenerateAccretor(double md, double ma, double fa, double MdDot, double a, double worb, double e, double dt, double JPrev, double jloss, double Ra, double thermalRateAccretor, double thermalRateDonor, double stellarTypeAccretor, double stellarTypeDonor, const programOptions & options);
			}
            else{
                afinal = a;
			}

			aChosenTemp	= afinal;
			RLfinal  = afinal*rocheLobeRadius(Mdfinal, Mafinal);
			RLfinalRsolTemp = RLfinal*AUToRsol;

        // Reduce mass of star by percentageMassLoss and recalculate the radius
        starcopyDonorTemp.m_Mass = starcopyDonorTemp.m_Mass * (1.0 - percentageMassLoss);
        starcopyDonorTemp.m_Mass0 = starcopyDonorTemp.m_Mass0 * (1.0 - percentageMassLoss);
        massAfterMassLossTemp = starcopyDonorTemp.m_Mass;

        // Modify Age
        starcopyDonorTemp.m_Age = timeChangeAfterMassLoss(starcopyDonorTemp.m_an_coefficients, starcopyDonorTemp.m_Mass0, starcopyDonorTemp.m_Mass, starcopyDonorTemp.m_coreMass, starcopyDonorTemp.m_Age, starcopyDonorTemp.m_timescales, starcopyDonorTemp.m_stellarType, starcopyDonorTemp.m_logMetallicityXi, false);

        // Recalculate radius of star
        starcopyDonorTemp.evolveOneTimestep(zero_time, useFakeTimestep, options, r);
        radiusAfterMassLossTemp = starcopyDonorTemp.m_Radius;

        if(radiusAfterMassLoss > radiusBeforeMassLoss){
			if(m_aFastPhaseRocheLobeOverFlowErrorPrinting== false){
				m_aFastPhaseRocheLobeOverFlowErrorPrinting = true;
				std::cerr << m_randomSeed << "\tError in adaptive RL function, fast phase. Star radii increases after removing mass. Stellartype, Mass:\t" << starcopyDonorTemp.m_stellarType << ", " << massAfterMassLossTemp << std::endl;
			}
		}


            absoluteErrorTemp = fabs(radiusAfterMassLossTemp-RLfinalRsolTemp);
            relativeErrorTemp = fabs(radiusAfterMassLossTemp-RLfinalRsolTemp)/RLfinalRsolTemp;

            if(i==1){
                absoluteError = absoluteErrorTemp;
                relativeError = relativeErrorTemp;

				aChosen		= aChosenTemp;
                RLfinalRsol = RLfinalRsolTemp;
                massAfterMassLoss = massAfterMassLossTemp;
                radiusAfterMassLoss = radiusAfterMassLossTemp;
            }

//            if(absoluteErrorTemp <= absoluteErrorThreshold){
//                if(absoluteErrorTemp < absoluteError){
			if(relativeErrorTemp <= relativeErrorThreshold){
                        if(relativeErrorTemp < relativeError){
							absoluteError = absoluteErrorTemp;
							relativeError = relativeErrorTemp;
							
							aChosen		= aChosenTemp;
							RLfinalRsol = RLfinalRsolTemp;
							massAfterMassLoss = massAfterMassLossTemp;
							radiusAfterMassLoss = radiusAfterMassLossTemp;
                }
            }

        percentageMassLoss+=(maxPercentageMassLoss/niter);

        if(debugging){
			std::cout << "m', r', RL',a', age, tMS, tau, rET:\t" << massAfterMassLossTemp<< "\t" << radiusAfterMassLossTemp << "\t" <<RLfinalRsolTemp << "\t" << afinal*AUToRsol << "\t" << starcopyDonorTemp.m_Age << "\t" << starcopyDonorTemp.m_timescales[0] << "\t" << starcopyDonorTemp.m_Age/starcopyDonorTemp.m_timescales[0] << "\t" << relativeError  << std::endl;
        }
    }


    if(debugging){
			std::cout << "m', r', RL',a', age, tMS, tau, rET:\t" << massAfterMassLossTemp<< "\t" << radiusAfterMassLossTemp << "\t" <<RLfinalRsolTemp << "\t" << afinal << "\t" << starcopyDonorTemp.m_Age << "\t" << starcopyDonorTemp.m_timescales[0] << "\t" << starcopyDonorTemp.m_Age/starcopyDonorTemp.m_timescales[0] << "\t" << relativeError  << std::endl;    
	}

    // If mass and radius are both sensible, take the log of them
    if(radiusAfterMassLoss <= 0.0){std::cerr << m_randomSeed << "\tError in massTransferFastCaseA, radius after fake mass loss < 0.0" << std::endl;}
    if(massAfterMassLoss <= 0.0){std::cerr << m_randomSeed << "\tError in massTransferFastCaseA, mass after fake mass loss < 0.0" << std::endl;}
	// ALEJANDRO - 19/01/2017 - Following error message may not classify as an error per se, but want to check it in the error file to see if there is any correlation with the error in the slow phase.
	if(massAfterMassLoss == massBeforeMassLoss*maxPercentageMassLoss){std::cerr << m_randomSeed << "\tMaximum mass allowed being lost during fast phase, maxPercentageMassLoss = " << maxPercentageMassLoss << std::endl;}
	
    deltaM = massAfterMassLoss - massBeforeMassLoss;
    deltaR = radiusAfterMassLoss - radiusBeforeMassLoss;

    if(debugging){
        std::cout << "a, a' [AU]:\t" << a << " " << aChosen << std::endl;
        std::cout << "a, a' [Rsol]:\t" << a*AUToRsol << " " << aChosen*AUToRsol << std::endl;
        std::cout << "RL, RL' [AU]:\t" << RL << " " << RLfinalRsol*RsolToAU << std::endl;
        std::cout << "RL, RL' [Rsol]:\t" << RLRsol << " " << RLfinalRsol << std::endl;
        std::cout << "R, R' [Rsol]:\t" << radiusBeforeMassLoss << " " << radiusAfterMassLoss << std::endl;
        std::cout << "M, M' [Msol]:\t" << massBeforeMassLoss << " " << massAfterMassLoss << std::endl;
        std::cout << "dM [Msol]:\t" << deltaM << std::endl;
   }

    // Return mass
    return deltaM;
}


bool isMassRatioUnstable(int stellarTypeDonor, double massDonor, int stellarTypeAccretor, double massAccretor, const programOptions &options, unsigned long m_randomSeed){
    // Function for determining if mass transfer produces a wet merger according to the mass ratio limit discussed by de Mink et al. 2013 and Claeys et al. 2014
    bool    debugging = false;
    double  q = massAccretor/massDonor;

//    debugging = true;
    if(debugging){
        std::cout << "CriticalMassRatio function" << std::endl;
        std::cout << "Stellar type donor: " << stellarTypeDonor << std::endl;
        std::cout << "Stellar type accretor: " << stellarTypeAccretor << std::endl;
        std::cout << "Mass donor: " << massDonor << std::endl;
        std::cout << "Mass donor: " << massAccretor << std::endl;
        std::cout << "q: " << q << std::endl;
    }

    if(stellarTypeAccretor <= HELIUM_WHITE_DWARF){  // Non-degenerate accretors

        if(stellarTypeDonor == MS_LESS_THAN_07){

				if(debugging)
					std::cout << "q_crit:\t" << options.massTransferCriticalMassRatioMSLowMassNonDegenerateAccretor << std::endl;

            if(options.massTransferCriticalMassRatioMSLowMass){					
                if(q < options.massTransferCriticalMassRatioMSLowMassNonDegenerateAccretor){
                    return true;
                }
                else{
                    return false;
                }
            }
            else{
                return false;
            }
        }

        if(stellarTypeDonor == MS_MORE_THAN_07){
			if(debugging){std::cout << "q_crit:\t" << options.massTransferCriticalMassRatioMSHighMassNonDegenerateAccretor << std::endl;}
            
			if(options.massTransferCriticalMassRatioMSHighMass){
                if(q < options.massTransferCriticalMassRatioMSHighMassNonDegenerateAccretor){
                    return true;
                }
                else{
                    return false;
                }
            }
            else{
                return false;
            }
        }
		
        if(stellarTypeDonor == HERTZSPRUNG_GAP){
			if(debugging){std::cout << "q_crit:\t" << options.massTransferCriticalMassRatioHGNonDegenerateAccretor << std::endl;}

            if(options.massTransferCriticalMassRatioHG){
                if(q < options.massTransferCriticalMassRatioHGNonDegenerateAccretor){
                    return true;
                }
                else{
                    return false;
                }
            }
            else{
                return false;
            }
        }

        if(stellarTypeDonor >= FIRST_GIANT_BRANCH and stellarTypeDonor <= THERMALLY_PULSING_ASYMPTOTIC_GIANT_BRANCH){
			if(debugging){std::cout << "q_crit:\t" << options.massTransferCriticalMassRatioGiantNonDegenerateAccretor << std::endl;}
					
            if(options.massTransferCriticalMassRatioGiant){
                if(q <options.massTransferCriticalMassRatioGiantNonDegenerateAccretor){
                    return true;
                }
                else{
                    return false;
                }
            }
            else{
                return false;
            }
        }

        if(stellarTypeDonor == NAKED_HELIUM_STAR_MS){
			if(debugging){std::cout << "q_crit:\t" << options.massTransferCriticalMassRatioHeliumMSNonDegenerateAccretor << std::endl;}
            
			if(options.massTransferCriticalMassRatioHeliumMS){
                if(q < options.massTransferCriticalMassRatioHeliumMSNonDegenerateAccretor){
                    return true;
                }
                else{
                    return false;
                }
            }
            else{
                return false;
            }
        }
		
        if(stellarTypeDonor == NAKED_HELIUM_STAR_HETZSPRUNG_GAP){
			if(debugging){std::cout << "q_crit:\t" << options.massTransferCriticalMassRatioHeliumHGNonDegenerateAccretor << std::endl;}

            if(options.massTransferCriticalMassRatioHeliumHG){
                if(q < options.massTransferCriticalMassRatioHeliumHGNonDegenerateAccretor){
                    return true;
                }
                else{
                    return false;
                }
            }
            else{
                return false;
            }
        }
		
        if(stellarTypeDonor == NAKED_HELIUM_STAR_GIANT_BRANCH){
			if(debugging){std::cout << "q_crit:\t" << options.massTransferCriticalMassRatioHeliumGiantNonDegenerateAccretor << std::endl;}

            if(options.massTransferCriticalMassRatioHeliumGiant){
                if(q < options.massTransferCriticalMassRatioHeliumGiantNonDegenerateAccretor){
                    return true;
                }
                else{
                    return false;
                }
            }
            else{
                return false;
            }
        }


    }
    else if( stellarTypeAccretor <= BLACK_HOLE){    // Degenerate accretors


        if(stellarTypeDonor == MS_LESS_THAN_07){
			if(debugging){std::cout << "q_crit:\t" << options.massTransferCriticalMassRatioMSLowMassDegenerateAccretor << std::endl;}
			
            if(options.massTransferCriticalMassRatioMSLowMass){
                if(q < options.massTransferCriticalMassRatioMSLowMassDegenerateAccretor){
                    return true;
                }
                else{
                    return false;
                }
            }
            else{
                return false;
            }
        }
		
        if(stellarTypeDonor == MS_MORE_THAN_07){
			if(debugging){std::cout << "q_crit:\t" << options.massTransferCriticalMassRatioMSHighMassDegenerateAccretor << std::endl;}

            if(options.massTransferCriticalMassRatioMSHighMass){
                if(q < options.massTransferCriticalMassRatioMSHighMassDegenerateAccretor){
                    return true;
                }
                else{
                    return false;
                }
            }
            else{
                return false;
            }
        }

        if(stellarTypeDonor == HERTZSPRUNG_GAP){
			if(debugging){std::cout << "q_crit:\t" << options.massTransferCriticalMassRatioHGDegenerateAccretor << std::endl;}

            if(options.massTransferCriticalMassRatioHG){
                if(q < options.massTransferCriticalMassRatioHGDegenerateAccretor){
                    return true;
                }
                else{
                    return false;
                }
            }
            else{
                return false;
            }
        }

        if(stellarTypeDonor >= FIRST_GIANT_BRANCH and stellarTypeDonor <= THERMALLY_PULSING_ASYMPTOTIC_GIANT_BRANCH){
			if(debugging){std::cout << "q_crit:\t" << options.massTransferCriticalMassRatioGiantDegenerateAccretor << std::endl;}
			
            if(options.massTransferCriticalMassRatioGiant){
                if(q <options.massTransferCriticalMassRatioGiantDegenerateAccretor){
                    return true;
                }
                else{
                    return false;
                }
            }
            else{
                return false;
            }
        }

        if(stellarTypeDonor == NAKED_HELIUM_STAR_MS){
			if(debugging){std::cout << "q_crit:\t" << options.massTransferCriticalMassRatioHeliumMSDegenerateAccretor << std::endl;}
			
            if(options.massTransferCriticalMassRatioHeliumMS){

                if(q < options.massTransferCriticalMassRatioHeliumMSDegenerateAccretor){
                    return true;
                }
                else{
                    return false;
                }
            }
            else{
                return false;
            }
        }
		
        if(stellarTypeDonor == NAKED_HELIUM_STAR_HETZSPRUNG_GAP){
			if(debugging){std::cout << "q_crit:\t" << options.massTransferCriticalMassRatioHeliumHGDegenerateAccretor << std::endl;}
			
            if(options.massTransferCriticalMassRatioHeliumHG){
                if(q < options.massTransferCriticalMassRatioHeliumHGDegenerateAccretor){
                    return true;
                }
                else{
                    return false;
                }
            }
            else{
                return false;
            }
        }


        if(stellarTypeDonor == NAKED_HELIUM_STAR_GIANT_BRANCH){
			if(debugging){std::cout << "q_crit:\t" << options.massTransferCriticalMassRatioHeliumGiantDegenerateAccretor << std::endl;}
			
            if(options.massTransferCriticalMassRatioHeliumGiant){
                if(q < options.massTransferCriticalMassRatioHeliumGiantDegenerateAccretor){
                    return true;
                }
                else{
                    return false;
                }
            }
            else{
                return false;
            }
        }


    }
    else{
        if(debugging){
            std::cerr << m_randomSeed << "\tError in criticalMassRatio function. Shouldn't get here." << std::endl;
        }
        return false;
    }

    return false;
}
