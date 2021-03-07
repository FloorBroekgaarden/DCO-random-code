//
//  programOptionsSorter.cpp

#include "programOptionsSorter.h"

void splashScreen(){
    // Print a nice splash screen when the program starts
    std::cout << "\nCOMPAS: Compact Object Mergers: Population Astrophysics and Statistics \nby Team COMPAS (http://www.sr.bham.ac.uk/compas/)\nA binary star simulator\n" << std::endl;
}

void commandLineSorter(int argc, char * argv[], int & programStatus, programOptions & options){
    /*
     Function to read in and process the command line arguments using BOOST

     Parameters
     -----------
     argc : int
        argument count : argc will be the number of strings pointed to by argv. This will (in practice) be 1 plus the number of arguments, as virtually all implementations will prepend the name of the program to the array.
     argv : char array
        argument vector : array containing arguments
     programStatus : int
        Status of the program (whether to exit due to error, continue etc.)
     options : programOptions
        Class to contain the program options.
     */

	splashScreen();

    // Create an alias for the program options namespace and the filesystem namespace
    namespace po = boost::program_options;
    namespace fs = boost::filesystem;

    // Try rewriting the simplest version of this.
    try {

        std::string pathTemp = "fred"; // This temporary string is rewritten later

        // Create program options object
        po::options_description desc("Options");
        desc.add_options()
		    ("help,h", "Print this help message")
		    ("quiet", "Suppress printing")
		    ("random-seed", po::value<unsigned long>(&options.randomSeed), "Random seed to use (default = 0)")
		    ("debugging", "Include print statements to aid in debugging")
		    ("RejectedSamplesPrinting", "print the samples that are rejected by COMPAS in a seperate output file")
		    ("only-double-compact-objects", "Only evolve binaries which may form double compact objects")
			("evolve-unbound-systems", "Keeps evolving stars even if the binary is disrupted")
			("number-of-binaries,n", po::value<int>(&options.nBinaries), "Specify the number of binaries to simulate (default = 1000000)")
		    ("outputPath,o", po::value<std::string>(&pathTemp), "Directory for output (default = CWD)")
		    ("maximum-evolution-time", po::value<double>(&options.maxEvolutionTime), "Maximum time to evolve binaries in Myrs (default = 5.0)")
		    ("maximum-number-iterations", po::value<int>(&options.maxNumberOfTimestepIterations), "Maximum number of timesteps to evolve binary before giving up (default = 99999)")
		    ("detailedOutput", "Print detailed output to file")
            ("populationDataPrinting", "Print details of population")

            ("single-star", "Evolve single star(s)")
		    ("individual-system", "Run an individual system")
		    ("individual-initial-primary-mass", po::value<double>(&options.primaryMass), "Initial mass for primary in Msol (default = 20.0)")
		    ("individual-initial-secondary-mass", po::value<double>(&options.secondaryMass), "Initial mass for secondary in Msol (default = 10.0)")
		    ("individual-initial-primary-metallicity", po::value<double>(&options.initialPrimaryMetallicity), "Initial metallicity for primary (default = 0.02)")
		    ("individual-initial-secondary-metallicity", po::value<double>(&options.initialSecondaryMetallicity), "Initial metallicity for secondary (default = 0.02)")
		    ("individual-initial-primary-type", po::value<int>(&options.primaryStellarType), "Initial stellar type for primary (not yet implemented) (default ZAMS from mass)")
		    ("individual-initial-secondary-type", po::value<int>(&options.secondaryStellarType), "Initial stellar type for secondary (not yet implemented) (default ZAMS from mass)")
		    ("individual-initial-primary-rotational-velocity", po::value<double>(&options.primaryRotationalVelocity), "Initial rotational velocity for primary (not yet implemented) (default = 0.0)")
		    ("individual-initial-secondary-rotational-velocity", po::value<double>(&options.secondaryRotationalVelocity), "Initial rotational velocity for secondary (not yet implemented) (default = 0.0)")
		    ("individual-initial-primary-core-mass", po::value<double>(&options.primaryCoreMass), "Initial core mass for primary in Msol (default = 0.0)")
		    ("individual-initial-secondary-core-mass", po::value<double>(&options.secondaryCoreMass), "Initial core mass for secondary in Msol (default = 0.0)")
		    ("individual-effective-initial-primary-mass", po::value<double>(&options.primaryEffectiveInitialMass), "Effective initial mass for primary in Msol (default = Mass)")
		    ("individual-effective-initial-secondary-mass", po::value<double>(&options.secondaryEffectiveInitialMass), "Effective initial mass for secondary in Msol (default = Mass)")
		    ("individual-initial-primary-age", po::value<double>(&options.primaryAge), "Initial age for primary in Myrs (default = 0.0)")
		    ("individual-initial-secondary-age", po::value<double>(&options.secondaryAge), "Initial age for secondary in Myrs (default = 0.0)")
		    ("individual-initial-orbital-separation", po::value<double>(&options.binarySeparation), "Initial orbital separation in AU (default = 1.0)")
		    ("individual-initial-orbital-period", po::value<double>(&options.binaryOrbitalPeriod), "Initial orbital period in days (default from masses and separation)")
		    ("individual-initial-orbital-eccentricity", po::value<double>(&options.binaryEccentricity), "Initial orbital eccentricity (default = 0.0)")

		    ("kick-velocity-distribution", po::value<std::string>(&options.kickVelocityDistribution), "Natal kick velocity distribution (options: ZERO, FLAT, MAXWELLIAN, MUELLER2016, MUELLER2016MAXWELLIAN, BRAYELDRIDGE. Default = MAXWELLIAN)")

		    //("kick-velocity-sigma", po::value<double>(&options.kickVelocityDistributionSigma), "Sigma for chosen kick velocity distribution (default = 250 km s^-1 )")
		    //kickVelocityDistributionSigmaCCSN_NS
		    ("kick-velocity-sigma-CCSN-NS", po::value<double>(&options.kickVelocityDistributionSigmaCCSN_NS), "Sigma for chosen kick velocity distribution for neutron stars (default = 250 km s^-1 )")
		    ("kick-velocity-sigma-CCSN-BH", po::value<double>(&options.kickVelocityDistributionSigmaCCSN_BH), "Sigma for chosen kick velocity distribution for black holes (default = 250 km s^-1 )")
			("kick-velocity-sigma-ECSN", po::value<double>(&options.kickVelocityDistributionSigmaForECSN), "Sigma for chosen kick velocity distribution for ECSN (default = 0 km s^-1 )")
			("kick-velocity-sigma-USSN", po::value<double>(&options.kickVelocityDistributionSigmaForUSSN), "Sigma for chosen kick velocity distribution for USSN (default = 20 km s^-1 )")

		    ("kick-velocity-max", po::value<double>(&options.kickVelocityDistributionMaximum), "Maximum drawn kick velocity in km s^-1. Ignored if < 0. Must be > 0 if using kick-velocity-distribution=FLAT")
			("kick-scaling-factor", po::value<double>(&options.kickScalingFactor), "Arbitrary factor used to scale kicks (default = 1.0 )")
		    ("fix-dimensionless-kick-velocity", po::value<double>(&options.fixedUK), "Fix dimensionless kick velocity uk to this value (default = -1, -ve values false, +ve values true)")
		    ("kick-direction", po::value<std::string>(&options.kickDirection), "Distribution for natal kick direction (options: ISOTROPIC, INPLANE, PERPENDICULAR, POWERLAW. Default = ISOTROPIC)")
		    ("kick-direction-power", po::value<double>(&options.kickDirectionPower), "Power for power law kick direction distribution (default = 0.0 = isotropic, +ve = polar, -ve = in plane)")
		  	("black-hole-kicks", po::value<std::string>(&options.blackHoleKicks), "Black hole kicks relative to NS kicks (options: FULL, REDUCED, ZERO, FALLBACK. Default = FALLBACK)")
		    ("remnant-mass-prescription", po::value<std::string>(&options.remnantMassPrescription), "Choose remnant mass prescription (options: hurley2000, belczynski2002, fryer2012, muller2016, muller2016Maxwellian. Default = fryer2012)")
		    ("fryer-supernova-engine", po::value<std::string>(&options.fryerSupernovaEngineString), "If using Fryer et al 2012 fallback prescription, select between 'delayed' and 'rapid' engines (default = 'delayed')")
		    ("neutrino-mass-loss-bh-formation", po::value<std::string>(&options.neutrinoMassLossAssumptionBHString), "Assumption about neutrino mass loss during BH formation (options: FIXED_FRACTION, FIXED_MASS. default = 'FIXED_FRACTION')")
		    ("neutrino-mass-loss-bh-formation-value", po::value<double>(&options.neutrinoMassLossValueBH), "Value corresponding to neutrino mass loss assumption (default = 0.1)")

		    ("pair-instability-supernovae", "Enable pair instability supernovae (PISN)")
		    ("PISN-upper-limit", po::value<double>(&options.pairInstabilityUpperLimit), "Maximum core mass for PISN (default = 135.0)")
		    ("PISN-lower-limit", po::value<double>(&options.pairInstabilityLowerLimit), "Minimum core mass for PISN (default = 60.0)")

		    ("pulsational-pair-instability", "Enable mass loss due to pulsational-pair-instability (PPI)")
		    ("PPI-lower-limit", po::value<double>(&options.pulsationalPairInstabilityLowerLimit), "Minimum core mass for PPI (default = 35.0)")
		    ("PPI-upper-limit", po::value<double>(&options.pulsationalPairInstabilityUpperLimit), "Maximum core mass for PPI (default = 60.0)")

		    ("pulsational-pair-instability-prescription", po::value<std::string>(&options.pulsationalPairInstabilityPrescriptionString), "Specify which prescription to use for pulsational pair instability (options: COMPAS, STARTRACK, MARCHANT default = COMPAS)")

			("maximum-neutron-star-mass", po::value<double>(&options.maximumNeutronStarMass), "Maximum mass of a neutron star (default = 3.0, as in StarTrack)")
			("nbatches-used", po::value<int>(&options.nBatchesUsed), "nr of batches used, only needed for STROOPWAFEL (AIS) (default = -1, not needed) sets itself automatically in pythonSubmit")

		    ("initial-mass-function,i", po::value<std::string>(&options.initialMassFunction), "Specify initial mass function to use (options: SALPETER, POWERLAW, UNIFORM, KROUPA. default = KROUPA)")
		    ("initial-mass-min", po::value<double>(&options.initialMassFunctionMin), "Minimum mass (in Msol) to generate using given IMF (default = 8)")
		    ("initial-mass-max", po::value<double>(&options.initialMassFunctionMax), "Maximum mass (in Msol) to generate using given IMF (default = 100)")
		    ("initial-mass-power", po::value<double>(&options.initialMassFunctionPower), "Single power law power to generate primary mass using given IMF")

		    ("spin-assumption", po::value<std::string>(&options.spinAssumption), "Which assumption of misalignedments to use (default = bothAligned)")
		    ("spin-distribution", po::value<std::string>(&options.spinDistribution), "Which distribution of spins to use (default = 0)")
		    ("spin-mag-min", po::value<double>(&options.spinDistributionMin), "Minimum magnitude of spin (default = )")
		    ("spin-mag-max", po::value<double>(&options.spinDistributionMax), "Maximum magnitude of spin (default = )")

		    ("semi-major-axis-distribution,a", po::value<std::string>(&options.semiMajorAxisDistribution), "Initial semi-major axis distribution, a (options: FLATINLOG, DuquennoyMayor1991, SANA2012. default = FLATINLOG)")
		    ("semi-major-axis-min", po::value<double>(&options.semiMajorAxisDistributionMin), "Minimum semi major axis in AU to generate (default = 0.1)")
		    ("semi-major-axis-max", po::value<double>(&options.semiMajorAxisDistributionMax), "Maximum semi major axis in AU to generate (default = 1000)")
		   	("orbital-period-min", po::value<double>(&options.periodDistributionMin), "Minimum period in days to generate (default = 1.1)")
		    ("orbital-period-max", po::value<double>(&options.periodDistributionMax), "Maximum period in days to generate (default = 1000)")

		    ("eccentricity-distribution,e", po::value<std::string>(&options.eccentricityDistribution), "Initial eccentricity distribution, e (options: ZERO, FIXED, FLAT, THERMALISED, GELLER+2013. Default = ZERO)")
		    ("eccentricity-min", po::value<double>(&options.eccentricityDistributionMin), "Minimum eccentricity to generate (default = 0.0)")
		    ("eccentricity-max", po::value<double>(&options.eccentricityDistributionMax), "Maximum eccentricity to generate (default = 1.0)")

		    ("mass-ratio-distribution,q", po::value<std::string>(&options.massRatioDistribution), "Initial mass ratio distribution for q=m2/m1 (options: FLAT, DuquennoyMayor1991, SANA2012. Default = FLAT)")
		    ("mass-ratio-min", po::value<double>(&options.massRatioDistributionMin), "Minimum mass ratio m2/m1 to generate (default = 0.0)")
		    ("mass-ratio-max", po::value<double>(&options.massRatioDistributionMax), "Maximum mass ratio m2/m1 to generate (default = 1.0)")

		    ("minimum-secondary-mass", po::value<double>(&options.minimumMassSecondary), "Minimum mass of secondary to generate in Msol (default = 0.0)")

		    ("rotational-velocity-distribution", po::value<std::string>(&options.rotationalVelocityString), "Initial rotational velocity distribution (options: ZERO, HURLEY, VLTFLAMES. Default = ZERO)")

		    ("pulsar-birth-magnetic-field-distribution", po::value<std::string>(&options.pulsarBirthMagneticFieldDistributionString), "Distribution of (log10 of) pulsar birth magnetic field in G (options: ZERO, FIXED, FLATINLOG, UNIFORM, LOGNORMAL. Default = ZERO)")
		    ("pulsar-birth-magnetic-field-distribution-min", po::value<double>(&options.pulsarBirthMagneticFieldDistributionMin), "Minimum (log10) pulsar birth magnetic field) (default = 11.0)")
		    ("pulsar-birth-magnetic-field-distribution-max", po::value<double>(&options.pulsarBirthMagneticFieldDistributionMax), "Maximum (log10) pulsar birth magnetic field (default = 13.0)")

		    ("pulsar-birth-spin-period-distribution", po::value<std::string>(&options.pulsarBirthSpinPeriodDistributionString), "Distribution of pulsar birth spin period in ms (options: ZERO, FIXED, UNIFORM, NORMAL. Default = ZERO)")
		    ("pulsar-birth-spin-period-distribution-min", po::value<double>(&options.pulsarBirthSpinPeriodDistributionMin), "Minimum pulsar birth spin period in ms (default = 0.0)")
		    ("pulsar-birth-spin-period-distribution-max", po::value<double>(&options.pulsarBirthSpinPeriodDistributionMax), "Maximum pulsar birth spin period in ms (default = 100.0)")

		    ("pulsar-magnetic-field-decay-timescale", po::value<double>(&options.pulsarMagneticFieldDecayTimescale), "Timescale on which magnetic field decays in Myrs (default = 1000.0)")
		    ("pulsar-magnetic-field-decay-massscale", po::value<double>(&options.pulsarMagneticFieldDecayMassscale), "Mass scale on which magnetic field decays during accretion in solar masses (default = 0.025)")
		    ("pulsar-minimum-magnetic-field", po::value<double>(&options.pulsarLog10MinimumMagneticField), "log10 of the minimum pulsar magnetic field in Gauss (default = 8.0)")

		    ("evolve-pulsars", "Whether to evolve pulsars")
		    
		    // NS EOS parameters
		    ("neutron-star-equation-of-state", po::value<std::string>(&options.neutronStarEquationOfStateString), "Specify which neutron star equation of state to use (options: SSE, ARP3 default = SSE)")

			("PN", "Enable post-newtonian evolution")
		    ("use-mass-loss", "Enable mass loss")
		    ("tides-prescription", po::value<std::string>(&options.tidesPrescriptionString), "Tides prescription to use (options: default = None)")
		    ("mass-loss-prescription", po::value<std::string>(&options.massLossPrescriptionString), "Mass loss prescription to use (options: NONE, HURLEY, VINK. Default = NONE)")
		    ("luminous-blue-variable-multiplier", po::value<double>(&options.luminousBlueVariableFactor), "Multiplicitive constant for LBV mass loss (default = 1.5, use 10 for Mennekens & Vanbeveren 2014)")
		    ("wolf-rayet-multiplier", po::value<double>(&options.wolfRayetFactor), "Multiplicitive constant for WR winds (default = 1.0)")

		   	("massTransfer", "Enable mass transfer")
			("circulariseBinaryDuringMassTransfer", "Circularise binary when it enters a Mass Transfer episode (default = False)")
			("forceCaseBBBCStabilityFlag", "Force case BB/BC mass transfer to be only stable or unstable (default = True)")
			("alwaysStableCaseBBBCFlag", "Choose case BB/BC mass transfer to be always stable (default = True)")
			("angularMomentumConservationDuringCircularisation", "Conserve angular momentum when binary is circularised when entering a Mass Transfer episode (default = False)")
		    ("mass-transfer-prescription", po::value<std::string>(&options.massTransferPrescriptionString), "Mass Transfer prescription to use (default = DEMINK)")
		    ("mass-transfer-angular-momentum-loss-prescription", po::value<std::string>(&options.massTransferAngularMomentumLossPrescriptionString), "Mass Transfer Angular Momentum Loss prescription to use (options: JEANS, ISOTROPIC, CIRCUMBINARY, ARBITRARY. Default = ISOTROPIC)")
		    ("mass-transfer-accretion-efficiency-prescription", po::value<std::string>(&options.massTransferAccretionEfficiencyPrescriptionString), "Mass Transfer Accretion Efficiency prescription to use (options: THERMAL, FIXED, CENTRIFUGAL. Default = THERMAL)")
			("mass-transfer-thermal-limit-accretor", po::value<std::string>(&options.massTransferThermallyLimitedVariationString), "Mass Transfer Thermal Accretion limit to use (default = CFACTOR)")
		    ("mass-transfer-rejuvenation-prescription", po::value<std::string>(&options.massTransferRejuvenationPrescriptionString), "Mass Transfer Rejuvenation prescription to use (options: NONE, STARTRACK. Default = NONE)")
		    ("mass-transfer-fa", po::value<double>(&options.massTransferFractionAccreted), "Mass Transfer fraction accreted (default = 1.0, fully conservative)")
		    ("mass-transfer-jloss", po::value<double>(&options.massTransferJloss), "Specific angular momentum with which the non-accreted system leaves the system (default = 1.0)")
			("mass-transfer-thermal-limit-C", po::value<double>(&options.massTransferCParameter), "Mass Transfer Thermal rate factor fo the accretor (default = 10.0, Hurley+2002)")
			("eddington-accretion-factor", po::value<double>(&options.eddingtonAccretionFactor), "Multiplication factor for eddington accretion for NS & BH, i.e. >1 is super-eddington and 0. is no accretion")

            ("critical-mass-ratio-MS-low-mass-non-degenerate-accretor", po::value<double>(&options.massTransferCriticalMassRatioMSLowMassNonDegenerateAccretor), "Critical mass ratio for MT from a MS star (default = 1.44, Claeys+ 2014). Specify both MS low mass flags to use. 0 is always stable, <0 is disabled")
            ("critical-mass-ratio-MS-low-mass-degenerate-accretor", po::value<double>(&options.massTransferCriticalMassRatioMSLowMassDegenerateAccretor), "Critical mass ratio for MT from a MS star to a degenerate accretor (default = 1.0 from Claeys+ 2014) Specify both MS low mass flags to use. 0 is always stable, <0 is disabled")
            ("critical-mass-ratio-MS-high-mass-non-degenerate-accretor", po::value<double>(&options.massTransferCriticalMassRatioMSHighMassNonDegenerateAccretor), "Critical mass ratio for MT from a MS star (default = 0.625, Claeys+ 2014). Specify both MS high mass flags to use. 0 is always stable, <0 is disabled")
            ("critical-mass-ratio-MS-high-mass-degenerate-accretor", po::value<double>(&options.massTransferCriticalMassRatioMSHighMassDegenerateAccretor), "Critical mass ratio for MT from a MS star to a degenerate accretor (default = 0.0 from Claeys+ 2014) Specify both MS high mass flags to use. 0 is always stable, <0 is disabled")			
            ("critical-mass-ratio-HG-non-degenerate-accretor", po::value<double>(&options.massTransferCriticalMassRatioHGNonDegenerateAccretor), "Critical mass ratio for MT from a HG star (default = 0.40 from de Mink+ 2013) Specify both HG flags to use. 0 is always stable, <0 is disabled")
            ("critical-mass-ratio-HG-degenerate-accretor", po::value<double>(&options.massTransferCriticalMassRatioHGDegenerateAccretor), "Critical mass ratio for MT from a HG star (default = 0.21 from Claeys+ 2014) Specify both HG flags to use. 0 is always stable, <0 is disabled")
            ("critical-mass-ratio-giant-non-degenerate-accretor", po::value<double>(&options.massTransferCriticalMassRatioGiantNonDegenerateAccretor), "Critical mass ratio for MT from a giant star (default = not implemented from Claeys+ 2014) Specify both giant flags to use. 0 is always stable, <0 is disabled")
            ("critical-mass-ratio-giant-degenerate-accretor", po::value<double>(&options.massTransferCriticalMassRatioGiantDegenerateAccretor), "Critical mass ratio for MT from a giant star (default = 0.87 from Claeys+ 2014) Specify both giant flags to use. 0 is always stable, <0 is disabled")
            ("critical-mass-ratio-helium-MS-non-degenerate-accretor", po::value<double>(&options.massTransferCriticalMassRatioHeliumMSNonDegenerateAccretor), "Critical mass ratio for MT from a helium MS star (default = 0.625) Specify both helium MS flags to use. 0 is always stable, <0 is disabled")
            ("critical-mass-ratio-helium-MS-degenerate-accretor", po::value<double>(&options.massTransferCriticalMassRatioHeliumMSDegenerateAccretor), "Critical mass ratio for MT from a helium MS star (default = 0.0) Specify both helium MS flags to use. 0 is always stable, <0 is disabled")
            ("critical-mass-ratio-helium-HG-non-degenerate-accretor", po::value<double>(&options.massTransferCriticalMassRatioHeliumHGNonDegenerateAccretor), "Critical mass ratio for MT from a helium HG star (default = 0.25 from Claeys+ 2014) Specify both helium HG flags to use. 0 is always stable, <0 is disabled")
            ("critical-mass-ratio-helium-HG-degenerate-accretor", po::value<double>(&options.massTransferCriticalMassRatioHeliumHGDegenerateAccretor), "Critical mass ratio for MT from a helium HG star (default = 0.21 from Claeys+ 2014) Specify both helium HG flags to use. 0 is always stable, <0 is disabled")
            ("critical-mass-ratio-helium-giant-non-degenerate-accretor", po::value<double>(&options.massTransferCriticalMassRatioHeliumGiantNonDegenerateAccretor), "Critical mass ratio for MT from a helium giant star (default = 1.28 from Claeys+ 2014) Specify both helium giant flags to use. 0 is always stable, <0 is disabled")
            ("critical-mass-ratio-helium-giant-degenerate-accretor", po::value<double>(&options.massTransferCriticalMassRatioHeliumGiantDegenerateAccretor), "Critical mass ratio for MT from a helium giant star (default = 0.87 from Claeys+ 2014) Specify both helium giant flags to use. 0 is always stable, <0 is disabled")
            ("critical-mass-ratio-white-dwarf-non-degenerate-accretor", po::value<double>(&options.massTransferCriticalMassRatioWhiteDwarfNonDegenerateAccretor), "Critical mass ratio for MT from a white dwarf (default = 0.0) Specify both white dwarf flags to use. 0 is always stable, <0 is disabled")
            ("critical-mass-ratio-white-dwarf-degenerate-accretor", po::value<double>(&options.massTransferCriticalMassRatioWhiteDwarfDegenerateAccretor), "Critical mass ratio for MT from a white dwarf (default = 1.6 from Claeys+ 2014) Specify both white dwarf flags to use. 0 is always stable, <0 is disabled")

            ("metallicity,z", po::value<double>(&options.metallicity), "Metallicity to use (default 0.02 is Zsol)")

		    ("common-envelope-alpha", po::value<double>(&options.commonEnvelopeAlpha), "Common Envelope efficiency alpha (default = 1.0)")
		    ("common-envelope-lambda", po::value<double>(&options.commonEnvelopeLambda), "Common Envelope lambda (default = 0.1)")
			("common-envelope-hertzsprung-gap-assumption", po::value<std::string>(&options.commonEnvelopeHertzsprungGapDonorString), "Assumption to make about HG stars in CE (default = OPTIMISTIC_HG_CE)")
			("common-envelope-lambda-prescription", po::value<std::string>(&options.commonEnvelopeLambdaPrescriptionString), "Prescription for CE lambda (options: LAMBDA_FIXED, LAMBDA_LOVERIDGE, LAMBDA_NANJING, LAMBDA_KRUCKOW, LAMBDA_DEWI. Default = LAMBDA_FIXED)")
		    ("common-envelope-slope-Kruckow", po::value<double>(&options.commonEnvelopeSlopeKruckow), "Common Envelope slope for Kruckow lambda (default = -4/5)")
		    ("common-envelope-alpha-thermal", po::value<double>(&options.commonEnvelopeAlphaThermal), "Defined such that lambda = alpha_th * lambda_b + (1.0 - alpha_th) * lambda_g (default = 1.0)")
		    ("common-envelope-lambda-multiplier", po::value<double>(&options.commonEnvelopeLambdaMultiplier), "Multiply lambda by some constant (default = 1.0)")
		    ("common-envelope-allow-main-sequence-survive", "Allow main sequence stars to survive common envelope evolution")

		    ("common-envelope-mass-accretion-prescription", po::value<std::string>(&options.commonEnvelopeMassAccretionPrescriptionString), "Assumption about whether NS/BHs can accrete mass during common envelope evolution (options: ZERO, CONSTANT, UNIFORM, MACLEOD+2014 . Default = ZERO)")
		    ("common-envelope-mass-accretion-min", po::value<double>(&options.commonEnvelopeMassAccretionMin), "Minimum amount of mass accreted by NS/BHs during common envelope evolution in solar masses (default = 0.04)")
		    ("common-envelope-mass-accretion-max", po::value<double>(&options.commonEnvelopeMassAccretionMax), "Maximum amount of mass accreted by NS/BHs during common envelope evolution in solar masses (default = 0.1)")
	    	("common-envelope-mass-accretion-constant", po::value<double>(&options.commonEnvelopeMassAccretionConstant), "Value of mass accreted by NS/BH during common envelope evolution if assuming all NS/BH accrete same amount of mass (common-envelope-mass-accretion-prescription CONSTANT). Ignored otherwise (default = 0)")

			("lambda-calculation-every-timeStep", "Calculate all values of lambda at each timestep")
			
			("common-envelope-zeta-prescription", po::value<std::string>(&options.commonEnvelopeZetaPrescriptionString), "Prescription for CE zeta (default = STARTRACK)")
		    ("zeta-adiabatic-arbitrary", po::value<double>(&options.zetaAdiabaticArbitrary), "Value of mass-radius exponent zeta adiabatic")
		    ("zeta-main-sequence", po::value<double>(&options.zetaMainSequence), "Value of mass-radius exponent zeta on the main sequence (default = 2.0)")
		    ("zeta-hertzsprung-gap", po::value<double>(&options.zetaHertzsprungGap), "Value of mass-radius exponent zeta on the hertzstrpung gap (default = 6.5)")
			("zeta-Calculation-Every-Time-Step", "Enable calculate MT zetas every timestep")
						
			("revised-energy-formalism-Nandez-Ivanova", "Enable revised energy formalism")
		    ("maximum-mass-donor-Nandez-Ivanova", po::value<double>(&options.maximumMassDonorNandezIvanova), "Maximum donor mass allowed for the revised common envelope formalism in Msol (default = 2.0)")
		    ("common-envelope-recombination-energy-density", po::value<double>(&options.commonEnvelopeRecombinationEnergyDensity), "Recombination energy density in ergs/g (default = 1.5x10^13)")

		   	
			
			("BeBinaries", "Enable Be Binaries study")
            ("CHEvolution","Enable output Chemically Homogeneous Evolution")

            ("AIS-exploratory-phase", "Run exploratory phase of STROOPWAFEL") // Floor 
            ("AIS-DCOtype", po::value<std::string>(&options.AISDCOtypeString), "DCO type selection in exploratory phase of STROOPWAFEL, select ALL, BBH, BNS or BHNS")
		    ("AIS-Hubble", "Excluding not in Hubble time mergers selection in exploratory phase of STROOPWAFEL")
		    ("AIS-RLOF", "RLOFSecondaryZAMS selection in exploratory phase of STROOPWAFEL")
		    ("AIS-Pessimistic", "Optimistic or Pessimistic selection in exploratory phase of STROOPWAFEL")
		    ("AIS-refinement-phase", "If true: run main sampling phase (step2) of STROOPWAFEL")
		    ("kappa-gaussians", po::value<double>(&options.KappaGaussians), "Scaling factor for the width of the Gaussian distributions in STROOPWAFEL main sampling phase" )

            ("RLOFPrinting","Enable output parameters before/after RLOF ")

			("sample-kick-velocity-sigma","Sample over Kick Velocity Sigma")
			("sample-kick-velocity-sigma-min",po::value<double>(&options.sampleKickVelocitySigmaMin),"Minimum for Uniform sampling over kick velocity sigma")
			("sample-kick-velocity-sigma-max",po::value<double>(&options.sampleKickVelocitySigmaMax),"Maximum for Uniform sampling over kick velocity sigma")

			("sample-kick-direction-power","Sample over kick direction powerlaw exponent")
			("sample-kick-direction-power-min",po::value<double>(&options.sampleKickDirectionPowerMin),"Minimum for Uniform sampling over kick direction powerlaw exponent")
			("sample-kick-direction-power-max",po::value<double>(&options.sampleKickDirectionPowerMax),"Maximum for Uniform sampling over kick direction powerlaw exponent")

			("sample-common-envelope-alpha","Sample over common envelope alpha")
			("sample-common-envelope-alpha-min",po::value<double>(&options.sampleCommonEnvelopeAlphaMin),"Minimum for Uniform sampling over common envelope alpha")
			("sample-common-envelope-alpha-max",po::value<double>(&options.sampleCommonEnvelopeAlphaMax),"Maximum for Uniform sampling over common envelope alpha")

			("sample-wolf-rayet-multiplier","Sample over WR winds multiplicative constant")
			("sample-wolf-rayet-multiplier-min",po::value<double>(&options.sampleWolfRayetMultiplierMin),"Minimum for Uniform sampling over multiplicative constant for WR winds")
			("sample-wolf-rayet-multiplier-max",po::value<double>(&options.sampleWolfRayetMultiplierMax),"Maximum for Uniform sampling over multiplicative constant for WR winds")

			("sample-luminous-blue-variable-multiplier","Sample over multiplicative constant from LBV mass loss")
			("sample-luminous-blue-variable-multiplier-min",po::value<double>(&options.sampleLuminousBlueVariableMultiplierMin),"Minimum for Uniform sampling over multiplicative constant for LBV mass loss")
			("sample-luminous-blue-variable-multiplier-max",po::value<double>(&options.sampleLuminousBlueVariableMultiplierMax),"Maximum for Uniform sampling over multiplicative constant for LBV mass loss")

			//("importance-sampling", "Use importance sampling")
			("mcmc","Use MCMC sampling (Not yet implemented. default = false)");
		;

        // Variables map
        po::variables_map vm;

        try {

            po::store(po::parse_command_line(argc, argv, desc), vm); // Can throw

            // --help option
            if (vm.count("help")){
                std::cout << "\nCOMPAS: Compact Object Mergers: Population Astrophysics and Statistics \nby Team COMPAS (http://www.sr.bham.ac.uk/compas/)\nA binary star simulator\n" << std::endl << desc << std::endl;
                programStatus = SUCCESS;
            }
            else{
                programStatus = CONTINUE;
            }

            // Invoke notify to assign user-input values to variables
            po::notify(vm); // Throws an error, so do after help just in case there are any problems.

            // Should do type/constraint checking on each variable as you enter it. e.g. negative number of binaries
            if(options.nBinaries <= 0){
                std::cout << "Asked for zero or a negative number of binaries - exiting" << std::endl;
                programStatus = ERROR_IN_COMMAND_LINE;
            }

            // Check if the output option was specified, check if user specified a valid path and assign options variable
            if (vm.count("outputPath") or vm.count("o")){

            	//std::cout << "Got here" << std::endl;

	            //std::cout << "test pathTemp" << pathTemp << std::endl;

                fs::path pathTempPath;              // Generate a temporary path from

                pathTempPath = pathTemp;           // If user specified a path, set temp path equal to that path.

                //std::cout << "pathTempPath " << pathTempPath << std::endl;

                // Now check that the path you just constructed is valid: if so, assign to the options path variable. Otherwise do nothing. (or cout/cerr an error saying it wasn't a valid path, using default path.)

                if(fs::is_directory(pathTempPath)){
                    options.outputPath = pathTempPath;
                }
                else{
                    std::cerr << "You did not enter a valid path. Using default (CWD)." << std::endl;
                    options.outputPath = options.defaultOutputPath;
                }

            }

            if(vm.count("random-seed")){
                options.fixedRandomSeed = true;
            }

            if (vm.count("quiet"))
            {
                options.quiet = true;
            }

            if(vm.count("debugging")){
                options.debugging = true;
            }

            if(vm.count("only-double-compact-objects")){
                options.onlyDoubleCompactObjects = true;
            }

            if(vm.count("evolve-unbound-systems")){
                options.evolveUnboundSystems = true;
            }

            if(vm.count("common-envelope-allow-main-sequence-survive")){
            	options.allowMainSequenceStarToSurviveCommonEnvelope = true;
            }

            // Check if user want to fix the dimensionless kick velocity
            if(vm.count("fix-dimensionless-kick-velocity")){
                if(options.fixedUK >= 0.0){
                    options.useFixedUK = true;
                }
            }

		    if(options.initialMassFunctionMin < 0.0){

		    	std::cerr << "User specified --initial-mass-min < 0" << std::endl;
		    	programStatus = ERROR_IN_COMMAND_LINE;

		    }

		    if(options.initialMassFunctionMax < 0.0){

		    	std::cerr << "User specified --initial-mass-max < 0" << std::endl;
		    	programStatus = ERROR_IN_COMMAND_LINE;

		    }

		    if(options.initialMassFunctionMax < options.initialMassFunctionMin){

		    	std::cerr << "User specified --initial-mass-max < --initial-mass-min" << std::endl;
		    	programStatus = ERROR_IN_COMMAND_LINE;

		    }

            // Check that user specified a maximum kick velocity if using the FLAT distribution
            if(boost::iequals(options.kickVelocityDistribution, "Flat") and (options.kickVelocityDistributionMaximum <= 0.0)){

            	std::cerr << "User specified --kick-velocity-distribution=Flat but did not specify --kick-velocity-max > 0. Exiting." << std::endl;
            	programStatus = ERROR_IN_COMMAND_LINE;

            }

            if(options.semiMajorAxisDistributionMin < 0.0){

            	std::cerr << "User specified --semi-major-axis-min < 0" << std::endl;
            	programStatus = ERROR_IN_COMMAND_LINE;

            }

            if(options.semiMajorAxisDistributionMax < 0.0){

            	std::cerr << "User specified --semi-major-axis-max < 0" << std::endl;
            	programStatus = ERROR_IN_COMMAND_LINE;

            }

            if(options.periodDistributionMin < 0.0){

            	std::cerr << "User specified --orbital-period-min < 0" << std::endl;
            	programStatus = ERROR_IN_COMMAND_LINE;

            }

            if(options.periodDistributionMax < 0.0){

            	std::cerr << "User specified --orbital-period-max < 0" << std::endl;
            	programStatus = ERROR_IN_COMMAND_LINE;

            }

            if(options.eccentricityDistributionMin < 0.0 or options.eccentricityDistributionMin > 1.0){

            	std::cerr << "--eccentricity-min must be between 0 and 1" << std::endl;
            	programStatus = ERROR_IN_COMMAND_LINE;

            } 

            if(options.eccentricityDistributionMax < 0.0 or options.eccentricityDistributionMax > 1.0){

            	std::cerr << "--eccentricity-max must be between 0 and 1" << std::endl;
            	programStatus = ERROR_IN_COMMAND_LINE;

            } 

            if(options.eccentricityDistributionMax < options.eccentricityDistributionMin){

            	std::cerr << "--eccentricity-max must be > --eccentricity-min" << std::endl;
            	programStatus = ERROR_IN_COMMAND_LINE;

            }

            if(options.massRatioDistributionMin < 0.0 or options.massRatioDistributionMin > 1.0){

            	std::cerr << "--mass-ratio-min must be between 0 and 1, was " << options.massRatioDistributionMin << std::endl;
            	programStatus = ERROR_IN_COMMAND_LINE;

            }

            if(options.massRatioDistributionMax < 0.0 or options.massRatioDistributionMax > 1.0){

            	std::cerr << "--mass-ratio-max must be between 0 and 1, was " << options.massRatioDistributionMax << std::endl;
            	programStatus = ERROR_IN_COMMAND_LINE;

            }

            if(options.massRatioDistributionMax < options.massRatioDistributionMin){

            	std::cerr << "--mass-ratio-max must be > --mass-ratio-min" << std::endl; 
            	programStatus = ERROR_IN_COMMAND_LINE;

            }

            if(options.minimumMassSecondary < 0.0){
            	std::cerr << "--minimum-secondary-mass must be >= 0, was " << options.minimumMassSecondary << " . Exiting." << std::endl;
            	programStatus = ERROR_IN_COMMAND_LINE; 
            }

            if(options.minimumMassSecondary > options.initialMassFunctionMax){
            	std::cerr << "--minimum-secondary-mass must be < --initial-mass-max, was " << options.minimumMassSecondary << " . Exiting." << std::endl;
            	programStatus = ERROR_IN_COMMAND_LINE; 
            }

            // Check user selected black hole kicks option
            if(vm.count("black-hole-kicks")){

                if(boost::iequals(options.blackHoleKicks, "FULL")){
                    options.blackHoleKicksOption = FULL;
                }
                else if(boost::iequals(options.blackHoleKicks, "REDUCED")){
                    options.blackHoleKicksOption = REDUCED;
                }
                else if(boost::iequals(options.blackHoleKicks, "ZERO")){
                    options.blackHoleKicksOption = ZERO;
                }
                else if(boost::iequals(options.blackHoleKicks, "FALLBACK")){
                    options.blackHoleKicksOption = FALLBACK;
                }
                else{
                    std::cerr << "Invalid choice for black hole kicks." << std::endl;
                    options.blackHoleKicksOption = FULL;
                }

            }

            // If user specifies which Fryer et al 2012 supernova engine to use, set the relevant variable.
            if(vm.count("fryer-supernova-engine")){
                if(boost::iequals(options.fryerSupernovaEngineString, "DELAYED")){
                    options.fryerSupernovaEngine = SN_DELAYED;
                }
                else if(boost::iequals(options.fryerSupernovaEngineString, "RAPID")){
                    options.fryerSupernovaEngine = SN_RAPID;
                }
                else{
                    std::cout << "Incorrect Fryer et al supernova engine. Using default DELAYED." << std::endl;
                    options.fryerSupernovaEngine = SN_DELAYED;
                    // Could also just exit the program here
                    // programStatus = ERROR_IN_COMMAND_LINE;
                }
            }

            // Set neutrino mass loss assumption variables
            if(vm.count("neutrino-mass-loss-bh-formation")){
                if(boost::iequals(options.neutrinoMassLossAssumptionBHString, "FIXED_FRACTION")){
                    options.neutrinoMassLossAssumptionBH = NEUTRINO_MASS_LOSS_BH_FIXED_FRACTION;
                }
                else if(boost::iequals(options.neutrinoMassLossAssumptionBHString, "FIXED_MASS")){
                    options.neutrinoMassLossAssumptionBH = NEUTRINO_MASS_LOSS_BH_FIXED_MASS;
                }
                else{
                    std::cout << "Incorrect neutrino mass loss assumption. Using default FIXED_FRACTION." << std::endl;
                    options.neutrinoMassLossAssumptionBH = NEUTRINO_MASS_LOSS_BH_FIXED_FRACTION;
                    // Could also just exit the program here
                    // programStatus = ERROR_IN_COMMAND_LINE;
                }
            }

            if(options.neutrinoMassLossAssumptionBH == NEUTRINO_MASS_LOSS_BH_FIXED_MASS){
            	if(options.neutrinoMassLossValueBH < 0.0){
            		std::cout << "Invalid neutrino mass loss value " << options.neutrinoMassLossValueBH << " requested. Must be > 0.0. Exiting" << std::endl;
                    programStatus = ERROR_IN_COMMAND_LINE;
            	}
            }

            if(options.neutrinoMassLossAssumptionBH == NEUTRINO_MASS_LOSS_BH_FIXED_FRACTION){
            	if(options.neutrinoMassLossValueBH < 0.0 or options.neutrinoMassLossValueBH > 1.0){
            		std::cout << "Invalid neutrino mass loss value " << options.neutrinoMassLossValueBH << " requested. Must be between 0 and 1. Exiting" << std::endl;
                    programStatus = ERROR_IN_COMMAND_LINE;
            	}
            }
            
            // Specify which tides prescription to use
            if(vm.count("tides-prescription")){
                if(boost::iequals(options.tidesPrescriptionString, "NONE")){
                    options.tidesPrescription = TIDES_PRESCRIPTION_NONE;
                }
                else if(boost::iequals(options.tidesPrescriptionString, "LOCKED_ENERGY")){
                    options.tidesPrescription = TIDES_PRESCRIPTION_LOCKED_ENERGY;
                }
                else if(boost::iequals(options.tidesPrescriptionString, "LOCKED_ANGMOMENTUM")){
                    options.tidesPrescription = TIDES_PRESCRIPTION_LOCKED_ANG_MOMENTUM;
                }
                else if(boost::iequals(options.tidesPrescriptionString, "HUT")){
                    options.tidesPrescription = TIDES_PRESCRIPTION_HUT;
                }
                else{
                    std::cout << "Invalid tides prescription requested. Exiting" << std::endl;
                    programStatus = ERROR_IN_COMMAND_LINE;
                }
            }

            // Check if user asked for PN evolution to be on
            if(vm.count("PN")){
                options.PNevolution = true;
            }

            if(vm.count("massTransfer")){
                options.useMassTransfer = true;
            }

            if(vm.count("critical-mass-ratio-MS-low-mass-non-degenerate-accretor") and vm.count("critical-mass-ratio-MS-low-mass-degenerate-accretor") and (options.massTransferCriticalMassRatioMSLowMassNonDegenerateAccretor >= 0.0) and (options.massTransferCriticalMassRatioMSLowMassDegenerateAccretor >= 0.0)){
                options.massTransferCriticalMassRatioMSLowMass = true;
            }
            else{
                options.massTransferCriticalMassRatioMSLowMass = false;
			}

            if(vm.count("critical-mass-ratio-MS-high-mass-non-degenerate-accretor") and vm.count("critical-mass-ratio-MS-high-mass-degenerate-accretor") and (options.massTransferCriticalMassRatioMSHighMassNonDegenerateAccretor >= 0.0) and (options.massTransferCriticalMassRatioMSHighMassDegenerateAccretor >= 0.0)){
                options.massTransferCriticalMassRatioMSHighMass = true;
            }
            else{
                options.massTransferCriticalMassRatioMSHighMass = false;
			}
	
            if(vm.count("critical-mass-ratio-HG-non-degenerate-accretor") and vm.count("critical-mass-ratio-HG-degenerate-accretor") and (options.massTransferCriticalMassRatioHGNonDegenerateAccretor >= 0.0) and (options.massTransferCriticalMassRatioHGDegenerateAccretor >= 0.0)){
                options.massTransferCriticalMassRatioHG = true;
            }
            else{
                options.massTransferCriticalMassRatioHG = false;
            }

            if(vm.count("critical-mass-ratio-giant-non-degenerate-accretor") and vm.count("critical-mass-ratio-giant-degenerate-accretor") and (options.massTransferCriticalMassRatioGiantNonDegenerateAccretor >= 0.0) and (options.massTransferCriticalMassRatioGiantDegenerateAccretor >= 0.0)){
                options.massTransferCriticalMassRatioGiant = true;
            }
            else{
                options.massTransferCriticalMassRatioGiant = false;
            }

            if(vm.count("critical-mass-ratio-helium-MS-non-degenerate-accretor") and vm.count("critical-mass-ratio-helium-MS-degenerate-accretor") and (options.massTransferCriticalMassRatioHeliumMSNonDegenerateAccretor >= 0.0) and (options.massTransferCriticalMassRatioHeliumMSDegenerateAccretor >= 0.0)){
                options.massTransferCriticalMassRatioHeliumMS = true;
            }
            else{
                options.massTransferCriticalMassRatioHeliumMS = false;
            }

            if(vm.count("critical-mass-ratio-helium-HG-non-degenerate-accretor") and vm.count("critical-mass-ratio-helium-HG-degenerate-accretor") and (options.massTransferCriticalMassRatioHeliumHGNonDegenerateAccretor >= 0.0) and (options.massTransferCriticalMassRatioHeliumHGDegenerateAccretor >= 0.0)){
                options.massTransferCriticalMassRatioHeliumHG = true;
            }
            else{
                options.massTransferCriticalMassRatioHeliumHG = false;
            }
			
            if(vm.count("critical-mass-ratio-helium-giant-non-degenerate-accretor") and vm.count("critical-mass-ratio-helium-giant-degenerate-accretor") and (options.massTransferCriticalMassRatioHeliumGiantNonDegenerateAccretor >= 0.0) and (options.massTransferCriticalMassRatioHeliumGiantDegenerateAccretor >= 0.0)){
                options.massTransferCriticalMassRatioHeliumGiant = true;
            }
            else{
                options.massTransferCriticalMassRatioHeliumGiant = false;
            }

            if(vm.count("critical-mass-ratio-white-dwarf-non-degenerate-accretor") and vm.count("critical-mass-ratio-white-dwarf-degenerate-accretor") and (options.massTransferCriticalMassRatioWhiteDwarfNonDegenerateAccretor >= 0.0) and (options.massTransferCriticalMassRatioWhiteDwarfDegenerateAccretor >= 0.0)){
                options.massTransferCriticalMassRatioWhiteDwarf = true;
            }
            else{
                options.massTransferCriticalMassRatioWhiteDwarf = false;
            }

            if(options.maxEvolutionTime < 0.0){
                std::cerr << "Error: Must evolve for a positive amount of time in Myrs (i.e. 15000). Exiting" << std::endl;
                programStatus = ERROR_IN_COMMAND_LINE;
            }
            
            if(vm.count("use-mass-loss")){
                options.useMassLoss = true;
            }
            else{
                options.useMassLoss = false;
            }
			
            // Pair instability supernovae and pulsational pair instability mass loss
            if(vm.count("pair-instability-supernovae")){
            	options.usePairInstabilitySupernovae = true;
            }
            else{
            	options.usePairInstabilitySupernovae = false;
            }

            if(vm.count("pulsational-pair-instability")){
            	options.usePulsationalPairInstability = true;
            }
            else{
            	options.usePulsationalPairInstability = false;
            }

            // Specify whether to run an individual system or do population synthesis
            if(vm.count("individual-system")){
                options.individualSystem = true;
            }
            else{
                options.individualSystem = false;
            }

            if(vm.count("neutron-star-equation-of-state")){
            	if(boost::iequals(options.neutronStarEquationOfStateString, "SSE")){
            		options.neutronStarEquationOfState = NS_EOS_SSE;
            	}
            	else if(boost::iequals(options.neutronStarEquationOfStateString, "ARP3")){
            		options.neutronStarEquationOfState = NS_EOS_ARP3;
            	}
            	else{
            		std::cerr << "Unrecognised neutron-star-equation-of-state " << options.neutronStarEquationOfStateString << ". Exiting" << std::endl;
            		programStatus = ERROR_IN_COMMAND_LINE;
            	}
            }
            
            // Specify whether to evolve a single star or a binary
            if(vm.count("single-star")){
                options.singleStar = true;
            }
            else{
                options.singleStar = false;
            }

            // Specify which mass loss prescription to use
            if(vm.count("mass-loss-prescription")){
                if(boost::iequals(options.massLossPrescriptionString, "NONE")){
                    options.massLossPrescription = MASS_LOSS_PRESCRIPTION_NONE;
                }
                else if(boost::iequals(options.massLossPrescriptionString, "Hurley")){
                    options.massLossPrescription = MASS_LOSS_PRESCRIPTION_HURLEY;
                }
                else if(boost::iequals(options.massLossPrescriptionString, "Vink")){
                    options.massLossPrescription = MASS_LOSS_PRESCRIPTION_VINK;
                }
                else{
                    std::cout << "Invalid mass loss prescription requested. Exiting" << std::endl;
                    programStatus = ERROR_IN_COMMAND_LINE;
                }
            }

            // Specify which mass transfer prescription to use
            if(options.useMassTransfer){
                if(vm.count("mass-transfer-prescription")){
                    if(boost::iequals(options.massTransferPrescriptionString, "DEMINK")){
                        options.massTransferPrescription = DEMINK_MASS_TRANSFER;
                    }
                    else if(boost::iequals(options.massTransferPrescriptionString, "BELCZYNSKI")){
                        options.massTransferPrescription = BELCZYNSKI_MASS_TRANSFER;
                    }
                    else{
                        std::cout << "Invalid mass transfer prescription requested. Exiting" << std::endl;
                        programStatus = ERROR_IN_COMMAND_LINE;
                    }
                }
            }

            // Specify which mass transfer angular momentum loss prescription to use
            if(options.useMassTransfer){
                if(vm.count("mass-transfer-angular-momentum-loss-prescription")){
                    if(boost::iequals(options.massTransferAngularMomentumLossPrescriptionString, "JEANS")){
                        options.massTransferAngularMomentumLossPrescription = JEANS_MASS_TRANSFER_ANGULAR_MOMENTUM_LOSS;
                    }
                    else if(boost::iequals(options.massTransferAngularMomentumLossPrescriptionString, "ISOTROPIC")){
                        options.massTransferAngularMomentumLossPrescription = ISOTROPIC_RE_EMISSION_MASS_TRANSFER_ANGULAR_MOMENTUM_LOSS;
                    }
                    else if(boost::iequals(options.massTransferAngularMomentumLossPrescriptionString, "CIRCUMBINARY")){
                        options.massTransferAngularMomentumLossPrescription = CIRCUMBINARY_RING_MASS_TRANSFER_ANGULAR_MOMENTUM_LOSS;
                    }
                    else if(boost::iequals(options.massTransferAngularMomentumLossPrescriptionString, "ARBITRARY")){
                        options.massTransferAngularMomentumLossPrescription = ARBITRARY_MASS_TRANSFER_ANGULAR_MOMENTUM_LOSS;
                    }
                    else{
                        std::cout << "Invalid mass transfer angular momentum loss prescription requested. Exiting" << std::endl;
                        programStatus = ERROR_IN_COMMAND_LINE;
                    }
                }
            }

            // Specify which mass transfer accretion efficiency prescription to use
            if(options.useMassTransfer){
                if(vm.count("mass-transfer-accretion-efficiency-prescription")){
                    if(boost::iequals(options.massTransferAccretionEfficiencyPrescriptionString, "THERMAL")){
                        options.massTransferAccretionEfficiencyPrescription = THERMALLY_LIMITED_MASS_TRANSFER;
                    }
                    else if(boost::iequals(options.massTransferAccretionEfficiencyPrescriptionString, "FIXED")){
                        options.massTransferAccretionEfficiencyPrescription = FIXED_FRACTION_MASS_TRANSFER;
                    }
                    else if(boost::iequals(options.massTransferAccretionEfficiencyPrescriptionString, "CENTRIFUGAL")){
                        options.massTransferAccretionEfficiencyPrescription = CENTRIFUGALLY_LIMITED_MASS_TRANSFER;
                    }
                    else{
                        std::cout << "Invalid mass transfer accretion efficiency prescription requested. Exiting" << std::endl;
                        programStatus = ERROR_IN_COMMAND_LINE;
                    }
                }
            }

            // Specify which limit to use for the accretor
            if(options.useMassTransfer){
                if(vm.count("mass-transfer-thermal-limit-accretor")){
                    if(boost::iequals(options.massTransferThermallyLimitedVariationString, "CFACTOR")){

                        options.massTransferThermallyLimitedVariation = THERMAL_C_FACTOR;

                        // If user didn't specify choice of C factor, use default based on choice of thermally limited variation
                        if(!vm.count("mass-transfer-thermal-limit-C")){
                        	options.massTransferCParameter = 10.0;
                        }

                    }
                    else if(boost::iequals(options.massTransferThermallyLimitedVariationString, "ROCHELOBE")){

                        options.massTransferThermallyLimitedVariation = THERMAL_RADIUS_TO_ROCHELOBE;

                        // If user didn't specify choice of C factor, use default based on choice of thermally limited variation
                        if(!vm.count("mass-transfer-thermal-limit-C")){
                        	options.massTransferCParameter = 1.0;
                        }

                    }
                    else{
                        std::cout << "Invalid mass transfer accretor thermal limit requested. Exiting" << std::endl;
                        programStatus = ERROR_IN_COMMAND_LINE;
                    }
                }
            }


            // Specify which mass transfer rejuvenation prescription to use
            if(options.useMassTransfer){
                if(vm.count("mass-transfer-rejuvenation-prescription")){
                    if(boost::iequals(options.massTransferRejuvenationPrescriptionString, "NONE")){
                        options.massTransferRejuvenationPrescription = MASS_TRANSFER_REJUVENATION_NONE;
                    }
                    else if(boost::iequals(options.massTransferRejuvenationPrescriptionString, "STARTRACK")){
                        options.massTransferRejuvenationPrescription = MASS_TRANSFER_REJUVENATION_STARTRACK;
                    }
                    else{
                        std::cout << "Invalid mass transfer rejuvenation prescription requested. Exiting" << std::endl;
                        programStatus = ERROR_IN_COMMAND_LINE;
                    }
                }
            }

			if(options.useMassTransfer){
                if(vm.count("circulariseBinaryDuringMassTransfer")){
					options.circulariseBinaryDuringMassTransfer = true;
					
					if(vm.count("angularMomentumConservationDuringCircularisation")){
						options.angularMomentumConservationDuringCircularisation = true;
					}
					else{
						options.angularMomentumConservationDuringCircularisation = false;
					}
                }
				else{
					options.circulariseBinaryDuringMassTransfer = false;
				}
            }

			if(options.useMassTransfer){
                if(vm.count("forceCaseBBBCStabilityFlag")){
					options.forceCaseBBBCStabilityFlag = true;

	                if(vm.count("alwaysStableCaseBBBCFlag")){
						options.alwaysStableCaseBBBCFlag = true;
					}
					else{
						options.alwaysStableCaseBBBCFlag = false;	
					}
                }
				else{
					options.forceCaseBBBCStabilityFlag = false;
				}
            }

            // Check common envelope parameters are sensible
            if(vm.count("common-envelope-alpha")){
            	if(options.commonEnvelopeAlpha < 0){
            		std::cerr << "Error: common-envelope-alpha < 0 requested. Exiting" << std::endl;
            		programStatus = ERROR_IN_COMMAND_LINE;
            	}
            }

            if(vm.count("common-envelope-lambda-multiplier")){
            	if(options.commonEnvelopeLambdaMultiplier < 0){
            		std::cerr << "Error: common-envelope-lambda-multiplier < 0 requested. Exiting" << std::endl;
            		programStatus = ERROR_IN_COMMAND_LINE;
            	}
            }

            if(vm.count("common-envelope-mass-accretion-prescription")){
            	if(boost::iequals(options.commonEnvelopeMassAccretionPrescriptionString, "ZERO")){
                    options.commonEnvelopeMassAccretionPrescription = COMMON_ENVELOPE_ACCRETION_ZERO;
                }
                else if(boost::iequals(options.commonEnvelopeMassAccretionPrescriptionString, "CONSTANT")){
                    options.commonEnvelopeMassAccretionPrescription = COMMON_ENVELOPE_ACCRETION_CONSTANT;
                }
                else if(boost::iequals(options.commonEnvelopeMassAccretionPrescriptionString, "UNIFORM")){
                    options.commonEnvelopeMassAccretionPrescription = COMMON_ENVELOPE_ACCRETION_UNIFORM;
                }
                else if(boost::iequals(options.commonEnvelopeMassAccretionPrescriptionString, "MACLEOD+2014")){
                    options.commonEnvelopeMassAccretionPrescription = COMMON_ENVELOPE_ACCRETION_MACLEOD;
                }
                else{
                    std::cout << "Unrecognised CE option." << options.commonEnvelopeMassAccretionPrescriptionString <<  "Exiting." << std::endl;
                    programStatus = ERROR_IN_COMMAND_LINE;
                }
            }

            if(vm.count("common-envelope-mass-accretion-min")){
            	if(options.commonEnvelopeMassAccretionMin < 0.0){
            		std::cout << "Minimum accreted mass should be > 0, got " << options.commonEnvelopeMassAccretionMin << std::endl;
            		programStatus = ERROR_IN_COMMAND_LINE;
            	}
            }

            if(vm.count("common-envelope-mass-accretion-max")){
            	if(options.commonEnvelopeMassAccretionMax < 0.0){
            		std::cout << "Maximum accreted mass should be > 0, got " << options.commonEnvelopeMassAccretionMax << std::endl;
            		programStatus = ERROR_IN_COMMAND_LINE;
            	}
            }

            if(vm.count("common-envelope-mass-accretion-constant")){
            	if(options.commonEnvelopeMassAccretionConstant < 0.0){
            		std::cout << "commonEnvelopeMassAccretionConstant should be >= 0 " << std::endl;
            		programStatus = ERROR_IN_COMMAND_LINE;
            	}
            }

            if(vm.count("common-envelope-alpha-thermal")){
            	if(options.commonEnvelopeAlphaThermal < 0 or options.commonEnvelopeAlphaThermal > 1){
            		std::cerr << "Error: common-envelope-alpha-thermal should be between 0 and 1. Exiting" << std::endl;
            		programStatus = ERROR_IN_COMMAND_LINE;
            	}
            }
			
            // Specify which HG CE assumption to make
            if(vm.count("common-envelope-hertzsprung-gap-assumption")){
                if(boost::iequals(options.commonEnvelopeHertzsprungGapDonorString, "OPTIMISTIC_HG_CE")){
                    options.commonEnvelopeHertzsprungGapDonor = OPTIMISTIC_HG_CE;
                }
                else if(boost::iequals(options.commonEnvelopeHertzsprungGapDonorString, "PESSIMISTIC_HG_CE")){
                    options.commonEnvelopeHertzsprungGapDonor = PESSIMISTIC_HG_CE;
                }
                else if(boost::iequals(options.commonEnvelopeHertzsprungGapDonorString, "STABLE_HG_CE")){
                    options.commonEnvelopeHertzsprungGapDonor = STABLE_HG_CE;
                }
                else{
                    std::cout << "Unrecognised CE option. Exiting." << std::endl;
                    programStatus = ERROR_IN_COMMAND_LINE;
                }
            }

			// Specify which CE lambda to use
            if(vm.count("common-envelope-lambda-prescription")){
                if(boost::iequals(options.commonEnvelopeLambdaPrescriptionString, "LAMBDA_FIXED")){
                    options.commonEnvelopeLambdaPrescription = LAMBDA_FIXED;
                }
                else if(boost::iequals(options.commonEnvelopeLambdaPrescriptionString, "LAMBDA_LOVERIDGE")){
                    options.commonEnvelopeLambdaPrescription = LAMBDA_LOVERIDGE;
                }
				else if(boost::iequals(options.commonEnvelopeLambdaPrescriptionString, "LAMBDA_NANJING")){
                    options.commonEnvelopeLambdaPrescription = LAMBDA_NANJING;
                }
				else if(boost::iequals(options.commonEnvelopeLambdaPrescriptionString, "LAMBDA_KRUCKOW")){
                    options.commonEnvelopeLambdaPrescription = LAMBDA_KRUCKOW;

					// If user didn't specify choice of the slope, use default
                    if(!vm.count("common-envelope-slope-Kruckow")){
						options.commonEnvelopeSlopeKruckow = -4.0/5.0;
					}
                }
                else if(boost::iequals(options.commonEnvelopeLambdaPrescriptionString, "LAMBDA_DEWI")){
                    options.commonEnvelopeLambdaPrescription = LAMBDA_DEWI;
                }
                else{
                    std::cout << "Unrecognised CE lambda prescription. Exiting." << std::endl;
                    programStatus = ERROR_IN_COMMAND_LINE;
                }
            }

            // Specify which CE zeta prescription to use
            if(vm.count("common-envelope-zeta-prescription")){
                if(boost::iequals(options.commonEnvelopeZetaPrescriptionString, "STARTRACK")){
                    options.commonEnvelopeZetaPrescription = ZETA_STARTRACK;
                }
                else if(boost::iequals(options.commonEnvelopeZetaPrescriptionString, "SOBERMAN")){
                    options.commonEnvelopeZetaPrescription = ZETA_SOBERMAN;
                }
				else if(boost::iequals(options.commonEnvelopeZetaPrescriptionString, "HURLEY")){
                    options.commonEnvelopeZetaPrescription = ZETA_HURLEY;
                }
				else if(boost::iequals(options.commonEnvelopeZetaPrescriptionString, "ARBITRARY")){
					options.commonEnvelopeZetaPrescription = ZETA_ARBITRARY;

					// If used didn't specify zetas, asign large value, which will favour stable MT
                    if(!vm.count("zeta-adiabatic-arbitrary")){
						options.zetaAdiabaticArbitrary = 10000.0;
					}

                    if(!vm.count("zeta-thermal-arbitrary")){
						options.zetaThermalArbitrary = 10000.0;
					}					
				}
                else{
                    std::cout << "Unrecognised CE zeta prescription. Exiting." << std::endl;
                    programStatus = ERROR_IN_COMMAND_LINE;
                }
            }

			if(vm.count("lambda-calculation-every-timeStep")){
                options.lambdaCalculationEveryTimeStep = true;
            }
            else{
                options.lambdaCalculationEveryTimeStep = false;
            }

            if(vm.count("revised-energy-formalism-Nandez-Ivanova")){
                options.revisedEnergyFormalismNandezIvanova = true;
            }
            else{
                options.revisedEnergyFormalismNandezIvanova = false;
            }
			
			if(vm.count("zeta-Calculation-Every-Time-Step")){
                options.zetaCalculationEveryTimeStep = true;
            }
            else{
                options.zetaCalculationEveryTimeStep = false;
            }

			if(vm.count("BeBinaries")){
				options.BeBinaries = true;
			}
            if(vm.count("CHEvolution")){
				options.CHEvolution = true;
			}


			// Floor 24-04-2018 - add options for Exploratory phase of Adaptive Importance Sampling
            if(vm.count("AIS-exploratory-phase")){
				options.AISexploratoryphase = true;  
			}
			else{
				options.AISexploratoryphase = false; 
			}

             // Specify DCO type to select for
            if(vm.count("AIS-DCOtype")){
            	// Check Exploratory phase settings for Adaptive Importance Sampling
                if(boost::iequals(options.AISDCOtypeString, "ALL")){
                    options.AISDCOtype = ALL;
                }
                else if(boost::iequals(options.AISDCOtypeString, "BBH")){
                    options.AISDCOtype = BBH;
                }
                else if(boost::iequals(options.AISDCOtypeString, "BNS")){
                    options.AISDCOtype = BNS; 
                }               
                else if(boost::iequals(options.AISDCOtypeString, "BHNS")){
                    options.AISDCOtype = BHNS;
                }
                else{
                    std::cout << "Unrecognised AISDCOtype option. Exiting." << std::endl;
                    programStatus = ERROR_IN_COMMAND_LINE;
                }
            }

            // Check if excluding binaries that merge outside Hubble time (exploratory phase AIS)
            if(vm.count("AIS-Hubble")){
				options.AISHubble = true;  
			}
			else{
				options.AISHubble = false; 
			}

            // Check if excluding binaries that RLOFSecondaryZAMS  (exploratory phase AIS)
            if(vm.count("AIS-RLOF")){
				options.AISRLOF = true;  
			}
			else{
				options.AISRLOF = false; 
			}

			// Check if excluding binaries that are Optimistic  (exploratory phase AIS)
            if(vm.count("AIS-Pessimistic")){
				options.AISPessimistic = true;  
			}
			else{
				options.AISPessimistic = false; 
			}

            if(vm.count("AIS-refinement-phase")){
				options.AISrefinementPhase = true;  
			}
			else{
				options.AISrefinementPhase = false; 
			}

            if(vm.count("RLOFPrinting")){
				options.RLOFPrinting = true;
			}

            // Check if user set the metallicity to use
            if(vm.count("metallicity")){
            	if(options.metallicity < 0 or options.metallicity > 1){
            		std::cerr << "Error: metallicity should be absolute metallicity and should be between 0 and 1. Exiting" << std::endl;
            		programStatus = ERROR_IN_COMMAND_LINE;
            	}
                options.fixedMetallicity = true;
            }

            if(vm.count("pulsar-magnetic-field-decay-timescale")){
            	if(options.pulsarMagneticFieldDecayTimescale < 0.0){
            		std::cerr << "Error: pulsar-magnetic-field-decay-timescale should be > 0" << std::endl;
            		programStatus = ERROR_IN_COMMAND_LINE;
            	}
            }

            if(vm.count("pulsar-magnetic-field-decay-massscale")){
            	if(options.pulsarMagneticFieldDecayMassscale < 0.0){
            		std::cerr << "Error: pulsar-magnetic-field-decay-massscale should be > 0" << std::endl;
            		programStatus = ERROR_IN_COMMAND_LINE;
            	}
            }

            // Whether to evolve pulsars
            if(vm.count("evolve-pulsars")){
            	options.evolvePulsars = true;       
            }

            // Check if user wants a detailed output of each simualted system
            if(vm.count("detailedOutput")){
                options.detailedOutput = true;
            }
            else{
                options.detailedOutput = false;
            }
            // Check if user wants the output of each rejected system
            if(vm.count("RejectedSamplesPrinting")){
                options.RejectedSamplesPrinting = true;
            }
            else{
                options.RejectedSamplesPrinting = false;
            }
            // Check if user wants to print certain values while running a population
            if(vm.count("populationDataPrinting")){
                options.populationDataPrinting = true;
            }
            else{
                options.populationDataPrinting = false;
            }

			//JIM BARRETT -- 06/07/2016 -- adding options to sample over some hyperparameters
			if(vm.count("sample-kick-velocity-sigma")){
				options.sampleKickVelocitySigma = true;
			}
			else{
				options.sampleKickVelocitySigma = false;
			}

			if(vm.count("sample-kick-direction-power")){
				options.sampleKickDirectionPower = true;
			}
			else{
				options.sampleKickDirectionPower = false;
			}

			if(vm.count("sample-common-envelope-alpha")){
				options.sampleCommonEnvelopeAlpha = true;
			}
			else{
				options.sampleCommonEnvelopeAlpha = false;
			}

			if(vm.count("sample-wolf-rayet-multiplier")){
				options.sampleWolfRayetMultiplier = true;
			}
			else{
				options.sampleWolfRayetMultiplier = false;
			}

			if(vm.count("sample-luminous-blue-variable-multiplier")){
				options.sampleLuminousBlueVariableMultiplier = true;
			}
			else{
				options.sampleLuminousBlueVariableMultiplier = false;
			}

			// Simon Stevenson - 15/03/2018 - adding MCMC functionality
			if(vm.count("mcmc")){
				options.useMCMC = true;
			}
			else{
				options.useMCMC = false;
			}

			if(vm.count("importance-sampling")){
				options.useImportanceSampling = true;
			}
			else{
				options.useImportanceSampling = false;
			}			

			if(vm.count("pulsational-pair-instability-prescription")){
				if(boost::iequals(options.pulsationalPairInstabilityPrescriptionString, "COMPAS")){
					options.pulsationalPairInstabilityPrescription = COMPASPPI;
				}
				else if(boost::iequals(options.pulsationalPairInstabilityPrescriptionString, "STARTRACK")){
					options.pulsationalPairInstabilityPrescription = STARTRACKPPI;
				}
				else if(boost::iequals(options.pulsationalPairInstabilityPrescriptionString, "MARCHANT")){
					options.pulsationalPairInstabilityPrescription = MARCHANTPPI;
				}
				else{
					std::cerr << "Error: Unrecognised pulsational-pair-instability-prescription . Exiting" << std::endl;
					programStatus = ERROR_IN_COMMAND_LINE;
				}
			}
        }
        catch(po::error& e){
            std::cerr << "Error: " << e.what() << std::endl;
            std::cerr << desc << std::endl;
            programStatus = ERROR_IN_COMMAND_LINE;
        }

    } catch (std::exception& e) {

        std::cerr << "Unhandled exception reached top of main : " << e.what() << ", application will now exit." << std::endl;
        programStatus = ERROR_UNHANDLED_EXCEPTION;
    }


}
