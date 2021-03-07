//
//  programOptions.cpp

#include "programOptions.h"

// Default constructor
programOptions::programOptions(){
    
    // This sets all of the program options to their default values -- can be modified via the command line.
    
    // Individual system run? default = true
//    individualSystem = true;
    individualSystem = false;                                   // Flag to evolve a specific individual system which you can specify initial parameters of
    
    singleStar = false;                                         // Flag to evolve a single star

    // Debugging? default = false
    debugging = false;
    
    //suppress some of the printing
    quiet = false;                                          // Flag to activate debugging statements
 
    RejectedSamplesPrinting = false;                           // Flag to turn on printing of the samples that are rejected by COMPAS    
    onlyDoubleCompactObjects = false;                           // Flag to turn on some shortcuts to only evolve systems which may form double compact objects

	evolveUnboundSystems = false;
    
	// Individual system variables
    primaryMass = 20.0;                                         // Initial primary mass in solar masses
    secondaryMass = 10.0;                                       // Initial secondary mass in solar masses
    
    initialPrimaryMetallicity = 0.02;                           // Initial metallicity of the primary
    initialSecondaryMetallicity = 0.02;                         // Initial metallicity of the secondary
    
    // Variables required to restart a binary/star halfway through
    primaryStellarType = -1;                                    // Initial primary stellar type (not yet implemented)
    secondaryStellarType = -1;                                  // Initial secondary stellar type (not yet implemented)
    
    primaryEffectiveInitialMass = primaryMass;                  // Effective initial mass for the primary in solar masses
    secondaryEffectiveInitialMass = secondaryMass;              // Effective initial mass for the secondary in solar masses
    
    primaryCoreMass = 0.0;                                      // Initial primary core mass in solar masses
    secondaryCoreMass = 0.0;                                    // Initial secondary core mass in solar masses
    
    primaryAge = 0.0;                                           // Effective age for the primary star in Myrs
    secondaryAge = 0.0;                                         // Effective age for the secondary star in Myrs
    
    binarySeparation = -1.0;                                    // Initial separation in AU
    binaryOrbitalPeriod = -1.0;                                 // Initial orbital period in day
    binaryEccentricity = 0.0;                                   // Initial eccentricity
    
    primaryRotationalVelocity = 0.0;                            // Initial rotational velocity of the primary
    secondaryRotationalVelocity = 0.0;                          // Initial rotational velocity of the secondary
    
    // Public population synthesis variables
    nBinaries = 10; 
    
    fixedRandomSeed = false;                                    // Whether to use a fixed random seed given by options.randomSeed (set to true if --random-seed is passed on command line) (default = false)
    randomSeed = 0;                                             // Random seed to use (default = 0)
    
    // Specify how long to evolve binaries for
    maxEvolutionTime = 13700.0;                                 // Maximum evolution time in Myrs
    maxNumberOfTimestepIterations = 99999;                      // Maximum number of timesteps to evolve binary for before giving up
    
    // Initial mass options
    initialMassFunction = "Kroupa";
    initialMassFunctionMin   = 8.0;
    initialMassFunctionMax   = 100.0;
    initialMassFunctionPower = -2.3;
    
    // Initial mass ratios
    massRatioDistribution = "Flat";                             // Most likely want Flat or SANA2012
    massRatioDistributionMin = 0.0;
    massRatioDistributionMax = 1.0;

    minimumMassSecondary = 0.0;                                 // Minimum mass of secondary to draw (in Msol)
    
    // Initial orbit options
    semiMajorAxisDistribution   = "FlatInLog";                  // Most likely want FlatInLog or SANA2012
    semiMajorAxisDistributionMin     = 0.1;
    semiMajorAxisDistributionMax     = 1000;
    semiMajorAxisDistributionPower   = -1;
    
    // Initial orbital period
    periodDistributionMin = 1.1;                                // Minimum initial period in days
    periodDistributionMax = 1000.0;                             // Maximum initial period in days
    
    // Eccentricity
    eccentricityDistribution = "Zero";                          // Most likely want ZERO, THERMAL or SANA2012
    eccentricityDistributionMin = 0.0;                          // Minimum initial eccentricity to sample
    eccentricityDistributionMax = 1.0;                          // Maximum initial eccentricity to sample
    
    // Kick options
    kickVelocityDistribution				= "Maxwellian";		// Which kick velocity distribution to use
    //kickVelocityDistributionSigma   		= 250;              // Characteristic kick velocity in km s^-1. SIMON : REPLACED
    kickVelocityDistributionSigmaCCSN_NS    = 250;              // Kick velocity sigma in km s^-1 for neutron stars (default = "250" )
    kickVelocityDistributionSigmaCCSN_BH    = 250;              // Kick velocity sigma in km s^-1 for black holes (default = "250" )
    kickVelocityDistributionMaximum         = -1.0;             // Maximum kick velocity to draw in km s^-1. Ignored if < 0
    kickVelocityDistributionSigmaForECSN   	= 30.0;              // Characteristic kick velocity for an ECSN in km s^-1
    kickVelocityDistributionSigmaForUSSN   	= 30.0;             // Characteristic kick velocity for an USSN in km s^-1
	kickScalingFactor						= 1.0;				// Arbitrary factor for scaling kicks
    
    // Black hole kicks
    blackHoleKicks = "FALLBACK";
    blackHoleKicksOption = FALLBACK;
    
    // Supernova remnant mass prescription options
    // 
    // postitNote
    // hurley2000
    // belczynski2002
    // fryer2012
    //
    
    remnantMassPrescription = "fryer2012";
    fryerSupernovaEngineString = "Delayed";
    fryerSupernovaEngine = SN_DELAYED;

    neutrinoMassLossAssumptionBHString = "FIXED_FRACTION";
    neutrinoMassLossAssumptionBH = NEUTRINO_MASS_LOSS_BH_FIXED_FRACTION;  // Assumption to make about neutrino mass loss for BH formation
    neutrinoMassLossValueBH = 0.1;                           // Value (corresponding to assumption) for neutrino mass loss for BH formation

    // Fixed uk options
    useFixedUK = false;
    fixedUK = -1.0;

    // Pair instability and pulsational pair instability mass loss
    usePairInstabilitySupernovae = true;                           // Whether to use pair instability supernovae (PISN)
    pairInstabilityUpperLimit = 135.0;                             // Maximum core mass leading to PISN (default = 135, Value in Belczynski+ 2016 is 135 Msol)
    pairInstabilityLowerLimit = 60.0;                              // Minimum core mass leading to PISN (default = 65,  Value in Belczynski+ 2016 is 65 Msol)

    usePulsationalPairInstability = true;                          // Whether to use pulsational pair instability (PPI)
    pulsationalPairInstabilityLowerLimit = 35.0;                   // Minimum core mass leading to PPI, default = 40, Value in Belczynski+ 2016 is 45 Msol
    pulsationalPairInstabilityUpperLimit = 60.0;                   // Maximum core mass leading to PPI, default = 65, Value in Belczynski+ 2016 is 65 Msol
    
    pulsationalPairInstabilityPrescriptionString = "COMPAS";        // String for which PPI prescription to use
    pulsationalPairInstabilityPrescription = COMPASPPI;             // Prescription for PPI to use

	maximumNeutronStarMass = 3.0;									// Maximum mass of a neutron star allowed, set to default in StarTrack
    nBatchesUsed = -1;                                              // nr of batches used, only needed for STROOPWAFEL (AIS) (default = -1, not needed)		
    // Kick direction option
    kickDirection = "Isotropic";                                    // Which assumption for SN kicks: Possibilities: Isotropic, poles, plane, power
    kickDirectionPower = 0.0;                                       // Power law power for the "power" SN kick direction choice
    
    // Get default output path
    outputPathString = ".";                                         // String to hold the output directory
    defaultOutputPath = boost::filesystem::current_path();          // Default output location
    outputPath = defaultOutputPath;                                 // Desired output location (default = CWD)

    // Spin options
    spinDistributionMin = 0.60;
    spinDistributionMax = 0.98;
    spinDistribution = "Fixed";
    spinAssumption = "bothAligned";
    
    // PN evolution options
    PNevolution = false;
    
    // Tides options
    tidesPrescriptionString = "None";                                   // String containing which tides prescription to use (default = "None")
    tidesPrescription = TIDES_PRESCRIPTION_NONE;                        // Tides prescription that will be used by the code
    // Tides prescription options
    //
    // NONE TIDES_PRESCRIPTION_NONE
    // LOCKED TIDES_PRESCRIPTION_LOCKED
    // HUT TIDES_PRESCRIPTION_HUT
    //
    
    // Magnetic Braking options
    //    useMagneticBraking = false;                                         // Wether to use Magnetic Braking
    
    // Mass loss options
    useMassLoss = true;                                                // Whether to use mass loss
    
    massLossPrescriptionString = "VINK";
    massLossPrescription = MASS_LOSS_PRESCRIPTION_VINK;
    // Mass loss prescription options
    // NONE     MASS_LOSS_PRESCRIPTION_NONE
    // HURLEY   MASS_LOSS_PRESCRIPTION_HURLEY
    // VINK     MASS_LOSS_PRESCRIPTION_VINK
    
    // Wind mass loss multiplicitive constants
    luminousBlueVariableFactor = 1.5;                         // Luminous blue variable mass loss enhancement factor
    wolfRayetFactor = 1.0;                                    // WR winds factor
    
    // Mass transfer options
    useMassTransfer = true;                                             // Whether to use mass transfer (default = true)
	circulariseBinaryDuringMassTransfer	=	false;						// Whether to circularise binary when it starts (default = false)
	forceCaseBBBCStabilityFlag = true;									// Whether if all case BB/BC systems are forced to be stable or unstable
	alwaysStableCaseBBBCFlag = true;									// Whether if case BB/BC is always stable
    
	massTransferPrescriptionString = "DEMINK";
    massTransferPrescription = DEMINK_MASS_TRANSFER;
    // Mass transfer prescription options
    // DEMINK       DEMINK_MASS_TRANSFER
    // BELCZYNSKI   BELCZYNSKI_MASS_TRANSFER

    // Options adaptive Roche Lobe Overflow prescription
    massTransferAdaptiveAlphaParameter = 0.5;
    maxPercentageAdaptiveMassTransfer = 0.01;
    // Options for mass transfer accretion efficiency
    massTransferFractionAccreted = 1.0;
    massTransferCParameter = 10.0;
    massTransferAccretionEfficiencyPrescriptionString = "THERMAL";
    massTransferAccretionEfficiencyPrescription = THERMALLY_LIMITED_MASS_TRANSFER;
    // Mass transfer accretion efficiency prescription options
    // THERMAL          THERMALLY_LIMITED_MASS_TRANSFER
    // FIXED            FIXED_FRACTION_MASS_TRANSFER
    // CENTRIFUGAL      CENTRIFUGALLY_LIMITED_MASS_TRANSFER
	
	massTransferThermallyLimitedVariationString = "CFACTOR";
	massTransferThermallyLimitedVariation = THERMAL_C_FACTOR;
//	massTransferThermallyLimitedVariation = THERMAL_RADIUS_TO_ROCHELOBE;
	// CFACTOR			THERMAL_C_FACTOR
	// ROCHELOBE		THERMAL_RADIUS_TO_ROCHELOBE

    eddingtonAccretionFactor = 1;      // Multiplication factor for eddington accretion for NS & BH
                                         //i.e. >1 is super-eddington 
                                         //0. is no accretion   
	
    massTransferJloss= 1.0;
    massTransferAngularMomentumLossPrescriptionString = "ISOTROPIC";
    massTransferAngularMomentumLossPrescription = ISOTROPIC_RE_EMISSION_MASS_TRANSFER_ANGULAR_MOMENTUM_LOSS;
    // Mass transfer angular momentum loss prescription options
    // JEANS           JEANS_MASS_TRANSFER_ANGULAR_MOMENTUM_LOSS
    // ISOTROPIC       ISOTROPIC_RE_EMISSION_MASS_TRANSFER_ANGULAR_MOMENTUM_LOSS
    // CIRCUMBINARY    CIRCUMBINARY_RING_MASS_TRANSFER_ANGULAR_MOMENTUM_LOSS
    // ARBITRARY       ARBITRARY_MASS_TRANSFER_ANGULAR_MOMENTUM_LOSS
    
    massTransferRejuvenationPrescriptionString = "NONE";
    massTransferRejuvenationPrescription = MASS_TRANSFER_REJUVENATION_NONE;
    // Mass transfer rejuvenation prescriptions
    // "Belczynski weird stuff no 77"
    // Which prescription to use for rejuvenating stars when they accrete material during mass transfer
    // NONE             MASS_TRANSFER_REJUVENATION_NONE                 // The BSE prescription
    // STARTRACK        MASS_TRANSFER_REJUVENATION_STARTRACK            // The StarTrack prescription described in section 5.6 of Belczynski 2008 http://arxiv.org/pdf/astro-ph/0511811v3.pdf
    //
    
    // Mass transfer critical mass ratios
    massTransferCriticalMassRatioMSLowMass = false;                    			// Whether to use critical mass ratios
    massTransferCriticalMassRatioMSLowMassNonDegenerateAccretor = 1.44;         // Critical mass ratio for MT from a MS low mass star (default = 1.44, Claeys+ 2014)
    massTransferCriticalMassRatioMSLowMassDegenerateAccretor = 1.0;             // Critical mass ratio for MT from a MS low mass star on to a degenerate accretor (default = 1.0, Claeys+ 2014)

    massTransferCriticalMassRatioMSHighMass = false;							// Whether to use critical mass ratios
    massTransferCriticalMassRatioMSHighMassNonDegenerateAccretor = 0.625;       // Critical mass ratio for MT from a MS high mass star (default = 0.625, Claeys+ 2014)
    massTransferCriticalMassRatioMSHighMassDegenerateAccretor = 0.0;            // Critical mass ratio for MT from a MS high mass star on to a degenerate accretor (default = 0)
    
    massTransferCriticalMassRatioHG = false;                                    // Whether to use critical mass ratios
    massTransferCriticalMassRatioHGNonDegenerateAccretor = 0.40;                // Critical mass ratio for MT from a HG star (default = 0.25, Claeys+ 2014)
    massTransferCriticalMassRatioHGDegenerateAccretor = 0.21;                   // Critical mass ratio for MT from a HG star on to a degenerate accretor (default = 0.21, Claeys+ 2014)
    
    massTransferCriticalMassRatioGiant = false;                                 // Whether to use critical mass ratios
    massTransferCriticalMassRatioGiantNonDegenerateAccretor = 0.0;              // Critical mass ratio for MT from a giant (default = 0.0)
    massTransferCriticalMassRatioGiantDegenerateAccretor = 0.87;                // Critical mass ratio for MT from a giant on to a degenerate accretor (default = 0.81, Claeys+ 2014)

    massTransferCriticalMassRatioHeliumMS = false;								// Whether to use critical mass ratios
    massTransferCriticalMassRatioHeliumMSNonDegenerateAccretor = 0.625;       	// Critical mass ratio for MT from a Helium MS star (default = 0.625)
    massTransferCriticalMassRatioHeliumMSDegenerateAccretor = 0.0;            	// Critical mass ratio for MT from a Helium MS star on to a degenerate accretor (default = 0)
	
    massTransferCriticalMassRatioHeliumHG = false;                              // Whether to use critical mass ratios
    massTransferCriticalMassRatioHeliumHGNonDegenerateAccretor = 0.25;			// Critical mass ratio for MT from a Helium HG star (default = 0.25, de Claeys+ 2014)
    massTransferCriticalMassRatioHeliumHGDegenerateAccretor = 0.21;				// Critical mass ratio for MT from a Helium HG star on to a degenerate accretor (default = 0.21, Claeys+ 2014)
    
    massTransferCriticalMassRatioHeliumGiant = false;                           // Whether to use critical mass ratios
    massTransferCriticalMassRatioHeliumGiantNonDegenerateAccretor = 1.28;       // Critical mass ratio for MT from a Helium giant (default = 0.25, Claeys+ 2014)
    massTransferCriticalMassRatioHeliumGiantDegenerateAccretor = 0.87;          // Critical mass ratio for MT from a Helium giant on to a degenerate accretor
	
    massTransferCriticalMassRatioWhiteDwarf = false;                           // Whether to use critical mass ratios
	massTransferCriticalMassRatioWhiteDwarfNonDegenerateAccretor = 0.0;       // Critical mass ratio for MT from a White Dwarf (default = 0.0)
    massTransferCriticalMassRatioWhiteDwarfDegenerateAccretor = 1.6;          // Critical mass ratio for MT from a White Dwarf on to a degenerate accretor (default = 1.6, Claeys+ 2014)
    
	
    // Common Envelope parameters
    commonEnvelopePrescriptionFlag = WEBBINK_COMMON_ENVELOPE;          // Which common envelope prescription to use
    commonEnvelopeAlpha = 1.0;                                         // Common envelope efficiency alpha parameter (default = 1.0)
    commonEnvelopeLambda = 0.1;                                        // Common envelope Lambda parameter (default = 0.1)
    commonEnvelopeHertzsprungGapDonorString = "OPTIMISTIC_HG_CE";      // String containing which prescription to use for Hertzsprung gap donors in a CE (default = "OPTIMISTIC_HG_CE")
    commonEnvelopeHertzsprungGapDonor = OPTIMISTIC_HG_CE;              // Which prescription to use for Hertzsprung gap donors in a CE (default = OPTIMISTIC_HG_CE)
	commonEnvelopeAlphaThermal = 1.0;                                  // lambda = (alpha_th * lambda_b) + (1-alpha_th) * lambda_g
    commonEnvelopeLambdaMultiplier = 1.0;                              // Multiply common envelope lambda by some constant
    allowMainSequenceStarToSurviveCommonEnvelope = false;              // Whether or not to allow a main sequence star to survive a common envelope event

    // Accretion during common envelope
    commonEnvelopeMassAccretionPrescriptionString = "ZERO";
    commonEnvelopeMassAccretionPrescription = COMMON_ENVELOPE_ACCRETION_ZERO;
    commonEnvelopeMassAccretionMin = 0.04; // Minimum amount of mass accreted during CE in solar masses
    commonEnvelopeMassAccretionMax = 0.1; // Maximum amount of mass accreted during CE in solar masses
    commonEnvelopeMassAccretionConstant = 0.0; // Constant value

	// Common envelope lambda prescription
	commonEnvelopeLambdaPrescriptionString = "LAMBDA_NANJING";			// String containing which prescription to use for CE lambda (default = "LAMBDA_NANJING")
	commonEnvelopeLambdaPrescription = LAMBDA_NANJING;					// Which prescription to use for CE lambda (default = LAMBDA_NANJING)
	lambdaCalculationEveryTimeStep = false;
	
	// Common envelope Nandez and Ivanova energy formalism
	revisedEnergyFormalismNandezIvanova	= false;						// Use the revised energy formalism from Nandez & Ivanova 2016 (default = false)
	maximumMassDonorNandezIvanova = 2.0;								// Maximum mass allowed to use the revised energy formalism in Msol (default = 2.0)
	commonEnvelopeRecombinationEnergyDensity = 1.5E13;					// Factor using to calculate the binding energy depending on the mass of the envelope. (default = 1.5x10^13 ergs/g)
	
	// Common envelope power factor for Kruckow fit
	commonEnvelopeSlopeKruckow = -4.0/5.0;								// Common envelope power factor for Kruckow fit normalized according to Kruckow+2016, Fig. 1
		
	// Which prescription to use for calculating zetas
	commonEnvelopeZetaPrescriptionString = "STARTRACK";						// String containing which prescription to use for calculating CE zetas (default = STARTRACK)
	commonEnvelopeZetaPrescription = ZETA_STARTRACK;					// Which prescription to use for calculating CE zetas (default = ZETA_ADIABATIC)
	zetaCalculationEveryTimeStep = false;
	zetaAdiabaticArbitrary = NEVER_SET;
	zetaThermalArbitrary = NEVER_SET;	
    zetaMainSequence 	= 2.0;
	zetaHertzsprungGap	= 6.5;
	
	// Which prescription to use for calculating zetas
	commonEnvelopeZetaPrescriptionString = "STARTRACK";						// String containing which prescription to use for calculating CE zetas (default = STARTRACK)
	commonEnvelopeZetaPrescription = ZETA_STARTRACK;					// Which prescription to use for calculating CE zetas (default = ZETA_ADIABATIC)
	
	BeBinaries = false;	
    CHEvolution = false;                                                        // Flag for printing the specific output for Chemically Homogenous evolution
    RLOFPrinting = false;
    
    // Afaptive Importance Sampling options
    AISexploratoryphase = false;                                            // Flag for whether to run the AIS exploratory phase
    AISDCOtypeString = "ALL";                                               // String containing which type of DCOs to focus on (default = "ALL")
    AISDCOtype = ALL;                                                       // Which prescription to use for DCO type (default = ALL)                                            
    AISHubble = false;                                                      // Flag for excluding DCOs that do not merge in Hubble
    AISRLOF = false;                                                        // Flag for excluding DCOs that RLOFSecondaryZAMS
    AISPessimistic = false;                                                 // Flag for excluding DCOs that are Optmistic
    AISrefinementPhase = false;                                             // Flag for whether to run the AIS refinement phase (step 2)
    KappaGaussians = 2;                                                     // scaling factor for the width of the Gaussian distributions in AIS main sampling phase


    // Metallicity options
    metallicity = Zsol;
    fixedMetallicity = true;
    
    // Rotational velocity distribution options
    rotationalVelocityString = "ZERO";                      //
    
    // Pulsar birth magnetic field distribution
    pulsarBirthMagneticFieldDistributionString = "ZERO";  // Which birth magnetic field distribution to use for pulsars
    pulsarBirthMagneticFieldDistributionMin = 11.0;       // Minimum pulsar birth magnetic field distribution
    pulsarBirthMagneticFieldDistributionMax = 13.0;       // Maximum pulsar birth magnetic field distribution
    
    // Pulsar birth spin period distribution string
    pulsarBirthSpinPeriodDistributionString = "ZERO";     // Which birth spin period distribution to use for pulsars
    pulsarBirthSpinPeriodDistributionMin = 0.0;           // Minimum birth spin period (ms)
    pulsarBirthSpinPeriodDistributionMax = 100.0;         // Maximum birth spin period (ms)

    pulsarMagneticFieldDecayTimescale = 1000.0;           // Timescale on which magnetic field decays (Myrs)
    pulsarMagneticFieldDecayMassscale = 0.025;            // Mass scale on which magnetic field decays during accretion (solar masses)
    pulsarLog10MinimumMagneticField = 8.0;                // log10 of the minimum pulsar magnetic field in Gauss

    evolvePulsars = false;                                // Whether to evolve pulsars

    // Detailed output
    detailedOutput = false;
    
    // Print certain data for small populations, but not for larger one
    populationDataPrinting = false;

	//JIM BARRETT -- 06/07/2016 -- adding options to sample over some hyperparameters
	sampleKickVelocitySigma = false;
	sampleKickVelocitySigmaMin = 0.0;
	sampleKickVelocitySigmaMax = 400.0;
	
	sampleKickDirectionPower = false;
	sampleKickDirectionPowerMin = -10.0;
	sampleKickDirectionPowerMax = 10.0;
	
	sampleCommonEnvelopeAlpha = false;
	sampleCommonEnvelopeAlphaMin = 0.0;
	sampleCommonEnvelopeAlphaMax = 5.0;

	sampleWolfRayetMultiplier = false;
	sampleWolfRayetMultiplierMin = 0.0;
	sampleWolfRayetMultiplierMax = 5.0;

	sampleLuminousBlueVariableMultiplier = false;
	sampleLuminousBlueVariableMultiplierMin = 1.0;
	sampleLuminousBlueVariableMultiplierMax = 12.0;

    useMCMC = false;

    useImportanceSampling = false;

    // Neutron star equation of state
    neutronStarEquationOfStateString = "SSE";
    neutronStarEquationOfState = 0;

}

void programOptions::setToFiducialValues()
{
	individualSystem = false;                                   // Flag to evolve a specific individual system which you can specify initial parameters of
    
    singleStar = false;                                         // Flag to evolve a single star

    // Debugging? default = false
    debugging = false;                                          // Flag to activate debugging statements
    RejectedSamplesPrinting = false;                           // Flag to turn on printing of the samples that are rejected by COMPAS      
    onlyDoubleCompactObjects = false;                           // Flag to turn on some shortcuts to only evolve systems which may form double compact objects

    // Individual system variables
    primaryMass = 96.2;                                         // Initial primary mass in solar masses
    secondaryMass = 60.2;                                       // Initial secondary mass in solar masses
    
    initialPrimaryMetallicity = 0.02;                           // Initial metallicity of the primary
    initialSecondaryMetallicity = 0.02;                         // Initial metallicity of the secondary
    
    // Variables required to restart a binary/star halfway through
    primaryStellarType = 1;                                    // Initial primary stellar type (not yet implemented)
    secondaryStellarType = 1;                                  // Initial secondary stellar type (not yet implemented)
    
    primaryEffectiveInitialMass = 0.;                  // Effective initial mass for the primary in solar masses
    secondaryEffectiveInitialMass = 0.;              // Effective initial mass for the secondary in solar masses
    
    primaryCoreMass = 0.0;                                      // Initial primary core mass in solar masses
    secondaryCoreMass = 0.0;                                    // Initial secondary core mass in solar masses
    
    primaryAge = 0.0;                                           // Effective age for the primary star in Myrs
    secondaryAge = 0.0;                                         // Effective age for the secondary star in Myrs
    
    binarySeparation = 11.5;                                    // Initial separation in AU
    binaryOrbitalPeriod = -1.0;                                 // Initial orbital period in day
    binaryEccentricity = 0.0;                                   // Initial eccentricity
    
    primaryRotationalVelocity = 0.0;                            // Initial rotational velocity of the primary
    secondaryRotationalVelocity = 0.0;                          // Initial rotational velocity of the secondary
    
    // Public population synthesis variables
    nBinaries = 10; 
    
    fixedRandomSeed = true;                                    // Whether to use a fixed random seed given by options.randomSeed (set to true if --random-seed is passed on command line) (default = false)
    randomSeed = 0;                                             // Random seed to use (default = 0)
    
    // Specify how long to evolve binaries for
    maxEvolutionTime = 13700.0;                                 // Maximum evolution time in Myrs
    maxNumberOfTimestepIterations = 99999;                      // Maximum number of timesteps to evolve binary for before giving up
    
    // Initial mass options
    initialMassFunction = "Kroupa";
    initialMassFunctionMin   = 8.;
    initialMassFunctionMax   = 100.;
    initialMassFunctionPower = 0.;
    
    // Initial mass ratios
    massRatioDistribution = "Flat";                             // Most likely want Flat or SANA2012
    massRatioDistributionMin = 0.0;
    massRatioDistributionMax = 1.0;
    
    minimumMassSecondary = 0.1;                                 // Minimum mass of secondary to draw (in Msol) (default = 0.1, brown dwarf limit)
    
    // Initial orbit options
    semiMajorAxisDistribution       = "FlatInLog";              // Most likely want FlatInLog or SANA2012
    semiMajorAxisDistributionMin    = 0.1;
    semiMajorAxisDistributionMax    = 1000.;
    semiMajorAxisDistributionPower  = -1;
    
    // Initial orbital period
    periodDistributionMin = 1.1;                                // Minimum initial period in days
    periodDistributionMax = 1000.0;                             // Maximum initial period in days
    
    // Eccentricity
    eccentricityDistribution = "Zero";                          // Most likely want ZERO, THERMAL or SANA2012
    eccentricityDistributionMin = 0.0;                          // Minimum initial eccentricity to sample
    eccentricityDistributionMax = 1.0;                          // Maximum initial eccentricity to sample
    
    // Kick options
    kickVelocityDistribution        		= "Maxwellian";     // Which kick velocity distribution to use
    //kickVelocityDistributionSigma   		= 250;				// Characteristic kick velocity in km s^-1 SIMON : REPLACED
    kickVelocityDistributionSigmaCCSN_NS    = 250;              // Kick velocity sigma in km s^-1 for neutron stars (default = "250" )
    kickVelocityDistributionSigmaCCSN_BH    = 250;              // Kick velocity sigma in km s^-1 for black holes (default = "250" )
    kickVelocityDistributionMaximum         = -1;               // Maximum kick velocity to draw in km s^-1. Ignored if < 0
    kickVelocityDistributionSigmaForECSN   	= 30.0;             // Characteristic kick velocity for an ECSN in km s^-1 (default = "30")
    kickVelocityDistributionSigmaForUSSN   	= 30.0;             // Characteristic kick velocity for an USSN in km s^-1 (default = "30")
	kickScalingFactor						= 1.0;				// Arbitrary factor for scaling kicks
    
    // Black hole kicks
    blackHoleKicks = "FALLBACK";
    blackHoleKicksOption = FALLBACK;
    
    // Supernova remnant mass prescription options
    // 
    // postitNote
    // hurley2000
    // belczynski2002
    // fryer2012
    //
    
    remnantMassPrescription = "FRYER2012";
    fryerSupernovaEngineString = "Delayed";
    fryerSupernovaEngine = SN_DELAYED;

    // Fixed uk options
    useFixedUK = false;
    fixedUK = -1.0;

    // Pair instability and pulsational pair instability mass loss
    usePairInstabilitySupernovae = true;                           // Whether to use pair instability supernovae (PISN)
    pairInstabilityUpperLimit = 135.0;                             // Maximum core mass leading to PISN (default = 135, Value in Belczynski+ 2016 is 135 Msol)
    pairInstabilityLowerLimit = 60.0;                              // Minimum core mass leading to PISN (default = 60, Value in Belczynski+ 2016 is 65 Msol)

    usePulsationalPairInstability = true;                          // Whether to use pulsational pair instability (PPI)
    pulsationalPairInstabilityLowerLimit = 35.0;                   // Minimum core mass leading to PPI (default = 40, Value in Belczynski+ 2016 is 45 Msol)
    pulsationalPairInstabilityUpperLimit = 60.0;                   // Maximum core mass leading to PPI (default = 60, Value in Belczynski+ 2016 is 65 Msol)

    pulsationalPairInstabilityPrescriptionString = "COMPAS";        // String for which PPI prescription to use
    pulsationalPairInstabilityPrescription = COMPASPPI;             // Prescription for PPI to use
    
	maximumNeutronStarMass = 3.0;								   // Maximum mass of a neutron star allowed, set to default in StarTrack
    nBatchesUsed = -1;                                              // nr of batches used, only needed for STROOPWAFEL (AIS) (default = -1, not needed) 
	
    // Kick direction option
    kickDirection = "Isotropic";                                    // Which assumption for SN kicks: Possibilities: Isotropic, poles, plane, power
    kickDirectionPower = 0.0;                                       // Power law power for the "power" SN kick direction choice
    
    // Get default output path
    outputPathString = ".";                                         // String to hold the output directory
    defaultOutputPath = boost::filesystem::current_path();          // Default output location
    outputPath = defaultOutputPath;                                 // Desired output location (default = CWD)

    // Spin options
    spinDistributionMin = 0.;
    spinDistributionMax = 1.;
    spinDistribution = "ZERO";
    spinAssumption = "bothAligned";
    
    // PN evolution options
    PNevolution = false;
    
    // Tides options
    tidesPrescriptionString = "None";                                   // String containing which tides prescription to use (default = "None")
    tidesPrescription = TIDES_PRESCRIPTION_NONE;                        // Tides prescription that will be used by the code
    // Tides prescription options
    //
    // NONE TIDES_PRESCRIPTION_NONE
    // LOCKED TIDES_PRESCRIPTION_LOCKED
    // HUT TIDES_PRESCRIPTION_HUT
    //
    
    // Magnetic Braking options
    //    useMagneticBraking = false;                                         // Wether to use Magnetic Braking
    
    // Mass loss options
    useMassLoss = true;                                                // Whether to use mass loss
    
    massLossPrescriptionString = "VINK";
    massLossPrescription = MASS_LOSS_PRESCRIPTION_VINK;
    // Mass loss prescription options
    // NONE     MASS_LOSS_PRESCRIPTION_NONE
    // HURLEY   MASS_LOSS_PRESCRIPTION_HURLEY
    // VINK     MASS_LOSS_PRESCRIPTION_VINK
    
    // Wind mass loss multiplicitive constants
    luminousBlueVariableFactor = 1.5;                                   // Luminous blue variable mass loss enhancement factor
    wolfRayetFactor = 1.0;                                              // WR winds factor
    
    // Mass transfer options
    useMassTransfer = true;											// Whether to use mass transfer (default = true)
	circulariseBinaryDuringMassTransfer	= true;						// Whether to circularise binary when it starts (default = true)
	forceCaseBBBCStabilityFlag = true;									// Whether if all case BB/BC systems are forced to be stable or unstable
	alwaysStableCaseBBBCFlag = true;									// Whether if case BB/BC is always stable
	angularMomentumConservationDuringCircularisation = true;		// Whether to conserve angular momentum while circularising or circularise to periastron (default = true)
    massTransferPrescriptionString = "DEMINK";
    massTransferPrescription = DEMINK_MASS_TRANSFER;
    // Mass transfer prescription options
    // DEMINK       DEMINK_MASS_TRANSFER
    // BELCZYNSKI   BELCZYNSKI_MASS_TRANSFER

    // Options adaptive Roche Lobe Overflow prescription
    massTransferAdaptiveAlphaParameter = 0.5;
    maxPercentageAdaptiveMassTransfer = 0.01;
    // Options for mass transfer accretion efficiency
    massTransferFractionAccreted = 1.0;
    massTransferCParameter = 10.0;
    massTransferAccretionEfficiencyPrescriptionString = "THERMAL";
    massTransferAccretionEfficiencyPrescription = THERMALLY_LIMITED_MASS_TRANSFER;
    // Mass transfer accretion efficiency prescription options
    // THERMAL          THERMALLY_LIMITED_MASS_TRANSFER
    // FIXED            FIXED_FRACTION_MASS_TRANSFER
    // CENTRIFUGAL      CENTRIFUGALLY_LIMITED_MASS_TRANSFER
	
	massTransferThermallyLimitedVariationString = "CFACTOR";
	massTransferThermallyLimitedVariation = THERMAL_C_FACTOR;
//	massTransferThermallyLimitedVariation = THERMAL_RADIUS_TO_ROCHELOBE;
	// CFACTOR			THERMAL_C_FACTOR
	// ROCHELOBE		THERMAL_RADIUS_TO_ROCHELOBE
	
    eddingtonAccretionFactor = 1;      // Multiplication factor for eddington accretion for NS & BH
                                         //i.e. >1 is super-eddington 
                                         //0. is no accretion   
    massTransferJloss= 1.0;
    massTransferAngularMomentumLossPrescriptionString = "ISOTROPIC";
    massTransferAngularMomentumLossPrescription = ISOTROPIC_RE_EMISSION_MASS_TRANSFER_ANGULAR_MOMENTUM_LOSS;
    // Mass transfer angular momentum loss prescription options
    // JEANS           JEANS_MASS_TRANSFER_ANGULAR_MOMENTUM_LOSS
    // ISOTROPIC       ISOTROPIC_RE_EMISSION_MASS_TRANSFER_ANGULAR_MOMENTUM_LOSS
    // CIRCUMBINARY    CIRCUMBINARY_RING_MASS_TRANSFER_ANGULAR_MOMENTUM_LOSS
    // ARBITRARY       ARBITRARY_MASS_TRANSFER_ANGULAR_MOMENTUM_LOSS
    
    massTransferRejuvenationPrescriptionString = "STARTRACK";
    massTransferRejuvenationPrescription = MASS_TRANSFER_REJUVENATION_STARTRACK;
    // Mass transfer rejuvenation prescriptions
    // "Belczynski weird stuff no 77"
    // Which prescription to use for rejuvenating stars when they accrete material during mass transfer
    // NONE             MASS_TRANSFER_REJUVENATION_NONE                 // The BSE prescription
    // STARTRACK        MASS_TRANSFER_REJUVENATION_STARTRACK            // The StarTrack prescription described in section 5.6 of Belczynski 2008 http://arxiv.org/pdf/astro-ph/0511811v3.pdf
    //
    
	    // Mass transfer critical mass ratios
    massTransferCriticalMassRatioMSLowMass = false;                    			// Whether to use critical mass ratios
    massTransferCriticalMassRatioMSLowMassNonDegenerateAccretor = 1.44;         // Critical mass ratio for MT from a MS low mass star (default = 1.44, Claeys+ 2014)
    massTransferCriticalMassRatioMSLowMassDegenerateAccretor = 1.0;             // Critical mass ratio for MT from a MS low mass star on to a degenerate accretor (default = 1.0, Claeys+ 2014)

    massTransferCriticalMassRatioMSHighMass = false;							// Whether to use critical mass ratios
    massTransferCriticalMassRatioMSHighMassNonDegenerateAccretor = 0.625;       // Critical mass ratio for MT from a MS high mass star (default = 0.625, Claeys+ 2014)
    massTransferCriticalMassRatioMSHighMassDegenerateAccretor = 0.0;            // Critical mass ratio for MT from a MS high mass star on to a degenerate accretor (default = 0)
    
    massTransferCriticalMassRatioHG = false;                                    // Whether to use critical mass ratios
    massTransferCriticalMassRatioHGNonDegenerateAccretor = 0.40;                // Critical mass ratio for MT from a HG star (default = 0.25, Claeys+ 2014)
    massTransferCriticalMassRatioHGDegenerateAccretor = 0.21;                   // Critical mass ratio for MT from a HG star on to a degenerate accretor (default = 0.21, Claeys+ 2014)
    
    massTransferCriticalMassRatioGiant = false;                                 // Whether to use critical mass ratios
    massTransferCriticalMassRatioGiantNonDegenerateAccretor = 0.0;              // Critical mass ratio for MT from a giant (default = 0.0)
    massTransferCriticalMassRatioGiantDegenerateAccretor = 0.87;                // Critical mass ratio for MT from a giant on to a degenerate accretor (default = 0.81, Claeys+ 2014)

    massTransferCriticalMassRatioHeliumMS = false;								// Whether to use critical mass ratios
    massTransferCriticalMassRatioHeliumMSNonDegenerateAccretor = 0.625;       	// Critical mass ratio for MT from a Helium MS star (default = 0.625)
    massTransferCriticalMassRatioHeliumMSDegenerateAccretor = 0.0;            	// Critical mass ratio for MT from a Helium MS star on to a degenerate accretor (default = 0)
	
    massTransferCriticalMassRatioHeliumHG = false;                              // Whether to use critical mass ratios
    massTransferCriticalMassRatioHeliumHGNonDegenerateAccretor = 0.25;			// Critical mass ratio for MT from a Helium HG star (default = 0.25, de Claeys+ 2014)
    massTransferCriticalMassRatioHeliumHGDegenerateAccretor = 0.21;				// Critical mass ratio for MT from a Helium HG star on to a degenerate accretor (default = 0.21, Claeys+ 2014)
    
    massTransferCriticalMassRatioHeliumGiant = false;                           // Whether to use critical mass ratios
    massTransferCriticalMassRatioHeliumGiantNonDegenerateAccretor = 1.28;       // Critical mass ratio for MT from a Helium giant (default = 0.25, Claeys+ 2014)
    massTransferCriticalMassRatioHeliumGiantDegenerateAccretor = 0.87;          // Critical mass ratio for MT from a Helium giant on to a degenerate accretor
	
    massTransferCriticalMassRatioWhiteDwarf = false;                           // Whether to use critical mass ratios
	massTransferCriticalMassRatioWhiteDwarfNonDegenerateAccretor = 0.0;       // Critical mass ratio for MT from a White Dwarf (default = 0.0)
    massTransferCriticalMassRatioWhiteDwarfDegenerateAccretor = 1.6;          // Critical mass ratio for MT from a White Dwarf on to a degenerate accretor (default = 1.6, Claeys+ 2014)
    
    // Common Envelope parameters
    commonEnvelopePrescriptionFlag = WEBBINK_COMMON_ENVELOPE;           // Which common envelope prescription to use
    commonEnvelopeAlpha = 1.0;                                          // Common envelope efficiency alpha parameter (default = 1.0)
    commonEnvelopeLambda = 0.1;                                         // Common envelope Lambda parameter (default = 0.1)
    commonEnvelopeHertzsprungGapDonorString = "PESSIMISTIC_HG_CE";      // String containing which prescription to use for Hertzsprung gap donors in a CE (default = "OPTIMISTIC_HG_CE")
    commonEnvelopeHertzsprungGapDonor = PESSIMISTIC_HG_CE;              // Which prescription to use for Hertzsprung gap donors in a CE (default = OPTIMISTIC_HG_CE)
    commonEnvelopeAlphaThermal = 1.0;                                   // lambda = (alpha_th * lambda_b) + (1-alpha_th) * lambda_g
    commonEnvelopeLambdaMultiplier = 1.0;                               // Multiply common envelope lambda by some constant
    allowMainSequenceStarToSurviveCommonEnvelope = false;               // Whether or not to allow a main sequence star to survive a common envelope event

    // Accretion during common envelope
    commonEnvelopeMassAccretionPrescriptionString = "ZERO";
    commonEnvelopeMassAccretionPrescription = COMMON_ENVELOPE_ACCRETION_ZERO;
    commonEnvelopeMassAccretionMin = 0.04; // Minimum amount of mass accreted during CE in solar masses
    commonEnvelopeMassAccretionMax = 0.1; // Maximum amount of mass accreted during CE in solar masses
    commonEnvelopeMassAccretionConstant = 0.0; // Constant value

	// Common envelope lambda prescription
	commonEnvelopeLambdaPrescriptionString = "LAMBDA_NANJING";					// String containing which prescription to use for CE lambda (default = "LAMBDA_NANJING")
	commonEnvelopeLambdaPrescription = LAMBDA_NANJING;								 // Which prescription to use for CE lambda (default = LAMBDA_NANJING)
	lambdaCalculationEveryTimeStep = false;
	
	// Common envelope Nandez and Ivanova energy formalism
	revisedEnergyFormalismNandezIvanova	= false;						// Use the revised energy formalism from Nandez & Ivanova 2016 (default = false)
	maximumMassDonorNandezIvanova = 2.0;								// Maximum mass allowed to use the revised energy formalism in Msol (default = 2.0)
	commonEnvelopeRecombinationEnergyDensity = 1.5E13;					// Factor using to calculate the bin
	
	// Common envelope power factor for Kruckow fit
	commonEnvelopeSlopeKruckow = -2.0/3.0;								// Common envelope power factor for Kruckow fit normalized according to Kruckow+2016, Fig. 1
		
	// Which prescription to use for calculating zetas
	commonEnvelopeZetaPrescriptionString = "STARTRACK";						// String containing which prescription to use for calculating CE zetas (default = STARTRACK)
	commonEnvelopeZetaPrescription = ZETA_STARTRACK;					// Which prescription to use for calculating CE zetas (default = ZETA_ADIABATIC)
	zetaCalculationEveryTimeStep = false;
	zetaAdiabaticArbitrary = NEVER_SET;
	zetaThermalArbitrary = NEVER_SET;
    zetaMainSequence 	= 6.5;
	zetaHertzsprungGap	= 2.0;
	
	BeBinaries = false;															// Flag if we want to print BeBinaries (main.cpp)	
	CHEvolution = false;
    RLOFPrinting = false;

    // Adaptive Importance Sampling Exploratory phase
    AISexploratoryphase = false;  // Floor 
    AISDCOtypeString = "ALL";                                               // String containing which type of DCOs to focus on (default = "ALL")
    AISDCOtype = ALL;                                                       // Which prescription to use for DCO type (default = ALL)                                            
    AISHubble = false;                                                      // Flag for excluding DCOs that do not merge in Hubble
    AISRLOF = false;                                                        // Flag for excluding DCOs that RLOFSecondaryZAMS
    AISPessimistic = false;                                                 // Flag for excluding DCOs that are Optmistic
    AISrefinementPhase = false;                                             // Flag for whether to run the AIS refinement phase (step 2)
    KappaGaussians = 2;                                                     // scaling factor for the width of the Gaussian distributions in AIS main sampling phase


    // Metallicity options
    metallicity = Zsol;
    fixedMetallicity = true;
    
    // Rotational velocity distribution options
    rotationalVelocityString = "ZERO";                  //
    
    // Pulsar birth magnetic field distribution
    pulsarBirthMagneticFieldDistributionString = "ZERO";        // Which birth magnetic field distribution to use for pulsars
    pulsarBirthMagneticFieldDistributionMin = 11.0;             // Minimum pulsar birth magnetic field distribution
    pulsarBirthMagneticFieldDistributionMax = 13.0;             // Maximum pulsar birth magnetic field distribution
    
    // Pulsar birth spin period distribution string
    pulsarBirthSpinPeriodDistributionString = "ZERO";           // Which birth spin period distribution to use for pulsars
    pulsarBirthSpinPeriodDistributionMin = 0.0;                 // Minimum birth spin period (ms)
    pulsarBirthSpinPeriodDistributionMax = 100.0;               // Maximum birth spin period (ms)
    
    pulsarMagneticFieldDecayTimescale = 1000.0;           // Timescale on which magnetic field decays (Myrs)
    pulsarMagneticFieldDecayMassscale = 0.025;            // Mass scale on which magnetic field decays during accretion (solar masses)  
    pulsarLog10MinimumMagneticField = 8.0;                // log10 of the minimum pulsar magnetic field in Gauss
    
    evolvePulsars = false;                                // Whether to evolve pulsars

    // Detailed output
    detailedOutput = false;
    
    // Print certain data for small populations, but not for larger one
    populationDataPrinting = false;

	//JIM BARRETT -- 06/07/2016 -- adding options to sample over some hyperparameters
	sampleKickVelocitySigma = false;
	sampleKickVelocitySigmaMin = 0.0;
	sampleKickVelocitySigmaMax = 400.0;
	
	sampleKickDirectionPower = false;
	sampleKickDirectionPowerMin = -10.0;
	sampleKickDirectionPowerMax = 10.0;
	
	sampleCommonEnvelopeAlpha = false;
	sampleCommonEnvelopeAlphaMin = 0.0;
	sampleCommonEnvelopeAlphaMax = 5.0;

	sampleWolfRayetMultiplier = false;
	sampleWolfRayetMultiplierMin = 0.0;
	sampleWolfRayetMultiplierMax = 5.0;

	sampleLuminousBlueVariableMultiplier = false;
	sampleLuminousBlueVariableMultiplierMin = 1.0;
	sampleLuminousBlueVariableMultiplierMax = 12.0;
    
    useMCMC = false;

    useImportanceSampling = false;
    
    // Neutron star equation of state
    neutronStarEquationOfStateString = "SSE";
    neutronStarEquationOfState = 0;
    
}
