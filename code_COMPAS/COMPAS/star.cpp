//
//  star.cpp
//  singleStarEvolution

#include "star.h"

// Static arrays containing the numbers for ZAMS, SSE

// L coefficients
// Table 1 in Tout et al 1996
static const int nrowsL    = 7;
static const int ncolsL    = 5;
static const double Lcoeff[nrowsL][ncolsL] = {
    {0.39704170,	-0.32913574,	0.34776688,     0.37470851,     0.09011915},
    {8.52762600,	-24.41225973,	56.43597107,	37.06152575,	5.45624060},
    {0.00025546,	-0.00123461,	-0.00023246,	0.00045519,     0.00016176},
    {5.43288900,	-8.62157806,	13.44202049,	14.51584135,	3.39793084},
    {5.56357900,	-10.32345224,	19.44322980,	18.97361347,	4.16903097},
    {0.78866060,	-2.90870942,	6.54713531,     4.05606657,     0.53287322},
    {0.00586685,	-0.01704237,	0.03872348,     0.02570041,     0.00383376}
};


// R coefficients
// Table 2 in Tout et al 1996
static const int nrowsR     = 9;
static const int ncolsR     = ncolsL;
static const double Rcoeff[nrowsR][ncolsR] = {
    {1.71535900,	0.62246212,     -0.92557761,	-1.16996966,	-0.30631491},
    {6.59778800,	-0.42450044,	-12.13339427,	-10.73509484,	-2.51487077},
    {10.08855000,	-7.11727086,	-31.67119479,	-24.24848322,	-5.33608972},
    {1.01249500,	0.32699690,     -0.00923418,	-0.03876858,	-0.00412750},
    {0.07490166,	0.02410413,     0.07233664,     0.03040467,     0.00197741},
    {0.01077422,	0.00000000,     0.00000000,     0.00000000,     0.00000000},
    {3.08223400,	0.94472050,     -2.15200882,	-2.49219496,	-0.63848738},
    {17.84778000,	-7.45345690,	-48.96066856,	-40.05386135,	-9.09331816},
    {0.00022582,	-0.00186899,	0.00388783,     0.00142402,     -0.00007671}
};


// Default constructor -- probably not used a lot except for testing
Star::Star(){

	// m_randomSeed = NEVER_SET;

    // Initialise age to 0
    m_time                  = 0.0;                          // Physical time the star has been evolved for
    m_Age                   = 0.0;                          // Effective age of the star, differs for different phases of evolution
    m_dt                    = 0.0;                          // Current timestep
    m_tau                   = 0.0;                          // Relative time
    m_nuclearBurningTime    = 0.0;                          // Nuclear burning time
    m_dtWind                = 0.0;                          // Timestep to check for wind-mass-loss cap
    
    // Calculate the metallicity dependent numbers
    m_Metallicity           = 0.01;
    m_logMetallicityXi      = log10(m_Metallicity / Zsol);
    m_logMetallicitySigma   = log10(m_Metallicity);
    m_logMetallicityRho     = m_logMetallicityXi + 1.0;
    
    // Store metallicities in an array as well?

    // Calculate metallicity dependent radius coefficients
    calculateAllRadiusCoefficients(m_logMetallicityXi, radius_coefficient_theta, radius_coefficient_iota, radius_coefficient_kappa, radius_coefficient_lamda, radius_coefficient_mu, radius_coefficient_nu, radius_coefficient_xi, radius_coefficient_omicron, radius_coefficient_Pi);

    // Calculate metallicity dependent luminosity coefficients
    calculateAllLuminosityCoefficients(m_logMetallicityXi, luminosity_coefficient_alpha, luminosity_coefficient_beta, luminosity_coefficient_gamma, luminosity_coefficient_delta, luminosity_coefficient_epsilon, luminosity_coefficient_zeta, luminosity_coefficient_eta);
    
    // Calculate the metallicity dependent numbers for ZAMS -- HOW DO YOU PLAN TO GIVE THEM ACCESS TO THE DATA FILES?
    // import them somewhere else and pass them around or
    // have them hardcoded as static arrays that just get generated once and are nice and easy to use via indexes?
    m_MHeF          = calcMHeF(m_logMetallicityXi);
    m_Mhook         = calcMhook(m_logMetallicityXi);
    m_MFGB          = calcMFGB(m_Metallicity);
    
    m_massCutoffs[0] = m_Mhook;
    m_massCutoffs[1] = m_MHeF;
    m_massCutoffs[2] = m_MFGB;
    
    // Initialise an, bn, GB arrays
    for(int i = 0; i < nrowsa; i++){
        m_an_coefficients[i] = 0.0;
    }
    
    for(int i = 0; i < nrowsb; i++){
        m_bn_coefficients[i] = 0.0;
    }
    
    for(int i = 0; i < nGBParams; i++){
        m_GBParams[i] = 0.0;
    }
    
    // Calculate metallicity dependent coefficients
    calculate_aCoefficients(m_an_coefficients, m_Metallicity);  // also needs logMetallicity right?
    calculate_bCoefficients(m_bn_coefficients, m_Metallicity, m_massCutoffs);  // Same here
    
    // Calculate the timescales
    m_t_BGB             = lifetimeToBGB(m_an_coefficients, m_MZAMS);
    m_t_hook            = lifetimeHook(m_an_coefficients, m_MZAMS, m_t_BGB);    // calculate t_hook
    m_t_MS              = lifetimeMainSequence(m_an_coefficients, m_MZAMS, m_logMetallicityXi, m_t_BGB, m_t_hook);  // MS Lifetime
 //   m_t_MS1			= m_t_MS; //for public access

    // Can you just calculate these once, or do they have to be recalculated due to mass loss
    // What about the other timescales?
    // Should be done in an update timescales function?
    
    // Initialise timescales
    for(int i = 0; i < nTimescales; i++){
        m_timescales[i] = 0.0;
    }
    
    // Again, store these in an easily accessible list/array that I have to remember/write the indexes for somewhere (python dictionaries would be wonderful for this :()
    
    // Calculate Zero Age Main Sequence (ZAMS) Radius and Luminosity using Pols et al formulae
    m_MZAMS             = 1.0;
    // m_RZAMS             = RadiusZAMS(m_MZAMS, m_logMetallicityXi);
    m_RZAMS             = RadiusZAMS(m_MZAMS, radius_coefficient_theta, radius_coefficient_iota, radius_coefficient_kappa, radius_coefficient_lamda, radius_coefficient_mu, radius_coefficient_nu, radius_coefficient_xi, radius_coefficient_omicron, radius_coefficient_Pi);
    //m_LZAMS             = LuminosityZAMS(m_MZAMS, m_logMetallicityXi);
    m_LZAMS             = LuminosityZAMS(m_MZAMS, luminosity_coefficient_alpha, luminosity_coefficient_beta, luminosity_coefficient_gamma, luminosity_coefficient_delta, luminosity_coefficient_epsilon, luminosity_coefficient_zeta, luminosity_coefficient_eta);
    m_TZAMS             = calculateTemperature(m_LZAMS, m_RZAMS);
    m_omegaZAMS         = 0.0;

    // m_omegaZAMS
    //rotationalAngularVelocityDistribution(m_MZAMS, m_RZAMS); 
    //rotationalAngularVelocityDistribution(m_MZAMS, m_RZAMS, options, r); 
    //0.0;
    //angularFrequencyZAMS(m_MZAMS,m_RZAMS);
    
    // Variables for the effective initial Zero Age Main Sequence parameters corresponding to Mass0
    m_RZAMS0            = m_RZAMS;
    m_LZAMS0            = m_LZAMS;
    m_TZAMS0            = m_TZAMS;
    m_omegaZAMS0        = m_omegaZAMS; 

    //rotationalAngularVelocityDistribution(m_MZAMS, m_RZAMS, options, r); 
    //0.0;
    //angularFrequencyZAMS(m_MZAMS,m_RZAMS);

    // Current timestep variables
    m_Mass              = m_MZAMS;      // Current mass
    m_Mass0             = m_MZAMS;      // Current Effective initial mass
    m_massLossDiff      = 0.0;          // Amount of mass lost due to winds
    m_Luminosity        = m_LZAMS;      // Current Luminosity
    m_Radius            = m_RZAMS;      // Current Radius
    m_Temperature       = m_TZAMS;      // Current Temperature
    m_coreMass          = 0.0;          // Current Core Mass
    m_HeCoreMass        = 0.0;          // Current He Core Mass
    m_COCoreMass        = 0.0;          // Current CO Core Mass
    m_envMass           = initialEnvelopeMass(m_Mass);       // Current envelope mass
    m_Mu                = 0.0;          // Current Mu
    m_coreRadius        = 0.0;          // Current core radius

    m_omegaBreak        = omegaBreak(m_Mass, m_Radius); // sqrt(pow(2.0*pi,2.0)*(m_MZAMS)/pow(RsolToAU*m_RZAMS,3.0)); // SIMON: Changed to use function
    m_omega             = m_omegaZAMS;

    m_omegaPrev         = m_omegaZAMS;

    m_momentInertia     = 0.0; //momentOfInertia3(m_Mass,  m_coreMass,  m_Radius, 0.0,
                                 //          m_stellarType);
    m_angularMomentum   = 0.0;//rotationalAngularMomentum(m_momentInertia, m_omega);

	m_radialExpansionTimescale = 0.0;
	
    // Previous timestep variables
    m_MassPrev              = m_MZAMS;      // Current mass
    m_Mass0Prev             = m_MZAMS;      // Current Effective initial mass
    m_massLossDiffPrev      = 0.0;          // Amount of mass lost due to winds
    m_LuminosityPrev        = m_LZAMS;      // Current Luminosity
    m_RadiusPrev            = m_RZAMS;      // Current Radius
    m_TemperaturePrev       = m_TZAMS;      // Current Temperature
    m_coreMassPrev          = 0.0;          // Current Core Mass
    m_HeCoreMassPrev        = 0.0;          // Current He Core Mass
    m_COCoreMassPrev        = 0.0;          // Current CO Core Mass
    m_envMassPrev           = initialEnvelopeMass(m_Mass);       // Current envelope mass
    m_MuPrev                = 0.0;          // Current Mu
    m_coreRadiusPrev        = 0.0;          // Current core radius

    //
    m_omegaTidesIndividualDiff = 0.0;
    m_omegaMagneticBrakingDiff = 0.0;
    m_MassTransferDiff = 0.0;
    m_massLossDiff      = 0.0;
    m_omegaTidesIndividualDiffPrev = 0.0;
    m_omegaMagneticBrakingDiffPrev = 0.0;
    m_MassTransferDiffPrev = 0.0;
    m_massLossDiffPrev      = 0.0;
    
	// Binding energies
	m_NanjingLambda			= 0.0;
	m_LoveridgeLambda		= 0.0;
	m_LoveridgeWindsLambda	= 0.0;
	m_KruckowTopLambda		= 0.0;
	m_KruckowMidLambda		= 0.0;
	m_KruckowBotLambda		= 0.0;
	m_DewiLambda			= 0.0;
	
	// Logarithmic mass-radius derivatives: dlogR/dlogM
	m_zetaThermal	= 0.0;
	m_zetaNuclear	= 0.0;
	m_zetaSoberman	= 0.0;
	m_zetaSobermanHelium = 0.0;
	m_zetaHurley	= 0.0;
	m_zetaHurleyHelium = 0.0;
	m_zetaSimple	= 0.0;	
	
	// Stellar timescales
	m_dynamicalTimescale	= 0.0;
	m_thermalTimescale		= 0.0;
	m_nuclearTimescale		= 0.0;
	
    // Initialise pulsar parameters
    m_pulsarMagneticField = NEVER_SET;                    // Pulsar magnetic field strength (G)
    m_pulsarSpinPeriod = NEVER_SET;                       // Pulsar spin period (ms)
    m_pulsarSpinFrequency = NEVER_SET;		          // Pulsar spin frequency (rads per second)
    m_pulsarSpinDownRate = NEVER_SET;                     // Pulsar spin down rate (Pdot, dimensionless)

	// Supernova Flags
	m_flagHrichSN		 = false;
	m_flagHpoorSN		 = false;
	flagUSSN			 = false;
    flagSNPrev  		 = false;
    flagExperiencedCCSN  = false;
	flagExperiencedECSN  = false;
    flagExperiencedPISN  = false;                            // Flag to know if at any point in time the star experienced a pair instability supernova
    flagExperiencedPPISN = false;                           // Flag to know if at any point in time the star experienced a pulsational pair instability supernova
	m_runawayFlag		 = false;
    m_fallback			 = 0.0;
	
	// Flags
    flagRLOF    	= false;
	experiencedRLOF = false;
    flagInitiateMT 	= false;
    flagSN      	= false;
    flagPISN        = false;
    flagPPISN       = false;
    flagECSN    	= false;
    m_error     	= false;
	m_LBVphaseFlag	= false;
    flagRLOFPrev	= false;
    m_errorPrev 	= false;
    m_fastPhaseCaseA    = false;
    m_firstMassTransferEpisode = false;
    m_initialMassTransferCase = NO_MASS_TRANSFER;
	flagRecycledNS	= false;
	flagRLOFontoaNS	= false;
    flagCHEvolution = false;
    
    m_Period = NEVER_SET;
    m_preSNeOrbitalEnergy=NEVER_SET;
    m_postSNeOrbitalEnergy=NEVER_SET;
    
    
    // Assign stellar stype
    if(m_MZAMS <= 0.7){
        m_stellarType = MS_LESS_THAN_07;
    }
    else{
        m_stellarType = MS_MORE_THAN_07;
    }

    m_stellarTypePrev = m_stellarType;
    
    // Options
    m_remnantMassPrescription = "fryer2012"; // Since no options have been provided.
    m_SNengine = SN_RAPID;
    
    m_totalMassAtCompactObjectFormation		= NEVER_SET;
    m_HeCoreMassAtCompactObjectFormation 	= NEVER_SET;
    m_COCoreMassAtCompactObjectFormation 	= NEVER_SET;
    m_coreMassAtCompactObjectFormation 		= NEVER_SET;        // What was the core mass of this star when it formed a compact object
    m_HeCoreMassAtCommonEnvelope 			= NEVER_SET;		// He mass at CE
	m_COCoreMassAtCommonEnvelope 			= NEVER_SET;		// CO mass at CE
    m_coreMassAtCommonEnvelope 				= NEVER_SET;        // Core mass at CE
	m_lambdaAtCommonEnvelope				= NEVER_SET;		// Lambda value used at CE
	m_bindingEnergyAtCommonEnvelope			= NEVER_SET;		// Binding energy value used during CE
    m_drawnKickVelocity 					= NEVER_SET;		// Kick velocity the system received during the supernova (km s^-1)
    m_kickVelocity 							= NEVER_SET;        // Kick velocity the system received during the supernova (km s^-1)
    m_supernovaPsi							= NEVER_SET;        // True anomaly (rad)
	
    m_supernovaTheta                        = 0.0;        // Angle between the orbital plane and the 'z' axis of supernovae vector (rad)
    m_supernovaPhi                          = 0.0;        // Angle between 'x' and 'y', both in the orbital plane of supernovae vector (rad)    

    //kickDirectionDistribution(options, r, m_supernovaTheta, m_supernovaPhi, options.kickDirectionPower);

    // Random number U(0,1) for choosing the supernova kick velocity magnitude
    //m_supernovaKickVelocityMagnitudeRandomNumber = uniformDistribution(r, 0.0, 1.0);        
    m_supernovaKickVelocityMagnitudeRandomNumber = 0.0;

    // Orbital anomalies
    //m_MeanAnomaly       = uniformDistribution(r, 0, 2*M_PI);    // Mean anomaly at instantaneous time of the SN - uniform in [0, 2pi]
    m_MeanAnomaly       = 0;                                    // Mean anomaly at instantaneous time of the SN - uniform in [0, 2pi]
    m_EccentricAnomaly  = 0;                                    // Eccentric anomaly at instataneous time of the SN
    m_TrueAnomaly       = 0;                                    // True anomaly at instantaneous time of the SN

    // booleans for fiddling with bits of the code whilst debugging.
    debugging = false;                                          // Whether to print debugging statements
    isPerturbation = true;                                      // Whether to use perturbations (disable for testing)
    isMassLoss = true;                                          // Whether to use mass loss (disable for testing)
    massLossPrescription = MASS_LOSS_PRESCRIPTION_HURLEY;       // Which mass loss prescription to use
    writeOutput = false;                                        // Whether to write parameters after timestep to standard output
    
}

// Default constructor -- probably not used a lot except for testing
Star::Star(double MZAMS, double Metallicity){
	
	// ALEJANDRO - 18/11/2016 - Adding random seed to keep tracks of errors in star.cpp while running a population. This seed should be the same as the currentBinary.m_randomSeed seed.
	// m_randomSeed = NEVER_SET;
    
    // Initialise age to 0
    m_time                  = 0.0;                          // Physical time the star has been evolved for
    m_Age                   = 0.0;                          // Effective age of the star, differs for different phases of evolution
    m_dt                    = 0.0;                          // Current timestep
    m_tau                   = 0.0;                          // Relative time
    m_nuclearBurningTime    = 0.0;                          // Nuclear burning time
    m_dtWind                = 0.0;                          // Timestep to check for wind-mass-loss cap    

    // Calculate the metallicity dependent numbers
    m_Metallicity           = Metallicity;
    m_logMetallicityXi      = log10(m_Metallicity / Zsol);
    m_logMetallicitySigma   = log10(m_Metallicity);
    m_logMetallicityRho     = m_logMetallicityXi + 1.0;
    
    // Store metallicities in an array as well?

    // Calculate metallicity dependent radius coefficients
    calculateAllRadiusCoefficients(m_logMetallicityXi, radius_coefficient_theta, radius_coefficient_iota, radius_coefficient_kappa, radius_coefficient_lamda, radius_coefficient_mu, radius_coefficient_nu, radius_coefficient_xi, radius_coefficient_omicron, radius_coefficient_Pi);

    // Calculate metallicity dependent luminosity coefficients
    calculateAllLuminosityCoefficients(m_logMetallicityXi, luminosity_coefficient_alpha, luminosity_coefficient_beta, luminosity_coefficient_gamma, luminosity_coefficient_delta, luminosity_coefficient_epsilon, luminosity_coefficient_zeta, luminosity_coefficient_eta);
    
    // Calculate the metallicity dependent numbers for ZAMS -- HOW DO YOU PLAN TO GIVE THEM ACCESS TO THE DATA FILES?
    // import them somewhere else and pass them around or
    // have them hardcoded as static arrays that just get generated once and are nice and easy to use via indexes?
    m_MHeF          = calcMHeF(m_logMetallicityXi);
    m_Mhook         = calcMhook(m_logMetallicityXi);
    m_MFGB          = calcMFGB(m_Metallicity);
    
    //m_massCutoffs   = {m_MHeF, m_Mhook, m_MFGB}; // difficult to assign to already initialised array -- use std::vector?
    //or
    m_massCutoffs[0] = m_Mhook;
    m_massCutoffs[1] = m_MHeF;
    m_massCutoffs[2] = m_MFGB;
    
    // Initialise an, bn, GB arrays
    for(int i = 0; i < nrowsa; i++){
        m_an_coefficients[i] = 0.0;
    }
    
    for(int i = 0; i < nrowsb; i++){
        m_bn_coefficients[i] = 0.0;
    }
    
    for(int i = 0; i < nGBParams; i++){
        m_GBParams[i] = 0.0;
    }
    
    // Calculate metallicity dependent coefficients
    calculate_aCoefficients(m_an_coefficients, m_Metallicity);  // also needs logMetallicity right?
    calculate_bCoefficients(m_bn_coefficients, m_Metallicity, m_massCutoffs);  // Same here
    
    // Calculate the timescales
    m_t_BGB             = lifetimeToBGB(m_an_coefficients, m_MZAMS);
    m_t_hook            = lifetimeHook(m_an_coefficients, m_MZAMS, m_t_BGB);    // calculate t_hook
    m_t_MS              = lifetimeMainSequence(m_an_coefficients, m_MZAMS, m_logMetallicityXi, m_t_BGB, m_t_hook);  // MS Lifetime
    // Can you just calculate these once, or do they have to be recalculated due to mass loss
    // What about the other timescales?
    // Should be done in an update timescales function?
    
    // Initialise timescales
    for(int i = 0; i < nTimescales; i++){
        m_timescales[i] = 0.0;
    }
    
    // Again, store these in an easily accessible list/array that I have to remember/write the indexes for somewhere (python dictionaries would be wonderful for this :()
    
    // Calculate Zero Age Main Sequence (ZAMS) Radius and Luminosity using Pols et al formulae
    m_MZAMS             = MZAMS;
    // m_RZAMS             = RadiusZAMS(m_MZAMS, m_logMetallicityXi);
    m_RZAMS             = RadiusZAMS(m_MZAMS, radius_coefficient_theta, radius_coefficient_iota, radius_coefficient_kappa, radius_coefficient_lamda, radius_coefficient_mu, radius_coefficient_nu, radius_coefficient_xi, radius_coefficient_omicron, radius_coefficient_Pi);
    //m_LZAMS             = LuminosityZAMS(m_MZAMS, m_logMetallicityXi);
    m_LZAMS             = LuminosityZAMS(m_MZAMS, luminosity_coefficient_alpha, luminosity_coefficient_beta, luminosity_coefficient_gamma, luminosity_coefficient_delta, luminosity_coefficient_epsilon, luminosity_coefficient_zeta, luminosity_coefficient_eta);
    m_TZAMS             = calculateTemperature(m_LZAMS, m_RZAMS);
    m_omegaZAMS         = 0.0; 

    //rotationalAngularVelocityDistribution(m_MZAMS, m_RZAMS); 
    //rotationalAngularVelocityDistribution(m_MZAMS, m_RZAMS, options, r); 
    //0.0; 
    //angularFrequencyZAMS(m_MZAMS,m_RZAMS);

    // Variables for the effective initial Zero Age Main Sequence parameters corresponding to Mass0
    m_RZAMS0            = m_RZAMS;
    m_LZAMS0            = m_LZAMS;
    m_TZAMS0            = m_TZAMS;
    m_omegaZAMS0        = m_omegaZAMS; 

    //rotationalAngularVelocityDistribution(m_MZAMS, m_RZAMS, options, r); 
    //0.0;
    //angularFrequencyZAMS(m_MZAMS,m_RZAMS);

    // Current timestep variables
    m_Mass              = m_MZAMS;      // Current mass
    m_Mass0             = m_MZAMS;      // Current Effective initial mass
    m_Luminosity        = m_LZAMS;      // Current Luminosity
    m_Radius            = m_RZAMS;      // Current Radius
    m_Temperature       = m_TZAMS;      // Current Temperature
    m_coreMass          = 0.0;          // Current Core Mass
    m_HeCoreMass        = 0.0;          // Current He Core Mass
    m_COCoreMass        = 0.0;          // Current CO Core Mass
    m_envMass           = initialEnvelopeMass(m_Mass);       // Current envelope mass
    m_Mu                = 0.0;          // Current Mu
    m_coreRadius        = 0.0;          // Current core radius

    m_omegaBreak        = omegaBreak(m_Mass, m_Radius); //sqrt(pow(2.0*pi,2.0)*(m_MZAMS)/pow(m_RZAMS,3.0)); // SIMON: Changed to use function
    m_omega             = m_omegaZAMS;
    m_omegaPrev         = m_omegaZAMS;

    m_momentInertia     = 0.0; //momentOfInertia3(m_Mass,  m_coreMass,  m_Radius, 0.0,
                                 //          m_stellarType);
    m_angularMomentum   = 0.0; //rotationalAngularMomentum(m_momentInertia, m_omega);

	m_radialExpansionTimescale = 0.0;
	
    // Previous timestep variables
    m_MassPrev              = m_MZAMS;      // Current mass
    m_Mass0Prev             = m_MZAMS;      // Current Effective initial mass
    m_massLossDiffPrev      = 0.0;          // Amount of mass lost due to winds
    m_LuminosityPrev        = m_LZAMS;      // Current Luminosity
    m_RadiusPrev            = m_RZAMS;      // Current Radius
    m_TemperaturePrev       = m_TZAMS;      // Current Temperature
    m_coreMassPrev          = 0.0;          // Current Core Mass
    m_HeCoreMassPrev        = 0.0;          // Current He Core Mass
    m_COCoreMassPrev        = 0.0;          // Current CO Core Mass
    m_envMassPrev           = initialEnvelopeMass(m_Mass);       // Current envelope mass
    m_MuPrev                = 0.0;          // Current Mu
    m_coreRadiusPrev        = 0.0;          // Current core radius

    //
    m_omegaTidesIndividualDiff = 0.0;
    m_omegaMagneticBrakingDiff = 0.0;
    m_MassTransferDiff = 0.0;
    m_massLossDiff      = 0.0;
    m_omegaTidesIndividualDiffPrev = 0.0;
    m_omegaMagneticBrakingDiffPrev = 0.0;
    m_MassTransferDiffPrev = 0.0;
    m_massLossDiffPrev      = 0.0;
	
	// Binding energies
	m_NanjingLambda			= NEVER_SET;
	m_LoveridgeLambda		= NEVER_SET;
	m_LoveridgeWindsLambda	= NEVER_SET;
	m_KruckowTopLambda		= NEVER_SET;
	m_KruckowMidLambda		= NEVER_SET;
	m_KruckowBotLambda		= NEVER_SET;
	m_DewiLambda			= NEVER_SET;
	
	// Logarithmic mass-radius derivatives: dlogR/dlogM
	m_zetaThermal	= NEVER_SET;
	m_zetaNuclear	= NEVER_SET;
	m_zetaSoberman	= NEVER_SET;
	m_zetaSobermanHelium = NEVER_SET;
	m_zetaHurley	= NEVER_SET;
	m_zetaHurleyHelium = NEVER_SET;
	m_zetaSimple	= NEVER_SET;

	// Stellar timescales
	m_dynamicalTimescale	= NEVER_SET;
	m_thermalTimescale		= NEVER_SET;
	m_nuclearTimescale		= NEVER_SET;

    // Initialise pulsar parameters
    m_pulsarMagneticField = NEVER_SET;                    // Pulsar magnetic field strength (G)
    m_pulsarSpinPeriod = NEVER_SET;                       // Pulsar spin period (ms)
    m_pulsarSpinFrequency = NEVER_SET;                    // Pulsar spin frequency (rads per second)    
    m_pulsarSpinDownRate = NEVER_SET;                     // Pulsar spin down rate (Pdot, dimensionless)

	// Supernova Flags
	m_flagHrichSN		= false;
	m_flagHpoorSN		= false;
	flagUSSN			= false;
    flagSNPrev  		= false;
    flagExperiencedCCSN = false;
	flagExperiencedECSN = false;
    flagExperiencedPISN  = false;                            // Flag to know if at any point in time the star experienced a pair instability supernova
    flagExperiencedPPISN = false;                           // Flag to know if at any point in time the star experienced a pulsational pair instability supernova
	m_runawayFlag		= false;
    m_fallback          = NEVER_SET;
	
	// Flags
    flagRLOF    	= false;
	experiencedRLOF = false;
    flagInitiateMT 	= false;
    flagSN      	= false;
    flagPISN        = false;
    flagPPISN       = false;
    flagECSN    	= false;
    m_error     	= false;
	m_LBVphaseFlag	= false;
    flagRLOFPrev	= false;
	m_errorPrev 	= false;
    m_fastPhaseCaseA    = false;
    m_firstMassTransferEpisode = false;
    m_initialMassTransferCase = NO_MASS_TRANSFER;
	flagRecycledNS	= false;
	flagRLOFontoaNS	= false;
    flagCHEvolution = false;
	
    m_Period = NEVER_SET;
    m_preSNeOrbitalEnergy=NEVER_SET;
    m_postSNeOrbitalEnergy=NEVER_SET;

    
    // Assign stellar stype
    if(m_MZAMS <= 0.7){
        m_stellarType = MS_LESS_THAN_07;
    }
    else{
        m_stellarType = MS_MORE_THAN_07;
    }
    
    m_stellarTypePrev = m_stellarType;
    
    // options not provided
    m_remnantMassPrescription = "fryer2012"; // Since no options have been provided.
    m_SNengine = SN_RAPID;

    m_totalMassAtCompactObjectFormation		= NEVER_SET;
    m_HeCoreMassAtCompactObjectFormation 	= NEVER_SET;
    m_COCoreMassAtCompactObjectFormation 	= NEVER_SET;
    m_coreMassAtCompactObjectFormation 		= NEVER_SET;        // What was the core mass of this star when it formed a compact object
    m_HeCoreMassAtCommonEnvelope 			= NEVER_SET;		// He mass at CE
    m_COCoreMassAtCommonEnvelope 			= NEVER_SET;		// CO mass at CE
    m_coreMassAtCommonEnvelope 				= NEVER_SET;        // Core mass at CE
	m_lambdaAtCommonEnvelope				= NEVER_SET;		// Lambda value used at CE
	m_bindingEnergyAtCommonEnvelope			= NEVER_SET;		// Binding energy value used during CE
    m_drawnKickVelocity 					= NEVER_SET;        // Kick velocity the system received during the supernova (km s^-1)
    m_kickVelocity 							= NEVER_SET;        // Kick velocity the system received during the supernova (km s^-1)
    m_supernovaPsi							= NEVER_SET;        // True anomaly (rad)
	
    m_supernovaTheta                        = NEVER_SET;        // Angle between the orbital plane and the 'z' axis of supernovae vector (rad)
    m_supernovaPhi                          = NEVER_SET;        // Angle between 'x' and 'y', both in the orbital plane of supernovae vector (rad)    

    //kickDirectionDistribution(options, r, m_supernovaTheta, m_supernovaPhi, options.kickDirectionPower);

    // Random number U(0,1) for choosing the supernova kick velocity magnitude
    //m_supernovaKickVelocityMagnitudeRandomNumber = uniformDistribution(r, 0.0, 1.0);     
    m_supernovaKickVelocityMagnitudeRandomNumber = 0.0;    
    
    // Orbital anomalies
    //m_MeanAnomaly       = uniformDistribution(r, 0, 2*M_PI);    // Mean anomaly at instantaneous time of the SN - uniform in [0, 2pi]
    m_MeanAnomaly       = 0.0;                                  // Mean anomaly at instantaneous time of the SN - uniform in [0, 2pi]
    m_EccentricAnomaly  = 0;                                    // Eccentric anomaly at instataneous time of the SN
    m_TrueAnomaly       = 0;                                    // True anomaly at instantaneous time of the SN

    // booleans for fiddling with bits of the code whilst debugging.
    debugging = false;                                          // Whether to print debugging statements
    isPerturbation = true;                                      // Whether to use perturbations (disable for testing)
    isMassLoss = true;                                          // Whether to use mass loss (disable for testing)
    massLossPrescription = MASS_LOSS_PRESCRIPTION_HURLEY;       // Which mass loss prescription to use
    writeOutput = false;                                        // Whether to write parameters after timestep to standard output
    
}

// Have a constructor for a given mass, metallicity (most commonly used version) and options
Star::Star(double MZAMS, double Metallicity, const programOptions &options, const gsl_rng *r){ 
    
	// ALEJANDRO - 18/11/2016 - Adding random seed to keep tracks of errors in star.cpp while running a population. This seed should be the same as the currentBinary.m_randomSeed seed.
	// m_randomSeed = NEVER_SET;
	
    // Initialise age etc to 0
    m_time                  = 0.0;                  // Physical time the star has been evolved for
    m_Age                   = 0.0;                  // Effective age of the star, differs for different phases of evolution
    m_nuclearBurningTime    = 0.0;                  // Nuclear burning time
    m_dt                    = 0.0;                  // Current timestep
    m_tau                   = 0.0;                  // Relative time
    m_dtWind                = 0.0;                  // Timestep to check for wind-mass-loss cap    

    // Assign stellar stype
    if(MZAMS <= 0.7){
        m_stellarType = MS_LESS_THAN_07;
    }
    else{
        m_stellarType = MS_MORE_THAN_07;
    }
    
    m_stellarTypePrev = m_stellarType;
    
    // Calculate the metallicity dependent numbers
    m_Metallicity           = Metallicity;
    m_logMetallicityXi      = log10(m_Metallicity / Zsol);
    m_logMetallicitySigma   = log10(m_Metallicity);
    m_logMetallicityRho     = m_logMetallicityXi + 1.0;
    
    // Calculate metallicity dependent radius coefficients
    calculateAllRadiusCoefficients(m_logMetallicityXi, radius_coefficient_theta, radius_coefficient_iota, radius_coefficient_kappa, radius_coefficient_lamda, radius_coefficient_mu, radius_coefficient_nu, radius_coefficient_xi, radius_coefficient_omicron, radius_coefficient_Pi);

    // Calculate metallicity dependent luminosity coefficients
    calculateAllLuminosityCoefficients(m_logMetallicityXi, luminosity_coefficient_alpha, luminosity_coefficient_beta, luminosity_coefficient_gamma, luminosity_coefficient_delta, luminosity_coefficient_epsilon, luminosity_coefficient_zeta, luminosity_coefficient_eta);
    
    // Store metallicities in an array as well?
    m_Mhook   = calcMhook(m_logMetallicityXi);
    m_MHeF    = calcMHeF(m_logMetallicityXi);
    m_MFGB    = calcMFGB(m_Metallicity);
    
    // Do you really need variables, and then put them into the (makes it a bit clearer what they are I guess?)
    m_massCutoffs[0] = m_Mhook;
    m_massCutoffs[1] = m_MHeF;
    m_massCutoffs[2] = m_MFGB;
    
    // Initialise an, bn, GB arrays
    for(int i = 0; i < nrowsa; i++){
        m_an_coefficients[i] = 0.0;
    }
    
    for(int i = 0; i < nrowsb; i++){
        m_bn_coefficients[i] = 0.0;
    }
    
    for(int i = 0; i < nGBParams; i++){
        m_GBParams[i] = 0.0;
    }
    
    // Array for the various alpha/betaLR parameters -- delta parameters?? or just calculate them.
    calculate_aCoefficients(m_an_coefficients, m_Metallicity);  // also needs logMetallicity right? You'd think .. .
    calculate_bCoefficients(m_bn_coefficients, m_Metallicity, m_massCutoffs);  // Same here...
    
    // Zero age main sequence parameters
    m_MZAMS             = MZAMS;
    // m_RZAMS             = RadiusZAMS(m_MZAMS, m_logMetallicityXi);
    m_RZAMS             = RadiusZAMS(m_MZAMS, radius_coefficient_theta, radius_coefficient_iota, radius_coefficient_kappa, radius_coefficient_lamda, radius_coefficient_mu, radius_coefficient_nu, radius_coefficient_xi, radius_coefficient_omicron, radius_coefficient_Pi);
    //m_LZAMS             = LuminosityZAMS(m_MZAMS, m_logMetallicityXi);
    m_LZAMS             = LuminosityZAMS(m_MZAMS, luminosity_coefficient_alpha, luminosity_coefficient_beta, luminosity_coefficient_gamma, luminosity_coefficient_delta, luminosity_coefficient_epsilon, luminosity_coefficient_zeta, luminosity_coefficient_eta);
    m_TZAMS             = calculateTemperature(m_LZAMS, m_RZAMS);
    m_omegaZAMS         = rotationalAngularVelocityDistribution(m_MZAMS, m_RZAMS, options, r); //0.0; 

    // rotationalAngularVelocityDistribution(m_MZAMS, m_RZAMS, options, r);
    // rotationalAngularVelocityDistribution(m_MZAMS, m_RZAMS); 
    //0.0;//angularFrequencyZAMS(m_MZAMS,m_RZAMS);
    
    // Variables for the effective initial Zero Age Main Sequence parameters corresponding to Mass0
    m_RZAMS0            = m_RZAMS;
    m_LZAMS0            = m_LZAMS;
    m_TZAMS0            = m_TZAMS;
    m_omegaZAMS0        = m_omegaZAMS; 

    //rotationalAngularVelocityDistribution(m_MZAMS, m_RZAMS, options, r); 
    //0.0;//angularFrequencyZAMS(m_MZAMS,m_RZAMS);

    // Current timestep variables
    m_Mass              = m_MZAMS;          // Current mass
    m_Mass0             = m_MZAMS;          // Current effective initial mass
    m_Luminosity        = m_LZAMS;          // Current luminosity
    m_Radius            = m_RZAMS;          // Current radius
    m_Temperature       = m_TZAMS;          // Current temperature
    m_coreMass          = 0.0;              // Current Core Mass
    m_HeCoreMass        = 0.0;              // Current He core mass
    m_COCoreMass        = 0.0;              // Current CO core mass
    m_envMass           = initialEnvelopeMass(m_Mass);           // Current envelope mass
    m_Mu                = 0.0;              // Current Mu
    m_coreRadius        = 0.0;              // Current core radius
    
    m_omegaBreak        = omegaBreak(m_Mass, m_Radius); // sqrt(pow(2.0*pi,2.0)*(m_MZAMS)/pow(RsolToAU*m_RZAMS,3.0)); // Can we have the equation for this in a function? // SIMON: Changed to use function
    m_omega             = m_omegaZAMS;
    m_omegaPrev         = m_omegaZAMS;

    
    m_momentInertia     = 0.0; //momentOfInertia3(m_Mass,  m_coreMass,  m_Radius, 0.0,
                                        //   m_stellarType);
    m_angularMomentum   = 0.0; //rotationalAngularMomentum(m_momentInertia, m_omega);

	m_radialExpansionTimescale = 0.0;
	
    // Previous timestep variables
    m_MassPrev              = m_MZAMS;          // Current mass
    m_Mass0Prev             = m_MZAMS;          // Current effective initial mass
    m_LuminosityPrev        = m_LZAMS;          // Current luminosity
    m_RadiusPrev            = m_RZAMS;          // Current radius
    m_TemperaturePrev       = m_TZAMS;          // Current temperature
    m_coreMassPrev          = 0.0;              // Current Core Mass
    m_HeCoreMassPrev        = 0.0;              // Current He core mass
    m_COCoreMassPrev        = 0.0;              // Current CO core mass
    m_envMassPrev           = initialEnvelopeMass(m_Mass);           // Current envelope mass
    m_MuPrev                = 0.0;              // Current Mu
    m_coreRadiusPrev        = 0.0;              // Current core radius
    
    //
    m_omegaTidesIndividualDiff = 0.0;
    m_omegaMagneticBrakingDiff = 0.0;
    m_MassTransferDiff = 0.0;
    m_massLossDiff      = 0.0;
    m_omegaTidesIndividualDiffPrev = 0.0;
    m_omegaMagneticBrakingDiffPrev = 0.0;
    m_MassTransferDiffPrev = 0.0;
    m_massLossDiffPrev      = 0.0;
	
	
	// Binding energies
	m_NanjingLambda			= NEVER_SET;
	m_LoveridgeLambda		= NEVER_SET;
	m_LoveridgeWindsLambda	= NEVER_SET;
	m_KruckowTopLambda		= NEVER_SET;
	m_KruckowMidLambda		= NEVER_SET;
	m_KruckowBotLambda		= NEVER_SET;
	m_DewiLambda			= NEVER_SET;

	// Logarithmic mass-radius derivatives: dlogR/dlogM
	m_zetaThermal	= NEVER_SET;
	m_zetaNuclear	= NEVER_SET;
	m_zetaSoberman	= NEVER_SET;
	m_zetaSobermanHelium = NEVER_SET;
	m_zetaHurley	= NEVER_SET;
	m_zetaHurleyHelium = NEVER_SET;
	m_zetaSimple	= NEVER_SET;
	
	// Stellar timescales
	m_dynamicalTimescale	= NEVER_SET;
	m_thermalTimescale		= NEVER_SET;
	m_nuclearTimescale		= NEVER_SET;
    
    // Initialise pulsar parameters
    m_pulsarMagneticField = NEVER_SET;                    // Pulsar magnetic field strength (G)
    m_pulsarSpinPeriod = NEVER_SET;                       // Pulsar spin period (ms)
    m_pulsarSpinFrequency = NEVER_SET;                    // Pulsar spin frequency (rads per second)
    m_pulsarSpinDownRate = NEVER_SET;                     // Pulsar spin down rate (Pdot, dimensionless)
	
	// Supernova Flags
	m_flagHrichSN		 = false;
	m_flagHpoorSN		 = false;
	flagUSSN			 = false;
    flagSNPrev  		 = false;
    flagExperiencedCCSN  = false;
	flagExperiencedECSN  = false;
    flagExperiencedPISN  = false;                            // Flag to know if at any point in time the star experienced a pair instability supernova
    flagExperiencedPPISN = false;                           // Flag to know if at any point in time the star experienced a pulsational pair instability supernova
	m_runawayFlag		 = false;
    m_fallback           = NEVER_SET;

	
	// Flags
    flagRLOF    	= false;
	experiencedRLOF = false;
    flagInitiateMT 	= false;
    flagSN      	= false;
    flagPISN        = false;
    flagPPISN       = false;
    flagECSN    	= false;
    m_error     	= false;
	m_LBVphaseFlag	= false;
    flagRLOFPrev	= false;
    m_errorPrev 	= false;
    m_fastPhaseCaseA    = false;
    m_firstMassTransferEpisode = false;
    m_initialMassTransferCase = NO_MASS_TRANSFER;
	flagRecycledNS	= false;
	flagRLOFontoaNS	= false;
    flagCHEvolution = false;
    m_Period = NEVER_SET;
    m_preSNeOrbitalEnergy=NEVER_SET;
    m_postSNeOrbitalEnergy=NEVER_SET;
    
    // Calculate the timescales -- might not be necessary here due to timescales function
    m_t_BGB             = lifetimeToBGB(m_an_coefficients, m_MZAMS);            // Lifetime to BGB
    m_t_hook            = lifetimeHook(m_an_coefficients, m_MZAMS, m_t_BGB);    // calculate t_hook
    m_t_MS              = lifetimeMainSequence(m_an_coefficients, m_MZAMS, m_logMetallicityXi, m_t_BGB, m_t_hook);  // MS Lifetime
    
    // Initialise timescales
    for(int i = 0; i < nTimescales; i++){
        m_timescales[i] = 0.0;
    }
    
    // options
    m_remnantMassPrescription = options.remnantMassPrescription; // Since no options have been provided.
    m_SNengine = options.fryerSupernovaEngine;
    
    m_totalMassAtCompactObjectFormation		= NEVER_SET;
    m_HeCoreMassAtCompactObjectFormation 	= NEVER_SET;
    m_COCoreMassAtCompactObjectFormation 	= NEVER_SET;
    m_coreMassAtCompactObjectFormation 		= NEVER_SET;        // What was the core mass of this star when it formed a compact object
    m_HeCoreMassAtCommonEnvelope 			= NEVER_SET;		// He mass at CE
    m_COCoreMassAtCommonEnvelope 			= NEVER_SET;		// CO mass at CE
    m_coreMassAtCommonEnvelope 				= NEVER_SET;        // Core mass at CE
	m_lambdaAtCommonEnvelope				= NEVER_SET;		// Lambda value used at CE
	m_bindingEnergyAtCommonEnvelope			= NEVER_SET;		// Binding energy value used during CE
    m_drawnKickVelocity 					= NEVER_SET;        // Kick velocity the system received during the supernova (km s^-1)
    m_kickVelocity 							= NEVER_SET;        // Kick velocity the system received during the supernova (km s^-1)
    m_supernovaPsi							= NEVER_SET;        // True anomaly (rad)
	
    m_supernovaTheta						= NEVER_SET;        // Angle between the orbital plane and the 'z' axis of supernovae vector (rad)
	m_supernovaPhi							= NEVER_SET;        // Angle between 'x' and 'y', both in the orbital plane of supernovae vector (rad)    

    double kickDirectionPower = options.kickDirectionPower;

    kickDirectionDistribution(options, r, m_supernovaTheta, m_supernovaPhi, kickDirectionPower);

    // Random number U(0,1) for choosing the supernova kick velocity magnitude
    m_supernovaKickVelocityMagnitudeRandomNumber = uniformDistribution(r, 0.0, 1.0);   
    
    // Orbital anomalies
    m_MeanAnomaly       = uniformDistribution(r, 0, 2*M_PI);    // Mean anomaly at instantaneous time of the SN - uniform in [0, 2pi]
    m_EccentricAnomaly  = 0;                                    // Eccentric anomaly at instataneous time of the SN
    m_TrueAnomaly       = 0;                                    // True anomaly at instantaneous time of the SN

    // booleans for fiddling with bits of the code whilst debugging.
//    debugging = false;                                      // Whether to print debugging statements
    debugging = options.debugging;
    isPerturbation = true;                                  // Whether to use perturbations (disable for testing)
    isMassLoss = options.useMassLoss;                       // Whether to use mass loss (disable for testing)
    massLossPrescription = options.massLossPrescription;    // Which mass loss prescription to use
    writeOutput = false;                                    // Whether to write parameters after timestep to standard output
    
}

//// Could also have a constructor for a star with given age, stellar type etc
//Star::Star(int resumeStar, programOptions options, const gsl_rng *r){
//    
//    if(resumeStar == 1){
//        
//        std::cout << "Resuming star 1" << std::endl;
//        
//    }
//    else if(resumeStar == 2){
//        
//        std::cout << "Resuming star 2" << std::endl;
//        
//    }
//    else{
//        
//        std::cout << "Error in selecting which star to resume" << std::endl;
//        
//    }
//    
//    // Initialise age etc to 0
//    m_time                  = 0.0;                  // Physical time the star has been evolved for
//    m_Age                   = 0.0;                  // Effective age of the star, differs for different phases of evolution
//    m_nuclearBurningTime    = 0.0;                  // Nuclear burning time
//    m_dt                    = 0.0;                  // Current timestep
//    m_tau                   = 0.0;                  // Relative time
//    
//    // Assign stellar stype
//    if(MZAMS <= 0.7){
//        m_stellarType = MS_LESS_THAN_07;
//    }
//    else{
//        m_stellarType = MS_MORE_THAN_07;
//    }
//    
//    m_stellarTypePrev = m_stellarType;
//    
//    // Calculate the metallicity dependent numbers
//    m_Metallicity           = Metallicity;
//    m_logMetallicityXi      = log10(m_Metallicity / Zsol);
//    m_logMetallicitySigma   = log10(m_Metallicity);
//    m_logMetallicityRho     = m_logMetallicityXi + 1.0;
//    
//    // Store metallicities in an array as well?
//    m_Mhook   = calcMhook(m_logMetallicityXi);
//    m_MHeF    = calcMHeF(m_logMetallicityXi);
//    m_MFGB    = calcMFGB(m_Metallicity);
//    
//    // Do you really need variables, and then put them into the (makes it a bit clearer what they are I guess?)
//    m_massCutoffs[0] = m_Mhook;
//    m_massCutoffs[1] = m_MHeF;
//    m_massCutoffs[2] = m_MFGB;
//    
//    // Initialise an, bn, GB arrays
//    for(int i = 0; i < nrowsa; i++){
//        m_an_coefficients[i] = 0.0;
//    }
//    
//    for(int i = 0; i < nrowsb; i++){
//        m_bn_coefficients[i] = 0.0;
//    }
//    
//    for(int i = 0; i < nGBParams; i++){
//        m_GBParams[i] = 0.0;
//    }
//    
//    // Array for the various alpha/betaLR parameters -- delta parameters?? or just calculate them.
//    calculate_aCoefficients(m_an_coefficients, m_Metallicity);  // also needs logMetallicity right? You'd think .. .
//    calculate_bCoefficients(m_bn_coefficients, m_Metallicity, m_massCutoffs);  // Same here...
//    
//    // Zero age main sequence parameters
//    m_MZAMS             = MZAMS;
//    m_RZAMS             = RadiusZAMS(m_MZAMS, m_logMetallicityXi);
//    m_LZAMS             = LuminosityZAMS(m_MZAMS, m_logMetallicityXi);
//    m_TZAMS             = calculateTemperature(m_LZAMS, m_RZAMS);
//    m_omegaZAMS         = rotationalAngularVelocityDistribution(m_MZAMS, m_RZAMS, options, r); //0.0;//angularFrequencyZAMS(m_MZAMS,m_RZAMS);
//    
//    // Variables for the effective initial Zero Age Main Sequence parameters corresponding to Mass0
//    m_RZAMS0            = m_RZAMS;
//    m_LZAMS0            = m_LZAMS;
//    m_TZAMS0            = m_TZAMS;
//    m_omegaZAMS0        = m_omegaZAMS; //rotationalAngularVelocityDistribution(m_MZAMS, m_RZAMS, options, r); //0.0;//angularFrequencyZAMS(m_MZAMS,m_RZAMS);
//    
//    m_omegaBreak        = sqrt(pow(2.0*pi,2.0)*(m_MZAMS)/pow(RsolToAU*m_RZAMS,3.0)); // Can we have the equation for this in a function
//    m_omega             = m_omegaZAMS;
//    m_omegaPrev         = m_omegaZAMS;
//    m_angularMomentum   = 0;
//    
//    // Soultion for Tides interaction
//    m_omegaTidesIndividualDiff = 0.0;
//    
//    // Solution for Magnetic Braking interaction
//    m_omegaMagneticBrakingDiff = 0.0;
//    
//    // Solution for Mass Transfer
//    m_MassTransferDiff = 0.0;
//    flagRLOF = false;
//    
//    // Solution for SuperNova
//    flagSN   = false;
//    
//    m_Mass              = m_MZAMS;          // Current mass
//    m_Mass0             = m_MZAMS;          // Current effective initial mass
//    m_MassPrev          = m_MZAMS;          // Mass at the previous timestep
//    m_massLossDiff      = 0.0;          // Amount of mass lost due to winds
//    m_fallback          = 0.0;              // Initialise fallback to 0, record if there is any
//    m_Luminosity        = m_LZAMS;          // Current luminosity
//    m_LuminosityPrev    = m_Luminosity;
//    m_Radius            = m_RZAMS;          // Current radius
//    m_RadiusPrev        = m_RZAMS;          // Radius of previous timestes is same as RZAMS
//    m_Temperature       = m_TZAMS;          // Current temperature
//    m_coreMass          = 0.0;              // Current Core Mass
//    m_HeCoreMassPrev      = 0.0;              // Previous timestep Core Mass
//    m_HeCoreMass        = 0.0;              // Current He core mass
//    m_COCoreMass        = 0.0;              // Current CO core mass
//    m_envMass           = initialEnvelopeMass(m_Mass);           // Current envelope mass
//    m_Mu                = 0.0;              // Current Mu
//    m_coreRadius        = 0.0;              // Current core radius
//    
//    // Calculate the timescales -- might not be necessary here due to timescales function
//    m_t_BGB             = lifetimeToBGB(m_an_coefficients, m_MZAMS);            // Lifetime to BGB
//    m_t_hook            = lifetimeHook(m_an_coefficients, m_MZAMS, m_t_BGB);    // calculate t_hook
//    m_t_MS              = lifetimeMainSequence(m_an_coefficients, m_MZAMS, m_logMetallicityXi, m_t_BGB, m_t_hook);  // MS Lifetime
//    
//    // Initialise timescales
//    for(int i = 0; i < nTimescales; i++){
//        m_timescales[i] = 0.0;
//    }
//    
//    // options
//    m_remnantMassPrescription = options.remnantMassPrescription; // Since no options have been provided.
//    m_SNengine = options.fryerSupernovaEngine;
//    
//    // booleans for fiddling with bits of the code whilst debugging.
//    debugging = false;                                      // Whether to print debugging statements
//    isPerturbation = true;                                  // Whether to use perturbations (disable for testing)
//    isMassLoss = options.useMassLoss;                       // Whether to use mass loss (disable for testing)
//    massLossPrescription = options.massLossPrescription;    // Which mass loss prescription to use
//    writeOutput = false;                                    // Whether to write parameters after timestep to standard output
//    
//}

// Dummy function to learn how to do things.
void Star::Print(){
    
    std::cout << "age = "           << m_Age            << std::endl;
    std::cout << "dt = "            << m_dt             << std::endl;
    std::cout << "MZAMS = "         << m_MZAMS          << std::endl;
    std::cout << "LZAMS = "         << m_LZAMS          << std::endl;
    std::cout << "RZAMS = "         << m_RZAMS          << std::endl;
    std::cout << "TZAMS = "         << m_TZAMS          << std::endl;
    std::cout << "Z = "             << m_Metallicity    << std::endl;
    std::cout << "MHeF = "          << m_MHeF           << std::endl;
    std::cout << "MFGB = "          << m_MFGB           << std::endl;
    std::cout << "Mass0 = "         << m_Mass0          << std::endl;
    std::cout << "Mass = "          << m_Mass           << std::endl;
    std::cout << "Radius = "        << m_Radius         << std::endl;
    std::cout << "Luminosity = "    << m_Luminosity     << std::endl;
    std::cout << "Temperature = "   << m_Temperature    << std::endl;
    std::cout << "stellar type = "  << m_stellarType    << std::endl;
    std::cout << "Mdot = "          << m_Mdot           << std::endl;
    
}

double Star::gettMS(){return m_t_MS;}

double calculateLuminosityCoefficient(int whichCoefficient, double logMetallicityXi){
    /* Calculate an alpha-like metallicity dependent constant
     
     Equation 3 in Tout et al 1996
     
     Parameters
     ------------
     whichCoefficient :
     which coefficient to evaluate
     COEFF_ALPHA     = 0
     COEFF_BETA      = 1
     COEFF_GAMMA     = 2
     COEFF_DELTA     = 3
     COEFF_EPSILON   = 4
     COEFF_ZETA      = 5
     COEFF_ETA       = 6
     
     logMetallicityXi:
     log10(Metallicity / Zsol) where Metallicity = Z (Z = 0.02 = Zsol)
     
     Returns
     --------
     ALPHA :
     constant alpha
     */
    
    double a = Lcoeff[whichCoefficient][0];
    double b = Lcoeff[whichCoefficient][1];
    double c = Lcoeff[whichCoefficient][2];
    double d = Lcoeff[whichCoefficient][3];
    double e = Lcoeff[whichCoefficient][4];
    
    return a + b * logMetallicityXi + c * pow(logMetallicityXi , 2.0) + d * pow(logMetallicityXi, 3.0) + e * pow(logMetallicityXi, 4.0);
    
}

double calculateRadiusCoefficient(int whichCoefficient, double logMetallicityXi){
    /*
     Calculate an alpha-like metallicity dependent constant for the radius
     
     Equation 4 in Tout et al 1996
     
     Parameters
     -----------
     whichCoefficient : int
     Which coefficient to evaluate:
     COEFF_THETA     = 0
     COEFF_IOTA      = 1
     COEFF_KAPPA     = 2
     COEFF_LAMBDA    = 3
     COEFF_MU        = 4
     COEFF_NU        = 5
     COEFF_XI        = 6
     COEFF_OMICRON   = 7
     COEFF_PI        = 8
     
     logMetallicityXi:
     log10(Metallicity / Zsol) where Metallicity = Z (Z = 0.02 = Zsol)
     
     Returns
     ----------
     alpha : float
     Constant alpha

     */
    
    double a = Rcoeff[whichCoefficient][0];
    double b = Rcoeff[whichCoefficient][1];
    double c = Rcoeff[whichCoefficient][2];
    double d = Rcoeff[whichCoefficient][3];
    double e = Rcoeff[whichCoefficient][4];
    
    return a + b * logMetallicityXi + c * pow(logMetallicityXi, 2.0) + d * pow(logMetallicityXi, 3.0) + e * pow(logMetallicityXi, 4.0);
}

void calculateAllLuminosityCoefficients(const double &logMetallicityXi, double &luminosity_coefficient_alpha, double &luminosity_coefficient_beta, double &luminosity_coefficient_gamma, double &luminosity_coefficient_delta, double &luminosity_coefficient_epsilon, double &luminosity_coefficient_zeta, double &luminosity_coefficient_eta){
    /*

    Precalculate all metallicity dependent luminosity coefficients

    See calculateLuminosityCoefficient
    
    */
    luminosity_coefficient_alpha   = calculateLuminosityCoefficient(COEFF_ALPHA,   logMetallicityXi);
    luminosity_coefficient_beta    = calculateLuminosityCoefficient(COEFF_BETA,    logMetallicityXi);
    luminosity_coefficient_gamma   = calculateLuminosityCoefficient(COEFF_GAMMA,   logMetallicityXi);
    luminosity_coefficient_delta   = calculateLuminosityCoefficient(COEFF_DELTA,   logMetallicityXi);
    luminosity_coefficient_epsilon = calculateLuminosityCoefficient(COEFF_EPSILON, logMetallicityXi);
    luminosity_coefficient_zeta    = calculateLuminosityCoefficient(COEFF_ZETA,    logMetallicityXi);
    luminosity_coefficient_eta     = calculateLuminosityCoefficient(COEFF_ETA,     logMetallicityXi);
}

// do these need to be normal functions of class functions?
double LuminosityZAMS(const double &MZAMS, const double &luminosity_coefficient_alpha, const double &luminosity_coefficient_beta, const double &luminosity_coefficient_gamma, const double &luminosity_coefficient_delta, const double &luminosity_coefficient_epsilon, const double &luminosity_coefficient_zeta, const double &luminosity_coefficient_eta){
    /*
     Calculate luminosity at ZAMS in Lsol
     
     Equation 1 in Tout et al 1996
     
     Parameters
     ------------
     MZAMS : float
        Zero age main sequence mass in Msol
     
     Returns
     --------
     LZAMS : float
        Luminosity in Lsol

     */
    
    double top     = luminosity_coefficient_alpha * pow(MZAMS, (5.5)) + luminosity_coefficient_beta * pow(MZAMS, (11.0));
    double bottom  = luminosity_coefficient_gamma + pow(MZAMS, (3.0)) + luminosity_coefficient_delta * pow(MZAMS, (5.0)) + luminosity_coefficient_epsilon * pow(MZAMS, (7.0)) + luminosity_coefficient_zeta * pow(MZAMS,(8.0)) + luminosity_coefficient_eta * pow(MZAMS,(9.5));
    
    return top/bottom;
}

void calculateAllRadiusCoefficients(const double &logMetallicityXi, double &radius_coefficient_theta, double &radius_coefficient_iota, double &radius_coefficient_kappa, double &radius_coefficient_lamda, double &radius_coefficient_mu, double &radius_coefficient_nu, double &radius_coefficient_xi, double &radius_coefficient_omicron, double &radius_coefficient_Pi){
    /*
    
    Calculate metallicity dependent radius coefficients at start up.

    See calculateRadiusCoefficient

    Parameters
    -----------
    logMetallicityXi : double
        log of metallicity

    Returns
    --------

    */
    radius_coefficient_theta   = calculateRadiusCoefficient(COEFF_THETA,   logMetallicityXi);
    radius_coefficient_iota    = calculateRadiusCoefficient(COEFF_IOTA,    logMetallicityXi);
    radius_coefficient_kappa   = calculateRadiusCoefficient(COEFF_KAPPA,   logMetallicityXi);
    radius_coefficient_lamda   = calculateRadiusCoefficient(COEFF_LAMBDA,  logMetallicityXi);
    radius_coefficient_mu      = calculateRadiusCoefficient(COEFF_MU,      logMetallicityXi);
    radius_coefficient_nu      = calculateRadiusCoefficient(COEFF_NU,      logMetallicityXi);
    radius_coefficient_xi      = calculateRadiusCoefficient(COEFF_XI,      logMetallicityXi);
    radius_coefficient_omicron = calculateRadiusCoefficient(COEFF_OMICRON, logMetallicityXi);
    radius_coefficient_Pi      = calculateRadiusCoefficient(COEFF_PI,      logMetallicityXi);
}

double RadiusZAMS(const double &MZAMS, const double &radius_coefficient_theta, const double &radius_coefficient_iota, const double &radius_coefficient_kappa, const double &radius_coefficient_lamda, const double &radius_coefficient_mu, const double &radius_coefficient_nu, const double &radius_coefficient_xi, const double &radius_coefficient_omicron, const double &radius_coefficient_Pi){
    /*
     Calculate radius at ZAMS in units of Rsol
     
     Equation 2 in Tout et al 1996
     
     Parameters
     -----------
     MZAMS : float
        Zero age main sequence mass in Msol
     
     Returns
     --------
     RZAMS : float
     Radius in units of Rsol

     */
    
    double top     = radius_coefficient_theta * pow(MZAMS,(2.5)) + radius_coefficient_iota * pow(MZAMS,(6.5)) + radius_coefficient_kappa * pow(MZAMS,(11.0)) + radius_coefficient_lamda * pow(MZAMS,(19.0)) + radius_coefficient_mu * pow(MZAMS,(19.5));
    double bottom  = radius_coefficient_nu + radius_coefficient_xi * pow(MZAMS,(2.0)) + radius_coefficient_omicron * pow(MZAMS,(8.5)) + pow(MZAMS,(18.5)) + radius_coefficient_Pi * pow(MZAMS,(19.5));
    
    return top/bottom;
}

void calculate_aCoefficients(double an_coefficients[nrowsa], double Metallicity){
    /*
     
     Parameters
     -----------
     an_coefficients : double array
        Array containing the a_n coefficients
     Metallicity : double
        Metallicity Z (Z = 0.02 = Zsol)
     
     Returns
     --------

     */
    
    // Declare variables
    double alpha    = 0;
    double beta     = 0;
    double gamma    = 0;
    double eta      = 0;
    double mu       = 0;
    
    // Ensure that Metallicity is in the correct units
    double logMetallicityXi        = log10(Metallicity / Zsol);         // called xi in Hurley Appendix
    double logMetallicitySigma     = log10(Metallicity);                // called sigma in Hurley Appendix
    //double logMetallicityRho       = logMetallicityXi + 1.0;            // called rho in Hurley Appendix
    
    for(int i = 0; i < nrowsa; i++){
        
        // General equation
        alpha   = aCoefficients[i][0];
        beta    = aCoefficients[i][1];
        gamma   = aCoefficients[i][2];
        eta     = aCoefficients[i][3];
        mu      = aCoefficients[i][4];
        
        an_coefficients[i] = alpha + beta*logMetallicityXi + gamma*pow(logMetallicityXi, 2.0) + eta*pow(logMetallicityXi, 3.0) + mu*pow(logMetallicityXi, 4.0);
        
    }
    
    // Watch out - you were just about to make the mistake where you assume element n is a_n, but it is actually a_n+1
    // Hence why part of this function should be to fix this!!
    // Does this pass by value or reference atm?
    
    // Special Cases -- format is i = n - 1 where n is the a_n number given in the appendix
    double logan = 0;                                            // WOLVERINE!! log_10 of the a_n
    
    an_coefficients[11 - 1] *= an_coefficients[14 - 1];
    an_coefficients[12 - 1] *= an_coefficients[14 - 1];
    
    logan = std::max((0.097 - 0.1072*(logMetallicitySigma + 3.0)), std::max(0.097, std::min(0.1461, (0.1461 + 0.1237*(logMetallicitySigma + 2.0)))));
    an_coefficients[17 - 1] = pow(10.0, logan);
    
    an_coefficients[18 - 1] *= an_coefficients[20 - 1];
    an_coefficients[19 - 1] *= an_coefficients[20 - 1];
    
    an_coefficients[29 - 1] = pow(an_coefficients[29 - 1], (an_coefficients[32 - 1]));
    
    an_coefficients[33 - 1] = std::min(1.4, 1.5135 + 0.3769*logMetallicityXi);
    
    an_coefficients[42 - 1] = std::min(1.25, std::max(1.1, an_coefficients[42 - 1]));
    
    an_coefficients[44 - 1] = std::min(1.3, std::max(0.45, an_coefficients[44 - 1]));
    
    an_coefficients[49 - 1] = std::max(an_coefficients[49 - 1], 0.145);
    
    an_coefficients[50 - 1] = std::min(an_coefficients[50 - 1], (0.306 + 0.053*logMetallicityXi));
    
    an_coefficients[51 - 1] = std::min(an_coefficients[51 - 1], (0.3625 + 0.062*logMetallicityXi));
    
    if(Metallicity > 0.01){
        an_coefficients[52 - 1] = std::min(an_coefficients[52 - 1], 1.0);
        an_coefficients[53 - 1] = std::min(an_coefficients[53 - 1], 1.1);
    }
    else{
        an_coefficients[52 - 1] = std::max(an_coefficients[52 - 1], 0.9);
        an_coefficients[53 - 1] = std::max(an_coefficients[53 - 1], 1.0);
    }
    
    an_coefficients[57 - 1] = std::min(1.4, an_coefficients[57 - 1]);
    an_coefficients[57 - 1] = std::max((0.6355 - 0.4192*logMetallicityXi), std::max(1.25, an_coefficients[57 - 1]));
    
    an_coefficients[62 - 1] = std::max(0.065, an_coefficients[62 - 1]);
    
    if(Metallicity < 0.004){
        an_coefficients[63 - 1] = std::min(0.055, an_coefficients[63 - 1]);
    }
    
    an_coefficients[64 - 1] = std::max(0.091, std::min(0.121, an_coefficients[64 - 1]));
    
    an_coefficients[66 - 1] = std::max(an_coefficients[66 - 1], std::min(1.6, -0.308 - 1.046*logMetallicityXi));
    
    an_coefficients[66 - 1] = std::max(0.8, std::min(0.8 - 2.0*logMetallicityXi, an_coefficients[66 - 1]));
    
    an_coefficients[68 - 1] = std::max(0.9, std::min(an_coefficients[68 - 1], 1.0));
    
    if(an_coefficients[68 - 1] > an_coefficients[66 -1]){
        // Need to connect this to B = aR(M = a66) given by Equation
        // For now calculate this here, can probably precalculate
        double B_alphaR = calculate_B_alphaR(an_coefficients);
        an_coefficients[64 - 1] = B_alphaR;
    }
    
    an_coefficients[68 - 1] = std::min(an_coefficients[68 - 1], an_coefficients[66 - 1]);
    
    if(Metallicity > 0.01){
        an_coefficients[72 - 1] = std::max(an_coefficients[72 - 1], 0.95);
    }
    
    an_coefficients[74 - 1] = std::max(1.4, std::min(an_coefficients[74 - 1], 1.6));
    
    an_coefficients[75 - 1] = std::max(1.0, std::min(an_coefficients[75 - 1], 1.27));
    an_coefficients[75 - 1] = std::max(an_coefficients[75 - 1], 0.6355 - 0.4192*logMetallicityXi);
    
    an_coefficients[76 - 1] = std::max(an_coefficients[76 - 1], -0.1015564 - 0.2161264*logMetallicityXi - 0.05182516*pow(logMetallicityXi, 2.0));
    
    an_coefficients[77 - 1] = std::max((-0.3868776 - 0.5457078*logMetallicityXi - 0.1463472*pow(logMetallicityXi, 2.0)), std::min(0.0, an_coefficients[77 - 1]));
    
    an_coefficients[78 - 1] = std::max(0.0, std::min(an_coefficients[78 - 1], 7.454 + 9.046*logMetallicityXi));
    
    an_coefficients[79 - 1] = std::min(an_coefficients[79 - 1], std::max(2.0, -13.3 - 18.6*logMetallicityXi));
    
    an_coefficients[80 - 1] = std::max(0.0585542, an_coefficients[80 - 1]);
    
    an_coefficients[81 - 1] = std::min(1.5, std::max(0.4, an_coefficients[81 - 1]));
    
}

// Same for bCoefficients
void calculate_bCoefficients(double bn_coefficients[nrowsb], double Metallicity, double massCutoffs[nMassCutoffs]){
    /*
     
     Parameters
     -----------
     bn_coefficients : double array
     Array containing the b_n coefficients - should really have bumped everything up by 1 to avoid this index - 1 problem.
     Metallicity : double
     Metallicity Z (Z = 0.02 = Zsol)
     massCutoffs : double array
     Array containing Mhook, MHeF, MFGB for this star
     
     Returns
     --------
     
     
     */
    
    double alpha    = 0;
    double beta     = 0;
    double gamma    = 0;
    double eta      = 0;
    double mu       = 0;
    
    // Ensure that Metallicity is in the correct units
    double logMetallicityXi        = log10(Metallicity / Zsol);         // called xi in Hurley Appendix
    double logMetallicitySigma     = log10(Metallicity);                // called sigma in Hurley Appendix
    double logMetallicityRho       = logMetallicityXi + 1.0;            // called rho in Hurley Appendix
    
    for(int i = 0; i < nrowsb; i++){
        
        // General equation
        alpha   = bCoefficients[i][0];
        beta    = bCoefficients[i][1];
        gamma   = bCoefficients[i][2];
        eta     = bCoefficients[i][3];
        mu      = bCoefficients[i][4];
        
        bn_coefficients[i] = alpha + beta*logMetallicityXi + gamma*pow(logMetallicityXi, 2.0) + eta*pow(logMetallicityXi, 3.0) + mu*pow(logMetallicityXi, 4.0);
        
    }
    
    // Special cases
    bn_coefficients[1 - 1] = std::min(0.54, bn_coefficients[1 - 1]);
    
    bn_coefficients[2 - 1] = pow(10.0, (-4.6739 - 0.9394 * logMetallicitySigma));
    bn_coefficients[2 - 1] = std::min(std::max(bn_coefficients[2 - 1], (-0.04167 + 55.67*Metallicity)), (0.4771 - 9329.21*pow(Metallicity, 2.94)));
    
    bn_coefficients[3 - 1] = std::max(-0.1451, (-2.2794 - 1.5175*logMetallicitySigma - 0.254*pow(logMetallicitySigma, 2.0)));
    bn_coefficients[3 - 1] = pow(10.0, bn_coefficients[3 - 1]);
    if(Metallicity > 0.004){
        bn_coefficients[3 - 1] = std::max(bn_coefficients[3 - 1], 0.7307 + 14265.1*pow(Metallicity, 3.395));
    }
    
    bn_coefficients[4 - 1] += 0.1231572 * pow(logMetallicityXi, 5.0);
    
    bn_coefficients[6 - 1] += 0.01640687 * pow(logMetallicityXi, 5.0);
    
    bn_coefficients[11 - 1] = pow(bn_coefficients[11 - 1], 2.0);
    
    bn_coefficients[13 - 1] = pow(bn_coefficients[13 - 1], 2.0);
    
    bn_coefficients[14 - 1] = pow(bn_coefficients[14 - 1], bn_coefficients[15 - 1]);
    
    bn_coefficients[16 - 1] = pow(bn_coefficients[16 - 1], bn_coefficients[15 - 1]);
    
    bn_coefficients[17 - 1] = 1.0;
    if(logMetallicityXi > -1.0){
        bn_coefficients[17 - 1] = 1.0 - 0.3880523 * pow((logMetallicityXi + 1.0), 2.862149);
    }
    
    bn_coefficients[24 - 1] = pow(bn_coefficients[24 - 1], bn_coefficients[28 - 1]);
    
    bn_coefficients[26 - 1] = 5.0 - 0.09138012 * pow(Metallicity, -0.3671407);
    
    bn_coefficients[27 - 1] = pow(bn_coefficients[27 - 1], (2.0 * bn_coefficients[28 - 1]));
    
    bn_coefficients[31 - 1] = pow(bn_coefficients[31 - 1], bn_coefficients[33 - 1]);
    
    bn_coefficients[34 - 1] = pow(bn_coefficients[34 - 1], bn_coefficients[33 - 1]);
    
    bn_coefficients[36 - 1] = pow(bn_coefficients[36 - 1], 4.0);
    
    bn_coefficients[37 - 1] = 4.0 * bn_coefficients[37 - 1];
    
    bn_coefficients[38 - 1] = pow(bn_coefficients[38 - 1], 4.0);
    
    bn_coefficients[40 - 1] = std::max(bn_coefficients[40 - 1], 1.0);
    
    bn_coefficients[41 - 1] = pow(bn_coefficients[41 - 1], bn_coefficients[42 - 1]);
    
    bn_coefficients[44 - 1] = pow(bn_coefficients[44 - 1], 5.0);
    
    bn_coefficients[45 - 1] = 1.0 - (2.47162*logMetallicityRho - 5.401682*pow(logMetallicityRho, 2.0) + 3.247361*pow(logMetallicityRho, 3.0));
    if(logMetallicityRho <= 0.0){
        bn_coefficients[45 - 1] = 1.0;
    }
    
    double MHeF = massCutoffs[1];                                      // MHeF
    double MFGB = massCutoffs[2];                                      // MFGB
    bn_coefficients[46 - 1] = -1.0 * bn_coefficients[46 - 1] * log10(MHeF / MFGB);
    
    bn_coefficients[47 - 1] = 1.127733*logMetallicityRho + 0.2344416*pow(logMetallicityRho, 2.0) - 0.3793726*pow(logMetallicityRho, 3.0);
    
    bn_coefficients[51 - 1] -= 0.1343798 * pow(logMetallicityXi, 5.0);
    
    bn_coefficients[53 - 1] += 0.4426929 * pow(logMetallicityXi, 5.0);
    
    bn_coefficients[55 - 1] = std::min((0.99164 - 743.123 * pow(Metallicity, 2.83)), bn_coefficients[55 - 1]);
    
    bn_coefficients[56 - 1] += 0.1140142 * pow(logMetallicityXi, 5.0);
    
    bn_coefficients[57 - 1] -= 0.01308728 * pow(logMetallicityXi, 5.0);
    
    // So you can't do : double myMetallicity = m_Metallicity; unless you pass it to the function
    // How is the best way to use these variables then? pass in all three? make an array with them and pass that around?
    
}



//**********************************************************************************************************//

// HAVE DUPLICATE FUNCTIONS -- NEED TO MAKESURE I ONLY HAVE 1 -- used in calculating a_n / b_n above -- delete (, metallicity) argument


// Radius constants
// Once again, may want to precompute these since they are constant for a given star
double calculate_B_alphaR(double an_coefficients[nrowsa]){
    /*
     Calculate the constant B_alphaR = alphaR(M = a66)
     
     Given in Equation 21a of Hurley et al 2000
     
     Parameters
     -----------
     an_coefficients : double array
     Array containing metallicity dependent a_n coefficients for this star
     
     Returns
     --------
     B_alphaR : double
     Metallicity dependent constant
     
     
     */
    double top      = an_coefficients[58 - 1] * pow(an_coefficients[66 - 1], an_coefficients[60 - 1]);
    double bottom   = an_coefficients[59 - 1] + pow(an_coefficients[66 - 1], an_coefficients[61 - 1]); // Wrong in the arxiv version - says = a59*M**(a61)
    return top / bottom;
}

double calculate_C_alphaR(double an_coefficients[nrowsa]){
    /*
     Calculate the constant C_alphaR = alphaR(M = a67)
     
     Given by Equation 21a in Hurley et al 2000
     
     Parameters
     -----------
     an_coefficients : double array
     Array containing metallicity dependent a_n coefficients for this star
     
     Returns
     --------
     C_alphaR : double
     Metallicity dependent constant
     
     
     */
    double top      = an_coefficients[58 - 1] * pow(an_coefficients[67 - 1], an_coefficients[60 - 1]);
    double bottom   = an_coefficients[59 - 1] + pow(an_coefficients[67 - 1], an_coefficients[61 - 1]); // Wrong in the arxiv version
    return top / bottom;
}


double calculate_alphaR(double an_coefficients[nrowsa], double Mass, unsigned long m_randomSeed){
    /*
     Calculate the radius constant alpha_R
     
     Given by Equations 21, 21a in Hurley et al 2000
     
     Parameters
     -----------
     Mass : double
     Mass in Msol
     
     Returns
     --------
     alphaR : double
     Radius constant alphaR
     
     
     */
    if(Mass < 0.5){
        // a62
        return an_coefficients[62 - 1];
    }
    else if(Mass >= 0.5 and Mass < 0.65){
        // a62 + (a63 - a62)*(Mass - 0.5)/0.15
        return an_coefficients[62 - 1] + (an_coefficients[63 - 1] - an_coefficients[62 - 1])*(Mass - 0.5)/0.15;
    }
    else if(Mass >= 0.65 and Mass < an_coefficients[68 - 1]){
        // a63 + (a64 - a63)*(Mass - 0.65)/(a68 - 0.65)
        return an_coefficients[63 - 1] + (an_coefficients[64 - 1] - an_coefficients[63 - 1])*(Mass - 0.65)/(an_coefficients[68 - 1] - 0.65);
    }
    else if(Mass >= an_coefficients[68 - 1] and Mass < an_coefficients[66 - 1]){
        // a64 + (B_alphaR - a64)*(Mass - a68)/(a66 - a68)
        double B_alphaR = calculate_B_alphaR(an_coefficients);
        return an_coefficients[64 - 1] + (B_alphaR - an_coefficients[64 - 1])*(Mass - an_coefficients[68 - 1])/(an_coefficients[66 - 1] - an_coefficients[68 - 1]);
    }
    else if(Mass >= an_coefficients[66 - 1] and Mass <= an_coefficients[67 - 1]){
        double alphaRtop    = an_coefficients[58 - 1] * pow(Mass, an_coefficients[60 - 1]);    // a58*Mass**(a60)
        double alphaRbottom = an_coefficients[59 - 1] + pow(Mass, an_coefficients[61 - 1]);    // a59 + Mass**(a61)
        return alphaRtop / alphaRbottom;
    }
    else if(Mass > an_coefficients[67 - 1]){
        // C_alphaR + a65*(Mass - a67)
        double C_alphaR = calculate_C_alphaR(an_coefficients);
        return C_alphaR + an_coefficients[65 - 1]*(Mass - an_coefficients[67 - 1]);
    }
    else{
        std::cerr << m_randomSeed << "\tError in calculating alphaR" << std::endl;
        return 0;
    }
}

//*********************************************************************************************************//

double calcMhook(double logMetallicityXi){
    /*
     Calculate the metallicity dependent parameter Mhook above which a hook appears on the MS
     
     Given by Equation 1 in Hurley et al 2000
     
     Parameters
     -----------
     logMetallicityXi : double
     log10(Z / Zsol) (Z = 0.02 = Zsolar)
     
     Returns
     ---------
     Mhook : double
     Mass above which a hook appears on the MS (in Msol)
     
     
     */
    return 1.0185 + 0.16015*logMetallicityXi + 0.0892*logMetallicityXi*logMetallicityXi;
}

double calcMHeF(double logMetallicityXi){
    /*
     Calculate the metallicity dependent maximum initial mass for which He ignites degenerately in the He Flash
     
     Equation 2 in Hurley et al 2000
     
     Parameters
     ------------
     logMetallicityXi : double
     log10(Z / Zsol) (Z = 0.02 = Zsolar)
     
     Returns
     --------
     MHeF : double
     Mass above which a hook appears on the MS
     
     
     */
    return 1.995 + 0.25*logMetallicityXi + 0.087*logMetallicityXi*logMetallicityXi;
}

double calcMFGB(double Metallicity){
    /*
     Calculate the metallicity dependent maximum mass at which He ignites degenerately on the First Giant Branch (FGB)
     
     Equation 3 in Hurley et al 2000
     
     Parameters
     -----------
     Metallicity : double
     Metallicity Z (Z = 0.02 = Zsolar)
     
     Returns
     --------
     MFGB : double
     Mass at which He ignites degenerately on the First Giant Branch (FGB) in Msol
     
     
     */
    double top      = 13.048 * pow((Metallicity / Zsol),(0.06));
    double bottom   = 1.0 + 0.0012 * pow((Zsol / Metallicity), (1.27));
    return top / bottom;
}

double lifetimeToBGB(double an_coefficients[nrowsa], double Mass){
    /*
     Calculate the lifetime to the base of the giant branch (BGB) / end of the Hertzsprung Gap (HG). For high mass stars, t_BGB = t_HeI.
     
     Equation 4 in Hurley et al 2000, plotted in their Figure 5
     
     Parameters
     -----------
     an_coefficients : double array
     Array of length nrowsa containing the metallicity dependent a_n coefficients for this star
     Mass : double
     Stellar mass in Msol
     
     Returns
     --------
     t_BGB : double
     Lifetime to the base of the giant branch BGB in Myr
     
     
     */
    
    double top      = an_coefficients[1 - 1] + an_coefficients[2 - 1] * pow(Mass, 4.0) + an_coefficients[3 - 1] * pow(Mass, 5.5) + pow(Mass, 7.0);
    double bottom   = an_coefficients[4 - 1] * pow(Mass, 2.0) + an_coefficients[5 - 1] * pow(Mass, 7.0);
    
//    std::cout << "lifetimeToBGB: " << Mass << "\t" << top << "\t" << bottom << "\t" << top / bottom << std::endl;
    
    // put tHeI if statement here
    
    return top / bottom;
}

double lifetimeHook(double an_coefficients[nrowsa], double Mass, double t_BGB){
    /*
     Calculate time to hook
     
     Equations 5, 6 and 7 in Hurley et al 2000
     
     Parameters
     -----------
     an_coefficients : double array
     Array of length nrowsa containing the metallicity dependent a_n coefficients for this star
     Mass : double
     Stellar mass in Msol
     
     Returns
     --------
     t_hook : double
     Lifetime to hook in Myr
     
     
     */
    
    // Should I be putting these - 1 things here -- I don't like it -- just have index like before
    double mu = std::max(0.5, (1.0 - 0.01*std::max((an_coefficients[6 - 1]/pow(Mass, an_coefficients[7 - 1])), (an_coefficients[8 - 1] + an_coefficients[9 - 1]/pow(Mass, an_coefficients[10 - 1])))));
    return mu * t_BGB;
}

// Calculate in timescales section of constructor
double lifetimeMainSequence(double an_coefficients[nrowsa], double Mass, double logMetallicityXi, double t_BGB, double t_hook){
    /*
     Equation 5 in Hurley et al 2000
     
     Parameters
     -----------
     an_coefficients : double array
     Array of length nrowsa containing the metallicity dependent a_n coefficients for this star
     Mass : double
     Mass in Msol
     logMetallicityXi : double
     log10(Z / Zsol) (Z = 0.02 = Zsol)
     t_BGB : double
     Lifetime to BGB in Myr
     t_hook : double
     Lifetime to hook in Myr
     
     Returns
     --------
     t_MS : double
     Main sequence lifetime in Myr
     
     
     */
    
    // For mass < Mhook, x > mu (i.e. for stars without a hook)
    double x = std::max(0.95, std::min((0.95 - 0.03*(logMetallicityXi + 0.30103)), 0.99));
    return std::max(t_hook, (x * t_BGB));
}

double luminosityEndMainSequence(double an_coefficients[nrowsa], double Mass, unsigned long m_randomSeed){
    /*
     Calculate luminosity at the end of the Main Sequence (MS), L_TMS
     
     Equation 8 in Hurley et al 2000
     
     Parameters
     ------------
     Mass : double
     Mass in Msol
     
     Returns
     --------
     L_TMS : double
     Luminoisty at the end of the main sequence in Lsol
     
     
     */
    double top      = an_coefficients[11 - 1] * pow(Mass, 3.0) + an_coefficients[12 - 1] * pow(Mass, 4.0) + an_coefficients[13 - 1] * pow(Mass, (an_coefficients[16 - 1] + 1.8));
    double bottom   = an_coefficients[14 - 1] + an_coefficients[15 - 1] * pow(Mass, 5.0) + pow(Mass, an_coefficients[16 - 1]);
    return top / bottom;
}

double radiusEndMainSequence(double an_coefficients[nrowsa], double Mass, double RZAMS, unsigned long m_randomSeed){
    /*
     Calculate the radius at the end of the Main Sequence (MS), R_TMS
     
     Equation 9a, 9b in Hurley et al 2000
     
     Parameters
     -----------
     an_coefficients : double array
     Array of length nrowsa containing the metallicity dependent a_n coefficients for this star
     Mass : double
     Mass in Msol
     RZAMS : double
     Zero Age Main Sequence (ZAMS) radius in Rsol
     
     Returns
     --------
     R_TMS : double
     Radius of a star at the end of the Main Sequence (MS)
     
     
     */
    
    
    // Declare variables
    double Masterisk    = an_coefficients[17 - 1] + 0.1;
    double RTMS         = 0;                                    // Value we want to calculate and return.
    
    if(Mass <= an_coefficients[17 - 1]){
        double top      = an_coefficients[18 - 1] + an_coefficients[19 - 1] * pow(Mass, an_coefficients[21 - 1]);
        double bottom   = an_coefficients[20 - 1] + pow(Mass, an_coefficients[22 - 1]);
        RTMS = top / bottom;
        if(Mass < 0.5){
            return std::max(RTMS, 1.5 * RZAMS);
        }
        else{
            return RTMS;
        }
    }
    else if(Mass >= Masterisk){
        double top      = cCoefficients[1 - 1] * pow(Mass, 3.0) + an_coefficients[23 - 1] * pow(Mass, an_coefficients[26 - 1]) + an_coefficients[24 - 1] * pow(Mass, (an_coefficients[26 - 1] + 1.5));
        double bottom   = an_coefficients[25 - 1] + pow(Mass, 5.0);
        RTMS =  top / bottom;
        return RTMS;
    }
    else{
        // For stars with masses between a17, a17 + 0.1 interpolate between the end points (y = mx + c)
        double y2_top           = cCoefficients[1 - 1] * pow(Masterisk, 3.0) + an_coefficients[23 - 1] * pow(Masterisk, an_coefficients[26 - 1]) + an_coefficients[24 - 1] * pow(Masterisk, (an_coefficients[26 - 1] + 1.5));
        double y2_bottom        = an_coefficients[25 - 1] + pow(Masterisk, 5.0);
        double y2               = y2_top / y2_bottom;                       // RTMS(Masterisk)
        
        double y1_top           = an_coefficients[18 - 1] + an_coefficients[19 - 1] * pow(an_coefficients[17 - 1], an_coefficients[21 - 1]);
        double y1_bottom        = an_coefficients[20 - 1] + pow(an_coefficients[17 - 1], an_coefficients[22 - 1]);
        double y1               = y1_top / y1_bottom;                       // RTMS(a17)
        
        double gradient_top     = y2 - y1;                                      // y2 - y1 = RTMS(Masterisk) - RTMS(a17)
        double gradient_bottom  = 0.1;                                          // x2 - x1 = a17 + 0.1 - a17
        double gradient         = gradient_top / gradient_bottom;
        
        double intercept        = y1 - gradient * an_coefficients[17 - 1];      // c = y - mx
        return gradient*Mass + intercept;
    }
}

double luminosityBaseGiantBranch(double an_coefficients[nrowsa], double Mass, unsigned long m_randomSeed){
    /*
     Calculate Luminoisty at the base of the Giant Branch (GB), LBGB
     
     Equation 10 in Hurley et al 2000
     
     Parameters
     -----------
     an_coefficients : double array
     Array of length nrowsa containing metallicity dependent an coefficients
     Mass : double
     Mass in Msol
     Metallicity : double
     Metallicity Z (Z = 0.02 = Zsol)
     
     Returns
     --------
     L_BGB : double
     Luminosity at the base of the giant branch (BGB)
     
     
     */
    
    double top = an_coefficients[27 - 1] * pow(Mass, an_coefficients[31 - 1]) + an_coefficients[28 - 1] * pow(Mass, cCoefficients[2 - 1]);
    double bottom = an_coefficients[29 - 1] + an_coefficients[30 - 1] * pow(Mass, cCoefficients[3 - 1]) + pow(Mass, an_coefficients[32 - 1]);
    
    return top / bottom;
    
}

// Note that this is a constant for a given star, so if you see that it is being calculated repeatedly, try pre-calculating it
double calculate_B_alphaL(double an_coefficients[nrowsa]){
    /*
     Calculate
     
     Given by
     
     Parameters
     -----------
     an_coefficients : double array
     Array containing metallicity dependent an coefficients for this star
     Metallicity : double
     Metallicity Z (Z = 0.02 = Zsol)
     
     Returns
     --------
     B_alphaL : double
     Constant B(alphaL)
     
     
     */
    
    double M = 2.0;
    
    double top      = an_coefficients[45 - 1] + an_coefficients[46 - 1] * pow(M, an_coefficients[48 - 1]);
    double bottom   = pow(M, 0.4) + an_coefficients[47 - 1] * pow(M, 1.9);
    
    return top / bottom;
    
}

// You might find that functions you give Metallicity to as an argument don't actually need it -- metallicity dependence is in an/MFGB etc.
// Also, if this is likely to be called in a loop, try to precompute it (only depends on initial values of mass/metallicity right? Of course the only problem with this kind of thing is it makes it less flexible if you have to change one of those)
double calculate_alphaL(double an_coefficients[nrowsa], double Mass, bool &m_error, unsigned long m_randomSeed){
    /*
     Calcluate the luminosity evolution constant alpha_L
     
     Given by Equations 19a and 19b in Hurley et al 2000
     
     Parameters
     ------------
     an_coefficients : double array
     Array containing Metallicity dependent (metallicity Z (Z = 0.02 = Zsol) a_n coefficients for this star
     Mass : double
     Mass in Msol
     
     Returns
     --------
     alpha_L : double
     constant
     
     
     */
    
    double B_alphaL = calculate_B_alphaL(an_coefficients);                                    //calculate_B_alphaL
    
    if(Mass < 0.5){
        return an_coefficients[49 - 1];
    }
    else if(Mass >= 0.5 and Mass < 0.7){
        return an_coefficients[49 - 1] + 5.0*(0.3 - an_coefficients[49 - 1])*(Mass - 0.5);
    }
    else if(Mass >= 0.7 and Mass < an_coefficients[52 - 1]){
        return 0.3 + (an_coefficients[50 - 1] - 0.3)*(Mass - 0.7)/(an_coefficients[52 - 1] - 0.7);
    }
    else if(Mass >= an_coefficients[52 - 1] and Mass < an_coefficients[53 - 1]){
        return an_coefficients[50 - 1] + (an_coefficients[51 - 1] - an_coefficients[50 - 1])*(Mass - an_coefficients[52 - 1])/(an_coefficients[53 - 1] - an_coefficients[52 - 1]);
    }
    else if(Mass >= an_coefficients[53 - 1] and Mass < 2.0){
        return an_coefficients[51 - 1] + (B_alphaL - an_coefficients[51 - 1])*(Mass - an_coefficients[53 - 1])/(2.0 - an_coefficients[53 - 1]);
    }
    else if(Mass >= 2.0){
        double top      = an_coefficients[45 - 1] + an_coefficients[46 - 1]*pow(Mass, an_coefficients[48 - 1]);
        double bottom   = pow(Mass, 0.4) + an_coefficients[47 - 1]*pow(Mass, 1.9);
        return top / bottom;
    }
    else{
        std::cerr << m_randomSeed << "\tError calculating alpha_L" << std::endl;
        m_error = true;
        return 0;
    }
    
}

double radiusEndHertzsprungGap(double an_coefficients[nrowsa], double bn_coefficients[nrowsb], double Mass, double coreMass, double massCutoffs[nMassCutoffs], double logMetallicityXi, unsigned long m_randomSeed){
    /*
     Calculate the radius of a star at the end of the Hertzsprung Gap R_EHG
     
     Equations 7-8 in Hurley et al 2000
     
     Parameters
     -----------
     an_coefficients : double array
     Array containing metallicity dependent an coefficients
     bn_coefficients : double array
     Array containing metallicity dependent bn coefficients
     Mass : double
     Mass in Msol
     coreMass : double
     Mass of core in Msol
     massCutoffs : double array
     Array containing metallicity dependent mass cutoffs
     logMetallicityXi : double
        log of metallicity
     
     Returns
     --------
     R_EHG : double
     Radius of a star at the end of the Hertzsprung Gap (HG) in Rsol
     
     
     */
    
    // Declare Variables
    double MFGB         = massCutoffs[2];
    
    if(Mass < MFGB){
        double LBGB = luminosityBaseGiantBranch(an_coefficients, Mass, m_randomSeed);
        return radiusFirstGiantBranch(bn_coefficients, Mass, LBGB, logMetallicityXi, m_randomSeed);
    }
    else{
        return radiusHeliumIgnition(bn_coefficients, Mass, coreMass, massCutoffs, logMetallicityXi, m_randomSeed);
    }
}

double luminosityEndHertzsprungGap(double an_coefficients[nrowsa], double bn_coefficients[nrowsb], double Mass, double massCutoffs[nMassCutoffs], unsigned long m_randomSeed){
    /*
     Calculate the luminosity of a star at the end of the Hertzsprung Gap L_EHG
     
     Given above Equation 8 in Hurley et al 2000
     
     Parameters
     ------------
     an_coefficients : double array
     Array containing metallicity dependent an coefficients
     bn_coefficients : double array
     Array containing metallicity dependent bn coefficients
     Mass : double
     Mass in Msol
     massCutoffs : double array
     Array containing metallicity dependent mass cutoffs
     
     Returns
     --------
     LEHG : double
     Luminosity at the end of the Hertzsprung Gap (HG)
     
     
     */
    
    double MFGB = massCutoffs[2];
    
    if(Mass < MFGB){
        return luminosityBaseGiantBranch(an_coefficients, Mass, m_randomSeed);
    }
    else{
        return luminosityHeliumIgnition(bn_coefficients, Mass, massCutoffs, m_randomSeed);
    }
}

// luminosity and radius for Hertzsprung Gap
double luminosityHertzsprungGap(double an_coefficients[nrowsa], double bn_coefficients[nrowsb], double massCutoffs[nMassCutoffs], double Mass, double tau, unsigned long m_randomSeed){
    /*
     Calculate the luminosity whilst on the Hertzsprung Gap (HG)
     
     Given by Equation 26 in Hurley et al 2000
     
     Parameters
     -----------
     an_coefficients : double array
     Array containing metallicity dependent an coefficients
     bn_coefficients : double array
     Array containing metallicity dependent bn coefficients
     massCutoffs : double array
     Array containing metallicity dependent mass cutoffs
     Mass : double
     Mass in Msol
     tau : double
     Relative time on HG [0,1]
     
     Returns
     --------
     LHG : double
     Luminosity on the HG in Lsol
     
     
     */
    
    double LTMS = luminosityEndMainSequence(an_coefficients, Mass, m_randomSeed);
    double LEHG = luminosityEndHertzsprungGap(an_coefficients, bn_coefficients, Mass, massCutoffs, m_randomSeed);
    
    return LTMS * pow((LEHG/LTMS), tau);
    
}

double radiusHertzsprungGap(double an_coefficients[nrowsa], double bn_coefficients[nrowsb], double massCutoffs[nMassCutoffs], double Mass, double coreMass, double RZAMS, double tau, double logMetallicityXi, unsigned long m_randomSeed){
    /*
     Calculate the radius on the Hertzsprung Gap (HG)
     
     Given by Equation 27 in Hurley et al 2000
     
     Parameters
     -----------
     an_coefficients : double array
     Array containing metallicity dependent an coefficients
     bn_coefficients : double array
     Array containing metallicity dependent bn coefficients
     massCutoffs : double array
     Array containing metallicity dependent mass cutoffs
     Mass : double
     Mass in Msol
     coreMass : double
     Core mass in Msol
     RZAMS : double
     Zero Age Main Sequence (ZAMS) Radius
     tau : double
     Relative time [0,1]
     
     Returns
     --------
     RHG : double
     Radius whilst on the Hertzsprung Gap
     
     
     */
    
    double RTMS = radiusEndMainSequence(an_coefficients, Mass, RZAMS, m_randomSeed);
    double REHG = radiusEndHertzsprungGap(an_coefficients, bn_coefficients, Mass, coreMass, massCutoffs, logMetallicityXi, m_randomSeed);
    
    return RTMS * pow((REHG/RTMS), (tau));
    
}

double calculate_B_betaL(double an_coefficients[nrowsa]){
    /*
     Calculate the constant B = beta_L(M = a57)
     
     Given by Equation 20 in Hurley et al 2000
     
     Parameters
     -----------
     an_coefficients : double array
     Array containing Metallicity dependent a_n coefficients for this star
     
     Returns
     --------
     B_betaL : double
     Constant
     
     
     */
    
    return std::max(0.0, (an_coefficients[54 - 1] - an_coefficients[55 - 1] * pow(an_coefficients[57 - 1], an_coefficients[56 - 1])));
}

double calculate_betaL(double an_coefficients[nrowsa], double Mass){
    /*
     Calculate the constant beta_L for this star
     
     Given by Equation 20 in Hurley et al 2000
     
     Parameter
     ----------
     an_coefficients : double array
     Array containing Metallicity dependent a_n coefficients for this star
     
     Returns
     --------
     beta_l : double
     constant
     
     
     */
    
    // Get value of B_betaL
    double B_betaL = calculate_B_betaL(an_coefficients);
    
    double beta_L = std::max(0.0, (an_coefficients[54 - 1] - an_coefficients[55 - 1] * pow(Mass, an_coefficients[56 - 1])));
    
    if(Mass > an_coefficients[57 - 1] and beta_L > 0.0){
        beta_L = std::max(0.0, (B_betaL - 10.0*(Mass - an_coefficients[57 - 1])*B_betaL));
        return beta_L;
    }
    else{
        return beta_L;
    }
}

// Once again, since this is likely to be calculated in a loop, is there a way to precalculate this, stick it in an array with the other constants and just pass that to the function getting called in the loop.
double calculate_B_DeltaL(double an_coefficients[nrowsa]){
    /*
     Calculate the constant B for the luminosity perturbation DeltaL
     
     Given by Equation 16 in Hurley et al 2000
     
     Parameters
     -----------
     an_coefficients : double array
     Array containing Metallicity dependent a_n coefficients for this star
     
     Returns
     --------
     B_DeltaL : double
     Luminosity perturbation constant B
     
     
     */
    
    double M = an_coefficients[33 - 1];
    return std::min((an_coefficients[34 - 1]/pow(M, an_coefficients[35 - 1])), (an_coefficients[36 - 1]/pow(M, an_coefficients[37 - 1])));
}

double calculate_DeltaL(double an_coefficients[nrowsa], double massCutoffs[nMassCutoffs], double Mass, unsigned long m_randomSeed){
    /*
     Calculate the luminosity perturbation DeltaL
     
     Given by Equation 16 in Hurley et al 2000
     
     Parameters
     -----------
     an_coefficients : double array
     Array containing the metallicity dependent a_n coefficients for this star
     Mass : double
     Mass in Msol
     
     Returns
     ---------
     DeltaL : double
     Luminosity perturbation
     
     
     */
    
    // Get value of B_DeltaL
    double B_DeltaL = calculate_B_DeltaL(an_coefficients);
    
    double Mhook = massCutoffs[0];
    
    if(Mass < Mhook){
        return 0;   // This really is supposed to be zero
    }
    else if(Mass > Mhook and Mass < an_coefficients[33 - 1]){
        double top      = Mass - Mhook;
        double bottom   = an_coefficients[33 - 1] - Mhook;
        double brackets = top / bottom;
        return B_DeltaL * pow(brackets, 0.4);
    }
    else if(Mass >= an_coefficients[33 - 1]){
        return std::min((an_coefficients[34 - 1]/pow(Mass, an_coefficients[35 - 1])), (an_coefficients[36 - 1]/pow(Mass, an_coefficients[37 - 1])));
    }
    else{
        std::cerr << m_randomSeed << "\tError in DeltaL" << std::endl;
        return 0;
    }
}


double calculate_B_betaR(double an_coefficients[nrowsa]){
    /*
     Calculate constant B = betaRPrime(M = 2.0)
     
     Given by Equation 22a in Hurley et al 2000
     
     Parameters
     ------------
     an_coefficients : double array
     Array containing metallicity dependent a_n coefficients for this star
     
     Returns
     --------
     B_betaR : double
     constant
     
     
     */
    double M        = 2.0;
    double top      = an_coefficients[69 - 1] * pow(M, 3.5);
    double bottom   = an_coefficients[70 - 1] + pow(M, an_coefficients[71 - 1]);
    return top / bottom;
}


double calculate_C_betaR(double an_coefficients[nrowsa]){
    /*
     Calculate the constant C = betaRPrime(M = 16.0)
     
     Given by Equation
     
     Parameters
     ------------
     an_coefficients : double array
     Array containing metallicity dependent a_n coefficients for this star
     
     Returns
     --------
     C_betaR : double
     constant
     
     
     */
    double M        = 16.0;
    double top      = an_coefficients[69 - 1] * pow(M, 3.5);
    double bottom   = an_coefficients[70 - 1] + pow(M, an_coefficients[71 - 1]);
    return top / bottom;
}


double calculate_betaR(double an_coefficients[nrowsa], double Mass, unsigned long m_randomSeed){
    /*
     Calculate the radius constant beta_R
     
     Given by Equations 22a, 22b in Hurley et al 2000
     
     Parameters
     -----------
     an_coefficients : double array
     Array containing metallicity dependent a_n coefficients for this star
     Mass : double
     Mass in Msol
     
     Returns
     --------
     beta_R : double
     constant
     
     
     */
    
    // Get values of B_betaR and C_betaR
    double B_betaR      = calculate_B_betaR(an_coefficients);
    double C_betaR      = calculate_C_betaR(an_coefficients);
    double betaRPrime   = 0;
    
    if(Mass <= 1.0){
        betaRPrime = 1.06;
    }
    else if(Mass > 1.0 and Mass < an_coefficients[74 - 1]){
        betaRPrime = 1.06 + (an_coefficients[72 - 1] - 1.06)*(Mass - 1.0)/(an_coefficients[74 - 1] - 1.06);
    }
    else if(Mass > an_coefficients[74 - 1] and Mass < 2.0){
        betaRPrime = an_coefficients[72 - 1] + (B_betaR - an_coefficients[72 - 1])*(Mass - an_coefficients[74 - 1])/(2.0 - an_coefficients[74 - 1]);
    }
    else if(Mass >= 2.0 and Mass <= 16.0){
        double top      = an_coefficients[69 - 1]*pow(Mass, 3.5);
        double bottom   = an_coefficients[70 - 1] + pow(Mass, an_coefficients[71 - 1]);
        betaRPrime      = top / bottom;
    }
    else if(Mass > 16.0){
        betaRPrime = C_betaR + an_coefficients[73 - 1]*(Mass - 16.0);
    }
    else{
        std::cerr << m_randomSeed << "\tError in betaR" << std::endl;
    }
    
    return betaRPrime - 1.0;
}

double calculate_B_DeltaR(double an_coefficients[nrowsa]){
    /*
     Calculate the value of B = DeltaR(M = 2.0)
     
     Given by Equation 17 in Hurley et al 2000
     
     Parameters
     ------------
     an_coefficients : double array
     Array containing metallicity dependent a_n coefficients for this star
     
     Returns
     ---------
     B_DeltaR : double
     constant
     
     
     */
    double M = 2.0;
    
    double top      = an_coefficients[38 - 1] + an_coefficients[39 - 1] * pow(M, 3.5);
    double bottom   = an_coefficients[40 - 1] * pow(M, 3.0) + pow(M, an_coefficients[41 - 1]);
    
    return (top / bottom) - 1.0;
}

double calculate_DeltaR(double an_coefficients[nrowsa], double massCutoffs[nMassCutoffs], double Mass, unsigned long m_randomSeed){
    /*
     Calculate the value of the radius perturbation DeltaR
     
     Given by Equation 17 in Hurley et al 2000
     
     Parameters
     -----------
     an_coefficients : double array
     Array containing metallicity dependent a_n coefficients for this star
     Mass : double
     Mass in Msol
     
     Returns
     --------
     DeltaR : double
     Radius perturbation
     
     
     */
    
    // Get Mhook and B_DeltaR
    double Mhook    = massCutoffs[0];
    double B_DeltaR = calculate_B_DeltaR(an_coefficients);
    
    if(Mass <= Mhook){
        return 0;   // This really is supposed to be 0
    }
    else if(Mass > Mhook and Mass <= an_coefficients[42 - 1]){
        double top      = Mass - Mhook;
        double bottom   = an_coefficients[42 - 1] - Mhook;
        double brackets = top / bottom;
        return an_coefficients[43 - 1] * pow(brackets, 0.5);
    }
    else if(Mass > an_coefficients[42 - 1] and Mass < 2.0){
        double top      = Mass - an_coefficients[42 - 1];
        double bottom   = 2.0 - an_coefficients[42 - 1];
        double brackets = top / bottom;
        return an_coefficients[43 - 1] + (B_DeltaR - an_coefficients[43 - 1])*pow(brackets, an_coefficients[44 - 1]);
    }
    else if(Mass >= 2.0){
        double top      = an_coefficients[38 - 1] + an_coefficients[39 - 1] * pow(Mass, 3.5);
        double bottom   = an_coefficients[40 - 1] * pow(Mass, 3.0) + pow(Mass, an_coefficients[41 -1]);
        return (top / bottom) - 1.0;
    }
    else{
        std::cerr << m_randomSeed << "\tError in calculate_DeltaR" << std::endl;
        return 0;
    }
}

double calculate_B_gamma(double an_coefficients[nrowsa]){
    /*
     Calculate the constant B = gamma(M=1.0)
     
     Given by Equation 23 in Hurley et al 2000
     
     Parameters
     ------------
     an_coefficients : double array
     Array containing metallicity dependent a_n coefficients for this star
     
     Returns
     --------
     B_gamma : double
     constant
     
     
     */
    double M = 1.0;
    return an_coefficients[76 - 1] + an_coefficients[77 - 1]*pow((M - an_coefficients[78 - 1]), an_coefficients[79 - 1]);
}

double calculate_C_gamma(double an_coefficients[nrowsa]){
    /*
     Calculate the constant C
     
     Given after Equation 23 in Hurley et al 2000
     
     Parameters
     ------------
     an_coefficients : double array
     Array containing metallicity dependent a_n coefficients for this star
     
     Returns
     --------
     C_gamma : double
     constant
     
     
     */
    if(an_coefficients[75 - 1] == 1.0){
        return calculate_B_gamma(an_coefficients);
    }
    else{
        return an_coefficients[80 - 1];
    }
}

//double calculate_gamma(double an_coefficients[nrowsa], double Mass){
double calculate_gamma(double an_coefficients[nrowsa], double Mass, bool &m_error, unsigned long m_randomSeed){
    /*
     Calculate the constant gamma
     
     Given by Equation 23 in Hurley et al 2000
     
     Parameters
     -----------
     an_coefficients : double array
     Array containing metallicity dependent a_n coefficients for this star
     Mass: double
     Mass in Msol
     
     ReturnsDeltaR from function
     
     --------
     gamma : double
     constant
     
     
     */
    
    // Debugging flag should really be passed
    bool debugging = false;
//    bool debugging = true;

	// Alejandro - 16/01/17 - Experimental use of static variables in order to print errors only once. Seem to be working, but to be checked thoroughly.
	static bool errorPrinted = false;
	static unsigned long seed = m_randomSeed;
	
	if (m_randomSeed != seed)
	{
		if(debugging){
			std::cout << "random_seed:\t" << m_randomSeed << "\terrorPrinted:\t" << errorPrinted << std::endl;
		}
		seed = m_randomSeed;
		errorPrinted = false;
	}
    
    // Get B and C for gamma
    double B_gamma  = calculate_B_gamma(an_coefficients);
    double C_gamma  = calculate_C_gamma(an_coefficients);
    double gamma    = 0.0;
    
    if(debugging){
        std::cout << "Mass: \t" << Mass << std::endl;
    }
    if(Mass <= 1.0){
        gamma = an_coefficients[76 - 1] + (an_coefficients[77 - 1] * pow((Mass - an_coefficients[78 - 1]), an_coefficients[79 - 1]));
        if(debugging){
            std::cout << "1) gamma: \t" << gamma << std::endl;
        }

    }
    else if(Mass > 1.0 and Mass <= an_coefficients[75 - 1]){
        double top      = Mass - 1.0;
        double bottom   = an_coefficients[75 - 1] - 1.0;
        double brackets = top / bottom;
        gamma = B_gamma + (an_coefficients[80 - 1] - B_gamma) * pow(brackets, an_coefficients[81 - 1]);
        if(debugging){std::cout << "2) gamma: \t" << gamma << std::endl;}

    }
    else if(Mass > an_coefficients[75 - 1] and Mass <= an_coefficients[75] + 0.1){ // The end point here is wrong in the arxiv version
        gamma = C_gamma - 10.0*(Mass - an_coefficients[75 - 1])*C_gamma;
        if(debugging){std::cout << "3) gamma: \t" << gamma << std::endl;}

    }
    else if(Mass > an_coefficients[75] + 0.1){
        gamma = 0.0;            // This really is zero
        if(debugging){std::cout << "4) gamma: \t" << gamma << std::endl;}

    }
    else{
        std::cerr << m_randomSeed << "\tError calculating gamma in SSE. Mass: " << Mass << std::endl;
        if(debugging){std::cout << "5) gamma: error"<< std::endl;}

    }
    
    // Check that gamma >= 0.0
    if(gamma >= 0.0){
        if(debugging){std::cout << "6) gamma: \t" << gamma << std::endl;}
        return gamma;
    }
    else{
        if(debugging){

        	std::cout << "7) gamma: \t" << gamma << std::endl;
			if(errorPrinted == false){
				std::cerr << m_randomSeed << "\tError calculating gamma. Must be >= 0. gamma, Mass: " << gamma << " " << Mass << std::endl;
				errorPrinted = true;
			}
		}

        return 0;   // gamma can go below 0 for low mass stars -- this bounds it to be 0
    }
    
    return gamma;
}

double calculate_eta(double Mass, double Metallicity, unsigned long m_randomSeed){
    /*
     Calculate the constant eta
     
     Given by Equation 18 in the Hurley et al 2000
     
     Parameters
     ------------
     Mass : double
     Mass in Msol
     Metallicity : double
     Metallicity Z (Z = 0.02 = Zsol)
     
     Returns
     --------
     eta : double
     constant
     
     
     */
    if(Metallicity <= 0.0009){
        if(Mass <= 1.0){
            return 10.0;
        }
        else if(Mass > 1.0 and Mass < 1.1){
            // Linear interpolation between end points
            return 100.0*Mass - 90.0;
        }
        else if(Mass >= 1.1){
            return 20.0;
        }
        else{
            std::cerr << m_randomSeed << "\tError in calculating eta. Using default." << std::endl;
            return 10.0;
        }
    }
    else{
        return 10.0;
    }
}

double luminosityMainSequence(double an_coefficients[nrowsa], double massCutoffs[nMassCutoffs], double timescales[nTimescales], double LZAMS, double mass, double logMetallicityXi, double Metallicity, double time, bool &m_error, unsigned long m_randomSeed){
    /*
     Calculates the logarithm of the luminosity (in units of the ZAMS luminosity) as a function of time on the main sequence
     
     Given by Equations 11, 12, 14, 15 in Hurley et al 2000
     
     Parameters
     -----------
     an_coefficients : double array
     Array containing metallicity dependent a_n coefficients for this star
     Time : double
     Time (after ZAMS) in MYRs
     Mass : double
     Mass in Msol
     Metallicity : double
     Metallicity Z (Z = 0.02 = Zsol)
     LZAMS: double
     Zero Age Main Sequence (ZAMS) Radius. Need to recalculate for m0?
     
     Returns
     --------
     log(L(t)/LZAMS) : double
     Log luminosity as a function of time
     
     
     */
    double epsilon = 0.01;               // deltaT = epsilon * t_hook
    
    // Get things that can be precalculated (pass these as arguments to this function to save recalculating them in a for loop)
    // Things which depend on mass should be treated carefully as we will want mass loss -- changing Mass
    double tMS      = timescales[0];
    double t_BGB    = timescales[1];
    
    double LTMS     = luminosityEndMainSequence(an_coefficients, mass, m_randomSeed);
    double alphaL   = calculate_alphaL(an_coefficients, mass, m_error, m_randomSeed);
    double betaL    = calculate_betaL(an_coefficients, mass);
    double DeltaL   = calculate_DeltaL(an_coefficients, massCutoffs, mass, m_randomSeed);
    
    // calculate mu given by equation 7 in Hurley et al 2000
    double mu       = std::max(0.5, (1.0 - 0.01 * std::max((an_coefficients[6 - 1]/pow(mass, an_coefficients[7 - 1])), (an_coefficients[8 - 1] + an_coefficients[9 - 1]/pow(mass, an_coefficients[10 - 1])))));
    
    // Can be calculated from t_BGB
    double t_hook   = mu * t_BGB;
    double eta      = calculate_eta(mass, Metallicity, m_randomSeed);
    
    double bracketsTop      = time - (1.0 - epsilon)*t_hook;
    double bracketsBottom   = epsilon * t_hook;
    double brackets         = bracketsTop / bracketsBottom;
    
    double tau1     = std::min(1.0, (time/t_hook));
    double tau2     = std::max(0.0, std::min(1.0, brackets));
    double tau      = time / tMS;
    
    // Dummy terms
    double one   = alphaL*tau;
    double two   = betaL*pow(tau, eta);
    double three = (log10(LTMS/LZAMS) - alphaL - betaL)*tau*tau;
    double four  = DeltaL * (tau1*tau1 - tau2*tau2);
    
    double logLMSoverLZAMS = one + two + three - four;
    
    return LZAMS * pow(10.0, logLMSoverLZAMS);
    
}

// Radius main sequence ////////////////////////////////////////////////////////////

double radiusMainSequence(double an_coefficients[nrowsa], double massCutoffs[nMassCutoffs], double timescales[nTimescales], double RZAMS, double Mass, double time, bool &m_error, unsigned long m_randomSeed){
    /*
     Calculate the radius whilst on the Main Sequence (MS)
     
     Given by Equation 13 in Hurley et al 2000
     
     Parameters
     -----------
     an_coefficients : double array
     Array containing metallicity dependent an coefficients
     timescales : double array
     Array containing timescales for this star
     massCutoffs : double array
     Array containing metallicity dependent mass cutoffs
     RZAMS: double
     Zero Age Main Sequence (ZAMS) Radius. Need to recalculate for m0?
     
     Returns
     --------
     Radius : double
     Radius on the main sequence
     
     
     */
    
    double epsilon = 0.01;
    
    // Required timescales
    double tMS  = timescales[0];
    double tBGB = timescales[1];
    
    double RTMS   = radiusEndMainSequence(an_coefficients, Mass, RZAMS, m_randomSeed);
    double alphaR = calculate_alphaR(an_coefficients, Mass, m_randomSeed);
    double betaR  = calculate_betaR(an_coefficients, Mass, m_randomSeed);
    double DeltaR = calculate_DeltaR(an_coefficients, massCutoffs, Mass, m_randomSeed); // Radius perturbation
    double gamma  = calculate_gamma(an_coefficients, Mass, m_error, m_randomSeed);
    
    // calculate mu given by equation 7 in Hurley et al 2000
    double mu       = std::max(0.5, (1.0 - 0.01 * std::max((an_coefficients[6 - 1]/pow(Mass, an_coefficients[7 - 1])), (an_coefficients[8 - 1] + an_coefficients[9 - 1]/pow(Mass, an_coefficients[10 - 1])))));
    
    double t_hook = mu * tBGB;
    
    double secondTermTop    = time - (1.0 - epsilon)*t_hook; // Second term in Equation 15 of Hurley et al 2000
    double secondTermBottom = epsilon * t_hook;
    double secondTerm       = secondTermTop / secondTermBottom;
    
    double tau1 = std::min(1.0, time/t_hook);
    double tau2 = std::max(0.0, std::min(1.0, secondTerm));
    double tau  = time / tMS;
    
    // Dummy terms
    double one   = alphaR * tau;
    double two   = betaR * pow(tau, 10.0);
    double three = gamma * pow(tau, 40.0);
    double four  = (log10(RTMS/RZAMS) - alphaR - betaR - gamma) * tau * tau * tau;
    double five  = DeltaR * (pow(tau1, 3.0) - pow(tau2, 3.0));
    
    // Equation 13
    double logRMSoverRZAMS = one + two + three + four - five;
    
    return RZAMS * pow(10.0, logRMSoverRZAMS);
    
}

double calculateBaseGiantBranchC(double an_coefficients[nrowsa], double GBParams[nGBParams], double Mass, double massCutoffs[nMassCutoffs], double logMetallicityXi, unsigned long m_randomSeed){
    /*
     Calculate the constant C for the base of the giant branch
     
     Given just after Equation 44 in Hurley et al 2000
     
     Parameters
     -----------
     an_coefficients : double array
     Array containing metallicity dependent an coefficients
     GBParams : double array
     Array containing GB parameters
     Mass : double
     Mass in Msol
     massCutoffs : double array
     Array containing metallicity dependent mass cutoffs
     logMetallicityXi : double
     logMetallicityXi
     
     Returns
     --------
     C : double
     Constant
     
     
     */
    
    // Get parameters
    double MHeF         = massCutoffs[1];
    double LBGB_MHeF    = luminosityBaseGiantBranch(an_coefficients, MHeF, m_randomSeed);
    double Mc_LBGB_MHeF = luminosityCoreMassRelation(GBParams, LBGB_MHeF, m_randomSeed);
    
    double c1           = 9.20925E-5;
    double c2           = 5.402216;
    
    return pow(Mc_LBGB_MHeF, 4.0) - c1*pow(MHeF, c2);
}

double coreMassBaseGiantBranch(double an_coefficients[nrowsa], double GBParams[nGBParams], double MCBAGB, double Mass, double massCutoffs[nMassCutoffs], double logMetallicityXi, unsigned long m_randomSeed){
    /*
     Calculate the core mass at the base of the giant branch (BGB)
     
     Given by Equation 44 in Hurley et al 2000
     
     For large enough M, we have McBGB ~ 0.098*Mass**(1.35)
     
     Parameters
     -----------
     an_coefficients : double array
     Array containing metallicity dependent an coefficients
     GBParams : double array
     Array containing GB params
     MCBAGB : double
     Core mass at base of asymptotic giant branch (BAGB)
     Mass : double
     Mass in Msol
     massCutoffs : double array
     Array containing metallicity dependent mass cutoffs
     logMetallicityXi : double
     logMetallicityXi
     
     Returns
     --------
     MCBGB : double
     Core Mass at the base of the giant branch
     
     
     */
    double c1       = 9.20925E-5;
    double c2       = 5.402216;
    double C        = calculateBaseGiantBranchC(an_coefficients, GBParams, Mass, massCutoffs, logMetallicityXi, m_randomSeed);
    
    return std::min((0.95 * MCBAGB), pow((C + c1 * pow(Mass, c2)), (0.25)));
    
}

double coreMassBaseAsymptoticGiantBranch(double bn_coefficients[nrowsb], double Mass, unsigned long m_randomSeed){
    /*
     Calculate the core mass at the Base of the Asymptotic Giant Branch (BAGB)
     
     Given by Equation 66 in Hurley et al 2000
     
     Parameters
     -----------
     bn_coefficients : double array
     Array containing metallicity dependent b_n coefficients for this star
     Mass   : double
     Mass in Msol
     
     Returns
     --------
     MCBAGB : double
     Core mass at the Base of the Asymptotic Giant Branch (AGB)
     
     
     */
    return pow((bn_coefficients[36 - 1]*pow(Mass,bn_coefficients[37 - 1]) + bn_coefficients[38 - 1]),(0.25));
}

double coreMassEndHertzsprungGap(double an_coefficients[nrowsa], double bn_coefficients[nrowsb], double GBParams[nGBParams], double Mass, double massCutoffs[nMassCutoffs], double logMetallicityXi, unsigned long m_randomSeed){
    /*
     Calculate the core mass at the end of the Hertzsprung Gap
     
     Given by Equation 28 in Hurley et al 2000
     
     Parameters
     -----------
     an_coefficients : double array
     Array containing metallicity dependent an coefficients
     bn_coefficients : double array
     Array containing metallicity dependent bn coefficients
     GBParams : double array
     Array containing GB params
     Mass : double
     Mass in Msol
     massCutoffs : double array
     Array containing metallicity dependent mass cutoffs
     logMetallicityXi : double
     logMetallicityXi
     
     Returns
     -------
     MCEHG : double
     Core Mass at the end of the Hertzsprung Gap / base of the giant branch
     
     
     */
    
    double MHeF  = massCutoffs[1];
    double MFGB  = massCutoffs[2];
    double MCBGB = GBParams[8];
    //double MCBGB  = coreMassBaseGiantBranch(an_coefficients, bn_coefficients, GBParams, Mass, massCutoffs, logMetallicityXi);
    
    double LBGB   = luminosityBaseGiantBranch(an_coefficients, Mass, m_randomSeed);
    double MCGB   = luminosityCoreMassRelation(GBParams, LBGB, m_randomSeed);
    double MCHeI  = coreMassHeliumIgnition(bn_coefficients, GBParams, Mass, massCutoffs, logMetallicityXi, m_randomSeed);
    
    if(Mass < MHeF){
        return MCGB;
    }
    else if(Mass >= MHeF and Mass < MFGB){
        return MCBGB;
    }
    else if(Mass >= MFGB){
        return MCHeI;
    }
    else{
        std::cerr << m_randomSeed << "\tError in calculating core mass at end of hertzsprung gap." << std::endl;
        return 0;
    }
}

double calculateRhoHG(double Mass){
    /*
     Calculate the parameter rho for the Hertzsprung Gap (HG)
     
     Given by Equation 29 in Hurley et al 2000
     
     Parameters
     ------------
     Mass : double
     Mass in Msol
     
     Returns
     ---------
     Rho : double
     Constant such that MCTMS = Rho * MCEHG
     
     
     */
    double top      = 1.586 + pow(Mass,(5.25));
    double bottom   = 2.434 + 1.02*pow(Mass,(5.25));
    return top/bottom;
}

double coreMassHertzsprungGap(double an_coefficients[nrowsa], double bn_coefficients[nrowsb], double massCutoffs[nMassCutoffs], double GBParams[nGBParams], double tau, double Mass, double MCHG, double Metallicity, double logMetallicityXi, unsigned long m_randomSeed){
    /*
     Calculate the core mass on the Hertzsprung Gap (HG)
     
     Given by Equation 30 in Hurley et al 2000
     
     Parameters
     ------------
     timescales : double array
     Array containing timescales
     an_coefficients : double array
     Array containing metallicity dependent an coefficients
     bn_coefficients : double array
     Array containing metallicity dependent bn coefficients
     massCutoffs : double array
     Array containing metallicity dependent mass cutoffs
     GBParams : double array
     Array containing GB parameters
     tau : double
     Relative time in [0,1]
     Mass : double
     Mass in Msol
     MCHG : double
     Current core mass on HG
     Metallicity : double
     Metallicity Z (Z = 0.02 = Zsol)
     logMetallicityXi : double
     log10(Z)
     
     Returns
     --------
     MCHG : double
     New Core mass on HG
     
     
     */
    
    // Calculate core mass at the end of the Hertzsprung gap
    double MCEHG    = coreMassEndHertzsprungGap(an_coefficients, bn_coefficients, GBParams, Mass, massCutoffs, logMetallicityXi, m_randomSeed);
    
    double rhoHG    = calculateRhoHG(Mass);         // Given by Equation 29
    
    // If the star is losing mass, choose MCHG as the maximum of the core mass at the previous time-step (i.e. should take MCHG as an argument) and the value given by Equation 30
    MCHG = std::max(((1.0 - tau)*rhoHG + tau)*MCEHG, (MCHG));
    
    // Debugging -- looks fine
    //    std::cout << "MCEHG = " << MCEHG << std::endl;
    //    std::cout << "rhoHG = " << rhoHG << std::endl;
    //    std::cout << "MCHG = " << MCHG << std::endl;
    
    return MCHG;
}

// Mass dependent -- does it need calculating each time?
double calculateCoreMassParameter_p(double Mass, double massCutoffs[nMassCutoffs], unsigned long m_randomSeed){
    /*
     Calculate the core mass - luminosity relation parameter p
     
     Given by Equation 31 - 38 in Hurley et al 2000
     
     Parameters
     ------------
     Mass : double
     Mass in Msol
     massCutoffs : double array
     Array containing metallicity dependent mass cutoffs
     
     Returns
     --------
     p : double
     p parameter
     
     
     */
    
    // Get things that have already been calculated
    double MHeF = massCutoffs[1];
    
    if(Mass <= MHeF){
        return 6.0;
    }
    else if(Mass > MHeF and Mass < 2.5){
        // Linear interpolation between end points
        double gradient    = (6.0 - 5.0)/(MHeF - 2.5);  // will be negative
        double intercept   = 5.0 - 2.5*gradient;
        return gradient*Mass + intercept;
    }
    else if(Mass >= 2.5){
        return 5.0;
    }
    else{
        std::cerr << m_randomSeed << "\tError calculating core mass parameter p!" << std::endl;
        return 5.0;     // default.
    }
}

double calculateCoreMassParameter_q(double Mass, double massCutoffs[nMassCutoffs], unsigned long m_randomSeed){
    /*
     Calculate the core mass - luminosity relation parameter q
     
     Given by Equations 31 - 38 in Hurley et al 2000
     
     Parameters
     -----------
     Mass : double
     Mass in Msol
     massCutoffs : double array
     Array containing metallicity dependent mass cutoffs
     
     Returns
     --------
     q  : double
     core mass - luminosity relation parameter q
     
     
     */
    
    // Get things that have already been calculated
    double MHeF = massCutoffs[1];
    
    if(Mass <= MHeF){
        return 3.0;
    }
    else if(Mass > MHeF and Mass < 2.5){
        // Linear interpolation between end points
        double gradient    = (3.0 - 2.0)/(MHeF - 2.5);
        double intercept   = 2.0 - 2.5*gradient;
        return gradient*Mass + intercept;
    }
    else if(Mass >= 2.5){
        return 2.0;
    }
    else{
        std::cerr << m_randomSeed << "\tError calculating core mass parameter q!" << std::endl;
        return 2.0;
    }
}

double calculateCoreMassParameter_B(double Mass){
    /*
     Calculate the core mass - luminosity relation parameter B
     
     Given by Equations 31-38 in Hurley et al 2000
     
     Parameters
     ------------
     Mass : double
     Mass in Msol
     
     Returns
     --------
     B : double
     B parameter
     
     
     */
    // Doesn't depend on metallicity
    return std::max(3E4, (500.0 + 1.75E4*pow(Mass,0.6)));
}

double calculateCoreMassParameter_D(double Mass, double massCutoffs[nMassCutoffs], double logMetallicityXi, unsigned long m_randomSeed){
    /*
     Calculate the core mass - luminosity relation parameter D
     
     Given by Equations 31-38 in Hurley et al 2000
     
     Parameters
     ------------
     Mass : double
     Mass in Msol
     massCutoffs : double array
     Array containing metallicity dependent mass cutoffs
     logMetallicityXi : double
     log10(Z / Zsol)
     
     Returns
     ---------
     D : double
     D parameter
     
     
     */
    // Get things in an easy to work with form.
    double MHeF = massCutoffs[1];
    
    double D0   = 5.37 + 0.135*logMetallicityXi;
    double D1   = 0.975 * D0 - 0.18 * Mass;
    double D2   = 0.5 * D0 - 0.06 * Mass;
    double logD = 0;
    
    if(Mass <= MHeF){
        logD = D0;
    }
    else if(Mass > MHeF and Mass < 2.5){
        // Linear interpolation between end points
        double gradient    = (D0 - D1)/(MHeF - 2.5);
        double intercept   = D0 - MHeF*gradient;
        logD        = gradient*Mass + intercept;
    }
    else if(Mass >= 2.5){
        logD = std::max(std::max(-1.0, D1), D2);
    }
    else{
        std::cerr << m_randomSeed << "\tError calculating core mass parameter D!" << std::endl;
        logD = D0;
    }
    
    return pow(10.0, logD);
}

double lifetimeHeliumIgnition(double bn_coefficients[nrowsb], double GBParams[nGBParams], double timescales[nTimescales], double Mass, double massCutoffs[nMassCutoffs], unsigned long m_randomSeed){
    /*
     Calculate the time until Helium ignition t_HeI
     
     Given by Equation 43 in Hurley et al 2000
     
     Parameters
     ------------
     bn_coefficients : double array
     Array containing metallicity dependent bn coefficients
     GBParams : double array
     Array containing GB parameters
     timescales : double array
     Array containing timescales
     Mass : double
     Mass in Msol
     massCutoffs : double array
     Array containing metallicity dependent mass cutoffs for this star
     
     GBParams
     ---------
     List of GB Parameters:
     0 : A_H
     1 : A_H,He
     2 : A_He
     3 : B
     4 : D
     5 : p
     6 : q
     7 : Mx
     8 : McBGB
     
     timescales
     -----------
     List of timescales (all in units of Myrs)
     0 : tMS  - Main sequence
     1 : tBGB - Base Giant Branch
     2 : tHeI - Helium ignition
     3 : tHe  - Helium burning
     #-- First Giant Branch (FGB)
     4 : Giant tinf1
     5 : Giant tinf2
     6 : Giant t(Mx)
     #-- Early Asymptotic Giant Branch (EAGB)
     7 : FAGB tinf1
     8 : FAGB tinf2
     9 : FAGB t(Mx)
     #-- Thermally Pulsating Asymptotic Giant Branch (TPAGB)
     10 : SAGB tinf1
     11 : SAGB tinf2
     12 : SAGB t(Mx)
     13 : TP (tDU?)
     14 : T(McMax)
     #-- not sure what TP, t(McMax) are
     #-- what about tHeMS etc?
     15 : tHeMS
     
     Returns
     --------
     t_HeI : double
     Time until Helium ignition in MYRs
     
     
     */
    
    // Get things that have already been calculated
    double LHeI = luminosityHeliumIgnition(bn_coefficients, Mass, massCutoffs, m_randomSeed);
    
    double AH = GBParams[0];
    double B  = GBParams[3];
    double D  = GBParams[4];
    double p  = GBParams[5];
    double q  = GBParams[6];
    
    double Lx = calculateLuminosityLx(GBParams);
    
    double tinf1 = timescales[4];
    double tinf2 = timescales[5];
    
    if(LHeI <= Lx){
        return tinf1 - ((1.0)/((p-1.0)*AH*D))*pow((D/LHeI), ((p-1.0)/p));
    }
    else if(LHeI > Lx){
        return tinf2 - ((1.0)/((q-1.0)*AH*B))*pow((B/LHeI), ((q-1.0)/q));
    }
    else{
        std::cerr << m_randomSeed << "\tError in tHeI calculation" << std::endl;
        return 0;
    }
    
}

double calculateHydrogenRateConstant(double Mass){
    /*
     Calculate the Hydrogen rate constant AH'
     
     Given just after Equation 43/ before Equation 44 in Hurley et al 2000
     
     Parameters
     -----------
     Mass : double
     Mass in Msol
     
     Returns
     --------
     AH' : double
     Hydrogen rate constant (Msol Lsol^-1 Myr^-1)
     
     
     */
    double logAH = std::max(-4.8, std::min((-5.7 + 0.8*Mass), (-4.1 + 0.14*Mass)));
    return pow(10.0, logAH);
}

double calculateHeliumRateConstant(){
    /*
     Calculate the Helium rate constant AHe
     
     Given by Equation 68 of Hurley et al 2000
     
     Parameters
     -----------
     
     Returns
     --------
     AHe : double
     Helium rate constant (Msol Lsol^-1 Myr^-1)
     
     
     */
    return 7.66E-5;
}

double calculateMixedHydrogenHeliumRateConstant(){
    /*
     Calculate the effective combined rate constant for both hydrogen and helium
     shell burning
     
     Given by Equation 71 in Hurley et al 2000
     
     Parameters
     -----------
     
     Returns
     --------
     AH,He : double
     Combined Hydrogen/Helium rate constant (Msol Lsol^-1 Myr^-1)
     
     
     */
    return 1.27E-5;
}

double calculateCoreMassParameter_Mx(double GBParams[nGBParams]){
    /*
     Calculate point where the two approximations cross, Mx
     
     Given by Equation 38 in Hurley et al 2000
     
     Parameters
     ----------
     GBParams : double array
     Array containing GB Parameters
     
     GBParams
     ---------
     List of GB Parameters:
     0 : A_H
     1 : A_H,He
     2 : A_He
     3 : B
     4 : D
     5 : p
     6 : q
     7 : Mx
     8 : McBGB
     
     Returns
     --------
     Mx : double
     Mass where approximations cross
     
     
     */
    double B = GBParams[3];
    double D = GBParams[4];
    double p = GBParams[5];
    double q = GBParams[6];
    return pow((B / D), (1.0 / (p - q)));
}

double calculateLuminosityLx(double GBParams[nGBParams]){
    /*
     Calculate the luminosity parameter on the first giant branch (FGB) Lx as a
     function of the core mass.
     
     Given by Equation 37 and 38 in Hurley et al 2000
     
     Parameters
     -----------
     GBParams : double array
     Array containing mass and metalllicity dependent GB parameters
     
     GBParams
     ---------
     List of GB Parameters:
     0 : A_H
     1 : A_H,He
     2 : A_He
     3 : B
     4 : D
     5 : p
     6 : q
     7 : Mx
     8 : McBGB
     
     Returns
     --------
     Lx : double
     Luminosity Lx constant in Lsol
     
     
     */
    
    double B = GBParams[3];
    double D = GBParams[4];
    double p = GBParams[5];
    double q = GBParams[6];
    
    double Mx = calculateCoreMassParameter_Mx(GBParams);
    
    return std::min((B*pow(Mx,q)), (D*pow(Mx, p)));   // both should give same answer by definition at the the point where they cross
}

double coreMassHeliumIgnition(double bn_coefficients[nrowsb], double GBParmas[nGBParams], double Mass, double massCutoffs[nMassCutoffs], double logMetallicityXi, unsigned long m_randomSeed){
    /*
     Calculate the core mass at Helium Ignition (HeI)
     
     Given by core mass - Luminosity relation for low mass stars, Equation 37 in
     Hurley et al 2000. For M >= X use Equation 44, replacing Mc(LBGB(MHeF)) with
     Mc(LHeI(MHeF)) (see end of section 5.3 for this note)
     
     Parameters
     -----------
     bn_coefficients : double array
     Array containing metallicity dependent bn coefficients
     GBParams : double array
     Array containing GB params
     Mass : double
     Mass in Msol
     massCutoffs : double array
     Array containing metallicity dependent mass cutoffs
     logMetallicityXi : double
     logMetallicityXi
     
     Returns
     ---------
     McHeI : double
     Core mass at Helium Ignition (HeI)
     
     
     */
    
    double MHeF = massCutoffs[1];
    
    double LHeI         = luminosityHeliumIgnition(bn_coefficients, Mass, massCutoffs, m_randomSeed);
    double LHeI_MHeF    = luminosityHeliumIgnition(bn_coefficients, MHeF, massCutoffs, m_randomSeed);
    double MC_LHeI_MHeF = luminosityCoreMassRelation(GBParmas, LHeI_MHeF, m_randomSeed);
    double MCBAGB       = coreMassBaseAsymptoticGiantBranch(bn_coefficients, Mass, m_randomSeed);
    
    if(Mass < MHeF){
        return luminosityCoreMassRelation(GBParmas, LHeI, m_randomSeed);
    }
    else if(Mass >= MHeF){
        double c1  = 9.20925E-5;
        double c2  = 5.402216;
        double C   = pow(MC_LHeI_MHeF, 4.0) - c1*pow(MHeF, c2);
        return std::min((0.95*MCBAGB), pow((C + c1*pow(Mass, c2)), (1.0/4.0)));
    }
    else{
        std::cerr << m_randomSeed << "\tError in calculating core mass at helium ignition" << std::endl;
        return 0;
    }
    
}

double coreMassGiantBranch(double an_coefficients[nrowsa], double bn_coefficients[nrowsb], double GBParams[nGBParams], double timescales[nTimescales], double time, double Mass, double massCutoffs[nMassCutoffs], double logMetallicityXi, unsigned long m_randomSeed){
    /*
     Calculate the core mass on the giant branch as a function of time
     
     Equation 39 and 45 in Hurley et al 2000
     
     Parameters
     -----------
     an_coefficients : double array
     Array containing metallicity dependent an coefficients
     bn_coefficients : double array
     Array containing metallicity dependent bn coefficients
     GBParams : double array
     Array containing GB parameters
     time : double
     Time after ZAMS in MYRS (tBGB <= time <= tHeI)
     Mass : double
     Mass in Msol
     logMetallicityXi : double
     log10(Z / Zsol) (Z = 0.02 = Zsol)
     
     GBParams
     ---------
     List of GB Parameters:
     0 : A_H
     1 : A_H,He
     2 : A_He
     3 : B
     4 : D
     5 : p
     6 : q
     7 : Mx
     8 : McBGB
     
     timescales
     -----------
     List of timescales (all in units of Myrs)
     0 : tMS  - Main sequence
     1 : tBGB - Base Giant Branch
     2 : tHeI - Helium ignition
     3 : tHe  - Helium burning
     #-- First Giant Branch (FGB)
     4 : Giant tinf1
     5 : Giant tinf2
     6 : Giant t(Mx)
     #-- Early Asymptotic Giant Branch (EAGB)
     7 : FAGB tinf1
     8 : FAGB tinf2
     9 : FAGB t(Mx)
     #-- Thermally Pulsating Asymptotic Giant Branch (TPAGB)
     10 : SAGB tinf1
     11 : SAGB tinf2
     12 : SAGB t(Mx)
     13 : TP (tDU?)
     14 : T(McMax)
     #-- not sure what TP, t(McMax) are
     #-- what about tHeMS etc?
     15 : tHeMS
     
     Returns
     --------
     MCGB : double
     Core mass on the Giant Branch (GB) in Msol
     
     
     */
    double MHeF  = massCutoffs[1];
    
    double AH    = GBParams[0];
    double B     = GBParams[3];
    double D     = GBParams[4];
    double p     = GBParams[5];
    double q     = GBParams[6];
    double MCBGB = GBParams[8];
    
    double tBGB  = timescales[1];
    double tHeI  = timescales[2];
    
    double tinf1 = timescales[4]; // double tinf1    = tBGB + ((1.0)/((p - 1.0)*AH*D)) * pow((D/LBGB), ((p-1.0)/p));
    double tinf2 = timescales[5]; // double tinf2    = tx + ((1.0)/((q - 1.0)*AH*B)) * pow((B/Lx), ((q-1.0)/q));
    double tx    = timescales[6]; // double tx       = tinf1 - (tinf1 - tBGB) * pow((LBGB/Lx), ((p-1.0)/p));
    
    double MCHeI = coreMassHeliumIgnition(bn_coefficients, GBParams, Mass, massCutoffs, logMetallicityXi, m_randomSeed);
    double tau      = (time - tBGB)/(tHeI - tBGB);   // 0 <= tau <= 1
    
    if(tau < 0 or tau > 1){
        std::cout << "tau is outside [0,1] in calculating the core mass on the giant branch" << std::endl;
    }
    
    // Declare variable for core mass on giant branch
    double MCGB = 0;
    
    if(time <= tx){
        MCGB = pow(((p-1.0)*AH*D*(tinf1 - time)), (1.0/(1.0-p)));
    }
    else{
        MCGB = pow(((q-1.0)*AH*B*(tinf2 - time)), (1.0/(1.0-q)));
    }
    
    if(Mass < MHeF){
        return MCGB;
    }
    else{
        // Equation 45
        return MCBGB + (MCHeI - MCBGB)*tau;
    }
}

double luminosityFirstGiantBranch(double GBParams[nGBParams], double timescales[nTimescales], double time, unsigned long m_randomSeed){
    /*
     Calculate the luminosity on the first giant branch (FGB) as a function of
     the core mass.
     
     Given by Equation 37 in Hurley et al 2000
     
     Parameters
     ----------
     GBParams : double array
     Array containing GB parameters
     timescales : double array
     Array containing timescales
     time : double
     time in MYRs
     
     GBParams
     ---------
     List of GB Parameters:
     0 : A_H
     1 : A_H,He
     2 : A_He
     3 : B
     4 : D
     5 : p
     6 : q
     7 : Mx
     8 : McBGB
     
     timescales
     -----------
     List of timescales (all in units of Myrs)
     0 : tMS  - Main sequence
     1 : tBGB - Base Giant Branch
     2 : tHeI - Helium ignition
     3 : tHe  - Helium burning
     #-- First Giant Branch (FGB)
     4 : Giant tinf1
     5 : Giant tinf2
     6 : Giant t(Mx)
     #-- Early Asymptotic Giant Branch (EAGB)
     7 : FAGB tinf1
     8 : FAGB tinf2
     9 : FAGB t(Mx)
     #-- Thermally Pulsating Asymptotic Giant Branch (TPAGB)
     10 : SAGB tinf1
     11 : SAGB tinf2
     12 : SAGB t(Mx)
     13 : TP (tDU?)
     14 : T(McMax)
     #-- not sure what TP, t(McMax) are
     #-- what about tHeMS etc?
     15 : tHeMS
     
     Returns
     --------
     luminosity : double
     luminosity in Lsol
     
     
     */
    
    // Calculate the core mass according to Equation 39, regardless of whether it is the correct expression to use given the stars mass (see coreMassGiantBranch) -- use it as a dummy variable
    
    double AH = GBParams[0];
    double B  = GBParams[3];
    double D  = GBParams[4];
    double p  = GBParams[5];
    double q  = GBParams[6];
    
    double tBGB  = timescales[1];
    double tHeI  = timescales[2];
    double tinf1 = timescales[4];
    double tinf2 = timescales[5];
    double tx    = timescales[6];
    
    double tau  = (time - tBGB)/(tHeI - tBGB);                  // 0 <= tau <= 1
    
    if(tau < 0 or tau > 1){
        std::cout << "tau is outside [0,1] in calculating luminosity on giant branch" << std::endl;
    }
    
    // declare dummy variable for core mass to calculate luminosity
    double coreMass = 0;
    
    if(time <= tx){
        coreMass = pow(((p-1.0)*AH*D*(tinf1 - time)), (1.0/(1.0-p)));
    }
    else{
        coreMass = pow(((q-1.0)*AH*B*(tinf2 - time)), (1.0/(1.0-q)));
    }
    
    // Use the dummy core mass to calculate the luminosity according to Equation 37
    return std::min((B*pow(coreMass,q)), (D*pow(coreMass, p)));
}

double coreMassLuminosityRelation(double GBParams[nGBParams], double coreMass, unsigned long m_randomSeed){
    /*
     Calculate the luminosity as a function of the core mass, used for AGB stars
     
     Given by Equation 37 in Hurley et al 2000
     
     Parameters
     -----------
     GBParams : double array
     Array containing GB params
     coreMass : double
     Core mass in Msol
     
     GBParams
     ---------
     List of GB Parameters:
     0 : A_H
     1 : A_H,He
     2 : A_He
     3 : B
     4 : D
     5 : p
     6 : q
     7 : Mx
     8 : McBGB
     
     Returns
     --------
     Luminosity : double
     Luminosity in Lsol
     
     
     */
    
    double B = GBParams[3];
    double D = GBParams[4];
    double p = GBParams[5];
    double q = GBParams[6];
    
    return std::min((B * pow(coreMass, q)), (D * pow(coreMass, p)));
    
}

double luminosityCoreMassRelation(double GBParams[nGBParams], double Luminosity, unsigned long m_randomSeed){
    /*
     Calculate the core mass for a given luminosity using the Mc - L relation
     
     Given by Equations 37 and 38 in Hurley et al 2000
     
     Parameters
     -----------
     GBParams : double array
     Array containing GB params
     Luminosity : double
     Luminosity in Lsol
     
     GBParams
     ---------
     List of GB Parameters:
     0 : A_H
     1 : A_H,He
     2 : A_He
     3 : B
     4 : D
     5 : p
     6 : q
     7 : Mx
     8 : McBGB
     
     Returns
     --------
     coreMass : double
     Core mass in Msol
     
     
     */
    
    double Lx = calculateLuminosityLx(GBParams);
    
    double B = GBParams[3];
    double D = GBParams[4];
    double p = GBParams[5];
    double q = GBParams[6];
    
    if(Luminosity > Lx){
        return pow((Luminosity / B), (1.0 / q));
    }
    else{
        return pow((Luminosity / D), (1.0 / p));
    }
    
}

// Metallicity dependent constant -- precompute?
double calculateRadiusGiantBranchConstantX(double logMetallicityXi, unsigned long m_randomSeed){
    /*
     Calculate the parameter x for the Giant Branch, a combination of b5 and b7
     
     Given by Equation 47 in Hurley et al 2000
     
     Parameters
     -----------
     logMetallicityXi : double
     log10(Z / Zsol)
     
     Returns
     --------
     X : double
     Constant
     
     
     */
    return 0.30406 + 0.0805*logMetallicityXi + 0.0897*logMetallicityXi*logMetallicityXi + 0.0878*logMetallicityXi*logMetallicityXi*logMetallicityXi + 0.0222*logMetallicityXi*logMetallicityXi*logMetallicityXi*logMetallicityXi;
}

double calculateRadiusGiantBranchConstantA(double bn_coefficients[nrowsb], double Mass){
    /*
     Calculate the constant A for the giant branch radius
     
     Given by Equation 46 in Hurley et al 2000
     
     Parameters
     -----------
     bn_coefficients : double array
     Array containing metallicity dependent b_n coefficients for this star
     Mass : double
     Mass in Msol
     
     Returns
     --------
     A : double
     Constant A for the giant branch radius
     
     
     */
    return std::min((bn_coefficients[4 - 1] * pow(Mass, (-bn_coefficients[5 - 1]))), (bn_coefficients[6 - 1] * pow(Mass, (-bn_coefficients[7 - 1]))));
}

//double radiusFirstGiantBranch(double bn_coefficients[nrowsb], double Mass, double Luminosity){
//    /*
//     Calculate the radius of a star on the first giant branch (FGB) as a function
//     of the initial mass and the current core mass
//     
//     Given by Equation 46 in Hurley et al 2000
//     
//     Parameters
//     -----------
//     bn_coefficients : double array
//     Array containing metallicity dependent bn coefficients for this star
//     Mass : double
//     Mass in Msol
//     Luminosity : double
//     Luminosity in Lsol
//     
//     Returns
//     --------
//     RGB : double
//     Radius on the Giant Branch (GB) in Rsol
//     
//     
//     */
//    // can I just precompute and pass this?
//    double A = calculateRadiusGiantBranchConstantA(bn_coefficients, Mass);
//    
//    return A * (pow(Luminosity, bn_coefficients[1 - 1]) + bn_coefficients[2 - 1] * pow(Luminosity, bn_coefficients[3 - 1]));
//    
//}

double calculateRadiusGiantBranchConstantAFromX(double bn_coefficients[nrowsb], double x, double Mass){
    /*
     Calculate the giant branch radius
     
     Given just above Equation 47 in Hurley et al 2000
     
     Parameters
     -----------
     bn_coefficients : double array
        Array containing metallicity dependent b_n coefficients for this star
     x : double
        Radius mass Exponent x
     Mass : double
        Mass in Msol
     
     Returns
     --------
     
     
     */
    return std::min((bn_coefficients[4 - 1] * pow(Mass, (-bn_coefficients[5 - 1]))), (bn_coefficients[6 - 1] * pow(Mass, (-bn_coefficients[7 - 1]))));
}

double radiusFirstGiantBranch(double bn_coefficients[nrowsb], double Mass, double Luminosity, double logMetallicityXi, unsigned long m_randomSeed){
    /*
     Calculate the radius of a star on the first giant branch (FGB) as a function
     of the initial mass and the current core mass
     
     Given by Equation 46 in Hurley et al 2000
     
     Parameters
     -----------
     bn_coefficients : double array
     Array containing metallicity dependent bn coefficients for this star
     Mass : double
     Mass in Msol
     Luminosity : double
     Luminosity in Lsol
     
     Returns
     --------
     RGB : double
     Radius on the Giant Branch (GB) in Rsol
     
     
     */
//    double A  = calculateRadiusGiantBranchConstantA(bn_coefficients, Mass);
    double  x           = calculateRadiusGiantBranchConstantX(logMetallicityXi, m_randomSeed);               // I think this can be precomputed and stored somewhere
    double  A           = calculateRadiusGiantBranchConstantAFromX(bn_coefficients, x, Mass);
    double  RGB         = A * (pow(Luminosity, bn_coefficients[1 - 1]) + bn_coefficients[2 - 1] * pow(Luminosity, bn_coefficients[3 - 1]));
    bool    debugging   = false;
    
    // Debugging
    if(debugging){
    std::cout << "Radius Giant Branch" << std::endl;
    std::cout << "RGB : Mass = " << Mass << std::endl;
    std::cout << "RGB : Luminosity = " << Luminosity << std::endl;
    std::cout << "RGB : x = " << x << std::endl;
    std::cout << "RBG : A = " << A << std::endl;
    std::cout << "RGB : RGB = " << RGB << std::endl;
    }
    return RGB;
    
}

// Metallicity dependent constant -- precalculate?
double calculateLuminosityHeliumIgnitionAlpha1(double bn_coefficients[nrowsb], double massCutoffs[nMassCutoffs], unsigned long m_randomSeed){
    /*
     Calculate alpha1 parameter for calculating luminosity at helium ignition
     
     Given just after Equation 49 in Hurley et al 2000
     
     Parameters
     ----------
     bn_coefficients : double array
     Array containing metallicity dependent bn coefficients for this star
     
     Returns
     --------
     alpha1 : double
     alpha1 parmeter
     
     
     */
    double MHeF = massCutoffs[1];
    
    double LHeI_MHeF_top    = bn_coefficients[11 - 1] + bn_coefficients[12 - 1] * pow(MHeF, 3.8);
    double LHeI_MHeF_bottom = bn_coefficients[13 - 1] + (MHeF*MHeF);
    double LHeI_MHeF        = LHeI_MHeF_top / LHeI_MHeF_bottom;
    
    double top      = ((bn_coefficients[9 - 1] * pow(MHeF, bn_coefficients[10 - 1])) - LHeI_MHeF);
    double bottom   = LHeI_MHeF;
    
    return top / bottom;
}



double luminosityHeliumIgnition(double bn_coefficients[nrowsb], double Mass, double massCutoffs[nMassCutoffs], unsigned long m_randomSeed){
    /*
     Calculate the luminosity at Helium Ignition (HeI)
     
     Given by Equation 49 in Hurley et al 2000
     
     Parameters
     ------------
     bn_coefficients : double array
     Array containing metallicity dependent bn coefficients
     Mass : double
     Mass in Msol
     massCutoffs : double array
     Array containing metallicity dependent mass cutoffs for this star
     
     Returns
     --------
     LHeI : double
     Luminosity at Helium Ignition (HeI) in Lsol
     
     
     */
    
    double alpha1   = calculateLuminosityHeliumIgnitionAlpha1(bn_coefficients, massCutoffs, m_randomSeed);
    double MHeF     = massCutoffs[1];
    
    if(Mass < MHeF){
        double top      = bn_coefficients[9 - 1] * pow(Mass, bn_coefficients[10 - 1]);
        double bottom   = 1.0 + alpha1*exp(15.0*(Mass - MHeF));
        return top / bottom;
    }
    else if(Mass >= MHeF){
        double top      = bn_coefficients[11 - 1] + bn_coefficients[12 - 1] * pow(Mass, 3.8);
        double bottom   = bn_coefficients[13 - 1] + pow(Mass, 2.0);
        return top / bottom;
    }
    else{
        std::cerr << m_randomSeed << "\tError calculating luminosity at helium ignition!" << std::endl;
        return 0;
    }
}

double radiusHeliumIgnition(double bn_coefficients[nrowsb], double Mass, double coreMass, double massCutoffs[nMassCutoffs], double logMetallicityXi, unsigned long m_randomSeed){
    /*
     Calculate the radius of the star at Helium Ignition (HeI)
     
     Given by Equation 50 in Hurley et al 2000
     
     Parameters
     -----------
     bn_coefficients : double array
     Array containing metallicity dependent bn coefficients for this star
     Mass : double
     Mass in Msol
     massCutoffs : double array
     Array containing metallicity dependent mass cutoffs
     
     Returns
     --------
     RHeI : double
     Radius at Helium Ignition (HeI) in Rsol
     
     
     */
    
    //
    double MFGB = massCutoffs[2];
    double LHeI = luminosityHeliumIgnition(bn_coefficients, Mass, massCutoffs, m_randomSeed);
    double RmHe = minimumRadiusCoreHeliumBurning(bn_coefficients, Mass, coreMass, massCutoffs, logMetallicityXi, m_randomSeed);
    
    double RAGB_LHeI = radiusAsymptoticGiantBranch(bn_coefficients, Mass, LHeI, massCutoffs, m_randomSeed);
    double RGB_LHeI  = radiusFirstGiantBranch(bn_coefficients, Mass, LHeI, logMetallicityXi, m_randomSeed);
    
    if(Mass <= MFGB){
        return RGB_LHeI;
    }
    else if(Mass >= std::max(MFGB, 12.0)){
        // DEBUGGING Check that for lower masses ~16 RmHe < RAGB_LHeI
        //std::cout << "RHeI = min(RmHe, RAGB_LHeI) = min(" << RmHe << "," << RAGB_LHeI << ")" << std::endl;
        return std::min(RmHe, RAGB_LHeI); // Given by Equation 55
    }
    else if(Mass > MFGB and Mass < 12.0){
        double muTop    = log10(Mass/12.0);
        double muBottom = log10(MFGB / 12.0);
        double mu       = muTop / muBottom;
        double RHeI     = RmHe * pow((RGB_LHeI/RmHe), (mu));
        return RHeI;
    }
    else{
        std::cerr << m_randomSeed << "\tError calculating radius at Helium Ignition (HeI)" << std::endl;
        return 0;
    }
    
}

double minimumLuminosityCoreHeliumBurning(double bn_coefficients[nrowsb], double Mass, double massCutoffs[nMassCutoffs], double logMetallicityXi, unsigned long m_randomSeed){
    /*
     Calculate the minimum luminosity during core helium burning (CHeB) for
     intermediate mass (IM) stars
     
     Given by Equation 51 in Hurley et al 2000
     
     Parameters
     -----------
     bn_coefficients : double array
     Array containing metallicity dependent bn coefficients for this star
     Mass : double
     Mass in Msol
     
     Returns
     ---------
     L_min_He : double
     Minimum luminosity during CHeB in Lsol
     
     
     */
    
    double MFGB = massCutoffs[2];
    
    // calculate the fractions piece by piece
    double cFirstTop    = bn_coefficients[17 - 1];
    double cFirstBottom = pow(MFGB, 0.1);
    
    double cSecondTop       = (bn_coefficients[16 - 1] * bn_coefficients[17 - 1]) - bn_coefficients[14 - 1];
    double cSecondBottom    = pow(MFGB, (bn_coefficients[15 - 1] + 0.1));
    
    double c                = (cFirstTop/cFirstBottom) + (cSecondTop/cSecondBottom);
    
    double LHeI             = luminosityHeliumIgnition(bn_coefficients, Mass, massCutoffs, m_randomSeed);
    
    double top              = bn_coefficients[14 - 1] + c*pow(Mass, (bn_coefficients[15 - 1] + 0.1));
    double bottom           = bn_coefficients[16 - 1] + pow(Mass, bn_coefficients[15 - 1]);
    
    return LHeI*(top/bottom);
}

double luminosityZeroAgeHorizontalBranch(double bn_coefficients[nrowsb], double Mass, double coreMass, double massCutoffs[nMassCutoffs], double logMetallicityXi, unsigned long m_randomSeed){
    /*
     Calculate the luminosity for Low Mass (LM) stars at Zero Age Horizontal
     Branch (ZAHB)
     
     Given by Equations 52 and 53 in Hurley et al 2000
     
     Parameters
     -----------
     bn_coefficients : double array
     Array containing the metallicity dependent bn coefficients
     Mass : double
     Mass in Msol
     coreMass : double
     Core mass in Msol
     
     Returns
     --------
     LZAHB : double
     Luminosity at ZAHB in Lsol
     
     
     */
    
    double MHeF     = massCutoffs[1];
    
    double muTop    = Mass - coreMass;
    double muBottom = MHeF - coreMass;
    double mu       = muTop / muBottom;     // 0 <= mu <= 1
    
    double LZHe     = luminosityHeliumMainSequenceZAMS(coreMass);
    double LminHe   = minimumLuminosityCoreHeliumBurning(bn_coefficients, MHeF, massCutoffs, logMetallicityXi, m_randomSeed);
    
    double alpha2Top    = bn_coefficients[18 - 1] + LZHe - LminHe;
    double alpha2Bottom = LminHe - LZHe;
    double alpha2       = alpha2Top / alpha2Bottom;
    
    double firstTop     = 1.0 + bn_coefficients[20 - 1];
    double firstBottom  = 1.0 + (bn_coefficients[20 - 1] * pow(mu, 1.6479));
    
    double secondTop    = bn_coefficients[18 - 1] * pow(mu, bn_coefficients[19 - 1]);
    double secondBottom = 1.0 + alpha2 * exp(15.0*(Mass - MHeF));
    
    return LZHe + (firstTop/firstBottom)*(secondTop/secondBottom);
    
}

double radiusZeroAgeHorizontalBranch(double bn_coefficients[nrowsb], double Mass, double coreMass, double massCutoffs[nMassCutoffs], double logMetallicityXi, unsigned long m_randomSeed){
    /*
     Calculate the radius on the Zero Age Horizontal Branch (ZAHB)
     
     Given by Equation 54 in Hurley et al 2000
     
     Parameters
     -----------
     bn_coefficients : double array
     Array containing metallicity dependent bn coefficients
     Mass : double
     Mass in Msol
     
     Returns
     ---------
     RZAHB : double
     Radius at ZAHB in Rsol
     
     
     */
    
    double MHeF = massCutoffs[1];
    
    double muTop    = Mass - coreMass;
    double muBottom = MHeF - coreMass;
    double mu       = muTop / muBottom;
    
    double fTop     = (1.0 + bn_coefficients[21 - 1]) * pow(mu, bn_coefficients[22 - 1]);
    double fBottom  = 1.0 + bn_coefficients[21 - 1] * pow(mu, bn_coefficients[23 - 1]);
    double constF   = fTop / fBottom;
    
    double RZHe     = radiusHeliumMainSequenceZAMS(coreMass);
    double LZAHB    = luminosityZeroAgeHorizontalBranch(bn_coefficients, Mass, coreMass, massCutoffs, logMetallicityXi, m_randomSeed);
    double RGB      = radiusFirstGiantBranch(bn_coefficients, Mass, LZAHB, logMetallicityXi, m_randomSeed);
    
    return (1.0 - constF)*RZHe + constF*RGB;
}

// Move to be with other minimum
double minimumRadiusCoreHeliumBurning(double bn_coefficients[nrowsb], double Mass, double coreMass, double massCutoffs[nMassCutoffs], double logMetallicityXi, unsigned long m_randomSeed){
    /*
     Calculate minimum radius during the blue loop
     
     Given by Equation 55 in Hurley et al 2000
     
     Parameters
     -----------
     Mass : double
     Mass in Msol
     coreMass : double
     Core mass in Msol
     
     Returns
     -------
     R_min_He : double
     Minimum radius on the blue loop in Rsol
     
     
     */
    double MHeF = massCutoffs[1];
    
    double LZAHB_MHeF   = luminosityZeroAgeHorizontalBranch(bn_coefficients, MHeF, coreMass, massCutoffs, logMetallicityXi, m_randomSeed);
    double LZAHB        = luminosityZeroAgeHorizontalBranch(bn_coefficients, Mass, coreMass, massCutoffs, logMetallicityXi, m_randomSeed);
    double mu           = Mass / MHeF;
    
    if(Mass >= MHeF){
        double top      = (bn_coefficients[24 - 1] * Mass) + pow((bn_coefficients[25 - 1] * Mass), bn_coefficients[26 - 1]) * pow(Mass, bn_coefficients[28 - 1]);
        double bottom   = bn_coefficients[27 - 1] + pow(Mass, bn_coefficients[28 - 1]);
        
        return top / bottom;
        
    }
    else if(Mass < MHeF){
        double fractionTopTop       = bn_coefficients[24 - 1]*MHeF + pow((bn_coefficients[25 - 1] * Mass), bn_coefficients[26 - 1]) * pow(MHeF, bn_coefficients[28 - 1]);
        double fractionTopBottom    = bn_coefficients[27 - 1] + pow(MHeF, bn_coefficients[28 - 1]);
        double fractionTop          = fractionTopTop / fractionTopBottom;
        // still unsure about this line.
        double fractionBottom       = radiusFirstGiantBranch(bn_coefficients, MHeF, LZAHB_MHeF, logMetallicityXi, m_randomSeed); // Mass or MHeF
        return radiusFirstGiantBranch(bn_coefficients, Mass, LZAHB, logMetallicityXi, m_randomSeed) * pow((fractionTop/fractionBottom), mu); //radiusFirstGiantBranch(initialMass, LZAHB, Metallicity)*((fractionTop/fractionBottom)**(mu));
    }
    else{
        std::cerr << m_randomSeed << "\tError calculating minimum radius during core helium burning" << std::endl;
        return 0;
    }
}

// Metallicity dependent, precompute?
double calculateAlphaBAGBAlpha3(double bn_coefficients[nrowsb], double massCutoffs[nMassCutoffs]){
    /*
     Calculate the constant alpha 3 for the BAGB
     
     Given just after Equation 56 in Hurley et al 2000
     
     Parameters
     -----------
     bn_coefficients : double array
     Array containing metallicity dependent bn coefficients for
     massCutoffs : double
     Metallicity dependent mass cutoffs
     
     Returns
     --------
     alpha3 : double
     Metallicity dependent constant

     */
    double MHeF         = massCutoffs[1];
    
    double LBAGB_top    = bn_coefficients[31 - 1] + bn_coefficients[32 - 1] * pow(MHeF, (bn_coefficients[33 - 1] + 1.8));
    double LBAGB_bottom = bn_coefficients[34 - 1] + pow(MHeF, bn_coefficients[33 - 1]);
    double LBAGB        = LBAGB_top / LBAGB_bottom;
    
    double top          = bn_coefficients[29 - 1] * pow(MHeF, bn_coefficients[30 - 1]) - LBAGB;
    double bottom       = LBAGB;
    
    return top / bottom;
    
}

double luminosityBaseAsymptoticGiantBranch(double bn_coefficients[nrowsb], double Mass, double massCutoffs[nMassCutoffs], unsigned long m_randomSeed){
    /*
     Calculate luminosity at the base of the Asymptotic Giant Branch (AGB)
     
     Given by Equation 56 in Hurley et al 2000
     
     Parameters
     -----------
     Mass : double
     Mass in Msol
     
     Returns
     --------
     LBAGB : double
     Luminosity at BAGB in Lsol
     
     
     */
    // Get things that have already been computed
    double MHeF     = massCutoffs[1];
    double alpha3   = calculateAlphaBAGBAlpha3(bn_coefficients, massCutoffs);
    
    if(Mass < MHeF){
        double top      = bn_coefficients[29 - 1] * pow(Mass, bn_coefficients[30 - 1]);
        double bottom   = 1.0 + alpha3 * exp(15.0*(Mass - MHeF));
        return top / bottom;
    }
    else if(Mass >= MHeF){
        double top      = bn_coefficients[31 - 1] + bn_coefficients[32 - 1] * pow(Mass, (bn_coefficients[33 - 1] + 1.8));
        double bottom   = bn_coefficients[34 - 1] + pow(Mass, bn_coefficients[33 - 1]);
        return top / bottom;
    }
    else{
        std::cerr << m_randomSeed << "\tError calculating luminosity at the base of the asymptotic giant branch" << std::endl;
        return 0;
    }
}

double radiusBaseAsymptoticGiantBranch(double bn_coefficients[nrowsb], double Mass, double LBAGB, double massCutoffs[nMassCutoffs], unsigned long m_randomSeed){
    /*
     Calculate the radius at the base of the asymptotic giant branch (BAGB)
     
     As it says just after Equation 56, the radius at the BAGB is given by
     RAGB(M, LBAGB) where RAGB is given by Equation 74 in Hurley et al 2000
     
     Parameters
     -----------
     bn_coefficients : double array
     Array containing metallicity dependent bn coefficents for this star
     Mass : double
     Mass in Msol
     LBAGB : double
     Luminosity at the Base of the Asymptotic Giant Branch (BAGB)
     massCutoffs : double array
     Array containing metallicity dependent mass cutoffs
     
     Returns
     --------
     RBAGB : double
     Radius at the Base of the Asymptotic Giant Branch (BAGB)
     
     
     */
    return radiusAsymptoticGiantBranch(bn_coefficients, Mass, LBAGB, massCutoffs, m_randomSeed);
}

double calculateRadiusAsymptoticGiantBranchConstantA(double bn_coefficients[nrowsb], double Mass, double massCutoffs[nMassCutoffs], unsigned long m_randomSeed){
    /*
     Calculate the constant A for the Asymptotic Giant Branch (AGB) radius
     
     Given just after Equation 74 in Hurley et al 2000
     
     Parameters
     -----------
     bn_coefficients : double array
     Array containing metallicity dependent bn coefficients
     Mass : double
     Mass in Msol
     massCutoffs : double array
     Array containing metallicity dependent mass cutoffs
     
     Returns
     --------
     A : double
     Constant A for the Asymptotic Giant Branch radius
     
     
     */
    double MHeF = massCutoffs[1];
    
    if(Mass >= MHeF){
        return std::min((bn_coefficients[51 - 1] * pow(Mass, -bn_coefficients[52 - 1])), (bn_coefficients[53 - 1] * pow(Mass, -bn_coefficients[54 - 1])));
    }
    else if(Mass <= MHeF - 0.2){
        return bn_coefficients[56 - 1] + bn_coefficients[57 - 1] * Mass;
    }
    else if(Mass > MHeF - 0.2 and Mass < MHeF){
        // Linear interpolation between end points
        double x1           = MHeF - 0.2;
        double y1           = bn_coefficients[56 - 1] + (bn_coefficients[57 - 1] * x1);
        double x2           = MHeF;
        double y2           = bn_coefficients[51 - 1] * pow(x2, -bn_coefficients[52 - 1]);
        double gradient     = (y2 - y1)/(x2 - x1);
        double intercept    = y2 - (gradient * x2);
        return (gradient * Mass) + intercept;
    }
    else{
        std::cerr << m_randomSeed << "\tError calculating constant A for the AGB radius" << std::endl;
        return 0;
    }
}

double radiusAsymptoticGiantBranch(double bn_coefficients[nrowsb], double Mass, double Luminosity, double massCutoffs[nMassCutoffs], unsigned long m_randomSeed){
    /*
     Calculate the radius of a star on the asymptotic giant branch (AGB) as a function
     of the initial mass, luminosity and metallicity
     
     Given by Equation 74 in Hurley et al 2000
     
     Parameters
     -----------
     bn_coefficients : double array
     Array containing the metallicity dependent bn coefficients for this star
     Mass : double
     Mass in Msol
     Luminosity : double
     Luminosity in Lsol
     massCutoffs : double array
     Array containing metallicity dependent mass cutoffs
     
     Returns
     -------
     RAGB : double
     Radius on the Asymptotic Giant Branch (AGB) in Rsol
     
     
     */
    
    double A = calculateRadiusAsymptoticGiantBranchConstantA(bn_coefficients, Mass, massCutoffs, m_randomSeed);
    
    double MHeF = massCutoffs[1];
    double b50  = 0;
    
    if(Mass >= MHeF){
        b50 = bn_coefficients[55 - 1] * bn_coefficients[3 - 1];
    }
    else if(Mass <= MHeF - 0.2){
        b50 = bn_coefficients[3 - 1];
    }
    else if(Mass > MHeF - 0.2 and Mass < MHeF){
        // Linear interpolation between end points given just after Equation 74
        double x1           = MHeF - 0.2;
        double y1           = bn_coefficients[3 - 1];
        double x2           = MHeF;
        double y2           = bn_coefficients[55 - 1] * bn_coefficients[3 - 1];
        double gradient     = (y2 - y1)/(x2 - x1);
        double intercept    = y2 - (gradient * x2);
        b50 = (gradient * Mass) + intercept;
    }
    else{
        std::cerr << m_randomSeed << "\tError calculating the radius during the asymptotic giant branch" << std::endl;
        b50 = bn_coefficients[55 - 1] * bn_coefficients[3 - 1];
    }
    
    double RAGB = A * (pow(Luminosity, bn_coefficients[1 - 1]) + (bn_coefficients[2 - 1] * pow(Luminosity, b50)));
    
    // DEBUGGING
    //std::cout << "A, b1, b2, b50 = " << A << " " << bn_coefficients[1-1] << " " << bn_coefficients[2-1] << " " << b50 << std::endl;
    //std::cout << "Luminosity = " << log10(Luminosity) << " " << Luminosity << std::endl;
    //std::cout << "RAGB = " << RAGB << std::endl;
    
    return RAGB;
    
}

//double calculateAlpha4(const double (&an_coefficients)[nrowsa], const double (&bn_coefficients)[nrowsb], const double (&massCutoffs)[nMassCutoffs]){
//    /*
//     Calculate the constant alpha4
//     
//     Given just after Equation 57 in Hurley et al 2000
//     
//     Parameter
//     -----------
//     an_coefficients : double array
//     Array containing metallicity dependent an coefficents
//     bn_coefficients : double array
//     Array containing metallicity dependent bn coefficents
//     massCutoffs : double array
//     Array containing metallicity dependent mass cutoffss
//     
//     Returns
//     --------
//     alpha4 : double
//     Metallicity dependent constant
//     
//     
//     */
//    
//    // Get/pass things
//    double MHeF = massCutoffs[1];
//    double tBGB_MHeF = lifetimeToBGB(an_coefficients, MHeF); // calculate tBGB for mass M = MHeF
//    
//    double top      = bn_coefficients[41 - 1] * pow(MHeF, bn_coefficients[42 - 1]) + bn_coefficients[43 - 1] * pow(MHeF, 5.0);
//    double bottom   = bn_coefficients[44 - 1] + pow(MHeF, 5.0);
//    
//    return tBGB_MHeF * (top / bottom);
//    
//}

double calculateAlpha4(double an_coefficients[nrowsa], double bn_coefficients[nrowsb], double massCutoffs[nMassCutoffs]){
    /*
     Calculate the constant alpha4

     Given just after Equation 57 in Hurley et al 2000

     Parameter
     -----------
     an_coefficients : double array
     Array containing metallicity dependent an coefficents
     bn_coefficients : double array
     Array containing metallicity dependent bn coefficents
     massCutoffs : double array
     Array containing metallicity dependent mass cutoffss

     Returns
     --------
     alpha4 : double
     Metallicity dependent constant

     
     */

    // Get/pass things
    double MHeF = massCutoffs[1];
    double tBGB_MHeF = lifetimeToBGB(an_coefficients, MHeF); // calculate tBGB for mass M = MHeF

    double top      = bn_coefficients[41 - 1] * pow(MHeF, bn_coefficients[42 - 1]) + bn_coefficients[43 - 1] * pow(MHeF, 5.0);
    double bottom   = bn_coefficients[44 - 1] + pow(MHeF, 5.0);

    return tBGB_MHeF * (top / bottom);

}


double lifetimeCoreHeliumBurning(double an_coefficients[nrowsa], double bn_coefficients[nrowsb], double Mass, double coreMass, double massCutoffs[nMassCutoffs], unsigned long m_randomSeed){
    /*
     Calculate the lifetime of the star during core helium burning
     
     Given by Equation 57 in Hurley et al 2000
     
     Parameters
     -----------
     bn_coefficients : double array
     Array containing metallicity dependent bn coefficients
     Mass : double
     Mass in Msol
     coreMass : double
     Core mass in Msol
     
     Returns
     --------
     t_He : double
     Lifetime during core helium burning in MYRs
     
     
     */
    
    // Get precomputed things
    double MHeF     = massCutoffs[1];
    double alpha4   = calculateAlpha4(an_coefficients, bn_coefficients, massCutoffs);
    
    // FINISH ME rewrite to use timescales, although will be called in timescales as well.
    double tHeMS    = lifetimeHeliumMainSequence(coreMass);  // timescales[15?]
    double tBGB     = lifetimeToBGB(an_coefficients, Mass); // may already be computed, but mass changing -- could use timescales? not used?
    
    double mu       = Mass / MHeF;
    
    if(Mass <  MHeF){
        double t_He = (bn_coefficients[39 - 1] + (tHeMS - bn_coefficients[39 - 1]) * pow((1.0 - mu), (bn_coefficients[40 - 1]))) * (1.0 + alpha4 * exp(15.0*(Mass - MHeF)));
        //std::cout << "t_He = " << t_He << std::endl; // DEBUGGING
        return t_He;
    }
    else if(Mass >= MHeF){
        double top      = bn_coefficients[41 - 1] * pow(Mass, bn_coefficients[42 - 1]) + bn_coefficients[43 - 1] * pow(Mass, 5.0);
        double bottom   = bn_coefficients[44 - 1] + pow(Mass, 5.0);
        double t_He     = tBGB * (top / bottom);
        //std::cout << "t_He = " << t_He << std::endl; // DEBUGGING
        return t_He;
    }
    else{
        std::cerr << m_randomSeed << "\tError calculating lifetime during core helium burning (CHeB)" << std::endl;
        return 0;
    }
    
}

double lifetimeBaseAsymptoticGiantBranch(double timescales[nTimescales]){
    /*
     Calculate lifetime until Base of the Asymptotic Giant Branch (BAGB) given by tBAGB = tHeI + tHe
     
     Parameters
     -----------
     timescales : double array
     Array containing timescales
     
     timescales
     -----------
     0 : tMS  - Main sequence
     1 : tBGB - Base Giant Branch
     2 : tHeI - Helium ignition
     3 : tHe  - Helium burning
     #-- First Giant Branch (FGB)
     4 : Giant tinf1
     5 : Giant tinf2
     6 : Giant t(Mx)
     #-- Early Asymptotic Giant Branch (EAGB)
     7 : FAGB tinf1
     8 : FAGB tinf2
     9 : FAGB t(Mx)
     #-- Thermally Pulsating Asymptotic Giant Branch (TPAGB)
     10 : SAGB tinf1
     11 : SAGB tinf2
     12 : SAGB t(Mx)
     13 : TP (tDU?)
     14 : T(McMax)
     #-- not sure what TP, t(McMax) are
     #-- what about tHeMS etc?
     15 : tHeMS
     
     Returns
     --------
     tBAGB : double
     Lifetime until Base of the Asymptotic Giant Branch (BAGB)
     
     
     */
    return timescales[2] + timescales[3];
}

double calculateBluePhaseFBL(double bn_coefficients[nrowsb], double Mass, double coreMass, double massCutoffs[nMassCutoffs], double logMetallicityXi, unsigned long m_randomSeed){
    /*
     Calculate the parameter f_bl
     
     Given just after Equation 58 in Hurley et al 2000
     
     Parameters
     -----------
     bn_coefficients : double array
     Array containing metallicity dependent bn coefficients
     Mass : double
     Mass in Msol
     coreMass :double
     coreMass in Msol
     massCutoffs : double array
     Metallicity dependent mass cutoffs
     logMetallicityXi : double
     Metallciity
     
     Returns
     --------
     f_bl : double
     constant
     
     
     */
    
    bool debugging = false;
    //debugging = true;

    // Calculate RmHe for M > MFGB > MHeF
    double topTop    = bn_coefficients[24 - 1] * Mass + pow((bn_coefficients[25 - 1] * Mass), bn_coefficients[26 - 1]) * pow(Mass, bn_coefficients[28 - 1]);
    double topBottom = bn_coefficients[27 - 1] + pow(Mass, bn_coefficients[28 - 1]);
    double top       = topTop / topBottom;

    // Might be that we are supposed to use min(RmHe, Rx=RHeI)
    double RHeI = radiusHeliumIgnition(bn_coefficients, Mass, coreMass, massCutoffs, logMetallicityXi, m_randomSeed);
    // radiusCoreHeliumBurningX
    if(debugging){
        std::cout << "min(RmHe, Rx=RHeI) = min(" << top << " , " << RHeI << ")" << std::endl; 
    }
    top = std::min(top, RHeI);

    // Calculate RAGB(LHeI(M)) for M > MFGB > MHeF
    double bottom       = radiusAsymptoticGiantBranch(bn_coefficients, Mass, luminosityHeliumIgnition(bn_coefficients, Mass, massCutoffs, m_randomSeed), massCutoffs, m_randomSeed);
    double brackets     = 1.0 - (top / bottom);
	if(brackets ==0){brackets = 1e-12;}  //If zero gives R=NaN Coen Neijssel 10-01-2017
    double fbl = pow(Mass, bn_coefficients[48 - 1]) * pow(brackets, bn_coefficients[49 - 1]);
    
    if(debugging){
        std::cout << "f_bl brackets top/bottom = " << top << " / " << bottom << std::endl;
        std::cout << "f_bl = " << fbl << std::endl;
    }
	
    return fbl;
}

double lifetimeBluePhase(double bn_coefficients[nrowsb], double Mass, double coreMass, double massCutoffs[nMassCutoffs], double logMetallicityXi, unsigned long m_randomSeed){
    /*
     Calculate lifetime of the blue phase of core helium burning relative to tHe
     
     Given by Equation 58 in Hurley et al 2000
     
     Parameters
     -----------
     bn_coefficients : double array
     Array containing metallicity dependent bn coefficients
     Mass : double
     Mass in Msol
     coreMass : double
     coreMass in Msol
     massCutoffs : double array
     Metallicity dependent mass cutoffs
     logMetalliciityXi : double
     Metallicity
     
     Returns
     ---------
     tau_bl : double
     lifetime of blue phase relative to tHe : 0 <= tau_bl <= 1
     
     
     */
    
    bool debugging = false;
    //debugging = true;

    double MHeF = massCutoffs[1];
    double MFGB = massCutoffs[2];
    
    double f_bl_M       = calculateBluePhaseFBL(bn_coefficients, Mass, coreMass, massCutoffs, logMetallicityXi, m_randomSeed);
    double f_bl_MFGB    = calculateBluePhaseFBL(bn_coefficients, MFGB, coreMass, massCutoffs, logMetallicityXi, m_randomSeed);

    if(debugging){

        std::cout << "0) f_bl_M, f_bl_MFGB = " << f_bl_M << " " << f_bl_MFGB << std::endl;

    }

    double tbl  = 0;
    
    if(Mass < MHeF){
        tbl = 1.0;
    }
    else if(Mass >= MHeF and Mass <= MFGB){
        
        double firstTerm = bn_coefficients[45 - 1] * pow((Mass/MFGB), (0.414));
        
        double secondTerm = (1.0 - (bn_coefficients[45 - 1] * pow((Mass / MFGB),(0.414)))) * pow((log10(Mass/MFGB)/log10(MHeF/MFGB)), (bn_coefficients[46 - 1]));
        
        tbl = firstTerm + secondTerm;
        
        tbl = std::min(1.0, std::max(0.0, tbl));
        
        // DEBUGGING
        //std::cout << "tbl = " << tbl << std::endl;
        
    }
    else if(Mass > MFGB){
        tbl = (1.0 - bn_coefficients[47 - 1]) * (f_bl_M / f_bl_MFGB);

        if(debugging){

            std::cout << "1) Calculating lifetime of blue phase: " << tbl << std::endl;

        }

        tbl = std::min(1.0, std::max(0.0, tbl));
    }
    else{
        std::cerr << m_randomSeed << "\tError calculating blue phase lifetime" << std::endl;
        tbl = 1.0;
    }

    if(debugging){

        std::cout << "end) Calculating lifetime of blue phase: " << tbl << std::endl;

    }
    
    return tbl;
    
}

double luminosityCoreHeliumBurningX(double bn_coefficients[nrowsb], double Mass, double coreMass, double massCutoffs[nMassCutoffs], double logMetallicityXi, unsigned long m_randomSeed){
    /*
     Calculate the luminosity at the start of the blue phase of core helium burning
     
     Given by Equation 59 in Hurley et al 2000
     
     Parameters
     -----------
     Mass : double
     Mass in Msol
     coreMass : double
     Core mass in Msol
     
     Returns
     --------
     Lx : double
     Luminosity X in Lsol
     
     
     */
    
    double MHeF = massCutoffs[1];
    double MFGB = massCutoffs[2];
    
    if(Mass < MHeF){
        return luminosityZeroAgeHorizontalBranch(bn_coefficients, Mass, coreMass, massCutoffs, logMetallicityXi, m_randomSeed);
    }
    else if(Mass >= MHeF and Mass < MFGB){
        return minimumLuminosityCoreHeliumBurning(bn_coefficients, Mass, massCutoffs, logMetallicityXi, m_randomSeed);
    }
    else if(Mass >= MFGB){
        return luminosityHeliumIgnition(bn_coefficients, Mass, massCutoffs, m_randomSeed);
    }
    else{
        std::cerr << m_randomSeed << "\tError in calculating luminosity core helium burning X" << std::endl;
        return 0;
    }
}

double radiusCoreHeliumBurningX(double bn_coefficients[nrowsb], double Mass, double coreMass, double massCutoffs[nMassCutoffs], double logMetallicityXi, unsigned long m_randomSeed){
    /*
     Calculate the radius at the start of the blue phase of core helium burning
     
     Given by Equation 60 in Hurley et al 2000
     
     Parameters
     -----------
     bn_coefficients : double array
     Array containing metallicity dependent bn coefficients
     Mass : double
     Mass in Msol
     coreMass : double
     Core mass in Msol
     
     Returns
     --------
     Rx : double
     Radius X in Rsol
     
     
     */
    
    double MHeF = massCutoffs[1];
    double MFGB = massCutoffs[2];
    
    if(Mass < MHeF){
        return radiusZeroAgeHorizontalBranch(bn_coefficients, Mass, coreMass, massCutoffs, logMetallicityXi, m_randomSeed);
    }
    else if(Mass >= MHeF and Mass < MFGB){
        double LminHe = minimumLuminosityCoreHeliumBurning(bn_coefficients, Mass, massCutoffs, logMetallicityXi, m_randomSeed);
        return radiusFirstGiantBranch(bn_coefficients, Mass, LminHe, logMetallicityXi, m_randomSeed);
    }
    else if(Mass >= MFGB){
        return radiusHeliumIgnition(bn_coefficients, Mass, coreMass, massCutoffs, logMetallicityXi, m_randomSeed);
    }
    else{
        std::cerr << m_randomSeed << "\tError in calculating the radius core helium burning X" << std::endl;
        return 0;
    }
    
}

double relativeAgeStartBluePhase(double bn_coefficients[nrowsb], double Mass, double coreMass, double massCutoffs[nMassCutoffs], double logMetallicityXi, unsigned long m_randomSeed){
    /*
     Calculate the relative age at the start of the blue phase of Core Helium Burning
     
     Given just above Equation 59 in Hurley et al 2000
     
     Parameters
     -----------
     bn_coefficients : double array
     Array containing metallicity dependent bn coefficients
     Mass : double
     Mass in Msol
     coreMass : double
     coreMass in Msol
     massCutoffs : double array
     Array containing metallicity dependent mass cutoffs.
     logMetallicityXi : double
     Metallicity
     
     Returns
     --------
     tx : double
     Relative age at the start of the blue phase
     
     
     */
    double MHeF = massCutoffs[1];
    double MFGB = massCutoffs[2];
    
    if(Mass < MHeF){                            // For Low Mass (LM) stars
        return 0.0;                               // Really should be zero
    }
    else if(Mass >= MHeF and Mass < MFGB){      // For Intermediate Mass (IM) stars
        return 1.0 - lifetimeBluePhase(bn_coefficients, Mass, coreMass, massCutoffs, logMetallicityXi, m_randomSeed);
    }
    else if(Mass >= MFGB){                      // For High Mass (HM) stars
        return 0.0;                               // Really should be zero
    }
    else{
        std::cerr << m_randomSeed << "\tError in relativeAgeStartBluePhase" << std::endl;
        return 0.0;
    }
}

double relativeAgeEndBluePhase(double bn_coefficients[nrowsb], double Mass, double coreMass, double massCutoffs[nMassCutoffs], double logMetallicityXi, unsigned long m_randomSeed){
    /*
     Calculate the relative age at the end of the blue phase of core helium burning
     
     Given just above Equation 64 in Hurley et al 2000
     
     Parameters
     -----------
     bn_coefficients : double array
     Array containing metallicity dependent bn coefficients
     Mass : double
     Mass in Msol
     coreMass : double
     Core mass in Msol
     massCutoffs : double array
     Array containing metallicity dependent mass cutoffs
     logMetallicityXi : double
     Metallicity
     
     Returns
     --------
     ty : double
     Relative age at the end of the blue phase of core helium burning
     
     
     */
    
    double MHeF = massCutoffs[1];
    double MFGB = massCutoffs[2];
    
    if(Mass < MHeF){                                // For Low Mass (LM) stars
        return 1.0;
    }
    else if(Mass >= MHeF and Mass < MFGB){          // For Intermediate Mass (IM) stars
        return 1.0;
    }
    else if(Mass >= MFGB){                          // For High Mass (HM) stars
        return lifetimeBluePhase(bn_coefficients, Mass, coreMass, massCutoffs, logMetallicityXi, m_randomSeed);
    }
    else{
        std::cout << "" << std::endl;
        return 1.0;
    }
}

double luminosityCoreHeliumBurning(double bn_coefficients[nrowsb], double timescales[nTimescales], double massCutoffs[nMassCutoffs], double Mass, double coreMass, double tau, double logMetallicityXi, bool &m_error, unsigned long m_randomSeed){
    /*
     Calculate the luminosity during Core Helium Burning (CHeB)
     
     Given by Equations 61, 62 and 62 in Hurley et al 2000
     
     Parameters
     -----------
     bn_coefficients : double array
     Array containing metallicity dependent bn coefficients
     timescales : double array
     Array containing timescales
     massCutoffs : double array
     Array containing metallicity dependent mass cutoffs
     Mass : double
     Mass in Msol
     coreMass : double
     Core mass in Msol
     tau : double
     Relative time defined as tau = (time-tHeI)/tHe above Equation 59
     
     Returns
     --------
     L : double
     Luminosity during CHeB in Lsol
     
     
     */
    
    //double tx = relativeAgeStartBluePhase(bn_coefficients, Mass, massCutoffs);
    //double ty = relativeAgeEndBluePhase(bn_coefficients, Mass, massCutoffs);
    
    // Tolerance either side of start end in time
    double tolerance = 1E-2;
    
    // Get timescales
    double tx = timescales[20]; // Will be 0 for LM and HM stars, non-zero for IM stars
    
    // Calculate radii
    double RminHe   = minimumRadiusCoreHeliumBurning(bn_coefficients, Mass, coreMass, massCutoffs, logMetallicityXi, m_randomSeed);
    double Rx       = radiusCoreHeliumBurningX(bn_coefficients, Mass, coreMass, massCutoffs, logMetallicityXi, m_randomSeed);
    
    // Calculate luminosities
    double Lx       = luminosityCoreHeliumBurningX(bn_coefficients, Mass, coreMass, massCutoffs, logMetallicityXi, m_randomSeed);
    double LBAGB    = luminosityBaseAsymptoticGiantBranch(bn_coefficients, Mass, massCutoffs, m_randomSeed);
    double LHeI     = luminosityHeliumIgnition(bn_coefficients, Mass, massCutoffs, m_randomSeed);
    double RHeI     = radiusHeliumIgnition(bn_coefficients, Mass, coreMass, massCutoffs, logMetallicityXi, m_randomSeed);
    
    bool debugging = false;
    
    if(debugging){
        
        std::cout << "Lx = " << Lx << std::endl;
        std::cout << "RHeI = " << RHeI << std::endl;
        
    }
    
    if(tau >= tx and tau <= 1.0){
        double epsilon      = std::min(2.5, std::max(0.4, (RminHe/Rx)));
        double lambdaTop    = tau - tx;
        double lambdaBottom = 1.0 - tx;
        double lambda       = pow((lambdaTop/lambdaBottom), epsilon);
        double L = Lx * pow((LBAGB/Lx), lambda);
        //std::cout << "1) L = " << L << std::endl;
        return L;
    }
    else if(tau >= 0 and tau < tx){
        double lambdaPrimeTop       = tx - tau;
        double lambdaPrimeBottom    = tx; // can be zero, divide by zero may be a problem.
        if(lambdaPrimeBottom <= 0){
            if(m_error == false){
                std::cerr << m_randomSeed << "\tError: Divide by zero in luminosityCoreHeliumBurning" << std::endl;
            }
           std::cerr << m_randomSeed << "\tError: Divide by zero in luminosityCoreHeliumBurning" << std::endl;
           m_error = true;
        }
        double lambdaPrime          = pow((lambdaPrimeTop / lambdaPrimeBottom), 3.0);
        double L = Lx * pow((LHeI/Lx), lambdaPrime);
        //std::cout << "2) L = " << L << std::endl;
        return L;
    }
    else if(tau > 1.0 and tau < 1.0 + tolerance){
        
        // Should really be tau = 1.0
        tau = 1.0;
        
        double epsilon      = std::min(2.5, std::max(0.4, (RminHe/Rx)));
        double lambdaTop    = tau - tx;
        double lambdaBottom = 1.0 - tx;
        double lambda       = pow((lambdaTop/lambdaBottom), epsilon);
        double L = Lx * pow((LBAGB/Lx), lambda);
        std::cout << "3) L = " << L << std::endl;
        return L;
        
    }
    else if(tau < 0.0 and tau > 0.0 - tolerance){
        
        // Luminosity error fix
        // Should really be tau = 0.0
        tau = 0.0;
        
        double lambdaPrime          = 0.0;
        double lambdaPrimeTop       = tx - tau;
        double lambdaPrimeBottom    = tx; // can be zero, divide by zero may be a problem.
        
        if(tau == 0.0 and tx == 0.0){
            
            lambdaPrime = 0.0;
            
        }
        else{
            if(lambdaPrimeBottom <= 0){
                if(m_error == false){
                    std::cerr << m_randomSeed << "\tError: Divide by zero in luminosityCoreHeliumBurning" << std::endl;
                }
                std::cerr << m_randomSeed << "\tError: Divide by zero in luminosityCoreHeliumBurning" << std::endl;
                m_error = true;
            }
            lambdaPrime          = pow((lambdaPrimeTop / lambdaPrimeBottom), 3.0);
        }
        
        double L = Lx * pow((LHeI/Lx), lambdaPrime);
        
        if(debugging){
            std::cout << "4) L = " << L << std::endl;
            std::cout << "LHeI, Lx = " << LHeI << " " << Lx << std::endl;
            std::cout << "Parameters: tau, tx, Mass, coreMass = "  << tau << " " << tx << " " << Mass << " " << coreMass << std::endl;
        }
        
        return L;
        
    }
    else{
        double L = 0;
        if(m_error == false){
            std::cerr << m_randomSeed << "\tError in luminosityCoreHeliumBurning: tau, tx, mass, coreMass = " << tau << " " << tx << " " << Mass << " " << coreMass << std::endl;
        }
        std::cerr << m_randomSeed << "\tError in luminosityCoreHeliumBurning: tau, tx, mass, coreMass = " << tau << " " << tx << " " << Mass << " " << coreMass << std::endl;
        m_error = true;
        return L;
    }
    
}

double radiusCoreHeliumBurningRho(double bn_coefficients[nrowsb], double timescales[nTimescales], double tau, double Mass, double coreMass, double massCutoffs[nMassCutoffs], double logMetallicityXi, bool &m_error, unsigned long m_randomSeed, double Mass0){
    /*
     Calculate the parameter Rho for the core helium burning radius
     
     Given by Equation 65 in Hurley et al 2000
     
     Parameters
     -----------
     bn_coefficients : double array
     Array containing metallicity dependent bn coefficients
     timescales : double array
     Array containing metallicity dependent timescales
     tau : double
     Relative time in CHeB defined as tau = (time - tHeI)/tHe above Eq 59
     Mass : double
     Mass in Msol
     coreMas : double
     Core mass in Msol
     massCutoffs : double array
     Array containing metallicity dependent mass cutoffs
     
     Returns
     --------
     rho : double
     Rho parameter
     
     
     */
    bool debugging = false;
    //debugging=true;

    //std::cout << "new radiusCoreHeliumBurningRho" << std::endl; // DEBUGGING
    
    // aren't these calculated by the function that calls this?
    //double tx   = relativeAgeStartBluePhase(bn_coefficients, Mass, massCutoffs);
    //double ty   = relativeAgeEndBluePhase(bn_coefficients, Mass, massCutoffs);
    
    // Get timescales
    double tx = timescales[20];
    double ty = timescales[21];
    
    // Calculate radii
    double RmHe = minimumRadiusCoreHeliumBurning(bn_coefficients, Mass, coreMass, massCutoffs, logMetallicityXi, m_randomSeed);
    double Rx   = radiusCoreHeliumBurningX(bn_coefficients, Mass, coreMass, massCutoffs, logMetallicityXi, m_randomSeed);
    double Rmin = std::min(RmHe, Rx);
    
    //std::cout << "tx, ty,  tau = " << tx << " " << ty << " " << tau << std::endl;
    double Ly = luminosityCoreHeliumBurning(bn_coefficients, timescales, massCutoffs, Mass0, coreMass, tau, logMetallicityXi, m_error, m_randomSeed);
    // used to use Ly to calculate Ry
    double Ry = radiusCHeBY(bn_coefficients, timescales, Mass, coreMass, massCutoffs, logMetallicityXi, m_error, m_randomSeed, Mass0);
    //double Ry = radiusAsymptoticGiantBranch(bn_coefficients, Mass, Ly, massCutoffs, m_randomSeed);
    //std::cout << "Rx, Ry = " << Rx << " " << Ry << std::endl;
    
    // Define the four bracketed terms indivually for improved readability
    double one      = pow(log((Ry / Rmin)), (1.0/3.0));
    double two      = (tau - tx)/(ty - tx);
    double three    = pow(log((Rx / Rmin)), (1.0/3.0));
    double four     = (ty - tau)/(ty - tx);
    
    double rho = (one * two) - (three * four);
    
    if(debugging){
        std::cout << "one=" << one << " tw0=" << two << " three=" << three << " four=" << four << " Ry="<<Ry<<" Rx="<<Rx<<std::endl;
    }

    return rho;
    
}

double radiusCHeBY(double bn_coefficients[nrowsb], double timescales[nTimescales], double Mass, double coreMass, double massCutoffs[nMassCutoffs], double logMetallicityXi, bool &m_error, unsigned long m_randomSeed, double Mass0){
    /*
     Calculate the radius parameter Ry
     
     Given just above Equation 64 in Hurley et al 2000
     
     Parameters
     -----------
     bn_coefficients : double array
     Array containing metallicity dependent bn coefficients
     timescales : double array
     Array containing timescales
     Mass : double
     Mass in Msol
     coreMass : double
     Core mass in Msol
     massCutoffs : double array
     Array containing metallicity dependent mass cutoffs
     logMetallicty : double
     Metallicity
     
     Returns
     --------
     Ry : double
     Radius parameter Ry
     
     
     */
    
    // Initialise variables
    double tx  = 0.0;
    double ty  = 0.0;
    double tau = 0.0;
    
    double tbl = lifetimeBluePhase(bn_coefficients, Mass, coreMass, massCutoffs, logMetallicityXi, m_randomSeed);
    
    // For LM stars tx = 0, ty = 1 and Ry = R(ty)
    if(Mass <= massCutoffs[1]){ // M <= MFGB
        ty  = 1.0;
        tx  = 0.0;    // Actually zero
        tau = ty;
        
    }
    // For IM stars tx = 1 - tbl, ty = 1.0, Ry = R(ty)
    else if(Mass > massCutoffs[1] and Mass <= massCutoffs[2]){
        ty  = 1.0;
        tx  = 1.0 - tbl; // 1 - tbl
        tau = ty;
        
    }
    // For HM stars Ly given by Equation 61 (Ly = LBAGB for M > MFGB) and Ry = RAGB(Ly)
    else if(Mass > massCutoffs[2]){ // M > MFGB
        tx  = 0.0; // Actually zero
        ty  = tbl; //tbl;
        tau = ty;
        
    }
    else{
        tx  = 0.0;
        ty  = 0.0;
        tau = ty;
        std::cerr << m_randomSeed << "\tError in radiusCHeBY" << std::endl;
		std::cout << "Check error flag" << std::endl;
        m_error = false;
        
    }
    
    double RminHe = minimumRadiusCoreHeliumBurning(bn_coefficients, Mass, coreMass, massCutoffs, logMetallicityXi, m_randomSeed);
    double Rx     = radiusCoreHeliumBurningX(bn_coefficients, Mass, coreMass, massCutoffs, logMetallicityXi, m_randomSeed); // In GBParams FINISH ME
    double Rmin   = std::min(RminHe, Rx);
    
    double Ly     = luminosityCHeBY(bn_coefficients, timescales, Mass0, coreMass, massCutoffs, logMetallicityXi, m_error, m_randomSeed);
    double Ry     = 0.0;   // Initialise Ry variable
    
    if(tau >= 0 and tau < tx){
        Ry = radiusFirstGiantBranch(bn_coefficients, Mass, Ly, logMetallicityXi, m_randomSeed);
        //std::cout << "RCHEBY 1) Ry = " << Ry << std::endl;
    }
    else if(tau >= ty and tau <= 1.0){ // added an equals sign here to the condition to just use RAGB
        Ry = radiusAsymptoticGiantBranch(bn_coefficients, Mass, Ly, massCutoffs, m_randomSeed);
        //std::cout << "RCHEBY 2) Ry = " << Ry << std::endl;
    }
    else if(tau >= tx and tau <= ty){
        double rho = 0;
        Ry = Rmin * exp(pow(fabs(rho), 3.0));
        //std::cout << "RCHEBY 3) Ry requires rho - not defined Ry = " << Ry << std::endl;
    }
    else{
        Ry = 0;
        std::cerr << m_randomSeed << "\tError in radiusCoreHeliumBurning, Ry = " << Ry << std::endl;
    }

    // DEBUGGING
    bool debugging = false;
    // debugging = true;
    if(debugging){
        std::cout << "RCHEBY : tx, ty, tau, Ry = " << tx << " " << ty << " " << tau << " " << Ry << std::endl;
    }
    
    return Ry;
    
}

double luminosityCHeBY(double bn_coefficients[nrowsb], double timescales[nTimescales], double Mass, double coreMass, double massCutoffs[nMassCutoffs], double logMetallicityXi, bool &m_error, unsigned long m_randomSeed){
    /*
     Calculate the luminosity constant Ly = L(tau = tau_y) during core helium burning
     
     Given by Equation 61, 62, 63 in Hurley et al 2000 and just before Equation 64
     
     Parameters
     -----------
     bn_coefficients : double array
     Array containing metallicity dependent bn coefficients
     timescales : double array
     Array containing timescales
     Mass : double
     Mass in Msol
     coreMass : double
     Core mass in Msol
     massCutoffs : double array
     Array containing metallicity dependent mass cutoffs
     
     Returns
     --------
     Ly : double
     Luminosity Ly during core helium burning (CHeB)
     
     
     */
    
    double Ly = 0;
    
    //double tau_y = relativeAgeEndBluePhase(bn_coefficients, Mass, massCutoffs);
    //double tau_x = relativeAgeStartBluePhase(bn_coefficients, Mass, massCutoffs);
    
    // Get timescales
    double tau_x = timescales[20];
    double tau_y = timescales[21];
    
    //DEBUGGING
    //std::cout << "Calc LY tx, ty = " << tau_x << " " << tau_y << std::endl;
    
    // Calculate radii
    double RminHe = minimumRadiusCoreHeliumBurning(bn_coefficients, Mass, coreMass, massCutoffs, logMetallicityXi, m_randomSeed);
    double Rx = radiusCoreHeliumBurningX(bn_coefficients, Mass, coreMass, massCutoffs, logMetallicityXi, m_randomSeed);
    
    // Calculate luminosities (could do up front)
    double Lx    = luminosityCoreHeliumBurningX(bn_coefficients, Mass, coreMass, massCutoffs, logMetallicityXi, m_randomSeed);
    double LBAGB = luminosityBaseAsymptoticGiantBranch(bn_coefficients, Mass, massCutoffs, m_randomSeed);
    double LHeI  = luminosityHeliumIgnition(bn_coefficients, Mass, massCutoffs, m_randomSeed);
    
    //tau = tau_y
    if(tau_y >= tau_x and tau_y <= 1.0){
        double epsilon      = std::min(2.5, std::max(0.4, (RminHe/Rx)));
        double lambdaTop    = tau_y - tau_x;
        double lambdaBottom = 1.0 - tau_x;
        double lambda       = pow((lambdaTop/lambdaBottom), epsilon);
        Ly = Lx * pow((LBAGB/Lx), lambda);
        //std::cout << "1) Ly = " << Ly << std::endl;
    }
    else if(tau_y >= 0 and tau_y < tau_x){
        double lambdaPrimeTop       = tau_x - tau_y;
        double lambdaPrimeBottom    = tau_x; // can be zero, divide by zero may be a problem.
        if(lambdaPrimeBottom <= 0){
            std::cout << "Divide by zero in luminosityCoreHeliumBurning" << std::endl;
        }
        double lambdaPrime          = pow((lambdaPrimeTop / lambdaPrimeBottom), 3.0);
        Ly = Lx * pow((LHeI/Lx), lambdaPrime);
        //std::cout << "2) Ly = " << Ly << std::endl;
    }
    else{
        Ly = 0;
 
        if(m_error == false){
           std::cerr << m_randomSeed << "\tError in luminosityCoreHeliumBurning: Ly, tau_y, tau_x = " << Ly << " " << tau_y << " " << tau_x << std::endl;
        }
        std::cerr << m_randomSeed << "\tError in luminosityCoreHeliumBurning: Ly, tau_y, tau_x = " << Ly << " " << tau_y << " " << tau_x << std::endl;
        m_error = true;

    }
    
    return Ly;
    
}

double radiusCoreHeliumBurning(double bn_coefficients[nrowsb], double massCutoffs[nMassCutoffs], double timescales[nTimescales], double tau, double Mass, double coreMass, double Luminosity, double logMetallicityXi, bool &m_error, unsigned long m_randomSeed, double Mass0){
    /*
     Calculate the radius during core helium burning
     
     Given by Equation 64 in Hurley et al 2000
     
     Parameters
     -----------
     bn_coefficients : double array
     Array containing metallicity dependent bn coefficients
     massCutoffs : double array
     Array containing metallicity dependent mass cutoffs
     timescales : double array
     Array containing timescales
     tau : double
     Relative time in MYRs defined as tau = (time - tHeI)/tHe just above Equation 59
     Mass : double
     Mass in Msol
     coreMass : double
     Core Mass in Msol
     Luminosity : double
     Luminosity in Lsol
     
     Returns
     --------
     R_CHeB : double
     Radius during Core Helium Burning (CHeB)
     
     
     */
   	bool	debugging = false;
	//debugging = true;

    // Calculate radii
    double RmHe   = minimumRadiusCoreHeliumBurning(bn_coefficients, Mass, coreMass, massCutoffs, logMetallicityXi, m_randomSeed);
    double Rx     = radiusCoreHeliumBurningX(bn_coefficients, Mass, coreMass, massCutoffs, logMetallicityXi, m_randomSeed); // IN GBPARAMS FINISH ME
    double Rmin   = std::min(RmHe, Rx);
    double R_CHeB = 0.0;
    
	if(debugging){
        std::cout << "Radius CHeB" << std::endl;
		std::cout << "Tau: " << tau << std::endl;
		std::cout << "RCHEB: RmHe = " << RmHe << std::endl;
		std::cout << "RCHEB: Rx = " << Rx << std::endl;
		std::cout << "RCHEB: Rmin = " << Rmin << std::endl;
    }
    // These aren't stored in timescales (even in fortran code) even though they probably should be
    //double tx   = relativeAgeStartBluePhase(bn_coefficients, Mass, massCutoffs);
    //double ty   = relativeAgeEndBluePhase(bn_coefficients, Mass, massCutoffs);
    

    double tx = timescales[20];
    double ty = timescales[21];
    

	if(debugging){
		std::cout << "RCHEB: tx, ty, tau = " << tx << " " << ty << " " << tau << std::endl;
    }
	
    double rho  = 0.0; // Initialise variable for rho
    
    if(tau >= 0 and tau < tx){
        R_CHeB = radiusFirstGiantBranch(bn_coefficients, Mass, Luminosity, logMetallicityXi, m_randomSeed);
		if(debugging){
			std::cout << "1) tau, R_CHeB, log10(R_CHEB) = " << tau << " " << R_CHeB << " " << log10(R_CHeB) << std::endl;
		}
    }
    else if(tau > ty and tau <= 1.0){
        R_CHeB = radiusAsymptoticGiantBranch(bn_coefficients, Mass, Luminosity, massCutoffs, m_randomSeed);
		if(debugging){
			std::cout << "2) tau, R_CHeB, log10(R_CHEB) = " << tau << " " << R_CHeB << " " << log10(R_CHeB) << std::endl;
		}
    }
    else if(tau >= tx and tau <= ty){
        //double a = 48.3395;
        //double b = -1.18817;
        //std::cout << "test with Rmin, rho = " << Rmin << " " << rho << " = " << log10(a * exp(pow(fabs(b), 3.0))) << std::endl;
        rho = radiusCoreHeliumBurningRho(bn_coefficients, timescales, tau, Mass, coreMass, massCutoffs, logMetallicityXi, m_error, m_randomSeed, Mass0);
        R_CHeB =  Rmin * exp(pow(fabs(rho), 3.0));
        //std::cout << "rho = " << rho << std::endl;
        // Apparently, should use fabs instead of abs
		if(debugging){
			std::cout << "3) tau, R_CHeB, log10(R_CHEB), rho = " << tau << " " << R_CHeB << " " << log10(R_CHeB) << " " << rho << std::endl;
		}
        if(R_CHeB != R_CHeB){
            std::cerr << m_randomSeed << "\tR_CHeB is NAN" << std::endl;
            //R_CHeB = 0.0;
        }
    }
    else{
        R_CHeB = 0.0;
		if(debugging){
			std::cout << "4) Error in radiusCoreHeliumBurning" << std::endl;
		}
    }
    
    return R_CHeB;
    
}

double coreMassHeliumBurning(double bn_coefficients[nrowsb], double GBParams[nGBParams], double tau, double Mass, double massCutoffs[nMassCutoffs], double logMetallicityXi, unsigned long m_randomSeed){
    /*
     Calculate the core mass between helium ignition ans the base of the asymptotic
     giant branch
     
     Given by Equation 67 in Hurley et al 2000
     
     Parameters
     -----------
     bn_coefficients : double array
     Array containing metallicity dependent bn coefficients
     GBParams : double array
     Array containing GB parameters
     tau : double
     Relative age defined as tau = (time - tHeI)/tHe above Equation 59
     Mass : double
     Mass in Msol
     massCutoffs : double array
     Array containing metallicity dependent mass cutoffs
     logMetallicityXi : double
     logMetallicityXi
     
     Returns
     --------
     Mc : double
     Core mass during helium burning in Msol
     
     
     */
    
    double McHeI  = coreMassHeliumIgnition(bn_coefficients, GBParams, Mass, massCutoffs, logMetallicityXi, m_randomSeed);
    double McBAGB = coreMassBaseAsymptoticGiantBranch(bn_coefficients, Mass, m_randomSeed);
    
    return (1.0 - tau)*McHeI + tau*McBAGB;
}

double coreMassSecondDredgeUp(double GBParams[nGBParams], unsigned long m_randomSeed){
    /*
     Calculate the mass of the core during dredge up
     
     Given by Equation 69 in Hurley et al 2000
     
     Parameters
     ----------
     GBParams : double array
     Array containing giant branch parameters
     
     Returns
     --------
     McDU : double
     Core Mass at dredge up
     
     
     */
    double McBAGB = GBParams[9]; // Get McBAGB
    double McDU   = 0.44*McBAGB + 0.448; // Declare variable for McDU
    
    if(McBAGB <= 0.8){
        McDU = McBAGB;
    }
    
    return McDU;
}

double lifetimeSecondDredgeUp(double bn_coefficients[nrowsb], double timescales[nTimescales], double GBParams[nGBParams], double Mass, double COCoreMass, double massCutoffs[nMassCutoffs], double logMetallicityXi, unsigned long m_randomSeed){
    /*
     Calculate the lifetime til second dredge up (nb there is no explicit first)
     This is the time of transition between the EAGB and the TPAGB
     
     Given by Equation 70, 71, 72 in Hurley et al 2000
     
     Parameters
     -----------
     bn_coefficients : double array
     Array containing metallicity dependent bn coefficients
     timescales : double array
     Array containing timescales
     GBParams : double array
     Array containing GB params
     Mass : double
     Mass in Msol
     COCoreMass : double
     CO Core Mass in Msol
     massCutoffs : double array
     Array containing metallicity dependent mass cutoffs
     logMetallicityXi : double
     logMetallicityXi
     
     timescales
     ------------
     List of timescales (all in units of Myrs)
     0 : tMS  - Main sequence
     1 : tBGB - Base Giant Branch
     2 : tHeI - Helium ignition
     3 : tHe  - Helium burning
     #-- First Giant Branch (FGB)
     4 : Giant tinf1
     5 : Giant tinf2
     6 : Giant t(Mx)
     #-- Early Asymptotic Giant Branch (EAGB)
     7 : FAGB tinf1
     8 : FAGB tinf2
     9 : FAGB t(Mx)
     #-- Thermally Pulsating Asymptotic Giant Branch (TPAGB)
     10 : SAGB tinf1
     11 : SAGB tinf2
     12 : SAGB t(Mx)
     13 : TP (tDU?)
     14 : T(McMax)
     #-- not sure what TP, t(McMax) are
     #-- what about tHeMS etc?
     15 : tHeMS
     
     GBParams
     ---------
     List of GB Parameters:
     0 : A_H
     1 : A_H,He
     2 : A_He
     3 : B
     4 : D
     5 : p
     6 : q
     7 : Mx
     8 : McBGB
     9 : McBAGB
     10 : McDU
     
     Returns
     --------
     tDU : double
     lifetime til dredge up
     
     
     */
    
    // Declare variable for second dredge up time
    double tDU = 0;
    
    // Get parameters
    double AHe  = GBParams[2];
    double B    = GBParams[3];
    double D    = GBParams[4];
    double p    = GBParams[5];
    double q    = GBParams[6];
    //double Mx   = GBParams[7];
    double McDU = GBParams[10];
    
    double tinf1 = timescales[7]; // Given by Equation 40
    double tinf2 = timescales[8]; // Given by Equation 42
    
    // Calculate luminosities
    double LDU   = coreMassLuminosityRelation(GBParams, McDU, m_randomSeed);
    double Lx    = calculateLuminosityLx(GBParams);
    
    if(LDU <= Lx){
        tDU = tinf1 - ((1.0)/((p - 1.0) * AHe * D)) * pow((D / LDU), ((p - 1.0)/(p))); // Equation 70
    }
    else if(LDU > Lx){
        tDU = tinf2 - ((1.0)/((q - 1.0) * AHe * B)) * pow((B / LDU), ((q - 1.0)/(q))); // Equation 70
    }
    else{
        std::cerr << m_randomSeed << "\tError calculating the lifetime of second dredge up" << std::endl;
    }
    
    return tDU;
    
}

double coreMassAsymptoticGiantBranch(double GBParams[nGBParams], double timescales[nTimescales], double time, unsigned long m_randomSeed){
    /*
     Calculate the CO core mass as a function of time on the asymptotic giant branch
     
     Given by Equation 39 in Hurley et al 2000 modified as described in Section 5.4
     
     Parameters
     -----------
     GBParams : double array
     Array containing GB parameters
     timescales : double array
     Array containing timescales
     time : double
     Time in MYRs
     
     GBParams
     ---------
     List of GB Parameters:
     0 : A_H
     1 : A_H,He
     2 : A_He
     3 : B
     4 : D
     5 : p
     6 : q
     7 : Mx
     8 : McBGB
     9 : McBAGB
     10 : McDU
     
     Timescales
     -----------
     List of timescales (all in units of Myrs)
     0 : tMS  - Main sequence
     1 : tBGB - Base Giant Branch
     2 : tHeI - Helium ignition
     3 : tHe  - Helium burning
     #-- First Giant Branch (FGB)
     4 : Giant tinf1
     5 : Giant tinf2
     6 : Giant t(Mx)
     #-- Early Asymptotic Giant Branch (EAGB)
     7 : FAGB tinf1
     8 : FAGB tinf2
     9 : FAGB t(Mx)
     #-- Thermally Pulsating Asymptotic Giant Branch (TPAGB)
     10 : SAGB tinf1
     11 : SAGB tinf2
     12 : SAGB t(Mx)
     13 : TP (tDU?)
     14 : T(McMax)
     #-- not sure what TP, t(McMax) are
     #-- what about tHeMS etc?
     15 : tHeMS
     #-- Heium giant branch
     16 : HeGB tinf1
     17 : HeGB tinf2
     18 : HeGB tx
     19 : Relative duration of blue loop taubl
     20 : Relative start of blue loop taux
     21 : Relative end of blue loop tauy
     
     Returns
     --------
     McAGB : double // McCO?
     Core mass on the Asymptotic Giant Branch (AGB)
     
     
     */
    
    // Get GB Parameters
    double AHe = GBParams[2];
    double B   = GBParams[3];
    double D   = GBParams[4];
    double p   = GBParams[5];
    double q   = GBParams[6];
    
    // Get timescales
    double tinf1    = timescales[7];
    double tx       = timescales[9];
    double tinf2    = timescales[8];
    
    // DEBUGGING
    //std::cout << "McAGB: time, tx = " << time << " " << tx << std::endl;
    
    double McAGB = 0.0;
    
    // Given by Equation 39
    if(time <= tx){
        McAGB = pow(((p - 1.0) * AHe * D * (tinf1 - time)), (1.0/(1.0 - p)));
    }
    else{
        McAGB = pow(((q - 1.0) * AHe * B * (tinf2 - time)), (1.0/(1.0 - q)));
    }
    
    //std::cout << "McAGB: McAGB = " << McAGB << std::endl;
    
    return McAGB;
    
}

double luminosityAsymptoticGiantBranch(double GBParams[nGBParams], double coreMass){
    /*
     Calculate the luminosity on the Asymptotic Giant Branch (AGB) using the core mass - luminosity relation for AGB stars.
     
     Given by Equation 37 of Hurley et al 2000 modified as described at the begining of Section 5.4
     
     Parameters
     -----------
     GBParams : double array
     Array containing GB Params
     coreMass : double
     Core mass in Msol
     
     GBParams
     ---------
     List of GB Parameters:
     0 : A_H
     1 : A_H,He
     2 : A_He
     3 : B
     4 : D
     5 : p
     6 : q
     7 : Mx
     8 : McBGB
     9 : McBAGB
     10 : McDU
     
     Returns
     ---------
     L_AGB : double
     Luminosity on the AGB in Lsol
     
     
     */
    
    double B = GBParams[3];
    double D = GBParams[4];
    double p = GBParams[5];
    double q = GBParams[6];
    
    return std::min((B * pow(coreMass, q)), (D * pow(coreMass, p)));
}

double calculateTPAGBLambda(double mass){
    /*
     Calculate the fraction
     
     Given by Equation 73 in Hurley et al 2000
     
     Parameters
     -----------
     Mass : double
     Mass in Msol
     
     Returns
     --------
     lambda : double
     Fraction of core mass eaten during thermal pulse
     
     
     */
    return std::min(0.9, 0.3 + 0.001 * pow(mass, 5.0));
}

double coreMassPrimeThermallyPulsatingAsymptoticGiantBranch(double GBParams[nGBParams], double timescales[nTimescales], double time, unsigned long m_randomSeed){
    /*
     Parameters
     -----------
     GBParams : double array
     Array containing GB parameters
     timescales : double array
     Array containing timescales
     time : double
     
     
     Returns
     --------
     McPrime : double
     Dummy core mass used to calculate luminosity and actual core mass on TPAGB
     
     
     */
    
    double McPrime = 0;
    
    // Get parameters
    double AHHe   = GBParams[1];
    double B      = GBParams[3];
    double D      = GBParams[4];
    double p      = GBParams[5];
    double q      = GBParams[6];
    
    double tinf1  = timescales[10];
    double tinf2  = timescales[11];
    double tx     = timescales[12];
    
    // McPrime used as a dummy variable to calculate Mc (After Equation 73)
    if(time <= tx){
        //std::cout << "1) " << std::endl; // DEBUGGING
        McPrime = pow(((p - 1.0) * AHHe * D * (tinf1 - time)), (1.0/(1.0 - p))); // Equation 39
    }
    else{
        //std::cout << "2) " << std::endl; // DEBUGGING
        McPrime = pow(((q - 1.0) * AHHe * B * (tinf2 - time)), (1.0/(1.0 - q))); // Equation 39
    }
    
    return McPrime;
    
}

double coreMassThermallyPulsatingAsymptoticGiantBranch(double McPrime, double GBParams[nGBParams], double Mass, unsigned long m_randomSeed){
    /*
     Calculate core mass on TPAGB
     
     Parameters
     -----------
     McPrime : double
     Dummy core mass used for calculating luminosity and core mann on TPAGB in Msol
     Mass : double
     Mass in Msol
     
     Returns
     --------
     Mc : double
     Core mass in Msol
     
     
     */
    
    double lambda = calculateTPAGBLambda(Mass);
    
    double Mc = McPrime;
    
    double McDU = GBParams[10];
    
    Mc = McDU +  (1.0- lambda) * (McPrime - McDU);
    
    return Mc;
}

double luminosityThermallyPulsatingAsymptoticGiantBranch(double GBParams[nGBParams], double McPrime, unsigned long m_randomSeed){
    /*
     Calculate luminosity on TPAGB
     
     Parameters
     -----------
     GBParams : double array
     Array containing GB parameters
     McPrime : double
     Dummy core mass used to calculate luminosity in Msol
     
     Returns
     --------
     LTPAGB : double
     Luminosity on TPAGB
     
     
     */
    return coreMassLuminosityRelation(GBParams, McPrime, m_randomSeed);
}

double luminosityHeliumMainSequenceZAMS(double Mass){
    /*
     Calculate the ZAMS luminosity of a naked helium star with Z = 0.02
     
     Given by Equation 77 of Hurley et al 2000
     
     Parameters
     ------------
     Mass : double
     Mass in Msol
     
     Returns
     --------
     LZHe : double
     Luminosity of naked helium star in Lsol
     
     
     */
    
    double top      = 15262.0 * pow(Mass, 10.25);
    double bottom   = pow(Mass, 9.0) + 29.54 * pow(Mass, 7.5) + 31.18 * pow(Mass, 6.0) + 0.0469;
    
    return top / bottom;
}

double radiusHeliumMainSequenceZAMS(double Mass){
    /*
     Calculate the ZAMS radius of a naked helium star with Z = 0.02
     
     Given by Equation 78 of Hurley et al 2000
     
     Parameters
     -----------
     Mass : double
     Mass in Msol
     
     Returns
     --------
     RZHe : double
     Radius of naked helium star in Rsol
     
     
     */
    
    double top      = 0.2391 * pow(Mass, 4.6);
    double bottom   = pow(Mass, 4.0) + 0.162 * pow(Mass, 3.0) + 0.0065;
    
    return top / bottom;
}

double lifetimeHeliumMainSequence(double Mass){
    /*
     Calculate the Helium Main Sequence (HeMS) lifetime of a naked helium star with Z = 0.02
     
     Given by Equation 79 in Hurley et al 2000
     
     Parameters
     -----------
     Mass : double
     Mass in Msol
     
     Returns
     --------
     tHeMS : double
     lifetime on HeMS
     
     
     */
    double top      = 0.4129 + 18.81 * pow(Mass, 4.0) + 1.853 * pow(Mass, 6.0);
    double bottom   = pow(Mass, 6.5);
    
    return top / bottom;
}

double luminosityHeliumMainSequence(double tau, double Mass){
    /*
     Calculate the luminosity as a function of time during central He burning
     
     Given by Equations 80 and 82 in Hurley et al 2000
     
     Parameters
     -----------
     tau : double
     Relative age wrt tHeMS should be in [0,1]
     Mass : double
     Mass in Msol
     
     Returns
     --------
     LHeMS : double
     Luminosity during central He burning in Lsol
     
     
     */
    
    double LZHe     = luminosityHeliumMainSequenceZAMS(Mass);
    
    double alpha    = std::max(0.0, 0.85 - 0.08 * Mass);
    
    return LZHe * (1.0 + 0.45 * tau + alpha * tau * tau);
}

double radiusHeliumMainSequence(double tau, double Mass){
    /*
     Parameters
     ------------
     tau : double
     Relative age wrt tHeMS should be [0,1]
     Mass : double
     Mass in Msol
     
     Returns
     --------
     RHeMS : double
     Radius during central He burning in Rsol
     
     
     */
    
    double RZHe     = radiusHeliumMainSequenceZAMS(Mass);
    
    double beta     = std::max(0.0, 0.4 - 0.22 * log10(Mass));
    
    return RZHe * (1.0 + beta * tau - beta * pow(tau, 6.0));
}

double luminosityEndHeliumMainSequence(double Mass){
    /*
     Calculate luminosity at the end of the helium main sequence
     
     Given by Equation 80 at tau = 1
     
     Parameters
     -----------
     Mass : double
     Mass in Msol
     
     Returns
     --------
     LTHe : double
     Luminosity at the end of the helium main sequence in Lsol
     
     
     */
    double alpha    = std::max(0.0, 0.85 - 0.08 * Mass);
    
    double LZHe     = luminosityHeliumMainSequenceZAMS(Mass);
    
    return LZHe * (1.0 + 0.45 + alpha);
}




double radiusEndHeliumMainSequence(double Mass){
    /*
     Calculate the radius at the end of the helium main sequence
     
     Given by Equation 81 at tau = 1
     
     Parameters
     ------------
     Mass : double
     Mass in Msol
     
     Returns
     ---------
     RTHe : double
     Radius at the end of the helium main sequence
     
     
     */
    double RZHe = radiusHeliumMainSequenceZAMS(Mass);
    return RZHe;
}

double ageEvolvedHeliumGiantBranch(double GBParams[nGBParams], double timescales[nTimescales], double Mass, double coreMass){
    /*
     Calculate the age of an evolved Helium Giant branch star
     
     Given by core mass - t relation (Equation 39) modified as described after Equation 84
     
     Parameters
     -----------
     GBParams : double array
     Array containing GB parameters
     Mass : double
     Mass in Msol
     coreMass : double
     Core mass in Msol
     
     GBParams
     ---------
     0  : A_H
     1  : A_H,He
     2  : A_He
     3  : B
     4  : D
     5  : p
     6  : q
     7  : Mx
     8  : McBGB
     9  : McBAGB
     10 : McDU
     
     Returns
     --------
     t : double
     Age in Myrs
     
     
     */
    
    // Get parameters
    double AHe = GBParams[2];
    double B   = GBParams[3];
    double D   = GBParams[4];
    double p   = GBParams[5];
    double q   = GBParams[6];
    double Mx  = GBParams[7];
    
    // timescales
    double tHeMS = timescales[15];
    double LTHe  = luminosityEndHeliumMainSequence(Mass);
    double Lx    = calculateLuminosityLx(GBParams);
    
    double tinf1 = tHeMS + ((1.0)/((p-1)*AHe*D)) * pow((D/LTHe), ((p-1)/p));
    double tx    = tinf1 - (tinf1 - tHeMS) * pow((LTHe/Lx), ((p-1)/p));
    double tinf2 = tx + ((1.0)/((q-1)*AHe*B)) * pow((B/Lx), ((q-1)/q));
    
    // Oh dear, you really are stupid sometimes Mx = Mc(t=tx)
    // Should really have a function to do these since you use them so much . . .
    //double Mx1 = pow(((p-1)*AHe*D*(tinf1 - tx)), (1.0/(1-p)));
    //double Mx2 = pow(((q-1)*AHe*B*(tinf2 - tx)), (1.0/(1-q)));
    //double Mx = pow((B/D), (1.0/(p-q))); // Equation 38 should be in GBparams
    
    // DEBUGGING
    //std::cout << "Check these are the same : " << Mx << " " << GBParams[7] << std::endl; // they are!
    
    // Return appropriate age
    if(coreMass <= Mx){
        //std::cout << "check this Mc > Mx : " << coreMass << " " << Mx << std::endl;
        return tinf1 - pow(coreMass, (1-p))/((p-1)*AHe*D);
    }
    else{
        //std::cout << "check this Mc < Mx : " << coreMass << " " << Mx << std::endl;
        return tinf2 - pow(coreMass, (1-q))/((q-1)*AHe*B);
    }
    
    // still confused as to why t <= tx translates to Mc >= Mx ?
    
}

// Need functions for evolution of core mass, radius, luminosity of helium giant branch
double coreMassHeliumGiantBranch(double GBParams[nGBParams], double timescales[nTimescales], double Mass, double time){
    /*
     Calculate the core mass on the Helium Giant branch
     
     Given by Equation 39 modified as described after Equation 84
     
     Parameters
     -----------
     GBParams : double array
     Array containing GB parameters
     timescales : double array
     Array containing timescales
     Mass : double
     Mass in Msol
     time : double
     Time in Myrs
     
     GBParams
     ---------
     0  : A_H
     1  : A_H,He
     2  : A_He
     3  : B
     4  : D
     5  : p
     6  : q
     7  : Mx
     8  : McBGB
     9  : McBAGB
     10 : McDU
     
     timescales
     -----------
     
     
     
     Returns
     --------
     Mc : double
     Core mass in Msol
     
     
     */
    
    // Get parameters
    double AHe = GBParams[2];
    double B   = GBParams[3];
    double D   = GBParams[4];
    double p   = GBParams[5];
    double q   = GBParams[6];
    
    // timescales
    double tHeMS = timescales[15];
    
    // Required luminosities
    double Lx    = calculateLuminosityLx(GBParams); // should just be able to use L < Lx rather than t < tx
    double LTHe  = luminosityEndHeliumMainSequence(Mass);
    
    // Where do these come from -- check -- should be computed by GBParams and passed to this function FINISH ME
    double tinf1 = tHeMS + (1.0/((p-1.0)*AHe*D)) * pow((D/LTHe), ((p-1.0)/p));
    double tx    = tinf1 - (tinf1 - tHeMS) * pow((LTHe / Lx), ((p-1.0)/p));
    double tinf2 = tx + (1.0/((q-1.0)*AHe*B)) * pow((B/Lx), ((q-1.0)/q));
    
    //std::cout << "time = " << time << std::endl;
    //std::cout << "tx, tinf1, tinf2 = " << tx << " " << tinf1 << " " << tinf2 << std::endl;
    
    if(time <= tx){
        double Mc = pow(((p-1.0)*AHe*D*(tinf1 - time)), (1.0/(1.0-p)));
        //std::cout << "1) McHeGB = " << Mc << std::endl;
        return Mc;
    }
    else{
        double Mc = pow(((q-1.0)*AHe*B*(tinf2 - time)), (1.0/(1.0-q)));
        //std::cout << "2) McHeGB = " << Mc << std::endl;
        return Mc;
    }
    
}

double luminosityHeliumGiantBranch(double GBParams[nGBParams], double Mass, double coreMass){
    /*
     Calculate the giant branch luminosity for a helium star
     
     Given by Equation 84 of Hurley et al 2000
     
     Parameters
     -----------
     GBParams : double array
     Array containing GB parameters
     Mass : double
     Mass in Msol
     coreMass : double
     Core Mass in Msol
     
     Returns
     --------
     LHeGB : double
     post-HeMS luminosity in Lsol
     
     
     */
    double B = GBParams[3]; // Should be 4.1E4
    double D = GBParams[4]; // Should be (5.5E4)/(1.0 + 0.4 * M^4), see just after Equation 84
    return std::min((B * pow(coreMass, 3.0)), (D * pow(coreMass, 5.0)));
}

double radiusHeliumGiantBranch(double Mass, double Luminosity, int &STELLAR_TYPE){
    /*
     Calculate the giant branch radius for a helium star
     
     Given by Equations 85, 86, 87, 88 in Hurley et al 2000
     
     Parameters
     ------------
     Mass : double
     Mass in Msol
     Luminosity : double
     Luminosity in Lsol
     
     Returns
     --------
     RHeGB : double
     Radius on the helium giant branch / post-HeMs
     
     
     */
    
    double RZHe = radiusHeliumMainSequenceZAMS(Mass);
    
    double LTHe = luminosityEndHeliumMainSequence(Mass);
    
    double lamda = 500.0 * (2.0 + pow(Mass, 5.0))/(pow(Mass, 2.5));
    
    // If R1 gives the radius, the star is on the naked helium HG
    double R1 = RZHe * pow((Luminosity/LTHe), 0.2) + (0.02 * (exp(Luminosity/lamda) - exp(LTHe/lamda)));
    
    // If R2 gives the radius, the star is on helium GB
    double R2 = 0.08 * pow(Luminosity, 0.75);                   // Mimics the Hayashi track
    
    bool debugging = false;
//    debugging = true;
    
    // May need to change this according to section 2.1.2 of http://iopscience.iop.org/article/10.1086/340304/pdf which talks about updated helium star evolution
    // or replace with helium star tracks from binary_c
    // Especially important for low mass helium stars, BNS progenitors
    // This paper suggest mass below which envelope is convective is Mconv = 4.5 Msol, they leave it as an uncertain model parameterﰗ-- here it is set in constants (as usual, may want to edit from commandLine)
    // Rapid expansion given by the second term in R1

    if(R1 < R2){
        
        // DEBUGGING - seems fine
        if(debugging){
            std::cout << "1) HeHG, R1 = " << R1 << std::endl;
            std::cout << "M = " << Mass << std::endl;
        }
        
        STELLAR_TYPE = NAKED_HELIUM_STAR_HETZSPRUNG_GAP;
        
        return R1;
        
    }
    else{
        
        // DEBUGGING - seems fine
        if(debugging){
            std::cout << "2) HeGB, R2 = " << R2 << std::endl;
            std::cout << "M = " << Mass << std::endl;
        }
        
        STELLAR_TYPE = NAKED_HELIUM_STAR_GIANT_BRANCH; // Has a deep convective envelope
        return R2;
    }
    
}

// Evolution of L, Mc with time on Helium Giant Branch? see just after equation 84

double maximumCoreMass(double Mass){
    /*
     Calculate the maximum core mass
     
     Given by Equation 89 in Hurley et al 2000
     
     Parameters
     -----------
     Mass : double
     Mass in Msol
     
     Returns
     --------
     McMax : double
     Maximum core mass in Msol
     
     
     */
    return std::min((1.45 * Mass - 0.31), (Mass));
}

double maximumCoreMassSN(double GBParams[nGBParams]){
    /*
     Calculate the core mass at which the AGB phase is terminated in a SN/loss of envelope
     
     Given by Equation 75 in Hurley et al 2000
     
     Parameters
     -----------
     GBParams : double array
     Array containing GB parameters
     
     Returns
     --------
     McSN : double
     Maximum core mass before supernova
     
     
     */
    double McBAGB = GBParams[9];
    double McSN   = std::max((Mch), (0.773 * McBAGB - 0.35));
    return McSN;
}

// Metallicity dependent, can precompute? put into mass cutoffs.
double calculateMup(double bn_coefficients[nrowsb]){
    /*
     Metallicity dependent value Mup, the mass below (?) which a SN will result in a massless remnant
     
     Given by inverting Equation 66 in Hurley et al for McBAGB = 1.6
     For example values, see Table 1 in Pols et al 1998
     
     Parameters
     -----------
     Metallicity : double
     Metallicity Z (Z = 0.02 = Zsol)
     
     Returns
     --------
     Mup : double
     Boundary mass for massless remnant after SN
     
     
     */
    
    double MCBAGB = 1.6;
    
    return pow(((pow(MCBAGB, 4.0) - bn_coefficients[38 - 1])/bn_coefficients[36 - 1]), (1.0/bn_coefficients[37 - 1]));
    
}

double calculateMec(double bn_coefficients[nrowsb]){
    /*
     Metallicity dependent value Mec, the mass below which a SN will be an electron capture (EC) supernova, resulting in a NS
     
     Given by inverting Equation 66 in Hurley et al for MCBAGB = 1.6
     
     For example values see Table 1 in Pols et al 1998
     
     Parameters
     ----------
     bn_coefficients : double array
     Array containing metallicity dependent bn coefficients
     
     Returns
     -------
     Mec : double
     Boundary mass for electron capture supernovae
     
     
     */
    
    double MCBAGB = 2.25;
    
    return pow(((pow(MCBAGB, 4.0) - bn_coefficients[38 - 1])/bn_coefficients[36 - 1]), (1.0/bn_coefficients[37 - 1]));
    
}

double calculateMetalicityX(double Metallicity){
    /*
     Calculate Hydrogen fraction X
     
     Parameters
     -----------
     Metallicity : double
     Metallicity Z (Z = 0.02 = Zsol)
     
     Returns
     ---------
     X : double
     Hydrogen fraction X (X = 0.70 = Xsol)
     
     
     */
    return 0.76 - 3.0 * Metallicity;
}

double calcluateMetallicityY(double Metallicity){
    /*
     Parameters
     ------------
     Metallicity : double
     Metallicity Z (Z = 0.02 = Zsol)
     Returns
     --------
     Y : double
     Helium fraction Y (Y = 0.28 = Ysol)
     
     
     */
    return 0.24 + 2.0 * Metallicity;
}

double perturbation_mu(double Mass, double coreMass, double Luminosity, int STELLAR_TYPE){
    /*
     Calculate the perturbation parameter mu
     
     Given by Equations 97 and 98 in Hurley et al 2000
     
     Parameters
     -----------
     Mass : double
     Mass in Msol
     coreMass : double
     coreMass in Msol
     Luminosity : double
     Luminosity in Lsol
     STELLAR_TYPE : int
     Stellar type
     
     Returns
     --------
     mu : double
     perturbation mu
     
     
     */
    // Main sequence without core envelope separation
    if(STELLAR_TYPE <= MS_MORE_THAN_07){
        return 5.0;     // mu > 1.0
    }
    // Any star with well-defined core/envelope structure (i.e. post-MS)
    else if(STELLAR_TYPE > MS_MORE_THAN_07 and STELLAR_TYPE <= THERMALLY_PULSING_ASYMPTOTIC_GIANT_BRANCH){
        double kappa = -0.5;
        double L0    = 7.0E4;
        double firstTermTop = Mass - coreMass;
        double firstTermBottom = Mass;
        double firstTerm = firstTermTop / firstTermBottom;
        double secondTerm = std::min(5.0, std::max(1.2, pow((Luminosity/L0), kappa)));
//        std::cout << "1) Mu: " << firstTerm*secondTerm << std::endl;
        return firstTerm * secondTerm;
    }
    // Helium Main Sequence
    else if(STELLAR_TYPE == NAKED_HELIUM_STAR_MS){
//        std::cout << "2) Mu: " << 5.0 << std::endl;
        return 5.0; // mu > 1.0
    }
    // Helium Giant Branch
    else if(STELLAR_TYPE == NAKED_HELIUM_STAR_GIANT_BRANCH or STELLAR_TYPE == NAKED_HELIUM_STAR_HETZSPRUNG_GAP){
        double McMax = maximumCoreMass(Mass);
        return 5.0 * ((McMax - coreMass)/McMax);
    }
    else{
        std::cout << "Error in perturbationMu" << std::endl;
        return 5.0;     // mu > 1.0
    }
}

double perturbation_b(double Mass){
    /*
     Calculate the perturbation parameter b
     
     Given by Equation 103 in Hurley et al 2000
     
     Parameters
     -----------
     Mass : double
     Mass in Msol
     
     Returns
     ---------
     b : double
     Perturbation parameter b
     
     
     */
    return 0.002 * std::max(1.0, (2.5 / Mass));
}

double perturbation_s(double mu, double Mass){
    /*
     Calculate the perturbation parameter s
     
     Given by Equation 101 in Hurley et al 2000
     
     Parameters
     -----------
     mu : double
     Decides when to perturb luminosity/radius due to small envelope
     Mass : double
     Mass in Msol
     
     
     Returns
     --------
     perturbation_s : double
     Perturbation parameter s
     
     
     */
    
    double b = perturbation_b(Mass);
    
    double mu_over_b = mu / b;
    
    double mu_over_b_cubed = mu_over_b * mu_over_b * mu_over_b;
    
    double top = (1.0 + (b * b * b)) * mu_over_b_cubed;
    double bottom = 1.0 + mu_over_b_cubed;
    
    double s = top / bottom;
    
    // Debugging
    //std::cout << "b, s = " << b << " " << s << std::endl;
    
    return s;
}

double perturbation_c(double Mass){
    /*
     Calculate the perturbation parameter c
     
     Given by Equation 104 in Hurley et al 2000
     
     Parameters
     -----------
     Mass : double
     Mass in Msol
     
     Returns
     --------
     c : double
     Perturbation parameter c
     
     
     */
    return 0.006 * std::max(1.0, (2.5 / Mass));
}

double perturbation_q(double Radius, double Rc){
    /*
     Calculate the perturbation parameter q
     
     Given by Equation 105 in Hurley et al 2000
     
     Parameters
     -----------
     Radius : double
     Radius in Rsol
     Rc : double
     Radius that the remnant would have if the star immediately lost its envelope (in Rsol)
     
     Returns
     --------
     q : double
     Perturbation parameter q
     
     
     */
    return log(Radius / Rc); // really is natural log
}

double perturbation_r(double mu, double Mass, double Radius, double Rc){
    /*
     Calculate the perturbation parameter r
     
     Given by Equation 102 in Hurley et al 2000
     
     Parameters
     -----------
     mu : double
     Determines whether the envelope is small
     Mass : double
     Mass in Msol
     Radius : double
     Radius in Rsol
     Rc : double
     Radius remnant would have if star immediately lost its envelope (in Rsol)
     
     Returns
     --------
     s : double
     Perturbation parameter
     
     
     */
    
    double r = 0.0;
    
    if(mu < 0.0){
        r = 0.0;
        
        //std::cout << "r = " << r << std::endl;
    }
    else{
        double c = perturbation_c(Mass);
        double q = perturbation_q(Radius, Rc);
        
        double mu_over_c = mu / c;
        double mu_over_c_cubed = mu_over_c * mu_over_c * mu_over_c;
        
        double power = 0.1 / q;
        double powerMax = -14.0 / log10(mu);
        
        power = std::min(power, powerMax);
        
        double top = (1.0 + (c * c * c))*mu_over_c_cubed*pow((mu), (power));
        double bottom = (1.0 + mu_over_c_cubed);
        
        r = top / bottom;
        
        //std::cout << "c, r = " << c << " " << r << std::endl;
    }
    
    return r;
}

void perturbation_Lc_Rc(double &Lc, double &Rc, double an_coefficients[nrowsa], double bn_coefficients[nrowsb], double GBParams[nGBParams], double massCutoffs[nMassCutoffs], double timescales[nTimescales], double time, double Mass, double coreMass, double HeCoreMass, double COCoreMass, double Metallicity, double logMetallicityXi, int STELLAR_TYPE, unsigned long m_randomSeed){
    /*
     Calculate the luminosity and radius of the remnant the star would become if it lost all of its envelope immediately.
     
     i.e. M = Mc, coreMass
     
     Given after Equation 105 in Hurley et al 2000
     
     Parameters
     -----------
     &Lc : double
     Luminosity that remnant would have if star lost all of its envelope
     &Rc : double
     Radius that remnant would have if star lost all of its envelope
     an_coefficients : double array
     Array containing metallicity dependent an coefficients
     bn_coefficients : double array
     Array containing metallicity dependent bn coefficients
     GBParams : double array
     Array containing GB parameters
     massCutoffs : double array
     Array of metallicity dependent mass cutoffs
     timescales : double
     Array containing timescales
     time : double
     Time in Myrs
     Mass : double
     Mass in Msol
     coreMass : double
     Core mass in Msol
     HeCoreMass : double
     He Core mass in Msol
     COCoreMass : double
     CO Core mass in Msol
     Metallicity : double
     Metallicity Z (Z = 0.02 = Zsol)
     logMetallicityXi : double
     logMetallicityXi
     STELLAR_TYPE : int
     Stellar type
     
     timescales
     -----------
     List of timescales (all in units of Myrs)
     0 : tMS  - Main sequence
     1 : tBGB - Base Giant Branch
     2 : tHeI - Helium ignition
     3 : tHe  - Helium burning
     #-- First Giant Branch (FGB)
     4 : Giant tinf1
     5 : Giant tinf2
     6 : Giant t(Mx)
     #-- Early Asymptotic Giant Branch (EAGB)
     7 : FAGB tinf1
     8 : FAGB tinf2
     9 : FAGB t(Mx)
     #-- Thermally Pulsating Asymptotic Giant Branch (TPAGB)
     10 : SAGB tinf1
     11 : SAGB tinf2
     12 : SAGB t(Mx)
     13 : TP (tDU?)
     14 : T(McMax)
     #-- not sure what TP, t(McMax) are
     #-- what about tHeMS etc?
     15 : tHeMS
     
     Returns
     --------
     Lc : double
     Luminosity of core in Lsol
     Rc : double
     Radius of core in Rsol
     
     
     */
    
    double MHeF = massCutoffs[1];
    
    double LZHe = luminosityHeliumMainSequenceZAMS(coreMass);
    double RZHe = radiusHeliumMainSequenceZAMS(coreMass);
    
    int currentStellarType = STELLAR_TYPE;

    bool debugging = false;
    //debugging = true;
    
    if(STELLAR_TYPE == HERTZSPRUNG_GAP or STELLAR_TYPE == FIRST_GIANT_BRANCH){
        if(Mass < MHeF){
            Lc = LZHe;
            Rc = RZHe;
        }
        else{
            Lc = luminosityWhiteDwarf(0, coreMass, Metallicity, HELIUM_WHITE_DWARF);
            Rc = radiusWhiteDwarf(coreMass);
        }

        if(debugging){

            std::cout << "Lc, Rc = " << Lc << " " << Rc << std::endl;

        }
    }
    else if(STELLAR_TYPE == CORE_HELIUM_BURNING){
        // Remnant will be an evolved helium MS star with M = Mc  and tau = (t - tHeI)/tHe so use appropriate functions:
        // Will these two be the same? can I just pass the timescales array?
        double tHeI = timescales[2];
        double tHe  = timescales[3];
        double tau  = (time - tHeI)/tHe; // Should be 0 <= tau <= 1
        
        if(tau < 0.0){
            tau = 0.0;
        }

        if(tau > 1.0){
            tau = 1.0;
        }

        Lc = luminosityHeliumMainSequence(tau, coreMass);
        Rc = radiusHeliumMainSequence(tau, coreMass);
        
        if(debugging){

            std::cout << "Lc, Rc = " << Lc << " " << Rc << std::endl;

        }

    }
    else if(STELLAR_TYPE == EARLY_ASYMPTOTIC_GIANT_BRANCH){
        
        // Remnant will be a helium HG or GB star with M = McHe, Mc = McCO
        
        Lc = luminosityHeliumGiantBranch(GBParams, HeCoreMass, COCoreMass);  // Comes from HeGB MC-L relation with Mc = McCO
        Rc = radiusHeliumGiantBranch(HeCoreMass, Lc, currentStellarType);    // Comes from HeGB Radius (McHe, Lc)
        
    }
    else if(STELLAR_TYPE == THERMALLY_PULSING_ASYMPTOTIC_GIANT_BRANCH or STELLAR_TYPE == NAKED_HELIUM_STAR_HETZSPRUNG_GAP or STELLAR_TYPE == NAKED_HELIUM_STAR_GIANT_BRANCH){
        
        Lc = luminosityWhiteDwarf(0.0, coreMass, Metallicity, CARBON_OXYGEN_WHITE_DWARF);
        Rc = radiusWhiteDwarf(coreMass);
        
    }
    else{
        std::cerr << "Error in perturbation_Lc_Rc. Unsupported stellar type in perturbation_Lc_Rc" << std::endl;
        Lc = 0.0;
        Rc = 0.0;
    }
}

void perturbedLuminosityRadius(double &Luminosity, double &Radius, double mu, double an_coefficients[nrowsa], double bn_coefficients[nrowsb], double GBParams[nGBParams], double massCutoffs[nMassCutoffs], double timescales[nTimescales], double time, double Mass, double coreMass, double HeCoreMass, double COCoreMass, double Metallicity, double logMetallicityXi, int STELLAR_TYPE, unsigned long m_randomSeed){
    /*
     Calculate the perturbed luminosity and radius
     
     Given by Equation 99 and 100 in Hurley et al 2000
     
     Parameters
     -----------
     &Luminosity : double
        Unperturbed Luminosity in Lsol to perturb if mu < 1
     &Radius : double
        Unperturbed Radius in Rsol to perturb if mu < 1
     mu : double
        Parameter to determine whether to peturb luminosity/radius for small envelopes
     an_coefficients : double array
        Array containing metallicity dependent an coefficients
     bn_coefficients : double array
        Array containing metallicity dependent bn coefficients
     GBParams : double array
        Array containing GB parameters
     massCutoffs : double array
        Array containing metallicity dependent mass cutoffs
     timescales : double array
        Array containing timescales
     time : double
        Age in Myrs
     Mass : double
        Mass in Msol
     coreMass : double
        Core Mass Mc in Msol
     HeCoreMass : double
        HeCoreMass in Msol
     COCoreMass : double
        COCoreMass in Msol
     Metallicity : double
        Metallicity Z (Z = 0.02 = Zsol)
     logMetallicityXi : double
        logMetallicityXi
     STELLAR_TYPE : int
        Stellar type
     
     Returns
     --------
     LPrime : double
        Perturbed luminosity
     RPrime : double
        Perturbed Radius
     
     
     */
    
    bool debugging = false;
    //debugging = true;

    if(debugging){

        std::cout << "Perturbed luminosity and radius function" << std::endl;
        std::cout << "Unperturbed L, R: " << Luminosity << " " << Radius << std::endl;

    }

    if(STELLAR_TYPE == MS_LESS_THAN_07 or STELLAR_TYPE == MS_MORE_THAN_07 or STELLAR_TYPE == NAKED_HELIUM_STAR_MS){
        
        Luminosity = Luminosity;
        Radius = Radius;
        
    }
    else if (STELLAR_TYPE > MS_MORE_THAN_07 and STELLAR_TYPE <= THERMALLY_PULSING_ASYMPTOTIC_GIANT_BRANCH){
        
        // If mu < 1.0 then do this.

        if(debugging){

            std::cout << "perturbation mu = " << mu << std::endl;

        }

        if(mu < 1.0){
            
            double Lc = 0.0;
            double Rc = 0.0;
            
            perturbation_Lc_Rc(Lc, Rc, an_coefficients, bn_coefficients, GBParams, massCutoffs, timescales, time, Mass, coreMass, HeCoreMass, COCoreMass, Metallicity, logMetallicityXi, STELLAR_TYPE, m_randomSeed);
            
            //double s = perturbation_s(Mass, HeCoreMass, Luminosity, STELLAR_TYPE);
            //double r = perturbation_r(Mass, HeCoreMass, Luminosity, Radius, Rc, Metallicity, STELLAR_TYPE);
            double s = perturbation_s(mu, Mass);
            double r = perturbation_r(mu, Mass, Radius, Rc);

            Luminosity = Lc * pow((Luminosity / Lc), s);
            Radius     = Rc * pow((Radius / Rc), r);

            if(debugging){

                std::cout << "perturbed L, R: " << Luminosity << " " << Radius << std::endl;

                // DEBUGGING -- Why does perturbation give Rc = 0, RWD = 0? should mu not be < 1?
                std::cout << "Mc = " << coreMass << std::endl;
                // std::cout << "LC, RC = " << Lc << " " << Rc << std::endl;
                //std::cout << "log10(Lc), log10(Rc) = " << log10(Lc) << " " << log10(Rc) << std::endl;

                // big difference between mu = 0.01 and mu = 0.001 which means we really have to sort out the mass loss
                //std::cout << "if mu = 0.001, r = " << perturbation_r(1E-3, Mass, Radius, Rc) << std::endl;

                // std::cout << "RWD = " << radiusWhiteDwarf(coreMass) << std::endl;
                // std::cout << "R_ZHe = " << radiusHeliumMainSequenceZAMS(coreMass) << std::endl;
                std::cout << "s, r = " << s << " " << r << std::endl;

            }
            
        }
        
    }
    else if(STELLAR_TYPE == HELIUM_WHITE_DWARF or STELLAR_TYPE == CARBON_OXYGEN_WHITE_DWARF or STELLAR_TYPE == OXYGEN_NEON_WHITE_DWARF){
        
        // L, R FINISH ME
        
        Luminosity = Luminosity;
        Radius = Radius;
        
    }
    else if(STELLAR_TYPE == NAKED_HELIUM_STAR_GIANT_BRANCH or STELLAR_TYPE == NAKED_HELIUM_STAR_HETZSPRUNG_GAP){
        
        if(mu < 1.0){
            
            double Lc = 0.0;
            double Rc = 0.0;
            
            perturbation_Lc_Rc(Lc, Rc, an_coefficients, bn_coefficients, GBParams, massCutoffs, timescales, time, Mass, coreMass, HeCoreMass, COCoreMass, Metallicity, logMetallicityXi, STELLAR_TYPE, m_randomSeed);
            
            //double s = perturbation_s(Mass, HeCoreMass, Luminosity, STELLAR_TYPE);
            //double r = perturbation_r(Mass, HeCoreMass, Luminosity, Radius, Rc, Metallicity, STELLAR_TYPE);
            double s = perturbation_s(mu, Mass);
            double r = perturbation_r(mu, Mass, Radius, Rc);
            
            //std::cout << "&Luminosity = " << Luminosity << std::endl;
            //std::cout << "&Radius = " << Radius << std::endl;
            
            // DEBUGGING -- Why does perturbation give Rc = 0, RWD = 0? should mu not be < 1?
            //std::cout << "LC, RC = " << Lc << " " << Rc << std::endl;
            //std::cout << "log10(Lc), log10(Rc) = " << log10(Lc) << " " << log10(Rc) << std::endl;
            //std::cout << "Mc = " << coreMass << std::endl;
            //std::cout << "RWD = " << radiusWhiteDwarf(coreMass) << std::endl;
            //std::cout << "s, r = " << s << " " << r << std::endl;
            
            Luminosity = Lc * pow((Luminosity / Lc), s);
            Radius     = Rc * pow((Radius / Rc), r);
            
        }
        
    }
    else{
        
        std::cout << "Invalid Stellar type for small envelope perturbations" << std::endl;
        
    }
    
}

double initialEnvelopeMass(double Mass){
    /*
     Calculate the initial convective envelope mass M_env0
     
     Given just after Equation 111 in Hurley et al 2000
     
     Parameters
     -----------
     Mass : double
     Mass in Msol
     
     Returns
     --------
     Menv0 : double
     ZAMS envelope mass
     
     
     */
    
    double envMass = 0.0;
    
    // Initial envelope mass (ZAMS)
    if(Mass < 0.35){ // Stars fully convective so Menv = M
        envMass = Mass;
    }
    else if(Mass >= 0.35 and Mass <= 1.25){
        envMass = 0.35 * pow(((1.25 - Mass)/0.9), (2.0)); // Given after Equation 111
    }
    else{
        envMass = 0.0; // Really 0
    }
    
    return envMass;
}

double envelopeMass(double Mass, double coreMass, double tau, int STELLAR_TYPE){
    /*
     Calculate envelope mass as a function of time
     
     Given after Equation 111 in Hurley et al 2000
     
     Parameters
     -----------
     Mass : double
     Mass in Msol
     coreMass : double
     Core mass in Msol
     tau : double
     Relative lifetime
     STELLAR_TYPE : int
     Stellar type
     
     
     */
    
    double Menv = 0;
    double Menv0 = initialEnvelopeMass(Mass);
    
    if(STELLAR_TYPE == MS_LESS_THAN_07 or STELLAR_TYPE == MS_MORE_THAN_07){
        
        Menv = Menv0 * pow((1.0 - tau), (0.25));
        
    }
    else if(STELLAR_TYPE == HERTZSPRUNG_GAP){
        
        Menv = tau * (Mass - coreMass);
        
    }
    else{ // For most phases with core envelope separation
        Menv = Mass - coreMass;
    }
    
    return Menv;
}

double Star::chooseTimestep(int STELLAR_TYPE, double time, double timescales[nTimescales],const programOptions &options ){
    /*
     Choose timestep for evolution based on current stellar type
     
     Can obviously do this your own way
     
     Given in the discussion in Hurley et al 2000
     
     Parameters
     -----------
     STELLAR_TYPE : int
     Current stellar type
     time : double
     Current age of star in Myrs
     timescales : array
     Timescales in MYrs
     
     List of timescales (all in units of Myrs)
     0 : tMS  - Main sequence
     1 : tBGB - Base Giant Branch
     2 : tHeI - Helium ignition
     3 : tHe  - Helium burning
     #-- First Giant Branch (FGB)
     4 : Giant tinf1
     5 : Giant tinf2
     6 : Giant t(Mx)
     #-- Early Asymptotic Giant Branch (EAGB)
     7 : FAGB tinf1
     8 : FAGB tinf2
     9 : FAGB t(Mx)
     #-- Thermally Pulsating Asymptotic Giant Branch (TPAGB)
     10 : SAGB tinf1
     11 : SAGB tinf2
     12 : SAGB t(Mx)
     13 : TP (tDU?)
     14 : T(McMax)
     #-- not sure what TP, t(McMax) are
     #-- what about tHeMS etc?
     15 : tHeMS
     #-- HeGB
     16 : HeGB tinf1
     17 : HeGB tinf2
     18 : HeGB tx
     
     Returns
     --------
     dt : double
     Suggested timestep
     
     
     */
    
    double dtk = 0.0;
    double dte = 0.0;
    double timestep = 0.0;
    double timeCap = 1E-4; // 100 years
    bool    debugging = false;
//    debugging = true;
    
    if(STELLAR_TYPE == MS_LESS_THAN_07 or STELLAR_TYPE == MS_MORE_THAN_07){
        dtk = 1E-2 * timescales[0];
        dte = timescales[0] - time; // tMS - t
        //std::cout << "dtk, dte = " << dtk << " " << dte << std::endl;
        // This timestep is not short enough to resolve the hook at the end of the main sequence for HM stars so choose shorter
        if(dte < dtk){
            dtk = 1E-3 * timescales[0];
        }
    }
    else if(STELLAR_TYPE == HERTZSPRUNG_GAP){
        dtk = 0.05 * (timescales[1] - timescales[0]);
        dte = timescales[1] - time;         // tBGB - t
		
		if(debugging){
			std::cout << "dtk: " << dtk << std::endl;
			std::cout << "dte: " << dte << std::endl;
		}
    }
    else if(STELLAR_TYPE == FIRST_GIANT_BRANCH){
        if(time <= timescales[6]){ // ah because timescales[4,5,6] are not calculated yet
            dtk = 0.02 * (timescales[4] - time); // tinf1 - t
        }
        else{
            dtk = 0.02 * (timescales[5] - time); // tinf2 - t
        }
        dte = timescales[2] - time; // tHeI - t
    }
    else if(STELLAR_TYPE == CORE_HELIUM_BURNING){
        // Messing around with step size
        //dtk = 0.02 * timescales[3]; // tHe
        //std::cout << "suggested dtk = " << dtk << std::endl;
        dtk = 2E-3 * timescales[3];
        dte = timescales[2] + timescales[3] - time; // tHeI + tHe - t
        //std::cout << "current dtk, dte = " << dtk << " " << dte << std::endl;
        if(dtk < dte){
            dtk = 2E-3 * timescales[3];
        }
        // Issues with timestep in CHeB?
    }
    else if(STELLAR_TYPE == EARLY_ASYMPTOTIC_GIANT_BRANCH){
        if(time <= timescales[9]){ // t <= tx
            dtk = 0.02 * (timescales[7] - time); // tinf1 - t
        }
        else{
            dtk = 0.02 * (timescales[8] - time); // tinf2 - t
        }
        //dte = timescales[13] - time; // tDU - t? // huh but once you get to tDU, dt = 0? I guess this has to be based on core mass
        dte = dtk; // FINISH ME
    }
    else if(STELLAR_TYPE == THERMALLY_PULSING_ASYMPTOTIC_GIANT_BRANCH){
        if(time <= timescales[12]){
            dtk = 0.02 * (timescales[10] - time); // tinf1 - t
        }
        else{
            dtk = 0.02 * (timescales[11] - time); // tinf2 - t
        }
        dte = 5E-3;
    }
    else if(STELLAR_TYPE == NAKED_HELIUM_STAR_MS){
        dtk = 0.05 * timescales[15]; // tHeMS
        dte = timescales[15] - time; // tHeMS - time
    }
    else if(STELLAR_TYPE == NAKED_HELIUM_STAR_HETZSPRUNG_GAP or STELLAR_TYPE == NAKED_HELIUM_STAR_GIANT_BRANCH){
        // Having problems here. Very crowded points.
        if(debugging){
            std::cout << "time, 'tx': " << time << " " << timescales[18] << std::endl;
        }
        if(time <= timescales[18]){ // if t <= tx
            if(debugging){
                std::cout << "If 1) timescales[16]: " << timescales[16] << std::endl;
                std::cout << "timescales[16] - time: " << timescales[16] - time << std::endl;
            }
            dtk = 0.02 * (timescales[16] - time); // tinf1 - t
        }
        else{
            if(debugging){
                std::cout << "If 2) timescales[17]: " << timescales[17] << std::endl;
                std::cout << "timescales[17] - time: " << timescales[17] - time << std::endl;
            }
            dtk = 0.02 * (timescales[17] - time); // tinf2 - t
        }
        if(debugging){
            std::cout << "dtk: " << dtk << std::endl;
        }
        dte = dtk; // how do you define time till end of HeGB? Depends on time till McCO = McSN right?
    }
    else if(STELLAR_TYPE >= HELIUM_WHITE_DWARF){
        dtk = std::max(0.1, 10.0 * time);
        dtk = std::max(0.1, 10.0 * dtk);
        dtk = std::min(dtk, 5E2);
        dte = dtk;
    }
    else{
        if(debugging){
            std::cout << "Unsupported stellar type" << std::endl;
        }
        dtk = std::max(0.1, 10.0 * time);
        dtk = std::max(0.1, 10.0 * dtk);
        dtk = std::min(dtk, 5E2);
        dte = dtk;
    }
    
    if(std::min(dtk, dte) < 0){
        if(debugging)
        std::cerr << m_randomSeed << "\tError - Negative timestep dt = " << std::min(dtk, dte) << std::endl;
        
    }



    
    timestep=std::max(std::min(dtk, dte), timeCap);
    
    if(debugging){
        std::cout << "timestep: " << timestep << std::endl;
        if(std::min(dtk, dte) <= timeCap){
            std::cerr << m_randomSeed << "\tUse 100 years as minimum timestep for nuclear evolution. " << std::endl;
        }
    }

    return timestep;
    
}
void calculateGiantBranchParameters(double GBParams[nGBParams], double an_coefficients[nrowsa], double bn_coefficients[nrowsb], double massCutoffs[nMassCutoffs], double Mass, double logMetallicityXi, int STELLAR_TYPE, unsigned long m_randomSeed){
    /*
     Set the Giant Branch (GB) parameters required by several functions
     
     Parameters
     -----------
     GBParams : double array
        Array containing the GB parameters which you want to set
     an_coefficients : double array
        Array containing metallicity dependent an coefficients
     bn_coefficients : double array
        Array containing metallicity dependent bn coefficients
     massCutoffs : double array
        Array containing metallicity dependent mass cutoffs
     logMetallicityXi : double
        logMetallicityXi
     STELLAR_TYPE : int
        Stellar evolution type (MS, WD, NS, BH etc)
     
     Returns
     --------
     updates GBParams
     
     List of GB Parameters:
     0  : A_H
     1  : A_H,He
     2  : A_He
     3  : B
     4  : D
     5  : p
     6  : q
     7  : Mx
     8  : McBGB
     9  : McBAGB
     10 : McDU
     11 : McCOBAGB
     
     
     */
    
	if(STELLAR_TYPE == MS_LESS_THAN_07 or STELLAR_TYPE == MS_MORE_THAN_07){
		//no GBParams needed for these types
	}
	else{

		GBParams[0]  = calculateHydrogenRateConstant(Mass);                                  // A_H
		GBParams[3]  = calculateCoreMassParameter_B(Mass);                                   // B
		GBParams[4]  = calculateCoreMassParameter_D(Mass, massCutoffs, logMetallicityXi, m_randomSeed);    // D

		if(STELLAR_TYPE == NAKED_HELIUM_STAR_GIANT_BRANCH or STELLAR_TYPE == NAKED_HELIUM_STAR_HETZSPRUNG_GAP){
			GBParams[3] = 4.1E4;
			GBParams[4] = 5.5E4 / (1.0 + 0.4 * pow(Mass, 4.0));
		}

		GBParams[5]  = calculateCoreMassParameter_p(Mass, massCutoffs, m_randomSeed);                      // p
		GBParams[6]  = calculateCoreMassParameter_q(Mass, massCutoffs, m_randomSeed);                      // q


		//Parameters only needed in specific cases
		if(STELLAR_TYPE == THERMALLY_PULSING_ASYMPTOTIC_GIANT_BRANCH){
			GBParams[1]  = calculateMixedHydrogenHeliumRateConstant();                           // A_H,He
		}
		if(STELLAR_TYPE >= CORE_HELIUM_BURNING){
			GBParams[2]  = calculateHeliumRateConstant();                                        // A_He
			GBParams[10] = coreMassSecondDredgeUp(GBParams, m_randomSeed); // McDU
		}
		

		
		GBParams[7]  = calculateCoreMassParameter_Mx(GBParams);                              // Mx
		GBParams[9]  = coreMassBaseAsymptoticGiantBranch(bn_coefficients, Mass, m_randomSeed);             // McBAGB

		if(STELLAR_TYPE == NAKED_HELIUM_STAR_GIANT_BRANCH or STELLAR_TYPE == NAKED_HELIUM_STAR_HETZSPRUNG_GAP){
			GBParams[9] = Mass; // Replace McBAGB with M in Equation 75 and discussion for final stages
		}

		double MCBAGB = GBParams[9];
		GBParams[8]  = coreMassBaseGiantBranch(an_coefficients, GBParams, MCBAGB, Mass, massCutoffs, logMetallicityXi, m_randomSeed);
		
    }
}


double angularFrequencyZAMS(double M, double RZAMS){
/* Function for calculating the angular frequency of the star at ZAMS according to eqs (107) and (108) from Hurley et al. 2000
Alejandro Vigna-Gomez - 10/2015 */
	double vrot = 330*(pow(M,3.3))/(15+pow(M,3.45));	// [vrot] 	= km s-1
	return 45.35*vrot/RZAMS;                            // [omega] 	= yr -1
}


void calculateTimescales(double timescales[nTimescales], double GBParams[nGBParams], double an_coefficients[nrowsa], double bn_coefficients[nrowsb], double massCutoffs[nMassCutoffs], double Mass0, double coreMass, double COCoreMass, double Metallicity, double logMetallicityXi, int STELLAR_TYPE, unsigned long m_randomSeed){
    /*
     Calculate timescales each loop
     
     Uses M0 as described toward the end of Section 7.1
     
     Parameters
     ------------
     timescales : double array
        Array containing timescales
     GBParams : double array
        Array containing GBParams
     an_coefficients : double array
        Array containing metallicity dependent an coefficients
     bn_coefficients : double array
        Array containing metallicity dependent bn coefficients
     massCutoffs : double array
        Array containing metallicity dependent mass cutoffs
     Mass0 : double
        Effective Initial Mass in Msol
     coreMass : double
        coreMass in Msol
     COCoreMass : double
        CO core mass in Msol
     Metallicity : double
        Metallicity Z (Z = 0.02 = Zsol)
     logMetallicityXi : double
        logMetallicityXi
     STELLAR_TYPE : int
        Stellar type
     
     GBParams
     ---------
     List of GB Parameters:
     0 : A_H
     1 : A_H,He
     2 : A_He
     3 : B
     4 : D
     5 : p
     6 : q
     7 : Mx
     8 : McBGB
     
     Returns
     --------
     timescales : double array
     
     List of timescales (all in units of Myrs)
     0 : tMS  - Main sequence
     1 : tBGB - Base Giant Branch
     2 : tHeI - Helium ignition
     3 : tHe  - Helium burning
     #-- First Giant Branch (FGB)
     4 : Giant tinf1
     5 : Giant tinf2
     6 : Giant t(Mx)
     #-- Early Asymptotic Giant Branch (EAGB)
     7 : FAGB tinf1
     8 : FAGB tinf2
     9 : FAGB t(Mx)
     #-- Thermally Pulsating Asymptotic Giant Branch (TPAGB)
     10 : SAGB tinf1
     11 : SAGB tinf2
     12 : SAGB t(Mx)
     13 : TP (tDU?)
     14 : T(McMax)
     #-- not sure what TP, t(McMax) are
     #-- what about tHeMS etc?
     15 : tHeMS
     #-- Heium giant branch
     16 : HeGB tinf1
     17 : HeGB tinf2
     18 : HeGB tx
     19 : Relative duration of blue loop taubl
     20 : Relative start of blue loop taux
     21 : Relative end of blue loop tauy
     
     
     */
    
    // Start off by just calculating those you know you can -- maybe add in if statements based on STELLAR_TYPE
    // Odd order due to things relying on each other

	double AH   = GBParams[0];
    double AHHe = GBParams[1];
    double AHe  = GBParams[2];
    double B    = GBParams[3];
    double D    = GBParams[4];
    double p    = GBParams[5];
    double q    = GBParams[6];
    double McDU = GBParams[10];    
    double MFGB = massCutoffs[2];
    
	if(STELLAR_TYPE <=  MS_MORE_THAN_07){
		timescales[1] = lifetimeToBGB(an_coefficients, Mass0); // tBGB
		double tHook = lifetimeHook(an_coefficients, Mass0, timescales[1]); // tHook
		timescales[0] = lifetimeMainSequence(an_coefficients, Mass0, logMetallicityXi, timescales[1], tHook); // tMS
	}


	if((STELLAR_TYPE == HERTZSPRUNG_GAP) or (STELLAR_TYPE == FIRST_GIANT_BRANCH)){
		// Requires
		double Lx = calculateLuminosityLx(GBParams);
		double LBGB = luminosityBaseGiantBranch(an_coefficients, Mass0, m_randomSeed);
		timescales[1] = lifetimeToBGB(an_coefficients, Mass0); // tBGB
		double tHook = lifetimeHook(an_coefficients, Mass0, timescales[1]); // tHook
		timescales[0] = lifetimeMainSequence(an_coefficients, Mass0, logMetallicityXi, timescales[1], tHook); // tMS
		// double tinf1    = tBGB + ((1.0)/((p - 1.0)*AH*D)) * pow((D/LBGB), ((p-1.0)/p));
		timescales[4] = timescales[1] + ((1.0)/((p - 1.0)*AH*D)) * pow((D/LBGB), ((p-1.0)/p));          // tinf1
		
		// double tx       = tinf1 - (tinf1 - tBGB) * pow((LBGB/Lx), ((p-1.0)/p));
		timescales[6] = timescales[4] - (timescales[4] - timescales[1]) * pow((LBGB/Lx), ((p-1.0)/p));  // tx
		
		// double tinf2    = tx + ((1.0)/((q - 1.0)*AH*B)) * pow((B/Lx), ((q-1.0)/q));
		timescales[5] = timescales[6] + ((1.0)/((q - 1.0)*AH*B)) * pow((B/Lx), ((q-1.0)/q));            // tinf2
		
		// Depends on timescales[4], timescales[5] should move it after them? But timescales[4] depends on timescales[1] argh!
		timescales[2] = lifetimeHeliumIgnition(bn_coefficients, GBParams, timescales, Mass0, massCutoffs, m_randomSeed); // tHeI
		timescales[15] = lifetimeHeliumMainSequence(Mass0); //Need HeMS timescale when we strip
    }


	if(STELLAR_TYPE == CORE_HELIUM_BURNING){
		// Requires
		double Lx = calculateLuminosityLx(GBParams);
		double LBGB = luminosityBaseGiantBranch(an_coefficients, Mass0, m_randomSeed);
		timescales[1] = lifetimeToBGB(an_coefficients, Mass0); // tBGB
		double tHook = lifetimeHook(an_coefficients, Mass0, timescales[1]); // tHook
		timescales[0] = lifetimeMainSequence(an_coefficients, Mass0, logMetallicityXi, timescales[1], tHook); // tMS
		// double tinf1    = tBGB + ((1.0)/((p - 1.0)*AH*D)) * pow((D/LBGB), ((p-1.0)/p));
		timescales[4] = timescales[1] + ((1.0)/((p - 1.0)*AH*D)) * pow((D/LBGB), ((p-1.0)/p));          // tinf1
		// double tx       = tinf1 - (tinf1 - tBGB) * pow((LBGB/Lx), ((p-1.0)/p));
		timescales[6] = timescales[4] - (timescales[4] - timescales[1]) * pow((LBGB/Lx), ((p-1.0)/p));  // tx
		// double tinf2    = tx + ((1.0)/((q - 1.0)*AH*B)) * pow((B/Lx), ((q-1.0)/q));
		timescales[5] = timescales[6] + ((1.0)/((q - 1.0)*AH*B)) * pow((B/Lx), ((q-1.0)/q));            // tinf2
		// Depends on timescales[4], timescales[5] should move it after them? But timescales[4] depends on timescales[1] argh!
		timescales[2] = lifetimeHeliumIgnition(bn_coefficients, GBParams, timescales, Mass0, massCutoffs, m_randomSeed); // tHeI
		timescales[3] = lifetimeCoreHeliumBurning(an_coefficients, bn_coefficients, Mass0, coreMass, massCutoffs, m_randomSeed);
		    // Blue loop -- may need to reorganise if required by other timescales/GBParams
		timescales[19] = lifetimeBluePhase(bn_coefficients, Mass0, coreMass, massCutoffs, logMetallicityXi, m_randomSeed);            // taubl
		timescales[20] = relativeAgeStartBluePhase(bn_coefficients, Mass0, coreMass, massCutoffs, logMetallicityXi, m_randomSeed);    // taux
		timescales[21] = relativeAgeEndBluePhase(bn_coefficients, Mass0, coreMass, massCutoffs, logMetallicityXi, m_randomSeed);      // tauy

		timescales[15] = lifetimeHeliumMainSequence(Mass0); //Need HeMS timescale when we strip
	}


    if(STELLAR_TYPE >= EARLY_ASYMPTOTIC_GIANT_BRANCH){
        // Requires
        double tBAGB = lifetimeBaseAsymptoticGiantBranch(timescales);
        double LBAGB = luminosityBaseAsymptoticGiantBranch(bn_coefficients, Mass0, massCutoffs, m_randomSeed);
        double Lx    = calculateLuminosityLx(GBParams);
        // Getting these equations is described in the paragraph just above Equation 68
        // double tinf1    = tBAGB + ((1.0)/((p - 1.0) * AHe * D)) * pow((D / LBAGB), ((p - 1.0)/(p)));// Given by Equation 40
        timescales[7] = tBAGB + ((1.0)/((p - 1.0) * AHe * D)) * pow((D / LBAGB), ((p - 1.0)/(p)));
        // double tx       = tinf1 - (tinf1 - tBAGB) * pow((LBAGB / Lx), ((p - 1.0)/(p)));      // Given by Equation 41
        timescales[9] = timescales[7] - (timescales[7] - tBAGB) * pow((LBAGB / Lx), ((p - 1.0)/(p)));
        // double tinf2    = tx + ((1.0)/((q - 1.0) * AHe * B)) * pow((B / Lx), ((q - 1.0)/(q)));          // Given by Equation 42
        timescales[8] = timescales[9] + ((1.0)/((q - 1.0) * AHe * B)) * pow((B / Lx), ((q - 1.0)/(q)));
    }
    


    timescales[13] = lifetimeSecondDredgeUp(bn_coefficients, timescales, GBParams, Mass0, COCoreMass, massCutoffs, logMetallicityXi, m_randomSeed); // tDU given by Equation 70
    
    if(STELLAR_TYPE >= THERMALLY_PULSING_ASYMPTOTIC_GIANT_BRANCH){
        // Calculate LDU, Lx
        double LDU   = coreMassLuminosityRelation(GBParams, McDU, m_randomSeed);
        double Lx    = calculateLuminosityLx(GBParams);
        // Used the fortran code to check this bit
        // Requires
        double tBAGB = lifetimeBaseAsymptoticGiantBranch(timescales);
        double LBAGB = luminosityBaseAsymptoticGiantBranch(bn_coefficients, Mass0, massCutoffs, m_randomSeed);
        // Getting these equations is described in the paragraph just above Equation 68
        // double tinf1    = tBAGB + ((1.0)/((p - 1.0) * AHe * D)) * pow((D / LBAGB), ((p - 1.0)/(p)));// Given by Equation 40
        timescales[7] = tBAGB + ((1.0)/((p - 1.0) * AHe * D)) * pow((D / LBAGB), ((p - 1.0)/(p)));
        // double tx       = tinf1 - (tinf1 - tBAGB) * pow((LBAGB / Lx), ((p - 1.0)/(p)));      // Given by Equation 41
        timescales[9] = timescales[7] - (timescales[7] - tBAGB) * pow((LBAGB / Lx), ((p - 1.0)/(p)));
        // double tinf2    = tx + ((1.0)/((q - 1.0) * AHe * B)) * pow((B / Lx), ((q - 1.0)/(q)));          // Given by Equation 42
        timescales[8] = timescales[9] + ((1.0)/((q - 1.0) * AHe * B)) * pow((B / Lx), ((q - 1.0)/(q)));
        if(LDU > Lx){
            timescales[10] = timescales[7]; // tinf1
            timescales[12] = timescales[9]; // tx
            // tinf2 = tDU + ((1.0)/((q-1)*AHHe*B)) * pow((B/LDU), ((q-1)/q));
            timescales[11] = timescales[13] + ((1.0)/((q-1)*AHHe*B)) * pow((B/LDU), ((q-1)/q)); // tinf2 from Equation 72
        }
        else{ // LDU <= Lx, Equations 40-42 modified as described after Equation 70
            // tinf1 = tDU + ((1.0)/((p-1)*AHHe*D)) * pow((D/LDU), ((p-1)/p));
            timescales[10] = timescales[13] + ((1.0)/((p-1)*AHHe*D)) * pow((D/LDU), ((p-1)/p));
            // tx = tinf1 - (tinf1 - tDU) * pow((LDU/Lx), ((p-1)/p));
            timescales[12] = timescales[10] - (timescales[10] - timescales[13])*pow((LDU/Lx), ((p-1)/p));
            // tinf2 = tx + ((1.0)/((q-1)*AHHe*B)) * pow((B/Lx), ((q-1)/q));
            timescales[11] = timescales[12] + ((1.0)/((q-1)*AHHe*B)) * pow((B/Lx), ((q-1)/q));
        }
    }
    
    // t(McMax) from Equation ?
    timescales[14] = 0.0;
    //if(STELLAR_TYPE == NAKED_HELIUM_STAR_MS){
    	// tHeMS
	    
	 if(STELLAR_TYPE == NAKED_HELIUM_STAR_MS){
		timescales[15] = lifetimeHeliumMainSequence(Mass0);
		}
    //}
    // Helium Giant branch timescales
    if(STELLAR_TYPE == NAKED_HELIUM_STAR_HETZSPRUNG_GAP or STELLAR_TYPE == NAKED_HELIUM_STAR_GIANT_BRANCH){
        timescales[15] = lifetimeHeliumMainSequence(Mass0);
        // Required luminosities
        double LTHe = luminosityEndHeliumMainSequence(Mass0);
        double Lx   = calculateLuminosityLx(GBParams);
        
        // tinf1 = tHeMS + ((1.0)/((p-1.0)*AHe*D)) * (D/LTHe)**((p-1.0)/p)
        timescales[16] = timescales[15] + ((1.0)/((p-1.0)*AHe*D)) * pow((D/LTHe), ((p-1.0)/p)); // tinf1
        
        // tx    = tinf1 - (tinf1 - tHeMS) * (LTHe/Lx)**((p-1.0)/p)
        timescales[18] = timescales[16] - (timescales[16] - timescales[15]) * pow((LTHe/Lx), ((p-1.0)/p)); // tx
        
        // tinf2 = tx + ((1.0)/((q-1.0)*AHe*B)) * (B/Lx)**((q-1.0)/q)
        timescales[17] = timescales[18] + ((1.0)/((q-1.0)*AHe*B)) * pow((B/Lx), ((q-1.0)/q));  // tinf2
        
    }
    

    
    // DEBUGGING
    //std::cout << "tbl, tx, ty = " << timescales[19] << " " << timescales[20] << " " << timescales[21] << std::endl;
    
}

void Star::test(double time){
    
    calculateGiantBranchParameters(m_GBParams, m_an_coefficients, m_bn_coefficients, m_massCutoffs, m_Mass0, m_logMetallicityXi, m_stellarType, m_randomSeed);
    
    calculateTimescales(m_timescales, m_GBParams, m_an_coefficients, m_bn_coefficients, m_massCutoffs, m_Mass0, m_coreMass, m_COCoreMass, m_Metallicity, m_logMetallicityXi, m_stellarType, m_randomSeed);
    
    // To use set precision:
    // double f = 3.145;
    // std::cout << std::setprecision(9) << f << '\n';
    
    // testing
    
    std::cout << "testing = " << lifetimeSecondDredgeUp(m_bn_coefficients, m_timescales, m_GBParams, m_Mass, m_COCoreMass, m_massCutoffs, m_logMetallicityXi, m_randomSeed) << std::endl;
    
    //double tt = ageEvolvedHeliumGiantBranch(m_GBParams, m_timescales, m_Mass, 0.5);
    
    // Seems to not agree
    //std::cout << time << "\t" << luminosityMainSequence(m_an_coefficients, m_massCutoffs, m_timescales, m_LZAMS, m_Mass, m_logMetallicityXi, m_Metallicity, time) << "\t" << radiusMainSequence(m_an_coefficients, m_massCutoffs, m_timescales, m_RZAMS, m_Mass, time) << std::endl;
    
    // Seems to agree
    //     std::cout << "testing = " << calculate_B_alphaL(m_an_coefficients) << std::endl;
    // betaL, deltaL
    // LZAMS, RZAMS
    //std::cout << "tMS, tBGB = " << m_timescales[0] << "\t" << m_timescales[1] << std::endl;
    
    //std::cout << "LTMS, RTMS = " << luminosityEndMainSequence(m_an_coefficients, m_Mass) << "\t" << radiusEndMainSequence(m_an_coefficients, m_Mass, m_RZAMS) << std::endl;
    //std::cout << "alpha, beta, deltaL = " << calculate_alphaL(m_an_coefficients, m_Mass) << "\t" << calculate_betaL(m_an_coefficients, m_Mass) << "\t" << calculate_DeltaL(m_an_coefficients, m_massCutoffs, m_Mass) << std::endl;
    // alphaL
    //std::cout << "testing = " << std::setprecision(5) << m_an_coefficients[50 - 1] << " " << m_an_coefficients[52 - 1] << std::endl;
    
    // calculate_eta(m_Mass, m_Metallicity);
    
    //    double mu = std::max(0.5, (1.0 - 0.01 * std::max((m_an_coefficients[6 - 1]/pow(m_Mass, m_an_coefficients[7 - 1])), (m_an_coefficients[8 - 1] + m_an_coefficients[9 - 1]/pow(m_Mass, m_an_coefficients[10 - 1]))))) ;
    //    double t_hook = mu * m_timescales[1];
    //
    //    double epsilon = 0.01;
    //
    //    double bracketsTop      = time - (1.0 - epsilon)*t_hook;
    //    double bracketsBottom   = epsilon * t_hook;
    //    double brackets         = bracketsTop / bracketsBottom;
    //
    //    double tau1     = std::min(1.0, (time/t_hook));
    //    double tau2     = std::max(0.0, std::min(1.0, brackets));
    //    double tau      = time / m_timescales[0];
    
}

double Star::omegaBreak(double m_Mass, double m_Radius){
    /* 
    Break up frequency of star in yr-1 units, where [G] = 4*pi^2 AU^3 yr-2 Msol-1

    Parameters
    -----------
    m_Mass : double
        Mass in Msol
    m_Radius : double
        Radius in Rsol

    Returns
    --------
    omega : double
        Break up frequency
    */
	return sqrt(pow(2.0*pi,2.0)*(m_Mass)/pow(RsolToAU*m_Radius,3.0));	// 	Change pow(2.0*pi,2.0) to a constant Gsomething
}

double Star::limitMaximumChangeStar(const programOptions &options, const gsl_rng *r){
	//This function tests if the m_dt is small enough to resolve a maximum
	//change of 10% in the Radius of the star. To do this we first copy m_dt
	//because we do not want to update m_dt during the loop.
	double timestep = m_dt;
	
	//We will not go smaller than dynamical timescales
	//models need to be in hydrostatic equilibrium
	double t_dynamic = sqrt( pow(m_Radius,3) / (2*m_Mass*G ))/(31556926.0*1e6); //unitconversion seconds to Myr
	//To test the radial change we evolve a copy of the star and compare
	//its radius to the original star
	Star starCopy = *this;
	starCopy.evolveOneTimestep(timestep, true, options, r);

	//Calculate absolute radial change between original and evolved copy
	double radial_change = fabs(m_Radius - starCopy.m_Radius)/m_Radius;

	//If radial change too big we need to find a smaller timestep such
	//that the radius of the copy changes less than 10% compared to original
	//We do not want to get in an infinite loop so
	//define and count the number of iterations
	int iteration_counter = 0;
	int max_iterations =30;
	while (radial_change >= 0.01){
		if (timestep < t_dynamic){
			break;
			}
		
		iteration_counter +=1;
		timestep = 0.5*timestep;
		starCopy = *this;
		starCopy.evolveOneTimestep(timestep, false, options, r);
		radial_change = fabs(m_Radius-starCopy.m_Radius)/m_Radius;
		
		if (iteration_counter > max_iterations){
			if (debugging == true){
			std::cout << m_randomSeed << "\tMax number of iterations to resolve Radius reached (makeTimesep() in star.cpp), radial change (%) =" <<TAB<<radial_change<<TAB<<m_stellarType<< std::endl;}
			break;
			}
		if (m_stellarTypePrev != m_stellarType){
			if (debugging == true){
			std::cout << m_randomSeed << "\tStellar type changed when resolving R (maketimestep() star.cpp" << std::endl;}
			}
	
	}//close while
	
	return timestep;
}

double Star::makeTimestep(const programOptions &options, const gsl_rng *r){
    // Since Giant Branch (GB) parameters generally used to calculate core mass on the giant branch, should use M0
    calculateGiantBranchParameters(m_GBParams, m_an_coefficients, m_bn_coefficients, m_massCutoffs, m_Mass0, m_logMetallicityXi, m_stellarType, m_randomSeed);
    
    // Calculate timescales for the star
    calculateTimescales(m_timescales, m_GBParams, m_an_coefficients, m_bn_coefficients, m_massCutoffs, m_Mass0, m_coreMass, m_COCoreMass, m_Metallicity, m_logMetallicityXi, m_stellarType, m_randomSeed);


    //choosing timestep based on timescales
    m_dt = chooseTimestep(m_stellarType, m_Age, m_timescales, options);	

	//Coen Neijssel 10-01-2016
	//Want to resolve change star with a certain accurancy
	// Future work incorporate winds below in this function
    m_dt = limitMaximumChangeStar(options, r);



    //COEN NEIJSSEL 16/11/2016 --- added this function to cap timestep according to max 1% change in mass star.
    m_dtWind = m_dt;
    bool dtWind_boolean = false; //fudge to use masslosscaller
    double dtWind_mass = 
    massLossCaller(m_Mass, m_Mass0, m_coreMass, m_Luminosity, m_Radius, m_Mu, m_Metallicity, m_Temperature, m_logMetallicityXi, m_Mdot, m_dtWind, m_Age, options.luminousBlueVariableFactor, options.wolfRayetFactor, m_stellarType, isMassLoss, dtWind_boolean, massLossPrescription, debugging, m_an_coefficients, m_timescales); //dtWind_mass is not used. We need to run massLossCaller to update m_dtWind
    
    
    // Choose suitable timestep to return  
    
    return std::min(m_dt, m_dtWind);
}


void Star::evolveOneTimestep(double dt, bool useFakeMassLoss, const programOptions &options, const gsl_rng *r){
    /*
     Evolve the star just one timestep
     
     Should set radius, luminosity, core mass, temperature, mass loss etc as a function of time
     
     Should be a way to merge this with the function above since all it does is check whether to keep going or not, evolution is same
     
     Means you should put the bulk of the above into a function for one timestep?
     
     Maybe this is the function that is called by the evolveMaxTime so that you only have to fix things once.
     
     
     */
    
    // Update both the physical simulation time and the age of the star
    /*
    m_time    += dt;
    m_Age     += dt;
    m_dt       = dt;
*/

    
    // Record properties of star at previous timestep (test of general version)
    m_stellarTypePrev           = m_stellarType;
    m_timePrev                  = m_time;
    m_dtPrev                    = m_dt;
    m_tauPrev                   = m_tau;
    m_AgePrev                   = m_Age;
    m_nuclearBurningTimePrev    = m_nuclearBurningTime;
    m_MassPrev                  = m_Mass;
    m_Mass0Prev                 = m_Mass0;
    m_LuminosityPrev            = m_Luminosity;
    m_RadiusPrev                = m_Radius;
    m_TemperaturePrev           = m_Temperature;
    m_envMassPrev               = m_envMass;
    m_coreMassPrev              = m_coreMass;
    m_HeCoreMassPrev            = m_HeCoreMass;
    m_COCoreMassPrev            = m_COCoreMass;
    m_MuPrev                    = m_Mu;
    m_coreRadiusPrev            = m_coreRadius;
    m_MdotPrev                  = m_Mdot;
    m_omegaPrev                 = m_omega;
    m_omegaBreakPrev            = m_omegaBreak;
    m_angularMomentumPrev       = m_angularMomentum;

    //
    m_time    += dt;
    m_Age     += dt;
    m_dt       = dt;

    debugging = false;

	   
//    // booleans for fiddling with bits of the code whilst debugging.
//    bool debugging = false;                     // Whether to print debugging statements
//    bool isPerturbation = true;                 // Whether to use perturbations (disable for testing)
//    bool isMassLoss = true;                     // Whether to use mass loss (disable for testing)
//    bool writeOutput = false;                   // Whether to write parameters after timestep to standard output
    
    // Set up variables/output
    double tau          = 0.0;                          // Age after adjusting for mass loss
    double timePrime    = 0.0;                          // Relative time
//    double dt           = 0.0;                          // Declare variable for timestep
    
    // Problem declaring things in the switch/case
    bool TPAGBexists = true; // Declare a boolean for whether the TPAGB exists
    double McDU      = 0.0;  // Declare variable for McDU during dredge up
    double McBAGB    = 0.0;  // Declare variable for core mass at base of asymptotic giant branch (AGB)
    double McMax     = 0.0;  // Declare variable for maximum core mass on AGB
    double McPrime   = 0.0;  // Declare variable for dummy core mass on AGB
    double McSN      = 0.0;  // Declare variable for supernova core mass
    double Mdot      = 0.0;  // Declare variable for mass loss rate : should this be an intrinsic variable of the star i.e. m_Mdot ? probably since mass gain of second star will depend upon this rate.
    double McCOBAGB  = 0.0;  // Declare variable for CO core mass at base of the asymptotic giant branch (BAGB)
	double tDU       = 0.0;  // Declare variable for time of dredfe up in EAGB
    
    // At each timestep we want to calculate:
    // luminosity of star
    // radius of star
    // temperature of star
    // core mass of star
    // stellar type of star
    // mass of star due to mass loss
	
	// ALEJANDRO - 02/12/2016 - Attempt to fix updating the star if it lost all of it's mass
	if(m_Mass <= 0.0){
		m_stellarType = MASSLESS_REMNANT;
	}
	else{
    // Since Giant Branch (GB) parameters generally used to calculate core mass on the giant branch, should use M0
    calculateGiantBranchParameters(m_GBParams, m_an_coefficients, m_bn_coefficients, m_massCutoffs, m_Mass0, m_logMetallicityXi, m_stellarType, m_randomSeed);
    
    // Calculate timescales for the star
    calculateTimescales(m_timescales, m_GBParams, m_an_coefficients, m_bn_coefficients, m_massCutoffs, m_Mass0, m_coreMass, m_COCoreMass, m_Metallicity, m_logMetallicityXi, m_stellarType, m_randomSeed);
    }
	
    switch (m_stellarType) {
        case MS_LESS_THAN_07:
            
            // Calculate effective ZAMS parameters
            //m_RZAMS0 = RadiusZAMS(m_Mass0, m_logMetallicityXi);
            m_RZAMS0 = RadiusZAMS(m_Mass0, radius_coefficient_theta, radius_coefficient_iota, radius_coefficient_kappa, radius_coefficient_lamda, radius_coefficient_mu, radius_coefficient_nu, radius_coefficient_xi, radius_coefficient_omicron, radius_coefficient_Pi);

            //m_LZAMS0 = LuminosityZAMS(m_Mass0, m_logMetallicityXi);
            m_LZAMS0 = LuminosityZAMS(m_Mass0, luminosity_coefficient_alpha, luminosity_coefficient_beta, luminosity_coefficient_gamma, luminosity_coefficient_delta, luminosity_coefficient_epsilon, luminosity_coefficient_zeta, luminosity_coefficient_eta);
            
            if(m_Age < m_timescales[0]){
                
                // Do main sequence evolution
                
                tau = m_Age / m_timescales[0];  // tau = t / tMS
                
                m_Luminosity    = luminosityMainSequence(m_an_coefficients, m_massCutoffs, m_timescales, m_LZAMS0, m_Mass0, m_logMetallicityXi, m_Metallicity, m_Age, m_error, m_randomSeed);
                m_Radius        = radiusMainSequence(m_an_coefficients, m_massCutoffs, m_timescales, m_RZAMS0, m_Mass, m_Age, m_error, m_randomSeed);
                m_coreMass      = 0.0;
                m_Temperature   = calculateTemperature(m_Luminosity, m_Radius);
                
                // Convective envelope mass diminishes as a function of time
                m_envMass       = envelopeMass(m_Mass, m_coreMass, tau, m_stellarType);
                m_Mu            = perturbation_mu(m_Mass, m_coreMass, m_Luminosity, m_stellarType);
                
                if(debugging){
                    
                    std::cout << "MS LUMINOISTY: " << m_Luminosity << std::endl;
                    std::cout << "MS RADIUS: " << m_Radius << std::endl;
                    std::cout << "MS MASS: " << m_Mass << std::endl;
                    std::cout << "MS MASS0: " << m_Mass0 << std::endl;
                    std::cout << "MS AGE: " << m_Age << std::endl;
                    
                }
                
            }
            else{
                
                // Reset the values to those at the end of the main sequence and evolve star onto HG
                tau           = 1.0;
                m_Age         = m_timescales[0];
                
                m_Luminosity  = luminosityEndMainSequence(m_an_coefficients, m_Mass0, m_randomSeed);
                m_Radius      = radiusEndMainSequence(m_an_coefficients, m_Mass, m_RZAMS0, m_randomSeed);
                m_coreMass    = 0.0;
                m_Temperature = calculateTemperature(m_Luminosity, m_Radius);
                
                // Convective envelope mass vanishes
                m_envMass     = envelopeMass(m_Mass, m_coreMass, tau, m_stellarType);
                m_Mu          = perturbation_mu(m_Mass, m_coreMass, m_Luminosity, m_stellarType);
                
                m_stellarType = HERTZSPRUNG_GAP;
                
                if(debugging){
                    // Tag to jump to end of MS
                    std::cout << "End of MS" << std::endl;
                }
            }
            
            break;
            
        case MS_MORE_THAN_07:
            
            // Calculate effective ZAMS parameters
            //m_RZAMS0 = RadiusZAMS(m_Mass0, m_logMetallicityXi);
            m_RZAMS0 = RadiusZAMS(m_Mass0, radius_coefficient_theta, radius_coefficient_iota, radius_coefficient_kappa, radius_coefficient_lamda, radius_coefficient_mu, radius_coefficient_nu, radius_coefficient_xi, radius_coefficient_omicron, radius_coefficient_Pi);

            //m_LZAMS0 = LuminosityZAMS(m_Mass0, m_logMetallicityXi);
            m_LZAMS0 = LuminosityZAMS(m_Mass0, luminosity_coefficient_alpha, luminosity_coefficient_beta, luminosity_coefficient_gamma, luminosity_coefficient_delta, luminosity_coefficient_epsilon, luminosity_coefficient_zeta, luminosity_coefficient_eta);
            
            if(m_Age < m_timescales[0]){
                
                // Do main sequence evolution
                
                tau = m_Age / m_timescales[0]; // tau = t / tMS
                
                debugging = false;
                
                if(debugging){
                    std::cout << "Mass, Mass0 = " << m_Mass << " " << m_Mass0 << std::endl;
                    std::cout << "tMS = " << m_timescales[0] << std::endl;
                    std::cout << "time, tau = " << m_Age << " " << tau << std::endl;
					std::cout << "Luminosity, Temperature, Radius:\t" << m_Luminosity << "\t" << m_Temperature << "\t" << m_Radius << std::endl;
                }
                
                m_Luminosity    = luminosityMainSequence(m_an_coefficients, m_massCutoffs, m_timescales, m_LZAMS0, m_Mass0, m_logMetallicityXi, m_Metallicity, m_Age, m_error, m_randomSeed);
                m_Radius        = radiusMainSequence(m_an_coefficients, m_massCutoffs, m_timescales, m_RZAMS0, m_Mass, m_Age, m_error, m_randomSeed);
                m_coreMass      = 0.0;
                m_Temperature   = calculateTemperature(m_Luminosity, m_Radius);
                
                // Convective envelope mass diminishes with time
                m_envMass       = envelopeMass(m_Mass, m_coreMass, tau, m_stellarType);
                m_Mu            = perturbation_mu(m_Mass, m_coreMass, m_Luminosity, m_stellarType);
                
                if(debugging){
                    
                    std::cout << "MS LUMINOISTY: " << m_Luminosity << std::endl;
                    std::cout << "MS log10(LUMINOISTY): " << log10(m_Luminosity) << std::endl;
                    std::cout << "MS RADIUS: " << m_Radius << std::endl;
                    std::cout << "MS MASS: " << m_Mass << std::endl;
                    std::cout << "MS MASS0: " << m_Mass0 << std::endl;
                    std::cout << "MS AGE: " << m_Age << std::endl;
                    
                }
                
            }
            else{
                
                // Reset the values to those at the end of the main sequence and evolve star to HG
                tau           = 1.0;
                m_Age         = m_timescales[0];
                
                m_Luminosity  = luminosityEndMainSequence(m_an_coefficients, m_Mass0, m_randomSeed);
                m_Radius      = radiusEndMainSequence(m_an_coefficients, m_Mass, m_RZAMS0, m_randomSeed);
                m_coreMass    = 0.0;
                m_Temperature = calculateTemperature(m_Luminosity, m_Radius);
                
                // Convective envelope mass vanishes
                m_envMass     = envelopeMass(m_Mass, m_coreMass, tau, m_stellarType);
                m_Mu          = perturbation_mu(m_Mass, m_coreMass, m_Luminosity, m_stellarType);
                
                m_stellarType = HERTZSPRUNG_GAP;
                
                // TAG to jump to end of MS -- Why does it not end at the same time as the fortran code?
                if(debugging){
                    std::cout << "END OF MS" << std::endl;
                }
                
                if(debugging){
                    
                    std::cout << "MS LUMINOISTY: " << m_Luminosity << std::endl;
                    std::cout << "MS log10(LUMINOISTY): " << log10(m_Luminosity) << std::endl;
                    std::cout << "MS RADIUS: " << m_Radius << std::endl;
                    std::cout << "MS MASS: " << m_Mass << std::endl;
                    std::cout << "MS MASS0: " << m_Mass0 << std::endl;
                    std::cout << "MS AGE: " << m_Age << std::endl;
                    
                }
                
            }
            
            break;
            
        case HERTZSPRUNG_GAP:
            
            if(m_Age < m_timescales[1]){ // t < tBGB
                // Do Hertzsprung Gap (HG) evolution
                
                tau = (m_Age - m_timescales[0])/(m_timescales[1] - m_timescales[0]);  //tau = (t - tMS)/(tBGB - tMS);
                
                if(debugging){
                    std::cout << "tauHG = " << tau << std::endl;
                    std::cout << "age = " << m_Age << std::endl;
                    std::cout << "tMS, tBGB = " << m_timescales[0] << " " << m_timescales[1] << std::endl;
                }
                
                m_Luminosity  = luminosityHertzsprungGap(m_an_coefficients, m_bn_coefficients, m_massCutoffs, m_Mass0, tau, m_randomSeed);
                
                m_coreMass    = coreMassHertzsprungGap(m_an_coefficients, m_bn_coefficients, m_massCutoffs, m_GBParams, tau, m_Mass0, m_coreMass, m_Metallicity, m_logMetallicityXi, m_randomSeed);
                
                m_HeCoreMass  = m_coreMass;
                
                m_Radius      = radiusHertzsprungGap(m_an_coefficients, m_bn_coefficients, m_massCutoffs, m_Mass, m_coreMass, m_RZAMS0, tau, m_logMetallicityXi, m_randomSeed);
                
                m_Mu          = perturbation_mu(m_Mass, m_coreMass, m_Luminosity, m_stellarType);
                
                if(isPerturbation){
                    perturbedLuminosityRadius(m_Luminosity, m_Radius, m_Mu, m_an_coefficients, m_bn_coefficients, m_GBParams, m_massCutoffs, m_timescales, m_Age, m_Mass, m_coreMass, m_HeCoreMass, m_COCoreMass, m_Metallicity, m_logMetallicityXi, m_stellarType, m_randomSeed);
                }
                
                m_Temperature = calculateTemperature(m_Luminosity, m_Radius);
                
                // Convective envelope gradually establishes itself
                m_envMass  = envelopeMass(m_Mass, m_coreMass, tau, m_stellarType);
                
                if(m_coreMass > m_Mass){ // EnvelopeLoss
                    
                    // Star has lost its envelope, forms a HeWD if has degenerate core (M < MHeF) or zero-age naked helium star
                    // Described just above Equation 76
                    // Testing envelope loss function
                    /*
                    if(m_Mass < m_massCutoffs[1]){
                        
                        // Reset mass/age parameters
                        m_coreMass    = m_Mass;
                        
                        m_Age         = 0;  // Check eq. 76, Hurley et al. 2000
                        
                        m_stellarType = HELIUM_WHITE_DWARF;
                        
                    }
                    else{
                        
                        // Zero age Helium MS star // see section 6.1
                        
                        m_coreMass    = m_Mass;
                        m_Mass0       = m_Mass;
                        
                        m_Age         = 0;  // Check eq. 76, Hurley et al. 2000
                        
                        m_stellarType = NAKED_HELIUM_STAR_MS;
                        
                        // Something for calculate radius and luminosity here?
                        
                        m_Radius      = radiusHeliumMainSequenceZAMS(m_Mass);
                        m_Luminosity  = luminosityHeliumMainSequenceZAMS(m_Mass);
                        std::cout << "Maximun success\t"<<m_Radius<<"\t"<<m_Mass<<std::endl;
                    }
                     */
                    modifyStarAfterLosingEnvelope(m_stellarType, m_Mass);
                    
                }
                
            }
            else{
                // time = end of HG, evolves on FGB
                tau   = 1.0;
                
                
                // Reset parameters to those at end of HG
                
                m_Luminosity  = luminosityEndHertzsprungGap(m_an_coefficients, m_bn_coefficients, m_Mass0, m_massCutoffs, m_randomSeed);
                m_coreMass    = coreMassEndHertzsprungGap(m_an_coefficients, m_bn_coefficients, m_GBParams, m_Mass0, m_massCutoffs, m_logMetallicityXi, m_randomSeed);
                m_HeCoreMass  = m_coreMass;
                m_COCoreMass  = 0.0;
                m_Radius      = radiusEndHertzsprungGap(m_an_coefficients, m_bn_coefficients, m_Mass, m_coreMass, m_massCutoffs, m_logMetallicityXi, m_randomSeed);
                
                m_Mu          = perturbation_mu(m_Mass, m_coreMass, m_Luminosity, m_stellarType);
                
                if(isPerturbation){
                    perturbedLuminosityRadius(m_Luminosity, m_Radius, m_Mu, m_an_coefficients, m_bn_coefficients, m_GBParams, m_massCutoffs, m_timescales, m_Age, m_Mass, m_coreMass, m_HeCoreMass, m_COCoreMass, m_Metallicity, m_logMetallicityXi, m_stellarType, m_randomSeed);
                }
                
                m_Temperature = calculateTemperature(m_Luminosity, m_Radius);
                
                if(m_Mass0 < m_massCutoffs[2]){ // M < MFGB
                    
                    m_stellarType = FIRST_GIANT_BRANCH;
                    m_Age = m_timescales[1];
                    if(debugging){
                        std::cout << "M < MFGB, Starting FGB" << std::endl;
                    }
                }
                else{
                    
                    m_stellarType = CORE_HELIUM_BURNING;
                    m_Age = m_timescales[2];
                    if(debugging){
                        std::cout << "M > MFGB, Starting CHeB" << std::endl;
                    }
                }
                
                m_envMass     = envelopeMass(m_Mass, m_coreMass, tau, m_stellarType);
                
                if(debugging){
                    // DEBUGGING
                    std::cout << "END OF HG" << std::endl;
                    
                    std::cout << "RminHe = " << log10(minimumRadiusCoreHeliumBurning(m_bn_coefficients, m_Mass, m_coreMass, m_massCutoffs, m_logMetallicityXi, m_randomSeed)) << std::endl;
                    std::cout << "RHeI = "   << log10(radiusHeliumIgnition(m_bn_coefficients, m_Mass, m_coreMass, m_massCutoffs, m_logMetallicityXi, m_randomSeed)) << std::endl;
                    std::cout << "Rx = "     << log10(radiusCoreHeliumBurningX(m_bn_coefficients, m_Mass, m_coreMass, m_massCutoffs, m_logMetallicityXi, m_randomSeed)) << std::endl;
                }
                
            }
            
            break;
            
        case FIRST_GIANT_BRANCH:
            
            // M >~ 16 Msol will skip FGB
            if(m_Mass0 < m_massCutoffs[2]){ // M < MFGB
                
                if(m_Age < m_timescales[2]){ // t < tHeI get from timescales
                    
                    // Do FGB evolution
                    
                    m_Luminosity  = luminosityFirstGiantBranch(m_GBParams, m_timescales, m_Age, m_randomSeed);
                    m_coreMass    = coreMassGiantBranch(m_an_coefficients, m_bn_coefficients, m_GBParams, m_timescales, m_Age, m_Mass0, m_massCutoffs, m_logMetallicityXi, m_randomSeed);
                    m_HeCoreMass  = m_coreMass;
                    m_COCoreMass  = 0.0;
                    m_Radius      = radiusFirstGiantBranch(m_bn_coefficients, m_Mass, m_Luminosity, m_logMetallicityXi, m_randomSeed);
                    
                    m_Mu          = perturbation_mu(m_Mass, m_coreMass, m_Luminosity, m_stellarType);
                    
                    if(isPerturbation){
                        perturbedLuminosityRadius(m_Luminosity, m_Radius, m_Mu, m_an_coefficients, m_bn_coefficients, m_GBParams, m_massCutoffs, m_timescales, m_Age, m_Mass, m_coreMass, m_HeCoreMass, m_COCoreMass, m_Metallicity, m_logMetallicityXi, m_stellarType, m_randomSeed);
                    }
                    
                    m_Temperature = calculateTemperature(m_Luminosity, m_Radius);
                    
                    if(m_coreMass > m_Mass){
                     modifyStarAfterLosingEnvelope(m_stellarType, m_Mass);
                    }
                    
                }
                else{
                    // Reset to values at end of FGB
                    
                    m_Age = m_timescales[2];
                    
                    m_Luminosity  = luminosityFirstGiantBranch(m_GBParams, m_timescales, m_Age, m_randomSeed);
                    m_coreMass    = coreMassGiantBranch(m_an_coefficients, m_bn_coefficients, m_GBParams, m_timescales, m_Age, m_Mass0, m_massCutoffs, m_logMetallicityXi, m_randomSeed);
                    m_HeCoreMass  = m_coreMass;
                    m_COCoreMass  = 0.0;
                    m_Radius      = radiusFirstGiantBranch(m_bn_coefficients, m_Mass, m_Luminosity, m_logMetallicityXi, m_randomSeed);
                    
                    m_Mu          = perturbation_mu(m_Mass, m_coreMass, m_Luminosity, m_stellarType);
                    
                    if(isPerturbation){
                        perturbedLuminosityRadius(m_Luminosity, m_Radius, m_Mu, m_an_coefficients, m_bn_coefficients, m_GBParams, m_massCutoffs, m_timescales, m_Age, m_Mass, m_coreMass, m_HeCoreMass, m_COCoreMass, m_Metallicity, m_logMetallicityXi, m_stellarType, m_randomSeed);
                    }
                    
                    m_Temperature = calculateTemperature(m_Luminosity, m_Radius);
                    
                    m_stellarType = CORE_HELIUM_BURNING;
                    
                    if(m_Mass0 < m_massCutoffs[1]){
                        
                        if(debugging){
                            // For LM star at ZAHB (end of GB/begining of CHeB) due to helium flash when doing mass loss
                            std::cout << "End of GB - LM M0 = " << m_Mass0 << " < " << m_massCutoffs[1] << std::endl;
                            std::cout << "At ZAHB, reset m_Mass0 = " << m_Mass << " and t = " << lifetimeHeliumIgnition(m_bn_coefficients, m_GBParams, m_timescales, m_Mass0, m_massCutoffs, m_randomSeed) << std::endl;
                        }
                        
                        m_Mass0 = m_Mass;
                        
                        calculateGiantBranchParameters(m_GBParams, m_an_coefficients, m_bn_coefficients, m_massCutoffs, m_Mass0, m_logMetallicityXi, m_stellarType, m_randomSeed);
                        calculateTimescales(m_timescales, m_GBParams, m_an_coefficients, m_bn_coefficients, m_massCutoffs, m_Mass0, m_coreMass, m_COCoreMass, m_Metallicity, m_logMetallicityXi, m_stellarType, m_randomSeed);
                        
                        if(debugging){
                            // DEBUGGING
                            std::cout << "tHeI(M), tHeI(M0) = " << lifetimeHeliumIgnition(m_bn_coefficients, m_GBParams, m_timescales, m_Mass, m_massCutoffs, m_randomSeed) << " " << lifetimeHeliumIgnition(m_bn_coefficients, m_GBParams, m_timescales, m_Mass, m_massCutoffs, m_randomSeed) << std::endl;
                        }
                        
                        m_Age = m_timescales[2]; //lifetimeHeliumIgnition(m_bn_coefficients, m_GBParams, m_timescales, m_Mass0, m_massCutoffs);
                        
                    }
                    
                    // Also Note that for M>zpars(3) MFGB there is no GB as the star goes from HG -> CHeB -> AGB. So in effect tscls(1) tBGB refers to the time of Helium ignition and not the BGB.
                    
                    if(debugging){
                        std::cout << "End of FGB" << std::endl;
                    }
                }
            }
            else{
                
                // time = tHeI?
                //m_Age = m_timescales[];
                
                m_stellarType = CORE_HELIUM_BURNING;
                
                if(debugging){
                    std::cout << "Start of CHEB" << std::endl;
                }
                
            }
            
            break;
            
        case CORE_HELIUM_BURNING:
			//debugging = true;
			
			if(debugging){
				std::cout << "Radius at CHeB: " << m_Radius << std::endl;
			}
            if(m_Age < m_timescales[2] + m_timescales[3]){ // t < t_HeI + tHe
                // Do Core helium burning evolution
                
                tau = (m_Age - m_timescales[2])/m_timescales[3]; // (t - tHeI)/tHe
                
                // DEBUGGING
                if(debugging){
                    std::cout << "tau, age, tHeI, tHe, tHeAlt = " << tau << " " << m_Age << " " << m_timescales[2] << " " << m_timescales[3] << " " << lifetimeCoreHeliumBurning(m_an_coefficients, m_bn_coefficients, m_Mass0, m_coreMass, m_massCutoffs, m_randomSeed) << std::endl;
                    
                    //std::cout << "R, Ralt = " << m_Radius << " " << radiusCoreHeliumBurningTauAlt(tau, m_Mass, m_coreMass, m_timescales, m_bn_coefficients, m_massCutoffs) << std::endl;
                }
                m_coreMass    = coreMassHeliumBurning(m_bn_coefficients, m_GBParams, tau, m_Mass0, m_massCutoffs, m_logMetallicityXi, m_randomSeed);
                m_HeCoreMass  = m_coreMass;
                m_COCoreMass  = 0.0;
                m_Luminosity  = luminosityCoreHeliumBurning(m_bn_coefficients, m_timescales, m_massCutoffs, m_Mass0, m_coreMass, tau, m_logMetallicityXi, m_error, m_randomSeed);
 
                // KLUDGE for now
                if(m_coreMass > m_Mass){
                    m_coreMass   = m_Mass;
                    m_HeCoreMass = m_coreMass;
                }
                
                m_Radius      = radiusCoreHeliumBurning(m_bn_coefficients, m_massCutoffs, m_timescales, tau, m_Mass, m_coreMass, m_Luminosity, m_logMetallicityXi, m_error, m_randomSeed, m_Mass0);

                //debugging = true;
				if(debugging){
					std::cout << "After eval radius at CHeB: " << std::endl;
					std::cout << "Radius: " << m_Radius << std::endl;
					std::cout << "tau: " << m_tau << std::endl;
					std::cout << "Mass: " << m_Mass << std::endl;
					std::cout << "coreMass: " << m_coreMass << std::endl;
					std::cout << "Luminosity: " << m_Luminosity << std::endl;
                }
				
                m_Mu          = perturbation_mu(m_Mass, m_coreMass, m_Luminosity, m_stellarType);
                m_envMass     = envelopeMass(m_Mass , m_coreMass, tau, m_stellarType);
                
                if(isPerturbation){
                    perturbedLuminosityRadius(m_Luminosity, m_Radius, m_Mu, m_an_coefficients, m_bn_coefficients, m_GBParams, m_massCutoffs, m_timescales, m_Age, m_Mass, m_coreMass, m_HeCoreMass, m_COCoreMass, m_Metallicity, m_logMetallicityXi, m_stellarType, m_randomSeed);
                }
                
                if(debugging){

                    std::cout << "After perturbation" << std::endl;
                    std::cout << "R_CHeB: " << m_Radius << std::endl;
                    std::cout << "Luminosity: " << m_Luminosity << std::endl;

                }

                m_Temperature = calculateTemperature(m_Luminosity, m_Radius);
                
                // if core mass > total mass then need to set stellar type to an evolved naked helium star and evolve that
                if(m_coreMass >= m_Mass){// EnvelopeLoss
					if(debugging){
						std::cout << "Inside EnvelopeLoss in CHeB" << std::endl;
					}
                    modifyStarAfterLosingEnvelope(m_stellarType, m_Mass);
                }
				
				if(debugging){
				std::cout << "Radius at end CHeB: " << m_Radius << std::endl;
				}
            }
            else{
                
                // Set things to values at end of CHeB/BAGB
                m_Age = m_timescales[2] + m_timescales[3];
                
                tau = (m_Age - m_timescales[2])/m_timescales[3]; // (t - tHeI)/tHe
                
                m_coreMass    = coreMassBaseAsymptoticGiantBranch(m_bn_coefficients, m_Mass0, m_randomSeed);
                m_HeCoreMass  = m_coreMass;
                m_COCoreMass  = 0.0;
                m_Luminosity  = luminosityBaseAsymptoticGiantBranch(m_bn_coefficients, m_Mass0, m_massCutoffs, m_randomSeed);
                m_Radius      = radiusBaseAsymptoticGiantBranch(m_bn_coefficients, m_Mass, m_Luminosity, m_massCutoffs, m_randomSeed);
                
                m_Mu          = perturbation_mu(m_Mass, m_coreMass, m_Luminosity, m_stellarType);
                
                if(debugging){
                    std::cout << "L, R, Mc = " << log10(m_Luminosity) << " " << log10(m_Radius) << " " << m_coreMass << std::endl;
                    
                    // DEBUGGING
                    std::cout << "END OF CHeB, START OF EAGB i.e. BAGB" << std::endl;
                    
                }
                
                if(isPerturbation){
                    perturbedLuminosityRadius(m_Luminosity, m_Radius, m_Mu, m_an_coefficients, m_bn_coefficients, m_GBParams, m_massCutoffs, m_timescales, m_Age, m_Mass, m_coreMass, m_HeCoreMass, m_COCoreMass, m_Metallicity, m_logMetallicityXi, m_stellarType, m_randomSeed);
                    
                    if(debugging){
                        std::cout << "L', R' = " << log10(m_Luminosity) << " " << log10(m_Radius) << std::endl;
                    }
                }
                
                m_Temperature = calculateTemperature(m_Luminosity, m_Radius);
                
                if(debugging){
                    std::cout << "T = " << log10(m_Temperature * Tsol) << std::endl;
                }
                
                if(debugging){
                    std::cout << "End of CHeB, start of EAGB" << std::endl;
                }
                
                // would move this up but NAN
                m_stellarType = EARLY_ASYMPTOTIC_GIANT_BRANCH;
                
                // Need to assign variables at Base of the Asymptotic Giant Branch (BAGB)
                m_HeCoreMass  = m_coreMass;  // Remains constant throughout EAGB.
            }

            // Why no supernova for CHeB stars?

            break;
            
        case EARLY_ASYMPTOTIC_GIANT_BRANCH:
            
            // Get parameters
            McBAGB   = m_GBParams[9];     // Value that McHe remains constant at
			if(McBAGB <= 0.8){
				McDU = McBAGB;
			}
			else if(McBAGB > 0.8 and McBAGB < 2.25){
				//instantaneous dredge up (some part of the core mixes with envelope this reducing core mass)
				McDU = 0.44*McBAGB + 0.448 ; //eq(69) 
			}
			else if(McBAGB >= 2.25){
                TPAGBexists = false;
			}
			tDU = lifetimeSecondDredgeUp(m_bn_coefficients, m_timescales, m_GBParams, m_Mass0, McDU, m_massCutoffs, m_logMetallicityXi, m_randomSeed);
            McCOBAGB = coreMassAsymptoticGiantBranch(m_GBParams, m_timescales, m_timescales[2] + m_timescales[3], m_randomSeed); // tBAGB
            McSN     = maximumCoreMassSN(m_GBParams);
            McSN     = std::max(McSN, 1.05*McCOBAGB); // Hack from Hurley fortran code, doesn't seem to be in the paper
            if(debugging){
                // DEBUGGING
                std::cout << "McBAGB = " << McBAGB << std::endl;
                std::cout << "McCOBAGB = " << McCOBAGB << std::endl;
                std::cout << "McDU = " << McDU << std::endl;
                std::cout << "McSN = " << McSN << std::endl;
                std::cout << "McCO = " << m_COCoreMass << std::endl;
                std::cout << "MendHack = " << 1.05*McSN << std::endl;
                // Why is the time step here so small, is essentially 0? Why is core mass not growing.?
            }
            
			//Check if at any moment in time the core is heavy enough to go supernova 
			//if not evolve star else => supernova remnant
			if(m_COCoreMass < McSN){
		        // Regular EAGB evolution (the CO-core has not consumed the entire HeCore yet)
					
					if (m_Age >= tDU and TPAGBexists){
							m_Age = m_timescales[13];
							m_COCoreMass  = McDU;
				            m_coreMass    = m_COCoreMass;
				            m_Luminosity  = coreMassLuminosityRelation(m_GBParams, m_COCoreMass, m_randomSeed);
				            m_Radius      = radiusAsymptoticGiantBranch(m_bn_coefficients, m_Mass, m_Luminosity, m_massCutoffs, m_randomSeed);
				            m_Mu          = perturbation_mu(m_Mass, m_coreMass, m_Luminosity, m_stellarType);
				            
				            if(isPerturbation){perturbedLuminosityRadius(m_Luminosity, m_Radius, m_Mu, m_an_coefficients,
											 m_bn_coefficients, m_GBParams, m_massCutoffs, m_timescales, m_Age, m_Mass,
											 m_coreMass, m_HeCoreMass, m_COCoreMass, m_Metallicity, m_logMetallicityXi, m_stellarType, m_randomSeed);
				            }
				            
				            m_Temperature = calculateTemperature(m_Luminosity, m_Radius);
		                	m_stellarType = THERMALLY_PULSING_ASYMPTOTIC_GIANT_BRANCH;
					}
					else{
						if(m_coreMass > m_Mass){ // Total mass now equal to core mass i.e. envelope is lost
		               		modifyStarAfterLosingEnvelope(m_stellarType, m_Mass);
		            	}
						else{
							m_coreMass    = McBAGB;
						    m_COCoreMass  = coreMassAsymptoticGiantBranch(m_GBParams , m_timescales, m_Age, m_randomSeed);
							m_Luminosity  = luminosityAsymptoticGiantBranch(m_GBParams, m_COCoreMass);
						    m_Radius      = radiusAsymptoticGiantBranch(m_bn_coefficients, m_Mass, m_Luminosity, m_massCutoffs, m_randomSeed);
						    m_Mu          = perturbation_mu(m_Mass, m_coreMass, m_Luminosity, m_stellarType);
						    if(isPerturbation){perturbedLuminosityRadius(m_Luminosity, m_Radius, m_Mu, m_an_coefficients, 
											 m_bn_coefficients, m_GBParams, m_massCutoffs, m_timescales, m_Age, m_Mass, 
											m_coreMass, m_HeCoreMass, m_COCoreMass, m_Metallicity, m_logMetallicityXi, m_stellarType, m_randomSeed);
				        	}
				                
				        	m_Temperature = calculateTemperature(m_Luminosity, m_Radius);
				        	tau = 0.0; // dummy tau to calculate envelope mass
				        	m_envMass     = envelopeMass(m_Mass, m_coreMass, tau, m_stellarType);
						}
					}
					
				}//close regular EAGB evolution



		    if(m_COCoreMass > McSN){ // Supernova

                //debugging = true;

                //m_coreMass    = McBAGB;
                //m_HeCoreMass  = McBAGB;
                //m_COCoreMass  = McBAGB;

		        if(debugging){
		            std::cout << "End of EAGB - Supernova" << std::endl;
                    std::cout << "McBAGB = " << McBAGB << std::endl;
                    std::cout << "HeCoreMass = " << m_HeCoreMass << std::endl;
                    std::cout << "coreMass = " << m_coreMass << std::endl;
                    std::cout << "COCoreMass = " << m_COCoreMass << std::endl;
		        }

		        supernovaMain(&m_totalMassAtCompactObjectFormation, &m_HeCoreMassAtCompactObjectFormation,
                   &m_COCoreMassAtCompactObjectFormation,&m_coreMassAtCompactObjectFormation, &m_Mass0,
                   &m_Age, &m_envMass, &m_coreMass, &m_HeCoreMass, &m_COCoreMass, &m_coreRadius, &m_Mass, m_MZAMS, m_GBParams[9],
                   &m_fallback, &m_stellarType, &m_Radius, &m_Luminosity, &m_Temperature, &flagSN,
                   &flagECSN, &flagExperiencedECSN, &flagExperiencedCCSN, &flagExperiencedPISN, &flagExperiencedPPISN,
                   &m_drawnKickVelocity, &m_kickVelocity, 
                   &m_error, &m_flagHrichSN, &m_flagHpoorSN, m_SNengine, m_randomSeed, &m_angularMomentum, &m_pulsarSpinPeriod, &m_pulsarSpinFrequency, &m_pulsarMagneticField, &m_pulsarSpinDownRate, &m_momentInertia, options, r);
		    }//closing Supernova 

               

		break;





       case THERMALLY_PULSING_ASYMPTOTIC_GIANT_BRANCH:
            
            if(TPAGBexists){ // Check TPAGB exists
                
                McBAGB   = m_GBParams[9];     // Value that McHe remains constant at
                McSN  = maximumCoreMassSN(m_GBParams);
                McMax = maximumCoreMass(m_Mass);
                
                // debugging
                if(debugging){
                    
                    std::cout << "TPAGB: mass, core mass = " << m_Mass << " " << m_coreMass << std::endl;
                    std::cout << "TPAGB: McSN, McMax = " << McSN << " " << McMax << std::endl;
                    
                }
                
                if(m_COCoreMass < std::min(McSN, m_Mass)){


                    
                    // Do TPAGB evolution
                    
                    McPrime       = coreMassPrimeThermallyPulsatingAsymptoticGiantBranch(m_GBParams, m_timescales, m_Age, m_randomSeed);
                    m_Luminosity  = luminosityThermallyPulsatingAsymptoticGiantBranch(m_GBParams, McPrime, m_randomSeed);
                    m_COCoreMass  = coreMassThermallyPulsatingAsymptoticGiantBranch(McPrime, m_GBParams, m_Mass0, m_randomSeed);
                    m_coreMass    = m_COCoreMass;
                    m_Radius      = radiusAsymptoticGiantBranch(m_bn_coefficients, m_Mass, m_Luminosity, m_massCutoffs, m_randomSeed);
                    m_Mu          = perturbation_mu(m_Mass, m_coreMass, m_Luminosity, m_stellarType);
                    
                    if(isPerturbation){
                        perturbedLuminosityRadius(m_Luminosity, m_Radius, m_Mu, m_an_coefficients, m_bn_coefficients, m_GBParams, m_massCutoffs, m_timescales, m_Age, m_Mass, m_coreMass, m_HeCoreMass, m_COCoreMass, m_Metallicity, m_logMetallicityXi, m_stellarType, m_randomSeed);
                    }
                    
                    m_Temperature = calculateTemperature(m_Luminosity, m_Radius);                  
                    // DEBUGGING
                    if(debugging){
                        std::cout << "TEST TPAGB: McSN, McMax, McPrime, McBAGB, L, Mc = " << McSN << " " << McMax << " " << McPrime << " " << McBAGB << " " << m_Luminosity << " " << m_COCoreMass << std::endl;
                    }
                    
                }
                if(m_COCoreMass > std::min(McSN, m_Mass)){ // End of TPAGB
                    // End of TPAGB either when core mass == stellar mass (loss of envelope)
                    // or when core mass == McSn
                    // so while McCO < min(McSN, M)
                    // Continuation condition FINISH ME to do with core mass = stellar mass becomes COWD
                    // 1.6 < MCBAGB < 2.25 (Mup, Mec) SN leaves NS via EC SN
                    // if McBAGB > 2.25 SN results in a NS being formed.
                    // if McSN > 7.0 (equivalent to McBAGB > 9.52) collapses to BH -- will be modified by Fryer stuff
                    // work out whether end result is a supernova or white dwarf or helium star.
                    // Make sure Mc !> M
                    
                    if(debugging){
                        
                        std::cout << "TEST END TPAGB: McSN, McMax, McPrime, McBAGB, L, Mc = " << McSN << " " << McMax << " " << McPrime << " " << McBAGB << " " << m_Luminosity << " " << m_COCoreMass << std::endl;
                        
                    }
                    
                    if(m_COCoreMass >= McSN and m_COCoreMass < m_Mass){

                        //debugging = true;

                        if(debugging){
                            std::cout << "end of TPAGB" << std::endl;
                        }

                        supernovaMain(&m_totalMassAtCompactObjectFormation, &m_HeCoreMassAtCompactObjectFormation,
                   &m_COCoreMassAtCompactObjectFormation,&m_coreMassAtCompactObjectFormation, &m_Mass0,
                   &m_Age, &m_envMass, &m_coreMass, &m_HeCoreMass, &m_COCoreMass, &m_coreRadius, &m_Mass, m_MZAMS, m_GBParams[9],
                   &m_fallback, &m_stellarType, &m_Radius, &m_Luminosity, &m_Temperature, &flagSN,
                   &flagECSN, &flagExperiencedECSN, &flagExperiencedCCSN, &flagExperiencedPISN, &flagExperiencedPPISN, &m_drawnKickVelocity, &m_kickVelocity, 
                   &m_error,  &m_flagHrichSN, &m_flagHpoorSN, m_SNengine, m_randomSeed, &m_angularMomentum, &m_pulsarSpinPeriod, &m_pulsarSpinFrequency, &m_pulsarMagneticField, &m_pulsarSpinDownRate, &m_momentInertia, options, r);
                    }
                    else{
                        m_COCoreMass = m_Mass;
                        modifyStarAfterLosingEnvelope(m_stellarType, m_Mass);
                        if(debugging){
                            std::cout << "Loss of envelope: M, Mc, McBAGB = " << m_Mass << " " << m_coreMass << " " << McBAGB << std::endl;
                        }
                    }
                    
                }
                
            }
            else{
                
                std::cout << "TPAGB shouldn't exist" << std::endl;
                
            }
            
            break;
            
        case NAKED_HELIUM_STAR_MS:
            
            // Do Naked HeMS evolution
            tau = m_Age / m_timescales[15]; // work out tau = t/tHeMS where t is counted from 0 at He ZAMS
            
            if(debugging){
                std::cout << "HeMS : tau = " << tau << std::endl;
            }
            
            if(tau == 0.0){
                
                if(debugging){
                    std::cout << "If 1) tau = 0" << std::endl;
                }
                
                // He ZAMS
                m_Luminosity  = luminosityHeliumMainSequenceZAMS(m_Mass0);
                m_Radius      = radiusHeliumMainSequenceZAMS(m_Mass);
                
                m_HeCoreMass  = m_Mass;
                m_COCoreMass  = 0.0;
                m_coreMass    = m_COCoreMass;
                
                m_Mu          = perturbation_mu(m_Mass, m_coreMass, m_Luminosity, m_stellarType);
                
                if(isPerturbation){
                    perturbedLuminosityRadius(m_Luminosity, m_Radius, m_Mu, m_an_coefficients, m_bn_coefficients, m_GBParams, m_massCutoffs, m_timescales, m_Age, m_Mass, m_coreMass, m_HeCoreMass, m_COCoreMass, m_Metallicity, m_logMetallicityXi, m_stellarType, m_randomSeed);
                }
                
                m_Temperature = calculateTemperature(m_Luminosity, m_Radius);
            }
            else if(tau > 0.0 and tau < 1.0){
                
                if(debugging){
                    std::cout << "If 2) 0 < tau < 1" << std::endl;
                }
                
                // core mass constant
                m_Luminosity  = luminosityHeliumMainSequence(tau, m_Mass0);
                m_Radius      = radiusHeliumMainSequence(tau, m_Mass);
                
                m_HeCoreMass  = m_Mass;
                m_COCoreMass  = 0.0;
                m_coreMass    = m_COCoreMass;
                
                m_Mu          = perturbation_mu(m_Mass, m_coreMass, m_Luminosity, m_stellarType);
                
                if(isPerturbation){
                    perturbedLuminosityRadius(m_Luminosity, m_Radius, m_Mu, m_an_coefficients, m_bn_coefficients, m_GBParams, m_massCutoffs, m_timescales, m_Age, m_Mass, m_coreMass, m_HeCoreMass, m_COCoreMass, m_Metallicity, m_logMetallicityXi, m_stellarType, m_randomSeed);
                }
                
                m_Temperature = calculateTemperature(m_Luminosity, m_Radius);
            }
            else{
                if(debugging){
                    std::cout << "If 3) tau > 1 & tau < 0)" << std::endl;
                }
                
                // End of the He MS
                tau           = 1.0;
                
                if(debugging){
                    std::cout << "End of HeMS, start of HeGB" << std::endl;
                }
                
                m_Age         = m_timescales[15]; // t = tHeMS
                //dt          = 0.0;
                //m_Luminosity  = luminosityHeliumMainSequence(tau, m_Mass);
                //m_Radius      = radiusHeliumMainSequence(tau, m_Mass);
                m_Luminosity  = luminosityEndHeliumMainSequence(m_Mass0);
                m_Radius      = radiusEndHeliumMainSequence(m_Mass);
                
                m_HeCoreMass  = m_Mass;
                m_COCoreMass  = 0.0;
                m_coreMass    = m_COCoreMass;
                
                m_Mu          = perturbation_mu(m_Mass, m_coreMass, m_Luminosity, m_stellarType);
                
                if(isPerturbation){
                    perturbedLuminosityRadius(m_Luminosity, m_Radius, m_Mu, m_an_coefficients, m_bn_coefficients, m_GBParams, m_massCutoffs, m_timescales, m_Age, m_Mass, m_coreMass, m_HeCoreMass, m_COCoreMass, m_Metallicity, m_logMetallicityXi, m_stellarType, m_randomSeed);
                }
                
                m_Temperature = calculateTemperature(m_Luminosity, m_Radius);
                
                // How do you know if HeGB or HeHG? Comes from whether using R1, R2, should be set by radiusHeliumGiantBranch
                m_stellarType = NAKED_HELIUM_STAR_HETZSPRUNG_GAP;
                
                // DEBUGGING -- what should CO core mass be at this point?
                // should the core mass shrink with mass loss, well yeah since no longer core, is just the mass of the star
                if(debugging){
                    std::cout << "Core Mass, He Core Mass, CO Core Mass = " << m_coreMass  << " " << m_HeCoreMass << " " <<  m_COCoreMass << std::endl;
                }
                
                // Do I need to recalculate timescales here?
                
            }
            
            // work out final stages of evolution and conditions for evolving to WD i.e. if due to mass loss mc > m// FINISH ME
            if(m_COCoreMass > m_Mass){// EnvelopeLoss
                if(debugging){
                    std::cout << "If 4) Envelope function" << std::endl;
                }
                // Testing envelope loss function
                /*
                m_stellarType = CARBON_OXYGEN_WHITE_DWARF;
                
                m_COCoreMass  = m_Mass;
                m_coreMass    = m_Mass;
                m_Mass        = m_Mass;
                m_Mass0       = m_Mass;
                m_envMass     = 0.0;
                m_Age         = 0.0;
                */
                modifyStarAfterLosingEnvelope(m_stellarType, m_Mass);
                
            }
            
            break;
            
        case NAKED_HELIUM_STAR_HETZSPRUNG_GAP:
            // Do Naked HeHG evolution
            
//            debugging = true;
            
            m_COCoreMass  = coreMassHeliumGiantBranch(m_GBParams, m_timescales, m_Mass0, m_Age);
            m_coreMass    = m_COCoreMass;
            m_Luminosity  = luminosityHeliumGiantBranch(m_GBParams, m_Mass0, m_COCoreMass);
            m_Radius      = radiusHeliumGiantBranch(m_Mass, m_Luminosity, m_stellarType);
            
            m_tau = 0.0;
            
            m_envMass     = envelopeMass(m_Mass, m_coreMass, m_tau, m_stellarType);
            
            // DEBUGGING -- somehow getting negative CO core mass. :( check t > tHeMS
            if(debugging){
                std::cout << "age, tHeMS = " << m_Age << " " << m_timescales[15] << std::endl;
                std::cout << "Mc, L, R (HeHG)= " << m_COCoreMass << " " << log10(m_Luminosity) << " " << log10(m_Radius) << std::endl;
                std::cout << "mass, type (HeHG) = " << m_Mass << " " << m_stellarType << std::endl;
            }
            
            m_Mu          = perturbation_mu(m_Mass, m_coreMass, m_Luminosity, m_stellarType);
            
            if(debugging){std::cout << "HeHG Mu = " << m_Mu << std::endl;}
            
            if(isPerturbation){
                perturbedLuminosityRadius(m_Luminosity, m_Radius, m_Mu, m_an_coefficients, m_bn_coefficients, m_GBParams, m_massCutoffs, m_timescales, m_Age, m_Mass, m_coreMass, m_HeCoreMass, m_COCoreMass, m_Metallicity, m_logMetallicityXi, m_stellarType, m_randomSeed);
            }
            
            // DEBUGGING radius
            if(debugging){
                std:: cout << "AP, log10(R) = " << log10(m_Radius) << std::endl;
            }
            
            m_Temperature = calculateTemperature(m_Luminosity, m_Radius);
            
            // work out final stages of evolution and conditions for evolving to WD
            McMax = maximumCoreMass(m_Mass);
            McSN  = maximumCoreMassSN(m_GBParams);
            
            // DEBUGGING -- check core mass
            if(debugging){
                std::cout << "HeHG core mass = " << m_COCoreMass << " " << m_coreMass << std::endl;
                std::cout << "McMax, McSN = " << McMax << " " << McSN << std::endl;
            }
            
            // Evolves to COWD when MC > MCMAX OR SN WHEN MC > MCSN NEAR EQUATION 89
            if(m_COCoreMass > McMax and McMax < McSN){
                
                if(debugging){
                    
                    std::cout << "Naked He HG star evolving to COWD" << std::endl;
                    
                }
                
                m_stellarType = CARBON_OXYGEN_WHITE_DWARF;
                
            }
            else if(m_COCoreMass > McSN){
                
                // debugging = true;

                if(debugging){
                    std::cout << "HeHG SN" << std::endl;
                }

                // Supernova, classify by M, initial Helium star mass
                supernovaMain(&m_totalMassAtCompactObjectFormation, &m_HeCoreMassAtCompactObjectFormation,
                   &m_COCoreMassAtCompactObjectFormation,&m_coreMassAtCompactObjectFormation, &m_Mass0,
                   &m_Age, &m_envMass, &m_coreMass, &m_HeCoreMass, &m_COCoreMass, &m_coreRadius, &m_Mass, m_MZAMS, m_GBParams[9],
                   &m_fallback, &m_stellarType, &m_Radius, &m_Luminosity, &m_Temperature, &flagSN,
                   &flagECSN, &flagExperiencedECSN, &flagExperiencedCCSN, &flagExperiencedPISN, &flagExperiencedPPISN, &m_drawnKickVelocity, &m_kickVelocity, 
                   &m_error,  &m_flagHrichSN, &m_flagHpoorSN, m_SNengine, m_randomSeed, &m_angularMomentum, &m_pulsarSpinPeriod, &m_pulsarSpinFrequency, &m_pulsarMagneticField, &m_pulsarSpinDownRate, &m_momentInertia, options, r);
                
            }
            else if(m_COCoreMass > m_Mass){
                modifyStarAfterLosingEnvelope(m_stellarType, m_Mass);
            }
    
            
            break;
            
        case NAKED_HELIUM_STAR_GIANT_BRANCH:
            
//            debugging = true;
            
            // Do Naked HeGB evolution
            
            m_COCoreMass  = coreMassHeliumGiantBranch(m_GBParams, m_timescales, m_Mass0, m_Age);
            m_coreMass    = m_COCoreMass;
            m_Luminosity  = luminosityHeliumGiantBranch(m_GBParams, m_Mass0, m_COCoreMass);
            m_Radius      = radiusHeliumGiantBranch(m_Mass, m_Luminosity, m_stellarType);
            
            m_Mu          = perturbation_mu(m_Mass, m_coreMass, m_Luminosity, m_stellarType);
            
            m_envMass     = envelopeMass(m_Mass, m_coreMass, m_tau, m_stellarType);
            
            if(isPerturbation){
                perturbedLuminosityRadius(m_Luminosity, m_Radius, m_Mu, m_an_coefficients, m_bn_coefficients, m_GBParams, m_massCutoffs, m_timescales, m_Age, m_Mass, m_coreMass, m_HeCoreMass, m_COCoreMass, m_Metallicity, m_logMetallicityXi, m_stellarType, m_randomSeed);
            }
            
            m_Temperature = calculateTemperature(m_Luminosity, m_Radius);
            
            // Work out final stages of evolution and conditions for evolving to WD // FINISH ME
            McMax = maximumCoreMass(m_Mass);
            McSN  = maximumCoreMassSN(m_GBParams);
            
            // DEBUGGING -- check core mass
            if(debugging){
                std::cout << "HeGB core mass = " << m_COCoreMass << " " << m_coreMass << std::endl;
                std::cout << "McMax, McSN = " << McMax << " " << McSN << std::endl;
                std::cout << "mass, type (HeHG) = " << m_Mass << " " << m_stellarType << std::endl;
            }
            
            // Evolves to COWD when MC > MCMAX OR SN WHEN MC > MCSN NEAR EQUATION 89
            if(m_COCoreMass > McMax and McMax < McSN){
                
                if(debugging){
                    
                    std::cout << "Naked He GB star evolving to COWD" << std::endl;
                    
                }
                
                m_stellarType = CARBON_OXYGEN_WHITE_DWARF;
                
            }
            // Evolves to COWD when MC > MCMAX OR SN WHEN MC > MCSN NEAR EQUATION 89
            else if(m_COCoreMass > McSN){

                // debugging = true;

                if(debugging){
                    std::cout << "HeGB SN" << std::endl;
                }

                supernovaMain(&m_totalMassAtCompactObjectFormation, &m_HeCoreMassAtCompactObjectFormation,
                   &m_COCoreMassAtCompactObjectFormation,&m_coreMassAtCompactObjectFormation, &m_Mass0,
                   &m_Age, &m_envMass, &m_coreMass, &m_HeCoreMass, &m_COCoreMass, &m_coreRadius, &m_Mass, m_MZAMS, m_GBParams[9],
                   &m_fallback, &m_stellarType, &m_Radius, &m_Luminosity, &m_Temperature, &flagSN,
                   &flagECSN, &flagExperiencedECSN, &flagExperiencedCCSN, &flagExperiencedPISN, &flagExperiencedPPISN, &m_drawnKickVelocity, &m_kickVelocity, 
                   &m_error,  &m_flagHrichSN, &m_flagHpoorSN, m_SNengine, m_randomSeed, &m_angularMomentum, &m_pulsarSpinPeriod, &m_pulsarSpinFrequency, &m_pulsarMagneticField, &m_pulsarSpinDownRate, &m_momentInertia, options, r);
            }


            else if(m_COCoreMass > m_Mass){
                modifyStarAfterLosingEnvelope(m_stellarType, m_Mass);
            }
            // Not sure this should lead to an error
            //                else{
            //
            //                    std::cout << "Error in NAKED_HELIUM_STAR_GIANT_BRANCH evolution" << std::endl;
            //
            //                }
            
            
            break;
            
        case HELIUM_WHITE_DWARF:
            
            // Do HeWD evolution (only expected in binaries where GB star loses envelope)
            
            m_Luminosity    = luminosityWhiteDwarf(m_Age, m_Mass, m_Metallicity, m_stellarType);
            m_Radius        = radiusWhiteDwarf(m_Mass);
            m_Temperature   = calculateTemperature(m_Luminosity, m_Radius);
            m_coreMass      = m_Mass;
            m_envMass       = 0.0;
            
            break;
            
        case CARBON_OXYGEN_WHITE_DWARF:
            
            // Do COWD evolution
            // Expected due to envelope loss of TPAGB star with M < Mup
            
            m_Luminosity    = luminosityWhiteDwarf(m_Age, m_Mass, m_Metallicity, m_stellarType);
            m_Radius        = radiusWhiteDwarf(m_Mass);
            m_Temperature   = calculateTemperature(m_Luminosity, m_Radius);
            m_coreMass      = m_Mass;
            m_envMass       = 0.0;
            
            if(m_Mass > Mch){
                
                // Supernova leaving a massless remnant
                
                if(debugging){
                    std::cout << "Massless remnant supernova from COWD 1)" << std::endl;
                }

                m_Mass = 0.0;
                m_Radius = 0.0;
                m_Luminosity = 0.0;
                m_Age = 0.0;
                m_stellarType = MASSLESS_REMNANT;
            }
            
            break;
            
        case OXYGEN_NEON_WHITE_DWARF:
            
            // Do ONeWD evolution
            // Expected from loss of envelope of TPAGB star with Mup <= M <= Mec
            
            m_Luminosity    = luminosityWhiteDwarf(m_Age, m_Mass, m_Metallicity, m_stellarType);
            m_Radius        = radiusWhiteDwarf(m_Mass);
            m_Temperature   = calculateTemperature(m_Luminosity, m_Radius);
            m_coreMass      = m_Mass;
            m_envMass       = 0.0;
            
            if(m_Mass > Mch){
                supernovaMain(&m_totalMassAtCompactObjectFormation, &m_HeCoreMassAtCompactObjectFormation,
                   &m_COCoreMassAtCompactObjectFormation,&m_coreMassAtCompactObjectFormation, &m_Mass0,
                   &m_Age, &m_envMass, &m_coreMass, &m_HeCoreMass, &m_COCoreMass, &m_coreRadius, &m_Mass, m_MZAMS, m_GBParams[9],
                   &m_fallback, &m_stellarType, &m_Radius, &m_Luminosity, &m_Temperature, &flagSN,
                   &flagECSN, &flagExperiencedECSN, &flagExperiencedCCSN, &flagExperiencedPISN, &flagExperiencedPPISN, &m_drawnKickVelocity, &m_kickVelocity, 
                   &m_error,  &m_flagHrichSN, &m_flagHpoorSN, m_SNengine, m_randomSeed, &m_angularMomentum, &m_pulsarSpinPeriod, &m_pulsarSpinFrequency, &m_pulsarMagneticField, &m_pulsarSpinDownRate, &m_momentInertia, options, r);
            }
            
            break;
            
        case NEUTRON_STAR:
            
            // Do NS evolution
            
            m_coreMass    = m_Mass;
            m_COCoreMass  = m_Mass;
            m_HeCoreMass  = m_Mass;
            m_Luminosity  = luminosityNeutronStar(m_Age, m_Mass);
            //m_Radius      = radiusNeutronStar(); // Radius of NS set to 10km
            m_Radius      = neutronStarRadius(m_Mass, options) * km_in_rsol; // neutronStarRadius is in km
            m_Temperature = calculateTemperature(m_Luminosity, m_Radius);
            
            break;
            
        case BLACK_HOLE:
            
            // Do Black hole evolution
            
            m_coreMass    = m_Mass;
            m_COCoreMass  = m_Mass;
            m_HeCoreMass  = m_Mass;
            m_Luminosity  = luminosityBlackHole(m_Mass);
            m_Radius      = SchwarzschildRadiusBlackHole(m_Mass); // Schwarzschild radius (not correct for rotating BH)
            m_Temperature = calculateTemperature(m_Luminosity, m_Radius);
            
            break;
            
        case MASSLESS_REMNANT:
            
            // Massless remnant. Stellar type : 15
            
            if(debugging){
                std::cout << "Star is a massless remnant" << std::endl;
            }
            
            break;
                
        default:
            
            
            std::cerr << m_randomSeed << "\tError in evolveOneTimestep function. Unknown stellar type : " << m_stellarType << std::endl;
            
            if(m_error == false){
                m_error = true;
            }
            
            break;
            
    }


//    massLossCaller(m_Mass, m_Mass0, m_coreMass, m_Luminosity, m_Radius, m_Mu, m_Metallicity, m_Temperature, m_logMetallicityXi, m_Mdot, m_dt, m_Age, m_stellarType, isMassLoss, useFakeMassLoss, massLossPrescription, debugging, m_an_coefficients, m_timescales);
    
    // Make sure to check this doesn't modify anything
    if(useFakeMassLoss){
        massLossPrescriptionFake(m_Mass, m_Mass0, m_coreMass, m_Luminosity, m_Radius, m_Mu, m_Metallicity, m_Temperature, m_logMetallicityXi, m_Mdot, m_dt, m_Age, m_stellarType, isMassLoss, massLossPrescription, debugging, m_an_coefficients, m_timescales);
    }

    // Could you just recalculate timescales here? Don't really want to have to do this twice per loop, shouldn't need to
    calculateGiantBranchParameters(m_GBParams, m_an_coefficients, m_bn_coefficients, m_massCutoffs, m_Mass0, m_logMetallicityXi, m_stellarType, m_randomSeed);
    calculateTimescales(m_timescales, m_GBParams, m_an_coefficients, m_bn_coefficients, m_massCutoffs, m_Mass0, m_coreMass, m_COCoreMass, m_Metallicity, m_logMetallicityXi, m_stellarType, m_randomSeed);
    // I guess at least then you would know that this was the problem -- yup
 
/* 28/10/2015 -- I think this lines are not needed anymore, dt and m_time, m_Age and  m_dt are calculated in another function. I still dont know about m_tau and m_MassPrev, though.
    // Choose suitable timestep -- shouldn't this be at the begining? I still don't know?
    dt = chooseTimestep(m_stellarType, m_Age, m_timescales); // I guess this will need to be settable or a property of the star
    
    // Update both the physical simulation time and the age of the star
    m_time    += dt;
    m_Age     += dt;
    m_tau      = tau;
    m_dt       = dt;
 */
//    m_Age     += dt;
    m_tau      = tau;
	
	// Calculate timescales
	m_dynamicalTimescale	= dynamicalTimescale(m_Mass, m_Radius);
	m_thermalTimescale 		= thermalTimescale(m_Mass, m_envMass, m_Radius, m_Luminosity, m_stellarType); 
	m_nuclearTimescale		= nuclearTimescale(m_Mass, m_Luminosity);

    // ALEJANDRO - 03/09/2019 - If statement to avoid recalculating m_radialExpansionTimescale when m_dt = 0.0. Approximated to m_dt >= m_dynamicalTimescale.
    if (m_dt >= m_dynamicalTimescale) {
        m_radialExpansionTimescale = calculateRadialExpansionTimescale(m_Radius,m_RadiusPrev,m_stellarTypePrev,m_stellarType,dt,m_time,m_dynamicalTimescale);
    }


}

/*
 MASS LOSS FUNCTIONS
 */
// Changed Mass and Age from pointer to value, because they will be updated outside winds. Is this a problem? Everything else should change according to winds, only this quantities have to do with Mass Transfer.
double  Star::massLossCaller(double Mass, double Mass0, double coreMass, double & Luminosity, double & Radius, double & Mu, double & Metallicity, double & Teff, double logMetallicityXi, double & Mdot, double & dt, double Age, double flbv, double wolfRayetFactor, int & stellarType, bool isMassLoss, bool useFakeMassLoss, int massLossPrescription, bool debugging, double an_coefficients[nrowsa], double timescales[nTimescales]){


    /*
     Pass parameters to relevant mass loss function
     
     Parameters
     -----------
     massLossPrescription : int
        Which mass loss prescription to use:
        Hurley et al 2000 mass loss prescription
        Vink et al 2001 mass loss prescription
     isMassLoss : bool
        Whether to include mass loss
     debugging : bool
        Whether we are debugging -- triggers more print statements.
     
     Returns
     --------
     
     */


    // Some debugging statements
    debugging = false;
    if(debugging){
        std::cout << "Using " << massLossPrescription << " mass loss prescription." << std::endl;
        
        if(isMassLoss){
            std::cout << "Mass loss is enabled" << std::endl;
        }
        else{
            std::cout << "Mass loss is disabled" << std::endl;
        }
    }
    // AVG: Modifiy so I can get the amount of mass loss/accreted from MT into this function. Calculate the amount mass due to winds, which will be returned by the function for modyfing the orbit, and then modify the mass loss rate within the function so it matches the mass loss with mass transfer. In case of turning off mass transfer, it should do the same as regular winds.
    
    if(isMassLoss){
        
        double timePrime=0;
        
        // Apply Mass Loss -- should this be a property of the star, probably.
        if(massLossPrescription == MASS_LOSS_PRESCRIPTION_HURLEY){
            Mdot = massLossRate(Mass, Luminosity, Radius, Mu, Metallicity, stellarType);         // Note also that mass loss rate given per year, times are in Myrs so will need to multiply by 10^6
        }
        else if (massLossPrescription == MASS_LOSS_PRESCRIPTION_VINK){
            Mdot = massLossRateVink(Mass, Luminosity, Radius, Mu, Metallicity, Teff, flbv, wolfRayetFactor, stellarType, m_LBVphaseFlag);   // Note also that mass loss rate given per year, times are in Myrs so will need to multiply by 10^6 when subtracting lost mass
        }
        else{
            std::cerr << m_randomSeed << "\tError in mass loss caller. No mass prescription selected.\n";
        }
        
        if(debugging){
            std::cout << "stellar type = " << stellarType << std::endl;
            std::cout << "Mdot, dt, Mdot * dt * 1E6 = " << Mdot << " " <<  dt << " " << Mdot * dt * 1E6 << std::endl;
            std::cout << "Before mass loss, M = " << Mass << "Msol" << std::endl;
            std::cout << "Before mass loss, Age = " << Age << "Myr" << std::endl;
        }
        // Limit mass loss rate to 1%, warn me otherwise -- is this accurate enough, or do you need to use trapezium etc?
        if(Mdot * dt * 1E6 < 0.01 * Mass){
                Mass -= Mdot * dt * 1E6;
        }
        else{

            if(debugging){
                std::cout << "Mass loss >1% = " << Mdot * dt * 1E6 << " in " << dt << " Myrs" << std::endl;
                std::cout << "Mass loss limited to " << 0.01 * Mass << std::endl;
            }
            
            // reset Mdot, dt to limit mass loss to 1% and then update this stars variables.
            dt = (0.01 * Mass) / (Mdot * 1E6); // Change timestep to time over which mass loss is 1% (assuming Mdot constant)
            
            Mdot = 0.01 * Mass / dt; // limited mass loss in Msol per year
            
            if(debugging){
                std::cout << "Mdot, dt limited to " << Mdot << " in " << dt << " Myrs" << std::endl;
            }
            
                Mass *= 0.99; // limit mass loss to 1% and update mass
        }
        
        // Mass loss should only effect the age of the star, not the time of the simulation
//        timePrime = timeChangeAfterMassLoss(an_coefficients, Mass0, Mass, coreMass, Age, timescales, stellarType, logMetallicityXi);
        
        if(debugging){
            std::cout << "Time before mass loss = " << Age << std::endl;
            std::cout << "After mass loss, M = " << Mass << " after losing " << Mdot << "Msol" << std::endl;
//            std::cout << "After mass loss, Age = " << timePrime << "Myr\n" << std::endl;
        }
//        Age = timePrime; // Update age due to mass loss.
        
        return  Mass;
    }
    else{
        std::cerr << m_randomSeed << "\tError in mass loss caller: Specify mass loss prescription\n";
        return Mass;
    }
}

//double Star::timeChangeAfterMassLoss(double an_coefficients[nrowsa], double &Mass0, double Mass, double coreMass, double time, double timescales[nTimescales], int STELLAR_TYPE, double logMetallicityXi){
//    /*
//     Calculate the new time for the new mass of the star after mass loss by assuming its fractional age remains constant
//     
//     Given towards the end of section 7.1 in Hurley et al 2000
//     
//     Parameters
//     -----------
//     an_coefficients : double array
//     Array containing metallicity dependent an coefficients
//     &Mass0 : double
//     Initial Mass in Msol
//     Mass : double
//     Mass in Msol after mass loss
//     coreMass : double
//     Core mass in Msol
//     time : double
//     Time/age prior to mass loss
//     timescales : double array
//     Array containing timescales for the star
//     STELLAR_TYPE : int
//     Stellar evolution type
//     logMetallicityXi : double
//     logMetallicityXi
//     
//     timescales
//     -----------
//     0 : tMS  - Main sequence
//     1 : tBGB - Base Giant Branch
//     2 : tHeI - Helium ignition
//     3 : tHe  - Helium burning
//     #-- First Giant Branch (FGB)
//     4 : Giant tinf1
//     5 : Giant tinf2
//     6 : Giant t(Mx)
//     #-- Early Asymptotic Giant Branch (EAGB)
//     7 : FAGB tinf1
//     8 : FAGB tinf2
//     9 : FAGB t(Mx)
//     #-- Thermally Pulsating Asymptotic Giant Branch (TPAGB)
//     10 : SAGB tinf1
//     11 : SAGB tinf2
//     12 : SAGB t(Mx)
//     13 : TP (tDU?)
//     14 : T(McMax)
//     #-- not sure what TP, t(McMax) are
//     #-- what about tHeMS etc?
//     15 : tHeMS
//     
//     Returns
//     --------
//     timePrime : double
//     Mass loss adjusted time in Myrs
//     
//     
//     */
//    // Calculate Mass and metallicity dependent mu (Equation 7)
//    double mu        = std::max(0.5, 1.0 - 0.01 * std::max((an_coefficients[6 - 1]/pow(Mass, an_coefficients[7 - 1])), (an_coefficients[8 - 1] + an_coefficients[9 - 1]/pow(Mass, an_coefficients[10 - 1]))));
//    
//    // Timescales without mass loss
//    double tMS      = timescales[0];
//    double tBGB     = timescales[1];
//    //double thook    = mu * tBGB;
//    double tHeMS    = timescales[15];
//    
//    // Declare variable for adjusted time/age
//    double timePrime = 0;
//    
//    // Calculate primed timescale variables
//    if(debugging)
//        std::cout << "timechangemassloss" << std::endl;
//    
//    // Check functions: lifetimeToBGB, lifetimeMainSequence, lifetimeHeliumMainSequence
//    double tBGBPrime  = lifetimeToBGB(an_coefficients, Mass);
//    double thookPrime = mu * tBGBPrime;
//    double tMSPrime   = lifetimeMainSequence(an_coefficients, Mass, logMetallicityXi, tBGBPrime, thookPrime);
//    double tHeMSPrime = lifetimeHeliumMainSequence(Mass);
//    
//    if(debugging){
//        std::cout << "\nMass: \t" << Mass << std::endl;
//        std::cout << "mu: \t" << mu << std::endl;
//        std::cout << "tBGB: \t" << tBGB << std::endl;
//        std::cout << "tBGBPrime: \t" << tBGBPrime << std::endl;
//        std::cout << "thookPrime: \t" << thookPrime << std::endl;
//        std::cout << "tMS: \t" << tMS << std::endl;
//        std::cout << "tMSPrime: \t" << tMSPrime << std::endl;
//    }
//    if(STELLAR_TYPE == MS_LESS_THAN_07 or STELLAR_TYPE == MS_MORE_THAN_07){
//        
//        // On the main sequence, M0 = Mt and the star is effectively aged by: t' = (tMS'/tMS) * t
//        Mass0 = Mass;
//        
//        if(debugging){
//            std::cout << "tMS: \t" << tMS << std::endl;
//            std::cout << "tMSPrime: \t" << tMSPrime << std::endl;
//        }
//        
//        timePrime = (tMSPrime / tMS) * time; // Given towards the end of section 7.1
//        
//    }
//    else if(STELLAR_TYPE == NAKED_HELIUM_STAR_MS){
//        
//        // Same for Naked Helium Main Sequence stars
//        Mass0 = Mass;
//        if(debugging){
//            std::cout << "tMS: \t" << tMS << std::endl;
//            std::cout << "tMSPrime: \t" << tMSPrime << std::endl;
//        }
//        
//        timePrime = (tHeMSPrime / tHeMS) * time; // Given towards the end of section 7.1
//        
//    }
//    // On the giant branches, age determines core mass, unaffected by change in envelope mass, so do not need to alter age/inital mass
//    else if(STELLAR_TYPE == HERTZSPRUNG_GAP){
//        
//        // M0 = Mt so long as M0 > Mc so as to avoid unphysical decrease in core mass.
//        if(Mass0 > coreMass){
//            if(debugging){
//                std::cout << "tMS: \t" << tMS << std::endl;
//                std::cout << "tMSPrime: \t" << tMSPrime << std::endl;
//            }
//            Mass0 = Mass;
//        }
//        else{
//            if(debugging){
//                std::cout << "tMS: \t" << tMS << std::endl;
//                std::cout << "tMSPrime: \t" << tMSPrime << std::endl;
//            }
//            Mass0 = Mass0;
//        }
//        
//        timePrime = tMSPrime + ((tBGBPrime - tMSPrime)/(tBGB - tMS)) * (time - tMS); // Given towards the end of section 7.1
//        
//    }
//    // Not sure if this is correct for NAKED HELIUM HG
//    //else if(STELLAR_TYPE == NAKED_HELIUM_STAR_HETZSPRUNG_GAP){
//    //
//    //    timePrime = tHeMSPrime + ((tBGBPrime - tHeMSPrime)/(tBGB - tHeMS)) * (time - tHeMS); // Given towards the end of section 7.1
//    //
//    //
//    //}
//    else{
//        
//        timePrime = time;
//        
//    }
//    // WHAT ABOUT FOR REMNANTS? NOT REALLY ANY MASS LOSS (MAY BE MASS GAIN IF ACCRETION) SO DO M0 = M ELSEWHERE FINISH ME
//    
//    
//    return timePrime;
//    
//}

double timeChangeAfterMassLoss(double an_coefficients[nrowsa], double &Mass0, double Mass, double coreMass, double age, double timescales[nTimescales], int STELLAR_TYPE, double logMetallicityXi, bool debugging){
    /*
     Calculate the new time for the new mass of the star after mass loss by assuming its fractional age remains constant
     
     Given towards the end of section 7.1 in Hurley et al 2000
     
     Parameters
     -----------
     an_coefficients : double array
     Array containing metallicity dependent an coefficients
     &Mass0 : double
     Initial Mass in Msol
     Mass : double
     Mass in Msol after mass loss
     coreMass : double
     Core mass in Msol
     age : double
        age prior to mass loss
     timescales : double array
     Array containing timescales for the star
     STELLAR_TYPE : int
     Stellar evolution type
     logMetallicityXi : double
     logMetallicityXi
     
     timescales
     -----------
     0 : tMS  - Main sequence
     1 : tBGB - Base Giant Branch
     2 : tHeI - Helium ignition
     3 : tHe  - Helium burning
     #-- First Giant Branch (FGB)
     4 : Giant tinf1
     5 : Giant tinf2
     6 : Giant t(Mx)
     #-- Early Asymptotic Giant Branch (EAGB)
     7 : FAGB tinf1
     8 : FAGB tinf2
     9 : FAGB t(Mx)
     #-- Thermally Pulsating Asymptotic Giant Branch (TPAGB)
     10 : SAGB tinf1
     11 : SAGB tinf2
     12 : SAGB t(Mx)
     13 : TP (tDU?)
     14 : T(McMax)
     #-- not sure what TP, t(McMax) are
     #-- what about tHeMS etc?
     15 : tHeMS
     
     Returns
     --------
     timePrime : double
     Mass loss adjusted time in Myrs
     
     
     */
    // Calculate Mass and metallicity dependent mu (Equation 7)
    double mu        = std::max(0.5, 1.0 - 0.01 * std::max((an_coefficients[6 - 1]/pow(Mass, an_coefficients[7 - 1])), (an_coefficients[8 - 1] + an_coefficients[9 - 1]/pow(Mass, an_coefficients[10 - 1]))));
    
    // Timescales without mass loss
    double tMS      = timescales[0];
    double tBGB     = timescales[1];
    //double thook    = mu * tBGB;
    double tHeMS    = timescales[15];
    
    // Declare variable for adjusted time/age
    double agePrime = 0;
    
    // Calculate primed timescale variables
    debugging = false;
    
    // Check functions: lifetimeToBGB, lifetimeMainSequence, lifetimeHeliumMainSequence
    double tBGBPrime  = lifetimeToBGB(an_coefficients, Mass);
    double thookPrime = mu * tBGBPrime;
    double tMSPrime   = lifetimeMainSequence(an_coefficients, Mass, logMetallicityXi, tBGBPrime, thookPrime);
    double tHeMSPrime = lifetimeHeliumMainSequence(Mass);
    if(debugging){
        std::cout << "timechangemassloss" << std::endl;
        std::cout << "\nMass: \t" << Mass << std::endl;
        std::cout << "mu: \t" << mu << std::endl;
        std::cout << "tBGB: \t" << tBGB << std::endl;
        std::cout << "tBGBPrime: \t" << tBGBPrime << std::endl;
        std::cout << "thookPrime: \t" << thookPrime << std::endl;
        std::cout << "tMS: \t" << tMS << std::endl;
        std::cout << "tMSPrime: \t" << tMSPrime << std::endl;
    }
    
    if(STELLAR_TYPE == MS_LESS_THAN_07 or STELLAR_TYPE == MS_MORE_THAN_07){
        
        // On the main sequence, M0 = Mt and the star is effectively aged by: t' = (tMS'/tMS) * t
        Mass0 = Mass;
        
        agePrime = (tMSPrime / tMS) * age; // Given towards the end of section 7.1
        
        if(debugging){
            std::cout << "tMS: \t" << tMS << std::endl;
            std::cout << "tMSPrime: \t" << tMSPrime << std::endl;
            std::cout << "time = \t" << age << std::endl;
            std::cout << "timePrime = \t" << agePrime << std::endl;
        }
        
        
    }
    else if(STELLAR_TYPE == NAKED_HELIUM_STAR_MS){
        
        // Same for Naked Helium Main Sequence stars
        Mass0 = Mass;
        
        if(debugging){
            std::cout << "tMS: \t" << tMS << std::endl;
            std::cout << "tMSPrime: \t" << tMSPrime << std::endl;
        }
        
        agePrime = (tHeMSPrime / tHeMS) * age; // Given towards the end of section 7.1
        
    }
    // On the giant branches, age determines core mass, unaffected by change in envelope mass, so do not need to alter age/inital mass
    else if(STELLAR_TYPE == HERTZSPRUNG_GAP){
        
        // M0 = Mt so long as M0 > Mc so as to avoid unphysical decrease in core mass.
        if(Mass0 > coreMass){
            if(debugging){
                std::cout << "tMS: \t" << tMS << std::endl;
                std::cout << "tMSPrime: \t" << tMSPrime << std::endl;
            }
            Mass0 = Mass;
        }
        else{
            if(debugging){
                std::cout << "tMS: \t" << tMS << std::endl;
                std::cout << "tMSPrime: \t" << tMSPrime << std::endl;
            }
            Mass0 = Mass0;
        }
        
        agePrime = tMSPrime + ((tBGBPrime - tMSPrime)/(tBGB - tMS)) * (age - tMS); // Given towards the end of section 7.1
        
    }
    // Not sure if this is correct for NAKED HELIUM HG
    //else if(STELLAR_TYPE == NAKED_HELIUM_STAR_HETZSPRUNG_GAP){
    //
    //    timePrime = tHeMSPrime + ((tBGBPrime - tHeMSPrime)/(tBGB - tHeMS)) * (time - tHeMS); // Given towards the end of section 7.1
    //
    //
    //}
    else{
        
        agePrime = age;
        
    }
    // WHAT ABOUT FOR REMNANTS? NOT REALLY ANY MASS LOSS (MAY BE MASS GAIN IF ACCRETION) SO DO M0 = M ELSEWHERE FINISH ME
    
    if(debugging){std::cout<< "age, agePrime: " << age << " " << agePrime << std::endl;}
    
    return agePrime;
    
}

void ageChangeAfterMassChange(Star & star, const programOptions & options){
    /*
     Update the age of a star after its mass changes, either due to mass loss through winds or mass loss or gain through mass transfer
     
     Parameters
     -----------
     star : Star
        Star to update age of
     options : programOptions
     
     Returns
     --------
     
     
     */
    
    // rejuvination factor
    double  f_rej = 1.0;
    bool    debugging = false;
    // debugging
    if(debugging){
        std::cout << "Mprev = " << star.m_MassPrev << std::endl;
        std::cout << "M = " << star.m_Mass << std::endl;
    }
    if(options.massTransferRejuvenationPrescription == MASS_TRANSFER_REJUVENATION_NONE){
        // Hurley et al 2000 prescription
        f_rej = 1.0;
    }
    else if (options.massTransferRejuvenationPrescription == MASS_TRANSFER_REJUVENATION_STARTRACK){
        // StarTrack 2008 prescription
        // See section 5.6 of http://arxiv.org/pdf/astro-ph/0511811v3.pdf
        
        if(star.m_Mass <= star.m_MassPrev){
            // Rejuvenation factor is unity for mass losing stars
            f_rej = 1.0;
        }
        else if (star.m_stellarType == MS_LESS_THAN_07){
            // Rejuvenation factor is unity for convective main sequence stars
            f_rej = 1.0;
        }
        else if (star.m_stellarType == MS_MORE_THAN_07 or star.m_stellarType == NAKED_HELIUM_STAR_MS){
            f_rej = star.m_MassPrev / star.m_Mass;
        }
        else if (star.m_stellarType == HERTZSPRUNG_GAP){
            f_rej = 1.0;
        }
        else{
            std::cout << "Shouldn't get here 1) ageChangeAfterMassChange" << std::endl;
        }
        
    }
    else{
        // Shouldn't get here
        std::cout << "Shouldn't get here 2) ageChangeAfterMassChange" << std::endl;
        f_rej = 1.0;
    }
    
    double agePrime = timeChangeAfterMassLoss(star.m_an_coefficients, star.m_Mass0, star.m_Mass, star.m_coreMass, star.m_Age, star.m_timescales, star.m_stellarType, star.m_logMetallicityXi, false);
    
    // Debugging
    if(debugging){
        std::cout << "f_rej = " << f_rej << std::endl;
        std::cout << "age = " << star.m_Age << std::endl;
        std::cout << "agePrime = " << agePrime << std::endl;
    }
    star.m_Age = agePrime * f_rej;
    
}

double massTransferRejuvenationFactor(Star & star, const programOptions & options){
    /*
     Calculate rejuvenation factor for stellar age based on mass lost/gained during mass transfer
     
     Parameters
     ----------
     star : Star
        Star to calculate rejuvenation factor for
     options : programOptions
        User specified program options
     
     Returns
     --------
     f_rej : double
        Rejuvenation factor => 1.0
     
     
    */
    // rejuvination factor
    double  f_rej = 1.0;
    bool    debugging = false;
    // debugging
    if(debugging){
    std::cout << "Mprev = " << star.m_MassPrev << std::endl;
    std::cout << "M = " << star.m_Mass << std::endl;
    }
    
    if(options.massTransferRejuvenationPrescription == MASS_TRANSFER_REJUVENATION_NONE){
        // Hurley et al 2000 prescription
        f_rej = 1.0;
    }
    else if (options.massTransferRejuvenationPrescription == MASS_TRANSFER_REJUVENATION_STARTRACK){
        // StarTrack 2008 prescription
        // See section 5.6 of http://arxiv.org/pdf/astro-ph/0511811v3.pdf
        
        if(star.m_Mass <= star.m_MassPrev){
            // Rejuvenation factor is unity for mass losing stars
            f_rej = 1.0;
        }
        else if (star.m_stellarType == MS_LESS_THAN_07){
            // Rejuvenation factor is unity for convective main sequence stars
            f_rej = 1.0;
        }
        else if (star.m_stellarType == MS_MORE_THAN_07 or star.m_stellarType == NAKED_HELIUM_STAR_MS){
            f_rej = star.m_MassPrev / star.m_Mass;
        }
        else if (star.m_stellarType == HERTZSPRUNG_GAP){
            f_rej = 1.0;
        }
        else{
            // Alejandro: I think there shouldn't be a problem for getting here. Commented error message, delete after discussing.
            // Case accretor into WD's still missing, though. Should f_rej = 1.0 here?
            // std::cout << "Shouldn't get here ageChangeAfterMassChange" << std::endl;
            if(debugging){
                std::cout << "Rejuvenation into GB, CHeB, EAGB, TPAGB, HeHG, HeGB and WDs not implemented. f_rej = 1.0" << std::endl;
            }
            f_rej = 1.0;

        }
        
    }
    else{
        // Shouldn't get here
//        std::cout << "Shouldn't get here 2) ageChangeAfterMassChange" << std::endl;
        std::cerr << star.m_randomSeed << "\tError: No Rejuvenation Prescription chosen." << std::endl;
        f_rej = 1.0;
    }
    
    return f_rej;
}

void Star::massLossPrescriptionFake(double & Mass, double Mass0, double coreMass, double & Luminosity, double & Radius, double & Mu, double & Metallicity, double & Teff, double logMetallicityXi, double & Mdot, double & dt, double Age, int & stellarType, bool isMassLoss, int massLossPrescription, bool debugging, double an_coefficients[nrowsa], double timescales[nTimescales]){
    /*
     Fake mass loss prescription used for determining how star responds to mass loss
     
     Parameters
     -----------
     isMassLoss : bool
     Whether to include mass loss
     debugging : bool
     Whether we are debugging -- triggers more print statements.
     
     Returns
     ----------
     
     
     */
    double timePrime=0;
    
    // Apply Mass Loss -- should this be a property of the star, probably.
    //Mdot = massLossRateVink(Mass, Luminosity, Radius, Mu, Metallicity, Teff, stellarType);   // Note also that mass loss rate given per year, times are in Myrs so will need to multiply by 10^6 when subtracting lost mass
    
    //        if(debugging){
    //            std::cout << "Mdot, dt, Mdot * dt * 1E6 = " << Mdot << " " <<  dt << " " << Mdot * dt * 1E6 << std::endl;
    //        }
    //        // Limit mass loss rate to 1%, warn me otherwise -- is this accurate enough, or do you need to use trapezium etc?
    //        if(Mdot * dt * 1E6 < 0.01 * Mass){
    //            Mass -= Mdot * dt * 1E6;
    //        }
    //        else{
    //
    //            if(debugging){
    //                std::cout << "Mass loss >1% = " << Mdot * dt * 1E6 << " in " << dt << " Myrs" << std::endl;
    //                std::cout << "Mass loss limited to " << 0.01 * Mass << std::endl;
    //            }
    //
    //            // reset Mdot, dt to limit mass loss to 1% and then update this stars variables.
    //            dt = (0.01 * Mass) / (Mdot * 1E6); // Change timestep to time over which mass loss is 1% (assuming Mdot constant)
    //
    //            Mdot = 0.01 * Mass / dt; // limited mass loss in Msol per year
    //
    //            if(debugging){
    //                std::cout << "Mdot, dt limited to " << Mdot << " in " << dt << " Myrs" << std::endl;
    //            }
    //
    //            Mass *= 0.99; // limit mass loss to 1% and update mass
    //
    //        }
    
    if(debugging){
        std::cout << "Mass before mass loss = " << Mass << std::endl;
    }
    
    double percentageMassLoss = 0.01;//0.0;//0.001;              // Use 1/0.1% mass loss and update mass
    double dM = Mass * (percentageMassLoss * 0.01);                      // Change in mass
    Mass *= (1.0 - percentageMassLoss * 0.01);                  // Subtract (could add) this amount of mass and see how star responds (adding avoids some issues)
    
    if(debugging){
        std::cout << "Mass after mass loss = " << Mass << std::endl;
    }
    
    // Mass loss should only effect the age of the star, not the time of the simulation
    timePrime = timeChangeAfterMassLoss(an_coefficients, Mass0, Mass, coreMass, Age, timescales, stellarType, logMetallicityXi, debugging);
    
    if(debugging){
        std::cout << "Time before mass loss = " << Age << std::endl;
        std::cout << "After mass loss, M = " << Mass << " after losing " << dM << "Msol" << std::endl;
    }
    Age = timePrime; // Update age due to mass loss.
        
}

/*
 ALEJANDROS FUNCTIONS
 */

void Star::modifyStarAfterLosingEnvelope(int stellarType, double mass){
    double  tau = 0.0;
    double  McBAGB = 0.0; // It is initialized like this in evolveOneTimestep function. Look for this in detail. AVG.
    bool    debugging = false;
    
    if(mass > 0.0){
        switch (m_stellarType) {
            case MS_LESS_THAN_07:
                
                // At this point, star is only envelope. Therefore, end? massless remnant? error?
                break;
                
            case MS_MORE_THAN_07:
                
                // At this point, star is only envelope. Therefore, end? massless remnant? error?
                break;
                
            case HERTZSPRUNG_GAP:
                
                // Star has lost its envelope, forms a HeWD if has degenerate core (M < MHeF) or zero-age naked helium star
                // Described just above Equation 76
                
                if(m_Mass < m_massCutoffs[1]){
                    
                    // Reset mass/age parameters
                    m_Mass        = m_coreMass;
                    m_Age         = 0;  // Check eq. 76, Hurley et al. 2000
                    m_stellarType = HELIUM_WHITE_DWARF;
                    m_Radius        = radiusWhiteDwarf(m_Mass);
                    
                }
                else{
                    
                    // Zero age Helium MS star // see section 6.1
                    
                    m_Mass    = m_coreMass;
                    m_Mass0       = m_Mass;
                    m_Age         = 0; // Check eq. 76, Hurley et al. 2000
                    m_stellarType = NAKED_HELIUM_STAR_MS;
                    
                    // Something for calculate radius and luminosity here?
                    m_Radius      = radiusHeliumMainSequenceZAMS(m_Mass);
                    m_Luminosity  = luminosityHeliumMainSequenceZAMS(m_Mass);
                }
                
                break;
                
            case FIRST_GIANT_BRANCH:
                
                
                // Star has lost its envelope, forms a HeWD if has degenerate core (M < MHeF) or zero-age naked helium star
                // see section 6.1
                
                if(m_Mass < m_massCutoffs[1]){ // M < MHeF
                    
                    // Star evolves to Helium White Dwarf - Reset parameters
                    
                    m_coreMass    = m_HeCoreMass;
                    m_Mass        = m_coreMass;
                    m_COCoreMass  = 0.0;
                    m_Age         = 0.0;
                    m_stellarType = HELIUM_WHITE_DWARF;
                    m_Radius      = radiusWhiteDwarf(m_Mass);
                }
                else{
                    
                    // Star evolves to Zero age Naked Helium Main Star and reset parameters
                    
                    m_coreMass    = m_HeCoreMass;
                    m_Mass        = m_coreMass;
                    m_Mass0       = m_Mass;
                    m_COCoreMass  = 0.0;
                    m_Age         = 0.0;
                    m_stellarType = NAKED_HELIUM_STAR_MS;
                    m_Radius      = radiusHeliumMainSequenceZAMS(m_Mass);
                }
                
                break;
                
            case CORE_HELIUM_BURNING:
                
                // if core mass > total mass then need to set stellar type to an evolved naked helium star and evolve that
                if(m_coreMass >= m_Mass){// EnvelopeLoss
                    
                    // Star becomes an evolved helium star
                    m_stellarType = NAKED_HELIUM_STAR_MS;
                    
                    // Reset toal mass to be  core mass
                    m_coreMass   = m_HeCoreMass;
                    m_Mass       = m_coreMass;
                    m_Mass0      = m_Mass;
                    m_COCoreMass = 0.0;
                    
 					// Set evolved time for naked helium star since already has some core mass.
					//Coen 10-01-2016 added prime parameters 
					double tPrime = m_Age;
					double tHeIPrime = m_timescales[2];
					double tHePrime  = m_timescales[3];
					calculateTimescales(m_timescales, m_GBParams, m_an_coefficients, m_bn_coefficients, m_massCutoffs, m_Mass0, m_coreMass, m_COCoreMass, m_Metallicity, m_logMetallicityXi, m_stellarType, m_randomSeed);
                    m_Age = ((tPrime - tHeIPrime)/tHePrime) * m_timescales[15];  // ((time - tHeI)/tHe) * tHeMS; Equation 76
                    
                    // Mass or Mass0 for GBParams?
                    calculateGiantBranchParameters(m_GBParams, m_an_coefficients, m_bn_coefficients, m_massCutoffs, m_Mass0, m_logMetallicityXi, m_stellarType, m_randomSeed);
                    
                    
                    tau = m_Age / m_timescales[15]; // tau = t/tHeMS , Just after equation 81
                    
                    m_Luminosity = luminosityHeliumMainSequence(tau, m_Mass);
                    m_Radius = radiusHeliumMainSequence(tau, m_Mass);
                    
                }
                break;
                
            case EARLY_ASYMPTOTIC_GIANT_BRANCH:
                
                // Envelope lost, form an evolved naked helium giant
                m_stellarType = NAKED_HELIUM_STAR_GIANT_BRANCH;
                
                // Set helium giant mass = He core mass
                m_Mass = m_HeCoreMass;
                m_Mass0 = m_Mass;
                // Set helium giant core mass = CO core mass
                m_coreMass = m_COCoreMass;
				
                calculateTimescales(m_timescales, m_GBParams, m_an_coefficients, m_bn_coefficients, m_massCutoffs, m_Mass0, m_coreMass, m_COCoreMass, m_Metallicity, m_logMetallicityXi, m_stellarType, m_randomSeed);
				calculateGiantBranchParameters(m_GBParams, m_an_coefficients, m_bn_coefficients, m_massCutoffs, m_Mass0, m_logMetallicityXi, m_stellarType, m_randomSeed);
                m_Age = ageEvolvedHeliumGiantBranch(m_GBParams, m_timescales, m_Mass, m_COCoreMass);
				m_Luminosity = luminosityHeliumGiantBranch(m_GBParams, m_Mass, m_COCoreMass);
				m_Radius      = radiusHeliumGiantBranch(m_Mass, m_Luminosity, m_stellarType);
                 // Given by HeGB Mc-t relation described after Equation 72
                
                break;
                
                
            case THERMALLY_PULSING_ASYMPTOTIC_GIANT_BRANCH:
                // What happens in this star type if there is a CEE? Following lines are for COCoreMass, more like a SNe process. AVG
                // Double check this
                
                McBAGB = m_GBParams[9];
                if(McBAGB < 1.6){
                    
                    // Remnant will be a Carbon-Oxygen White Dwarf (COWD)
                    
                    m_stellarType = CARBON_OXYGEN_WHITE_DWARF;
                    
                    //m_Mass     = m_Mass;
                    m_Mass     = m_coreMass;
                    m_Mass0    = m_Mass;
                    m_Age      = 0;
                    // m_envMass = 0.0; Should I add this line?
                    m_Radius        = radiusWhiteDwarf(m_Mass);
                }
                else{
                    
                    // Remnant will be an Oxygen-Neon White Dwarf (ONeWD)
                    
                    m_stellarType = OXYGEN_NEON_WHITE_DWARF;
                    
                    //m_Mass     = m_Mass;
                    m_Mass = m_coreMass;
                    m_Mass0    = m_Mass;
                    m_Age      = 0;
                    // m_envMass = 0.0; Should I add this line?
                    m_Radius        = radiusWhiteDwarf(m_Mass);
                }
                break;
                
            case NAKED_HELIUM_STAR_MS:
                // Double check this
                m_stellarType = CARBON_OXYGEN_WHITE_DWARF;
                //stripped of helium COCore remains
                m_coreMass    = m_COCoreMass;
                m_Mass        = m_coreMass;
                m_Mass0       = m_Mass;
                m_envMass     = 0.0;
                m_Age         = 0.0;
                m_Radius        = radiusWhiteDwarf(m_Mass);
                break;
                
            case NAKED_HELIUM_STAR_HETZSPRUNG_GAP:
                    m_coreMass    = m_COCoreMass;
                    m_Mass        = m_coreMass;
                    m_Mass0       = m_Mass;
                    m_envMass     = 0.0;
                    m_Age         = 0.0;                   
                
                // Remnant core becomes a white dwarf
                
                if(m_Mass < 1.6){
                    m_stellarType = CARBON_OXYGEN_WHITE_DWARF;
                    m_Radius        = radiusWhiteDwarf(m_Mass);
                }
                else{
                    m_stellarType = OXYGEN_NEON_WHITE_DWARF;
                    m_Radius        = radiusWhiteDwarf(m_Mass);
                }
                
                break;
                
            case NAKED_HELIUM_STAR_GIANT_BRANCH:
                    m_coreMass    = m_COCoreMass;
                    m_Mass        = m_coreMass;
                    m_Mass0       = m_Mass;
                    m_envMass     = 0.0;
                    m_Age         = 0.0;                   
                
                // Remnant core becomes a white dwarf
                
                if(m_Mass < 1.6){
                    m_stellarType = CARBON_OXYGEN_WHITE_DWARF;
                    m_Radius        = radiusWhiteDwarf(m_Mass);
                }
                else{
                    m_stellarType = OXYGEN_NEON_WHITE_DWARF;
                    m_Radius        = radiusWhiteDwarf(m_Mass);
                }
                
                break;
				
			// ALEJANDRO - 09/02/2017 
			// Following objects do not have an envelope to be stripped off, at least not in SSE. They keep their previous radius. 

			case HELIUM_WHITE_DWARF:
				break;
                 
			case CARBON_OXYGEN_WHITE_DWARF:
				break;
           
			case OXYGEN_NEON_WHITE_DWARF:
                break;
				
			case NEUTRON_STAR:
				break;
				
			case BLACK_HOLE:
				break;
				
			case MASSLESS_REMNANT:

                if(debugging)
                    std::cout << "Star is a massless remnant" << std::endl;

                break;
                
            default:
                
                    std::cerr << m_randomSeed << "\tError in modifyStarAfterLosingEnvelope function. Unknown stellar type : " << m_stellarType << std::endl;
                
                if(m_error == false){
                    m_error = true;
                }
                
                break;
            
        }
    }
    else{
        m_COCoreMass  = 0.0;
        m_HeCoreMass  = 0.0;
        m_coreMass    = 0.0;
        m_Mass0       = 0.0;
        m_envMass     = 0.0;
        m_coreRadius  = 0.0;
        m_Luminosity  = 0.0;
        m_Radius      = 0.0;
        m_Temperature = 0.0;
        m_stellarType = MASSLESS_REMNANT;
        
        if(debugging){
            std::cout << "Mass = " << m_Mass << std::endl;
            std::cout << "Star is a massless remnant" << std::endl;
        }
    }
 
    m_envMass = 0.0;
}


double Star::radiusRemnantStarAfterLosingEnvelope(){
// ALEJANDRO - 14/02/2017 - Function to determine the "core radius" or the core after the star losses its whole envelope, as explained after eq. 105 Hurley+2000. 	
    double  tau = 0.0;
    double  McBAGB = 0.0; // It is initialized like this in evolveOneTimestep function. Look for this in detail. AVG.
    bool    debugging = false;
    
    if(m_Mass > 0.0){
        switch (m_stellarType) {
            case MS_LESS_THAN_07:
				// No clear core-envelope structure, therefore not modify radius.
                break;
                
            case MS_MORE_THAN_07:
                // No clear core-envelope structure, therefore not modify radius.
                break;
                
            case HERTZSPRUNG_GAP:
                
                // Star has lost its envelope, forms a HeWD if has degenerate core (M < MHeF) or zero-age naked helium star
                // Described just above Equation 76
                
                if(m_Mass < m_massCutoffs[1]){                    
                    // Reset mass/age parameters
                    m_Radius        = radiusWhiteDwarf(m_Mass);
                    
                }
                else{
                    // Zero age Helium MS star // see section 6.1                    
                    m_Radius      = radiusHeliumMainSequenceZAMS(m_Mass);
                }
                
                break;
                
            case FIRST_GIANT_BRANCH:
                
                // Star has lost its envelope, forms a HeWD if has degenerate core (M < MHeF) or zero-age naked helium star
                // see section 6.1
                if(m_Mass < m_massCutoffs[1]){ // M < MHeF
                    
                    // Star evolves to Helium White Dwarf - Reset parameters
                    m_Radius      = radiusWhiteDwarf(m_Mass);
                }
                else{
                    
                    // Star evolves to Zero age Naked Helium Main Star and reset parameters
                    m_Radius      = radiusHeliumMainSequenceZAMS(m_Mass);
                }
                
                break;
                
            case CORE_HELIUM_BURNING:
                
                // if core mass > total mass then need to set stellar type to an evolved naked helium star and evolve that
                if(m_coreMass >= m_Mass){// EnvelopeLoss
                    
                    // Star becomes an evolved helium star
                    m_stellarType = NAKED_HELIUM_STAR_MS;
                    
                    // Reset core mass to be total mass
                    m_coreMass   = m_Mass;
                    m_Mass0      = m_Mass;
                    m_HeCoreMass = m_coreMass;
                    m_COCoreMass = 0.0;
                    
 					// Set evolved time for naked helium star since already has some core mass.
					//Coen 10-01-2016 added prime parameters 
					double tPrime = m_Age;
					double tHeIPrime = m_timescales[2];
					double tHePrime  = m_timescales[3];
					calculateTimescales(m_timescales, m_GBParams, m_an_coefficients, m_bn_coefficients, m_massCutoffs, m_Mass0, m_coreMass, m_COCoreMass, m_Metallicity, m_logMetallicityXi, m_stellarType, m_randomSeed);
                    m_Age = ((tPrime - tHeIPrime)/tHePrime) * m_timescales[15];  // ((time - tHeI)/tHe) * tHeMS; Equation 76
                    
                    // Mass or Mass0 for GBParams?
                    calculateGiantBranchParameters(m_GBParams, m_an_coefficients, m_bn_coefficients, m_massCutoffs, m_Mass0, m_logMetallicityXi, m_stellarType, m_randomSeed);
                    
                    
                    tau = m_Age / m_timescales[15]; // tau = t/tHeMS , Just after equation 81
                    
                    m_Luminosity = luminosityHeliumMainSequence(tau, m_Mass);
                    m_Radius = radiusHeliumMainSequence(tau, m_Mass);
                    
                }
                break;
                
            case EARLY_ASYMPTOTIC_GIANT_BRANCH:
                
                // Envelope lost, form an evolved naked helium giant
                m_stellarType = NAKED_HELIUM_STAR_GIANT_BRANCH;
                
                // Set helium giant mass = He core mass
                m_Mass = m_HeCoreMass;
                m_Mass0 = m_Mass;
                // Set helium giant core mass = CO core mass
                m_coreMass = m_COCoreMass;
				m_Age = ageEvolvedHeliumGiantBranch(m_GBParams, m_timescales, m_Mass, m_COCoreMass);
				calculateGiantBranchParameters(m_GBParams, m_an_coefficients, m_bn_coefficients, m_massCutoffs, m_Mass0, m_logMetallicityXi, m_stellarType, m_randomSeed);
				m_Luminosity = luminosityHeliumGiantBranch(m_GBParams, m_Mass, m_COCoreMass);
				m_Radius      = radiusHeliumGiantBranch(m_Mass, m_Luminosity, m_stellarType);
                 // Given by HeGB Mc-t relation described after Equation 72
                
                break;
                
                
            case THERMALLY_PULSING_ASYMPTOTIC_GIANT_BRANCH:
                // What happens in this star type if there is a CEE? Following lines are for COCoreMass, more like a SNe process. AVG
                // Double check this
                if(McBAGB < 1.6){
                    
                    // Remnant will be a Carbon-Oxygen White Dwarf (COWD)
                    m_Radius        = radiusWhiteDwarf(m_Mass);
                }
                else{
                    
                    // Remnant will be an Oxygen-Neon White Dwarf (ONeWD)
                    m_Radius        = radiusWhiteDwarf(m_Mass);
                }
                break;
                
            case NAKED_HELIUM_STAR_MS:
                m_Radius        = radiusWhiteDwarf(m_Mass);
                break;
                
            case NAKED_HELIUM_STAR_HETZSPRUNG_GAP:              
                
                // Remnant core becomes a white dwarf
                
                if(m_Mass < 1.6){
                    
                    m_Radius        = radiusWhiteDwarf(m_Mass);
                }
                else{
                    m_Radius        = radiusWhiteDwarf(m_Mass);
                }
                
                break;
                
            case NAKED_HELIUM_STAR_GIANT_BRANCH:
                
                if(m_Mass < 1.6){
                    m_Radius        = radiusWhiteDwarf(m_Mass);
                }
                else{
                    m_Radius        = radiusWhiteDwarf(m_Mass);
                }
                break;
				
			// ALEJANDRO - 09/02/2017 
			// Following objects do not have an envelope to be stripped off, at least not in SSE. They keep their previous radius. 

			case HELIUM_WHITE_DWARF:
				break;
                 
			case CARBON_OXYGEN_WHITE_DWARF:
				break;
           
			case OXYGEN_NEON_WHITE_DWARF:
                break;
				
			case NEUTRON_STAR:
				break;
				
			case BLACK_HOLE:
				break;
				
			case MASSLESS_REMNANT:

                if(debugging)
                    std::cout << "Star is a massless remnant" << std::endl;

                break;
                
            default:
                
                    std::cerr << m_randomSeed << "\tError in modifyStarAfterLosingEnvelope function. Unknown stellar type : " << m_stellarType << std::endl;
                
                if(m_error == false){
                    m_error = true;
                }
                
                break;
            
        }
    }
    else{
        m_Radius      = 0.0;
        
        if(debugging){
            std::cout << "Mass = " << m_Mass << std::endl;
            std::cout << "Star is a massless remnant" << std::endl;
        }
    }
 
    return	m_Radius;
}

double calculateRadialExtentConvectiveEnvelope(Star & star){
	// Function to calculate the radial extent of the convective envelope of stars with type 0,1,3,5,6,8,9 
	// For details see sec. 2.3, particularly eqs. 36-40 of subsec. 2.3.1 of Hurley+2002
	// [radialExtentConvectiveEnvelope] == Rsol
	
	bool	debugging = false;
	double	Renv = NEVER_SET;
	
	Star	starCopy= star;
	double	Rc	=	starCopy.radiusRemnantStarAfterLosingEnvelope();	// Radius of remnant of primary after stripping its envelope, in Rsol

	if(debugging){std::cout<< star.m_randomSeed << "\tInside radialExtentConvectiveEnvelope() function." << std::endl;}
	if(debugging){std::cout<< star.m_randomSeed << "\tStar with convective envelope" << std::endl;}
	
	if(star.m_stellarType == MS_LESS_THAN_07 || star.m_stellarType == MS_MORE_THAN_07){
		double Renv0	= NEVER_SET;
		
		if (star.m_Mass <= 0.35){
			Renv0 = star.m_Radius;
			if(debugging){std::cout<< star.m_randomSeed << "\tm_Mass, Renv_0:\t" << star.m_Mass << "\t" << Renv0 << std::endl;}
		}
		else if (star.m_Mass < 1.25){
			if(debugging){std::cout<< star.m_randomSeed << "\tCurrenly R' (eq. 36) not implemented. Solve and delete this comment." << std::endl;}
            // m_RZAMS0 = RadiusZAMS(m_Mass0, radius_coefficient_theta, radius_coefficient_iota, radius_coefficient_kappa, radius_coefficient_lamda, radius_coefficient_mu, radius_coefficient_nu, radius_coefficient_xi, radius_coefficient_omicron, radius_coefficient_Pi);
            // m_Radius        = radiusMainSequence(m_an_coefficients, m_massCutoffs, m_timescales, m_RZAMS0, m_Mass, m_Age, m_error, m_randomSeed);
			// Probably best way to calculate Rprime is  to make a table of M=0.35 Msol at different metallicities during different metallicities and fractional ages.
			double	Rprime = NEVER_SET;
			//Renv0 = Rprime*pow((1.25-m_Mass)/0.9,1.0/2.0);
			Renv0 = NEVER_SET;
			if(debugging){std::cout<< star.m_randomSeed << "\tm_Mass, R' and Renv_0:\t" << star.m_Mass << "\t" << Rprime << "\t" << Renv0 << std::endl;}
		}
		else{
			Renv0=0.0;
			if(debugging){std::cout<< star.m_randomSeed << "\tm_Mass, Renv_0:\t" << star.m_Mass << "\t" << Renv0 << std::endl;}
		}
		
		Renv = Renv0*pow(1-star.m_tau,1.0/4.0);
	}
	else if(star.m_stellarType == HERTZSPRUNG_GAP){
		if(debugging){std::cout<< star.m_randomSeed << "\ttau, Radius, Rc:\t" << star.m_tau << "\t" << star.m_Radius << "\t" << Rc << std::endl;}
		Renv=pow(star.m_tau,1.0/2.0)*(star.m_Radius-Rc);
	}
	else if(star.m_stellarType == FIRST_GIANT_BRANCH || star.m_stellarType == EARLY_ASYMPTOTIC_GIANT_BRANCH ||
			 star.m_stellarType == THERMALLY_PULSING_ASYMPTOTIC_GIANT_BRANCH || star.m_stellarType == NAKED_HELIUM_STAR_HETZSPRUNG_GAP ||
			 star.m_stellarType == NAKED_HELIUM_STAR_GIANT_BRANCH){
		Renv=star.m_Radius-Rc; 
			 }
	else{
		std::cerr<< star.m_randomSeed << "Shouldn't get here, stellar type doesn't have a convective envelope." << std::endl;
		Renv=NEVER_SET; 
	}

	if(debugging){std::cout<< star.m_randomSeed << "\tStellarType, Renv:\t" << star.m_stellarType << "\t" << Renv << std::endl;}		
			
	return	Renv;
}

double calculateEddyTurnoverTimescale(Star & star){
	// Function to calculate the eddy turnover timescale 
	// For details see sec. 2.3, particularly eq. 31 of subsec. 2.3.1 of Hurley+2002
	// [eddyTurnoverTimescale] = yr

	bool	debugging = false;

	if(debugging){std::cout<< star.m_randomSeed << "\tInside eddyTurnoverTimescale function" << std::endl;}
	
	double	Renv 		= calculateRadialExtentConvectiveEnvelope(star);
	double	numerator 	= star.m_Mass*Renv*(star.m_Radius-0.5*Renv);
	double	denominator = 3*star.m_Luminosity;
	double	tau_conv	=0.4311*pow(numerator/denominator,1.0/3.0);
	
	if(debugging){std::cout<< star.m_randomSeed << "\tRenv,tau_conv:\t" << Renv << "\t" << tau_conv << std::endl;}

	return	tau_conv;
}

double calculatekOverTc(Star & star){
    // Function to calculate factor (f/T)_c for stars with convective envelope
    // For details see section 2.3.1, particularly eq. 30 of Hurley+2002
    // [(k/T)_c] = yr
    bool    debugging = false;
    double  f_conv = 1.0;
    double  tau_conv    = calculateEddyTurnoverTimescale(star);
    double  kOverT_c    = (2.0/21.0)*(f_conv/tau_conv)*((star.m_Mass-star.m_coreMass)/star.m_Mass);

    if(debugging){
        std::cout << star.m_randomSeed << "\ttau_conv:\t" << tau_conv << std::endl;
        std::cout << star.m_randomSeed << "\tf_conv:\t" << f_conv << std::endl;
        std::cout << "kOverT_c:\t" << kOverT_c << "yr"<< std::endl;
    }

    return  kOverT_c;
}

double calculateCircularisationTimescale(Star & star, double Porb, double a, double mCompanion){
	// Function to calculate the circularization timescale
	// For details see sec. 2.3 of Hurley+2002
	// [circularisationTimescale] = yr^{-1}

	bool	debugging = false;
	bool	mainSequenceConvectiveFlag = false;
	double	circularisationTimescale = NEVER_SET;
	double	denominator = NEVER_SET;
	double	q2			= mCompanion/star.m_Mass;
    double  kOverT_c    = NEVER_SET;
	
	if(debugging){std::cout << star.m_randomSeed << "\tInside calculateCircularisationTimescale." << std::endl;}
	
	if(debugging){
			std::cout << star.m_randomSeed << "\tMass, envMass, coreMass:\t" << star.m_Mass << "\t" << star.m_envMass << "\t" << star.m_coreMass << std::endl;	
			std::cout << star.m_randomSeed << "\tRadius, separation:\t" << star.m_Radius << "\t" << a << std::endl;
			std::cout << star.m_randomSeed << "\tStellartype:\t" << star.m_stellarType << std::endl;
			std::cout << star.m_randomSeed << "\tq2:\t" << q2 << std::endl;
	}

	
	// Determine if MS star should be treated as with convective or radiative envelope; see subsection 2.3.1
	if(star.m_Mass < 1.25 and star.m_stellarType < HERTZSPRUNG_GAP)
		mainSequenceConvectiveFlag = true;
	
	if(	star.m_stellarType == HERTZSPRUNG_GAP ||  star.m_stellarType == FIRST_GIANT_BRANCH ||  
		star.m_stellarType == EARLY_ASYMPTOTIC_GIANT_BRANCH || star.m_stellarType == THERMALLY_PULSING_ASYMPTOTIC_GIANT_BRANCH || 
		star.m_stellarType == NAKED_HELIUM_STAR_HETZSPRUNG_GAP || star.m_stellarType == NAKED_HELIUM_STAR_GIANT_BRANCH ||
		mainSequenceConvectiveFlag == true){
		// Solve for stars with convective envelope, according to tides section
		// For more details, see sec. 2.3.1 of Hurley+2002
		if(debugging){std::cout << star.m_randomSeed << "\tThe star has a convective envelope according to Hurley tides section. This include all MS stars." << std::endl;}	
	
		double	tau_conv	= calculateEddyTurnoverTimescale(star);
		double	f_conv = 1.0;	// Currently, as COMPAS doesn't has rotating stars tested, we set f_conv = 1 always. This should be solved, as shown commented below, once rotation is included.
		//	double	P_tid		= NEVER_SET;
		//	double	Pspin = 2*pi/star.m_omega;
		//	if(Pspin <= 0.0){
		//		f_conv = 1.0;
		//	}
		//	else{
		//
		//		P_tid 	= 1.0/std::abs((1.0/Porb)-(1.0/Pspin));
		//		f_conv	= std::min(1.0,pow(P_tid/(2.0*tau_conv),2));
		//	}

        kOverT_c = calculatekOverTc(star);      
		denominator = (f_conv/tau_conv)*((star.m_Mass-star.m_coreMass)/star.m_Mass)*q2*(1.0+q2)*pow(star.m_Radius/a,8.0);
		circularisationTimescale = 1.0/denominator;

		if(debugging)
			std::cout << star.m_randomSeed << "\tq2:\t" << q2 << std::endl;
		

	}
	else if(star.m_stellarType == MS_LESS_THAN_07 ||  star.m_stellarType == MS_MORE_THAN_07 ||
			star.m_stellarType == CORE_HELIUM_BURNING || star.m_stellarType == NAKED_HELIUM_STAR_MS){
		// Solve for stars with radiative envelope
		// For more details, see sec. 2.3.2 of Hurley+2002
		if(debugging){std::cout << star.m_randomSeed << "\tThe star has a radiative envelope according to Hurley tides section." << std::endl;}	
	
		double	secondOrderTidalCoefficient = 1.592*pow(10.0,-9.0)*pow(star.m_Mass,2.84);	// Also known as E_2.
		double	freeFallFactor = sqrt(G1*star.m_Mass/pow(star.m_Radius*RsolToAU,3.0));
		
		denominator = (21.0/2.0)*freeFallFactor*q2*pow(1.0+q2,11.0/6.0)*secondOrderTidalCoefficient*pow(star.m_Radius/a,21.0/2.0);
		circularisationTimescale = 1.0/denominator;

		if(debugging){
			std::cout << star.m_randomSeed << "\tsecondOrderTidalCoefficient:\t" << secondOrderTidalCoefficient << std::endl;
			std::cout << star.m_randomSeed << "\tfreeFallFactor:\t" << freeFallFactor << std::endl;
		}
	}
	else{
		if(debugging){std::cout << star.m_randomSeed << "\tCircularization timescale being calculate for a compact object." << std::endl;}
		// Circularization timescale for a compact object currently set to 0. 
		// Can be calculate, e.g. using the Peters formula for circularisation gravitational radiation from 2 point particles.
		circularisationTimescale = 0.0;	
	}
	
	if(debugging){std::cout << star.m_randomSeed << "\tcircularisationTimescale:\t" << circularisationTimescale << std::endl << std::endl;}

	return	circularisationTimescale;
}
	
double calculateSyncrhonizationTimescale(Star & star, double Porb, double a, double mCompanion){
	// Function to calculate the synchronization timescale
	// For details see sec. 2.3 of Hurley+2002
	// [syncrhonizationTimescale] = yr

	bool	debugging = false;
	bool	mainSequenceConvectiveFlag = false;
	double	syncrhonizationTimescale = NEVER_SET;
	double	denominator = NEVER_SET;
	double	q2			= mCompanion/star.m_Mass;
	double	gyrationRadiusSquared = k_definition(star.m_RZAMS, star.m_Radius, star.m_Mass, star.m_stellarType);
    double  kOverT_c    = NEVER_SET;
	
	if(debugging){std::cout << star.m_randomSeed << "\tInside calculateSyncrhonizationTimescale." << std::endl;}
	
	if(debugging){
			std::cout << star.m_randomSeed << "\tMass, envMass, coreMass:\t" << star.m_Mass << "\t" << star.m_envMass << "\t" << star.m_coreMass << std::endl;	
			std::cout << star.m_randomSeed << "\tRadius, separation:\t" << star.m_Radius << "\t" << a << std::endl;
			std::cout << star.m_randomSeed << "\tStellartype:\t" << star.m_stellarType << std::endl;
			std::cout << star.m_randomSeed << "\tgyrationRadiusSquared:\t" << gyrationRadiusSquared << std::endl;
	}

	// Determine if MS star should be treated as with convective or radiative envelope; see subsection 2.3.1
	if(star.m_Mass < 1.25 and star.m_stellarType < HERTZSPRUNG_GAP)
		mainSequenceConvectiveFlag = true;
	
	if(	star.m_stellarType == HERTZSPRUNG_GAP ||  star.m_stellarType == FIRST_GIANT_BRANCH ||  
		star.m_stellarType == EARLY_ASYMPTOTIC_GIANT_BRANCH || star.m_stellarType == THERMALLY_PULSING_ASYMPTOTIC_GIANT_BRANCH || 
		star.m_stellarType == NAKED_HELIUM_STAR_HETZSPRUNG_GAP || star.m_stellarType == NAKED_HELIUM_STAR_GIANT_BRANCH ||
		mainSequenceConvectiveFlag == true){
		// Solve for stars with convective envelope, according to tides section
		// For more details, see sec. 2.3.1 of Hurley+2002
		if(debugging){std::cout << star.m_randomSeed << "\tThe star has a convective envelope according to Hurley tides section. This include all MS stars." << std::endl;}	

		double	tau_conv	= calculateEddyTurnoverTimescale(star);
		double	f_conv = 1.0;	// Currently, as COMPAS doesn't has rotating stars tested, we set f_conv = 1 always. This should be solved, as shown commented below, once rotation is included.
		//	double	P_tid		= NEVER_SET;
		//	double	Pspin = 2*pi/star.m_omega;
		//	if(Pspin <= 0.0){
		//		f_conv = 1.0;
		//	}
		//	else{
		//
		//		P_tid 	= 1.0/std::abs((1.0/Porb)-(1.0/Pspin));
		//		f_conv	= std::min(1.0,pow(P_tid/(2.0*tau_conv),2));
		//	}

        kOverT_c = calculatekOverTc(star);      		
		denominator = 3*kOverT_c*q2*pow(gyrationRadiusSquared,-1.0)*pow(star.m_Radius/a,6.0);

		if(debugging)
			std::cout << star.m_randomSeed << "\tq2:\t" << q2 << std::endl;

	}
	else if(star.m_stellarType == MS_LESS_THAN_07 ||  star.m_stellarType == MS_MORE_THAN_07 ||
			star.m_stellarType == CORE_HELIUM_BURNING || star.m_stellarType == NAKED_HELIUM_STAR_MS){
		// Solve for stars with radiative envelope
		// For more details, see sec. 2.3.2 of Hurley+2002
		if(debugging){std::cout << star.m_randomSeed << "\tThe star has a radiative envelope according to Hurley tides section." << std::endl;}	

		double	secondOrderTidalCoefficient = 1.592*pow(10.0,-9.0)*pow(star.m_Mass,2.84);	// Also known as E_2.
		double	freeFallFactor = sqrt(G1*star.m_Mass/pow(star.m_Radius*RsolToAU,3.0));
		
		denominator = pow(52.0,5.0/3.0)*freeFallFactor*pow(gyrationRadiusSquared,-1.0)*q2*q2*pow(1.0+q2,5.0/6.0)*secondOrderTidalCoefficient*pow(star.m_Radius/a,17.0/2.0);

		if(debugging){
			std::cout << star.m_randomSeed << "\tsecondOrderTidalCoefficient:\t" << secondOrderTidalCoefficient << std::endl;
			std::cout << star.m_randomSeed << "\tfreeFallFactor:\t" << freeFallFactor << std::endl;
		}
	}
	else{
		if(debugging){std::cout << star.m_randomSeed << "\tSynchronization timescale being calculate for a compact object." << std::endl;}
		
		denominator = (1.0/(1.3*pow(10.0,7.0)))*pow(star.m_Luminosity/star.m_Mass,5.0/7.0)*pow(star.m_Radius/a,6.0);
	}
	
	syncrhonizationTimescale = 1.0/denominator;
	if(debugging){std::cout << star.m_randomSeed << "\tsyncrhonizationTimescale:\t" << syncrhonizationTimescale << std::endl << std::endl;}
		
	return	syncrhonizationTimescale;
}

double k_definition(double RZAMS, double Radius, double Mass, int stellarType){
    // Define gyration radius 'k=r_g^2' using fit from de Mink et al. 2013, calling k_definition function
    // Created by Alejandro Vigna-Gomez on 11/2015.
	
	// ALEJANDRO - 14/02/2017 - This should a function of the star class, not the binary class. Ilya changed that in his own version of COMPAS.
	// The original fits from de Mink+2013 where made for MS stars a Z=0.02.
    double 	k0,c,C;
	double	gyrationRadius = 1.0;
    double 	mass=log10(Mass);
    double 	ratio=Radius/RZAMS;
	bool	debugging = false; //true;

	if((stellarType==MS_LESS_THAN_07)or(stellarType==MS_MORE_THAN_07)){

		if(mass<=1.3){
			c=0;
		}
		else{
			c=-0.055*pow(mass-1.3,2.0);
		}

		if(mass<=0){
			C=-2.5;
		}
		else if((mass>0)&&(mass<=0.2)){
			C=-2.5+5*mass;
		}
		else{
			C=-1.5;
		}

		k0=c+std::min(0.21,std::max(0.09-0.27*mass,0.037+0.033*mass));
		gyrationRadius = ((k0-0.025)*pow(ratio,C))+(0.025*pow(ratio,-0.1));
		
		if(debugging){
			std::cout << "Inside k_definition() function" << std::endl;
			std::cout << "log10(mass):\t" << mass << std::endl;
			std::cout << "ratio:\t" << ratio << std::endl;
			std::cout << "c:\t" << c << std::endl;
			std::cout << "c:\t" << c << std::endl;
			std::cout << "C:\t" << C << std::endl;
			std::cout << "k0:\t" << k0 << std::endl;
			std::cout << "gyrationRadius:\t" << gyrationRadius << std::endl;
 		}
	}
	else if(stellarType == NAKED_HELIUM_STAR_MS){
		// Test for Naked Helium Star MS. Not sure about the moment. 
		gyrationRadius = 0.1;		
	}
	else if(stellarType == FIRST_GIANT_BRANCH or stellarType == EARLY_ASYMPTOTIC_GIANT_BRANCH or stellarType == THERMALLY_PULSING_ASYMPTOTIC_GIANT_BRANCH or stellarType==NAKED_HELIUM_STAR_GIANT_BRANCH){
		// As defined after eq. 109 of Hurley+2000 for giants. Single number approximation.
		gyrationRadius = 0.1;
	}
	else if(stellarType == HERTZSPRUNG_GAP or stellarType == CORE_HELIUM_BURNING or stellarType == NAKED_HELIUM_STAR_HETZSPRUNG_GAP or stellarType == HELIUM_WHITE_DWARF or stellarType == CARBON_OXYGEN_WHITE_DWARF or stellarType == OXYGEN_NEON_WHITE_DWARF or stellarType == NEUTRON_STAR){
		// As defined after eq. 109 of Hurley+2000 for n=3/2 polytrope or dense convective core. Single number approximation.
		gyrationRadius = 0.21;
	}
	else if(stellarType == BLACK_HOLE){
		if(debugging){
			std::cout << "Not sure about what to do with the gyration radius of a BH. It is point particle but you can define it as a solid sphere (k=2/5) with a R=R_{Schwarzschild}." << std::endl;
		} 
		gyrationRadius = 2.0/5.0;
	}
	else{
		// For massless remnant
		gyrationRadius = 0.0;
	}

	return gyrationRadius;	
}


double momentOfInertia(Star starCopy, double Mass, double coreMass, double envMass, double Radius, double Rc, double RZAMS, int stellarType){
	// ALEJANDRO - 15/02/2017 - Function to calculate the moment of inertia of the star
	// I = k*M*R^2
	// [I] = Msol * AU^2
	
	double  I = 0.0;													// Moment of inertia. Value to return.
	bool	debugging = false;
	
	if(debugging){
		std::cout << "Inside momentOfInertia() function" << std::endl;
		std::cout << "Stellartype:\t" << stellarType << std::endl;
		std::cout << "Mass [Msol]:\t" << Mass << std::endl;
		std::cout << "coreMass [Msol]:\t" << coreMass << std::endl;
		std::cout << "envMass [Msol]:\t" << envMass << std::endl;
		std::cout << "Radius [AU]:\t" << Radius << std::endl;
		std::cout << "Rc [AU]:\t" << Rc << std::endl;
	}
	
	if((stellarType==MS_LESS_THAN_07)or(stellarType==MS_MORE_THAN_07)or(stellarType==NAKED_HELIUM_STAR_MS)){
		double	gyrationRadius = k_definition(RZAMS, Radius, Mass, stellarType);
		I = gyrationRadius*Mass*Radius*Radius;
		if(debugging){
			std::cout << "Moment of inertia of a MS star, I:\t" << I << std::endl;
			std::cout << "gyration radius:\t" << gyrationRadius << std::endl;
		}
	}
	else if((stellarType==FIRST_GIANT_BRANCH)or(stellarType==EARLY_ASYMPTOTIC_GIANT_BRANCH)or(stellarType==THERMALLY_PULSING_ASYMPTOTIC_GIANT_BRANCH)or(stellarType==NAKED_HELIUM_STAR_GIANT_BRANCH)or(stellarType==HERTZSPRUNG_GAP)or(stellarType==CORE_HELIUM_BURNING)or(stellarType==NAKED_HELIUM_STAR_HETZSPRUNG_GAP)or(stellarType==HELIUM_WHITE_DWARF)or(stellarType==CARBON_OXYGEN_WHITE_DWARF)or(stellarType==OXYGEN_NEON_WHITE_DWARF)or(stellarType==NEUTRON_STAR)){
		// k2=0.1 and k3=0.21 as defined after eq. 109 Hurley+2000
		double	momentOfInertiaCore = 0.1*envMass*Radius*Radius;
		double	momentOfInertiaEnvelope = 0.21*coreMass*Rc*Rc;
		I = momentOfInertiaCore + momentOfInertiaEnvelope;
		
		if(debugging){
			std::cout << "Moment of inertia of an evolved star, Icore, Ienv, Itotal:\t" << momentOfInertiaCore << "\t" << momentOfInertiaEnvelope << "\t" << I << std::endl;
		}
	}
	else if(stellarType==BLACK_HOLE){
		I = (2.0/5.0)*Mass*Radius*Radius;
		if(debugging){
			std::cout << "Not sure about what to do with the moment of inertia of a BH. It is point particle but you can define it as a solid sphere (k=2/5) with a R=R_{Schwarzschild}." << std::endl;
		} 
	}
	else{
		if(debugging){
			std::cout << "Massless remnants don't have a momment of inertia. Setting to I_{MR} = 0.0" << std::endl;
		} 
	}
	
	return I;
}

double	calculateRadialExpansionTimescale(double Radius, double RadiusPrev, double stellarTypePrev, double stellarType, double dt, double time, double tauDynamical){
	// Function to calculate the radial expansion timescale
	// [radialExpansionTimescale] = yr
	
	bool	debugging = false;
	double	tauRadialExpansion = NEVER_SET;
    double  timePrev = time-dt;
    double  tauDynamicalYears = tauDynamical*MyearToyear;

	
	if (stellarTypePrev==stellarType and RadiusPrev!=Radius){
		tauRadialExpansion = dt*MyearToyear*RadiusPrev/(Radius-RadiusPrev);
	}
	else{
        tauRadialExpansion = tauDynamicalYears;
	}

    if (tauRadialExpansion < 0.0 and debugging){
        std::cout << "Radial expansion timescale is negative. Star is contracting." << std::endl;
    }

    if(debugging){
        std::cout << "StellartypePrev, Stellartype:\t" << stellarTypePrev << " " << stellarType << std::endl;
        std::cout << "RadiusPrev, Radius [Rsol]:\t" << RadiusPrev << " " << Radius << std::endl;
        std::cout << "dtPrev [yr]:\t" << dt*MyearToyear << std::endl;
        std::cout << "timePrev, time [Myr]\t" << timePrev << " ," << time << std::endl;
        std::cout << "tauRadialExpansion [yr]:\t" << tauRadialExpansion << std::endl;
        std::cout << "tauDynamical [yr]:\t" << tauDynamicalYears << std::endl << std::endl;        
    }

	return tauRadialExpansion;
}

double bindingEnergy(double coreMass, double envMass, double Radius, double Lambda, unsigned long randomSeed){
	//	ALEJANDRO - 08/03/2017 - Function for calculating the absolute value of the binding energy core to the envelope of the star
	//	|E_{bind}| = G*M_{donor}*M_{env}/(\lambda*R_{donor})
	//	[E_{bind}] = ergs
	
	double	coreMassCGS		= coreMass*MsolTog;
	double	envMassCGS		= envMass*MsolTog;
	double	totalMassCGS	= coreMassCGS+envMassCGS;
	double	radiusCGS		= Radius*RsolTocm;
	double	bindingEnergyCGS= 0.0;
	bool	debugging		= false;
	
	if(radiusCGS <= 0.0){
			std::cerr << randomSeed << "\tRadius of star <= 0. Binding energy not set. Shouldn't get here.";
	}
	else if(Lambda <= 0.0){
		// Not necesarily zero as sometimes lambda is made 0, or maybe weird values for certain parameters of the fit. Not sure about the latter.
		if(debugging){
			std::cout << "Lambda of star <= 0. Binding energy not set.";
		}		
	}
	else{
		bindingEnergyCGS	= Gcgs*totalMassCGS*envMassCGS/(Lambda*radiusCGS);
	}
	
	if(debugging){
		std::cout << "Inside binding energy function" << std::endl;
		std::cout << "coreMass [g]:\t" << coreMassCGS << std::endl;
		std::cout << "envMass [g]:\t" << envMassCGS << std::endl;
		std::cout << "totalMass [g]:\t" << totalMassCGS << std::endl;
		std::cout << "Radius [m]:\t" << radiusCGS << std::endl;
		std::cout << "lambda:\t" << Lambda << std::endl;
		std::cout << "Binding energy [ergs]:\t" << bindingEnergyCGS << std::endl;
	}
	
	return	bindingEnergyCGS;
}

double ZadiabaticSPH(double md, double mdHEcore){
    // Calculation of Zadiabatic according to the fit from Soberman, Phinney, vdHeuvel (1997)
    double  m   = mdHEcore/md;          // Eq (57) Soberman, Phinney, vdHeuvel (1997)
    return      ((2.0/3.0)*m/(1.0-m))-((1.0/3.0)*((1.0-m)/(1.0+2*m)))-(0.03*m)+(0.2*m/(1.0+pow(1.0-m,-6.0)));                        // Eq (61) Soberman, Phinney, vdHeuvel (1997)
}

double ZadiabaticHurley2002(double md, double mdHEcore){
    // From Hurley et al 2002
    double m = mdHEcore/md;
    double x = -0.3; // Depends on composition, should use x from Hurley et al 2000
    return -x + (2.0*m*m*m*m*m);
}

double dynamicalTimescale(double Mass, double Radius){
    // Calculation of dynamical timescale, eq. (1), Kalogera & Webbink 1996
    // [Mass] = Msol
    // [Radius] = Rsol
    // [timescale] = Myr

    bool    debugging = false;
    double  dynamicalTimescaleYears= 5.0*pow(10.0,-5.0)*pow(Mass,-0.5)*pow(Radius,1.5);
    double  dynamicalTimescaleMYears = dynamicalTimescaleYears*yearToMyear;

    if(debugging){
        std::cout << "Mass [Msol]:\t" << Mass << std::endl;
        std::cout << "Radius [Rsol]:\t" << Radius << std::endl;
        std::cout << "Tdyn [yrs]: " << dynamicalTimescaleYears << std::endl;
        std::cout << "Tdyn [Myrs]: " << dynamicalTimescaleMYears << std::endl;
    }
    return dynamicalTimescaleMYears;
}

double thermalTimescale(double Mass, double Menv, double Radius, double Luminosity, int stellarType){
    // Calculation of mass transfer thermal timescale, eq. (2), Kalogera & Webbink 1996
    // Check if the star has an envelope
    // G*Msol^2/(Lsol*Rsol) ~ 30
    // [Mass] = Msol
    // [Radois] = Rsol
    // [Luminosity] = Lsol
    // [timescale] = Myr

    bool    debugging = false;

    if(stellarType == MS_LESS_THAN_07 or stellarType == MS_MORE_THAN_07 or stellarType == NAKED_HELIUM_STAR_MS){
        return 30*(pow(Mass,2.0))/(Radius*Luminosity);
    }
    else if(stellarType == HERTZSPRUNG_GAP or stellarType == FIRST_GIANT_BRANCH or stellarType == CORE_HELIUM_BURNING or stellarType ==EARLY_ASYMPTOTIC_GIANT_BRANCH or stellarType == THERMALLY_PULSING_ASYMPTOTIC_GIANT_BRANCH or stellarType ==NAKED_HELIUM_STAR_HETZSPRUNG_GAP or stellarType == NAKED_HELIUM_STAR_GIANT_BRANCH){
        return 30*Mass*Menv/(Radius*Luminosity);
    }
    else{
        if(debugging){std::cout << "MassTransferTimescaleThermal. Non-valid star type. Return 0." << std::endl;}
        return 0;
    }

}

double nuclearTimescale(double Mass, double Luminosity){
    // Calculation of nuclear timescale, eq. (3), Kalogera & Webbink 1996
    // [Mass] = Msol
    // [Luminosity] = Lsol
    // [timescale] = Myr
    return pow(10.0,10.0)*Mass/Luminosity*yearToMyear;
}

int envelopeType(int stellarType, const programOptions & options, unsigned long m_randomSeed){
    /* Function to determine if the star has a radiative or convective envelope. Some calculations on this can be found in sec. 2.3.4 of Belczynski et al. 2008. For now, we will only do the calculation using stellarType.
        EnvelopeType = 0 for a Radiative Envelope.
        EnvelopeType = 1 for a Convective Envelope.
        EnvelopeType = 2 for a Remmant.
     */

    if(options.commonEnvelopeHertzsprungGapDonor == STABLE_HG_CE){

        if(stellarType == MS_MORE_THAN_07 or stellarType == NAKED_HELIUM_STAR_MS or stellarType == HERTZSPRUNG_GAP){
            return RADIATIVE_ENVELOPE;
        }
		// ALEJANDRO - 12/09/2018 - Rechecking this function, don't thins CHeB type should be classified as "convective" envelope.
		// Tried to quickly fix this and failed. See: http://gitlab.sr.bham.ac.uk/COMPAS/COMPAS/issues/135
        else if(stellarType == MS_LESS_THAN_07 or stellarType == FIRST_GIANT_BRANCH or stellarType == CORE_HELIUM_BURNING or stellarType == EARLY_ASYMPTOTIC_GIANT_BRANCH or stellarType == THERMALLY_PULSING_ASYMPTOTIC_GIANT_BRANCH or stellarType == NAKED_HELIUM_STAR_HETZSPRUNG_GAP or stellarType == NAKED_HELIUM_STAR_GIANT_BRANCH or stellarType == HELIUM_WHITE_DWARF or stellarType == CARBON_OXYGEN_WHITE_DWARF or stellarType == OXYGEN_NEON_WHITE_DWARF){
            return CONVECTIVE_ENVELOPE;
        }
        else{
            // For BH, NS and massless remmnants
            return REMNANT_ENVELOPE;
        }

    }
    else if(options.commonEnvelopeHertzsprungGapDonor == OPTIMISTIC_HG_CE or options.commonEnvelopeHertzsprungGapDonor == PESSIMISTIC_HG_CE){

        if(stellarType == MS_MORE_THAN_07 or stellarType == NAKED_HELIUM_STAR_MS){
            return RADIATIVE_ENVELOPE;
        }
        else if(stellarType == MS_LESS_THAN_07 or stellarType == HERTZSPRUNG_GAP or stellarType == FIRST_GIANT_BRANCH or stellarType == CORE_HELIUM_BURNING or stellarType == EARLY_ASYMPTOTIC_GIANT_BRANCH or stellarType == THERMALLY_PULSING_ASYMPTOTIC_GIANT_BRANCH or stellarType == NAKED_HELIUM_STAR_HETZSPRUNG_GAP or stellarType == NAKED_HELIUM_STAR_GIANT_BRANCH or stellarType == HELIUM_WHITE_DWARF or stellarType == CARBON_OXYGEN_WHITE_DWARF or stellarType == OXYGEN_NEON_WHITE_DWARF){
            // HG may have radiative envelope. Check that later.
            return CONVECTIVE_ENVELOPE;
        }
        else{
            // For BH, NS and massless remmnants
            return REMNANT_ENVELOPE;
        }
    }
    else{

        std::cerr << m_randomSeed << "\tError in determing envelope type. Shouldn't get here." << std::endl;

        return 0;

    }
}

//
//double calculateZetaThermal(Star starcopy, double percentageMassLoss, bool addMass, const programOptions &options, bool & calculateZetaThermalErrorMassFlag, bool & calculateZetaThermalErrorRadiusFlag){
//    /*
//     Calculate the radius-mass exponent zeta, assuming the star has had time to recover its thermal equilibrium. This is calculated using a fake mass loss step and recomputing the stars radius.
//
//     Parameters
//     ----------
//     starcopy : Star
//     Star to calculate zeta thermal for
//     percentageMassLoss : double
//     Percentage of mass to artificially lose from the star while attempting to calcluate zeta_thermal
//     addMass : bool
//     Whether to add mass or remove it when calculating zeta_thermal
//
//     Returns
//     --------
//     Zeta_thermal : double
//     Radius-Mass exponent Zeta for the thermal timescale
//
//     */
//    bool debugging = false;
//    if(addMass){
//        percentageMassLoss *= -1.0;
//    }
//
//    double zero_time = 0.0;             // We want to evolve the star for zero time, but update its radius
//    bool useFakeTimestep = true;        // We want to use a fake timestep
//
//    double radiusBeforePreEvolve = starcopy.m_Radius;
//    double massBeforePreEvolve = starcopy.m_Mass;
//
//    if(debugging){
//        std::cout << "Radius before pre evolve (zeta_thermal) = " << radiusBeforePreEvolve << std::endl;
//        std::cout << "Mass before pre evolve (zeta_thermal) = " << massBeforePreEvolve << std::endl;
//    }
//    // Allow star to respond to previous mass loss changes
//    starcopy.evolveOneTimestep(zero_time, useFakeTimestep, options);
//    starcopy.evolveOneTimestep(zero_time, useFakeTimestep, options); // what happens if you do this again? Seems to fix some issues, should come back and look at order of updating radius
//
//    // Calculate properties of the star before fake mass loss
//    double radiusBeforeMassLoss = starcopy.m_Radius;
//    double massBeforeMassLoss = starcopy.m_MassPrev;                // Due to order of updating radius bug
//
//    // Initialise some variables for the mass and radius after mass loss
//    double radiusAfterMassLoss = 0.0;
//    double massAfterMassLoss = 0.0;
//
//    // To calculate radius-mass exponent need logarithmns of radii, masses
//    double logRadiusBeforeMassLoss = 0.0;
//    double logMassBeforeMassLoss = 0.0;
//    double logRadiusAfterMassLoss = 0.0;
//    double logMassAfterMassLoss = 0.0;
//
//    // Variables to calculate difference in mass, difference in radius
//    double deltaR = 0.0;
//    double deltaM = 0.0;
//
//    if(debugging){
//        std::cout << "Radius before mass loss (zeta_thermal) = " << radiusBeforeMassLoss << std::endl;
//        std::cout << "Mass before mass loss (zeta_thermal) = " << massBeforeMassLoss << std::endl;
//    }
//    // Check that radius, mass are sensible
//    if(radiusBeforeMassLoss > 0.0){
//        logRadiusBeforeMassLoss = log(radiusBeforeMassLoss);
//    }
//    else{
//		if(calculateZetaThermalErrorRadiusFlag==false){
//            std::cerr << starcopy.m_randomSeed << "\tError in calculateZetaThermal, radius before fake mass loss < 0.0" << std::endl;
//			calculateZetaThermalErrorRadiusFlag = true;
//		}
//    }
//    if(massBeforeMassLoss > 0.0){
//        logMassBeforeMassLoss = log(massBeforeMassLoss);
//    }
//    else{
//		if(calculateZetaThermalErrorMassFlag==false){
//			std::cerr << starcopy.m_randomSeed << "\tError in calculateZetaThermal, mass before fake mass loss < 0.0" << std::endl;
//			calculateZetaThermalErrorMassFlag = true;
//		}
//    }
//
//    // Reduce mass of star by percentageMassLoss and recalculate the radius
//    starcopy.m_Mass = starcopy.m_Mass * (1.0 - (percentageMassLoss/100.0));
//
//    massAfterMassLoss = starcopy.m_Mass;
//
//    if(debugging){std::cout << "Mass after mass loss (zeta_thermal) = " << massAfterMassLoss << std::endl;}
//
//    // Recalculate radius of star
//    starcopy.evolveOneTimestep(zero_time, useFakeTimestep, options);
//
//    radiusAfterMassLoss = starcopy.m_Radius;
//
//    if(debugging)
//        std::cout << "Radius after mass loss (zeta_thermal) = " << radiusAfterMassLoss << std::endl;
//
//    // If mass and radius are both sensible, take the log of them
//    if(radiusAfterMassLoss > 0.0){
//        logRadiusAfterMassLoss = log(radiusAfterMassLoss);
//    }
//    else{
//		if(calculateZetaThermalErrorRadiusFlag==false){
//            std::cerr << starcopy.m_randomSeed << "\tError in calculateZetaThermal, radius after fake mass loss < 0.0" << std::endl;
//			calculateZetaThermalErrorRadiusFlag = true;
//		}
//    }
//
//    if(massAfterMassLoss > 0.0){
//        logMassAfterMassLoss = log(massAfterMassLoss);
//    }
//    else{
//		if(calculateZetaThermalErrorMassFlag==false){
//            std::cerr << starcopy.m_randomSeed << "\tError in calculateZetaThermal, mass after fake mass loss < 0.0" << std::endl;
//			calculateZetaThermalErrorMassFlag = true;
//		}
//    }
//
//    double Zeta_Thermal = (logRadiusAfterMassLoss - logRadiusBeforeMassLoss)/(logMassAfterMassLoss - logMassBeforeMassLoss);
//
//    deltaM = massAfterMassLoss - massBeforeMassLoss;
//    deltaR = radiusAfterMassLoss - radiusBeforeMassLoss;
//
//    if(debugging){
//		// ALEJANDRO - 31/01/2017 - Following statements where displayed as 'std::cerr' when they should be 'std::cout'. Not sure why, but changed it.
//        std::cout << "Final Zeta_thermal from function = " << Zeta_Thermal << std::endl;
//        std::cout << "DeltaR from function = " << deltaR << std::endl;
//        std::cout << "DeltaM from function = " << deltaM << std::endl;
//    }
//
//    // Return zeta
//    return Zeta_Thermal;
//}
//
//double determineConvergedMassStepZetaThermal(const Star & star_to_calculate_zeta_thermal, const programOptions &options, bool & calculateZetaThermalErrorMassFlag, bool & calculateZetaThermalErrorRadiusFlag){
//    /*
//     Calculate the radius response of a star to mass loss on the thermal timescale as characterised by the mass-radius exponent zeta_thermal
//
//     Paramters
//     ----------
//     star_to_calculate_zeta_thermal : Star
//     Star to calculate zeta_thermal for
//
//     Returns
//     --------
//     zeta_thermal : double
//     Mass-radius exponent zeta_thermal = dlnR/dlnM
//
//     */
//    int niter = 0;                                                          // Current mass iteration
//    int niter_max = 10;                                                     // Maximum number of mass step iterations to make before failing loudly
//    bool debugging = false;
//
//    if(star_to_calculate_zeta_thermal.m_stellarType == NEUTRON_STAR or star_to_calculate_zeta_thermal.m_stellarType == BLACK_HOLE){
//        // I guess this is not really right for NS but is correct for BHs.
//        return 1.0;
//    }
//    else if(star_to_calculate_zeta_thermal.m_stellarType == FIRST_GIANT_BRANCH or star_to_calculate_zeta_thermal.m_stellarType == CORE_HELIUM_BURNING or star_to_calculate_zeta_thermal.m_stellarType == EARLY_ASYMPTOTIC_GIANT_BRANCH or star_to_calculate_zeta_thermal.m_stellarType == THERMALLY_PULSING_ASYMPTOTIC_GIANT_BRANCH or star_to_calculate_zeta_thermal.m_stellarType == NAKED_HELIUM_STAR_GIANT_BRANCH or star_to_calculate_zeta_thermal.m_stellarType == NAKED_HELIUM_STAR_HETZSPRUNG_GAP){
//        // As described in BSE paper (Hurley et al 2002) section 2.6.1, the equation is Equaion 47 of Hurley et al 2000
//        // Then modified according to Equation 56 of Hurley et al 2002
//        double x = calculateRadiusGiantBranchConstantX(star_to_calculate_zeta_thermal.m_logMetallicityXi, star_to_calculate_zeta_thermal.m_randomSeed);
//        double exponent = -x + (2.0 * pow((star_to_calculate_zeta_thermal.m_coreMass/star_to_calculate_zeta_thermal.m_Mass), (5.0)));
//
//        // Debugging
//        if(debugging){
//        std::cout << "Zeta_thermal_giant : x = " << x << std::endl;
//        std::cout << "Zeta_thermal_giant : exponent = " << exponent << std::endl;
//        }
//        return exponent;
//    }
//    else{
//
//        double Zeta_Thermal = 0.0;                                              // Stores final value of mass-radius exponent
//        double Zeta_Thermal_Previous = 0.0;                                     // Stores value of mass-radius exponent at previous mass step
//
//        //    double initialMassLossPercentage = 1.0;                                 // Initial percentage of mass to try removing from the star to calculate zeta_thermal
//        double initialMassLossPercentage = 1E-3;            // try initially smaller value?
//        double currentMassLossPercentage = initialMassLossPercentage;           // Current percentage of mass to try removing from the star to calculate zeta_thermal
//        //    double currentMassLossPercentage = 0.0;                                 // Check what happens at end of EAGB for 0 mass loss
//        // Yeah same thing happens, you remove some (0) mass and it collapses to a black hole. Why? This is what the pre-evolution step was supposed to solve?
//
//        double Zeta_Thermal_tolerance = 1E-3;                                   // Tolerance between mass steps
//
//        //    if(star_to_calculate_zeta_thermal.m_stellarType == EARLY_ASYMPTOTIC_GIANT_BRANCH){
//        //        std::cout << "Broken near the end of EAGB?. Constant deltaR and decreasing deltaM" << std::endl;
//        //    }
//
//        if(debugging){
//            std::cout << "Previous, current stellar type (thermal) = " << star_to_calculate_zeta_thermal.m_stellarTypePrev << " " << star_to_calculate_zeta_thermal.m_stellarType << std::endl;
//
//            if(star_to_calculate_zeta_thermal.m_stellarType == HERTZSPRUNG_GAP){
//                std::cout << "HERTZSPRUNG GAP" << std::endl;
//            }
//
//            if(star_to_calculate_zeta_thermal.m_stellarType == CORE_HELIUM_BURNING){
//                std::cout << "CORE HELIUM BURNING" << std::endl;
//            }
//        }
//
//        // Seems to get dR/RM > 0 which means shrinks with mass loss which doesn't seem right
//
//        while(niter < niter_max){
//
//            // Calculate zeta using this functon
//            Zeta_Thermal = calculateZetaThermal(star_to_calculate_zeta_thermal, currentMassLossPercentage, true, options, calculateZetaThermalErrorMassFlag, calculateZetaThermalErrorRadiusFlag);
//
//            if((1.0 - fabs(Zeta_Thermal_Previous/ Zeta_Thermal)) < Zeta_Thermal_tolerance){
//                if(debugging){
//                std::cout << "Zeta_Thermal converged with:" << std::endl;
//                std::cout << "Zeta_Thermal_Previous = " << Zeta_Thermal_Previous << std::endl;
//                std::cout << "Zeta_Thermal = " << Zeta_Thermal << std::endl;
//                std::cout << "test < tol = " << (1.0 - fabs(Zeta_Thermal_Previous/ Zeta_Thermal)) << " < " << Zeta_Thermal_tolerance << std::endl;
//                }
//                break;
//            }
//            else{
//                // Reduce mass step and increment loop counter
//                if(debugging){
//                std::cout << "Zeta_Thermal_Previous = " << Zeta_Thermal_Previous << std::endl;
//                std::cout << "Zeta_Thermal = " << Zeta_Thermal << std::endl;
//                std::cout << "Mass loss percentage = " << currentMassLossPercentage << std::endl;
//                std::cout << "niter = " << niter << std::endl;
//                }
//                currentMassLossPercentage *= 0.1;
//                niter += 1;
//                Zeta_Thermal_Previous = Zeta_Thermal;
//                if(debugging){
//                std::cout << "Zeta_Thermal_Previous = " << Zeta_Thermal_Previous << std::endl;
//                std::cout << "Zeta_Thermal = " << Zeta_Thermal << std::endl;
//                std::cout << "Mass loss percentage = " << currentMassLossPercentage << std::endl;
//                std::cout << "niter = " << niter << std::endl;
//                }
//            }
//        }
//
//        if(niter > niter_max){
//            std::cerr << star_to_calculate_zeta_thermal.m_randomSeed << "\tFailed to converge calculating Zeta_Thermal" << std::endl;
//            if(debugging){
//			// std::cout << "Failed to converge calculating Zeta_Thermal" << std::endl;
//            std::cout << "Zeta_Thermal_Previous = " << Zeta_Thermal_Previous << std::endl;
//            std::cout << "Zeta_Thermal = " << Zeta_Thermal << std::endl;
//            std::cout << "Mass loss percentage = " << currentMassLossPercentage << std::endl;
//            std::cout << "niter = " << niter << std::endl;
//            }
//        }
//
//        return Zeta_Thermal;
//    }
//}
//
//
//double calculateZetaNuclear(Star starcopy, double timestep, const programOptions &options, bool & calculateZetaNuclearErrorRadiusFlag, bool & calculateZetaNuclearErrorAgeFlag){
//    /*
//     Calculate the radius-time exponent zeta, which is derived by the nuclear evolution of the star according to SSE. This is calculated using a fake time step and recomputing the stars radius.
//
//     Parameters
//     ----------
//     starcopy : Star
//     Copy of the star to calculate Zeta Nuclear for
//
//     Returns
//     --------
//     Zeta_nuclear : double
//     Radius-Time exponent Zeta for the nuclear timescale
//
//     */
//
//    // We now want to use fake timesteps to estimate zeta rather than actually evolve the star
//    bool debugging = false;
//    bool useFakeTimestep = false;
//    int niter = 0;                              // Initialise a counter
//    double zero_time = 0.0;                     // Time = 0, not sure if evolution routine will like that
//
//    double Zeta_Nuclear = 0;                    // dlnR/dt | nuclear -- radius-time exponent for star on the nuclear timescale
//
//    double radiusBeforeTimeStep = starcopy.m_Radius;
//    double ageBeforeTimeStep    = starcopy.m_Age;
//
//    double radiusAfterTimeStep = 0.0;
//    double ageAfterTimeStep = 0.0;
//
//    // To calculate radius-time exponent need logarithmns of radii
//    double logRadiusBeforeTimeStep  = 0.0;
//    double logRadiusAfterTimeStep   = 0.0;
//
//    double deltaR = 0.0;
//    double deltaT = 0.0;
//
//    // Check that radius, timestep are sensible
//    if(radiusBeforeTimeStep > 0.0){
//        logRadiusBeforeTimeStep = log(radiusBeforeTimeStep);
//    }
//    else{
//		if(calculateZetaNuclearErrorRadiusFlag==false){
//            std::cerr << starcopy.m_randomSeed << "\tError in calculateZetaNuclear function, radius before fake mass loss < 0.0" << std::endl;
//			calculateZetaNuclearErrorRadiusFlag=true;
//		}
//    }
//    if(ageBeforeTimeStep < 0.0){
//		if(calculateZetaNuclearErrorAgeFlag==false){
//            std::cerr << starcopy.m_randomSeed << "\tError in calculateZetaNuclear function, age < 0.0" << std::endl;
//			calculateZetaNuclearErrorAgeFlag = true;
//		}
//    }
//
//    double initialRadiusBeforeTimeStep = radiusBeforeTimeStep;
//    double initialAgeBeforeTimeStep = ageBeforeTimeStep;
//
//    starcopy.evolveOneTimestep(0.0, useFakeTimestep, options);
//
//    // With the current way of doing things, this star has already lost a bunch of mass due to the mass loss caller, need to update radius of star to respond to this, before we evolve it for a short amount of time assuming no mass loss. When Alejandro moves the massLossCaller function, we will need to do things differently, and this will need fixing.
//    radiusBeforeTimeStep = starcopy.m_Radius;
//    ageBeforeTimeStep = starcopy.m_Age;
//
//    starcopy.evolveOneTimestep(timestep, useFakeTimestep, options);
//
//    radiusAfterTimeStep = starcopy.m_Radius;
//    ageAfterTimeStep = starcopy.m_Age;
//
//    if(radiusAfterTimeStep > 0.0){
//        logRadiusAfterTimeStep = log(radiusAfterTimeStep);
//    }
//    else{
//		if(calculateZetaNuclearErrorRadiusFlag==false){
//			std::cerr << starcopy.m_randomSeed << "\tError in calculateZetaNuclear, radius after timestep < 0.0" << std::endl;
//			calculateZetaNuclearErrorRadiusFlag=true;
//		}
//    }
//
//    if(ageAfterTimeStep < 0.0){
//		if(calculateZetaNuclearErrorAgeFlag==false){
//			std::cerr << starcopy.m_randomSeed << "\tError in calculateZetaNuclear, age after timestep < 0.0" << std::endl;
//			calculateZetaNuclearErrorAgeFlag=true;
//		}
//    }
//
//    deltaR = radiusAfterTimeStep - initialRadiusBeforeTimeStep;
//    deltaT = ageAfterTimeStep - ageBeforeTimeStep;
//
//    Zeta_Nuclear = (logRadiusAfterTimeStep - logRadiusBeforeTimeStep)/(ageAfterTimeStep - ageBeforeTimeStep);
//
//    // Debugging
//    //std::cout << "Radius, age after timestep while attempting to calculate zeta_nuclear = " << starcopy.m_Radius << " " << starcopy.m_Age << std::endl;
//    //std::cout << "DeltaR from function = " << deltaR << std::endl;
//    //std::cout << "DeltaT from function = " << deltaT << std::endl;
//    //std::cout << "Zeta nuclear         = " << Zeta_Nuclear << std::endl;
//
//    // Return zeta
//    return Zeta_Nuclear;
//
//}
//
//double determineConvergedTimestepZetaNuclear(const Star & star_to_calculate_zeta_nuclear, const programOptions &options, bool & calculateZetaNuclearErrorRadiusFlag, bool & calculateZetaNuclearErrorAgeFlag){
//    /*
//     Calculate the radius-time exponent zeta, which is derived by the nuclear evolution of the star according to SSE. This is calculated using a fake time step and recomputing the stars radius.
//
//     Iterate evolving the star for smaller timesteps until the value of Zeta nuclear doesn't change by more than some tolerance
//
//     Parameters
//     -----------
//     star_to_calculate_zeta_nuclear : Star
//     Star to calculate zeta_nuclear for
//
//     Returns
//     --------
//     Zeta_Nuclear : double
//     dlnR/dt | nuc
//
//     */
//	 
//	 bool	debugging = false;
//	 
//    if(star_to_calculate_zeta_nuclear.m_stellarType == BLACK_HOLE){
//        return 0.0;
//    }
//    else{
//        // We now want to use fake timesteps to estimate zeta rather than actually evolve the star
//        int niter = 0;                              // Initialise a counter
//        int niter_max = 10;                         // Maximum number of times to try
//
//        double Zeta_Nuclear = 0;                    // dlnR/dt | nuclear -- radius-time exponent for star on the nuclear timescale
//        double Zeta_Nuclear_Previous = 0.0;         // dlbR/dt | nuclear at previous size timestep
//
//        double timestep_initial = 1E-3;         // Initial timestep to use
//        //    double timestep_initial = 1.0;         // Initial timestep to use
//        double timestep = timestep_initial;     // Current timestep to use
//        double zeta_nuclear_tolerance = 1E-3;   // Fractional difference to tolerate when reducing stepsize calculating zeta_nuclear
//
//        // Check convergence of Zeta nuclear
//        // Want to make it so that the function works out what timestep makes sense to use itself
//
//        while(niter < niter_max){
//
//            // Calculate zeta using this functon
//            Zeta_Nuclear = calculateZetaNuclear(star_to_calculate_zeta_nuclear, timestep, options, calculateZetaNuclearErrorRadiusFlag, calculateZetaNuclearErrorAgeFlag);
//
//            if((1.0 - fabs(Zeta_Nuclear_Previous / Zeta_Nuclear)) < zeta_nuclear_tolerance){
//                //            std::cout << "Zeta_nuclear converged with:" << std::endl;
//                //            std::cout << "Zeta_Nuclear_Previous = " << Zeta_Nuclear_Previous << std::endl;
//                //            std::cout << "Zeta_Nuclear = " << Zeta_Nuclear << std::endl;
//                //            std::cout << "test < tol = " << (1.0 - fabs(Zeta_Nuclear_Previous / Zeta_Nuclear)) << " < " << zeta_nuclear_tolerance << std::endl;
//                break;
//            }
//            else{
//                // Reduce timestep and increment loop counter
//                //            std::cout << "zeta_nuclear_previous = " << Zeta_Nuclear_Previous << std::endl;
//                //            std::cout << "timestep = " << timestep << std::endl;
//                //            std::cout << "niter = " << niter << std::endl;
//                timestep *= 0.1;
//                niter += 1;
//                Zeta_Nuclear_Previous = Zeta_Nuclear;
//                //            std::cout << "zeta_nuclear_previous = " << Zeta_Nuclear_Previous << std::endl;
//                //            std::cout << "timestep = " << timestep << std::endl;
//                //            std::cout << "niter = " << niter << std::endl;
//            }
//        }
//
//        if(niter > niter_max){
//            std::cerr << star_to_calculate_zeta_nuclear.m_randomSeed << "\tFailed to converge calculating Zeta_Nuclear" << std::endl;
//			if(debugging){
//				std::cout << "zeta_nuclear_previous = " << Zeta_Nuclear_Previous << std::endl;
//				std::cout << "timestep = " << timestep << std::endl;
//				std::cout << "niter = " << niter << std::endl;
//			}
//        }
//
//        // Return zeta
//        return Zeta_Nuclear;
//    }
//}
