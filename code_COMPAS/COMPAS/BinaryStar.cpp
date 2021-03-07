//
//  BinaryStar.cpp

#include "BinaryStar.h"
#include <fstream>

// Initialise s_nIDGenerator
int BinaryStar::s_nIDGenerator = 0;

//bool BinaryStar::s_firstSupernova = true;
//std::ofstream BinaryStar::s_supernovaeOutput("./supernovae.txt",std::ios_base::app);

//bool BinaryStar::s_commonEnvelope = true;
//std::ofstream BinaryStar::s_commonEnvelopeOutput("./commonEnvelopes.txt");

bool BinaryStar::s_firstBinary = true;
std::ofstream BinaryStar::s_runtimeOutput("./runtimes.txt");       // Static as used for all binaries

// Default Constructor -- shouldn't be needed.
BinaryStar::BinaryStar(){

    //std::cout << "Calling default BinaryStar constructor 1)" << std::endl;

    // Initial Binary Parameters -- kept constant as a record of the initial parameters of the binary
    m_SemiMajorAxisInitial = 0.0;
    m_EccentricityInitial = 0.0;

    // Initialise variables to hold parameters prior to 2nd supernova
    m_SemiMajorAxisPre2ndSN = 0.0;
    m_EccentricityPre2ndSN = 0.0;

    // Initialise variables to hold parameters at DCO formation
    m_EccentricityAtDCOFormation = 0.0;
    m_SemiMajorAxisAtDCOFormation = 0.0;

    // Initialise runtime
    m_runtime = 0.0;
    m_clockStart = clock();

    // Current Binary parameters
    // CPLB: Specify units!
    //m_ID                = s_nIDGenerator++;                 // Generate unique ID
    m_Mass1             = 1;                                // Default value - can't initialise to 0 since q = m2/m1 (in Msol)
    m_Mass2             = 1;                                // Default value - can't initialise to 0 since q = m2/m1 (in Msol)
    m_MassRatio         = m_Mass1 / m_Mass2;                // Calculated from above
    m_SemiMajorAxis     = 1;                                // Default value (in AU)
    m_SemiMajorAxisPrev = m_SemiMajorAxis;
    m_SemiMajorAxisPrime= m_SemiMajorAxis;
    m_RadiusRL1         = 1;
    m_RadiusRLprev1     = 0;
    m_RadiusRL2         = 1;
    m_RadiusRLprev2     = 0;
    m_Eccentricity      = 0;                                // Default value
    m_EccentricityPrev = 0.0;
    m_EccentricityPrime = 0.0;
    m_Inclinaton        = 0;                                // Initialise to 0
    m_delayTime         = 0;                                // Initialise to 0 (UNITS?)
    m_Metallicity1      = Zsol;                             // Default solar metallicity
    m_Metallicity2      = Zsol;                             // Default solar metallicity

    // Conserved quantities
    m_TotalEnergy = 1;                                  // Default total energy
    m_TotalAngularMomentum = 1;                         // Default total angular momentum
    m_TotalEnergyPrime = 1;
    m_TotalAngularMomentumPrime = 1;
    m_TotalAngularMomentumPrev = 1;

	// Quantities for calculating orbital energy and angular momentum
	m_totalMass 		= 1;
	m_totalMassPrime 	= 1;
	m_totalMassPrev		= 1;
	m_reducedMass 		= 1;
	m_reducedMassPrime	= 1;
	m_reducedMassPrev	= 1;
	m_totalOrbitalEnergyPrime = 1;
	m_totalOrbitalEnergyPrev = 1;
	m_totalOrbitalAngularMomentumPrime = 1;
	m_totalOrbitalAngularMomentumPrev = 1;

    m_time = 0.0;
	m_dt   = 0.0;
    m_timePrev = 0.0;
    m_flagFirstTimeStep = true;

    m_aTidesDiff = 0.0;
    m_omegaTidesDiff = 0.0;
    m_omegaTides = 0.0;

    m_aMagneticBrakingDiff = 0.0;
    m_angularMomentumMagneticBrakingDiff = 0.0;

    m_aMassLossDiff = 0.0;
    m_omegaMassLossDiff = 0.0;

    m_aMassTransferDiff = 0.0;
    m_omegaMassTransferDiff = 0.0;

    massTransferTracker = 0.0;
    rocheLobeTracker1 = 0.0;
    rocheLobeTracker2 = 0.0;
	m_massTransferTrackerHistory = NO_MASS_TRANSFER;

    m_massTransferFlag = false;
    m_massTransferZeroFlag=false;
    m_initiateMassTransfer = false;
    m_insideMassTransfer = false;
    m_fastPhaseCaseA=false;
	m_aRocheLobeOverFlowErrorPrinting = false;
	m_aFastPhaseRocheLobeOverFlowErrorPrinting = false;
    m_dynamicalFlag = false;
    m_massTransferType = 0;
    m_massTransferRateThermal = 0;
    m_massTransferThermalFlag = false;
	m_massTransferFractionAccreted = 1.0;
	m_massTransferJloss = 1.0;
	m_stableRLOFafterCEE = false;

	// Mass transfer error flags
	calculateZetaThermalErrorMassFlag = false;
	calculateZetaThermalErrorRadiusFlag = false;
	calculateZetaThermalErrorAgeFlag = false;
	calculateZetaNuclearErrorMassFlag = false;
	calculateZetaNuclearErrorRadiusFlag = false;
	calculateZetaNuclearErrorAgeFlag = false;
	adaptiveRocheLobeOverFlowErrorMassFlag = false;
	adaptiveRocheLobeOverFlowErrorRadiusFlag = false;
	adaptiveRocheLobeOverFlowErrorAgeFlag = false;
	massTransferFastPhaseCaseAErrorMassFlag = false;
	massTransferFastPhaseCaseAErrorRadiusFlag = false;
	massTransferFastPhaseCaseAErrorAgeFlag = false;
    
	// Common Envelope
    m_commonEnvelopeFlag = false;
    m_doubleCoreCommonEnvelopeFlag = false;
    m_stellarMerger = false;
	m_commonEnvelopeOccuredAtLeastOnce = false;
    m_RLOFPrimaryAfterCEE = false;
	m_RLOFSecondaryAfterCEE = false;
	m_optimisticCommonEnvelopeFlag = false;
	m_EccentricityPreCEE = NEVER_SET;
	m_EccentricityPostCEE = NEVER_SET;
	m_SemiMajorAxisPreCEE = NEVER_SET;
	m_SemiMajorAxisPostCEE = NEVER_SET;
	m_rocheLobe1to2PreCEE = NEVER_SET;
	m_rocheLobe1to2PostCEE = NEVER_SET;
	m_rocheLobe2to1PreCEE = NEVER_SET;
	m_rocheLobe2to1PostCEE = NEVER_SET;
    m_counterCEE = 0;
    m_immediateRLOFAfterCEE = false;
    m_simultaneousRLOFleadingToCEEFlag = false;
    m_mainSequenceAccretorDuringCEEFlag = false;
    m_EccentricityRLOF = NEVER_SET;

    //Addition by Coen 18-10-2017 for zeta study
    m_zetaRLOFanalytic	= NEVER_SET;
    m_zetaRLOFnumerical	= NEVER_SET;
	//Alejandro - 07/08/2018 for CE study
	m_zetaStarCompare	= NEVER_SET;
	
    // X-ray binary
    m_is_xray_binary = false;

    // Spins
    m_Chi1          = 0;                                // Default dimensionless magnitude for spin 1
    m_Chi2          = 0;                                // Default dimensionless magnitude for spin 2 (should only get magnitude after SN!)

    // Need to make it so binary star contains two instances of star, star1, star2, to hold masses, radii and luminosities.
    // star 1 initially more massive
    star1 = Star(m_Mass1, m_Metallicity1);
    star2 = Star(m_Mass2, m_Metallicity2);
    
    // Misalignments
    m_theta1_i      = 0;                                // Initial misalignment of star 1
    m_theta2_i      = 0;                                // Initial misalignment of star 2
    m_theta1        = 0;                                // By default, aligned
    m_theta2        = 0;                                // By default, aligned
    m_theta12       = 0;                                // By default, aligned
    m_deltaPhi      = 0;                                // By default, aligned
    m_deltaPhi_i    = 0;                                // By default, aligned

    m_S1x           = 0;                                // x component of spin1
    m_S1y           = 0;                                // y component of spin1
    m_S1z           = 0;                                // z component of spin1

    m_S2x           = 0;                                // x component of spin2
    m_S2y           = 0;                                // y component of spin2
    m_S2z           = 0;                                // z component of spin2

    m_nu            = 0;                                // Initialise orbital velocity
    m_iota          = 0;

    // Initialise other parameters to 0
    m_MSN           = 0;
    m_MSNPrime      = 0;
    m_MC            = 0;
    m_MCPrime       = 0;

    m_Theta         = 0;
    m_Phi           = 0;
    m_Psi           = 0;
    m_vRel          = 0;
    m_uK            = 0;
    m_Radius        = 0;
    //JIM BARRETT - 05/12/2016 - removing m_aPrime as it isn't used at all
    //m_aPrime        = 0;
    m_EPrime        = 0;
    m_ePrime        = 0;
    m_cosiPrime     = 0;
    m_iPrime        = 0;
    m_tc            = 0;
    m_beta          = 0;                // Angle between r and v, related to eccentricity (= pi/2 for circular e = 0)
    m_Survived      = true;
    m_Merged        = false;
    m_errorFlag     = false;
    m_mergesInHubbleTimeFlag = false;

    m_systemicVelocity = NEVER_SET;                  // Post supernova systemic velocity

    m_CHEvolutionFlag = false;
	//kIT x-RAY-BINARY parameters
    // m_HMXB          = false;
    // m_HMXBtime      = 0;//kit
    // m_XRayPower     = 0;//
    // m_XRayAveMass1  = 0;
    // m_XRayAveMass2  = 0;
    // m_XRayAveSep    = 0;
    // m_XRayBinaryType= 0;

    // redundant variables
    m_Mass1_Keep    = 0;
    m_Mass2_Keep    = 0;

    eventCounter = 0;
    same_RLOF_loop = false;
    current_stellar_pair =0;
    mt_primary_counter = 0;
    mt_secondary_counter =0;
    mt_primary_ep1 = 0;
    mt_primary_ep1_K1 = 0;
    mt_primary_ep1_K2 = 0;
    mt_primary_ep2 = 0;
    mt_primary_ep2_K1 = 0;
    mt_primary_ep2_K2 = 0;
    mt_primary_ep3 = 0;
    mt_primary_ep3_K1 = 0;
    mt_primary_ep3_K2 = 0;
    mt_secondary_ep1 =0;
    mt_secondary_ep1_K1 =0;
    mt_secondary_ep1_K2 =0;
    mt_secondary_ep2 =0;
    mt_secondary_ep2_K1 =0;
    mt_secondary_ep2_K2 =0;
    mt_secondary_ep3 =0;
    mt_secondary_ep3_K1 =0;
    mt_secondary_ep3_K2 =0;
    SN_primary_type_1 =0;
    SN_primary_type_2 =0;
    SN_primary_type_3 =0;
    SN_secondary_type_1 =0;
    SN_secondary_type_2 =0;
    SN_secondary_type_3 =0;
    CEE =0;
    CEE_instigator =0;
    CEE_failed =0;
    CEE_failed_instigator =0;
    CEE_wet =0;
    CEE_wet_instigator =0;
    stellar_type_K1 = 0;
    stellar_type_K2 = 0;
	binary_disbound = false;
	
    m_systemicVelocity = NEVER_SET;                  // Post supernova systemic velocity

	m_synchronizationTimescale = NEVER_SET;
	m_circularizationTimescale = NEVER_SET;

    weight = 1;
    samplingphase = -1;

}

// Constructor where binary is generated according to distributions specified in options
BinaryStar::BinaryStar(gsl_rng *r, programOptions &options, AISvariables &aisvariables){     

    //std::cout << "Calling default BinaryStar constructor 2)" << std::endl;

    // Initialise runtime
    m_runtime = 0.0;
    m_clockStart = clock();

    m_ID          = s_nIDGenerator++;                       // Unique ID
    m_randomSeed  = options.randomSeed + m_ID;              // Set the random seed

    // I want each binary to be reproducible; therefore I reseed the random number generator here -- this will also ensure that each binary gets a unique random seed
    if(options.fixedRandomSeed){

//        std::cout << "Reseeding binary with : " << m_randomSeed << std::endl;
        gsl_rng_set(r, m_randomSeed);                           // Reseed the random number generator

    }
    else{
//         This is for using a truly random random seed (based on system time)
        if(options.populationDataPrinting){
            std::cout << "Using default random seed " << gsl_rng_default_seed << std::endl;
        }
        m_randomSeed = gsl_rng_default_seed + m_ID;
        gsl_rng_set(r, m_randomSeed);
    }

    if(options.individualSystem){

        m_Mass1 = options.primaryMass;                             // Initial primary mass (in Msol)
        m_Mass2 = options.secondaryMass;                           // Initial secondary mass (in Msol)

        m_MassRatio = m_Mass2 / m_Mass1;                                  // Initial mass ratio m2/m1

        // Check if user specified separation, orbital period or both (in which case we print a warning and ignore period)
        if(options.binarySeparation > 0.0 and options.binaryOrbitalPeriod <= 0.0){
            m_SemiMajorAxis = options.binarySeparation;                      // Initial semi-major axis in AU
        }
        else if(options.binarySeparation <= 0.0 and options.binaryOrbitalPeriod > 0.0){
            m_SemiMajorAxis = convertPeriodInDaysToSemiMajorAxisInAU(m_Mass1, m_Mass2, options.binaryOrbitalPeriod);
        }
        else if(options.binarySeparation > 0.0 and options.binaryOrbitalPeriod > 0.0){
            std::cout << "Warning, both separation and period specified. Using separation by default" << std::endl;
            m_SemiMajorAxis = options.binarySeparation;
        }
        else{
            std::cerr << m_randomSeed <<  "\tError setting individual system semi major axis" << std::endl;
            m_SemiMajorAxis = 0.0;
        }

        m_SemiMajorAxisPrime = m_SemiMajorAxis;                             // Variable for semi-major axis
        m_SemiMajorAxisPrev  = m_SemiMajorAxis;                             // Variable for semi-major axis at previous timestep
        m_RadiusRL1          = rocheLobeRadius(m_Mass1, m_Mass2);           // Primary dimensionless roche lobe radius
        m_RadiusRLprev1      = m_RadiusRL1;                                 // Primary dimensionless roche lobe radius at previous timestep
        m_RadiusRL2          = rocheLobeRadius(m_Mass2, m_Mass1);           // Secondary dimensionless roche lobe radius
        m_RadiusRLprev2      = m_RadiusRL2;                                 // Secondary dimensionless roche lobe radius at previous timestep
        m_Eccentricity       = options.binaryEccentricity;                  // Initial eccentricity
        m_EccentricityPrime  = m_Eccentricity;                              // Initial eccentricity
        m_Inclinaton         = 0;                                           // Initial spin-tilt inclination is 0
        m_Metallicity1       = options.initialPrimaryMetallicity;           // Metallicity for first star
        m_Metallicity2       = options.initialSecondaryMetallicity;         // Metallicity for second star
        // m_totalAngularMomentum = 0.0;                                    // Initial total angular momentum

        // Initial Binary Parameters -- kept constant as a record of the initial parameters of the binary
        m_SemiMajorAxisInitial = m_SemiMajorAxis;
        m_EccentricityInitial  = m_Eccentricity;

        // Initialise variables to hold parameters prior to 2nd supernova
        m_SemiMajorAxisPre2ndSN = 0.0;
        m_EccentricityPre2ndSN = 0.0;

        // Initialise variables to hold parameters at DCO formation
        m_EccentricityAtDCOFormation = 0.0;
        m_SemiMajorAxisAtDCOFormation = 0.0;

        // Sampled parameters are class members, so also need to set them here
        //m_kickVelocityDistributionSigma = options.kickVelocityDistributionSigma;

        //double m_kickVelocityDistributionSigma;           // Kick velocity sigma in km s^-1 (default = "250" )
        m_kickVelocityDistributionSigmaCCSN_NS = options.kickVelocityDistributionSigmaCCSN_NS;      // Kick velocity sigma in km s^-1 for neutron stars (default = "250" )
        m_kickVelocityDistributionSigmaCCSN_BH = options.kickVelocityDistributionSigmaCCSN_BH;      // Kick velocity sigma in km s^-1 for black holes (default = "250" )
        m_kickVelocityDistributionSigmaForECSN = options.kickVelocityDistributionSigmaForECSN;      // Kick velocity sigma for ECSN in km s^-1 (default = "0" )
        m_kickVelocityDistributionSigmaForUSSN = options.kickVelocityDistributionSigmaForUSSN;      // Kick velocity sigma for USSN in km s^-1 (default = "20" )
        m_kickVelocityDistributionMaximum = options.kickVelocityDistributionMaximum;           // Maximum kick velocity to draw. If negative, no maximum

        m_kickDirectionPower = options.kickDirectionPower;
        m_commonEnvelopeAlpha = options.commonEnvelopeAlpha;
        m_wolfRayetFactor = options.wolfRayetFactor;
        m_luminousBlueVariableFactor = options.luminousBlueVariableFactor;

        // Binary star contains two instances of star, star1, star2, to hold masses, radii and luminosities.
        // star 1 initially more massive
        star1 = Star(m_Mass1, m_Metallicity1, options, r);
        star2 = Star(m_Mass2, m_Metallicity2, options, r);

    }
    else{
        // Simon : Sampling over NS and BH kick distiributions seperately is not yet implemented.
    	// if (options.sampleKickVelocitySigma){
    	// 	double diff = options.sampleKickVelocitySigmaMax - options.sampleKickVelocitySigmaMin;
    	// 	double rand = (gsl_rng_uniform(r) * diff) + options.sampleKickVelocitySigmaMin;
    	// 	m_kickVelocityDistributionSigma = rand;
    	// }
    	// else{
        //  		m_kickVelocityDistributionSigma = options.kickVelocityDistributionSigma;
    	// }

        m_kickVelocityDistributionSigmaCCSN_NS = options.kickVelocityDistributionSigmaCCSN_NS;
        m_kickVelocityDistributionSigmaCCSN_BH = options.kickVelocityDistributionSigmaCCSN_BH;

    	if (options.sampleKickDirectionPower){
    		double diff = options.sampleKickDirectionPowerMax - options.sampleKickDirectionPowerMin;
    		double rand = (gsl_rng_uniform(r) * diff) + options.sampleKickDirectionPowerMin;
    		m_kickDirectionPower = rand;
    	}
        else{
          	m_kickDirectionPower = options.kickDirectionPower;
        }

    	if (options.sampleCommonEnvelopeAlpha){
    		double diff = options.sampleCommonEnvelopeAlphaMax - options.sampleCommonEnvelopeAlphaMin;
    		double rand = (gsl_rng_uniform(r) * diff) + options.sampleCommonEnvelopeAlphaMin;
    		m_commonEnvelopeAlpha = rand;
    	}
    	else{
      		m_commonEnvelopeAlpha = options.commonEnvelopeAlpha;
    	}

    	if (options.sampleWolfRayetMultiplier){
    		double diff = options.sampleWolfRayetMultiplierMax - options.sampleWolfRayetMultiplierMin;
    		double rand = (gsl_rng_uniform(r) * diff) + options.sampleWolfRayetMultiplierMin;
    		m_wolfRayetFactor = rand;
    	}
    	else{
      		m_wolfRayetFactor = options.wolfRayetFactor;
    	}

    	if (options.sampleLuminousBlueVariableMultiplier){
    		double diff = options.sampleLuminousBlueVariableMultiplierMax - options.sampleLuminousBlueVariableMultiplierMin;
    		double rand = (gsl_rng_uniform(r) * diff) + options.sampleLuminousBlueVariableMultiplierMin;
    		m_luminousBlueVariableFactor = rand;
    	}
    	else{
      		m_luminousBlueVariableFactor = options.luminousBlueVariableFactor;
    	}

        m_binaryRLOFonZAMS = true;                  // Whether binary fills its Roche Lobe at ZAMS
        m_secondarySmallerThanMiniumMass = true;    // Whether m2 < m2min

        // Floor 01/05/2018 - the bool initialParametersOutsideParameterSpace is used to check if we need to resample, 
        // when using Adaptive Importance Sampling due to sampling outside of the parameter space. 
        bool initialParametersOutsideParameterSpace = true;
        // Random drawn Gaussian that we will sample from in Adaptive Importance Sampling


        // Generate initial properties of binary
        // We check that it does not overflow its Roche Lobe at ZAMS
        // We also check m2 > m2min
        // using STROOPWAFEL we also check whether it is inside the parameter space
        // draw new random binary using same randomSeed when one of above does not hold. 
        while(m_binaryRLOFonZAMS or m_secondarySmallerThanMiniumMass or initialParametersOutsideParameterSpace){
            
            // If AISrefinementPhase = true we will run AIS step 2 and sample from importance sampling distribution 
            if(options.AISrefinementPhase){
                // the next block decides whether to draw from gaussians or from priors. since FractionSamplingPriors% is always still prior sampled                
                // double rnd11 = gsl_rng_uniform(r);  // Draw a random number between 0 and 1
                aisvariables.DrawingFromAISdistributions = true;   //  we are drawing from AIS distributions in AIS phase 2 
                    // draw a random gaussian for Adaptive importance Sampling       
                int max_rnd = aisvariables.cov_loga.size()-1;  // the max int (i.e. index of gaussian) we want to sample. -1 since cov_loga also read in empty line 
                aisvariables.RandomGaussianDraw = gsl_rng_uniform_int(r, max_rnd); // draw a random integer between 0 and max_rnd-1  
            }
            
            m_Mass1                 = initialMassDistribution(options, r, aisvariables);            // Initial primary mass (in Msol)
            m_MassRatio             = qDistribution(options, r, aisvariables);                      // Should be between 0 or 1
            m_Mass2                 = m_MassRatio*m_Mass1;                                          // Should be less than or equal to m_Mass1 (in Msol)
            m_SemiMajorAxis         = aDistribution(options, r, m_Mass1, m_Mass2, aisvariables);    // Initial separation (in AU)
            m_Eccentricity          = eDistribution(options, r, m_randomSeed);                      // Initial eccentricity
            m_Metallicity1          = metallicityDistribution(options, r);                          // Metallicity for first star
            m_Metallicity2          = metallicityDistribution(options, r);                          // Metallicity for second star

            // Binary star contains two instances of star, star1, star2, to hold masses, radii and luminosities.
            // star 1 initially more massive
            star1 = Star(m_Mass1, m_Metallicity1, options, r);
            star2 = Star(m_Mass2, m_Metallicity2, options, r);
            // Calculate ZAMS rocheLobeRadius
            m_RadiusRL1             = rocheLobeRadius(m_Mass1, m_Mass2);
            m_RadiusRL2             = rocheLobeRadius(m_Mass2, m_Mass1);

            rocheLobeTracker1 = (star1.m_Radius*RsolToAU)/(m_SemiMajorAxis*(1.0 - m_Eccentricity)*m_RadiusRL1);
            rocheLobeTracker2 = (star2.m_Radius*RsolToAU)/(m_SemiMajorAxis*(1.0 - m_Eccentricity)*m_RadiusRL2);

            if(rocheLobeTracker1 > 1.0 or rocheLobeTracker2 > 1.0){
                m_binaryRLOFonZAMS = true;   
                }        
            else{
                m_binaryRLOFonZAMS = false;
            }

            // Check if secondary mass is greater than user specified minimum
            if(m_Mass2 > options.minimumMassSecondary){
                m_secondarySmallerThanMiniumMass = false;
            }
            else if(m_Mass2 < options.minimumMassSecondary){
                m_secondarySmallerThanMiniumMass = true;             
            }
            // When using Adaptive Importance Sampling (step 2) check if drawns from Gaussians are inside the COMPAS parameter space
            if(options.AISrefinementPhase)
            {
                // Resample if m_Mass1 is outside parameter space
                if (m_Mass1 < options.initialMassFunctionMin or m_Mass1 > options.initialMassFunctionMax){
                    initialParametersOutsideParameterSpace = true;  
                }
                // Resample if m_MassRatio is outside parameter space
                else if (m_MassRatio < options.massRatioDistributionMin or m_MassRatio > options.massRatioDistributionMax){
                    initialParametersOutsideParameterSpace = true;
                }
                // Resample if m_SemiMajorAxis is outside parameter space
                else if (m_SemiMajorAxis < options.semiMajorAxisDistributionMin or m_SemiMajorAxis > options.semiMajorAxisDistributionMax){
                    initialParametersOutsideParameterSpace = true;
                }
                // Otherwise continue
                else{
                    initialParametersOutsideParameterSpace = false;
                }
            }
            // If not using Adaptive Importance Sampling the drawn samples should by definition never be outside the parameter space
            else{
                initialParametersOutsideParameterSpace = false;
            }
            // print the randomSeed when the binary is rejected if requested (see RejectedSamplesPrinting in pythonSubmit)
            if(options.RejectedSamplesPrinting){
            	if (m_binaryRLOFonZAMS or m_secondarySmallerThanMiniumMass or initialParametersOutsideParameterSpace){
            		printRejectedSamplesParameters(options.outputPath, "rejectedSamples.txt", m_randomSeed);
            	}
            }
        }	
			
        

        
        // Initialise other parameters
        m_RadiusRLprev1         = m_RadiusRL1;
        m_RadiusRLprev2         = m_RadiusRL2;

        //    m_RadiusRLprev1     = 0.0;
        //    m_RadiusRLprev2     = 0.0;

        m_EccentricityPrime     = m_Eccentricity;                                   // Initial eccentricity
        m_Inclinaton            = 0;                                                // Initial spin-tilt inclination is 0

        m_SemiMajorAxisPrime    = m_SemiMajorAxis;
        m_SemiMajorAxisPrev     = m_SemiMajorAxis;

    //    m_totalAngularMomentum = 0.0;				// Initial total angular momentum

        // Initial Binary Parameters -- kept constant as a record of the initial parameters of the binary
        m_SemiMajorAxisInitial = m_SemiMajorAxis;
        m_EccentricityInitial  = m_Eccentricity;

        // Initialise variables to hold parameters prior to 2nd supernova
        m_SemiMajorAxisPre2ndSN = 0.0;
        m_EccentricityPre2ndSN = 0.0;

        // Initialise variables to hold parameters at DCO formation
        m_EccentricityAtDCOFormation = 0.0;
        m_SemiMajorAxisAtDCOFormation = 0.0;

        m_kickVelocityDistributionSigmaCCSN_NS = options.kickVelocityDistributionSigmaCCSN_NS;      // Kick velocity sigma in km s^-1 for neutron stars (default = "250" )
        m_kickVelocityDistributionSigmaCCSN_BH = options.kickVelocityDistributionSigmaCCSN_BH;      // Kick velocity sigma in km s^-1 for black holes (default = "250" )
        m_kickVelocityDistributionSigmaForECSN = options.kickVelocityDistributionSigmaForECSN;      // Kick velocity sigma for ECSN in km s^-1 (default = "0" )
        m_kickVelocityDistributionSigmaForUSSN = options.kickVelocityDistributionSigmaForUSSN;      // Kick velocity sigma for USSN in km s^-1 (default = "20" )
        m_kickVelocityDistributionMaximum = options.kickVelocityDistributionMaximum;           // Maximum kick velocity to draw. If negative, no maximum

        m_kickDirectionPower = options.kickDirectionPower;
        m_commonEnvelopeAlpha = options.commonEnvelopeAlpha;
        m_wolfRayetFactor = options.wolfRayetFactor;
        m_luminousBlueVariableFactor = options.luminousBlueVariableFactor;

    }

    m_orbitalVelocity	= sqrt(G1*(m_Mass1+m_Mass2)/pow(m_SemiMajorAxis,3.0));		// Initial orbital velocity according to previous parameters in yr-1 units, where [G] = 4*pi^2 AU^3 yr-2 Msol-1 
    m_orbitalVelocityPrime = m_orbitalVelocity;
    m_orbitalVelocityPrev = m_orbitalVelocity;

    // Conserved quantities
    m_TotalEnergy = calculateTotalEnergy(m_SemiMajorAxis,m_Mass1,m_Mass2,star1.m_RZAMS,star2.m_RZAMS,star1.m_omega,star2.m_omega,m_orbitalVelocity,k_definition(star1.m_RZAMS, star1.m_Radius, star1.m_Mass, star1.m_stellarType),k_definition(star2.m_RZAMS, star2.m_Radius, star2.m_Mass, star2.m_stellarType));    // Default total energy
    m_TotalAngularMomentum = calculateAngularMomentum(m_SemiMajorAxis,m_Eccentricity,m_Mass1,m_Mass2,star1.m_RZAMS,star2.m_RZAMS,star1.m_omega,star2.m_omega,m_orbitalVelocity,k_definition(star1.m_RZAMS, star1.m_Radius, star1.m_Mass, star1.m_stellarType),k_definition(star2.m_RZAMS, star2.m_Radius, star2.m_Mass, star2.m_stellarType));   // Default total angular momentum

	// Quantities for calculating orbital energy and angular momentum
	m_totalMass 						= m_Mass1 + m_Mass2;
	m_totalMassPrime 					= m_totalMass;
	m_totalMassPrev						= m_totalMass;
	m_reducedMass 						= (m_Mass1*m_Mass2)/m_totalMass;
	m_reducedMassPrime					= m_reducedMass;
	m_reducedMassPrev					= m_reducedMass;
	m_totalOrbitalEnergy 				= orbitalEnergy(m_reducedMass, m_totalMass, m_SemiMajorAxis);
	m_totalOrbitalEnergyPrime 			= m_totalOrbitalEnergy;
	m_totalOrbitalEnergyPrev 			= m_totalOrbitalEnergy;
	m_totalOrbitalAngularMomentum 		= orbitalAngularMomentum(m_reducedMass, m_totalMass, m_SemiMajorAxis);
	m_totalOrbitalAngularMomentumPrime 	= m_totalOrbitalAngularMomentum;
	m_totalOrbitalAngularMomentumPrev 	= m_totalOrbitalAngularMomentum;

    m_time = 0.0;
	m_dt   = 0.0;
    m_timePrev = 0.0;
    m_flagFirstTimeStep = true;

    // Differential quantities
    m_aTidesDiff = 0.0;
    m_omegaTidesDiff = 0.0;
    m_omegaTides = 0.0;

    m_aMagneticBrakingDiff = 0.0;
    m_angularMomentumMagneticBrakingDiff = 0.0;

    m_aMassLossDiff = 0.0;
    m_omegaMassLossDiff = 0.0;

    m_aMassTransferDiff = 0.0;
    m_omegaMassTransferDiff = 0.0;

    massTransferTracker = 0.0;
    rocheLobeTracker1 = 0.0;
    rocheLobeTracker2 = 0.0;
	m_massTransferTrackerHistory = NO_MASS_TRANSFER;

    m_massTransferFlag = false;
    m_massTransferZeroFlag=false;
    m_fastPhaseCaseA=false;
	m_aRocheLobeOverFlowErrorPrinting = false;
	m_aFastPhaseRocheLobeOverFlowErrorPrinting = false;
    m_initiateMassTransfer = false;
    m_insideMassTransfer = false;
    m_dynamicalFlag = false;
    m_massTransferType = 0;
    m_massTransferRateThermal = 0;
    m_massTransferThermalFlag = false;
	m_massTransferFractionAccreted = options.massTransferFractionAccreted;
	m_massTransferJloss = options.massTransferJloss;
	m_stableRLOFafterCEE = false;

	// Mass transfer error flags
	calculateZetaThermalErrorMassFlag = false;
	calculateZetaThermalErrorRadiusFlag = false;
	calculateZetaThermalErrorAgeFlag = false;
	calculateZetaNuclearErrorMassFlag = false;
	calculateZetaNuclearErrorRadiusFlag = false;
	calculateZetaNuclearErrorAgeFlag = false;
	adaptiveRocheLobeOverFlowErrorMassFlag = false;
	adaptiveRocheLobeOverFlowErrorRadiusFlag = false;
	adaptiveRocheLobeOverFlowErrorAgeFlag = false;
	massTransferFastPhaseCaseAErrorMassFlag = false;
	massTransferFastPhaseCaseAErrorRadiusFlag = false;
	massTransferFastPhaseCaseAErrorAgeFlag = false;

    // Common Envelope
    m_commonEnvelopeFlag = false;
    m_doubleCoreCommonEnvelopeFlag = false;
    m_stellarMerger = false;
	m_commonEnvelopeOccuredAtLeastOnce = false;
	m_RLOFPrimaryAfterCEE = false;
    m_RLOFSecondaryAfterCEE = false;
	m_optimisticCommonEnvelopeFlag = false;
	
	m_EccentricityPreCEE = NEVER_SET;
	m_EccentricityPostCEE = NEVER_SET;
	m_SemiMajorAxisPreCEE = NEVER_SET;
	m_SemiMajorAxisPostCEE = NEVER_SET;
	m_rocheLobe1to2PreCEE = NEVER_SET;
	m_rocheLobe1to2PostCEE = NEVER_SET;
	m_rocheLobe2to1PreCEE = NEVER_SET;
	m_rocheLobe2to1PostCEE = NEVER_SET;
    m_counterCEE = 0;
    m_immediateRLOFAfterCEE = false;
    m_simultaneousRLOFleadingToCEEFlag = false;
    m_mainSequenceAccretorDuringCEEFlag = false;
    m_EccentricityRLOF = NEVER_SET;

    //Addition by Coen 18-10-2017 for zeta study
    m_zetaRLOFanalytic= NEVER_SET;
    m_zetaRLOFnumerical= NEVER_SET;
	//Alejandro - 07/08/2018 for CE study
	m_zetaStarCompare	= NEVER_SET;
	
    // X-ray binary
    m_is_xray_binary = false;

    // Energy quantities (dynamical)
    m_TotalEnergyPrime = m_TotalEnergy;
    m_TotalAngularMomentumPrime = m_TotalAngularMomentum;
    m_TotalAngularMomentumPrev = m_TotalAngularMomentum;

    // Even though these are "binary" quantities, they are different for each star for each supernova, so we set them as properties of the star
    // // Orbital anomalies
    // m_MeanAnomaly       = 0;                 // Mean anomaly at instantaneous time of the SN - uniform in [0, 2pi]
    // m_EccentricAnomaly  = 0;                 // Eccentric anomaly at instataneous time of the SN
    // m_TrueAnomaly       = 0;                 // True anomaly at instantaneous time of the SN

    // Spin variables
    m_Chi1             = spinDistribution(options, r);         // Spin of star 1
    m_Chi2             = spinDistribution(options, r);         // Spin of star 2 (shouldn't get spin until supernova!)

    // Misalignments
    m_theta1_i          = 0;                                    // Initial misalignment of star 1
    m_theta2_i          = 0;                                    // Initial misalignment of star 2
    m_theta1            = 0;                                    // By default, aligned
    m_theta2            = 0;                                    // By default, aligned
    m_theta12           = 0;                                    // By default, aligned
    m_deltaPhi          = 0;                                    // By default, aligned - could it somehow be this causing the prob?
    m_deltaPhi_i        = 0;                                    // DEBUGGING

    m_S1x               = 0;                                    // x component of spin1
    m_S1y               = 0;                                    // y component of spin1
    m_S1z               = 0;                                    // z component of spin1

    m_S2x               = 0;                                    // x component of spin2
    m_S2y               = 0;                                    // y component of spin2
    m_S2z               = 0;                                    // z component of spin2

    m_nu                = 0;                                    // Initialise orbital velocity
    m_iota              = 0;

    // Initialise other parameters to 0
    m_MSN           = 0;
    m_MSNPrime      = 0;
    m_MC            = 0;
    m_MCPrime       = 0;

    m_Theta         = 0;
    m_Phi           = 0;
    m_Psi           = 0;
    m_vRel          = 0;
    m_uK            = 0;
    m_Radius        = 0;
    //JIM BARRETT - 05/12/2016 - removing m_aPrime as it isn't used at all
    //m_aPrime        = 0;
    m_EPrime        = 0;
    m_ePrime        = 0;
    m_cosiPrime     = 0;
    m_iPrime        = 0;
    m_tc            = 0;
    m_beta          = 0;                // Angle between r and v, related to eccentricity (= pi/2 for circular e = 0)
    m_Survived      = true;
    m_Merged        = false;
    m_errorFlag     = false;
    m_mergesInHubbleTimeFlag = false;



    m_CHEvolutionFlag = false;
	//kIT x-RAY-BINARY parameters
    // m_HMXB          = false;
    // m_HMXBtime      = 0;//kit
    // m_XRayPower     = 0;//
    // m_XRayAveMass1  = 0;
    // m_XRayAveMass2  = 0;
    // m_XRayAveSep    = 0;
    // m_XRayBinaryType= 0;

    // Other useful things

    m_Mass1_Keep        = m_Mass1;                              // For now I want a copy of the original drawn mass (isn't that what m1, m1prime are for?)
    m_Mass2_Keep        = m_Mass2;                              // copy original drawn mass

    eventCounter = 0;
    same_RLOF_loop = false;
    current_stellar_pair =0;
    mt_primary_counter = 0;
    mt_secondary_counter =0;
    mt_primary_ep1 = 0;
    mt_primary_ep1_K1 = 0;
    mt_primary_ep1_K2 = 0;
    mt_primary_ep2 = 0;
    mt_primary_ep2_K1 = 0;
    mt_primary_ep2_K2 = 0;
    mt_primary_ep3 = 0;
    mt_primary_ep3_K1 = 0;
    mt_primary_ep3_K2 = 0;
    mt_secondary_ep1 =0;
    mt_secondary_ep1_K1 =0;
    mt_secondary_ep1_K2 =0;
    mt_secondary_ep2 =0;
    mt_secondary_ep2_K1 =0;
    mt_secondary_ep2_K2 =0;
    mt_secondary_ep3 =0;
    mt_secondary_ep3_K1 =0;
    mt_secondary_ep3_K2 =0;
    SN_primary_type_1 =0;
    SN_primary_type_2 =0;
    SN_primary_type_3 =0;
    SN_secondary_type_1 =0;
    SN_secondary_type_2 =0;
    SN_secondary_type_3 =0;
    CEE =0;
    CEE_instigator =0;
    CEE_failed =0;
    CEE_failed_instigator =0;
    CEE_wet =0;
    CEE_wet_instigator =0;
    stellar_type_K1 = 0;
    stellar_type_K2 = 0;
	binary_disbound = false;

    m_systemicVelocity = NEVER_SET;                  // Post supernova systemic velocity

	m_synchronizationTimescale = NEVER_SET;
	m_circularizationTimescale = NEVER_SET;

    weight = 1;
    samplingphase = -1;

}

// Specific constructor - if you want to produce and study a binary with these properties -- should now be superseded by options.individualSystem
BinaryStar::BinaryStar(double Mass1_ZAMS, double Mass2_ZAMS, double Metallicity1, double Metallicity2, double SemiMajorAxis, double Eccentricity, double RotationalVelocity1, double RotationalVelocity2, double Inclination, const gsl_rng *r, programOptions const &options, AISvariables const &aisvariables){
    // DEBUG: Used to check we are using correct constructor
    //std::cout << "Calling specific BinaryStar constructor 3)" << std::endl;

    // Initialise runtime
    m_runtime = 0.0;
    m_clockStart = clock();

    // Binary parameters
    // CPLB: Specify units!
    m_ID            = s_nIDGenerator++;                          // Generate unique ID
    m_Mass1         = Mass1_ZAMS;                                // Default value - can't initialise to 0 since q = m2/m1 (in Msol)
    m_Mass2         = Mass2_ZAMS;                                // Default value - can't initialise to 0 since q = m2/m1 (in Msol)
    m_MassRatio     = m_Mass1 / m_Mass2;                         // Calculated from above
    m_SemiMajorAxis = SemiMajorAxis;                             // Default value (in AU)
    m_SemiMajorAxisPrime= m_SemiMajorAxis;
    m_SemiMajorAxisPrev = m_SemiMajorAxis;
    m_RadiusRL1         = rocheLobeRadius(m_Mass1, m_Mass2);
    m_RadiusRLprev1     = m_RadiusRL1;
    //m_RadiusRLprev1     = 0;
    m_RadiusRL2         = rocheLobeRadius(m_Mass2, m_Mass1);
    m_RadiusRLprev2     = m_RadiusRL2;
    //    m_RadiusRLprev2     = 0;
    m_Eccentricity  	= Eccentricity;                              // Default value
    m_EccentricityPrev 	= Eccentricity;
    m_EccentricityPrime = Eccentricity;
    m_Inclinaton    = Inclination;                               // Initialise to 0
    m_delayTime     = 0;                                         // Initialise to 0 (UNITS?)
    m_Metallicity1  = Metallicity1;                              // Default solar metallicity (Zsol = 0.02)
    m_Metallicity2  = Metallicity2;                              // Default solar metallicity (Zsol = 0.02)


    // Initialise variables to hold parameters prior to 2nd supernova
    m_SemiMajorAxisPre2ndSN = 0.0;
    m_EccentricityPre2ndSN = 0.0;

    // Initial Binary Parameters -- kept constant as a record of the initial parameters of the binary
    m_SemiMajorAxisInitial = m_SemiMajorAxis;
    m_EccentricityInitial  = m_Eccentricity;

    m_orbitalVelocity	= sqrt(G1*(m_Mass1+m_Mass2)/pow(m_SemiMajorAxis,3.0));		// Initial orbital velocity according to previous parameterin yr-1 units, where [G] = 4*pi^2 AU^3 yr-2 Msol-1
    m_orbitalVelocityPrime = m_orbitalVelocity;
    m_orbitalVelocityPrev = m_orbitalVelocity;

    m_time = 0.0;
	m_dt   = 0.0;
    m_timePrev = 0.0;
    m_flagFirstTimeStep = true;

    // Differential quantities
    m_aTidesDiff = 0.0;
    m_omegaTidesDiff = 0.0;
    m_omegaTides = 0.0;

    m_aMagneticBrakingDiff = 0.0;
    m_angularMomentumMagneticBrakingDiff = 0.0;

    m_aMassLossDiff = 0.0;
    m_omegaMassLossDiff = 0.0;

    m_aMassTransferDiff = 0.0;
    m_omegaMassTransferDiff = 0.0;

    massTransferTracker = 0.0;
    rocheLobeTracker1 = 0.0;
    rocheLobeTracker2 = 0.0;
	m_massTransferTrackerHistory = NO_MASS_TRANSFER;

    m_massTransferFlag = false;
    m_massTransferZeroFlag=false;
    m_fastPhaseCaseA=false;
	m_aRocheLobeOverFlowErrorPrinting = false;
	m_aFastPhaseRocheLobeOverFlowErrorPrinting = false;
    m_initiateMassTransfer = false;
    m_insideMassTransfer = false;
    // m_massTransferFlag = false;
    m_dynamicalFlag = false;
    m_massTransferType = 0;
    m_massTransferRateThermal = 0;
    m_massTransferThermalFlag = false;
	m_massTransferFractionAccreted = options.massTransferFractionAccreted;
	m_massTransferJloss = options.massTransferJloss;
	m_stableRLOFafterCEE = false;
	
	// Mass transfer error flags
	calculateZetaThermalErrorMassFlag = false;
	calculateZetaThermalErrorRadiusFlag = false;
	calculateZetaThermalErrorAgeFlag = false;
	calculateZetaNuclearErrorMassFlag = false;
	calculateZetaNuclearErrorRadiusFlag = false;
	calculateZetaNuclearErrorAgeFlag = false;
	adaptiveRocheLobeOverFlowErrorMassFlag = false;
	adaptiveRocheLobeOverFlowErrorRadiusFlag = false;
	adaptiveRocheLobeOverFlowErrorAgeFlag = false;
	massTransferFastPhaseCaseAErrorMassFlag = false;
	massTransferFastPhaseCaseAErrorRadiusFlag = false;
	massTransferFastPhaseCaseAErrorAgeFlag = false;
	
    // X-ray binary
    m_is_xray_binary = false;

    // Common Envelope
    m_commonEnvelopeFlag = false;
    m_doubleCoreCommonEnvelopeFlag = false;
    m_stellarMerger = false;
	m_commonEnvelopeOccuredAtLeastOnce = false;
	m_RLOFPrimaryAfterCEE = false;
	m_RLOFSecondaryAfterCEE = false;
	m_optimisticCommonEnvelopeFlag = false;
	
	m_EccentricityPreCEE = NEVER_SET;
	m_EccentricityPostCEE = NEVER_SET;
	m_SemiMajorAxisPreCEE = NEVER_SET;
	m_SemiMajorAxisPostCEE = NEVER_SET;
	m_rocheLobe1to2PreCEE = NEVER_SET;
	m_rocheLobe1to2PostCEE = NEVER_SET;
	m_rocheLobe2to1PreCEE = NEVER_SET;
	m_rocheLobe2to1PostCEE = NEVER_SET;
    m_counterCEE = 0;
    m_immediateRLOFAfterCEE = false;
    m_simultaneousRLOFleadingToCEEFlag = false;
    m_mainSequenceAccretorDuringCEEFlag = false;
    m_EccentricityRLOF = NEVER_SET;

    //Addition by Coen 18-10-2017 for zeta study
    m_zetaRLOFanalytic= NEVER_SET;
    m_zetaRLOFnumerical= NEVER_SET;
	//Alejandro - 07/08/2018 for CE study
	m_zetaStarCompare	= NEVER_SET;


    // Spins
//    m_RotationalVelocity1 = RotationalVelocity1;                 // Rotational velocity of primary star
//    m_RotationalVelocity2 = RotationalVelocity2;                 // Rotational velocity of secondary star


    // Conserved quantities
    m_TotalEnergy = calculateTotalEnergy(m_SemiMajorAxis,m_Mass1,m_Mass2,star1.m_RZAMS,star2.m_RZAMS,star1.m_omega,star2.m_omega,m_orbitalVelocity,k_definition(star1.m_RZAMS, star1.m_Radius, star1.m_Mass, star1.m_stellarType),k_definition(star2.m_RZAMS, star2.m_Radius, star2.m_Mass, star2.m_stellarType));    // Default total energy
    m_TotalAngularMomentum = calculateAngularMomentum(m_SemiMajorAxis,m_Eccentricity,m_Mass1,m_Mass2,star1.m_RZAMS,star2.m_RZAMS,star1.m_omega,star2.m_omega,m_orbitalVelocity,k_definition(star1.m_RZAMS, star1.m_Radius, star1.m_Mass, star1.m_stellarType),k_definition(star2.m_RZAMS, star2.m_Radius, star2.m_Mass, star2.m_stellarType));   // Default total angular momentum

    // Energy quantities (dynamical)
    m_TotalEnergyPrime = m_TotalEnergy;
    m_TotalAngularMomentumPrime = m_TotalAngularMomentum;
    m_TotalAngularMomentumPrev = m_TotalAngularMomentum;


    m_Chi1          = 0;                                        // Default dimensionless magnitude for spin 1
    m_Chi2          = 0;                                        // Default dimensionless magnitude for spin 2 (should only get magnitude after SN!)

    // Need to make it so binary star contains two instances of star, star1, star2, to hold masses, spins, radii and luminosities etc.
    // star 1 initially more massive
    star1 = Star(m_Mass1, m_Metallicity1, options, r);
    star2 = Star(m_Mass2, m_Metallicity2, options, r);

    star1.m_omegaZAMS = RotationalVelocity1;
    star2.m_omegaZAMS = RotationalVelocity2;
    star1.m_omegaPrev = RotationalVelocity1;
    star2.m_omegaPrev = RotationalVelocity2;
    star1.m_omega = RotationalVelocity1;
    star2.m_omega = RotationalVelocity2;

    // Misalignments
    m_theta1_i      = 0;                                // Initial misalignment of star 1
    m_theta2_i      = 0;                                // Initial misalignment of star 2
    m_theta1        = 0;                                // By default, aligned
    m_theta2        = 0;                                // By default, aligned
    m_theta12       = 0;                                // By default, aligned
    m_deltaPhi      = 0;                                // By default, aligned
    m_deltaPhi_i    = 0;                                // By default, aligned

    m_S1x           = 0;                                // x component of spin1
    m_S1y           = 0;                                // y component of spin1
    m_S1z           = 0;                                // z component of spin1

    m_S2x           = 0;                                // x component of spin2
    m_S2y           = 0;                                // y component of spin2
    m_S2z           = 0;                                // z component of spin2

    m_nu            = 0;                                // Initialise orbital velocity
    m_iota          = 0;

    // Initialise other parameters to 0
    m_MSN           = 0;
    m_MSNPrime      = 0;
    m_MC            = 0;
    m_MCPrime       = 0;

    m_Theta         = 0;
    m_Phi           = 0;
    m_Psi           = 0;
    m_vRel          = 0;
    m_uK            = 0;
    m_Radius        = 0;
    //JIM BARRETT - 05/12/2016 - removing m_aPrime as it isn't used at all
    //m_aPrime        = 0;
    m_EPrime        = 0;
    m_ePrime        = 0;
    m_cosiPrime     = 0;
    m_iPrime        = 0;
    m_tc            = 0;
    m_beta          = 0;                // Angle between r and v, related to eccentricity (= pi/2 for circular e = 0)
    m_Survived      = true;
    m_Merged        = false;
    m_errorFlag     = false;
    m_mergesInHubbleTimeFlag = false;



    m_CHEvolutionFlag = false;
	//kIT x-RAY-BINARY parameters
    // m_HMXB          = false;
    // m_HMXBtime      = 0;//kit
    // m_XRayPower     = 0;//
    // m_XRayAveMass1  = 0;
    // m_XRayAveMass2  = 0;
    // m_XRayAveSep    = 0;
    // m_XRayBinaryType= 0;

    // redundant variables
    m_Mass1_Keep    = 0;
    m_Mass2_Keep    = 0;

        //Coen -- 10/11/2016 -- Formation Channel Printing parameters
    eventCounter = 0;
    same_RLOF_loop = false;
    current_stellar_pair =0;
    mt_primary_counter = 0;
    mt_secondary_counter =0;
    mt_primary_ep1 = 0;
    mt_primary_ep1_K1 = 0;
    mt_primary_ep1_K2 = 0;
    mt_primary_ep2 = 0;
    mt_primary_ep2_K1 = 0;
    mt_primary_ep2_K2 = 0;
    mt_primary_ep3 = 0;
    mt_primary_ep3_K1 = 0;
    mt_primary_ep3_K2 = 0;
    mt_secondary_ep1 =0;
    mt_secondary_ep1_K1 =0;
    mt_secondary_ep1_K2 =0;
    mt_secondary_ep2 =0;
    mt_secondary_ep2_K1 =0;
    mt_secondary_ep2_K2 =0;
    mt_secondary_ep3 =0;
    mt_secondary_ep3_K1 =0;
    mt_secondary_ep3_K2 =0;
    SN_primary_type_1 =0;
    SN_primary_type_2 =0;
    SN_primary_type_3 =0;
    SN_secondary_type_1 =0;
    SN_secondary_type_2 =0;
    SN_secondary_type_3 =0;
    CEE =0;
    CEE_instigator =0;
    CEE_failed =0;
    CEE_failed_instigator =0;
    CEE_wet =0;
    CEE_wet_instigator = 0;
    stellar_type_K1 = 0;
    stellar_type_K2 = 0;
	binary_disbound = false;

    m_systemicVelocity = NEVER_SET;                  // Post supernova systemic velocity

	m_synchronizationTimescale = NEVER_SET;
	m_circularizationTimescale = NEVER_SET;

    weight = 1;
    samplingphase = -1;
}

// Default Destructor
BinaryStar::~BinaryStar(void){

    // Calculate runtime for this binary in s
    m_runtime = ( clock() - m_clockStart ) / (double) CLOCKS_PER_SEC;

    if(s_firstBinary){
        s_firstBinary = false;
        
        s_runtimeOutput << "#" << TAB << "seconds" << std::endl;
        s_runtimeOutput << "int" << TAB << "float" << std::endl;
        s_runtimeOutput << "SEED" << TAB << "runtime" << std::endl;
    }

    s_runtimeOutput << m_randomSeed << TAB << m_runtime << std::endl;
    //std::cout << m_randomSeed << TAB << m_runtime << std::endl;

    //std::cout << "Binary is being deleted" << std::endl;
}

// Getter Functions

double BinaryStar::getMass1(bool si){return si? star1.m_Mass * Msol : star1.m_Mass;}
double BinaryStar::getMass2(bool si){return si? star2.m_Mass * Msol : star2.m_Mass;}
double BinaryStar::getMassRatio(){return star2.m_Mass / star1.m_Mass;}
double BinaryStar::getTotalMass(bool si){return si? (star1.m_Mass + star2.m_Mass) * Msol : (star1.m_Mass + star2.m_Mass);}
double BinaryStar::getReducedMass(bool si){return si? (star2.m_Mass * star1.m_Mass * Msol)/(star1.m_Mass + star2.m_Mass) : (star1.m_Mass * star2.m_Mass)/(star1.m_Mass + star2.m_Mass);}

double BinaryStar::getMZAMS1(){return star1.m_MZAMS;}
double BinaryStar::getMZAMS2(){return star2.m_MZAMS;}
double BinaryStar::getRZAMS1(){return star1.m_RZAMS;}
double BinaryStar::getRZAMS2(){return star2.m_RZAMS;}
double BinaryStar::getRadius1(){return star1.m_Radius;}
double BinaryStar::getRadius2(){return star2.m_Radius;}
double BinaryStar::getOmega1(){return star1.m_omega;}
double BinaryStar::getOmega2(){return star2.m_omega;}
double BinaryStar::getOrbitalVelocity(){return m_orbitalVelocity;}
double BinaryStar::getTotalEnergy(){return m_TotalEnergy;}
double BinaryStar::getTotalAngularMomentum(){return m_TotalAngularMomentum;}


double BinaryStar::getSemiMajorAxis(bool si){return si? m_SemiMajorAxis * AU : m_SemiMajorAxis;}
double BinaryStar::getEccentricity(){return m_Eccentricity;}
double BinaryStar::getMetallicity1(){return m_Metallicity1;}
double BinaryStar::getMetallicity2(){return m_Metallicity2;}
double BinaryStar::getInclination(){return m_Inclinaton;}
double BinaryStar::getRadius(bool si){return si? m_Radius*AU : m_Radius;}
double BinaryStar::getOrbitalEnergy(){return orbitalEnergy(getReducedMass(), getTotalMass(), getSemiMajorAxis());}
double BinaryStar::getChi1(){return m_Chi1;}
double BinaryStar::getChi2(){return m_Chi2;}
double BinaryStar::getTheta1(){return m_theta1;}
double BinaryStar::getTheta2(){return m_theta2;}
double BinaryStar::getTheta12(){return m_theta12;}
double BinaryStar::getDeltaPhi(){return m_deltaPhi;}



// Get post-SN parameters
//JIM BARRETT - 05/12/2016 - removing m_aPrime as it isn't used at all. None of these getters are used anymore
// double BinaryStar::getSemiMajorAxisPrime(bool si){return si? m_aPrime * AU : m_aPrime;}
// double BinaryStar::getEccentricityPrime(){return m_ePrime;}
// double BinaryStar::getM1Prime(bool si){return si? m_MSNPrime * Msol : m_MSNPrime;}
// double BinaryStar::getM2Prime(bool si){return si? m_MCPrime * Msol : m_MC;}


// Get period of binary (can do a version in different units (days?). -- use the INSI bool)
double BinaryStar::getPeriod(){
    return sqrt((4*pi*pi)/(G*(m_Mass1 + m_Mass2)))*pow(m_SemiMajorAxis, (3.0/2.0)); // return period in SI units (s)
}

void BinaryStar::supernova(programOptions const &options, const gsl_rng *r, Star &starSupernova, const Star &starCompanion, unsigned long m_randomSeed, double which_star){
    /*
     One of the stars has gone supernova! 

     Assign a random supernova kick according to the user specified options and then update the orbit

     Parameters
     -----------
     options : programOptions
        User specified program options
     r : gsl_rng
        GSL random number generator
     starSupernova : Star
        Star that went supernova
     starCompanion : Star
        The other star
     m_randomSeed : unsigned long
        random seed
     which_star : double
        which star?

     Returns
     --------
     */
    
    bool debugging = false;

    if(debugging){
        std::cout << "Kick before: " << starSupernova.m_kickVelocity << std::endl;
    }

    if((starSupernova.flagSN == true) or (starSupernova.flagECSN == true) or (starSupernova.flagUSSN == true)){

		// Masses should already be correct, mass before SN given by star.m_MassPrev

		// Declare variables and initialise to zero
		//double vK = 0, theta = 0, phi = 0, psi = 0, uK = 0, vRel = 0, radius = 0, beta = 0, vsys=0;
		double vK = 0, uK = 0, vRel = 0, radius = 0, beta = 0, psi=0, vsys=0;

		double theta = starSupernova.m_supernovaTheta;
		double phi = starSupernova.m_supernovaPhi;

        double M1Prime = starSupernova.m_Mass;
        double deltaM1 = starSupernova.m_MassPrev - starSupernova.m_Mass;

		// Generate true anomaly - (for e=0, should be a flat distribution) - updates Eccentric anomaly and True anomaly automatically
		// ALEJANDRO - 09/05/2018 - If statement to avoid solving Kepler's equation for an unbound orbit; it may be of interest to have SN of unbound stars in the supernovae.txt file.
		if(m_SemiMajorAxisPrime > 0.0){
    		solveKeplersEquation(starSupernova.m_MeanAnomaly, m_Eccentricity, starSupernova.m_EccentricAnomaly, starSupernova.m_TrueAnomaly);
    		psi = starSupernova.m_TrueAnomaly;
		}
		else{
    		// ALEJANDRO - 09/05/2018 - Following 3 lines copied from else statement in the end.
    		eventCounter +=1;
    		binary_disbound = true;
            m_Survived = false;
		}
    

		if(debugging){
			std::cout << "Mean Eccentric True Anomaly " << std::endl;
			std::cout << starSupernova.m_MeanAnomaly << " " << starSupernova.m_EccentricAnomaly << " " << starSupernova.m_TrueAnomaly << std::endl;
		}

		if(starSupernova.flagECSN){
			// ALEJANDRO - 04/05/2017 - Allow for ECSN to have kicks different than zero. Still, should be low kicks. Default set to zero.
			double	sigmaECSN = options.kickVelocityDistributionSigmaForECSN;
			
			if(debugging){
				std::cout << "Maxwellian sigma for an ECSN:\t" << sigmaECSN << std::endl;
				std::cout << "This is an electron capture supernova" << std::endl;
			}
		
			vK = kickVelocityDistribution(options, r, sigmaECSN, starSupernova.m_COCoreMassAtCompactObjectFormation, m_randomSeed, starSupernova.m_supernovaKickVelocityMagnitudeRandomNumber, deltaM1, M1Prime);

			// Draw a theta and phi according to the user specified distribution
			//kickDirectionDistribution(options, r, theta, phi, m_kickDirectionPower, m_randomSeed); // TODO: Set as initial params
        
			starSupernova.m_fallback = 0.0;
			starSupernova.m_drawnKickVelocity = vK;
        
		}
		else if(starSupernova.flagUSSN){
			// ALEJANDRO - 25/08/2017 - Allow for USSN to have a separate kick.
			double	sigmaUSSN = options.kickVelocityDistributionSigmaForUSSN;
				
			if(debugging){
				std::cout << "Maxwellian sigma for an USSN:\t" << sigmaUSSN << std::endl;
				std::cout << "This is an ultra-stripped supernova" << std::endl;
			}
		
			vK = kickVelocityDistribution(options, r, sigmaUSSN, starSupernova.m_COCoreMassAtCompactObjectFormation, m_randomSeed, starSupernova.m_supernovaKickVelocityMagnitudeRandomNumber, deltaM1, M1Prime);

			// Draw a theta and phi according to the user specified distribution
			//kickDirectionDistribution(options, r, theta, phi, m_kickDirectionPower, m_randomSeed); // TODO: Set as initial params
        
			starSupernova.m_fallback = 0.0;
			starSupernova.m_drawnKickVelocity = vK;
        
		}
		else if(starSupernova.flagSN){
	
			// Draw a random kick velocity from the user selected distribution, with sigma based on whether compact object is a NS or BH
			if(starSupernova.m_stellarType == NEUTRON_STAR){
				vK = kickVelocityDistribution(options, r, options.kickVelocityDistributionSigmaCCSN_NS, starSupernova.m_COCoreMassAtCompactObjectFormation, m_randomSeed, starSupernova.m_supernovaKickVelocityMagnitudeRandomNumber, deltaM1, M1Prime);
			}
			else if(starSupernova.m_stellarType == BLACK_HOLE){
				vK = kickVelocityDistribution(options, r, options.kickVelocityDistributionSigmaCCSN_BH, starSupernova.m_COCoreMassAtCompactObjectFormation, m_randomSeed, starSupernova.m_supernovaKickVelocityMagnitudeRandomNumber, deltaM1, M1Prime);
			}

			// Save kick velocity drawn from Maxwellian
			starSupernova.m_drawnKickVelocity = vK;
            
			// Draw a theta and phi according to the user specified distribution
			//kickDirectionDistribution(options, r, theta, phi, m_kickDirectionPower, m_randomSeed); // TODO: Set as initial params

			if(debugging){
				std::cout << "This is a CCSN supernova" << std::endl;
				std::cout << "stellar type " << starSupernova.m_stellarType << std::endl;
				std::cout << "Drawn kick = " << starSupernova.m_drawnKickVelocity << std::endl;
			}

			// Re-weight kicks by mass of remnant accoridng to user specification
			blackHoleKicks(options, vK, starSupernova.m_fallback, starSupernova.m_Mass, starSupernova.m_stellarType);

			if(debugging){
				std::cout << "Kick after applying blackHoleKicks with fallback of " << starSupernova.m_fallback << " = " << vK << " km s^-1" << std::endl;
			}
		}
		else{
        std::cerr << m_randomSeed <<  "\tWhat kind of supernova gets here??" << std::endl;
		}
	
		if(debugging){
			std::cout << "Mass of supernova star, remnant = " << starSupernova.m_MassPrev << ", " << starSupernova.m_Mass << " Msol" << std::endl;
            std::cout << "Mass of companion star = " << starCompanion.m_Mass << " Msol" << std::endl;
            std::cout << "Type of supernova remnant = " << starSupernova.m_stellarType << std::endl;
            std::cout << "Type of companion = " << starCompanion.m_stellarType << std::endl;
            std::cout << "Kick direction theta, phi = " << theta << " " << phi << std::endl;
            std::cout << "Kick before applying blackHoleKicks = "<< vK << " km s^-1" << std::endl;
		}

		// Save kick velocity after fallback
		starSupernova.m_kickVelocity = vK;

		// Calculate radius of orbit at current time in AU as a function of the true anomaly psi
		radius = radiusFunctionOfPsi(m_SemiMajorAxisPrime, m_Eccentricity, psi);    

		// Assign radius to this binaries orbital radius value
		m_Radius = radius;

		if(debugging){
			std::cout << "Initial semi-major axis, eccentricity = " << m_SemiMajorAxisPrime << " " << m_Eccentricity << std::endl;
			std::cout << "radius = " << radius << std::endl;
		}

		// Calculate the angle between the position and velocity vectors
		double Mtot = starSupernova.m_MassPrev + starCompanion.m_MassPrev;
		double reducedMass = (starSupernova.m_MassPrev * starCompanion.m_MassPrev) / Mtot;

		// calculate mass combinations after SN
		double MtotPrime = starSupernova.m_Mass + starCompanion.m_Mass;
		double reducedMassPrime = (starSupernova.m_Mass * starCompanion.m_Mass) / MtotPrime;

		if(debugging){
			std::cout << "Mtot, MtotPrime = " << Mtot << " " << MtotPrime << std::endl;
			std::cout << "m_semiMajorAxisPrime, m_Eccentricity, m_Radius = " << m_SemiMajorAxisPrime << " " << m_Eccentricity << " " << m_Radius << std::endl;
		}
    
		// JIM BARRETT - 28/11/2016 - Passing through the random seed for error tracking purposes
		// beta = angleBetweenPositionAndVelocity(m_SemiMajorAxisPrime, m_Eccentricity, Mtot, m_Radius);
		// beta = angleBetweenPositionAndVelocity(m_SemiMajorAxisPrime, m_Eccentricity, Mtot, m_Radius, m_randomSeed);
		beta = angleBetweenPositionAndVelocity(m_SemiMajorAxisPrime, m_Eccentricity, MtotPrime, m_Radius, m_randomSeed);

		if(debugging){std::cout << "beta = " << beta << std::endl;}
    
		///////////////////////////////////////////////////////////////////////////////////////
		//              AT THE MOMENT, QUANTITIES BEYOND HERE ARE IN SI (NOT IDEAL)
		///////////////////////////////////////////////////////////////////////////////////////

		// Convert vK to m s^-1 // Would be nice to draw this in nicer units to avoid this secion
		vK  *= KM;

		// Equation for vRelSquared at this point in the orbit prior to the SN. (constants to convert to ms^-1)
		// Returns vRel
		// vRel = orbitalVelocity(Mtot * Msol, getRadius(INSI), getSemiMajorAxis(INSI));
		vRel = orbitalVelocity(Mtot * Msol, m_Radius * AU, m_SemiMajorAxisPrime * AU);

		m_orbitalVelocityPre2ndSN = vRel;

		// Since the kick velocity always occurs in equations as vk/vrel, we need to know vrel. We can also create a variable which is uk = vk/vrel (dimensionless - should be around 10) #-- not 0.1?

		uK = vK/vRel;

		if(debugging){std::cout << "uK before testing if fixed = " << uK << std::endl;}

		///////////////////////////////////////////////////////////////////////////////////////
		//                      SHOULD BE BACK TO NICE UNITS NOW
		///////////////////////////////////////////////////////////////////////////////////////

		// Check if user wants to fix uK to a certain value and if so change uk to that value.
		fixUK(options, uK);

		if(debugging){std::cout << "uK after testing if fixed = " << uK << std::endl;}

		// Check pre-SN orbital energy
		double E = 0;

	//    E = orbitalEnergy(reducedMass, Mtot, getSemiMajorAxis());    // Should be -ve by construction
		E = orbitalEnergy(reducedMass, Mtot, m_SemiMajorAxisPrime);    // Should be -ve by construction

		// Seemed to be getting into this loop occasionally with E > 0 but E ~ 0 (1e-37 for example) -- what's going on?
		if(E > 0){
			// Check that binaries are bound - E should be -ve by construction
			std::cerr << m_randomSeed <<  "\tERROR: E > 0. That's not good. E :\t" << E << std::endl;
		}

		// Now calculate post-SN orbital properties

		double ePrime = 0;
		double aPrime = 0;
		double EPrime = 0;
		double cosiPrime = 0;
		double iPrime    = 0;

		// Record the semi major axis and eccentricity just before each supernova
		m_SemiMajorAxisPre2ndSN = m_SemiMajorAxisPrime;
		m_EccentricityPre2ndSN = m_Eccentricity;

		//    aPrime = aPostSupernova(uK, Mtot, MtotPrime, radius, getSemiMajorAxis(), theta, phi, m_randomSeed);
		//    ePrime = ePostSupernova(uK, theta, phi, beta, Mtot, MtotPrime, radius, m_SemiMajorAxis, m_randomSeed);
		aPrime = aPostSupernova(uK, Mtot, MtotPrime, radius, m_SemiMajorAxisPrime, theta, phi, m_randomSeed);
		ePrime = ePostSupernova(uK, theta, phi, beta, Mtot, MtotPrime, radius, m_SemiMajorAxisPrime, m_randomSeed);

		if(debugging){

			double ePrimeTest = ePostSupernovaInitiallyCircular(Mtot, MtotPrime, theta, phi, uK);
			double aPrimeTest = aPostSupernovaInitiallyCircular(m_SemiMajorAxisPrime, Mtot, MtotPrime, theta, phi, uK);

			std::cout << "aPrime, ePrime = " << aPrime << " " << ePrime << std::endl;
			std::cout << "aPrime, ePrime assume circular  = " << aPrimeTest << " " << ePrimeTest << std::endl;
			std::cout << "test1 a: " << getSemiMajorAxis() << std::endl;
			std::cout << "test2 a: " << m_SemiMajorAxis << std::endl;
			std::cout << "test3 a: " << m_SemiMajorAxisPrev << std::endl;
			std::cout << "test4 a: " << m_SemiMajorAxisPrime << std::endl;
		}

		// Check post-SN orbital energy, check if still bound

		EPrime = orbitalEnergy(reducedMassPrime, MtotPrime, aPrime);

		// Caclulate dimensionless post-SN orbital energy

		double epsilon = 0;

		epsilon = -EPrime / E;

		// Do some debugging
		if(debugging){
			std::cout << "vk : " << vK << std::endl;
			std::cout << "uk : " << uK << std::endl;
			std::cout << "E : " << E << std::endl;
			std::cout << "EPrime : " << EPrime << std::endl;
			std::cout << "epsilon : " << epsilon << std::endl;
			std::cout << "a : " << m_SemiMajorAxisPrime << std::endl;
			std::cout << "e : " << m_Eccentricity << std::endl;
		}

		// Check if still bound
		if(epsilon < 0){

			// Confirm that it survived:
			m_Survived = true;

			// Calculate post-SN orbital inclination using the equation for arbitrary eccentricity orbits
			cosiPrime   = cosFinalPlaneTilt(uK, beta, theta, phi);
			iPrime      = acos(cosiPrime);

			// void function to assign the spins based on the assumption we are using
			// TODO: I think this function is currently broken for two supernovae -- check!!
			assignMisalignments(m_theta1_i, m_theta2_i, iPrime, options, r);

			// Variables to evolve
			m_theta1 = m_theta1_i;
			m_theta2 = m_theta2_i;

			// why semimajor axis and aPrime?

			// Save to output those binaries which survive
			//output << m_ID << TAB << starCompanion.m_MassPrev << TAB << starCompanion.m_Mass << TAB << starSupernova.m_MassPrev << TAB << starSupernova.m_Mass << TAB << m_SemiMajorAxis << TAB << aPrime << TAB << m_Eccentricity << TAB << ePrime << TAB << uK << TAB << theta << TAB << phi << TAB << m_TrueAnomaly << TAB << m_EccentricAnomaly << TAB << EPrime << TAB << cosiPrime << TAB << iPrime << NL;

			// Calculate post-SN systemic (center-of-mass) velocity
			M1Prime = starSupernova.m_Mass;
			deltaM1 = starSupernova.m_MassPrev - starSupernova.m_Mass;
			double M2 = starCompanion.m_Mass;

			vsys = postSNSystemicVelocity(uK, vRel, Mtot, MtotPrime, M1Prime, deltaM1, M2, theta, phi, beta); // in ms s^-1
			vsys /= KM;     // Convert to km s^-1
        
			if(debugging){std::cout << vRel << " " << vsys << std::endl;}

			}
    else{
		
        // Confirm that it did not survive
		eventCounter +=1;
		binary_disbound=true;
        m_Survived = false;
		
		if(starCompanion.m_stellarType <= NAKED_HELIUM_STAR_GIANT_BRANCH)
			star2.m_runawayFlag = true;
    }
	

    ///////////////////////////////////////////////////////////////////////////////////////
    //                      UPDATE BINARIES PARAMETERS AND SAVE TO OUTPUT
    ///////////////////////////////////////////////////////////////////////////////////////

    m_MSN       = starSupernova.m_MassPrev;                     // Exploding star pre-SN mass
    m_MSNPrime  = starSupernova.m_Mass;                         // Exploding star post-SN mass
    m_MC        = starCompanion.m_MassPrev;                     // Companion star pre-SN mass
    m_MCPrime   = starCompanion.m_Mass;                         // Companion star post-SN mass

    //JIM BARRETT - 05/12/2016 - removing m_aPrime as it isn't used at all
    //m_aPrime    = aPrime;
    m_ePrime    = ePrime;
    m_EPrime    = EPrime;

    m_cosiPrime = cosiPrime;
    m_iPrime    = iPrime;

    m_vRel      = vRel;
    m_uK        = uK;
	// Alejandro - 17/09/2017 - Added star variables to have them handy
    starSupernova.m_supernovaTheta     = theta;
    starSupernova.m_supernovaPhi       = phi;

    m_Radius    = radius;
    m_beta      = beta;

    m_SemiMajorAxisPrime = aPrime;
    m_Eccentricity = ePrime;
    m_EccentricityPrime = ePrime;

    starSupernova.m_preSNeOrbitalEnergy=E;
    starSupernova.m_postSNeOrbitalEnergy=EPrime;

    m_systemicVelocity = vsys;                  // Post supernova systemic velocity

    if(debugging){
        std::cout << "a' : " << m_SemiMajorAxisPrime << std::endl;
        std::cout << "e' : " << m_Eccentricity << std::endl;
    }
   
    // Should move this output to a printing.cpp function. Supernova printing
    std::ofstream supernovae((options.outputPath/"supernovae.txt").string(),std::ios_base::app);
    supernovae << 	m_randomSeed << TAB << 	
					starSupernova.m_drawnKickVelocity << TAB << 
					starSupernova.m_kickVelocity << TAB << 	
					starSupernova.m_fallback << TAB << 
					vRel << TAB << 	
					uK << TAB << 
					psi << TAB << 
					starSupernova.m_supernovaTheta << TAB << 
					starSupernova.m_supernovaPhi << TAB << 	
					starSupernova.flagECSN 	<< TAB << 
					starSupernova.flagSN 	<< TAB << 	
					starSupernova.flagUSSN 	<< TAB << 
                    starSupernova.flagExperiencedPISN    << TAB <<   
                    starSupernova.flagExperiencedPPISN  << TAB << 
					m_Survived	<< TAB << 	
					starSupernova.m_MZAMS  	<< TAB << 
					starCompanion.m_MZAMS	<< TAB <<	
					starSupernova.m_totalMassAtCompactObjectFormation << TAB << 
					starCompanion.m_Mass << TAB << 	
					starSupernova.m_COCoreMassAtCompactObjectFormation << TAB << 
					starSupernova.m_Mass         << TAB <<  
					starSupernova.experiencedRLOF << TAB <<
					starSupernova.m_stellarType  << TAB <<	
					which_star << TAB <<
					starSupernova.m_stellarTypePrev << TAB <<  
					starCompanion.m_stellarTypePrev << TAB <<
					starSupernova.m_coreMassAtCompactObjectFormation << TAB <<  
					starSupernova.m_Metallicity << TAB << 
					m_commonEnvelopeOccuredAtLeastOnce << TAB <<  
					m_stableRLOFafterCEE << TAB << 
					starSupernova.flagRLOFontoaNS << TAB <<  
					m_time << TAB <<
                    m_EccentricityPre2ndSN << TAB <<  
					m_Eccentricity <<  TAB << 
                    m_SemiMajorAxisPre2ndSN*AUToRsol  << TAB <<  
					m_SemiMajorAxisPrime*AUToRsol << TAB << 
                    m_systemicVelocity << TAB <<
					starSupernova.m_flagHrichSN << TAB <<
					starSupernova.m_flagHpoorSN << TAB << 
					star2.m_runawayFlag << std::endl;

    supernovae.close();
    // Update binary masses redundant since should use star.mass
    //m_Mass1     = MCPrime;                          // Primary mass
    //m_Mass2     = MSNPrime;                         // Secondary mass

    // Save everything regardless of surviving or not to examine full distributions.
    // ideally would want m1, m2 here to be the initial m1, m2 of those binaries which will go supernova (or alternatively, just add an m1init m2init columns?)
    // HERE -- m_Mass1 << TAB << m_Mass2

    //outputAll << m_ID << TAB << m_Mass1_Keep << TAB << m_Mass2_Keep << TAB << starCompanion.m_MassPrev << TAB << starCompanion.m_Mass << TAB << starSupernova.m_MassPrev << TAB << starSupernova.m_Mass << TAB << m_SemiMajorAxis << TAB << aPrime << TAB << m_Eccentricity << TAB << ePrime << vK << TAB << vRel << TAB << uK << TAB << theta << TAB << phi << TAB << m_TrueAnomaly << TAB << m_EccentricAnomaly << TAB << EPrime << TAB << cosiPrime << TAB << iPrime << NL;

    // why are m_SemiMajorAxis and aPrime different, shouldn't this update the orbital parameters?
        starSupernova.flagSN = false;
        starSupernova.flagECSN = false;
    }
    else{
        if(debugging){
            std::cout << "No Supernova. This function is called ofted to check for SNe." << std::endl;
        }
    }

    if(debugging){
        std::cout << "Kick after: " << starSupernova.m_kickVelocity << std::endl;
    }

}


void BinaryStar::coalesce(const programOptions &options, std::ofstream & output){
    /*
     Function which calculates time to coalescence, whether binary merges within hubble time, records parameters of those that do

     Parameters
     -----------
     options : programOptions
        User specified programOptions used for this run
     ouput : std::ofstream
        Output stream to use to record parameters of merging binaries

     Returns
     --------

     */
    // Calculate the time for the binary to coalesce due to emission of gravitational radiation.
    bool    debugging = false;
    double tc = 0;

    // Define DCO formation to be now
    m_SemiMajorAxisAtDCOFormation = m_SemiMajorAxisPrime;
    m_EccentricityAtDCOFormation = m_Eccentricity;

    tc = timeToCoalescenceUsingInterpolation(m_SemiMajorAxisPrime * AU, m_Eccentricity, star1.m_Mass * Msol, star2.m_Mass * Msol);

    // Record coalescence time in Myrs
    m_tc = (tc/year)*yearToMyear;
    
    // DEBUG
    if(debugging){
        std::cout << "a : " << m_SemiMajorAxisPrime << std::endl;
        std::cout << "e : " << m_Eccentricity << std::endl;
        std::cout << "a at DCO formation : " << m_SemiMajorAxisAtDCOFormation << std::endl;
        std::cout << "e at DCO formation : " << m_EccentricityAtDCOFormation << std::endl;
        std::cout << "tc : " << m_tc << std::endl;
        std::cout << "M1 : " << star1.m_Mass << std::endl;
        std::cout << "M2 : " << star2.m_Mass << std::endl;
        std::cout << "No error: tc:\t" << tc << std::endl;
    }
    // test if shorter than HubbleTime (will need to worry about time delays eventually and time when born)

    if(tc < HubbleTime){

        // It merged in a hubble time
        // Why do we have 2 flags that do the same thing?
        m_Merged = true;
        m_mergesInHubbleTimeFlag = true;

        if(!options.quiet){
            std::cout << m_randomSeed << "\tMerges in less than a Hubble time, tc: " << tc << std::endl;
        }

    }
    else{

        // It did not merge -- why do we have two flags that do the same thing?
        m_Merged = false;
        m_mergesInHubbleTimeFlag = false;

        if(!options.quiet){
            std::cout << m_randomSeed << "\tDoesn't merge in a Hubble time, tc (in s), tc/t_Hubble = " << tc << " " << tc / HubbleTime << std::endl;
        }

    }
    
    // Output parameters of merging binary. Should move to a function in printing.cpp
    output	<< m_ID 
			<< TAB << m_randomSeed 
            << TAB << weight 
            << TAB << samplingphase
			<< TAB << m_SemiMajorAxisInitial 
			<< TAB << m_EccentricityInitial 
			<< TAB << m_SemiMajorAxisPre2ndSN 
			<< TAB << m_EccentricityPre2ndSN 
			<< TAB << m_orbitalVelocityPre2ndSN 
			<< TAB << m_SemiMajorAxisAtDCOFormation 
			<< TAB << m_EccentricityAtDCOFormation 
			<< TAB << star1.m_Metallicity 
			<< TAB << star2.m_Metallicity 
			<< TAB << star1.m_MZAMS 
			<< TAB << star1.m_totalMassAtCompactObjectFormation 
			<< TAB << star1.m_HeCoreMassAtCompactObjectFormation 
			<< TAB << star1.m_COCoreMassAtCompactObjectFormation 
			<< TAB << star1.m_coreMassAtCompactObjectFormation 
			<< TAB << star1.m_HeCoreMassAtCommonEnvelope 
			<< TAB << star1.m_COCoreMassAtCommonEnvelope 
			<< TAB << star1.m_coreMassAtCommonEnvelope 
			<< TAB << star1.m_kickVelocity 
			<< TAB << star1.m_supernovaTheta
			<< TAB << star1.m_supernovaPhi
			<< TAB << star1.m_Mass 
			<< TAB << star1.m_stellarType 
			<< TAB << star2.m_MZAMS 
			<< TAB << star2.m_totalMassAtCompactObjectFormation 
			<< TAB << star2.m_HeCoreMassAtCompactObjectFormation 
			<< TAB << star2.m_COCoreMassAtCompactObjectFormation 
			<< TAB << star2.m_coreMassAtCompactObjectFormation 
			<< TAB << star2.m_HeCoreMassAtCommonEnvelope 
			<< TAB << star2.m_COCoreMassAtCommonEnvelope 
			<< TAB << star2.m_coreMassAtCommonEnvelope 
			<< TAB << star2.m_kickVelocity 
			<< TAB << star2.m_supernovaTheta
			<< TAB << star2.m_supernovaPhi			
			<< TAB << star2.m_Mass 
			<< TAB << star2.m_stellarType 
			<< TAB << m_tc 
			<< TAB << m_time 
			<< TAB << m_luminousBlueVariableFactor 
            << TAB << m_kickVelocityDistributionSigmaCCSN_NS 
            << TAB << m_kickVelocityDistributionSigmaCCSN_BH
            << TAB << m_commonEnvelopeAlpha 
			<< TAB << m_kickDirectionPower 
			<< TAB << m_wolfRayetFactor 
			<< TAB << m_RLOFSecondaryAfterCEE 
			<< TAB << star1.m_initialMassTransferCase 
			<< TAB << star2.m_initialMassTransferCase 
			<< TAB << star1.m_preSNeOrbitalEnergy 
			<< TAB << star1.m_postSNeOrbitalEnergy 
			<< TAB << star2.m_preSNeOrbitalEnergy 
			<< TAB << star2.m_postSNeOrbitalEnergy 
			<< TAB << m_commonEnvelopeOccuredAtLeastOnce 
			<< TAB << star1.m_lambdaAtCommonEnvelope 
			<< TAB << star2.m_lambdaAtCommonEnvelope 
			<< TAB << m_EccentricityPreCEE 
			<< TAB << m_EccentricityPostCEE 
			<< TAB << m_SemiMajorAxisPreCEE 
			<< TAB << m_SemiMajorAxisPostCEE 
			<< TAB << m_rocheLobe1to2PreCEE 
			<< TAB << m_rocheLobe1to2PostCEE 
			<< TAB << m_rocheLobe2to1PreCEE 
			<< TAB << m_rocheLobe2to1PostCEE 
			<< TAB << m_optimisticCommonEnvelopeFlag 
			<< TAB << m_mergesInHubbleTimeFlag 
			<< TAB << m_doubleCoreCommonEnvelopeFlag 
			<< TAB << star1.m_bindingEnergyAtCommonEnvelope 
			<< TAB << star2.m_bindingEnergyAtCommonEnvelope 
			<< TAB << star1.flagRecycledNS 
			<< TAB << star2.flagRecycledNS 
			<< TAB << star1.flagUSSN 
			<< TAB << star2.flagUSSN 
			<< TAB << star1.flagExperiencedECSN
            << TAB << star2.flagExperiencedECSN
            << TAB << star1.flagExperiencedPISN
            << TAB << star2.flagExperiencedPISN
            << TAB << star1.flagExperiencedPPISN 
			<< TAB << star2.flagExperiencedPPISN << std::endl;

}



// ODEINT definitions
// Define state_type
const int nVariables = 10;
typedef std::array<double,nVariables> state_type;                            // The "normal" way

int nCalc = 0;              // For interest (maybe useful later) count how many times the deriv function is called

struct spinPrecession{

    // Parameters
    double m_m1, m_m2, m_q;

    // Constructor
    spinPrecession(double m1, double m2){
        m_m1 = m1;
        m_m2 = m2;
        m_q  = m_m2/m_m1;           // needed?
    }

    void operator()(const state_type &x, state_type &dxdt, double t) const {
        // So there is an X state_type which has all of the variables in.
        // And then there is the dxdt state_type which contains the derivatives of these variables - one for each variable
        // You only need to assign values to each of the derivatives, the updating of the 'X' variables is handled by the stepper

        // State vectors:
        // X[0]  = S1x
        // X[1]  = S1y
        // X[2]  = S1z
        // X[3]  = S2x
        // X[4]  = S2y
        // X[5]  = S2z
        // X[6]  = Lhatx
        // X[7]  = Lhaty
        // X[8]  = Lhatz
        // X[9]  = nu

        // dxdt[0]  = dS1x / dt                          // dS1 / dt = Omega1 X S1
        // dxdt[1]  = dS1y / dt                          // dS1 / dt = Omega1 X S1
        // dxdt[2]  = dS1z / dt                          // dS1 / dt = Omega1 X S1
        // dxdt[3]  = dS2x / dt                          // dS2 / dt = Omega2 X S2
        // dxdt[4]  = dS2y / dt                          // dS2 / dt = Omega2 X S2
        // dxdt[5]  = dS2z / dt                          // dS2 / dt = Omega2 X S2
        // dxdt[6]  = dLhatx / dt
        // dxdt[7]  = dLhaty / dt
        // dxdt[8]  = dLhatz / dt
        // dxdt[9]  = dnu / dt

        dS1dt(dxdt[0], dxdt[1], dxdt[2], x[0], x[1], x[2], x[3], x[4], x[5], m_m1, m_m2, x[6], x[7], x[8], x[9]);
        dS2dt(dxdt[3], dxdt[4], dxdt[5], x[0], x[1], x[2], x[3], x[4], x[5], m_m1, m_m2, x[6], x[7], x[8], x[9]);
        dLhatdt(dxdt[6], dxdt[7], dxdt[8], m_m1, m_m2, x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], dxdt[0], dxdt[1], dxdt[2], dxdt[3], dxdt[4], dxdt[5]);
        dnudt(m_m1, m_m2, dxdt[9], x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9]);

        nCalc++;
    }
};
// JIM BARRETT - 05/12/2016 - This function isn't currently called anywhere, and has a bug with using the outdated
// getSemiMajorAxisPrime() getter method. Commenting out for now, this function will need thorough checking if we ever use it

// void BinaryStar::PNSpinEvolution(const gsl_rng *r){
//     // evolve spins using PN evolution.
//
//     double Mtot = getTotalMass() * Msol_nat;                            // Total mass in G = c = 1 units
//     double q    = getMassRatio();                                       // Mass ratio
//     double m1   = Mtot / (1.0 + q);                                     // Primary mass in G = c = 1 units
//     double m2   = Mtot / (1.0 + 1.0/q);                                 // Secondary mass in G = c = 1 units
//
//     double chi1 = getChi1();                                            // dimensionless spin 1
//     double chi2 = getChi2();                                            // dimensionless spin 2
//
// //    // DEBUGGING
// //    std::cout << "m1 : " << m1 << std::endl;
// //    std::cout << "m2 : " << m2 << std::endl;
// //    std::cout << "chi1 : " << chi1 << std::endl;
// //    std::cout << "chi2 : " << chi2 << std::endl;
//
//     // Do the masses match up with those output by the program? problem with getMassRatio()
//
//     // Spins
//     double theta1       = getTheta1();                                  // Angle between S1 and orbital angular momentum
//     double theta2       = getTheta2();                                  // Angle between S2 and orbital angular momentum
//     double deltaPhi     = gsl_rng_uniform(r) * 2.0 * M_PI;              // Generate deltaPhi uniform between 0, 2pi
//
//     // Initialise spin vectors:
//     double s1x = 0.0;                                           // Initialise spin 1
//     double s1y = 0.0;                                           // Initialise spin 1
//     double s1z = 0.0;                                           // Initialise spin 1
//
//     double s2x = 0.0;                                           // Initialise spin 2
//     double s2y = 0.0;                                           // Initialise spin 2
//     double s2z = 0.0;                                           // Initialise spin 2
//
//     s1x = chi1 * m1 * m1 * sin(theta1) * sin(deltaPhi);         // x component of spin 1
//     s1y = chi1 * m1 * m1 * sin(theta1) * cos(deltaPhi);         // y component of spin 1
//     s1z = chi1 * m1 * m1 * cos(theta1);                         // z component of spin 1
//
// //    // DEBUG check magnitude
// //    double norms1 = norm3(s1x, s1y, s1z);                       // magnitude of spin 1 (needed?)
// //    std::cout << "Initially:" << std::endl;
// //    std::cout << norms1 << std::endl;
//
//     s2x = chi2 * m2 * m2 * 0.0;                                 // x component of S2
//     s2y = chi2 * m2 * m2 * sin(theta2);                         // y component of S2
//     s2z = chi2 * m2 * m2 * cos(theta2);                         // z component of S2
//
// //    // DEBUG
// //    double norms2 = norm3(s2x, s2y, s2z);                       // magnitude of spin 2  (needed?)
// //    std::cout << norms2 << std::endl;
// //    //unit3(s2x, s2y, s2z, us2x, us2y, us2z); // just divides by norm3(s2x, s2y, s2z)
// //    // how do I check this is consistent?
//
//     // USE THESE INSTEAD?
//     // Need to fix the following 2 functions to deal with deltaPhi = 0 or deltaPhi = NAN? So just fix the calculate deltaPhi funciton? This should now be done.
//     //calculateS1(s1x, s1y, s1z, chi1, m1, theta1, deltaPhi);
//     //calculateS2(s2x, s2y, s2z, chi2, m2, theta2);
//
//     // Initialise orbital angular momentum vector
//     double Lhatx = 0.0;                                         // Angular momentum along x
//     double Lhaty = 0.0;                                         // Angular momentum along y
//     double Lhatz = 1.0;                                         // Angular momentum along z
//     //double norml = norm3(Lhatx, Lhaty, Lhatz);
//
//     // Calculate initial angle between spins
//     // All of these functions are in the PNEvolution header
//     // Don't really need to calculate this here as we are specifying it above.
//     //deltaPhiAlt = calculateDeltaPhi(s1x, s1y, s1z, s2x, s2y, s2z, Lhatx, Lhaty, Lhatz);
//     // Can now calculate the initial angle between the two spins
//     double Theta12 = calculateTheta12(s1x, s1y, s1z, s2x, s2y, s2z);
//
//     // similarly:
//     // check the function to do this
//     //double theta1check2 = calculateTheta(s1x, s1y, s1z, Lhatx, Lhaty, Lhatz);
//     //double theta2check2 = calculateTheta(s2x, s2y, s2z, Lhatx, Lhaty, Lhatz);
//
//     // Should check that the current semimajor axis a is not already less than 1000M; if it is start at its current value. Will presumably be more important when you have common envelope evolution as well.
//
//     // Initial orbit at a = 1000M
//     double currentSemiMajorAxis = getSemiMajorAxisPrime() * AU;                 // get current semiMajorAxis and convert to G = c = 1 units (m)
//
//     double a0 = 0;
//     if(currentSemiMajorAxis < 1000 * Mtot){
//         a0  = currentSemiMajorAxis;
//         std::cout << "this binary started PN evolution less than 1000M : " << getID() << std::endl;
//     }
//     else{
//         a0  = 1000*Mtot;                                        // Initial separation to start evolving from
//     }
//     //double e0       = 0;                                                // Initial eccentricity (should always be about 0)
//     //double L0       = angMomentumNat2(m1, m2, a0, e0);                  // Initial angular momentum
//     double nu0      = orbitalVelocity(Mtot, a0);                        // Initial orbital velocity
//     double omega0   = orbitalFrequencyFuncNu(Mtot, nu0);                // Initial orbital frequency
//     double omegaGW0 = gravitationalWaveFreqFuncNu(Mtot, nu0);           // 2.0 * omega0 -- Initial gravitational wave frequency
//
//     ///////////////////////////////////////////////////////////////////////////
//     //                      SOLVE THE ODES
//     ///////////////////////////////////////////////////////////////////////////
//
//     typedef odeint::runge_kutta_dopri5< state_type > dopri;         // Does this ever get used?
//
//     const double atol = 1E-8;               // Set absolute tolerance
//     const double rtol = atol;               // Set relative tolerance
//
//     // Create a version of the petersEvolutions structure that contains the derivatives - m1 m2 are the masses
//     spinPrecession p(m1, m2);
//
//     // Create state type for initial conditions - S1, S2, Lhat, nu
//     state_type x = {{ s1x, s1y, s1z, s2x, s2y, s2z, Lhatx, Lhaty, Lhatz, nu0 }};
//
//     // Set initial step size -- have to guess this -- needs to approximately match the precessional timescale
//     double dtinsi = 1;                                  // in s
//     //double dtinsi = 1e+5;                           // commonEnvelopein s - was used for the longer radiation reaction evolution
//     double dtinit = dtinsi * c;                         // multiply by c to convert to G = c = 1 units.
//     double dt = dtinit;                                 // Store the initial dt and copy to a working variable
//
//     // Make adaptive
//     auto adopri = make_dense_output(atol, rtol, odeint::runge_kutta_dopri5<state_type>());
//
//     // make sure to reset t
//     double t = 0;
//
//     adopri.initialize(x, t, dt);
//
//     size_t loopcount = 0;                    // counter for number of steps
//
//     // Variables for output - how many of these do you actually need/want at the moment?
// //    double a            = a0;
// //    double e            = e0;
// //    double L            = L0;
//     double nu           = nu0;                                          // Orbital velocity
//     double omega        = omega0;                                       // Orbital frequency
//     double omegaGW      = omegaGW0;                                     // Gravitational wave frequency
//     double theta1out    = theta1;
//     double theta2out    = theta2;
//     double theta12out   = Theta12;
//     double dPhiout      = deltaPhi;
//
//     // More temporary spin variables so that can output spin at each timestep (could use ones from above?)
//     double s1xt = 0.0;
//     double s1yt = 0.0;
//     double s1zt = 0.0;
//
//     double s2xt = 0.0;
//     double s2yt = 0.0;
//     double s2zt = 0.0;
//
//     double Lhatxt = 0.0;
//     double Lhatyt = 0.0;
//     double Lhatzt = 0.0;
//
//     // Continue to integrate spins whilst f_GW < 10Hz i.e. whilst nu < (2 pi 5M/c)^(1/3) = tenhertz
//     double tenhertz = pow((2.0 * M_PI * 5.0 * Mtot / c), (1.0/3.0));
//
//     // Should I maybe output some more numbers -- initial angles and dPhi and theta12?
//
//     // Integrate spins
//     while(adopri.current_state()[9] < tenhertz){
//
//         // Take a step with the step size still specified by the adaptive stepper (state and time managed internally)
//         adopri.do_step(p);
//
//         // See what the new time and timestep are (i.e. t_n+1 = t_n + dt)
//         t   = adopri.current_time();
//         dt  = adopri.current_time_step();
//
//         // Grab the current spins and angular momentum vectors
//         s1xt = adopri.current_state()[0];
//         s1yt = adopri.current_state()[1];
//         s1zt = adopri.current_state()[2];
//
//         s2xt = adopri.current_state()[3];
//         s2yt = adopri.current_state()[4];
//         s2zt = adopri.current_state()[5];
//
//         Lhatxt = adopri.current_state()[6];
//         Lhatyt = adopri.current_state()[7];
//         Lhatzt = adopri.current_state()[8];
//
//         // Grab the current orbital velocity, and calculate the orbital and GW frequencies
//         nu      = adopri.current_state()[9];
//         omega   = orbitalFrequencyFuncNu(Mtot, nu);
//         omegaGW = gravitationalWaveFreqFuncNu(Mtot, nu);
//
//         // Calculate theta1out, theta2out, theta12out and dPhiout at each time step
//         theta1out   = calculateTheta(s1xt, s1yt, s1zt, Lhatxt, Lhatyt, Lhatzt);
//         theta2out   = calculateTheta(s2xt, s2yt, s2zt, Lhatxt, Lhatyt, Lhatzt);
//         theta12out  = calculateTheta12(s1xt, s1yt, s1zt, s2xt, s2yt, s2zt);
//         dPhiout     = calculateDeltaPhi(s1xt, s1yt, s1zt, s2xt, s2yt, s2zt, Lhatxt, Lhatyt, Lhatzt);
//
//         loopcount++;
//
//     }
//
//     // Generate view direction iota isotropically
//     double cosiota = gsl_rng_uniform(r)*2 - 1.0;                       // generate iota uniform in cos theta
//
//     // Check  -1.0 < cosiota < 1.0
//     if(cosiota > 1.0){
//         cosiota = 1.0;
//     }
//     else if(cosiota < -1.0){
//         cosiota = -1.0;
//     }
//
//     double iota = acos(cosiota);
//
//     // Convert the final spins to the LAL radiation frame
//     convertToLAL(s1xt, s1yt, s1zt, s2xt, s2yt, s2zt, Lhatxt, Lhatyt, Lhatzt, iota);
//
//     // Return (update binary) the final spins, angles (which )and velocity
//     m_theta1    = theta1out;
//     m_theta2    = theta2out;
//     m_theta12   = theta12out;
//     m_deltaPhi  = dPhiout;
//     m_nu        = nu;                                               // Return the final orbital velocity at 10Hz
//     m_iota      = iota;
//
//     // Presumably these are the full vectors not unit vectors so may need to post-process
//     m_S1x = s1xt;                           // save x component of spin 1
//     m_S1y = s1yt;                           // save y component of spin 1
//     m_S1z = s1zt;                           // save z component of spin 1
//
//     m_S2x = s2xt;                           // save x component of spin 2
//     m_S2y = s2yt;                           // save y component of spin 2
//     m_S2z = s2zt;                           // save z component of spin 2
//
// //    // DEBUG: Check magnitudes after conversion
// //    std::cout << "After conversion:" << std::endl;
// //    std::cout << "norm S1 = " << norm3(s1xt, s1yt, s1zt) << std::endl;
// //    std::cout << "chi1 * m1 * m1 = " << chi1 * m1 * m1 << std::endl;
// //    std::cout << "norm S2 = " << norm3(s2xt, s2yt, s2zt) << std::endl;
// //    std::cout << "chi2 * m2 * m2 = " << chi2 * m2 * m2 << std::endl;
//
//     // These should still be able to be unit vectors though. What went wrong?
//     // Check the convertToLAL function and the randomSpinsFile.
//
// }

double BinaryStar::calculateAngularMomentum(double a, double e, 
	double m1, double m2, double R1, double R2, double w1, double w2, double w, double k_s1, double k_s2){

	double 	Is1 = k_s1*m1*pow(R1*RsolToAU,2.0);
	double 	Is2 = k_s2*m2*pow(R2*RsolToAU,2.0);
    double  Jorb = ((m1*m2)/(m1+m2))*sqrt(G1*(m1+m2)*a*(1.0-e*e));
	double 	totalAngularMomentum = (Is1*w1)+(Is2*w2)+Jorb;
	bool 	debugging = false;
	
	if(debugging){
		std::cout << "Is1, w1:\t" << Is1 << "\t" << w1 << std::endl; 
		std::cout << "Is2, w2:\t" << Is2 << "\t" << w2 << std::endl; 
		std::cout << "Jorb:\t" << Jorb << std::endl; 
		std::cout << "totalAngularMomentum:\t" << totalAngularMomentum << std::endl; 
	}
    	return	totalAngularMomentum;
}

//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
double BinaryStar::calculateTotalEnergy(double a,double m1,double m2, double R1, double R2, double w1, double w2, double w, double k_s1, double k_s2){
//  Created by Alejandro Vigna-Gomez on 11/2015.
//  Copyright (c) 2015 Alejandro Vigna-Gomez. All rights reserved.//
   // 	double const G1		= 4.0*pi*pi;				// graviational constant 'G' in [G] = AU^3 * yr-1 * Msol-1
//	double const RsolToAU	= 0.0046491; 				// Constant for going from Rsol to AU
	double 	Is1 = k_s1*m1*pow(R1*RsolToAU,2.0);
	double 	Is2 = k_s2*m2*pow(R2*RsolToAU,2.0);
	double	d1 = a*m2/(m1+m2);
	double	d2 = a*m1/(m1+m2);
	double	Itot = (m1*pow(d1,2.0))+(m2*pow(d2,2.0));
	double	totalEnergy = (0.5*Is1*w1*w1)+(0.5*Is2*w2*w2)+(0.5*Itot*w*w)-(G1*m1*m2/(a));
	bool	debugging = false;
	
	if(debugging){
		std::cout << "Is1, w1:\t" << Is1 << "\t" << w1 << std::endl; 
		std::cout << "Is2, w2:\t" << Is2 << "\t" << w2 << std::endl; 
		std::cout << "totalEnergy:\t" << totalEnergy << std::endl; 
	}
	
	return	totalEnergy;	// E
}

double BinaryStar::aCorrectionMassLoss(double a, double m1, double m2, double m1final, double m2final){
    // New semi-amjor axis for circuralized orbit
	// StarTrack assumes 	a*(m1+m2) = constant
	// vdSluys assumes		

    double  M   = m1+m2;
    double  Mfinal  = m1final + m2final;

    return a/(2.0-(M/Mfinal));

}
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//


//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
/*
double BinaryStar::k_definition(double Mass, double Radius, double RZAMS){
    //Define gyration radius 'k' using fit from de Mink et al. 2013, calling k_definition function
    //  Created by Alejandro Vigna-Gomez on 11/2015.
    //  Copyright (c) 2015 Alejandro Vigna-Gomez. 
	
	// ALEJANDRO - 14/02/2017 - This should a function of the star class, not the binary class. Ilya changed that in his own version of COMPAS.
    double k0,c,C;
    double mass=log10(Mass);
    double ratio=Radius/RZAMS;

    if(mass<=1.3){c=0;}
    else{c=-0.055*pow(mass-1.3,2.0);}

    if(mass<=0){C=-2.5;}
    else if((mass>0)&&(mass<=0.2)){C=-2.5+5*mass;}
    else{C=-1.5;}

    k0=c+std::min(0.21,std::max(0.09-0.27*mass,0.037+0.033*mass));
    return ((k0-0.025)*pow(ratio,C))+(0.025*pow(ratio,-0.1));
}
 * */
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//


double BinaryStar::orbitalEnergy(double mu, double M, double a){

    // Calculate the orbital energy (not in J, needs a factor of G) using the vis-viva equation
    // mu is the reduced mass in MSOL
    // M is the total mass in MSOL
    // a is the semi-major axis in AU

	// Added the G factor to have correct units. 
        // In most (if not all) functions that call this function, shouldn't affect as they check ratios.
	// G = 4*pi*pi
	// [G] = AU^3 yr^-2 Msol^-2
	// [E] = AU^2 yr^-2 Msol

    return - (G1*mu*M)/(2.0*a);

}

double BinaryStar::orbitalAngularMomentum(double mu, double M, double a){

    // Calculate the total angular momentum of the binary over sqrt(G)
    // mu is the reduced mass in MSOL
    // M is the total mass in MSOL
    // a is the semi-major axis in AU
	// J/sqrt(G) = mu*sqrt(M*a)

	// G1 = 4*pi*pi
	// [G] = AU^3 yr^-2 Msol^-2
	// [J] = AU^2 yr^-1 Msol^3/2

    return mu*sqrt(G1*M*a);

}

//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
// Begin Tides
void BinaryStar::tides(programOptions const &options){
//  Created by Alejandro Vigna-Gomez on 11/2015.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Tides function was initially developed comparing the previous timestep to the new timestep and evolving
// the system using our simple tides prescription (see COMPAS notes). We later noticed this created some
// problems when we had events such as CE or MT. The way to implement tides now is to solve them at the end of
// binary evolution, tidally locking them. This is still work in progress. ALEJANDRO - 04/10/2016
//Initial parameters
    double m1   	= star1.m_Mass;                                     // Primary mass in Msol
    double m2   	= star2.m_Mass;                                     // Secondary mass in Msol
    double Mtot 	= m1+m2;                                            // Total mass in Msol
    double mu 		= (m1 * m2)/(Mtot);                                 // Initial reduced mass in Msol, same timestep
    double m1prev  	= m1;                                  				// Primary mass in Msol, same timestep
    double m2prev  	= m2;                                  				// Primary mass in Msol, same timestep
    double Mtotprev = Mtot;                                    			// Total mass in Msol, same timestep
    double muprev 	= (m1prev * m2prev)/(Mtotprev);                     // Initial reduced mass in Msol, same timestep
    double R1		= star1.m_Radius * RsolToAU;                        // Primary radius in AU
    double R2		= star2.m_Radius * RsolToAU;                        // Secondary radius in AU
    double R1prev	= R1;			                    				// Primary radius in AU 1 previous timestep
    double R2prev	= R2; 			                                    // Secondary radius in AU 1 previous timestep

    double R1_ZAMS	= getRZAMS1() * RsolToAU;                           // Primary radius at ZAMS in AU, where 1 Rsol = 0.0046491 AU
    double R2_ZAMS	= getRZAMS2() * RsolToAU;                           // Secondary radius at ZAMS in AU, where 1 Rsol = 0.0046491 AU
    double w1		= star1.m_omega;                                    // Primary angular velocity in yr-1 units
    double w2		= star2.m_omega;                                    // Secondary angular velocity in yr-1 units

    double a		= m_SemiMajorAxisPrime;                             // Semi-major axis in default units, AU
    double d1prev	= a*m2prev/(m1prev+m2prev);                         // initial distance to CM from star 1 in AU
    double d2prev 	= a*m1prev/(m1prev+m2prev);                         // initial distance to CM from star 2 in AU
    double w_orb	= m_orbitalVelocityPrime;                           // initial orbital speed of the system in yr-1 units
    double e		= m_EccentricityPrime;

	double m1core	= star1.m_coreMass;
	double m2core	= star2.m_coreMass;	
	double m1env	= star1.m_envMass;
	double m2env	= star2.m_envMass;	
	int	   sType1	= star1.m_stellarType;
	int	   sType2	= star2.m_stellarType;

	Star	starCopy1=star1;
	Star	starCopy2=star2;
	
	double	Rc1	=	starCopy1.radiusRemnantStarAfterLosingEnvelope()*RsolToAU;	// Radius of remnant of primary after stripping its envelope, in AU
	double	Rc2	=	starCopy2.radiusRemnantStarAfterLosingEnvelope()*RsolToAU;	// Radius of remnant of secondary after stripping its envelope, in AU
	
//Define gyration radius 'k' using fit from de Mink et al. 2013, calling k_definition function
	double k_s1 = k_definition(R1_ZAMS,R1,m1,sType1);
	double k_s2 = k_definition(R2_ZAMS,R2,m2,sType2);
    double k_s1prev = k_s1;
    double k_s2prev = k_s2;

//Calculate moments of inertia
	double	I_s1prev = momentOfInertia(starCopy1, m1, m1core, m1env, R1, Rc1, R1_ZAMS, sType1);
	double	I_s2prev = momentOfInertia(starCopy2, m2, m2core, m2env, R2, Rc2, R2_ZAMS, sType2);
	double 	I_totprev=(m1prev*d1prev*d1prev)+(m2prev*d2prev*d2prev);

    double  I_s1 = I_s1prev;
    double  I_s2 = I_s2prev;

    double 	L = calculateAngularMomentum(a,e,m1prev,m2prev,R1prev,R2prev,w1,w2,w_orb,k_s1prev,k_s2prev);
    double  E = calculateTotalEnergy(a, m1prev, m2prev, R1prev, R2prev, w1, w2, w_orb, k_s1prev, k_s2prev);

    bool    debugging 		= false; // true;

    if(debugging){
		std::cout << "Beggining of tides function" << std::endl;
		printingBinaryVariables();
		/*
		std::cout << "ap (AU)\t\t"  << a << std::endl;
        std::cout << "worb (yr-1)\t" << w_orb << std::endl;
        std::cout << "w (yr-1)\t" << m_orbitalVelocityPrime << std::endl;
        std::cout << "w1 (yr-1)\t" << w1 << std::endl;
        std::cout << "w2 (yr-1)\t" << w2 << std::endl;
        std::cout << "wb1 (yr-1)\t" << star1.m_omegaBreak << std::endl;
        std::cout << "wb2 (yr-1)\t" << star2.m_omegaBreak << std::endl;
        std::cout << "Mtot (Msol)\t" << Mtot<< std::endl;
        std::cout << "MtotPrev (Msol)\t" << Mtotprev << std::endl;
        std::cout << "m1p (Msol)\t" << m1prev << std::endl;
        std::cout << "m1 (Msol)\t" << m1 << std::endl;
        std::cout << "m2p (Msol)\t" << m2prev << std::endl;
        std::cout << "m2 (Msol)\t" << m2 << std::endl;
        std::cout << "R1p (Rsol)\t" << R1prev << std::endl;
        std::cout << "R1 (Rsol)\t" << R1 << std::endl;
        std::cout << "R2p (Rsol)\t" << R2prev << std::endl;
        std::cout << "R2 (Rsol)\t" << R2 << std::endl;
        std::cout << "L \t" << L << std::endl;
        std::cout << "E \t" << E << std::endl;
		*/
		std::cout << "Rc1 [AU]:\t" << Rc1 << std::endl;
		std::cout << "Rc2 [AU]:\t" << Rc2 << std::endl;
        std::cout << "k1p \t" << k_s1prev << std::endl;
        std::cout << "k1 \t" << k_s1 << std::endl;
        std::cout << "k2p \t" << k_s2prev << std::endl;
        std::cout << "k2 \t" << k_s2 << std::endl;
        std::cout << "Is1p \t" << I_s1prev << std::endl;
        std::cout << "Is1 \t" << I_s1 << std::endl;
        std::cout << "Is2p \t" << I_s2prev << std::endl;
        std::cout << "Is2 \t" << I_s2 << std::endl;
        std::cout << "Itotp\t" << I_totprev << std::endl;

    }


	if(a <= 0.0){	// ALEJANDRO - 24/01/2017 - Checking for sensible quantities. Add at will.
		if(debugging){std::cout << m_randomSeed <<  "\tError in tides. Semi-major axis <= 0.0. a = " << a << std::endl;}
		// Do nothing
	}
	else if(m_stellarMerger == true){
		if(debugging){std::cout << m_randomSeed <<  "\tA failed common envelope just took place, leading to a stellar merger. Tides are not evolved in this case. a = " << a << std::endl;}
		// Do nothing		
	}
    else if(options.tidesPrescription == TIDES_PRESCRIPTION_NONE){
        // Do nothing
    }
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	else if(options.tidesPrescription == TIDES_PRESCRIPTION_LOCKED_ENERGY){
            std::cerr << m_randomSeed <<  "\tError in tides function TIDES_PRESCRIPTION_LOCKED_ENERGY. Prescription no longer supported nor implemented." << std::endl;
    }

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    else if(options.tidesPrescription == TIDES_PRESCRIPTION_LOCKED_ANG_MOMENTUM){


        // Number of coefficients (is order of polynomial + 1)
        const size_t polynomialOrder = 4;
        const size_t nCoefficients = polynomialOrder + 1;
        const size_t nSolutions = 2 * polynomialOrder; // real part + imaginary part for each solution

        // Calculate the coefficients of the polynomial
        double spin1 = I_s1prev*w1;                                     // Initial spin angular momenta of object 1
        double spin2 = I_s2prev*w2;                                     // Initial spin angular momenta of object 2

        double A = sqrt(G1 * Mtot) * mu;
        double B = (w_orb * a*a * muprev) + spin1 + spin2;
        double C = (I_s1 * sqrt(G1 * Mtot)) + (I_s2 * sqrt(G1 * Mtot));
		
		// Polynomial of the form:	c_4*x^4 + c_3*x^3 + c_2*x^2 + c_1*x + c_0 = 0 
        double c_4 = A;
        double c_3 = B;
        double c_2 = 0.0;
        double c_1 = 0.0;
        double c_0 = C;

		if(debugging)
			std::cout << "c_4, c_3, c_0:\t" << c_4 << "\t" << c_3 << "\t" << c_0 << std::endl;
		
        // Package these coefficients into an array
        double coefficients[nCoefficients] = {c_0, c_1, c_2, c_3, c_4};
        double solutions[nSolutions];
        double angularMomentumSolutions[polynomialOrder];
        double orbitalEnergySolutions[polynomialOrder];

        // Do GSL rootfinding magic
        gsl_poly_complex_workspace *workspace = gsl_poly_complex_workspace_alloc (nCoefficients);
        gsl_poly_complex_solve(coefficients, nCoefficients, workspace, solutions);
        gsl_poly_complex_workspace_free(workspace);

        double  aSolution[3];
        double  wsync[3];
        double  E1[3];
        double  L1[3];
        double  delta_E1[3];
        double  delta_L1[3];
		double	errorPermitted	= 0.5;
        bool    chooseValues = false;

        // DEBUG: Output the solutions
        for (int i = 0; i < polynomialOrder; i++){
            aSolution[i]=std::pow(solutions[2*i],2.0);
            wsync[i] = sqrt((m1+m2)*G1*pow(aSolution[i],-3.0));
            E1[i]	= calculateTotalEnergy(aSolution[i], m1, m2, R1, R2, wsync[i], wsync[i], wsync[i], k_s1, k_s2);
            L1[i]	= calculateAngularMomentum(aSolution[i], e, m1, m2, R1, R2, wsync[i], wsync[i], wsync[i], k_s1, k_s2);	// To verify for conservation of angular momentum
            delta_E1[i] = std::abs(E-E1[i]);
			delta_L1[i] = std::abs(L-L1[i]);
			
            if(debugging)
                std::cout << "delta_E, delta_L, a:\t" << delta_E1[i] << "\t" << delta_L1[i] << "\t" << aSolution[i] << std::endl;


            if((delta_E1[i] < errorPermitted*std::abs(E))&&(solutions[2*i+1] == 0)){
                    star1.m_omegaTidesIndividualDiff = - star1.m_omegaPrev + wsync[i];
                    star2.m_omegaTidesIndividualDiff = - star2.m_omegaPrev + wsync[i];
                    m_aTidesDiff        = - m_SemiMajorAxisPrev + aSolution[i];
                    m_omegaTidesDiff 	= - m_orbitalVelocityPrev + wsync[i];
                    m_omegaTides        = wsync[i];
                    chooseValues = true;
                    break;
            }
			
			

        }

        if(chooseValues == false){
			// Check for solutions which maybe not be exact, but close to instant tidal locking
	        for (int i = 0; i < polynomialOrder; i++){
				if((delta_L1[i] < errorPermitted*L)&&(solutions[2*i+1] == 0)){
                    star1.m_omegaTidesIndividualDiff = - star1.m_omegaPrev + wsync[i];
                    star2.m_omegaTidesIndividualDiff = - star2.m_omegaPrev + wsync[i];
                    m_aTidesDiff        = - m_SemiMajorAxisPrev + aSolution[i];
                    m_omegaTidesDiff 	= - m_orbitalVelocityPrev + wsync[i];
                    m_omegaTides        = wsync[i];
                    chooseValues = true;
                    break;
				}
			}
			
	        if(chooseValues == false){	
				m_errorFlag     = true;
				std::cerr << m_randomSeed <<  "\tError in tides function TIDES_PRESCRIPTION_LOCKED_ANG_MOMENTUM. No solution to the polynomial within the limits defined by the user." << std::endl;
			}
        }

    }
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    else{
        m_errorFlag     = true;
        std::cerr << m_randomSeed <<  "\tError: Invalid tides prescription! Shouldn't get here." << std::endl;
    }

}
// End Tides
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
// Begin Magnetic Breaking
void BinaryStar::evaluateMagneticBraking(double dt, double wsync){

    star1.m_omegaMagneticBrakingDiff = 0.0;
    star2.m_omegaMagneticBrakingDiff = 0.0;
    m_angularMomentumMagneticBrakingDiff = 0.0;
	m_aMagneticBrakingDiff = 0.0;
	// ALEJANDRO - 15/02/2017 - Re-write this function
	/*
    double m1prev  = star1.m_MassPrev;                                  // Primary mass in Msol 1 timestep before
    double m2prev  = star2.m_MassPrev;                                  // Primary mass in Msol 2 timestep before
    double R1prev	= star1.m_RadiusPrev * RsolToAU;                    // Primary radius in AU 1 previous timestep
    double R2prev	= star2.m_RadiusPrev * RsolToAU;                    // Primary radius in AU 1 previous timestep
    double R1_ZAMS	= getRZAMS1() * RsolToAU;                           // Primary radius at ZAMS in AU, where 1 Rsol = 0.0046491 AU
    double R2_ZAMS	= getRZAMS2() * RsolToAU;                           // Secondary radius at ZAMS in AU, where 1 Rsol = 0.0046491 AU
    double w1prev	= wsync*MyearToyear;                                    // Primary angular velocity in Myr-1 units
    double w2prev	= wsync*MyearToyear;                                    // Secondary angular velocity in Myr-1 units

    double a	= m_SemiMajorAxisPrev;                                  // Semi-major axis in default units, AU
    double k_s1prev = k_definition(m1prev,R1prev,R1_ZAMS);              //
    double k_s2prev = k_definition(m2prev,R2prev,R2_ZAMS);              //


    std::cout << "Begin MB"<< std::endl;
    std::cout << "a (AU)\t\t"  << a << std::endl;
    std::cout << "w1 (yr-1)\t" << w1prev << std::endl;
    std::cout << "w2 (yr-1)\t" << w2prev << std::endl;
    std::cout << "m1p (Msol)\t" << m1prev << std::endl;
    std::cout << "m2p (Msol)\t" << m2prev << std::endl;
    std::cout << "R1p (Rsol)\t" << R1prev << std::endl;
    std::cout << "R2p (Rsol)\t" << R2prev << std::endl;
    std::cout << "k1p \t" << k_s1prev << std::endl;
    std::cout << "k2p \t" << k_s2prev << std::endl;
    std::cout << "dt \t" << dt << std::endl;


    star1.m_omegaMagneticBrakingDiff = 0.0;
    star2.m_omegaMagneticBrakingDiff = 0.0;
    m_angularMomentumMagneticBrakingDiff = 0.0;
    m_aMagneticBrakingDiff = magneticBrakingSemiMajorAxisRepetto(k_s1prev, m1prev, R1prev, m2prev, a, dt) + magneticBrakingSemiMajorAxisRepetto(k_s2prev, m2prev, R2prev, m1prev, a, dt);
    std::cout << "dw1MB:\t" << star1.m_omegaMagneticBrakingDiff << std::endl;
    std::cout << "dw2MB:\t" << star2.m_omegaMagneticBrakingDiff << std::endl;
    std::cout << "daMB:\t" << m_aMagneticBrakingDiff << std::endl;
    std::cout << "dLMB:\t" << m_angularMomentumMagneticBrakingDiff << std::endl;
    std::cout << "End MB"<< std::endl;
    std::cout << std::endl;
	 * */
}
// End Magnetic Breaking
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//


//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
// Begin Timestep determination
double BinaryStar::chooseTimestep(programOptions const &options, double dt, double dtPrev){
// Function made to calculate the proper timestep, according to star nuclear timescales, 
// angular momentum loss timescales or the specific
// case of the shortest timescales with RLOF. Based on section 2.8 of Hurley et al. (2002)
    double  orbitalTimestep=0;  // From eq 88 Hurley et al. (2002)
    double  largeTimestep=dt;   // From nuclear evolution SSE
    double  chosenTimestep=dt;  // Returned value at the end
    double  dtRL1=0;            // Why do I have to initialize value?
    double  dtRL2=0;            // Why do I have to initialize value?
    double  km;                 // From eq 91 Hurley et al. (2002)
    bool  debugging = false;
    double	min_timestep = 1E-4;
    //debugging = true;
    if(debugging){std::cout << "Before timestep function: \t" << dt << std::endl;}

	// ALEJANDRO - 02/12/2016 - This function still have a lot of problems and flaws. The core is, to some extent, there, but doesn't deal well with some edge cases. Re-visit or re-write.
    //
    /*
        if(m_TotalAngularMomentumPrev!=m_TotalAngularMomentumPrime){
            orbitalTimestep= (0.02*m_TotalAngularMomentumPrev/std::abs(m_TotalAngularMomentumPrime-m_TotalAngularMomentumPrev))*(dtPrev*MyearToyear); // Eq 88 Hurley et al. (2002)
            largeTimestep=std::min(orbitalTimestep,dt);
        }
        else{
            largeTimestep=dt;  // There is no change in angular momentum, the decision is made acccording to star nuclear timescales.
        }


        // Timestep calculation for RLOF/MassTransfer
        //Calculate RLOF timestep for donor star 1
        if(m_SemiMajorAxisPrev*rocheLobeRadius(star1.m_Mass,star2.m_Mass) <= star1.m_Radius*RsolToAU){
        //    if(m_SemiMajorAxisPrev*rocheLobeRadius(star1.m_MassPrev,star2.m_MassPrev) <= star1.m_RadiusPrev*RsolToAU){
             std::cout<<"RLOF from star 1"<<std::endl;
            if(star1.flagRLOF==false){
                dtRL1=largeTimestep*pow(10.0,-3.0);
                star1.flagRLOF = true;
            }
            else{
                km=std::min(2*pow(10.0,-3.0)*largeTimestep*m_orbitalVelocityPrime/(2*pi),(0.005*star1.m_Mass)/std::abs(star1.m_Mass-star1.m_MassPrev)); // Eq. (92) Hurley et al. (2002)
                dtRL1=km*2*pi/m_orbitalVelocityPrime; // Eq. (90) Hurley et al. (2002)
            }
         }
        else{
            star1.flagRLOF=false;
        }

        //Calculate RLOF timestep for donor star 2
        if(m_SemiMajorAxisPrev*rocheLobeRadius(star2.m_Mass,star1.m_Mass) <= star2.m_Radius*RsolToAU){
        //    if(m_SemiMajorAxisPrev*rocheLobeRadius(star2.m_MassPrev,star1.m_MassPrev) <= star2.m_RadiusPrev*RsolToAU){
            std::cout<<"RLOF from star 2"<<std::endl;
            if(star2.flagRLOF==false){
                dtRL2=largeTimestep*pow(10.0,-3.0);
                star2.flagRLOF = true;
            }
            else{
                km=std::min(2*pow(10.0,-3.0)*largeTimestep*m_orbitalVelocityPrime/(2*pi),(0.005*star2.m_Mass)/std::abs(star2.m_Mass-star2.m_MassPrev)); // Eq. (92) Hurley et al. (2002)
                dtRL2=km*2*pi/m_orbitalVelocityPrime; // Eq. (90) Hurley et al. (2002)
            }
        }
        else{
            star2.flagRLOF=false;
        }

        // Set timestep according to analysis
        if((star1.flagRLOF==false)&&(star2.flagRLOF==false)){
             if(options.useMassLoss == true){
                 if(((star1.m_Mass <= star1.m_Mass0*(1.0/3.0)))||(star2.m_Mass <= star2.m_Mass*(1.0/3.0))){
                     chosenTimestep = largeTimestep/128;
                 }
                 else if(((star1.m_Mass <= star1.m_Mass0*(2.0/3.0)))||(star2.m_Mass <= star2.m_Mass*(2.0/3.0))){
                     chosenTimestep = largeTimestep/64;
                 }
                 else{
                     chosenTimestep = largeTimestep/8;
                 }
             }

        //         else{
                 chosenTimestep = largeTimestep;
        //         }

        }

        else if((star1.flagRLOF==true)&&(star2.flagRLOF==false)){
            chosenTimestep = dtRL1/2;
        }
        else if ((star1.flagRLOF==false)&&(star2.flagRLOF==true)){
            chosenTimestep = dtRL2/2;
        }
        else{
            std::cout<<"Contact system has been formed"<<std::endl;
        }
    */
    // Alternative chooseTimestep function for Binary evolution


    // Track ratio of stars size to the size of its roche lobe 
    // (computed at the point of closest approach, periapsis)
//    rocheLobeTracker1 = (star1.m_RadiusPrev*RsolToAU)/(m_SemiMajorAxisPrev*(1.0 - m_Eccentricity)*rocheLobeRadius(star1.m_MassPrev,star2.m_MassPrev));
//    rocheLobeTracker2 = (star2.m_RadiusPrev*RsolToAU)/(m_SemiMajorAxisPrev*(1.0 - m_Eccentricity)*rocheLobeRadius(star2.m_MassPrev,star1.m_MassPrev));
    rocheLobeTracker1 = (star1.m_Radius*RsolToAU)/(m_SemiMajorAxisPrime*(1.0 - m_Eccentricity)*rocheLobeRadius(star1.m_Mass,star2.m_Mass));
    rocheLobeTracker2 = (star2.m_Radius*RsolToAU)/(m_SemiMajorAxisPrime*(1.0 - m_Eccentricity)*rocheLobeRadius(star2.m_Mass,star1.m_Mass));

    if(debugging){
        std::cout << "R1prev, R2prev, M1prev, M2Prev, a: " << star1.m_RadiusPrev << " " << star2.m_RadiusPrev << " " << star1.m_MassPrev<< " " << star2.m_MassPrev << " " << m_SemiMajorAxisPrev << std::endl;
        if(m_SemiMajorAxisPrev <= 0.0){
            std::cout << "Semi-major axis <= 0. No need to calculate timestep. End simulation." << std::endl;
        }
    }

//        double  tau_th  = thermalTimescale(md, envelopeD, Rd, Ld, stellarTypeDonor);                            // Eq (47) Belczynsky et al. (2008)

    double RLlowerThreshold = 0.70;           // Threshold -- Star is getting close to overflowing its Roche lobe, decrease timesteps
    double RLthreshold = 0.9;           // Threshold -- Star is getting close to overflowing its Roche lobe, decrease timesteps
    double massTransferthreshold = 0.95;
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------//
	// ALEJANDRO - 02/12/2015 - Dealing with MASSLESS_REMNANT timescale by setting it to 0.0
	if((star1.m_stellarType == MASSLESS_REMNANT)or(star2.m_stellarType == MASSLESS_REMNANT)){
		if(debugging){std::cout << m_randomSeed << "\tMassless remant calculation of timestep in chooseTimestep function" << std::endl;}
		chosenTimestep = 0.0;
	}
	else{

// JIM BARRETT - 21/11/2016 - Error in Binary evolution can only come from these if statements 
// (or the equivalent Star2 conditionals)
// STAR1
    if(rocheLobeTracker1 < RLthreshold){
        if(debugging){
            std::cout<<"RLtracker1: \t" << rocheLobeTracker1 << std::endl;
        }
        dtRL1=dt;
    }
    else if(rocheLobeTracker1 >= RLlowerThreshold and rocheLobeTracker1 < RLthreshold){
        if(debugging){std::cout<<"RLtracker1: \t" << rocheLobeTracker1 << std::endl;}
        dtRL1=dt/5;
    }
    else if(rocheLobeTracker1 >= RLthreshold and rocheLobeTracker1 < 1.0){
        if(debugging){std::cout<<"RLtracker1: \t" << rocheLobeTracker1 << std::endl;}
        dtRL1=dt/10;
    }
    else if(rocheLobeTracker1 >= massTransferthreshold){
        // Should be done in this or in previous timestep values?

        if(envelopeType(star1.m_stellarType, options, m_randomSeed)==RADIATIVE_ENVELOPE){
            if(debugging){std::cout << "S1 Radiative Envelope\n";}
            if(m_fastPhaseCaseA==false){
                // Should only do this for fastPhase aRLOF, just one timestep.
                dtRL1=dynamicalTimescale(star1.m_Mass, star1.m_Radius);
                if(debugging){std::cout << "M1, R1, Tdyn1: " << star1.m_Mass << "\t" << star1.m_Radius << "\t" << dtRL1 << std::endl;}
            }
            else{
                dtRL1=thermalTimescale(star1.m_Mass, star1.m_envMass, star1.m_Radius, star1.m_Luminosity, star1.m_stellarType);
            }
        }
        else{
            if(debugging){std::cout << "S1 Convective Envelope\n";}
            dtRL1=thermalTimescale(star1.m_Mass, star1.m_envMass, star1.m_Radius, star1.m_Luminosity, star1.m_stellarType);
        }

        if(debugging){
        std::cout << "S1 RLOF\n";
        std::cout << "RLtracker1: \t" << rocheLobeTracker1 << std::endl;
        std::cout << "dtRL1: \t" << dtRL1 << std::endl;
        }
    }
    else{
        if(debugging){
            std::cout << "RLtracker1: \t" << rocheLobeTracker1 << std::endl;
            std::cout << "R1p: \t" << star1.m_RadiusPrev << std::endl;
            std::cout << "ap: \t" << m_SemiMajorAxisPrev << std::endl;
            std::cout << "M1p, M2p: \t" << star1.m_MassPrev << "\t" << star2.m_MassPrev << std::endl;
        }

        std::cerr << m_randomSeed <<  "\tError in ChooseTimestep function, for Star1." << std::endl;
        m_errorFlag = true;
        dtRL1=dt;
    }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------//
// STAR2
    if(rocheLobeTracker2 < RLthreshold){
        if(debugging){std::cout<<"RLtracker2: \t" << rocheLobeTracker2 << std::endl;}
        dtRL2=dt;
    }
    else if(rocheLobeTracker2 >= RLlowerThreshold and rocheLobeTracker2 < RLthreshold){
        if(debugging){std::cout<<"RLtracker2: \t" << rocheLobeTracker2 << std::endl;}
        dtRL2=dt/5;
    }
    else if(rocheLobeTracker2 >= RLthreshold and rocheLobeTracker2 < 1.0){
        if(debugging){std::cout<<"RLtracker2: \t" << rocheLobeTracker2 << std::endl;}
        dtRL2=dt/10;
    }
    else if(rocheLobeTracker2 >= massTransferthreshold){
        // Should be done in this or in previous timestep values?

        if(envelopeType(star2.m_stellarType, options, m_randomSeed)==RADIATIVE_ENVELOPE){
            if(debugging){std::cout << "S2 Radiative Envelope\n";}
            if(m_fastPhaseCaseA==false){
                // Should only do this for fastPhase aRLOF, just one timestep.
                dtRL2=dynamicalTimescale(star2.m_Mass, star2.m_Radius);
                if(debugging){std::cout << "M2, R2, Tdyn2: " << star2.m_Mass << "\t" << star2.m_Radius << "\t" << dtRL2 << std::endl;}
            }
            else{
                dtRL2=thermalTimescale(star2.m_Mass, star2.m_envMass, star2.m_Radius, star2.m_Luminosity, star2.m_stellarType);
            }
        }
        else{
            if(debugging){std::cout << "S2 Convective Envelope\n";}
            dtRL2=thermalTimescale(star2.m_Mass, star2.m_envMass, star2.m_Radius, star2.m_Luminosity, star2.m_stellarType);
        }

        if(debugging){std::cout<<"RLtracker2: \t" << rocheLobeTracker2 << std::endl;}
    }
    else{
        if(debugging){
            std::cout<<"RLtracker2: \t" << rocheLobeTracker2 << std::endl;
            std::cout << "R2p: \t" << star2.m_RadiusPrev << std::endl;
            std::cout << "ap: \t" << m_SemiMajorAxisPrev << std::endl;
            std::cout << "M1p, M2p: \t" << star1.m_MassPrev << "\t" << star2.m_MassPrev << std::endl;
        }

	std::cerr << m_randomSeed <<  "\tError in ChooseTimestep function, for Star2." << std::endl;
        m_errorFlag = true;
        dtRL2=dt;
    }

    /*
    // Check when to choose timestep.
    if(options.useMassTransfer == true){
        chosenTimestep = std::min(dtRL1,dtRL2);
    }
    else{
        chosenTimestep = dt;
    }

     */
	}

    chosenTimestep = std::max(std::min(dtRL1,dtRL2),min_timestep);

    if(debugging){
        std::cout << "Star1 RLOF flag: \t" << star1.flagRLOF << std::endl;
        std::cout << "Star2 RLOF flag: \t" << star2.flagRLOF << std::endl;
        std::cout << "dtRL1: " << dtRL1 << std::endl;
        std::cout << "dtRL2: " << dtRL2 << std::endl;
        std::cout << "dt: " << dt << std::endl;
        std::cout << "After timestep function: \t" << chosenTimestep << std::endl << std::endl;
    }


    return chosenTimestep;
}
// End Timestep determination
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//



//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
// Begin Mass Transfer
void BinaryStar::MassTransfer(double dt, double dtPrev, const programOptions &options, const gsl_rng *r){

//	printingBinaryVariables();
	double JPrev 	= m_TotalAngularMomentumPrime;

    double jloss    = m_massTransferJloss;                            		// Specific angular momentum.
    double fa       = m_massTransferFractionAccreted;                 		// Fraction of material transferred. 0 < fa <= 1
    double a        = m_SemiMajorAxisPrime;                                  // Semi-major axis in default units, AU, current timestep
    double afinal   = NEVER_SET;
    double w        = m_orbitalVelocityPrime;                                // initial orbital speed of the system in yr-1 units
    double wfinal   = NEVER_SET;
    double e	    = m_EccentricityPrime;

    double md       = star1.m_Mass;                                     // Donor star mass in Msol 1 current timestep
    double mdEnv    = star1.m_envMass;                                  // Donor star envelope mass
    double envelopeD= star1.m_Mass - star1.m_coreMass;              // Donor star efective envelope mass
    double mdCore   = star1.m_coreMass;                                 // Donor star core mass in Msol
    double mdHeCore = star1.m_HeCoreMass;                               // Donor star He-core mass in Msol
    double mdCOCore = star1.m_COCoreMass;                               // Donor star CO-core mass in Msol
    double Ld       = star1.m_Luminosity;                               // Donor star luminosity in Lsol units

    double ma       = star2.m_Mass;                                     // Accretor star mass in Msol 1 current timestep
    double maEnv    = star2.m_envMass;                                  // Accretor star envelope mass
    double envelopeA= star2.m_Mass - star2.m_coreMass;              // Accretor star efective envelope mass
    double maCore   = star2.m_coreMass;                                 // Accretor star core mass in Msol
    double maHeCore = star2.m_HeCoreMass;                               // Accretor star He-core mass in Msol
    double maCOCore = star2.m_COCoreMass;                               // Accretor star CO-core mass in Msol

    double La       = star2.m_Luminosity;                               // Accretor star luminosity in current timestep in default units
    double Rd       = star1.m_Radius * RsolToAU;                        // Donor star radius in AU 1 current timestep                      -- check units
    double rRLd     = a*rocheLobeRadius(md, ma);                            // Roche lobe radius in AU, seen by donor star, 1 current timestep
    double Ra       = star2.m_Radius * RsolToAU;                        // Accretor star radius in AU 1 current timestep
	double rRLa     = a*rocheLobeRadius(ma, md);                            // Roche lobe radius in AU, seen by accretor star, 1 current timestep

    int		stellarTypeDonor    = star1.m_stellarType;                      // Stellar type of the donor according to Hurley et al. (2002)
    int    	stellarTypeAccretor = star2.m_stellarType;                      // Stellar type of the accretor according to Hurley et al. (2002)
    bool   	primaryDonor        = true;                                     // Flag to know if the donor star is the primary star
    bool   	debugging           = false;
    bool   	addMass             = false;
	bool	CEEcurrentMTepisode	= false;									// Flag for determining if there was a CEE in the current MT episode. 

    double 	zetaThermal		= NEVER_SET;
    double 	zetaNuclear		= NEVER_SET;
	double	zetaAdiabatic	= NEVER_SET;
	double	Zlob 			= NEVER_SET;
	double	ZlobNum 		= NEVER_SET;
	double	ZlobAna 		= NEVER_SET;
	double	zetaMainSequence 	= options.zetaMainSequence;
	double	zetaHertzsprungGap 	= options.zetaHertzsprungGap;
    double zetaCompare		= NEVER_SET;

    //Addition Coen 18-10-2017
    calculatingZetas(options, r);


    if(star2.flagRLOF==true){           				// Case when donor mass is star 2 instead of star 1 (m1=md by default)
        // Update values
        md       = star2.m_Mass;                                    	// Donor star mass in Msol 1 current timestep
        mdEnv    = star2.m_envMass;                                  	// Donor star envelope mass
        envelopeD= star2.m_Mass - star2.m_coreMass;           			// Donor star efective envelope mass
        mdCore   = star2.m_coreMass;                                 	// Donor star core mass in Msol
        mdHeCore = star2.m_HeCoreMass;                               	// Donor star He-core mass in Msol
        mdCOCore = star2.m_COCoreMass;                               	// Donor star CO-core mass in Msol
        Ld       = star2.m_Luminosity;                              	// Donor star luminosity in Lsol units

        maEnv    = star1.m_envMass;                                  	// Accretor star envelope mass
        envelopeA= star1.m_Mass - star1.m_coreMass;           			// Donor star efective envelope mass
        ma       = star1.m_Mass;                                    	// Accretor star mass in Msol 1 Primeious timestep
        maCore   = star1.m_coreMass;                                 	// Accretor star core mass in Msol
        maHeCore = star1.m_HeCoreMass;                               	// Accretor star He-core mass in Msol
        maCOCore = star1.m_COCoreMass;                               	// Accretor star CO-core mass in Msol

        La       = star1.m_Luminosity;                              	// Accretor star luminosity in current timestep default units
        Rd       = star2.m_Radius * RsolToAU;                       	// Donor star radius in AU 1 current timestep
        rRLd     = a*rocheLobeRadius(md, ma);                           // Roche lobe radius in AU, seen by donor star, 1 current timestep
        Ra       = star1.m_Radius * RsolToAU;                       	// Accretor star radius in AU 1 current timestep
        rRLa     = a*rocheLobeRadius(ma, md);                           // Roche lobe radius in AU, seen by accretor star, 1 current timestep

		stellarTypeDonor    = star2.m_stellarType;                      // Stellar type of the donor according to Hurley et al. (2002)
        stellarTypeAccretor = star1.m_stellarType;                      // Stellar type of the accretor according to Hurley et al. (2002)
        primaryDonor = false;                                           // Secondary star is the donor
    }

	double	J        = (md*ma)*sqrt(G1*(md+ma)*a*(1.0-(e*e)))/(md+ma);

	int	caseMT = massTransferCase(stellarTypeDonor);

	// Check for stability
    bool qCritFlag = false;
	bool massRatioUnstable = isMassRatioUnstable(stellarTypeDonor, md, stellarTypeAccretor, ma, options, m_randomSeed); // Use when fixing stability
	
	if(	options.massTransferCriticalMassRatioMSLowMass == true or
		options.massTransferCriticalMassRatioMSHighMass == true or 
		options.massTransferCriticalMassRatioHG == true or 
		options.massTransferCriticalMassRatioGiant == true or 
		options.massTransferCriticalMassRatioHeliumGiant == true or
		options.massTransferCriticalMassRatioHeliumMS == true or
		options.massTransferCriticalMassRatioHeliumHG == true or
		options.massTransferCriticalMassRatioHeliumGiant == true or
		options.massTransferCriticalMassRatioWhiteDwarf == true){
			
        qCritFlag = true;
        if(debugging){
			std::cout << "massRatioUnstable:\t" << massRatioUnstable << std::endl;
            std::cout << "qCritFlag: " << qCritFlag << std::endl;
            std::cout << "Mdonor: " << md << std::endl;
            std::cout << "Maccretor: " << ma << std::endl;
            std::cout << "q: " << ma/md << std::endl;
            std::cout << "Kd: " << stellarTypeDonor << std::endl;
            std::cout << "Ka: " << stellarTypeAccretor << std::endl;
        }
    }
	
    if((massRatioUnstable == true)and(qCritFlag == true)){
        m_dynamicalFlag = true;
        m_commonEnvelopeFlag = true;
    }
    else{
    // Check and record if the first mass transfer event is case A, B or C
    //    if(m_initiateMassTransfer == false){

        if(primaryDonor){
            if( star1.m_firstMassTransferEpisode == false){
                star1.m_initialMassTransferCase = initialMassTransferCase(stellarTypeDonor, options);
                    if(debugging){
                        std::cout << "Primary donor, MT CASE: " << star1.m_initialMassTransferCase << std::endl;
                    }
            }
            star1.m_firstMassTransferEpisode = true;
        }
        else{
            if( star2.m_firstMassTransferEpisode == false){
                star2.m_initialMassTransferCase = initialMassTransferCase(stellarTypeDonor, options);
                    if(debugging){
                        std::cout << "Secondary donor, MT CASE: " << star2.m_initialMassTransferCase << std::endl;
                    }
            }
            star2.m_firstMassTransferEpisode = true;
        }


        double  tauThermalDonor     	= thermalTimescale(md, envelopeD, Rd/RsolToAU, Ld, stellarTypeDonor);           // Thermal timescale of donor star
        double  tauThermalAccretor  	= thermalTimescale(ma, envelopeA, Ra/RsolToAU, La, stellarTypeAccretor);        // Thermal timescale of accretor star
		double	tauDynamicalDonor 		= dynamicalTimescale(md, Rd*AUToRsol);											// Dynamical timescale of donor star
		double	tauDynamicalAccretor	= dynamicalTimescale(ma, Ra*AUToRsol);											// Dynamical timescale of accretor star
        double	tauNuclearDonor			= nuclearTimescale(md, Ld);														// Nuclear timescale of donor star
        double	tauNuclearAccretor		= nuclearTimescale(ma, La);														// Nuclear timescale of accretor star
        double  thermalRateDonor    	= md/tauThermalDonor;        													// Thermal timescale of donor star
		double  thermalRateAccretor 	= ma/tauThermalAccretor;    													// Thermal timescale of accretor star
		double	dynamicalRateDonor		= md/tauDynamicalDonor;																	// Mdot donor star on the dynamical timescale
		double	dynamicalRateAccretor	= ma/tauDynamicalAccretor;																// Mdot accretor star on the dynamical timescale
		double 	rateAccretor			= NEVER_SET;
			
		// ALEJANDRO - 19/01/2017 - Quick fix for re-defining thermal rate for a star with clear core-envelope separation. This should be re-written at the begginng of the function.
		if(stellarTypeDonor == HERTZSPRUNG_GAP or stellarTypeDonor == FIRST_GIANT_BRANCH or stellarTypeDonor == CORE_HELIUM_BURNING or stellarTypeDonor == EARLY_ASYMPTOTIC_GIANT_BRANCH or stellarTypeDonor == THERMALLY_PULSING_ASYMPTOTIC_GIANT_BRANCH or stellarTypeDonor ==NAKED_HELIUM_STAR_HETZSPRUNG_GAP or stellarTypeDonor == NAKED_HELIUM_STAR_GIANT_BRANCH){
			thermalRateDonor    	= envelopeD/tauThermalDonor;
		}
		
		if(stellarTypeAccretor == HERTZSPRUNG_GAP or stellarTypeAccretor == FIRST_GIANT_BRANCH or stellarTypeAccretor == CORE_HELIUM_BURNING or stellarTypeAccretor == EARLY_ASYMPTOTIC_GIANT_BRANCH or stellarTypeAccretor == THERMALLY_PULSING_ASYMPTOTIC_GIANT_BRANCH or stellarTypeAccretor ==NAKED_HELIUM_STAR_HETZSPRUNG_GAP or stellarTypeAccretor == NAKED_HELIUM_STAR_GIANT_BRANCH){
			thermalRateAccretor 	= envelopeA/tauThermalAccretor;
		}
		
		if(options.massTransferThermallyLimitedVariation == THERMAL_RADIUS_TO_ROCHELOBE){
			tauThermalAccretor = thermalTimescale(ma, envelopeA, rRLa/RsolToAU, La, stellarTypeAccretor);        	// Thermal timescale of accretor star if radius is assumed to be same as RL
			thermalRateAccretor = ma/tauThermalAccretor;        													// Thermal rate of accretor star if radius is assumed to be same as RL
		}
        
		if(debugging){
				std::cout << "dt\t\t\t[Myrs]\t"<< dt<<std::endl;
				std::cout << "dtprev\t\t\t[Myrs]\t"<< dtPrev<<std::endl;
				std::cout << "time\t\t\t[Myrs]\t" << m_time <<std::endl;
                std::cout << "semiMajorAxis\t\t[AU]\t" << m_SemiMajorAxisPrime << std::endl;
                std::cout << "eccentricity\t\t[]\t" << m_Eccentricity << std::endl;
				std::cout << "fa\t\t\t[]\t" << fa<< std::endl;
                std::cout << "jloss\t\t\t[]\t" << jloss << std::endl;
				std::cout << "CASE [A=1,B=2,C=3]\t\t" << caseMT << std::endl;

				std::cout << "DONOR" << std::endl;
                std::cout << "StellarTypeDonor\t[]\t" << stellarTypeDonor << std::endl;
                std::cout << "RadiusDonor\t\t[AU]\t" << Rd << std::endl;
				std::cout << "RadiusDonor\t\t[Rsol]\t" << Rd*AUToRsol << std::endl;
                std::cout << "RocheLobeDonor\t\t[AU]\t" << rRLd << std::endl;
				std::cout << "ThermalRateDonor\t[Msol Myrs^-1]\t" << thermalRateDonor << std::endl;
				std::cout << "ThermalTimescaleDonor\t[Myrs]\t" << tauThermalDonor << std::endl;
				std::cout << "DynamicalTimescaleDonor\t[Myrs]\t" << tauDynamicalDonor << std::endl;
				std::cout << "md\t\t\t[Msol]\t" << md << std::endl;
                std::cout << "mdEnvelope\t\t[Msol]\t" << envelopeD << std::endl;
                std::cout << "mdHeCore\t\t[Msol]\t" << mdHeCore << std::endl;
                std::cout << "mdCOCore\t\t[Msol]\t" << mdCOCore << std::endl;
                std::cout << "mdCore\t\t\t[Msol]\t" << mdCore << std::endl;
				std::cout << "Ld\t\t\t[Lsol]\t" << Ld << std::endl;

				std::cout << "ACCRETOR" << std::endl;
                std::cout << "StellarTypeAccretor\t[]\t" << stellarTypeAccretor << std::endl;
				std::cout << "RadiusAccretor\t\t[AU]\t" << Ra << std::endl;
                std::cout << "RadiusAccretor\t\t[Rsol]\t" << Ra*AUToRsol << std::endl;
				std::cout << "RocheLobeAccretor\t[AU]\t" << rRLa << std::endl;
				std::cout << "thermalRateAccretor\t[Msol Myrs^-1]\t" << thermalRateAccretor << std::endl;
				std::cout << "ThermalTimescaleAcc\t[Myrs]\t" << tauThermalAccretor << std::endl;
                std::cout << "DynamicalTimescaleAcc\t[Myrs]\t" << tauDynamicalAccretor << std::endl;
				std::cout << "ma\t\t\t[Msol]\t" << ma << std::endl;
				std::cout << "maEnvelope\t\t[Msol]\t" << envelopeA << std::endl;
                std::cout << "maHeCore\t\t[Msol]\t" << maHeCore << std::endl;
				std::cout << "maCOCore\t\t[Msol]\t" << maCOCore << std::endl;
                std::cout << "maCore\t\t\t[Msol]\t" << maCore << std::endl;
				std::cout << "La\t\t\t[Lsol]\t" << La << std::endl;
				
				
				std::cout << "ZETAS" << std::endl;
//				std::cout << "ZetaThermalDonor\t[]\t" << zetaThermal << std::endl;
//                std::cout << "ZetaNuclearDonor\t[]\t" << zetaNuclear << std::endl;
                std::cout << "ZetaRL\t\t\t[]\t" << Zlob << std::endl;
                std::cout << "ZetaRLnum\t\t[]\t" << ZlobNum << std::endl;
                std::cout << "ZetaRLana'\t\t[]\t" << ZlobAna << std::endl;
//                std::cout << "MdotEq\t\t\t[]\t" << MeqDot << std::endl;           // Initial guess at mass loss rate
				
				std::cout << std::endl;
		}

		// Begin Mass Transfer
        m_massTransferZeroFlag = false;
        int    massTransferType = NO_MASS_TRANSFER;

//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
        if(options.massTransferPrescription == BELCZYNSKI_MASS_TRANSFER){
			// REWRITE THIS FUNCTION AND TEST
					std::cerr << "BELCZYNSKI_MASS_TRANSFER in not currently supported." << std::endl;
        }
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
        else if(options.massTransferPrescription == DEMINK_MASS_TRANSFER){


            if(debugging){std::cout << "Using deMink MT prescription\n";}
            // Mass Transfer variables for deMink prescription
            ///////////////////////////////////////////////////////////////

            if((stellarTypeDonor >= HERTZSPRUNG_GAP and stellarTypeDonor <= NAKED_HELIUM_STAR_MS) or stellarTypeDonor == NAKED_HELIUM_STAR_GIANT_BRANCH or stellarTypeDonor == NAKED_HELIUM_STAR_HETZSPRUNG_GAP){
                double coreMass;
                //Decision to make HeMS star zeta for soberman infinite (always stable)
                if(stellarTypeDonor ==NAKED_HELIUM_STAR_MS){coreMass = mdHeCore;}
                else{coreMass = mdCore;}

                if(options.commonEnvelopeZetaPrescription == ZETA_STARTRACK){
                    zetaAdiabatic = zetaThermal;
					zetaCompare = zetaAdiabatic;
                }
                else if(options.commonEnvelopeZetaPrescription == ZETA_SOBERMAN){
                    zetaAdiabatic = ZadiabaticSPH(md, coreMass);
                    zetaThermal   = 0.0;
					zetaCompare = zetaAdiabatic;
                }
                else if(options.commonEnvelopeZetaPrescription == ZETA_HURLEY){
                    zetaAdiabatic = ZadiabaticHurley2002(md, coreMass);
                    zetaThermal   = zetaAdiabatic;
                }
                else if(options.commonEnvelopeZetaPrescription == ZETA_ARBITRARY){
                    zetaAdiabatic = options.zetaAdiabaticArbitrary;
                    zetaThermal   = options.zetaThermalArbitrary;
					zetaCompare = zetaThermal;
                }
                else{
					zetaCompare = 0.0;
                    std::cout << "Error in Mass transfer. No zeta prescription selected" << std::endl;
                }

            }
			
            if(debugging){
                std::cout << "zetaAdiabatic = " << zetaAdiabatic << std::endl;
                std::cout << "zetaThermal = " << zetaThermal << std::endl;
            }

            ///////////////////////////////////////////////////////////////
			// Here we check if we have case A, B or C mass transfer, in order to calculate the proper mass acceptance rate according to the relevant timescale.
			// Case A check still pending
			if((caseMT == CASE_A) or (caseMT == CASE_B)){
				rateAccretor = massAcceptanceRate(Ra*AUToRsol, thermalRateAccretor, thermalRateDonor, stellarTypeAccretor, stellarTypeDonor, options, fa, m_randomSeed);   // Determine the maximum mass that can be accepted by the accretor according to its stellar type.
			}
			else if(caseMT == CASE_C){
				rateAccretor = massAcceptanceRate(Ra*AUToRsol, dynamicalRateAccretor, dynamicalRateDonor, stellarTypeAccretor, stellarTypeDonor, options, fa, m_randomSeed);   // Determine the maximum mass that can be accepted by the accretor according to its stellar type.
			}
			else{
				std::cerr << m_randomSeed <<  "\tNot valid mass transfer case." << std::endl;
			}

			if(options.massTransferAngularMomentumLossPrescription != ARBITRARY_MASS_TRANSFER_ANGULAR_MOMENTUM_LOSS){
			jloss	= gammaAngularMomentumLoss(md, ma, options);
			}

			if(debugging){
				std::cout << "faAfter\t\t\t[]\t" << fa << std::endl;
				std::cout << "jloss\t\t\t[]\t" << jloss << std::endl;
				std::cout << "RateAccretor\t\t[]\t" << rateAccretor << std::endl;
				std::cout << "deltaMdonor\t\t[Msol]\t" << thermalRateDonor*tauThermalDonor << std::endl;
				std::cout << "deltaMaccretor\t\t[Msol]\t" << rateAccretor*tauThermalDonor << std::endl;
				std::cout << "fa*mEnv\t\t\t[Msol]\t" << fa*envelopeD << std::endl;
			}

            ZlobNum = numericalZRocheLobe(md, ma, a, fa, jloss, 0.01, Ra, thermalRateAccretor, thermalRateDonor, stellarTypeAccretor, stellarTypeDonor, options);  // Eq (41) Belczynsky et al. (2008), using Woods et al. (2012) formula
            ZlobAna = ZRocheLobe(md, ma, fa);                           // Eq (60) Sluys (2013)
            Zlob = ZlobNum;
			
			// Addition by Coen 18-10-2017 for zeta study
			// ALEJANDRO  - 04/10/2017 - Moved this to calculate it for deMink MT.
			m_zetaRLOFanalytic= ZlobAna;
			m_zetaRLOFnumerical= ZlobNum;
	
            if(debugging){
                std::cout << "IN DEMINK" << std::endl;
                std::cout << "ZlobNum, ZlobAna, Zlob = " << ZlobNum << " " << ZlobAna << " " << Zlob << std::endl;
            }

            if(envelopeType(stellarTypeDonor, options, m_randomSeed)==RADIATIVE_ENVELOPE){
			// Need to know which is the donor star, make it lose enough mass to stay within its Roche lobe
				
				// ALEJANDRO - 20/10/2017 - Code arbitrary zeta cut for case A mass transfer. To use in BNS paper. Should be properly coded.
				if(stellarTypeDonor == MS_MORE_THAN_07){
					if(debugging){
						std::cout << "Stellartype Donor:\t" << stellarTypeDonor << std::endl;
						std::cout << "Stellartype Accretor:\t" << stellarTypeAccretor << std::endl;
						std::cout << "Mass Donor, Mass Accretor, Mass Ratio:\t" << md << "\t" << ma << "\t" << md/ma <<std::endl;
						std::cout << "zetaMainSequence, ZlobAna:\t" << zetaMainSequence << "\t" << ZlobAna << std::endl;							
					}
					
					if(zetaMainSequence < ZlobAna){
                        m_dynamicalFlag = true;
                        m_stellarMerger = true;
                        CEE_wet = true;
						CEEcurrentMTepisode = true;
					}
				}
				
                if(debugging){
                    std::cout << "Adaptive Roche Lobe Overflow here" << std::endl;
                    std::cout << "Donor star has a radiative Envelope\n";
                    std::cout << "time: " << m_time << std::endl;
                    std::cout << "fastPhaseCaseA1: " << star1.m_fastPhaseCaseA << std::endl;
                    std::cout << "fastPhaseCaseA2: " << star2.m_fastPhaseCaseA << std::endl;
                }
				
                double dM = 0.0;
                double dMaccreted = 0.0;

                if(primaryDonor==true){

					if(debugging){std::cout << "1) Primary star is the donor" << std::endl;}					
					m_massTransferTrackerHistory = STABLE_FROM_1_TO_2;

					if(star1.m_fastPhaseCaseA==false){
					
						if(debugging){std::cout << "Fast phase from primary\n" << std::endl;}
						star1.m_fastPhaseCaseA=true;
						dM=massTransferFastPhaseCaseA(star1, ma, Ra, stellarTypeAccretor, a, w, e, options, fa, jloss, thermalRateDonor, thermalRateAccretor, tauThermalDonor, m_randomSeed, m_aFastPhaseRocheLobeOverFlowErrorPrinting, massTransferFastPhaseCaseAErrorMassFlag, massTransferFastPhaseCaseAErrorRadiusFlag, r);
						dMaccreted = -dM*fa;
					}
					else{
						dM = adaptiveRocheLobeOverFlow(star1, ma, a, options, fa, jloss, m_randomSeed, m_aRocheLobeOverFlowErrorPrinting, massTransferFastPhaseCaseAErrorMassFlag, massTransferFastPhaseCaseAErrorRadiusFlag, r);
						dMaccreted = -dM*fa;
					}

					star1.m_MassTransferDiff = dM;
					star2.m_MassTransferDiff = dMaccreted;
				}
                else{
					
                    if(debugging){std::cout << "2) Secondary star is the donor" << std::endl;}
					m_massTransferTrackerHistory = STABLE_FROM_2_TO_1;

                    if(star2.m_fastPhaseCaseA==false){

						if(debugging){std::cout << "Fast phase from secondary\n" << std::endl;}
                        star2.m_fastPhaseCaseA=true;
						dM=massTransferFastPhaseCaseA(star2, ma, Ra, stellarTypeAccretor, a, w, e, options, fa, jloss, thermalRateDonor, thermalRateAccretor, tauThermalDonor, m_randomSeed, m_aFastPhaseRocheLobeOverFlowErrorPrinting, massTransferFastPhaseCaseAErrorMassFlag, massTransferFastPhaseCaseAErrorRadiusFlag, r);
                        dMaccreted = -dM*fa;
                    }
                    else{
                        dM = adaptiveRocheLobeOverFlow(star2, ma, a, options, fa, jloss, m_randomSeed, m_aRocheLobeOverFlowErrorPrinting, massTransferFastPhaseCaseAErrorMassFlag, massTransferFastPhaseCaseAErrorRadiusFlag, r);
                        dMaccreted = -dM*fa;
                    }

                    star2.m_MassTransferDiff = dM;
                    star1.m_MassTransferDiff = dMaccreted;
				}

                if(debugging){
					std::cout << "fa: " << fa << std::endl;
                    std::cout << "aRL: dM = " << dM << std::endl;
                    std::cout << "aRL: dMacc = " << dMaccreted << std::endl;
                    std::cout << "aRL: ma = " << ma << std::endl;
                    std::cout << "star1.m_MassTransferDiff: " << star1.m_MassTransferDiff << std::endl;
                    std::cout << "star2.m_MassTransferDiff: " << star2.m_MassTransferDiff << std::endl;
				}
					
                double MdDot = dM/dt;

				if(debugging){
					std::cout << "Solving orbit in MT function" << std::endl;
					std::cout << "time [Myrs]:\t" << m_time << std::endl;
					std::cout << "md [Msol]:\t" << md << std::endl;
					std::cout << "ma [Msol]:\t" << ma << std::endl;
					std::cout << "MdDot:\t" << MdDot << std::endl;
					std::cout << "a [AU]:\t" << a << std::endl;
					std::cout << "e :\t" << e << std::endl;
					std::cout << "dt [Myr]:\t" << dt << std::endl;
					std::cout << "J :\t" << J << std::endl;
					std::cout << "fa :\t" << fa << std::endl;
					std::cout << "jloss :\t" << jloss << std::endl;
					std::cout << "Ra [Rsol]:\t" << Ra << std::endl;
					std::cout << "trateD :\t" << thermalRateDonor << std::endl;
					std::cout << "rrateA :\t" << thermalRateAccretor << std::endl;
					std::cout << "dM [Msol]:\t" << dM << std::endl;
				}

				// Before solving for orbit
                if(stellarTypeAccretor <= NAKED_HELIUM_STAR_GIANT_BRANCH){
					if(debugging){std::cout << "From MT to a NON-degenerate accretor, a':\t" << afinal << std::endl;}
                    afinal = calculateMassTransferOrbitNonDegenerateAccretor(md, envelopeD, ma, fa, MdDot, a, w, e, dt, tauThermalDonor, J, jloss, Ra, thermalRateAccretor, thermalRateDonor, stellarTypeAccretor, stellarTypeDonor, options, m_randomSeed);       // Eq (33) Belczynski et al. (2008)
                }
                else if(stellarTypeAccretor >= HELIUM_WHITE_DWARF and stellarTypeAccretor <= BLACK_HOLE){
					if(debugging){std::cout << "From MT to a degenerate accretor, a':\t" << afinal << std::endl;}
                    afinal = calculateMassTransferOrbitDegenerateAccretor(md, ma, fa, MdDot, a, w, e, dt, J, jloss, Ra, thermalRateAccretor, thermalRateDonor, stellarTypeAccretor, stellarTypeDonor, options, m_randomSeed);
                }
                else{
                    std::cerr << m_randomSeed <<  "\tError in mass transfer. Massless remnant accreting mass. Return 'a'." << std::endl;
                    afinal = a;
                }

				// After solving for orbit
				if(debugging){
					std::cout << "a [Rsol]:\t" << a*AUToRsol << std::endl;
					std::cout << "afinal [Rsol]:\t" << afinal*AUToRsol << std::endl;
					std::cout << "a [AU]:\t" << a << std::endl;
					std::cout << "afinal [AU]:\t" << afinal << std::endl;
				}

				wfinal = sqrt(G1*(md+ma)*pow(afinal,-3.0));
                m_aMassTransferDiff = afinal - a;
                m_omegaMassTransferDiff = wfinal - w;
                massTransferTracker = MdDot;
            }
            else if(envelopeType(stellarTypeDonor, options, m_randomSeed)==CONVECTIVE_ENVELOPE){ // Either Case B or Case C

				if(options.forceCaseBBBCStabilityFlag){
					if(debugging){std::cout << "forceCaseBBBCStabilityFlag:\t" << options.forceCaseBBBCStabilityFlag << std::endl;}

					// ALEJANDRO - 24/08/2017 - Check for case BB or BC mass transfer; particularly for BNS project
					if(stellarTypeDonor == NAKED_HELIUM_STAR_HETZSPRUNG_GAP or stellarTypeDonor == NAKED_HELIUM_STAR_GIANT_BRANCH){
//						
						if(stellarTypeAccretor == NEUTRON_STAR or stellarTypeAccretor == BLACK_HOLE){
//						if(stellarTypeAccretor == NEUTRON_STAR){
							if(primaryDonor){
								star1.flagUSSN = true;
							}
							else{
								star2.flagUSSN = true;
							}
						}
						
						// Hard code stability
						if(options.alwaysStableCaseBBBCFlag == true){
							zetaCompare = 1.0;
							ZlobAna		= 0.0;
						}
						else{
							zetaCompare = 0.0;
							ZlobAna		= 1.0;
						}

						if(debugging){
							std::cout << "forceCaseBBBCStabilityFlag:\t" << options.forceCaseBBBCStabilityFlag << std::endl;
							std::cout << "Always stable case BB/BC Flag:\t" << options.alwaysStableCaseBBBCFlag << std::endl;
							std::cout << "Stellartype Donor:\t" << stellarTypeDonor << std::endl;
							std::cout << "Stellartype Accretor:\t" << stellarTypeAccretor << std::endl;
							std::cout << "Zcompare, ZlobAna:\t" << zetaCompare << "\t" << ZlobAna << std::endl;
						}
					}
				}
				
				// ALEJANDRO - 20/10/2017 - Code arbitrary zeta cut for case B mass transfer. To use in BNS paper. Should be properly coded.
				if(stellarTypeDonor == HERTZSPRUNG_GAP){
						zetaCompare = zetaHertzsprungGap;
						if(debugging){
							std::cout << "Stellartype Donor:\t" << stellarTypeDonor << std::endl;
							std::cout << "Stellartype Accretor:\t" << stellarTypeAccretor << std::endl;
							std::cout << "Mass Donor, Mass Accretor, Mass Ratio:\t" << md << "\t" << ma << "\t" << ma/md <<std::endl;
							std::cout << "Zcompare, ZlobAna:\t" << zetaCompare << "\t" << ZlobAna << std::endl;							
						}
				}
						
				// ALEJANDRO - 07/08/2018 - Added m_zetaStarCompare for CE study.
				// Beware as this variable may have different values, e.g. fixed for MS, fixed for HG, different for Soberman, '1' for case BB.
				m_zetaStarCompare = zetaCompare;
				
				if(zetaCompare > ZlobAna){
					
                    // Stable Mass Transfer
                    if(debugging){std::cout << "STABLE CASE B/C MT\n";}

                    // Check for conservative or non-conservative MT
					double mdEnvAccreted = envelopeD*fa;
                    double MdDot = -envelopeD/dt;

                    if(primaryDonor==true){
						m_massTransferTrackerHistory = STABLE_FROM_1_TO_2;
                        // Primary star is the donor
                        star1.m_MassTransferDiff = -envelopeD;
                        star2.m_MassTransferDiff = mdEnvAccreted;
                        star1.modifyStarAfterLosingEnvelope(star1.m_stellarType, star1.m_Mass);   // Only other interaction that adds/removes mass is winds. So, think its safe to update star here.

                    }
                    else{
						m_massTransferTrackerHistory = STABLE_FROM_2_TO_1;
                        // Secondary star is the donor
                        star2.m_MassTransferDiff = -envelopeD;
                        star1.m_MassTransferDiff = mdEnvAccreted;
                        star2.modifyStarAfterLosingEnvelope(star2.m_stellarType, star2.m_Mass);   // Only other interaction that adds/removes mass is winds. So, think its safe to update star here.

                    }
					
					if(debugging){
						std::cout << "Solving orbit in MT function" << std::endl;
						std::cout << "time [Myrs]:\t" << m_time << std::endl;
						std::cout << "md [Msol]:\t" << md << std::endl;
						std::cout << "ma [Msol]:\t" << ma << std::endl;
						std::cout << "MdDot:\t" << MdDot << std::endl;
						std::cout << "a [AU]:\t" << a << std::endl;
						std::cout << "e :\t" << e << std::endl;
						std::cout << "dt [Myr]:\t" << dt << std::endl;
						std::cout << "J :\t" << J << std::endl;
						std::cout << "fa :\t" << fa << std::endl;
						std::cout << "jloss :\t" << jloss << std::endl;
						std::cout << "Ra [Rsol]:\t" << Ra << std::endl;
						std::cout << "trateD :\t" << thermalRateDonor << std::endl;
						std::cout << "rrateA :\t" << thermalRateAccretor << std::endl;
						std::cout << "envelopeD [Msol]:\t" << envelopeD << std::endl;
					}

					// Update and solve the orbit
					if(stellarTypeAccretor <= NAKED_HELIUM_STAR_GIANT_BRANCH){
						afinal = calculateMassTransferOrbitNonDegenerateAccretor(md, envelopeD, ma, fa, MdDot, a, w, e, dt, tauThermalDonor, J, jloss, Ra, thermalRateAccretor, thermalRateDonor, stellarTypeAccretor, stellarTypeDonor, options, m_randomSeed);       // Eq (33) Belczynski et al. (2008)
					}
                    else if(stellarTypeAccretor >= HELIUM_WHITE_DWARF and stellarTypeAccretor <= BLACK_HOLE){
						afinal = calculateMassTransferOrbitDegenerateAccretor(md, ma, fa, MdDot, a, w, e, dt, J, jloss, Ra, thermalRateAccretor, thermalRateDonor, stellarTypeAccretor, stellarTypeDonor, options, m_randomSeed);
                    }
                    else{
						std::cerr << m_randomSeed <<  "\tError in mass transfer. Massless remnant accreting mass. Return 'a'." << std::endl;
                        afinal = a;
					}

                    wfinal = sqrt(G1*(md+ma)*pow(afinal,-3.0));
                    m_aMassTransferDiff = afinal - a;
                    m_omegaMassTransferDiff = wfinal - w;

                    if(debugging){
                        std::cout << "semiMajorAxis\t\t[AU]\t" << a << std::endl;
                        std::cout << "semiMajorAxisFinal\t[AU]\t" << afinal << std::endl;
                    }

                    massTransferTracker = envelopeD/dt;
					
					// Check for stable mass transfer after any CEE
					if(m_commonEnvelopeOccuredAtLeastOnce==true){
						if((m_massTransferTrackerHistory == STABLE_FROM_2_TO_1)or(m_massTransferTrackerHistory == STABLE_FROM_1_TO_2)){
							m_stableRLOFafterCEE = true;
								if(debugging){
									std::cout << "m_commonEnvelopeOccuredAtLeastOnce:\t" << m_commonEnvelopeOccuredAtLeastOnce << std::endl;
									std::cout << "m_massTransferTrackerHistory:\t" << m_massTransferTrackerHistory << std::endl;
									std::cout << "m_stableRLOFafterCEE:\t" << m_stableRLOFafterCEE << std::endl;
								}
						}
					}
                }
                else{
                    // Unstable Mass Transfer
                    if(debugging){
                        std::cout << "UNSTABLE MT\n";
                        std::cout << "zetaThermal\t\t[]\t" << zetaThermal << std::endl;
                        std::cout << "ZlobAna\t\t\t[]\t" << ZlobAna << std::endl;
                        std::cout << "Zlobnum\t\t\t[]\t" << ZlobNum << std::endl;
                    }

                    if (stellarTypeDonor < HERTZSPRUNG_GAP) {    // How to deal with CEE here? Just worry about donor?
                        if(debugging){std::cout << "WET MERGER" << std::endl;}
                        m_dynamicalFlag = true;
                        m_stellarMerger = true; // ?
                        CEE_wet = true;
						CEEcurrentMTepisode = true;
                    }
                    else{
                        if(debugging){std::cout << "CE/MERGER" << std::endl;}
                        m_dynamicalFlag = true;
                        m_commonEnvelopeFlag = true;
						CEEcurrentMTepisode = true;						
                    }
                }
            }
            else{
                if(debugging){std::cout << "Time: " << m_time << std::endl;}
                std::cerr << m_randomSeed <<  "\tError: Mass transfer from BH, NS or massless remmant. Shouldn't get here." << std::endl;
                m_errorFlag = true;
            }
        }
        else{
            std::cerr << m_randomSeed <<  "\tError: No mass transfer prescription chosen. Error.";
        }
    }
	
	// Check for recycled pulsars. Not considering CEE as a way of recycling NSs.
	if(CEEcurrentMTepisode == false){
		if(stellarTypeAccretor == NEUTRON_STAR){
			if(primaryDonor == true){
				// The primary is the donor, therefore the secondary is the NS accretor and recycled pulsar
				star2.flagRecycledNS 	= true;
				star1.flagRLOFontoaNS	= true;
			}				
			else{
				// The secondary is the donor, therefore the secondary is the NS accretor and recycled pulsar
				star1.flagRecycledNS 	= true;
				star2.flagRLOFontoaNS	= true;
			}
		        
                        debugging = false;
			if(debugging){
				std::cout << "CEEcurrentMTepisode:\t" << CEEcurrentMTepisode << std::endl;
				std::cout << "Checking for recycled pulsars. StellarType Donor, stellartype accretor:\t" << stellarTypeDonor << ", " << stellarTypeAccretor << std::endl;
				std::cout << "PrimaryDonor:\t" << primaryDonor << std::endl;
				std::cout << "flagRecycledNS1, flagRecycledNS2:\t" << star1.flagRecycledNS << ", " << star2.flagRecycledNS << std::endl;
				std::cout << "flagRLOFontoaNS1, flagRLOFontoaNS2:\t" << star1.flagRLOFontoaNS << ", " << star2.flagRLOFontoaNS << std::endl;
				std::cout << "massTransferTracker1, massTransferTracker2:\t" << star1.flagRLOF  << " " << star2.flagRLOF << std::endl;
			}
		}
	}

	if(debugging){std::cout << "\n\n";}
}
// End Mass Transfer
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
// Begin Common Envelope
void BinaryStar::CommonEnvelopeEvent(const gsl_rng *r, programOptions const &options){
    /*
    The binary has entered a common envelope event. This function updates the binary parameters accordingly

    "Common-envelope evolution occurs either as a result of a collision between 
    a star with a dense core (k1 {2,3,4,5,6,8,9}) or at the onset of RLOF where mass 
    is transferred from a giant (k1 {2,3,4,5,6,8,9}) on a dynamical time-scale" Hurley et al. 2002, sec. 2.7.1.

    JIM BARRETT - 28/09/2016 - Need to use local value for hyperparameter sampling

    Parameters
    -----------
    options : programOptions
        User specified program options

    Returns
    ---------
    */
	Star	starCopy1=star1;
	Star	starCopy2=star2;
	
    double  alphaCE = m_commonEnvelopeAlpha;                            // alphaCE: CE efficiency parameter
	double  lambda1  = NEVER_SET;                             			// lambda1: parameter to measure "the central concentration of the donor"
	double  lambda2  = NEVER_SET;                             			// lambda2: parameter to measure "the central concentration of the donor"
	double	eta		= options.commonEnvelopeRecombinationEnergyDensity;	// eta: recombination energy density

    double  a       = m_SemiMajorAxisPrime;                             // Current semi-major axis in default units, AU
	double	e		= m_EccentricityPrime;								// Current eccentricity. 
	
    double  m1      = star1.m_Mass;                                     // Star 1 mass, in Msol units
    double  m2      = star2.m_Mass;                                     // Star 2 mass, in Msol units
    int 	type1   = star1.m_stellarType;                              // Star 1 stellar type
    int 	type2   = star2.m_stellarType;                              // Star 2 stellar type
	int 	finalType1 = NEVER_SET;
	int 	finalType2 = NEVER_SET;

    double  menv1   = m1 - star1.m_coreMass;                            // Star 1 envelope mass, in Msol units
    double  menv2   = m2 - star2.m_coreMass;                            // Star 2 envelope mass, in Msol units

    double mfin1 = 0.0;   // Star 1 mass after losing its envelope, in Msol units // NOTE: in this case, we asume it losses all of its envelope
    double mfin2 = 0.0;   // Star 2 mass after losing its envelope, in Msol units // NOTE: in this case, we asume it losses all of its envelope

    if(options.allowMainSequenceStarToSurviveCommonEnvelope){
        if(type1 == MS_LESS_THAN_07 or type1 == MS_MORE_THAN_07 or type1 == NAKED_HELIUM_STAR_MS){
            mfin1 = m1;
            menv1=0.0;
        }
        else{
            mfin1 = star1.m_coreMass;
        }

        if(type2 == MS_LESS_THAN_07 or type2 == MS_MORE_THAN_07 or type2 == NAKED_HELIUM_STAR_MS){
            mfin2 = m2;
            menv2=0.0;
        }
        else{
            mfin2 = star2.m_coreMass;
        }

    }
    else{
        mfin1   = star1.m_coreMass;                                 
        mfin2   = star2.m_coreMass;  
    }

    double  rRLd1   = a*rocheLobeRadius(m1, m2);                        // Roche lobe radius in AU at the moment where CEE begins, seen by star 1
    double  rRLd2   = a*rocheLobeRadius(m2, m1);                        // Roche lobe radius in AU at the moment where CEE begins, seen by star 2
    double  rRLdfin1= NEVER_SET;                        				// Roche lobe radius in AU after CEE, seen by star 1
    double  rRLdfin2= NEVER_SET;                        				// Roche lobe radius in AU after CEE, seen by star 2
	double	radius1 = star1.m_Radius;									// Radius of the primary at the onset of CEE, in Rsol
	double	radius2 = star2.m_Radius;									// Radius of the secondary at the onset of CEE, in Rsol
	double	radius1AfterStripping = radius1;							// Radius of the primary after (if) being stripped, in Rsol
	double	radius2AfterStripping = radius2;							// Radius of the secondary after (if) being stripped, in Rsol

    bool    envelopeFlag1 = false;
    bool    envelopeFlag2 = false;
    bool    debugging = false;

    // Alejandro - 10/05/2019 - Add flag for MS accretors for CEE catalogue paper. Can also be determined in post-processing.
    bool    primaryDonor = true;
    if(star1.flagRLOF == false and star2.flagRLOF == true ){primaryDonor = false;}

    if(primaryDonor == true){
        if (type2 == MS_LESS_THAN_07 or type2 == MS_MORE_THAN_07 or type2 == NAKED_HELIUM_STAR_MS){
            m_mainSequenceAccretorDuringCEEFlag = true;
        }
    }
    else{
        if (type1 == MS_LESS_THAN_07 or type1 == MS_MORE_THAN_07 or type1 == NAKED_HELIUM_STAR_MS){
            m_mainSequenceAccretorDuringCEEFlag = true;
        }
    }
        
    // Count how many times a CEE occurred 
    m_counterCEE++;

    // Reset double-core CEE flag
    m_doubleCoreCommonEnvelopeFlag = false;

    // ALEJANDRO - 29/01/2019 - Check for simultaneous RLOF
    if(star1.flagRLOF == true and star2.flagRLOF == true){
        m_simultaneousRLOFleadingToCEEFlag = true;
    }
    else{
        m_simultaneousRLOFleadingToCEEFlag = false;       
    }

  //   if(menv1 < 0.0){
		// std::cerr << m_randomSeed <<  "\tError in CE function. Mass of envelope of the primary has a negative value.\n";
  //   }
  //   if(menv2 < 0.0){
		// std::cerr << m_randomSeed <<  "\tError in CE function. Mass of envelope of the secondary has a negative value.\n";
  //   }


	if(debugging){
		std::cout << "Inside CEE function" << std::endl;	
		std::cout << "Printing values at the begginging of Common Envelope function." << std::endl;
		std::cout << "alpha:\t" << alphaCE << std::endl;
		std::cout << "eta:\t" << eta << std::endl;
		printingBinaryVariables();
  
		std::cout << "ORBIT" << std::endl;
		std::cout << "Semi-major axis [AU]:\t" << a << std::endl;
		std::cout << "Eccentricity:\t" << m_Eccentricity << std::endl;
		std::cout << "EccentricityPrime:\t" << e << std::endl;
		std::cout << "rRLd1 [AU]: " << rRLd1 << std::endl;
		std::cout << "rRLd1 [Rsol]: " << rRLd1*AUToRsol << std::endl;
		std::cout << "rRLd1/a: " << rRLd1/a << std::endl;
		std::cout << "rRLd2 [AU]: " << rRLd2 << std::endl;
		std::cout << "rRLd2 [Rsol]: " << rRLd2*AUToRsol << std::endl;
		std::cout << "rRLd2/a: " << rRLd2/a << std::endl;
    }

	bool	useWinds = false;

	// Lambdas
    // SIMON: Do we really need to calculate all of these? 
	// Alejandro: Not really, but I like to have them handy for the common envelopes file. It is only done once or twice per run, so I think it should be fine.
	double 	Nanjing1	= lambdaNanjing(star1.m_MZAMS, star1.m_Mass, star1.m_coreMass, star1.m_Radius, star1.m_stellarType, star1.m_Metallicity, options);
	double	Nanjing2	= lambdaNanjing(star2.m_MZAMS, star2.m_Mass, star2.m_coreMass, star2.m_Radius, star2.m_stellarType, star2.m_Metallicity, options);
	double	Loveridge1	= lambdaLoveridgeEnergyFormalism(star1.m_Metallicity, star1.m_MZAMS, m1, menv1, star1.m_COCoreMass, radius1, false, m_randomSeed);
	double	Loveridge2	= lambdaLoveridgeEnergyFormalism(star2.m_Metallicity, star2.m_MZAMS, m2, menv2, star2.m_COCoreMass, radius2, false, m_randomSeed);
	double	Loveridge1Winds	= lambdaLoveridgeEnergyFormalism(star1.m_Metallicity, star1.m_MZAMS, m1, menv1, star1.m_COCoreMass, radius1, true, m_randomSeed);
	double	Loveridge2Winds	= lambdaLoveridgeEnergyFormalism(star2.m_Metallicity, star2.m_MZAMS, m2, menv2, star2.m_COCoreMass, radius2, true, m_randomSeed);
	double	Kruckow1		= lambdaKruckow(star1.m_Radius, options.commonEnvelopeSlopeKruckow);
	double	Kruckow2		= lambdaKruckow(star2.m_Radius, options.commonEnvelopeSlopeKruckow);
	double	Dewi1		= lambdaDewi(star1.m_Mass, star1.m_coreMass, star1.m_Radius, star1.m_RZAMS, star1.m_Luminosity, star1.m_stellarType, m_randomSeed);
	double	Dewi2		= lambdaDewi(star2.m_Mass, star2.m_coreMass, star2.m_Radius, star2.m_RZAMS, star2.m_Luminosity, star2.m_stellarType, m_randomSeed);

	// Binding energies
	// m_bindingEnergyAtCommonEnvelope
	double	bindingEnergyFixed1 = bindingEnergy(mfin1, menv1, radius1, options.commonEnvelopeLambda, m_randomSeed);
	double	bindingEnergyFixed2 = bindingEnergy(mfin2, menv2, radius2, options.commonEnvelopeLambda, m_randomSeed);
	double	bindingEnergyNanjing1 = bindingEnergy(mfin1, menv1, radius1, Nanjing1, m_randomSeed);
	double	bindingEnergyNanjing2 = bindingEnergy(mfin2, menv2, radius2, Nanjing2, m_randomSeed);
	double	bindingEnergyLoveridge1 = bindingEnergy(mfin1, menv1, radius1, Loveridge1, m_randomSeed);
	double	bindingEnergyLoveridge2 = bindingEnergy(mfin2, menv2, radius2, Loveridge2, m_randomSeed);
	double	bindingEnergyLoveridge1Winds = bindingEnergy(mfin1, menv1, radius1, Loveridge1Winds, m_randomSeed);
	double	bindingEnergyLoveridge2Winds = bindingEnergy(mfin2, menv2, radius2, Loveridge2Winds, m_randomSeed);
	double	bindingEnergyKruckow1 =	 bindingEnergy(mfin1, menv1, radius1, Kruckow1, m_randomSeed);
	double	bindingEnergyKruckow2 =	 bindingEnergy(mfin2, menv2, radius2, Kruckow2, m_randomSeed);
	

	if(debugging){std::cout << "LAMBDAS" << std::endl;}

	if(options.commonEnvelopeLambdaPrescriptionString == "LAMBDA_FIXED"){
		lambda1  = options.commonEnvelopeLambda;
		lambda2  = options.commonEnvelopeLambda;
		star1.m_bindingEnergyAtCommonEnvelope = bindingEnergyFixed1;
		star2.m_bindingEnergyAtCommonEnvelope = bindingEnergyFixed2;
	}
	else if(options.commonEnvelopeLambdaPrescriptionString == "LAMBDA_LOVERIDGE"){
		lambda1 = Loveridge1;
		lambda2 = Loveridge2;
		star1.m_bindingEnergyAtCommonEnvelope = bindingEnergyLoveridge1;
		star2.m_bindingEnergyAtCommonEnvelope = bindingEnergyLoveridge2;
    }
	else if(options.commonEnvelopeLambdaPrescriptionString == "LAMBDA_NANJING"){
		lambda1 = Nanjing1;
		lambda2 = Nanjing2;
		star1.m_bindingEnergyAtCommonEnvelope = bindingEnergyNanjing1;
		star2.m_bindingEnergyAtCommonEnvelope = bindingEnergyNanjing2;
	}
	else if(options.commonEnvelopeLambdaPrescriptionString == "LAMBDA_KRUCKOW"){
		lambda1 = Kruckow1;
		lambda2 = Kruckow2;
		star1.m_bindingEnergyAtCommonEnvelope = bindingEnergyKruckow1;
		star2.m_bindingEnergyAtCommonEnvelope = bindingEnergyKruckow2;
	}
	else{
		std::cerr << m_randomSeed << "\tNo common envelope lambda prescription chosen. Shouldn't get here." << std::endl;
	}

	if(lambda1 < 0.00001){
		lambda1 = 0.0;
	}
	
	if(lambda2 < 0.00001){
		lambda2 = 0.0;
	}

    if(debugging){
        std::cout << "lambda 1,2: " << lambda1 << "\t" << lambda2 << std::endl;
    }

    // Multiply by constant (default = 1)
    lambda1 *= options.commonEnvelopeLambdaMultiplier;
    lambda2 *= options.commonEnvelopeLambdaMultiplier;

    // Record parameters
    star1.m_HeCoreMassAtCommonEnvelope = star1.m_HeCoreMass;
    star1.m_COCoreMassAtCommonEnvelope = star1.m_COCoreMass;
    star1.m_coreMassAtCommonEnvelope = star1.m_coreMass;

    star2.m_HeCoreMassAtCommonEnvelope = star2.m_HeCoreMass;
    star2.m_COCoreMassAtCommonEnvelope = star2.m_COCoreMass;
    star2.m_coreMassAtCommonEnvelope = star2.m_coreMass;

    if((menv1 > 0.0)and(mfin1 > 0.0)){
		envelopeFlag1 = true;
		if(debugging){
			std::cout << "Loveridge1:\t" << Loveridge1 << std::endl;
			std::cout << "Loveridge1W:\t" << Loveridge1Winds << std::endl;
			std::cout << "Nanjing1:\t" << Nanjing1 << std::endl;
			std::cout << "Kruckow1Bot:\t" << lambdaKruckow(star1.m_Radius, -1.0)<< std::endl;
			std::cout << "Kruckow1Mid:\t" << lambdaKruckow(star1.m_Radius, -4.0/5.0)<< std::endl;
			std::cout << "Kruckow1Top:\t" << lambdaKruckow(star1.m_Radius, -2.0/3.0)<< std::endl;
			std::cout << "Dewi1:\t\t" << Dewi1 << std::endl;
			std::cout << "rRLd1:\t" << rRLd1*AUToRsol << std::endl;
			std::cout << "RLd:\t" << rRLd1/a << std::endl;
			std::cout << "bEFixed1 [ergs]:\t" << bindingEnergyFixed1 << std::endl;			
			std::cout << "bENanjing1 [ergs]:\t" << bindingEnergyNanjing1 << std::endl;
			std::cout << "bELoveridge1 [ergs]:\t" << bindingEnergyLoveridge1 << std::endl;
			std::cout << "bELoveridge1Winds [ergs]:\t" << bindingEnergyLoveridge1Winds << std::endl;
			std::cout << "bEKruckow1 [ergs]:\t" << bindingEnergyKruckow1 << std::endl;
		}
	}
    if((menv2 > 0.0)and(mfin2 > 0.0)){
		envelopeFlag2 = true;
		if(debugging){
			std::cout << "Loveridge2:\t" << Loveridge2 << std::endl;
			std::cout << "Loveridge2W:\t" << Loveridge2Winds << std::endl;
			std::cout << "Nanjing2:\t" << Nanjing2 << std::endl;
			std::cout << "Kruckow2Bot:\t" << lambdaKruckow(star2.m_Radius, -1.0)<< std::endl;
			std::cout << "Kruckow2Mid:\t" << lambdaKruckow(star2.m_Radius, -4.0/5.0)<< std::endl;
			std::cout << "Kruckow2Top:\t" << lambdaKruckow(star2.m_Radius, -2.0/3.0)<< std::endl;
			std::cout << "Dewi2:\t\t" << Dewi1 << std::endl;
			std::cout << "rRLd2:\t" << rRLd2*AUToRsol << std::endl;
			std::cout << "RLd:\t" << rRLd2/a << std::endl;
			std::cout << "bEFixed2 [ergs]:\t" << bindingEnergyFixed2 << std::endl;			
			std::cout << "bENanjing2 [ergs]:\t" << bindingEnergyNanjing2 << std::endl;
			std::cout << "bELoveridge2 [ergs]:\t" << bindingEnergyLoveridge2 << std::endl;
			std::cout << "bELoveridge2Winds [ergs]:\t" << bindingEnergyLoveridge2Winds << std::endl;
			std::cout << "bEKruckow2 [ergs]:\t" << bindingEnergyKruckow2 << std::endl;
		}
	}

    double  afinal 		= NEVER_SET;           	// Initialise variable
    double  omegafinal 	= NEVER_SET;       		// Initialise variable

    // Check if we have an HG star, and if are allowing such a star to survive the CE

    if(type1 == HERTZSPRUNG_GAP or type2 == HERTZSPRUNG_GAP or type1 == NAKED_HELIUM_STAR_HETZSPRUNG_GAP or type2 == NAKED_HELIUM_STAR_HETZSPRUNG_GAP){
            if(debugging){std::cout << "Hertzsprung Gap star in a CE. Making optimistic assumption and applying energy balance" << std::endl;}
			// Set the optimistic CE flag to true
			m_optimisticCommonEnvelopeFlag = true;
	}
	
	// CommonEnvelopeEvent
    afinal = doubleCommonEnvelopePhase(alphaCE, lambda1, lambda2, a, m1, mfin1, menv1, rRLd1, m2, mfin2, menv2, rRLd2, type1, type2, m_doubleCoreCommonEnvelopeFlag);   // Double CE prescription becomes single CE if menv2 = 0.0. Therefore, use as default.

    omegafinal = sqrt((mfin1+mfin2)*G1*pow(afinal,-3.0)); // SIMON: Should be a function for this equation

	rRLdfin1 = afinal*rocheLobeRadius(mfin1, mfin2);
	rRLdfin2 = afinal*rocheLobeRadius(mfin2, mfin1);

    // Correct for stellar types, stellar mass, separation and period
    m_SemiMajorAxisPrime = afinal;
    m_orbitalVelocityPrime = omegafinal;

    // We assume that a common envelope event (CEE) circularises the binary
    m_Eccentricity = 0.0;
    m_EccentricityPrime = 0.0;

    // Update star mass
    if(star1.m_stellarType != NEUTRON_STAR){
        star1.m_Mass = mfin1;
	//std::cout << "star1.m_stellarType " << star1.m_stellarType << std::endl;
        //std::cout << "star1.m_Mass_notNS  " << star1.m_Mass << std::endl;
    }

    if(star2.m_stellarType != NEUTRON_STAR){
        star2.m_Mass = mfin2;
	//std::cout << "star2.m_stellarType " << star2.m_stellarType << std::endl;
        //std::cout << "star2.m_Mass_notNS  " << star2.m_Mass << std::endl;
    }

    // Accretion during common envelope onto neutron star
    if(star1.m_stellarType == NEUTRON_STAR){
	//std::cout << "star1.m_stellarType " << star1.m_stellarType << std::endl;
	double accretedNSMass1 = massAccretionNeutronStar(r, options, m2, radius2) ;
	//std::cout << "star1 mass NS before accretion = " << star1.m_Mass << std::endl;
        star1.m_Mass += accretedNSMass1 ;
	star1.m_MassTransferDiff = accretedNSMass1;
        //std::cout << "star1.m_Mass_NS after accretion " << star1.m_Mass << std::endl;
    }

    if(star2.m_stellarType == NEUTRON_STAR){
	//std::cout << "star2.m_stellarType " << star2.m_stellarType << std::endl;
	double accretedNSMass2 = massAccretionNeutronStar(r, options, m1, radius1) ;
	
        star2.m_Mass += accretedNSMass2;
	star2.m_MassTransferDiff = accretedNSMass2;
        //std::cout << "star2.m_Mass_NS" << star2.m_Mass << std::endl;
    }
    
    if(debugging){
		std::cout << "a\t[AU]=\t" << a << std::endl;
		std::cout << "afinal\t[AU]=\t" << afinal << std::endl;
		std::cout << "a\t[Rsol]=\t" << a*AUToRsol << std::endl;
		std::cout << "afinal\t[Rsol]=\t" << afinal*AUToRsol << std::endl;
		std::cout << "rRLd1\t[Rsol]=\t" << rRLd1*AUToRsol << std::endl;
        std::cout << "radius1\t[Rsol]=\t" << radius1 << std::endl;
		std::cout << "rRLdfin1[Rsol]=\t" << rRLdfin1*AUToRsol << std::endl;
        std::cout << "radius1AfterStripping[Rsol]=\t" << radius1AfterStripping << std::endl;
		std::cout << "rRLd2\t[Rsol]=\t" << rRLd2*AUToRsol << std::endl;
        std::cout << "radius2\t[Rsol]=\t" << radius2 << std::endl;
		std::cout << "rRLdfin2[Rsol]=\t" << rRLdfin2*AUToRsol << std::endl;
        std::cout << "radius2AfterStripping[Rsol]=\t" << radius2AfterStripping << std::endl;
		std::cout << "Eccentricity=\t" << e << std::endl;
		std::cout << "EccentricityFin=\t" << m_Eccentricity << std::endl;
	}

	// Alejandro - 18/02018 - Calculate tidal timescales for primary donor. Only change if secondary is donor (see below).
	double periastronRsol	= a*AUToRsol*(1-e);
    double separationRsol   = a*AUToRsol;
    double chosenSeparationTides = separationRsol;
	m_synchronizationTimescale	= calculateSyncrhonizationTimescale(starCopy1, m_orbitalVelocityPrime, chosenSeparationTides, m2);
	m_circularizationTimescale 	= calculateCircularisationTimescale(starCopy1, m_orbitalVelocityPrime, chosenSeparationTides, m2);

    // update stellar type after losing it's envelope. Star1, Star2 or both if double CEE.
    if((envelopeFlag1 == true)&&(envelopeFlag2 == false)){
        star1.modifyStarAfterLosingEnvelope(star1.m_stellarType, star1.m_Mass);
		radius1AfterStripping = star1.m_Radius;
		radius2AfterStripping = star2.m_Radius;
		m_massTransferTrackerHistory = CE_FROM_1_TO_2;
    }
    else if((envelopeFlag1 == false)&&(envelopeFlag2 == true)){
        star2.modifyStarAfterLosingEnvelope(star2.m_stellarType, star2.m_Mass);
		radius1AfterStripping = star1.m_Radius;
		radius2AfterStripping = star2.m_Radius;
		m_massTransferTrackerHistory = CE_FROM_2_TO_1;
		m_synchronizationTimescale	= calculateSyncrhonizationTimescale(starCopy2, m_orbitalVelocityPrime, chosenSeparationTides, m1);
		m_circularizationTimescale 	= calculateCircularisationTimescale(starCopy2, m_orbitalVelocityPrime, chosenSeparationTides, m1);
    }
    else if((envelopeFlag1 == true)&&(envelopeFlag2 == true)){
        // Case of double CEE
        star1.modifyStarAfterLosingEnvelope(star1.m_stellarType, star1.m_Mass);
        star2.modifyStarAfterLosingEnvelope(star2.m_stellarType, star2.m_Mass);
		radius1AfterStripping = star1.m_Radius;
		radius2AfterStripping = star2.m_Radius;
		m_massTransferTrackerHistory = CE_DOUBLE_CORE;
        m_doubleCoreCommonEnvelopeFlag = true;
    }
    else{
        if(debugging){
            std::cout << "menv1: " << menv1 << " and menv2: " << menv2 << ". No clear core-envelope separation.\n";
            std::cout << "type1: " << type1 << " and type2: " << type2 << ". \n";
            std::cout << "No clear core-envelope separation.\n";
        }


        if(type1<=NAKED_HELIUM_STAR_MS and type2<=NAKED_HELIUM_STAR_MS){
            if(debugging){
                std::cout << "Wet merger." << std::endl;
            }
            m_massTransferTrackerHistory = CE_BOTH_MS;            
        }
        else{
            if(debugging){
                std::cout << "CEE of a MS donor with a CO accretor." << std::endl;
            }
            m_massTransferTrackerHistory = CE_MS_WITH_CO;            
        }

        // Here MS-WD systems are flagged as CE_BOTH_MS
        m_stellarMerger = true;
    }
	
	// Assign final stellar type after CEE
	finalType1 = star1.m_stellarType;
	finalType2 = star2.m_stellarType;

    
    if(afinal <= 0.0 or ( (radius1AfterStripping+radius2AfterStripping) > afinal*AUToRsol)){
        if(debugging){
            std::cout << "R1+R2 [Rsol]:\t" << star1.m_Radius+star2.m_Radius << std::endl;
            std::cout << "a [Rsol]:\t" << a << std::endl;
			std::cout << "afinal [Rsol]:\t" << m_SemiMajorAxisPrime << std::endl;
            std::cout << "a [AU]:\t" << a << std::endl;
			std::cout << "afinal [AU]:\t" << m_SemiMajorAxisPrime << std::endl;
            printingBinaryVariables();

	       if(afinal <= 0.0)
		      std::cout << "Inside CEE. afinal <= 0." << std::endl;
	       else
		      std::cout << "Inside CEE. R1+R2 > afinal." << std::endl;

	   }

	m_stellarMerger = true;
    }

	// ALEJANDRO - 10/11/2016 - Record the value of lambda used in the CEE.  Double check how to deal with this if there are several CEE.
	star1.m_lambdaAtCommonEnvelope		= lambda1;
	star2.m_lambdaAtCommonEnvelope		= lambda2;
	m_commonEnvelopeOccuredAtLeastOnce 	= true;

	// ALEJANDRO - 06/12/2016 - Record values pre-post CE, for populations studies. All separations in Rsol.
	m_EccentricityPreCEE = e;
	m_EccentricityPostCEE = m_Eccentricity;
	m_SemiMajorAxisPreCEE = a*AUToRsol;
	m_SemiMajorAxisPostCEE = afinal*AUToRsol;
	m_rocheLobe1to2PreCEE = rRLd1*AUToRsol;
	m_rocheLobe1to2PostCEE = rRLdfin1*AUToRsol;
	m_rocheLobe2to1PreCEE = rRLd2*AUToRsol;
	m_rocheLobe2to1PostCEE = rRLdfin2*AUToRsol;

    // ALEJANDRO - 28/01/2019 - Check for RLOF immediatedly after the CEE. A check for it during the next timestep is done in evaluateBinary funtion.
    if((radius1AfterStripping >= m_rocheLobe1to2PostCEE)or(radius2AfterStripping >= m_rocheLobe2to1PostCEE)){
        m_immediateRLOFAfterCEE = true;
    }

	
	if(debugging){
		std::cout << "Printing values at the end of Common Envelope function." << std::endl;
		printingBinaryVariables();
	}


    // SIMON: Ideally this should be just printCommonEnvelope, a function defined in printing.cpp
    // However, I can't figure out how to pass this binary to such a function, so I bail on this for now
    // Headers printed by printCommonEnvelopeHeader in printing.cpp
    std::ofstream commonEnvelope((options.outputPath/"commonEnvelopes.txt").string(),std::ios_base::app);
	commonEnvelope 	<< m_randomSeed << TAB 
					<< m_time << TAB 
					<< alphaCE << TAB 
					<< lambda1 << TAB 
					<< lambda2 << TAB 
					<< star1.m_bindingEnergyAtCommonEnvelope << TAB 
					<< star2.m_bindingEnergyAtCommonEnvelope << TAB 
                    << m_EccentricityRLOF << TAB 
					<< m_EccentricityPreCEE << TAB 
					<< m_EccentricityPostCEE << TAB 
					<< m_SemiMajorAxisPreCEE << TAB 
					<< m_SemiMajorAxisPostCEE << TAB 
					<< m_rocheLobe1to2PreCEE << TAB 
					<< m_rocheLobe1to2PostCEE << TAB 
					<< m_rocheLobe2to1PreCEE << TAB 
					<< m_rocheLobe2to1PostCEE << TAB 
					<< star1.m_MZAMS << TAB 
					<< m1 << TAB 
					<< menv1 << TAB 
					<< mfin1 << TAB 
					<< radius1 << TAB 
					<< radius1AfterStripping << TAB 
					<< type1 << TAB 
					<< finalType1 << TAB 
					<< options.commonEnvelopeLambda << TAB << Nanjing1 << TAB 
					<< Loveridge1 << TAB 
					<< Loveridge1Winds << TAB 
					<< Kruckow1 << TAB 
					<< bindingEnergyFixed1 << TAB 
					<< bindingEnergyNanjing1 << TAB 
					<< bindingEnergyLoveridge1 << TAB 
					<< bindingEnergyLoveridge1Winds << TAB 
					<< bindingEnergyKruckow1 << TAB 
					<< star2.m_MZAMS  << TAB 
					<< m2 << TAB 
					<< menv2 << TAB 
					<< mfin2 << TAB 
					<< radius2 << TAB 
					<< radius2AfterStripping << TAB 
					<< type2 << TAB 
					<< finalType2 << TAB 
					<< options.commonEnvelopeLambda << TAB 
					<< Nanjing2 << TAB 
					<< Loveridge2 << TAB 
					<< Loveridge2Winds << TAB 
					<< Kruckow2 << TAB 
					<< bindingEnergyFixed2 << TAB 
					<< bindingEnergyNanjing2 << TAB 
					<< bindingEnergyLoveridge2 << TAB 
					<< bindingEnergyLoveridge2Winds << TAB 
					<< bindingEnergyKruckow2 << TAB 
					<< m_massTransferTrackerHistory << TAB 
					<< m_stellarMerger << TAB 
					<< m_optimisticCommonEnvelopeFlag << TAB 
					<< m_counterCEE << TAB
					<< m_doubleCoreCommonEnvelopeFlag << TAB					
					<< star1.flagRLOF << TAB
					<< star1.m_Luminosity << TAB
					<< star1.m_Temperature << TAB
					<< star1.m_dynamicalTimescale << TAB
					<< star1.m_thermalTimescale << TAB
					<< star1.m_nuclearTimescale << TAB
					<< star2.flagRLOF << TAB
					<< star2.m_Luminosity << TAB 
					<< star2.m_Temperature << TAB
					<< star2.m_dynamicalTimescale << TAB
					<< star2.m_thermalTimescale << TAB
					<< star2.m_nuclearTimescale << TAB
					<< m_zetaStarCompare << TAB
					<< m_zetaRLOFanalytic << TAB
					<< m_synchronizationTimescale << TAB
					<< m_circularizationTimescale << TAB
					<< star1.m_radialExpansionTimescale << TAB
					<< star2.m_radialExpansionTimescale << TAB 
                    << m_immediateRLOFAfterCEE << TAB
                    << m_simultaneousRLOFleadingToCEEFlag << TAB
                    << m_mainSequenceAccretorDuringCEEFlag
					<< std::endl;

    commonEnvelope.close();
}
// End Common Envelope
//-------------------------------------------------------------------------------------------------------------------------------------//
//-------------------------------------------------------------------------------------------------------------------------------------//
// Testing Loveridge Binding energy fits
void    BinaryStar::bindingEnergyFits(const programOptions &options){
    // not used by default

    std::cout << "Z0\tZ\tMzams(Mo)\tM(Mo)\tR(Ro)\tGB\tLogBE\n";
    double  CO = 0.0;
    double  Radius = 100.0;
    double  MassZams = 1.0;
    double  Mass = 1.0;

    for(double Z=0.0001; Z<=0.03; Z+=0.002){
        calculateLogBindingEnergyLoveridge(Z, 1.5, 1.00, 0.0, 100.0, false, m_randomSeed);
    }
	std::cout <<  std::endl;

	for(double Z=0.0001; Z<=0.03; Z+=0.002){
        calculateLogBindingEnergyLoveridge(Z, 1.5, 1.00, 0.0, 100.0, true, m_randomSeed);
    }

	std::cout << std::endl;


    for(double Z=0.0001; Z<=0.03; Z+=0.002){
        calculateLogBindingEnergyLoveridge(Z, 25.0, 20.0, 0.0, 100.0, false, m_randomSeed);
    }
	std::cout <<  std::endl;

	for(double Z=0.0001; Z<=0.03; Z+=0.002){
        calculateLogBindingEnergyLoveridge(Z, 25.0, 20.0, 0.0, 100.0, true, m_randomSeed);
    }

	std::cout << std::endl;
/*
     1:  Star of 1Mo, Z=0.02, no wind mass loss and different radii (RGB+AGB):
     #           Z  Mzams (Mo)      M (Mo)      R (Ro)      gb  branch        log(BE/erg)      BE (erg)      log(BE_r/erg)    BE_r (erg)
     1      0.0200        1.00        1.00       30.00       1     RGB             47.153      1.42E+47             46.311      2.05E+46
     2      0.0200        1.00        1.00       60.00       1     RGB             46.942      8.75E+46             46.271      1.86E+46
     3      0.0200        1.00        1.00       90.00       1     RGB             46.825      6.68E+46             46.226      1.68E+46
     4      0.0200        1.00        1.00      120.00       1     RGB             46.739      5.48E+46             46.180      1.51E+46
     5      0.0200        1.00        1.00      150.00       1     RGB             46.666      4.64E+46             46.136      1.37E+46
     6      0.0200        1.00        1.00      180.00       2     AGB             46.769      5.88E+46             46.093      1.24E+46
     7      0.0200        1.00        1.00      210.00       2     AGB             46.690      4.90E+46             46.052      1.13E+46
     8      0.0200        1.00        1.00      240.00       2     AGB             46.615      4.12E+46             46.012      1.03E+46
     9      0.0200        1.00        1.00      270.00       2     AGB             46.545      3.51E+46             45.974      9.41E+45
     10      0.0200        1.00        1.00      300.00       2     AGB             46.477      3.00E+46             45.937      8.64E+45
     */

    /*
     2:  Stars of different masses, for Z=0.02, R=25Ro (RGB):
     #           Z  Mzams (Mo)      M (Mo)      R (Ro)      gb  branch        log(BE/erg)      BE (erg)      log(BE_r/erg)    BE_r (erg)
     1      0.0200        1.11        1.00       25.00       1     RGB             47.246      1.76E+47             46.317      2.07E+46
     2      0.0200        2.22        2.00       25.00       1     RGB             47.701      5.02E+47             46.718      5.22E+46
     3      0.0200        3.33        3.00       25.00       1     RGB             48.184      1.53E+48             46.925      8.42E+46
     4      0.0200        4.44        4.00       25.00       1     RGB             48.675      4.73E+48             47.062      1.15E+47
     5      0.0200        5.56        5.00       25.00       1     RGB             49.062      1.15E+49             47.162      1.45E+47
     6      0.0200        6.67        6.00       25.00       1     RGB             49.189      1.54E+49             47.241      1.74E+47
     7      0.0200        7.78        7.00       25.00       1     RGB             49.295      1.97E+49             47.306      2.02E+47
     8      0.0200        8.89        8.00       25.00       1     RGB             49.388      2.44E+49             47.361      2.29E+47
     9      0.0200       10.00        9.00       25.00       1     RGB             49.470      2.95E+49             47.408      2.56E+47
     10      0.0200       11.11       10.00       25.00       1     RGB             49.543      3.49E+49             47.449      2.81E+47
     */

    /*
     3:  Stars of 1Mo and different metallicities (code chooses nearest Z), for R=100Ro (RGB):
     #           Z  Mzams (Mo)      M (Mo)      R (Ro)      gb  branch        log(BE/erg)      BE (erg)      log(BE_r/erg)    BE_r (erg)
     1      0.0001        1.00        1.00      100.00       1     RGB             46.865      7.32E+46             46.063      1.16E+46
     2      0.0021        1.00        1.00      100.00       1     RGB             46.823      6.66E+46             46.134      1.36E+46
     3      0.0041        1.00        1.00      100.00       1     RGB             46.810      6.46E+46             46.186      1.53E+46
     4      0.0061        1.00        1.00      100.00       1     RGB             46.810      6.46E+46             46.186      1.53E+46
     5      0.0081        1.00        1.00      100.00       1     RGB             46.810      6.46E+46             46.186      1.53E+46
     6      0.0101        1.00        1.00      100.00       1     RGB             46.810      6.46E+46             46.186      1.53E+46
     7      0.0121        1.00        1.00      100.00       1     RGB             46.810      6.46E+46             46.186      1.53E+46
     8      0.0141        1.00        1.00      100.00       1     RGB             46.800      6.31E+46             46.183      1.52E+46
     9      0.0161        1.00        1.00      100.00       1     RGB             46.800      6.31E+46             46.183      1.52E+46
     10      0.0181        1.00        1.00      100.00       1     RGB             46.794      6.22E+46             46.210      1.62E+46
     11      0.0201        1.00        1.00      100.00       1     RGB             46.794      6.22E+46             46.210      1.62E+46
     12      0.0221        1.00        1.00      100.00       1     RGB             46.794      6.22E+46             46.210      1.62E+46
     13      0.0241        1.00        1.00      100.00       1     RGB             46.794      6.22E+46             46.210      1.62E+46
     14      0.0261        1.00        1.00      100.00       1     RGB             46.790      6.17E+46             46.215      1.64E+46
     15      0.0281        1.00        1.00      100.00       1     RGB             46.790      6.17E+46             46.215      1.64E+46
     */

}
// End testing Loveridge Binding energy fits
//------------------------------------------------------------------------------------------------------------------------------------//
//-------------------------------------------------------------------------------------------------------------------------------------//
// Testing Loveridge Binding energy fits
void    BinaryStar::comparingLambdas(const programOptions &options){
    // not used by default

    double  a       = m_SemiMajorAxisPrime;                             // Current semi-major axis in default units, AU
    double  m1      = star1.m_Mass;                                     // Star 1 mass, in Msol units
    double  m2      = star2.m_Mass;                                     // Star 2 mass, in Msol units
    int type1       = star1.m_stellarType;                              // Star 1 stellar type
    int type2       = star2.m_stellarType;                              // Star 2 stellar type

    double  menv1   = m1 - star1.m_coreMass;                            // Star 1 envelope mass, in Msol units
    double  menv2   = m2 - star2.m_coreMass;                            // Star 2 envelope mass, in Msol units
    double  mfin1   = star1.m_coreMass;                                 // Star 1 mass after losing its envelope, in Msol units // NOTE: in this case, we asume it losses all of its envelope
    double  mfin2   = star2.m_coreMass;                                 // Star 2 mass after losing its envelope, in Msol units // NOTE: in this case, we asume it losses all of its envelope
    double  rRLd1   = a*rocheLobeRadius(m1, m2);                        // Roche lobe radius in AU at the moment where CEE begins, seen by star 1
    double  rRLd2   = a*rocheLobeRadius(m2, m1);                        // Roche lobe radius in AU at the moment where CEE begins, seen by star 2

	bool	useWinds = false;

	double	Loveridge1 = NEVER_SET;
	double	Loveridge2 = NEVER_SET;
	double	Nanjing1 = NEVER_SET;
	double	Nanjing2 = NEVER_SET;

  // JIM BARRETT - 28/11/2016 - Passing through the random seed for error tracking purposes
	Loveridge1	= lambdaLoveridgeEnergyFormalism(star1.m_Metallicity, star1.m_MZAMS, m1, menv1, star1.m_COCoreMass, star1.m_Radius,useWinds, m_randomSeed);
	Loveridge2	= lambdaLoveridgeEnergyFormalism(star2.m_Metallicity, star2.m_MZAMS, m2, menv2, star2.m_COCoreMass, star2.m_Radius,useWinds, m_randomSeed);
	Nanjing1	= lambdaNanjing(star1.m_MZAMS, star1.m_Mass, star1.m_coreMass, star1.m_Radius, star1.m_stellarType, star1.m_Metallicity, options);
	Nanjing2	= lambdaNanjing(star2.m_MZAMS, star2.m_Mass, star2.m_coreMass, star2.m_Radius, star2.m_stellarType, star2.m_Metallicity, options);

	std::cout << "Lambdas:" << std::endl;
	std::cout << "Loveridge1\tNanjing1\tLoveridge2\tNanjing2" << std::endl;
	std::cout << Loveridge1 << "\t" << Nanjing1 << "\t\t" << Loveridge2 << "\t" << Nanjing2 << std::endl << std::endl;
}
//------------------------------------------------------------------------------------------------------------------------------------//
void    BinaryStar::calculatingLambdas(const programOptions &options){
// ALEJANDRO - 04/11/2016 - Lambda calculations as tracker for binding energy; lambda is recalculated in the CommonEnvelopeEvent function. Maybe we can give it as an argument later, as it is already being calculated here.
	star1.m_NanjingLambda			= lambdaNanjing(star1.m_MZAMS, star1.m_Mass, star1.m_coreMass, star1.m_Radius, star1.m_stellarType, star1.m_Metallicity, options);
	star1.m_LoveridgeLambda			= lambdaLoveridgeEnergyFormalism(star1.m_Metallicity, star1.m_MZAMS, star1.m_Mass, star1.m_Mass-star1.m_coreMass, star1.m_COCoreMass, star1.m_Radius, false, m_randomSeed);
	star1.m_LoveridgeWindsLambda	= lambdaLoveridgeEnergyFormalism(star1.m_Metallicity, star1.m_MZAMS, star1.m_Mass, star1.m_Mass-star1.m_coreMass, star1.m_COCoreMass, star1.m_Radius, true, m_randomSeed);
	star1.m_KruckowTopLambda		= lambdaKruckow(star1.m_Radius, -2.0/3.0);
	star1.m_KruckowMidLambda		= lambdaKruckow(star1.m_Radius, -4.0/5.0);
	star1.m_KruckowBotLambda		= lambdaKruckow(star1.m_Radius, -1.0);
	star1.m_DewiLambda				= lambdaDewi(star1.m_Mass, star1.m_coreMass, star1.m_Radius, star1.m_RZAMS, star1.m_Luminosity, star1.m_stellarType, m_randomSeed);

	star2.m_NanjingLambda			= lambdaNanjing(star2.m_MZAMS, star2.m_Mass, star2.m_coreMass, star2.m_Radius, star2.m_stellarType, star2.m_Metallicity, options);
	star2.m_LoveridgeLambda			= lambdaLoveridgeEnergyFormalism(star2.m_Metallicity, star2.m_MZAMS, star2.m_Mass, star2.m_Mass-star2.m_coreMass, star2.m_COCoreMass, star2.m_Radius, false, m_randomSeed);
	star2.m_LoveridgeWindsLambda	= lambdaLoveridgeEnergyFormalism(star2.m_Metallicity, star2.m_MZAMS, star2.m_Mass, star2.m_Mass-star2.m_coreMass, star2.m_COCoreMass, star2.m_Radius, true, m_randomSeed);
	star2.m_KruckowTopLambda		= lambdaKruckow(star2.m_Radius, -2.0/3.0);
	star2.m_KruckowMidLambda		= lambdaKruckow(star2.m_Radius, -4.0/5.0);
	star2.m_KruckowBotLambda		= lambdaKruckow(star2.m_Radius, -1.0);
	star2.m_DewiLambda				= lambdaDewi(star2.m_Mass, star2.m_coreMass, star2.m_Radius, star2.m_RZAMS, star2.m_Luminosity, star2.m_stellarType, m_randomSeed);
}
//------------------------------------------------------------------------------------------------------------------------------------//
void    BinaryStar::calculatingZetas(const programOptions &options, const gsl_rng *r){
// ALEJANDRO - 16/10/2017 - Zeta calculations as a mass transfer tracker.
	bool	debugging = false;
	Star   starCopy1 = star1;
	Star   starCopy2 = star2;
	
	star1.m_zetaThermal 		= determineConvergedMassStepZetaThermal(starCopy1, options, star1.m_calculateZetaThermalErrorMassFlag, star1.m_calculateZetaThermalErrorRadiusFlag, r);
	star1.m_zetaNuclear			= determineConvergedTimestepZetaNuclear(starCopy1, options, star1.m_calculateZetaNuclearErrorRadiusFlag, star1.m_calculateZetaNuclearErrorAgeFlag, r);
	star1.m_zetaSoberman		= ZadiabaticSPH(star1.m_Mass, star1.m_coreMass);
	star1.m_zetaSobermanHelium	= ZadiabaticSPH(star1.m_Mass, star1.m_HeCoreMass);
	star1.m_zetaHurley			= ZadiabaticHurley2002(star1.m_Mass, star1.m_coreMass);
	star1.m_zetaHurleyHelium	= ZadiabaticHurley2002(star1.m_Mass, star1.m_HeCoreMass);

	star2.m_zetaThermal 		= determineConvergedMassStepZetaThermal(starCopy2, options, star2.m_calculateZetaThermalErrorMassFlag, star2.m_calculateZetaThermalErrorRadiusFlag, r);
	star2.m_zetaNuclear			= determineConvergedTimestepZetaNuclear(starCopy2, options, star2.m_calculateZetaNuclearErrorRadiusFlag, star2.m_calculateZetaNuclearErrorAgeFlag, r);
	star2.m_zetaSoberman		= ZadiabaticSPH(star2.m_Mass, star2.m_coreMass);
	star2.m_zetaSobermanHelium	= ZadiabaticSPH(star2.m_Mass, star2.m_HeCoreMass);
	star2.m_zetaHurley			= ZadiabaticHurley2002(star2.m_Mass, star2.m_coreMass);
	star2.m_zetaHurleyHelium	= ZadiabaticHurley2002(star2.m_Mass, star2.m_HeCoreMass);

	if(debugging){
		std::cout << "Star\tType\tzetaThermal\tzetaNuclear\tzetaSPH\tZetaSPHHe\tzetaHurley\tzetaHurleyHe\tzetaSimple" << std::endl;
		std::cout << "1)\t" << 	star1.m_stellarType << "\t" << star1.m_zetaThermal << "\t" << 
								star1.m_zetaNuclear << "\t" << star1.m_zetaSoberman << "\t" << 
								star1.m_zetaSobermanHelium << "\t" << star1.m_zetaHurley << "\t" << 
								star1.m_zetaHurleyHelium << "\t" << star1.m_zetaSimple << std::endl;

		std::cout << "2)\t" << 	star2.m_stellarType << "\t" << star2.m_zetaThermal << "\t" << 
								star2.m_zetaNuclear << "\t" << star2.m_zetaSoberman << "\t" << 
								star2.m_zetaSobermanHelium << "\t" << star2.m_zetaHurley << "\t" << 
								star2.m_zetaHurleyHelium << "\t" << star2.m_zetaSimple << std::endl  << std::endl;
	}

}
//------------------------------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------------------------------//
void BinaryStar::XRayBinaries(const programOptions &options, std::ofstream & XRayData, std::ofstream & XPopData){
    /*
	// ALEJANDRO - 02/12/2016
	// Function to calculate and display all the values needed for the study of HMXB and ULX.
    // Based on Belczynski et al 2008 (https://arxiv.org/pdf/astro-ph/0511811.pdf), in turn based on Hurley+ 2002
	// Kit
    */
    bool    debugging = false;
    debugging = false;//true;


    double R2 = star2.m_Radius;                     //Radius of Donor star                                          //Rsol
    double Ma = star1.m_Mass;                       //Mass of accretor, will always be the primary                  //Msol
    double Md = star2.m_Mass;                       //Mass of Donor, will always be the secondary                   //Msol
    double Sa = star1.m_stellarType;                //Stellar tpye of accretor NS or BH
    double Sd = star2.m_stellarType;                //Stellar type of Donor
    double a = m_SemiMajorAxisPrime*AUToRsol;       //Seperation                                                    //Rsol
    double dmwind = star2.m_massLossDiff;           //DM wind                                                       //Msol
    double dt=m_time-m_timePrev;                    //timestep                                                      //Myr
    double ecc = m_EccentricityPrime;               //eccentricity
    double Mddot = 0.0;                                                                                                   //Msol per Myr
    Mddot = dmwind/dt;                              //Mass transfer rate from donor via  wind transfer

    double Vaccorb = 0.0;                           //accretors orbital velocity squared
    Vaccorb = G*Md/a;                               // Note that this is the velocity *squared* in units of??

    double Bwind = 0.0;                             //escape velocity coefficent dependent on stellar type and mass

    if(star2.m_stellarType == MS_LESS_THAN_07 or star2.m_stellarType == MS_MORE_THAN_07){   // Main Sequence
        if(Md<=1.4){
            Bwind=0.5;
        }else if(Md>=120.0){
            Bwind=7;
		}else{
            Bwind=0.5+6.5*(Md-1.4)/118.6;
        }
    }else if(star2.m_stellarType >= HERTZSPRUNG_GAP and star2.m_stellarType <= THERMALLY_PULSING_ASYMPTOTIC_GIANT_BRANCH){     // giants (Kdon = 2, 3, 4, 5, 6)
        Bwind = 0.125;
    }
    else if(star2.m_stellarType >= NAKED_HELIUM_STAR_MS and star2.m_stellarType <= NAKED_HELIUM_STAR_GIANT_BRANCH){           // He Rich stars (7,8,9)
        if(Md<=10){
            Bwind=0.125;
        }else if(Md>=120.0){
            Bwind=7.0;
        }else{
            Bwind=0.125+6.875*(Md-10.0)/110.0;
        }
    }
    // lots of 11 v 11 binaries
    // else{
    //     std::cerr << "Unimplemented donor stellar type in xray binary S1, S2 = " << star2.m_stellarType << " " << star2.m_stellarType << " setting Beta Wind = 0" << std::endl;
    //     Bwind = 0.0;
    // }

    double Vwind;                                   //wind velocity squared in units of ??
    Vwind = 2.0*Bwind*G*Md/R2;

    double V;                                       //ratio of velocities                                    //Dimensionless
    V = Vaccorb/Vwind;                                //Remember each of these is squared

    double Awind = 1.5;                             //accretion efficency
    double Madot;                                   //mass accretion rate
    Madot = (-Mddot*Awind*pow((G*Ma/(Vwind)),2.0))/(sqrt(1.0-pow(ecc,2.0))*2.0*pow(a,2.0)*pow((1.0 + V),1.5));                 //Msol per Myr
                                                    //Fwind stops accretion rate exceeding 0.8 donation rate

    // Check for NAN
    if(Madot != Madot){
        std::cerr << "Madot is NAN: " << Madot << " " << Mddot << " " << star1.m_stellarType << " " << star2.m_stellarType << " " << ecc << " " << a << " " << V << " " << Vwind << " " << std::endl;
        Madot = 0.0;
    }

    // Apply Fwind to limit accreted mass to 80% of donated mass
    if(Madot >= (-0.8*Mddot)){
        Madot = -0.8*Mddot;
		if(debugging){
			std::cout<<"Fwind"<<std::endl;
		}
    }
    else if(Madot < (-0.8*Mddot)){
        Madot = Madot; // 
    }
    else{
        std::cerr << "Xray binaries -- shouldn't get here: " << Madot << " " << -0.8*Mddot << std::endl;
        Madot = 0.0;
    }


    double Lum = 0.0;                                     //Luminosity
    double eta = 0.0;                                     //Matter to energy conversion efficency based off GPE binding energy ~M/2R
    // not bolometric luminosity correction as in Belczynski? in that case eta = 0.15 for NSs, 0.8 for BHs
    // Lbol = ??
    // Conversion to xray luminosity?

    // Where are these numbers from?
    if(Sa==BLACK_HOLE){
        eta=0.083;                              //Black Hole
    }else if(Sa==NEUTRON_STAR){
        eta=0.07;                               //Neutron Star
    }else{
        eta=0;                                  //primary is not a CO
    }

    // Kits formalism
    Lum=eta*Madot*pow(c,2.0);                         //Luminosity proptional to mdot*c^2                      //Msol per Myr * (m/s)^2
    Lum=Lum*Msol/(MyearToyear*year);                //Converting from Madot in Msol per Myr to kg per second    //watts
    Lum=Lum*JtoErg;                                 //Converting into erg/s       //Erg/s

    m_XRayLuminosity = Lum;                         // Store x-ray luminosity

    // Eddington luminosity?
    // Limit mass accretion rate to eddington rate?
    // MdotEdd = 4.0 * PI * G * M * mp / epsilon c sigmaT  // epsilon is fraction of GPE, sigmaT is thompson cross section (6.6E-29 m^2), mp is mass of proton (~1.67E-27 kg) (ionized hydrogen) 
    // LEdd ~ 1.3 × 10^38 (M1/Msol) erg s−1
    // Super eddington accretion?

    // Belczynski eq 83
    // epsilon G Macc Mdotacc / Racc (epsilon = 1.0 for NS, 0.5 for BH)

    // What about transient XRBS?

    // Calculate X-ray luminosity based on Zuo et al 2014 https://arxiv.org/abs/1310.5424
    double EddingtonLuminosity = 1.3E38 * star1.m_Mass;         // Eddington luminosity (eq 4 of Zuo et al 2014) in erg s-1
    double bolometricLuminosity = 0.1 * Madot * c * c;          // Just after equation 5 of Zuo et al 2014

    bolometricLuminosity = (bolometricLuminosity * Msol)/(MyearToyear * year);  // Converting from Madot in Msol per Myr to kg per second    //watts
    bolometricLuminosity = bolometricLuminosity * JtoErg;                       // Converting into erg/s       //Erg/s

    double etaEddington = 100.0;                        // Zuo et al use 100 for BHs, 5 for NSs -- check/vary this
    if(star1.m_stellarType == NEUTRON_STAR){
        etaEddington = 5.0;
    }
    double etaBolometric = 0.4;                         // The bolometric correction factor (from Zuo et al 2014, after Equation 5) converting the bolometric luminosity (Lbol) to the 0.5−8 keV X-ray luminosity (Belczynski & Taam 2004) observed by e.g. Chandra, adopted as 0.4 though its range is ∼0.1−0.5 for different types of XRB

    double xRayLuminosity = etaBolometric * std::min(bolometricLuminosity, etaEddington * EddingtonLuminosity);

    m_XRayLuminosity = xRayLuminosity;

}

void BinaryStar::BeXRayBinaries(const programOptions &options){
    /*
	// ALEJANDRO - 13/10/2017
	// Function to calculate and display all the values needed for the study of Be/X-Ray Binaries (BeXs)
	// Based on the XRayBinaries() function coded by Kit Boyett, Simon Stevenson and Alejandro Vigna-Gomez
    */
	
    bool    debugging = false; //true;

	// San and Alice: feel free to change variable names and add as many as you need. Also print statements.
	// This function is a safe zone and you can do whatever you want. Please don't modify other parts of the code without asking Coen or me first.
    double stellarTypeAccretor	= star1.m_stellarType;                //Stellar tpye of accretor
    double stellarTypeDonor		= star2.m_stellarType;                //Stellar type of donor
	
	if(debugging){
		std::cout << "Inside BeXs function" << std::endl;
		std::cout << "randomSeed:\t" << m_randomSeed << std::endl;
		std::cout << "StellarType 1,2:\t" << stellarTypeAccretor << "\t" << stellarTypeDonor << std::endl << std::endl;
	}
		


}

void BinaryStar::printingBinaryVariables(){
	// ALEJANDRO - 22/11/2016 - Tidying up print statements in evaluateBinary. Adding and reorganizing quantities encouraged, specially for ecperienced users of COMPAS ;)
		std::cout << "\n|-------------------------------------------------------------------------------|" << std::endl;
		std::cout << "printingBinaryVariables" << std::endl;
		std::cout << "TimePrev [Myrs]:\t" << m_timePrev << std::endl;
		std::cout << "Time\t [Myrs]:\t" << m_time << std::endl;
		std::cout << "Timestep [Myrs]:\t" << m_time-m_timePrev << std::endl;

		std::cout << "\nENERGIES" << std::endl;
		std::cout << "Eprev\t\t[]:\t" << m_TotalEnergyPrev << std::endl;
		std::cout << "E\t\t[]:\t" << m_TotalEnergyPrime << std::endl;
        std::cout << "LPrev\t\t[]:\t" << m_TotalAngularMomentumPrev << std::endl;
        std::cout << "L\t\t[]:\t" << m_TotalAngularMomentumPrime << std::endl;
		std::cout << "J_orbPrev\t[]:\t" << m_totalOrbitalAngularMomentumPrev << std::endl;
		std::cout << "J_orb\t\t[]:\t" << m_totalOrbitalAngularMomentumPrime << std::endl;
		std::cout << "E_orb\t\t[]:\t" << m_totalOrbitalEnergyPrime << std::endl;
		std::cout << "E_orbPrev\t[]:\t" << m_totalOrbitalEnergyPrev << std::endl;
		std::cout << "DeltaE_{orb}\t[]:\t" << m_totalOrbitalEnergyPrime - m_totalOrbitalEnergyPrev << std::endl;

		// Following quantities were made to test tides. Commented for now, may be un-commented at will.
        std::cout << "wp\t[]:\t"<< m_orbitalVelocityPrev << std::endl;
        std::cout << "dwML\t[]:\t"<< m_omegaMassLossDiff << std::endl;
        std::cout << "dwT\t[]:\t"<< m_omegaTidesDiff << std::endl;
        std::cout << "w\t[]:\t" << m_orbitalVelocityPrime << std::endl;
        std::cout << "w1p\t[]:\t" << star1.m_omegaPrev << std::endl;
        std::cout << "w2p\t[]:\t" << star2.m_omegaPrev << std::endl;
        std::cout << "w1\t[]:\t" << star1.m_omega << std::endl;
        std::cout << "w2\t[]:\t" << star2.m_omega << std::endl;
		

		std::cout << "\nORBIT" << std::endl;
        std::cout << "aPrev\t[AU]:\t"<< m_SemiMajorAxisPrev << std::endl;
        std::cout << "aPrev\t[Rsol]:\t"<< m_SemiMajorAxisPrev/RsolToAU << std::endl;
        std::cout << "a\t[AU]:\t"<< m_SemiMajorAxisPrime << std::endl;
		std::cout << "a\t[Rsol]:\t" << m_SemiMajorAxisPrime/RsolToAU << std::endl;
        std::cout << "wPrev\t[yr-1]:\t" << m_orbitalVelocityPrev << std::endl;
        std::cout << "w\t[yr-1]:\t" << m_orbitalVelocityPrime << std::endl;
		std::cout << "ePrev:\t\t" << m_EccentricityPrev << std::endl;
		std::cout << "e:\t\t" << m_Eccentricity << std::endl;							// Delete this one when properly incoporated to other parts of the code
		std::cout << "ePrime:\t\t" << m_EccentricityPrime << std::endl;

		std::cout << "\nPRIMARY" << std::endl;
		std::cout << "K1prev:\t\t" << star1.m_stellarTypePrev << std::endl;
		std::cout << "K1:\t\t" << star1.m_stellarType << std::endl;
		std::cout << "R1prev\t[Rsol]:\t" << star1.m_RadiusPrev << std::endl;
		std::cout << "R1\t[Rsol]:\t" << star1.m_Radius << std::endl;
		std::cout << "M1prev\t[Msol]:\t" << star1.m_MassPrev << std::endl;
		std::cout << "M1\t[Msol]:\t" << star1.m_Mass << std::endl;
		std::cout << "Menv1\t[Msol]:\t" << star1.m_envMass << std::endl;
		std::cout << "Mcore1\t[Msol]:\t" << star1.m_coreMass << std::endl;
		std::cout << "MHecore1[Msol]:\t" << star1.m_HeCoreMass << std::endl;
		std::cout << "MCOcore1[Msol]:\t" << star1.m_COCoreMass << std::endl;
		std::cout << "Age1\t[Myrs]:\t" << star1.m_Age << std::endl;

		std::cout << "\nSECONDARY" << std::endl;
		std::cout << "K2prev:\t\t" << star2.m_stellarTypePrev << std::endl;
		std::cout << "K2:\t\t" << star2.m_stellarType << std::endl;
		std::cout << "R2prev\t[Rsol]:\t" << star2.m_RadiusPrev << std::endl;
		std::cout << "R2\t[Rsol]:\t" << star2.m_Radius << std::endl;
		std::cout << "M2prev\t[Msol]:\t" << star2.m_MassPrev << std::endl;
		std::cout << "M2\t[Msol]:\t" << star2.m_Mass << std::endl;
		std::cout << "Menv2\t[Msol]:\t" << star2.m_envMass << std::endl;
		std::cout << "Mcore2\t[Msol]:\t" << star2.m_coreMass << std::endl;
		std::cout << "MHecore2[Msol]:\t" << star2.m_HeCoreMass << std::endl;
		std::cout << "MCOcore2[Msol]:\t" << star2.m_COCoreMass << std::endl;
		std::cout << "Age2\t[Myrs]:\t" << star2.m_Age << std::endl;
		std::cout << "|-------------------------------------------------------------------------------|" << std::endl;

}
//------------------------------------------------------------------------------------------------------------------------------------//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void BinaryStar::evaluateBinary(const programOptions &options, const gsl_rng *r, double dt, double dtPrev){

    bool    debugging = false;
	
	// ALEJANDRO - 14/11/2016 - After some profiling done by Jim, seems like calculate the Loveridge lambda at each timestep takes a lot of time. Therefore, I moved it to this function and only calculated it if explicitely indicated.
	if(options.lambdaCalculationEveryTimeStep){calculatingLambdas(options);}
	// ALEJANDRO - 16/10/2017 - Calculate zetas if specified
	if(options.zetaCalculationEveryTimeStep){calculatingZetas(options, r);}

	// ALEJANDRO - 16/11/2016 - Initiating flag, every timestep, to NO_MASS_TRANSFER. If it undergoes to MT or CEE, it should change.
	m_massTransferTrackerHistory = NO_MASS_TRANSFER;

    star1.flagRLOFPrev = star1.flagRLOF;
    star2.flagRLOFPrev = star2.flagRLOF;

	// Mass Transfer
    star1.flagRLOF = false;
    star2.flagRLOF = false;

    //Coen --10/11/2016-- used to check for updating formation printing
    same_RLOF_loop = false;

	// Find ratio of stars size to its Roche lobe radius, calculated at the point of closest approach, periapsis
    rocheLobeTracker1 = (star1.m_Radius*RsolToAU)/(m_SemiMajorAxisPrime*(1.0 - m_Eccentricity)*rocheLobeRadius(star1.m_Mass,star2.m_Mass));
    rocheLobeTracker2 = (star2.m_Radius*RsolToAU)/(m_SemiMajorAxisPrime*(1.0 - m_Eccentricity)*rocheLobeRadius(star2.m_Mass,star1.m_Mass));

    m_TotalAngularMomentumPrev = calculateAngularMomentum(m_SemiMajorAxisPrev, m_EccentricityPrev, star1.m_MassPrev, star2.m_MassPrev, star1.m_RadiusPrev, star2.m_RadiusPrev, star1.m_omegaPrev,star2.m_omegaPrev, m_orbitalVelocityPrev,k_definition(star1.m_RZAMS, star1.m_Radius, star1.m_Mass, star1.m_stellarType),k_definition(star2.m_RZAMS, star2.m_Radius, star2.m_Mass, star2.m_stellarType));

    // If star is larger than its Roche lobe, it is Roche lobe overflowing
    if(rocheLobeTracker1 >= 1.0){
		star1.flagRLOF = true;
		star1.experiencedRLOF = true;
	}
    if(rocheLobeTracker2 >= 1.0){
		star2.flagRLOF = true;
		star2.experiencedRLOF= true;
	}

	// In following statemets, do a m_commonEnvelopeFlag instigator for each mass
    // Flag check for RLOF of primary just after the CEE
    if((star1.flagRLOF == true)and(m_commonEnvelopeFlag == true)){m_RLOFPrimaryAfterCEE = true;}

    // Flag check for RLOF of secondary just after the CEE
    if((star2.flagRLOF == true)and(m_commonEnvelopeFlag == true)){m_RLOFSecondaryAfterCEE = true;}

     //JIM BARRETT - 18/11/2016 - We shouldn't be using doubles as Booleans. Issue 50
    //if(rocheLobeTracker1 == true or rocheLobeTracker2 == true){
    m_massTransferFlag=false;
    m_EccentricityRLOF = m_EccentricityPrime;

    if(star1.flagRLOF == true or star2.flagRLOF == true){

		if(options.circulariseBinaryDuringMassTransfer){
        // Circularise binary to the periapsis separation
    		if(debugging){
				std::cout << "Circularising binary" << std::endl;
    			std::cout << "Time: " << m_time << std::endl;
    			std::cout << "Semi-major axis prime: " << m_SemiMajorAxisPrime << std::endl;
    			std::cout << "Eccentricity: " << m_Eccentricity << std::endl;
    		}

			if(options.angularMomentumConservationDuringCircularisation == true){
				// Circulariste de a'=a*(1-e^2), which conserves angular momentum
				if(debugging){std::cout << "Conserving angular momentum while circularising the orbit due to mass transfer." << std::endl;}
				m_SemiMajorAxisPrime = m_SemiMajorAxisPrime * (1.0 - m_Eccentricity*m_Eccentricity);				
			}
			else{
				// Circularise to periastron, which does not conserves angular momentum. As used in StarTrack.
				if(debugging){std::cout << "Circularising to periastron. Angular momentum is not conserved if e > 0." << std::endl;}
				m_SemiMajorAxisPrime = m_SemiMajorAxisPrime * (1.0 - m_Eccentricity);
			}
			
			// ALEJANDRO - 22/11/2016 - Think shouldn't use m_Eccentricity but m_EccentrictyPrime. Right now setting both. Check later.
			m_Eccentricity = 0.0;
			m_EccentricityPrime = 0.0;

			// ALEJANDRO - 23/11/2016 - Bux fix for systems which enter MT being eccentric.
			// Previous values have to be the ones for periastron as later orbit is modified according to previous values. If you don't do this, you end up modifying pre-MT pre-circularisation orbit
			m_SemiMajorAxisPrev		= m_SemiMajorAxisPrime;
			m_EccentricityPrev		= m_EccentricityPrime;
			m_orbitalVelocityPrev	= m_orbitalVelocityPrime;

    		if(debugging){
    			std::cout << "Semi-major axis prime: " << m_SemiMajorAxisPrime << std::endl;
    			std::cout << "Eccentricity: " << m_Eccentricity << std::endl;
    		}
		}
		else{
			if(debugging){
				std::cout << m_randomSeed << "\tcircularise_binary_during_mass_transfer flag set to false. Do nothing to separation." << std::endl;
			}
		}			

		m_massTransferFlag=true;
    }


    massTransferTracker 		= 0.0;
    star1.m_MassTransferDiff 	= 0.0;
    star2.m_MassTransferDiff 	= 0.0;
    m_aMassTransferDiff 		= 0.0;
    m_omegaMassTransferDiff 	= 0.0;

    m_commonEnvelopeFlag = false;

    if(options.useMassTransfer == true){
        if(star1.flagRLOF==true and star2.flagRLOF==true){
			m_commonEnvelopeFlag = true;
			if(debugging){std::cout << "Contact system. CE/Merger, both stars RLOFlowing" << std::endl;}
        }
        else if(star1.flagRLOF==true or star2.flagRLOF==true){

			// COEN - 10/11/2016 - call updateFormationHistory function
            updateFormationHistory(star1.flagRLOF, star2.flagRLOF,same_RLOF_loop, m_commonEnvelopeFlag, star1.flagExperiencedCCSN, star2.flagExperiencedCCSN, star1.flagExperiencedECSN, star2.flagExperiencedECSN, star1.m_stellarType, star2.m_stellarType,  eventCounter,  mt_primary_counter, mt_secondary_counter, mt_primary_ep1, mt_primary_ep1_K1,  mt_primary_ep1_K2, mt_primary_ep2, mt_primary_ep2_K1,  mt_primary_ep2_K2, mt_primary_ep3, mt_primary_ep3_K1, mt_primary_ep3_K2, mt_secondary_ep1, mt_secondary_ep1_K1, mt_secondary_ep1_K2, mt_secondary_ep2, mt_secondary_ep2_K1, mt_secondary_ep2_K2, mt_secondary_ep3, mt_secondary_ep3_K1, mt_secondary_ep3_K2, SN_primary_type_1, SN_primary_type_2, SN_primary_type_3, SN_secondary_type_1, SN_secondary_type_2, SN_secondary_type_3, CEE,    CEE_instigator, CEE_failed, CEE_failed_instigator, CEE_wet, CEE_wet_instigator, stellar_type_K1,  stellar_type_K2, m_stellarMerger, star1.m_stellarTypePrev, star2.m_stellarTypePrev);
            same_RLOF_loop = true;

			if(star1.flagRLOF==true){
                if(debugging){std::cout << "\nRLOF Star 1";}
            }
            else{
                if(debugging){std::cout << "\nRLOF Star 2";}
            }

            if(debugging){
				std::cout <<std::endl << "Begin Mass Transfer" <<std::endl;
				printingBinaryVariables();
			}
                
                if(options.RLOFPrinting){
                    //RLOFp.printRLOFParameters(options.outputPath, "BeforeRLOF.txt", *this);
                    }
                
                MassTransfer(dt, dtPrev, options, r); // think i can remove dt , dtPrev, m_TotalAngularMomentumPrev quantities from the function               
                if(options.RLOFPrinting){
                    //RLOFprinting.printRLOFParameters(options.outputPath, "AfterRLOF.txt", *this);
                    }
                
			if(debugging){
				std::cout << "After Mass Transer:" << std::endl;
				printingBinaryVariables();
			}
				if(debugging){std::cout << "End Mass Transfer" <<std::endl<<std::endl;}
        }
        else{
            if(debugging){std::cout << "No mass transfer" << std::endl;}
        }
    }


// End Mass Transfer
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Winds Mass loss evolution of binary
    m_aMassLossDiff = 0.0;
    m_omegaMassLossDiff = 0.0;
    double  newMassWinds1 = 0.0;
    double  newMassWinds2 = 0.0;

    if(options.useMassTransfer == true and m_massTransferFlag == true){
    // If statement used for halting winds when in mass transfer (first approach).
            star1.m_massLossDiff = 0.0;
            star2.m_massLossDiff = 0.0;
    }
    else{
        if(options.useMassLoss == true){
            bool    debugging = false;      // Value specified in "massLossPrescriptionHurley" function in star.cpp file. Set to false.
            // should really just pass this thing options rather than lots of bits of options
            //JIM BARRETT - 28/09/2016 - Changed the LBV and WR multipliers to use the BinaryStar member variable
            newMassWinds1 = star1.massLossCaller(star1.m_Mass, star1.m_Mass0, star1.m_coreMass, star1.m_Luminosity, star1.m_Radius, star1.m_Mu, star1.m_Metallicity, star1.m_Temperature, star1.m_logMetallicityXi, star1.m_Mdot, dt, star1.m_Age, m_luminousBlueVariableFactor, m_wolfRayetFactor, star1.m_stellarType, options.useMassLoss, false, options.massLossPrescription, debugging, star1.m_an_coefficients, star1.m_timescales);
            newMassWinds2 = star2.massLossCaller(star2.m_Mass, star2.m_Mass0, star2.m_coreMass, star2.m_Luminosity, star2.m_Radius, star2.m_Mu, star2.m_Metallicity, star2.m_Temperature, star2.m_logMetallicityXi, star2.m_Mdot, dt, star2.m_Age, m_luminousBlueVariableFactor, m_wolfRayetFactor, star2.m_stellarType, options.useMassLoss, false, options.massLossPrescription, debugging, star2.m_an_coefficients, star2.m_timescales);

            double  newSemiMajorAxisWinds = aCorrectionMassLoss(m_SemiMajorAxisPrev, star1.m_MassPrev, star2.m_MassPrev,newMassWinds1, newMassWinds2);
            star1.m_massLossDiff = newMassWinds1-star1.m_Mass;
            star2.m_massLossDiff = newMassWinds2-star2.m_Mass;
            m_aMassLossDiff = -m_SemiMajorAxisPrev + newSemiMajorAxisWinds;
            m_omegaMassLossDiff = - m_orbitalVelocityPrev + sqrt(G1*(newMassWinds1+newMassWinds2)/pow(newSemiMajorAxisWinds,3.0));
        }
    }
// Ends Winds Mass loss evolution of binary
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Common Envelope Event and Supernovas: both classified as immediate events
	//    debugging = true;
    if(m_commonEnvelopeFlag == true or star1.flagSN == true or star2.flagSN ==true or star1.flagECSN == true or star2.flagECSN == true or m_stellarMerger == true){

        //debugging = true;
		// ALEJANDRO - 08/05/2018 - Don't understand why which_star was made a double. It should be an int type. Ideally, a boolean: primaryExplodes or secondaryExplodes.
        double which_star = 0; //passed to supernova() to keep track of which star
                               // explodes 1=primary, 2=secondary, 3=both

		// Check is a CCE or a stellar merger occured; if so, get into the CEE function.
        if(m_commonEnvelopeFlag == true || m_stellarMerger == true){
            if(debugging){
                std::cout << "Before CEE" << std::endl;
				std::cout << "m_commonEnvelopeFlag:\t" << m_commonEnvelopeFlag << std::endl;
				std::cout << "m_stellarMerger:\t" << m_stellarMerger << std::endl;
				printingBinaryVariables();
            }
                       // std::cout << "common envelope event funtion was called " << std::endl;
			// Call the function that evaluates CEE
			CommonEnvelopeEvent(r, options);

			if(debugging){
                std::cout << "After CEE" << std::endl;
				printingBinaryVariables();
                debugging = false;
            }
        }
        else if((star1.flagSN == true or star1.flagECSN == true) and (star2.flagSN == true or star2.flagECSN == true) and m_SemiMajorAxisPrime > 0.0){
            which_star = 3;
            if(debugging){std::cout << "Both stars went supernova at the same time. (1st check)" << std::endl;}
            supernova(options, r, star1, star2, m_randomSeed, which_star); // Calculate one supernova first

			// ALEJANDRO - 08/05/2018 - As we want to be able to keep evolving runaway stars, not checking for bound/unbound in order to print to supernovae.txt file.
            supernova(options, r, star2, star1, m_randomSeed, which_star);
            /*
			if(m_EPrime > 0.0){
                if(debugging){std::cout << "Unbound, no need to calculate second supernova" << std::endl;}
                star2.flagSN = false;
                star2.flagECSN = false;
            }
			
            else{
               if(debugging){std::cout << "Survived first supernova, calculate second" << std::endl;}
                supernova(options, r, star2, star1, m_randomSeed, which_star);
             }
			*/
            if(debugging){std::cout << "End of double supernova" << std::endl;}
        }
		// else if((star1.flagSN == true or star1.flagECSN == true) and m_SemiMajorAxisPrime > 0.0){
        else if((star1.flagSN == true or star1.flagECSN == true)){
            which_star = 1;
            if(debugging){
                std::cout << "Supernova1 begins (1st check)\n";
                std::cout << "Kick1 before: " << star1.m_kickVelocity << std::endl;
                std::cout << "a) \n";
            }
            supernova(options, r, star1, star2, m_randomSeed, which_star);
            if(debugging){
                std::cout << "Supernova1 Ends\n";
            }
        }
        // else if((star2.flagSN == true or star2.flagECSN == true) and m_SemiMajorAxisPrime > 0.0){
        else if((star2.flagSN == true or star2.flagECSN == true)){
            which_star=2;

            if(debugging){
                std::cout << "Supernova2 begins (1st check)\n";
                std::cout << "Kick2 before: " << star2.m_kickVelocity << std::endl;
                std::cout << "b) \n";
            }
            supernova(options, r, star2, star1, m_randomSeed,which_star);
            if(debugging){
                std::cout << "Supernova2 ends\n";
            }
        }
        else{
            if(debugging){std::cout << "Possible several conditions between SN1, SN2 and CEE. Check.\n";}
        }
        debugging = false;
    }
	// End of inmediate events, such as Common Envelope Event and Supernovas
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    else{
	// No CEE nor SNe
	// Update mass of each star.
        if(debugging){
            std::cout << "\nPre Mass: \n";
			printingBinaryVariables();
        }

		star1.m_Mass  = star1.m_MassPrev + star1.m_massLossDiff + star1.m_MassTransferDiff;
		star2.m_Mass  = star2.m_MassPrev + star2.m_massLossDiff + star2.m_MassTransferDiff;

        if(debugging){
            std::cout << "After Mass: \n";
			printingBinaryVariables();
        }

		// Update the age of the star according to winds and/or mass transfer loss
		//    star1.m_Age = star1.timeChangeAfterMassLoss(star1.m_an_coefficients, star1.m_Mass0, star1.m_Mass, star1.m_coreMass, star1.m_Age, star1.m_timescales, star1.m_stellarType, star1.m_logMetallicityXi);
		//    star2.m_Age = star2.timeChangeAfterMassLoss(star2.m_an_coefficients, star2.m_Mass0, star2.m_Mass, star2.m_coreMass, star2.m_Age, star2.m_timescales, star2.m_stellarType, star2.m_logMetallicityXi);

        if(debugging){std::cout << "PRE Age1, Age2: " << star1.m_Age << " " << star2.m_Age << std::endl;}
        double agePrime1 = timeChangeAfterMassLoss(star1.m_an_coefficients, star1.m_Mass0, star1.m_Mass, star1.m_coreMass, star1.m_Age, star1.m_timescales, star1.m_stellarType, star1.m_logMetallicityXi, false);
        double agePrime2 = timeChangeAfterMassLoss(star2.m_an_coefficients, star2.m_Mass0, star2.m_Mass, star2.m_coreMass, star2.m_Age, star2.m_timescales, star2.m_stellarType, star2.m_logMetallicityXi, false);
        double f_rej_1 = massTransferRejuvenationFactor(star1, options);
        double f_rej_2 = massTransferRejuvenationFactor(star2, options);

        if(debugging){
            std::cout << "agePrime1 = " << agePrime1 << std::endl;
            std::cout << "agePrime2 = " << agePrime2 << std::endl;
            std::cout << "f_rej_1 = " << f_rej_1 << std::endl;
            std::cout << "f_rej_2 = " << f_rej_2 << std::endl;
            std::cout << "Updating stars ages" << std::endl;
        }
        star1.m_Age = agePrime1 * f_rej_1;
        star2.m_Age = agePrime2 * f_rej_2;

        if(debugging){std::cout << "AFF Age1, Age2: " << star1.m_Age << " " << star2.m_Age << std::endl;}

        if(debugging){
            std::cout << "Post Age: \n";
			printingBinaryVariables();
        }

		// Commenting star rotation and doing it at the end, after tides
		/*
		star1.m_angularMomentum = k_definition(star1.m_Mass,star1.m_Radius,star1.m_RZAMS)*pow(star1.m_Radius*RsolToAU,2.0)*star1.m_omegaPrev;
		star2.m_angularMomentum = k_definition(star2.m_Mass,star2.m_Radius,star2.m_RZAMS)*pow(star2.m_Radius*RsolToAU,2.0)*star2.m_omegaPrev;
		star1.m_omega = star1.m_omegaPrev + star1.m_omegaMagneticBrakingDiff + star1.m_omegaTidesIndividualDiff;
		star2.m_omega = star2.m_omegaPrev + star2.m_omegaMagneticBrakingDiff + star2.m_omegaTidesIndividualDiff;
		*/

		// Update binary
		if(debugging){
			std::cout << "\nPre-modify orbit" << std::endl;
			std::cout << "m_SemiMajorAxisPrev:\t" << m_SemiMajorAxisPrev << std::endl;
			std::cout << "m_SemiMajorAxisPrime:\t" << m_SemiMajorAxisPrime << std::endl;
			std::cout << "m_aMassLossDiff:\t" << m_aMassLossDiff << std::endl;
			std::cout << "m_aMassTransferDiff:\t" << m_aMassTransferDiff << std::endl;
		}
		m_orbitalVelocityPrime	= m_orbitalVelocityPrev + m_omegaMassLossDiff + m_omegaMassTransferDiff; //should here be a diff quantity because of MB?
		m_SemiMajorAxisPrime    = m_SemiMajorAxisPrev + m_aMassLossDiff + m_aMassTransferDiff;
		// m_EccentricityPrime
		if(debugging){
			std::cout << "Post-modify orbit" << std::endl;
			std::cout << "m_SemiMajorAxisPrev:\t" << m_SemiMajorAxisPrev << std::endl;
			std::cout << "m_SemiMajorAxisPrime:\t" << m_SemiMajorAxisPrime << std::endl;
		}

		// ALEJANDRO - 16/11/2016 - Calculating orbital energy and angular momentum
		m_totalMassPrev						= m_totalMassPrime;
		m_reducedMassPrev					= m_reducedMassPrime;
		m_totalOrbitalEnergyPrev			= m_totalOrbitalEnergyPrime;
		m_totalOrbitalAngularMomentumPrev	= m_totalOrbitalAngularMomentumPrime;

		m_totalMassPrime					= star1.m_Mass+star2.m_Mass;
		m_reducedMassPrime					= (star1.m_Mass*star2.m_Mass)/m_totalMassPrime;
		m_totalOrbitalEnergyPrime			= orbitalEnergy(m_reducedMassPrime, m_totalMassPrime, m_SemiMajorAxisPrime);
		m_totalOrbitalAngularMomentumPrime 	= orbitalAngularMomentum(m_reducedMassPrime, m_totalMassPrime, m_SemiMajorAxisPrime);

		// ALEJANDRO - 16/11/2016 - Calculating energy and angular momentum using regular conservation of energy, specially useful for checking tides and rotational effects.
		m_TotalEnergyPrime = calculateTotalEnergy(m_SemiMajorAxisPrime, star1.m_Mass, star2.m_Mass, star1.m_Radius, star2.m_Radius, star1.m_omega, star2.m_omega, m_orbitalVelocityPrime, k_definition(star1.m_RZAMS, star1.m_Radius, star1.m_Mass, star1.m_stellarType),k_definition(star2.m_RZAMS, star2.m_Radius, star2.m_Mass, star2.m_stellarType));
		m_TotalAngularMomentumPrime = calculateAngularMomentum(m_SemiMajorAxisPrime, m_EccentricityPrime, star1.m_Mass, star2.m_Mass, star1.m_Radius, star2.m_Radius, star1.m_omega,star2.m_omega, m_orbitalVelocityPrime, k_definition(star1.m_RZAMS, star1.m_Radius, star1.m_Mass, star1.m_stellarType),k_definition(star2.m_RZAMS, star2.m_Radius, star2.m_Mass, star2.m_stellarType));

        if(debugging){
			printingBinaryVariables();
        }
    }
	

	//
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Update stellar values
	//    debugging = true;
		if(debugging){
			std::cout << "\nPre evolve\n";
			printingBinaryVariables();
		}
        
		star1.evolveOneTimestep(0.0, false, options, r);
		star2.evolveOneTimestep(0.0, false, options, r);

		//debugging = true;
        double which_star = 0.0; //passed to supernova() to keep track of which star
                               // explodes 1=orimary, 2=secondary, 3=both
		if((star1.flagSN == true or star1.flagECSN == true) and (star2.flagSN == true or star2.flagECSN == true) and m_SemiMajorAxisPrime > 0.0){
            which_star=3.0;
			//std::cerr  << m_randomSeed << "\tError in both star Supernova, star core is being evolved in BinaryStar" << std::endl;
			if(debugging){std::cout << "Both stars went supernova at the same time. (2nd check)" << std::endl;}
			supernova(options, r, star1, star2, m_randomSeed, which_star); // Calculate one supernova first
				if(m_EPrime > 0.0){
				if(debugging){std::cout << "Unbound, no need to calculate second supernova" << std::endl;}
				star2.flagSN = false;
				star2.flagECSN = false;
			}
			else{
				if(debugging){std::cout << "Survived first supernova, calculate second" << std::endl;}
				supernova(options, r, star2, star1, m_randomSeed, which_star);
			}
			if(debugging){std::cout << "End of double supernova" << std::endl;}
		}
		else if((star1.flagSN == true or star1.flagECSN == true) and m_SemiMajorAxisPrime > 0.0){
             which_star=1.0;
			//std::cerr << m_randomSeed << "\tError in star1 Supernova, star core is being evolved in BinaryStar" << std::endl;
			if(debugging){std::cout << "Supernova1 begins (2nd check)\n";}
			supernova(options, r, star1, star2, m_randomSeed, which_star);
			if(debugging){std::cout << "Supernova1 Ends\n";}
		}
		else if((star2.flagSN == true or star2.flagECSN == true) and m_SemiMajorAxisPrime > 0.0){
            which_star=2.0;
			//std::cerr << m_randomSeed << "\tError in star2 Supernova, star core is being evolved in BinaryStar" << std::endl;
			if(debugging){std::cout << "Supernova2 begins (2nd check)\n";}
			supernova(options, r, star2, star1, m_randomSeed, which_star);
			if(debugging){std::cout << "Supernova2 Ends\n";}
		}

		if(debugging){
			std::cout << "Post evolve\n";
			printingBinaryVariables();
		}

		// Assign new values to "previous" values, for following timestep
		m_EccentricityPrev	= m_EccentricityPrime;
		m_SemiMajorAxisPrev = m_SemiMajorAxisPrime;
		m_orbitalVelocityPrev = m_orbitalVelocityPrime;

		// COEN - 10/11/2016 - call updateFormationHistory function
		updateFormationHistory(star1.flagRLOF, star2.flagRLOF,same_RLOF_loop, m_commonEnvelopeFlag, star1.flagExperiencedCCSN, star2.flagExperiencedCCSN, star1.flagExperiencedECSN, star2.flagExperiencedECSN,star1.m_stellarType, star2.m_stellarType,  eventCounter,  mt_primary_counter, mt_secondary_counter, mt_primary_ep1, mt_primary_ep1_K1,  mt_primary_ep1_K2, mt_primary_ep2, mt_primary_ep2_K1,  mt_primary_ep2_K2, mt_primary_ep3, mt_primary_ep3_K1, mt_primary_ep3_K2, mt_secondary_ep1, mt_secondary_ep1_K1, mt_secondary_ep1_K2, mt_secondary_ep2, mt_secondary_ep2_K1, mt_secondary_ep2_K2, mt_secondary_ep3, mt_secondary_ep3_K1, mt_secondary_ep3_K2, SN_primary_type_1, SN_primary_type_2, SN_primary_type_3, SN_secondary_type_1, SN_secondary_type_2, SN_secondary_type_3, CEE,    CEE_instigator, CEE_failed, CEE_failed_instigator, CEE_wet, CEE_wet_instigator, stellar_type_K1,  stellar_type_K2, m_stellarMerger, star1.m_stellarTypePrev, star2.m_stellarTypePrev);
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Tides evolution of binary
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		m_aTidesDiff = 0.0;
		m_omegaTidesDiff = 0.0;
		star1.m_omegaTidesIndividualDiff = 0.0;
		star2.m_omegaTidesIndividualDiff = 0.0;

		if(debugging){
			std::cout << "Before Tides" << std::endl;
			printingBinaryVariables();
		}

		if(m_orbitalVelocityPrime > 0.0){
			if(options.tidesPrescription != TIDES_PRESCRIPTION_NONE){
				if(debugging){std::cout << "Begin tides" << std::endl;}
				tides(options);
				if(debugging){std::cout << "End tides" << std::endl << std::endl;}
			}
		}
		else{
			if(debugging){std::cout << m_randomSeed <<  "\tm_orbitalVelocityPrime <= 0.0" << std::endl;}
		}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// End Tides evolution of binary
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


		if(star1.m_stellarType <= HELIUM_WHITE_DWARF){
			star1.m_angularMomentum = k_definition(star1.m_RZAMS, star1.m_Radius, star1.m_Mass, star1.m_stellarType)*pow(star1.m_Radius*RsolToAU,2.0)*star1.m_omega;
		}

		if(star2.m_stellarType <= HELIUM_WHITE_DWARF){
			star2.m_angularMomentum = k_definition(star2.m_RZAMS, star2.m_Radius, star2.m_Mass, star2.m_stellarType)*pow(star2.m_Radius*RsolToAU,2.0)*star2.m_omega;
		}
		

		star1.m_omega += star1.m_omegaTidesIndividualDiff;
		star2.m_omega += star2.m_omegaTidesIndividualDiff;

		// Update binary
		m_orbitalVelocityPrime	+= m_omegaTidesDiff; //should here be a diff quantity because of MB?
		m_SemiMajorAxisPrime    += m_aTidesDiff;

		// ALEJANDRO - 16/11/2016 - Calculating orbital energy and angular momentum
		m_totalMassPrev						= m_totalMassPrime;
		m_reducedMassPrev					= m_reducedMassPrime;
		m_totalOrbitalEnergyPrev			= m_totalOrbitalEnergyPrime;
		m_totalOrbitalAngularMomentumPrev	= m_totalOrbitalAngularMomentumPrime;

		m_totalMassPrime					= star1.m_Mass+star2.m_Mass;
		m_reducedMassPrime					= (star1.m_Mass*star2.m_Mass)/m_totalMassPrime;
		m_totalOrbitalEnergyPrime			= orbitalEnergy(m_reducedMassPrime, m_totalMassPrime, m_SemiMajorAxisPrime);
		m_totalOrbitalAngularMomentumPrime 	= orbitalAngularMomentum(m_reducedMassPrime, m_totalMassPrime, m_SemiMajorAxisPrime);

		// ALEJANDRO - 16/11/2016 - Calculating energy and angular momentum using regular conservation of energy, specially useful for checking tides and rotational effects
		m_TotalEnergyPrime = calculateTotalEnergy(m_SemiMajorAxisPrime, star1.m_Mass, star2.m_Mass, star1.m_Radius, star2.m_Radius, star1.m_omega, star2.m_omega, m_orbitalVelocityPrime, k_definition(star1.m_RZAMS, star1.m_Radius, star1.m_Mass, star1.m_stellarType),k_definition(star2.m_RZAMS, star2.m_Radius, star2.m_Mass, star2.m_stellarType));
		m_TotalAngularMomentumPrime = calculateAngularMomentum(m_SemiMajorAxisPrime, m_EccentricityPrime, star1.m_Mass, star2.m_Mass, star1.m_Radius, star2.m_Radius, star1.m_omega,star2.m_omega, m_orbitalVelocityPrime, k_definition(star1.m_RZAMS, star1.m_Radius, star1.m_Mass, star1.m_stellarType),k_definition(star2.m_RZAMS, star2.m_Radius, star2.m_Mass, star2.m_stellarType));

		if(debugging){
			std::cout << "After Tides" << std::endl;
			printingBinaryVariables();
			std::cout << std::endl;
		}
        
		 //std::cout << "\nEnd of EvaluateBinary:" << std::endl;
		// printingBinaryVariables();
		//
		
		// Update pulsar parameters
		// 
		
		double epsilonPulsar = 1.0;
		double alphaPulsar = 0.0;

		// For star 1
		if(star1.m_stellarType == NEUTRON_STAR){

		
			updateMagneticFieldAndSpin( star1.flagRecycledNS, m_commonEnvelopeFlag,  m_dt * MyearToyear * year, star1.m_momentInertia * g_to_kg * pow (cm_to_m,2.0), star1.m_Mass * Msol, star1.m_Radius * Rsol, star1.m_MassTransferDiff * Msol, options.pulsarMagneticFieldDecayMassscale * Msol , epsilonPulsar, options.pulsarMagneticFieldDecayTimescale * MyearToyear * year, pow(10.0, options.pulsarLog10MinimumMagneticField) * gauss_to_tesla,  star1.m_angularMomentum, star1.m_pulsarSpinFrequency, star1.m_pulsarMagneticField, star1.m_pulsarSpinDownRate, alphaPulsar);	
		}

		// For star 2
		if(star2.m_stellarType == NEUTRON_STAR){
                        

			updateMagneticFieldAndSpin( star2.flagRecycledNS, m_commonEnvelopeFlag, m_dt * MyearToyear * year, star2.m_momentInertia * g_to_kg * pow (cm_to_m,2.0), star2.m_Mass, star2.m_Radius * Rsol, star2.m_MassTransferDiff * Msol, options.pulsarMagneticFieldDecayMassscale * Msol ,  epsilonPulsar, options.pulsarMagneticFieldDecayTimescale * MyearToyear * year, pow(10.0, options.pulsarLog10MinimumMagneticField) * gauss_to_tesla,  star2.m_angularMomentum, star2.m_pulsarSpinFrequency,  star2.m_pulsarMagneticField, star2.m_pulsarSpinDownRate, alphaPulsar);
		}
}
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
