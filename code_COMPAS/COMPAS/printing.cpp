#include "printing.h"

void printDoubleCompactObjectsHeader(const boost::filesystem::path &outputPath, std::string filename){

    // the / operator does concatenation
    std::ofstream dcos((outputPath/filename).string());
    
	// ALEJANDRO - 22/05/2018 - Headers checked. Dimension match (74 columns) with units and variable type.
    dcos << "#" << TAB <<
    		"#" << TAB << 
			"#" << TAB <<
			"#" << TAB <<  
			"AU" << TAB << 
			"#" << TAB << 
			"AU" << TAB << 
			"#" << TAB << 
			"kms^-1" << TAB << 
			"AU" << TAB << 
			"#" << TAB << 
			"#" << TAB << 
			"#" << TAB << 
			"Msol" << TAB << 
			"Msol" << TAB << 
			"Msol" << TAB << 
			"Msol" << TAB << 
			"Msol" << TAB << 
			"Msol" << TAB << 
			"Msol" << TAB << 
			"Msol" << TAB << 
			"kms^-1" << TAB << 
			"#" << TAB << 
			"#" << TAB << 
			"Msol" << TAB << 
			"#" << TAB << 
			"Msol" << TAB << 
			"Msol" << TAB << 
			"Msol" << TAB << 
			"Msol" << TAB << 
			"Msol" << TAB << 
			"Msol" << TAB << 
			"Msol" << TAB << 
			"Msol" << TAB << 
			"kms^-1" << TAB << 
			"#" << TAB << 
			"#" << TAB << 
			"Msol" << TAB << 
			"#" << TAB << 
			"Myr" << TAB << 
			"Myr" << TAB << 
			"#" << TAB << 
			"kms^-1" << TAB << 
			"kms^-1" << TAB << 
			"#" << TAB << 
			"#" << TAB << 
			"#" << TAB << 
			"#"<< TAB << 
			"#" << TAB << 
			"#" << TAB << 
			"Msol^2AU^-1" << TAB << 
			"Msol^2AU^-1" << TAB << 
			"Msol^2AU^-1" << TAB << 
			"Msol^2AU^-1" << TAB << 
			"#" << TAB << 
			"#" << TAB << 
			"#" << TAB << 
			"#" << TAB << 
			"#" << TAB << 
			"Rsol" << TAB << 
			"Rsol" << TAB << 
			"Rsol" << TAB << 
			"Rsol" << TAB << 
			"Rsol" << TAB << 
			"Rsol" << TAB << 
			"#" << TAB << 
			"#" << TAB << 
			"#" << TAB << 
			"ergs" << TAB << 
			"ergs" << TAB << 
			"#" << TAB << 
			"#" << TAB << 
			"#" << TAB << 
			"#" << TAB << 
			"#" << TAB << 
			"#" << TAB <<
			"#" << TAB <<
			"#" << TAB <<
			"#" << TAB <<
			"#" << std::endl;

    dcos << "int" << TAB << 
			"int" << TAB << 
			"float" << TAB <<
			"int" << TAB <<
			"float" << TAB << 
			"float" << TAB << 
			"float" << TAB << 
			"float" << TAB << 
			"float" << TAB << 
			"float" << TAB << 
			"float" << TAB << 
			"float" << TAB << 
			"float" << TAB << 
			"float" << TAB << 
			"float" << TAB << 
			"float" << TAB << 
			"float" << TAB << 
			"float" << TAB << 
			"float" << TAB << 
			"float" << TAB << 
			"float" << TAB << 
			"float" << TAB << 
			"float" << TAB << 
			"float" << TAB << 
			"float" << TAB << 
			"int" << TAB << 
			"float" << TAB << 
			"float" << TAB << 
			"float" << TAB << 
			"float" << TAB << 
			"float" << TAB << 
			"float" << TAB << 
			"float" << TAB << 
			"float" << TAB << 
			"float" << TAB << 
			"float" << TAB << 
			"float" << TAB << 
			"float" << TAB << 
			"int" << TAB << 
			"float" << TAB << 
			"float" << TAB << 
			"float" << TAB << 
			"float" << TAB << 
			"float" << TAB << 
			"float" << TAB << 
			"float" << TAB << 
			"float" << TAB << 
			"bool" << TAB << 
			"int" << TAB << 
			"int" << TAB << 
			"float" << TAB << 
			"float" << TAB << 
			"float" << TAB << 
			"float" << TAB << 
			"bool" << TAB << 
			"float" << TAB << 
			"float" << TAB << 
			"float" << TAB << 
			"float" << TAB << 
			"float" << TAB << 
			"float" << TAB << 
			"float" << TAB << 
			"float" << TAB << 
			"float" << TAB << 
			"float" << TAB << 
			"bool" << TAB << 
			"bool" << TAB << 
			"bool" << TAB << 
			"float" << TAB << 
			"float" << TAB << 
			"bool" << TAB << 
			"bool" << TAB << 
			"bool" << TAB << 
			"bool" << TAB << 
			"bool" << TAB << 
			"bool" << TAB <<
			"bool" << TAB << 
			"bool" << TAB << 
			"bool" << TAB << 
			"bool" << std::endl;

    dcos << "ID" << TAB << 
			"seed" << TAB << 
			"weight" << TAB <<
			"samplingPhase" << TAB <<
			"separationInitial" << TAB << 
			"eccentricityInitial" << TAB << 
			"separationPrior2ndSN" << TAB << 
			"eccentricityPrior2ndSN" << TAB << 
			"relativeVelocity2ndSN" << TAB << 
			"separationDCOFormation" << TAB << 
			"eccentricityDCOFormation" << TAB << 
			"Metallicity1" << TAB << 
			"Metallicity2" << TAB << 
			"M1ZAMS" << TAB << 
			"totalMassDCOFormation1" << TAB << 
			"HeCoreMassDCOFormation1" << TAB << 
			"COCoreMassDCOFormation1" << TAB << 
			"coreMassDCOFormation1" << TAB << 
			"HeCoreMassCE1" << TAB << 
			"COcoreMassCE1" << TAB << 
			"coreMassCE1" << TAB << 
			"drawnKick1" << TAB << 
			"thetaSupernova1" << TAB << 
			"phiSupernova1" << TAB << 
			"M1" << TAB << 
			"stellarType1" << TAB << 
			"M2ZAMS" << TAB << 
			"totalMassDCOFormation2" << TAB << 
			"HeCoreMassDCOFormation2" << TAB << 
			"COCoreMassDCOFormation2" << TAB << 
			"coreMassDCOFormation2" << TAB << 
			"HeCoreMassCE2" << TAB << 
			"COCoreMassCE2" << TAB << 
			"coreMassCE2" << TAB << 
			"drawnKick2" << TAB << 
			"thetaSupernova2" << TAB << 
			"phiSupernova2" << TAB << 
			"M2" << TAB << 
			"stellarType2" << TAB << 
			"tc" << TAB << 
			"tform" << TAB << 
			"flbv" << TAB << 
			"sigmaKickNS" << TAB << 
			"sigmaKickBH" << TAB << 
			"CEalpha" << TAB << 
			"kickDirectionPower" << TAB << 
			"wolfRayetMultiplier" << TAB << 
			"RLOFSecondaryAfterCEE" << TAB << 
			"PrimaryMTCase" << TAB << 
			"SecondaryMTCase" << TAB << 
			"preSNeOrbitalEnergy1" << TAB << 
			"postSNeOrbitalEnergy1" << TAB << 
			"preSNeOrbitalEnergy2" << TAB << 
			"postSNeOrbitalEnergy2" << TAB << 
			"CEflag" << TAB << 
			"lambda1CE" << TAB << 
			"lambda2CE" << TAB << 
			"EccentricityPreCEE" << TAB << 
			"EccentricityPostCEE" << TAB << 
			"SemiMajorAxisPreCEE" << TAB << 
			"SemiMajorAxisPostCEE" << TAB << 
			"RL1to2PreCEE" << TAB << 
			"RL1to2PostCEE" << TAB << 
			"RL2to1PreCEE" << TAB << 
			"RL2to1PostCEE" << TAB << 
			"optimisticCEFlag" << TAB << 
			"mergesInHubbleTimeFlag" << TAB << 
			"doubleCommonEnvelopeFlag" << TAB << 
			"bindingEnergyCEEStar1" << TAB << 
			"bindingEnergyCEEStar2" << TAB << 
			"recycledPrimary" << TAB << 
			"recycledSecondary" << TAB << 
			"USSNPrimary" << TAB << 
			"USSNSecondary" << TAB << 
			"ECSNPrimary" << TAB << 
			"ECSNSecondary"  << TAB <<
			"PISNPrimary" << TAB <<
			"PISNSecondary" << TAB <<
			"PPISNPrimary" << TAB <<
			"PPISNSecondary" << std::endl;

    dcos.close();

}

void printSystemHeader(const boost::filesystem::path &outputPath, std::string filename){

    // the / operator does concatenation
	std::ofstream systems((outputPath/filename).string());
    
    systems << "#\t#\t#\t#\tMsol\tMsol\tAU\t#\t#\t#\t#\t#\t#\t#\t#\t#\tyr^-1\tyr^-1\tkms^-1\tkms^-1\tkms^-1\tkms^-1\t#\t#\t#\t#\t#\t#\t#\t#\t#\t#" << std::endl;
    systems << "int\tint\tfloat\tint\tfloat\tfloat\tfloat\tfloat\tfloat\tfloat\tfloat\tfloat\tfloat\tfloat\tfloat\tfloat\tfloat\tfloat\tfloat\tfloat\tfloat\tfloat\tfloat\tfloat\tfloat\tfloat\tfloat\tbool\tbool\tint\tint\tbool" << std::endl;
    systems << "ID\tSEED\tweight\tsamplingPhase\tmass1\tmass2\tseparation\teccentricity\trk1\ttheta1\tphi1\tmeanAnomaly1\trk2\ttheta2\tphi2\tmeanAnomaly2\tomega1\tomega2\tsigma_kick_CCSN_NS\tsigma_kick_CCSN_BH\tsigma_kick_ECSN\tsigma_kick_USSN\tLBV_multiplier\tWR_multiplier\tCE_Alpha\tMetallicity1\tMetallicity2\tdisbound\tstellar_merger\tS1type_final\tS2type_final\terrorFlag" << std::endl;

    systems.close();

}

void printSystemParameters(const boost::filesystem::path &outputPath, std::string filename, const BinaryStar &currentBinary){

    // Open output stream
    std::ofstream systems((outputPath/filename).string(),std::ios_base::app);

	systems << 	currentBinary.m_ID << TAB << 
				currentBinary.m_randomSeed << TAB << 
				currentBinary.weight << TAB << 
				currentBinary.samplingphase << TAB << 
				currentBinary.star1.m_MZAMS << TAB << 
				currentBinary.star2.m_MZAMS << TAB << 
				currentBinary.m_SemiMajorAxisInitial << TAB << 
				currentBinary.m_EccentricityInitial << TAB << 
				currentBinary.star1.m_supernovaKickVelocityMagnitudeRandomNumber << TAB << 
				currentBinary.star1.m_supernovaTheta << TAB << 
				currentBinary.star1.m_supernovaPhi << TAB <<
            	currentBinary.star1.m_MeanAnomaly << TAB << 
            	currentBinary.star2.m_supernovaKickVelocityMagnitudeRandomNumber << TAB << 
            	currentBinary.star2.m_supernovaTheta << TAB << 
            	currentBinary.star2.m_supernovaPhi << TAB <<
            	currentBinary.star2.m_MeanAnomaly << TAB << 
            	currentBinary.star1.m_omegaZAMS << TAB << 
            	currentBinary.star2.m_omegaZAMS << TAB << 
            	currentBinary.m_kickVelocityDistributionSigmaCCSN_NS << TAB <<
            	currentBinary.m_kickVelocityDistributionSigmaCCSN_BH << TAB << 
            	currentBinary.m_kickVelocityDistributionSigmaForECSN <<	TAB << 
            	currentBinary.m_kickVelocityDistributionSigmaForUSSN << TAB << 
            	currentBinary.m_luminousBlueVariableFactor << TAB << 
            	currentBinary.m_wolfRayetFactor << TAB << 
            	currentBinary.m_commonEnvelopeAlpha << TAB << 
            	currentBinary.star1.m_Metallicity << TAB << 
            	currentBinary.star2.m_Metallicity << TAB << 
            	currentBinary.binary_disbound << TAB << 
            	currentBinary.m_stellarMerger << TAB << 
            	currentBinary.star1.m_stellarType << TAB << 
            	currentBinary.star2.m_stellarType << TAB << 
            	currentBinary.m_errorFlag << std::endl; 

    // Close output stream
    systems.close();

}

void printSupernovaeHeader(const boost::filesystem::path &outputPath, std::string filename){

    // the / operator does concatenation
	std::ofstream supernovae((outputPath/filename).string());

    supernovae << 	"#\tkms^-1\tkms^-1\t#\tkms^-1\t#\t#\t#\t#\t#\t#\t#\t#\t#\t#\tMSol\tMSol\tMSol\tMSol\tMSol\tMsol\t#\t#\t#\t#\t#\tMSol\t#\t#\t#\t#\tMyr\t#\t#\tRsun\tRsun\tkms^-1\t#\t#\t#" <<std::endl;

    supernovae << 	"int\tfloat\tfloat\tfloat\tfloat\tfloat\tfloat\tfloat\tfloat\tbool\tbool\tbool\tbool\tbool\tbool\tfloat\tfloat\tfloat\tfloat\tfloat\tfloat\tbool\tint\tint\tint\tint\tfloat\tfloat\tbool\tbool\tbool\tfloat\tfloat\tfloat\tfloat\tfloat\tfloat\tfloat\tfloat\tfloat" << std::endl;

    supernovae << 	"randomSeed\tdrawnKickVelocity\tkickVelocity\tfallback\tvRel\tuK\tpsi\ttheta\tphi\tflagECSN\tflagSN\tflagUSSN\tflagPISN\tflagPPISN\tSurvived\tMSNZAMS\tMCompanionZAMS\tMassStarSN\tMassStarCompanion\tMassCOCoreSN\tMrem\texperiencedRLOF\tRemnantType\twhichStar\tpreviousStellarTypeSN\tpreviousStellarTypeCompanion\tMassCoreSN\tMetallicity\tcommonEnvelopeOccuredAtLeastOnce\tm_stableRLOFafterCEE\tflagRLOFontoaNS    \ttime\teccentricityBefore\teccentricityAfter\tseparationBefore\tseparationAfter\tsystemicVelocity\tflagHrichSN\tflagHpoorSN\trunawayFlag" << std::endl;

    supernovae.close();
}

//alternative way of formatting tabs 
void printCommonEnvelopeHeader(const boost::filesystem::path &outputPath, std::string filename){

    std::ofstream commonEnvelope((outputPath/filename).string());

    commonEnvelope 	<< "#" << TAB << "Myrs" << TAB << "#" << TAB << "#" << TAB <<"#" << TAB << "ergs" << TAB << "ergs" << TAB <<
								"#" << TAB << "#" << TAB << "#" << TAB << "Rsol" << TAB << "Rsol" << TAB << 
								 "Rsol" << TAB << "Rsol" << TAB << "Rsol" << TAB << "Rsol"<< TAB <<
								 "Msol" << TAB << "Msol" << TAB << "Msol" << TAB << "Msol" << TAB << "Rsol" << TAB << "Rsol" << TAB << "#" << TAB << "#" << TAB <<
								 "#" << TAB << "#" << TAB << "#" << TAB << "#" << TAB << "#" << TAB <<
								 "ergs" << TAB << "ergs" << TAB << "ergs" << TAB << "ergs" << TAB << "ergs" << TAB <<
								 "Msol" << TAB << "Msol" << TAB << "Msol" << TAB << "Msol" << TAB << "Rsol" << TAB << "Rsol" << TAB << "#" << TAB << "#" << TAB <<
								 "#" << TAB << "#" << TAB << "#" << TAB << "#" << TAB << "#" << TAB <<
								 "ergs" << TAB << "ergs" << TAB << "ergs" << TAB << "ergs" << TAB << "ergs" << TAB <<
								 "#" << TAB << "#" << TAB << "#" << TAB << "#" << TAB << "#" << TAB <<  
 								 "#" << TAB << "Lsol" << TAB << "Tsol" << TAB << "Myrs" << TAB << "Myrs" << TAB << "Myrs" << TAB << 
 								 "#" << TAB << "Lsol" << TAB << "Tsol" << TAB << "Myrs" << TAB << "Myrs" << TAB << "Myrs" << TAB << 
								 "#" << TAB << "#" << TAB << "yrs" << TAB << "yrs" << TAB << "yrs" << TAB << "yrs"  << TAB <<
								 "#" << TAB << "#" << TAB << "#"
								 << std::endl;
        
     commonEnvelope  << "int" << TAB << "float" << TAB << "float" << TAB << "float" << TAB << "float" << TAB << "float" << TAB << "float" << TAB << 
                                "float" << TAB << "float" << TAB << "float" << TAB << "float" << TAB << "float" << TAB << 
                                "float" << TAB << "float" << TAB << "float" << TAB << "float" << TAB << 
                                "float" << TAB << "float" << TAB << "float" << TAB << "float" << TAB << "float" << TAB << "float" << TAB << "int" << TAB << "int" << TAB << 
                                "float" << TAB << "float" << TAB << "float" << TAB << "float" << TAB << "float" << TAB << 
                                "float" << TAB << "float" << TAB << "float" << TAB << "float" << TAB << "float" << TAB << 
                                "float" << TAB << "float" << TAB << "float" << TAB << "float" << TAB << "float" << TAB << "float" << TAB << "int" << TAB << "int" << TAB << 
                                "float" << TAB << "float" << TAB << "float" << TAB << "float" << TAB << "float" << TAB << 
                                "float" << TAB << "float" << TAB << "float" << TAB << "float" << TAB << "float" <<  TAB << 
                                "float" << TAB << "bool" << TAB << "bool" << TAB << "int" << TAB <<  "bool" <<  TAB << 
								"float" << TAB << "float" << TAB << "float" << TAB << "float" << TAB << "float" << TAB << "float" << TAB <<
								"float" << TAB << "float" << TAB << "float" << TAB << "float" << TAB << "float" << TAB << "float" << TAB <<
								"float" << TAB << "float" << TAB << "float" << TAB << "float" << TAB << "float" << TAB << "float" << TAB <<
								"bool" << TAB << "bool" << TAB << "bool"
								<< std::endl;
        
    commonEnvelope 	<< "randomSeed" << TAB << "time" << TAB << "alphaCE" << TAB << "lambda1" << TAB << "lambda2" << TAB << "bindingEnergy1" << TAB << "bindingEnergy2" << TAB << 
								"m_EccentricityRLOF" << TAB << "EccentricityPreCEE" << TAB << "EccentricityPostCEE" << TAB << "SemiMajorAxisPreCEE" << TAB << "SemiMajorAxisPostCEE" << TAB << 
								"rocheLobe1to2PreCEE" << TAB << "rocheLobe1to2PostCEE" << TAB << "rocheLobe2to1PreCEE" << TAB << "rocheLobe2to1PostCEE" << TAB << 
								"mass1ZAMS" << TAB << "mass1" << TAB << "massEnvelope1" << TAB << "massCore1" << TAB << 	"radius1" << TAB << "radius1AfterStripping" << TAB << "stellarType1" << TAB << "finalStellarType1" << TAB << 
								"lambdaFixed1" << TAB << "lambdaNanjing1" << TAB << "lambdaLoveridge1" << TAB << "lambdaLoveridgeWinds1" << TAB << "lambdaKruckow1" << TAB << 
								"bEFixed1" << TAB << "bENanjing1" << TAB << "bELoveridge1" << TAB << "bELoveridgeWinds1" << TAB << "bEKruckow1" << TAB << 
								"mass2ZAMS" << TAB << "mass2" << TAB << "massEnvelope2" << TAB << "massCore2" << TAB << "radius2" << TAB << "radius2AfterStripping" << TAB << "stellarType2" << TAB << "finalStellarType2" << TAB << 
								"lambdaFixed2" << TAB << "lambdaNanjing2" << TAB << "lambdaLoveridge2" << TAB << "lambdaLoveridgeWinds2" << TAB << "lambdaKruckow2" << TAB << 
								"bEFixed2" << TAB << "bENanjing2" << TAB << "bELoveridge2" << TAB << "bELoveridgeWinds2" << TAB << "bEKruckow2" <<  TAB << 
								"massTransferTrackerHistory" << TAB << "stellarMerger" << TAB << "optimisticCommonEnvelopeFlag" << TAB << "counterCEE" <<  TAB << "doubleCoreCommonEnvelopeFlag" <<  TAB << 
								"flagRLOF1" << TAB << "luminosity1" << TAB << "Teff1" << TAB << "tauDynamical1" << TAB << "tauThermal1" << TAB << "tauNuclear1" << TAB <<
								"flagRLOF2" << TAB << "luminosity2" << TAB << "Teff2" << TAB << "tauDynamical2" << TAB << "tauThermal2" << TAB << "tauNuclear2" << TAB <<
								"zetaStarCompare" <<  TAB << "zetaRLOFanalytic" << TAB << "tauSync" << TAB << "tauCirc" << TAB << "tauRexp1" << TAB << "tauRexp2" << TAB <<
								"immediateRLOFAfterCEE" << TAB << "simultaneousRLOFleadingToCEEFlag" << TAB << "MSAccretorDuringCEEFlag"
								<< std::endl;
    commonEnvelope.close();
}




// void printCommonEnvelope(const boost::filesystem::path &outputPath, std::string filename){
    // Difficult since called from BinaryStar::CommonEnvelopeEvent so how can you pass through parameters/binary?
//     // std::ios_base::app appends to the end of the file after the headers.
//     //     std::ofstream commonEnvelope((options.outputPath/"commonEnvelopes.txt").string(),std::ios_base::app);
//     std::ofstream commonEnvelope((outputPath/filename).string(), std::ios_base::app);

//     // commonEnvelope  << m_randomSeed << TAB << m_time << TAB << alphaCE << TAB << lambda1 
//     //                 << TAB << lambda2 << TAB << star1.m_bindingEnergyAtCommonEnvelope 
//     //                 << TAB << star2.m_bindingEnergyAtCommonEnvelope << TAB <<  m_EccentricityPreCEE 
//     //                 << TAB << m_EccentricityPostCEE << TAB << m_SemiMajorAxisPreCEE 
//     //                 << TAB << m_SemiMajorAxisPostCEE << TAB <<  m_rocheLobe1to2PreCEE 
//     //                 << TAB << m_rocheLobe1to2PostCEE << TAB << m_rocheLobe2to1PreCEE 
//     //                 << TAB << m_rocheLobe2to1PostCEE << TAB <<  star1.m_MZAMS 
//     //                 << TAB << m1 << TAB << menv1 << TAB << mfin1 << TAB << radius1 
//     //                 << TAB << radius1AfterStripping << TAB << type1 << TAB << finalType1 
//     //                 << TAB << options.commonEnvelopeLambda << TAB << Nanjing1 
//     //                 << TAB << Loveridge1 << TAB << Loveridge1Winds << TAB << Kruckow1 
//     //                 << TAB <<  bindingEnergyFixed1 << TAB << bindingEnergyNanjing1 
//     //                 << TAB << bindingEnergyLoveridge1 << TAB << bindingEnergyLoveridge1Winds 
//     //                 << TAB << bindingEnergyKruckow1 << TAB << star2.m_MZAMS  
//     //                 << TAB << m2 << TAB << menv2 << TAB << mfin2 << TAB << radius2 
//     //                 << TAB << radius2AfterStripping << TAB << type2 << TAB << finalType2 
//     //                 << TAB << options.commonEnvelopeLambda << TAB << Nanjing2 
//     //                 << TAB << Loveridge2 << TAB << Loveridge2Winds << TAB << Kruckow2 
//     //                 << TAB <<  bindingEnergyFixed2 << TAB << bindingEnergyNanjing2 
//     //                 << TAB << bindingEnergyLoveridge2 << TAB << bindingEnergyLoveridge2Winds 
//     //                 << TAB << bindingEnergyKruckow2 << TAB <<  m_massTransferTrackerHistory 
//     //                 << TAB << m_stellarMerger << TAB << m_optimisticCommonEnvelopeFlag 
//     //                 << TAB << m_counterCEE << std::endl;
//     // commonEnvelope.close();
// }

void printFormationHeader(const boost::filesystem::path &outputPath, std::string filename){

    std::ofstream formationHistory;

    formationHistory.open((outputPath/filename).string());

    formationHistory << "#\t#\t#\t#\t#\t#\t#\t#\t#\t#\t#\t#\t#\t#\t#\t#\t#\t#\t#\t#\t#\t#\t#\t#\t#\t#\t#\t#\t#\t#\t#\t#\t#\t#\t#\t#" << std::endl;

    formationHistory << "int\tint\tint\tint\tint\tint\tint\tint\tint\tint\tint\tint\tint\tint\tint\tint\tint\tint\tint\tint\tint\tint\tint\tint\tint\tint\tint\tint\tint\tint\tint\tint\tint\tint\tbool\tbool" << std::endl;

    formationHistory << "m_randomSeed\teventCounter\tmt_primary_ep1\tmt_primary_ep1_K1\tmt_primary_ep1_K2\tmt_primary_ep2\t mt_primary_ep2_K1\tmt_primary_ep2_K2\tmt_primary_ep3\t mt_primary_ep3_K1\tmt_primary_ep3_K2\tmt_secondary_ep1\tmt_secondary_ep1_K1\t mt_secondary_ep1_K2\tmt_secondary_ep2\tmt_secondary_ep2_K1\tmt_secondary_ep2_K2\tmt_secondary_ep3\tmt_secondary_ep3_K1\tmt_secondary_ep3_K2\tSN_primary_type_1\tSN_primary_type_2\tSN_primary_type_3\tSN_secondary_type_1\tSN_secondary_type_2\t SN_secondary_type_3\tCEE\tCEE_instigator\tCEE_failed\tCEE_failed_instigator\tCEE_wet\tCEE_wet_instigator\tstellar_type_K1\tstellar_type_K2\tmerged_in_Hubble_time\tbinary_disbound" << std::endl;

    formationHistory.close();
}

void printFormationChannel(const boost::filesystem::path &outputPath, std::string filename, const BinaryStar &currentBinary){
    
    // std::ios_base::app appends to the end of the file after the headers.
    std::ofstream formationHistory((outputPath/filename).string(), std::ios_base::app);

    formationHistory <<  currentBinary.m_randomSeed  << TAB<<  currentBinary.eventCounter << TAB << currentBinary.mt_primary_ep1
        << TAB <<   currentBinary.mt_primary_ep1_K1   << TAB <<   currentBinary.mt_primary_ep1_K2
        << TAB <<   currentBinary.mt_primary_ep2      << TAB <<   currentBinary.mt_primary_ep2_K1
        << TAB <<   currentBinary.mt_primary_ep2_K2   << TAB <<   currentBinary.mt_primary_ep3
        << TAB <<   currentBinary.mt_primary_ep3_K1   << TAB <<   currentBinary.mt_primary_ep3_K2 
        << TAB <<   currentBinary.mt_secondary_ep1    << TAB <<   currentBinary.mt_secondary_ep1_K1
        << TAB <<   currentBinary.mt_secondary_ep1_K2 << TAB <<   currentBinary.mt_secondary_ep2
        << TAB <<   currentBinary.mt_secondary_ep2_K1 << TAB <<   currentBinary.mt_secondary_ep2_K2 
        << TAB <<   currentBinary.mt_secondary_ep3    << TAB <<   currentBinary.mt_secondary_ep3_K1
        << TAB <<   currentBinary.mt_secondary_ep3_K2 << TAB <<   currentBinary.SN_primary_type_1
        << TAB <<   currentBinary.SN_primary_type_2   << TAB <<   currentBinary.SN_primary_type_3
        << TAB <<   currentBinary.SN_secondary_type_1 << TAB <<   currentBinary.SN_secondary_type_2
        << TAB <<   currentBinary.SN_secondary_type_3 << TAB <<   currentBinary.CEE 
        << TAB <<   currentBinary.CEE_instigator      << TAB <<   currentBinary.CEE_failed
        << TAB<<    currentBinary.CEE_failed_instigator<< TAB <<   currentBinary.CEE_wet
        << TAB <<   currentBinary.CEE_wet_instigator   << TAB <<   currentBinary.stellar_type_K1
        << TAB <<   currentBinary.stellar_type_K2      << TAB <<   currentBinary.m_Merged 
		<< TAB <<   currentBinary.binary_disbound   << std::endl;

    formationHistory.close();
}

void printXRayBinaryHeader(const boost::filesystem::path &outputPath, std::string filename){
    
        std::ofstream XRayBinary((outputPath/filename).string());
        XRayBinary << "#\tAU\tMSol\tMSol\tMSol\tMSol\tAU\t#\t#\t#\tLSol\tLSol\tergss^-1\t#\t#\tMyr\tMyr"<< std::endl;
        XRayBinary << "int\tfloat\tfloat\tfloat\tfloat\tfloat\tfloat\tfloat\tint\tint\tfloat\tfloat\tfloat\tbool\tbool\tfloat\tfloat"<< std::endl;
        XRayBinary << "m_randomSeed\tseparationInitial\tM1ZAMS\tM2ZAMS\tMass1\tMass2\tseparation\teccentricity\tType1\tType2\tLum1\tLum2\tLx\tflagRLOF1\tflagRLOF2\tdt\ttime"<< std::endl;
        XRayBinary.close();
}


void printXRayBinaryParameters(const boost::filesystem::path &outputPath, std::string filename, const BinaryStar &currentBinary){
        std::ofstream XRayBinary((outputPath/filename).string(),std::ios_base::app);
        XRayBinary << currentBinary.m_randomSeed << TAB << currentBinary.m_SemiMajorAxisInitial << TAB <<
                    currentBinary.star1.m_MZAMS << TAB << currentBinary.star2.m_MZAMS << TAB <<
                    currentBinary.star1.m_Mass << TAB << currentBinary.star2.m_Mass << TAB <<
                    currentBinary.m_SemiMajorAxisPrime*AUToRsol << TAB <<
                    currentBinary.m_Eccentricity << TAB << currentBinary.star1.m_stellarType << TAB << currentBinary.star2.m_stellarType << TAB <<
                    currentBinary.star1.m_Luminosity << TAB << currentBinary.star2.m_Luminosity << TAB << currentBinary.m_XRayLuminosity << TAB << currentBinary.star1.flagRLOF
                    << TAB << currentBinary.star2.flagRLOF << TAB << std::setprecision(10) << currentBinary.m_dt << TAB << currentBinary.m_time << std::endl;
        XRayBinary.close();
}

void printCHEvolutionHeader(const boost::filesystem::path &outputPath, std::string filename){

	std::ofstream CHEvolution((outputPath/filename).string());
    CHEvolution << "#\t#\tMsol\tMsol\tRsol\t#\t#\t#\t!\t!" << std::endl;
    CHEvolution << "int\tint\tfloat\tfloat\tfloat\tbool\tbool\tbool\tfloat\tfloat" << std::endl; 
    CHEvolution << "ID\tSEED\tmass1\tmass2\tseparation\tflagCHE1\tflagCHE2\tflagCHEBinary\tperiod1\tperiod2" << std::endl;            
    CHEvolution.close();
}

void printCHEvolutionParameters(const boost::filesystem::path &outputPath, std::string filename, const BinaryStar &currentBinary){
    std::ofstream CHEvolution((outputPath/filename).string(),std::ios_base::app);
	CHEvolution << currentBinary.m_ID << TAB << currentBinary.m_randomSeed << TAB << currentBinary.star1.m_Mass 
			<< TAB << currentBinary.star2.m_Mass << TAB << currentBinary.m_SemiMajorAxisPrime*AUToRsol<< TAB << currentBinary.star1.flagCHEvolution<< TAB << currentBinary.star2.flagCHEvolution << TAB << currentBinary.m_CHEvolutionFlag << TAB << currentBinary.star1.m_Period << TAB << currentBinary.star2.m_Period  << std::endl; 
    CHEvolution.close();
}


void printpulsarEvolutionHeader(const boost::filesystem::path &outputPath, std::string filename){

        std::ofstream pulsarEvolution((outputPath/filename).string());
    pulsarEvolution << "#\t#\t#\tMsol\tMsol\t#\t#\tRsol\t#\tTesla\tTesla\trad/s^1\trad/s^1\trad/s^2\trad/s^2\tMyrs\tMyrs" << std::endl;
    pulsarEvolution << "int\tint\tbool\tfloat\tfloat\tint\tint\tfloat\tbool\tfloat\tfloat\tfloat\tfloat\tfloat\tfloat\tfloat\tfloat" << std::endl;
    pulsarEvolution << "ID\tSEED\tdisbound\tmass1\tmass2\ttype1\ttype2\tseparation\tmassTransferTracker\tmagneticField1\tmagneticField2\tspinFrequency1\tspinFrequency2\tspinDownRate1\tspinDownRate2\ttime\tdt" << std::endl;
    pulsarEvolution.close();
}

void printpulsarEvolutionParameters(const boost::filesystem::path &outputPath, std::string filename, const BinaryStar &currentBinary){
    std::ofstream pulsarEvolution((outputPath/filename).string(),std::ios_base::app);
        pulsarEvolution << currentBinary.m_ID << TAB << currentBinary.m_randomSeed << TAB <<  currentBinary.binary_disbound << TAB <<  currentBinary.star1.m_Mass
                        << TAB << currentBinary.star2.m_Mass << TAB << currentBinary.star1.m_stellarType << TAB <<  currentBinary.star2.m_stellarType << TAB <<  currentBinary.m_SemiMajorAxisPrime*AUToRsol<< TAB << currentBinary.m_massTransferTrackerHistory << TAB << currentBinary.star1.m_pulsarMagneticField<< TAB << currentBinary.star2.m_pulsarMagneticField << TAB << currentBinary.star1.m_pulsarSpinFrequency << TAB << currentBinary.star2.m_pulsarSpinFrequency<< TAB << currentBinary.star1.m_pulsarSpinDownRate << TAB << currentBinary.star2.m_pulsarSpinDownRate << TAB << currentBinary.m_time << TAB << std::setprecision(10) << currentBinary.m_dt << std::endl;
    pulsarEvolution.close();
}





void printRLOFHeader(const boost::filesystem::path &outputPath, std::string filename){

    std::ofstream RLOF((outputPath/filename).string());
    RLOF << "#\t#\tMsun\tMsun\tRsun\tRsun\ttype1\ttype2\tRsun\t#\tMyr\t#\t#\t#\tMsun\tMsun\tRsun\tRsun\t#\t#\tRsun\t#\tMyr\t#\t#\t#\t#\t#\t#\t#\t#\t#\t#\t#\t#\t#\t#\t#\t#\t#\t#" << std::endl;
    RLOF << "int\tint\tfloat\tfloat\tfloat\tfloat\tint\tint\tfloat\tint\tfloat\tbool\tbool\tbool\tfloat\tfloat\tfloat\tfloat\tint\tint\tfloat\tint\tfloat\tbool\tbool\tfloat\tfloat\tfloat\tfloat\tfloat\tfloat\tfloat\tfloat\tfloat\tfloat\tfloat\tfloat\tfloat\tfloat\tfloat\tfloat" << std::endl; 
    RLOF <<    "ID\trandomSeed\tmass1\tmass2\tradius1\tradius2\ttype1\ttype2\tseparation\teventCounter\ttime\tflagRLOF1\tflagRLOF2\tflagCEE\tmass1Prev\tmass2Prev\tradius1Prev\tradius2Prev\ttype1Prev\ttype2Prev\tseparationPrev\teventCounterPrev\ttimePrev\tflagRLOF1Prev\tflagRLOF2Prev\tzetaThermal1\tzetaNuclear1\tzetaSoberman1\tzetaSobermanHelium1\tzetaHurley1\tzetaHurleyHelium1\tzetaSimple1\tzetaThermal2\tzetaNuclear2\tzetaSoberman2\tzetaSobermanHelium2\tzetaHurley2\tzetaHurleyHelium2\tzetaSimple2\tzetaRLOFanalytical\tzetaRLOFnumerical" <<std::endl;
    RLOF.close(); 



}


void printBeBinariesHeader(const boost::filesystem::path &outputPath, std::string filename){
    
    // std::ios_base::app appends to the end of the file.
    std::ofstream BeB((outputPath/filename).string());
    
    BeB << "#\t#\tMyrs\tMyrs\tMsol\tTODO\tMsol\tLsol\tTsol\tRsol\t#\tRsol\t#" << std::endl;
    BeB << "int\tint\tfloat\tfloat\tfloat\tfloat\tfloat\tfloat\tfloat\tfloat\tbool\tfloat\tfloat" << std::endl; 
    BeB << "ID\trandomSeed\tdt\ttotalTime\tmassNS\tspinNS\tcompanionMass\tcompanionL\tcompanionT\tcompanionR\tcompanionAccreted\tseparation\teccentricity" << std::endl; 


    BeB.close();  

}




void printRejectedSamplesHeader(const boost::filesystem::path &outputPath, std::string filename){
    // printing header of rejected samples output by COMPAS 
    // std::ios_base::app appends to the end of the file.
    std::ofstream RejectedSamples((outputPath/filename).string());
    
    RejectedSamples << "#" << std::endl;
    RejectedSamples << "int" << std::endl; 
    RejectedSamples << "randomSeed" << std::endl; 

    RejectedSamples.close();  
}


void printRejectedSamplesParameters(const boost::filesystem::path &outputPath, std::string filename,  unsigned long randomSeed){
	    // printing Seeds of rejected samples by COMPAS 
    std::ofstream RejectedSamples((outputPath/filename).string(),std::ios_base::app);
        RejectedSamples <<  randomSeed << std::endl;
    RejectedSamples.close();
}




std::string detailedOutputFilename(int i){
    /*
    */
    // Detailed output
    std::string Number;                          // string which will contain the result
    std::ostringstream convert;                  // stream used for the conversion
    convert << i;                           // insert the textual representation of the i-th binary
    Number = convert.str();

    std::string dataOutputFile = "dataOutput_";
    dataOutputFile  += Number;
    dataOutputFile  += ".dat";

    return dataOutputFile;
}

void printDetailedOutputHeader(const boost::filesystem::path &outputPath, int i){
    /*
    */

    std::string filename = detailedOutputFilename(i);

    std::ofstream data((outputPath/filename).string());

    data << "%adim" << TAB << "Myr" << TAB << "Myr" << TAB << "Rsol" << TAB << "adim" << TAB << 
    "Msol" << TAB << "Msol" << TAB << "Msol" << TAB << "Msol" << TAB << "Msol" << TAB << "Msol" << TAB << 
    "Msol" << TAB << "Msol" << TAB << "Msol" << TAB << "Msol" << TAB << "Msol" << TAB << "Msol" << TAB << 
    "Rsol" << TAB << "Rsol" << TAB << "adim" << TAB << "adim" << TAB << "adim" << TAB << "adim" << TAB << 
    "yr-1" << TAB << "yr-1" << TAB << "yr-1" << TAB << "yr-1" << TAB << "adim" << TAB << "adim" << TAB << 
    "Myr"  << TAB << "Myr"  << TAB << "Lsol" << TAB << "Lsol" << TAB << "Tsol" << TAB << "Tsol" << TAB << 
    "Msol*AU^2*yr-1" << TAB << "Msol*AU^2*yr-1" << TAB << "Myr" << TAB << "Myr" << TAB << "Myr" << TAB << 
    "Myr"  << TAB << "Myr" << TAB << "Myr" << TAB << "adim" << TAB << "adim" << TAB << "adim" << TAB << 
    "adim" << TAB << "adim" << TAB << "adim" << TAB << "adim" << TAB << "adim" << TAB << "adim" << TAB << 
    "adim" << TAB << "adim" << TAB << "adim" << TAB << "adim" << TAB << "adim" << TAB << "Msol" << TAB << 
    "Msol" << TAB << "Msol" << TAB << "Msol" << TAB << "adim" << TAB << "Msol*AU^2*yr-1"<< TAB << 
    "Msol*AU^2*yr-2" << TAB << "adim"<< TAB << "adim"<< TAB << "adim"<< TAB << "adim"<< TAB << 
    "adim"<< TAB << "adim"<< TAB << "adim" << TAB << "Tesla" << TAB << "Tesla"<< TAB << "rad/s"<< TAB << 
    "rad/s"<< TAB << "rad/s^2" << TAB << "rad/s^2" << TAB << "yr" << TAB << "yr" << std::endl;

    data << "%SEED" << TAB << "dt" << TAB << "time" << TAB << "separation" << TAB << "eccentricity" << TAB << 
    "mass1_0" << TAB << "mass2_0" << TAB << "mass1" << TAB << "mass2" << TAB << "massEnv1" << TAB << 
    "massEnv2" << TAB << "massCore1" << TAB << "massCore2" << TAB << "massHeCore1" << TAB << "massHeCore2" << TAB << 
    "massCOCore1" << TAB << "massCOCore2" << TAB << "Radius1" << TAB << "Radius2" << TAB << "RocheLobe1/a" << TAB << 
    "RocheLobe2/a" << TAB << "Radius1/RL1" << TAB << "Radius2/RL2" << TAB << "omega1" << TAB << "omega2" << TAB << 
    "omegaBreak1" << TAB << "omegaBreak2" << TAB << "S1type" << TAB << "S2type" << TAB << "Age1" << TAB << 
    "Age2" << TAB << "Lum1" << TAB << "Lum2" << TAB << "Teff1" << TAB << "Teff2" << TAB << "AngMom1" << TAB << 
    "AngMom2" << TAB << "tauDynamical1" << TAB << "tauDynamical2" << TAB << "tauThermal1" << TAB << 
    "tauThermal2" << TAB << "tauNuclear1" << TAB << "tauNuclear2" << TAB << "ZThermal1" << TAB << 
    "ZThermal2" << TAB << "ZNuclear1" << TAB << "ZNuclear2" << TAB << "ZSPH1" << TAB << "ZSPH2" << TAB << 
    "ZSPHHe1" << TAB << "ZSPHHe2" << TAB << "ZHurley1" << TAB << "ZHurley2" << TAB << "ZHurleyHe1" << TAB <<
    "ZHurleyHe2" << TAB << "Zsimple1" << TAB << "Zsimple2" << TAB << "dmWinds1" << TAB << "dmWinds2" << TAB << 
    "dmMassTransfer1" << TAB << "dmMassTransfer2" << TAB << "MTtype" << TAB << "AngMomTotal"<< TAB << 
    "EnergyTotal" << TAB << "Nanjing1" << TAB << "Nanjing2" << TAB << "Loveridge1" << TAB << 
    "Loveridge2" << TAB << "Kruckow1Top" << TAB << "Kruckow2Top" << TAB << "Kruckow1Mid" << TAB << 
    "Kruckow2Mid" << TAB << "Kruckow1Bot" << TAB << "Kruckow2Bot" << TAB << "Metallicity1" << TAB << 
    "Metallicity2" << TAB << "massTransferTracker" << TAB <<  "pulsarMagneticField1" << TAB <<
     "pulsarMagneticField2" << TAB << "pulsarSpinFrequency1" << TAB << "pulsarSpinFrequency2" << TAB <<
     "pulsarSpinDownRate1" << TAB << "pulsarSpinDownRate2" << TAB << "tauRadial1" << TAB << "tauRadial2" <<std::endl ;

    data.close();

}

void printDetailedOutput(const boost::filesystem::path &outputPath, int i, const BinaryStar &currentBinary){
    /*
    */

    std::string filename = detailedOutputFilename(i);

    std::ofstream data((outputPath/filename).string(),std::ios_base::app);

    data << std::setprecision(10) << currentBinary.m_randomSeed << TAB << 
                    0.0 << TAB << 
                    currentBinary.m_time << TAB << 
                    currentBinary.m_SemiMajorAxisPrime*AUToRsol << TAB << 
                    currentBinary.m_EccentricityPrime << TAB << 
                    currentBinary.star1.m_Mass0 << TAB << 
                    currentBinary.star2.m_Mass0 << TAB << 
                    currentBinary.star1.m_Mass << TAB << 
                    currentBinary.star2.m_Mass << TAB << 
                    currentBinary.star1.m_envMass << TAB << 
                    currentBinary.star2.m_envMass << TAB << 
                    currentBinary.star1.m_coreMass << TAB << 
                    currentBinary.star2.m_coreMass << TAB << 
                    currentBinary.star1.m_HeCoreMass << TAB << 
                    currentBinary.star2.m_HeCoreMass << TAB << 
                    currentBinary.star1.m_COCoreMass << TAB << 
                    currentBinary.star2.m_COCoreMass << TAB << 
                    currentBinary.star1.m_Radius << TAB << 
                    currentBinary.star2.m_Radius << TAB << 
                    rocheLobeRadius(currentBinary.star1.m_Mass, currentBinary.star2.m_Mass) << TAB << 
                    rocheLobeRadius(currentBinary.star2.m_Mass, currentBinary.star1.m_Mass) << TAB << 
                    currentBinary.rocheLobeTracker1 << TAB << 
                    currentBinary.rocheLobeTracker2 << TAB << 
                    currentBinary.star1.m_omega << TAB << 
                    currentBinary.star2.m_omega << TAB << 
                    currentBinary.star1.m_omegaBreak << TAB << 
                    currentBinary.star2.m_omegaBreak << TAB << 
                    currentBinary.star1.m_stellarType << TAB << 
                    currentBinary.star2.m_stellarType << TAB << 
                    currentBinary.star1.m_Age << TAB << 
                    currentBinary.star2.m_Age << TAB << 
                    currentBinary.star1.m_Luminosity << TAB << 
                    currentBinary.star2.m_Luminosity << TAB << 
                    currentBinary.star1.m_Temperature << TAB << 
                    currentBinary.star2.m_Temperature << TAB << 
                    currentBinary.star1.m_angularMomentum << TAB << 
                    currentBinary.star2.m_angularMomentum << TAB << 
                    currentBinary.star1.m_dynamicalTimescale << TAB <<
                    currentBinary.star2.m_dynamicalTimescale << TAB <<
                    currentBinary.star1.m_thermalTimescale << TAB <<
                    currentBinary.star2.m_thermalTimescale << TAB <<
                    currentBinary.star1.m_nuclearTimescale << TAB <<
                    currentBinary.star2.m_nuclearTimescale << TAB <<
                    currentBinary.star1.m_zetaThermal << TAB <<  
                    currentBinary.star2.m_zetaThermal << TAB << 
                    currentBinary.star1.m_zetaNuclear << TAB << 
                    currentBinary.star2.m_zetaNuclear << TAB << 
                    currentBinary.star1.m_zetaSoberman << TAB << 
                    currentBinary.star2.m_zetaSoberman << TAB << 
                    currentBinary.star1.m_zetaSobermanHelium << TAB << 
                    currentBinary.star2.m_zetaSobermanHelium << TAB << 
                    currentBinary.star1.m_zetaHurley << TAB <<
                    currentBinary.star2.m_zetaHurley << TAB <<
                    currentBinary.star1.m_zetaHurleyHelium << TAB <<
                    currentBinary.star2.m_zetaHurleyHelium << TAB <<
                    currentBinary.star1.m_zetaSimple << TAB <<
                    currentBinary.star2.m_zetaSimple << TAB <<
                    currentBinary.star1.m_massLossDiff << TAB << 
                    currentBinary.star2.m_massLossDiff << TAB << 
                    currentBinary.star1.m_MassTransferDiff << TAB << 
                    currentBinary.star2.m_MassTransferDiff << TAB << 
                    currentBinary.m_massTransferType << TAB << 
                    currentBinary.m_TotalAngularMomentumPrime << TAB << 
                    currentBinary.m_TotalEnergyPrime << TAB << 
                    currentBinary.star1.m_NanjingLambda << TAB << 
                    currentBinary.star2.m_NanjingLambda << TAB << 
                    currentBinary.star1.m_LoveridgeLambda << TAB << 
                    currentBinary.star2.m_LoveridgeLambda << TAB << 
                    currentBinary.star1.m_KruckowTopLambda << TAB << 
                    currentBinary.star2.m_KruckowTopLambda << TAB << 
                    currentBinary.star1.m_KruckowMidLambda << TAB << 
                    currentBinary.star2.m_KruckowMidLambda << TAB << 
                    currentBinary.star1.m_KruckowBotLambda << TAB << 
                    currentBinary.star2.m_KruckowBotLambda << TAB << 
                    currentBinary.star1.m_Metallicity << TAB << 
                    currentBinary.star2.m_Metallicity << TAB << 
                    currentBinary.m_massTransferTrackerHistory << TAB << 
                    currentBinary.star1.m_pulsarMagneticField << TAB << 
                    currentBinary.star2.m_pulsarMagneticField << TAB <<
                    currentBinary.star1.m_pulsarSpinFrequency << TAB <<
                    currentBinary.star2.m_pulsarSpinFrequency << TAB << 
                    currentBinary.star1.m_pulsarSpinDownRate << TAB <<
                    currentBinary.star2.m_pulsarSpinDownRate << TAB << 
                    currentBinary.star1.m_radialExpansionTimescale << TAB << 
                    currentBinary.star2.m_radialExpansionTimescale << std::endl;

    data.close();

}


// ALEJANDRO - 02/12/2016 - Testing Coen's printing functions for Kit's project
/*
void printDetailedXRayBinariesHeader(){

	            // Detailed output
            string Number;                          // string which will contain the result
            ostringstream convert;                  // stream used for the conversion
            convert << i;                           // insert the textual representation of the i-th binary
            Number = convert.str();

            string dataOutputFile = "dataOutput_";
            dataOutputFile  += Number;
            dataOutputFile  += ".dat";

            data.open((options.outputPath/dataOutputFile).string());

			
    std::ofstream XRayBinariesDetailedOutput;
    XRayBinariesDetailedOutput.open("XRayBinariesDetailedOutput.dat");
    XRayBinariesDetailedOutput << "randomSeed" << std::endl;
    XRayBinariesDetailedOutput.close();
}
*/


