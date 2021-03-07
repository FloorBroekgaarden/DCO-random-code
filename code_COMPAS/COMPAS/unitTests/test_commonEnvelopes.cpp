//
//  test_kicks.cpp
//  
//  Created by Jim Barrett on 12/01/2017.
//
//	unit tests for the functionality of supernova kicks
//

#define BOOST_TEST_DYN_LINK

#include "../testingUtils.h"
#include "../supernovaFunctions.h"

#include <vector>
#include <cmath>
#include <boost/test/unit_test.hpp>

 
BOOST_AUTO_TEST_SUITE(CommonEnvelopes)

BOOST_AUTO_TEST_CASE(TestDoubleCommonEnvelope)
{
    //JWB - Draws kicks from our Maxwellian distribution and checks the RMS matches the chosen one
    //Draws 10000 samples per distribution for distributions with 10 randomly chosen sigmas
    
    gsl_rng* r = getGslRng();
    
    programOptions* options = getDefaultCompasOptions();
    
    options->detailedOutput = true;
    
    options->kickVelocityDistribution = "ZERO";
    
    for (int i=0; i<100; i++)
    {
        BinaryStar testBinary = generateRandomBinary();
    
        testBinary.star1.m_Mass = 10.01;
        testBinary.star2.m_Mass = 9.99;
        
        double a = 9. - (i*0.05);
        testBinary.m_SemiMajorAxisInitial = a;
        testBinary.setSemiMajorAxis(a);
        testBinary.m_SemiMajorAxisPrime = a;
        testBinary.m_SemiMajorAxisPrev = a;
        
        testBinary.star1.evolveOneTimestep(0.,false,*options);
        testBinary.star2.evolveOneTimestep(0.,false,*options);

        bool firstTimeFlag = true;
        
        double dt = 0., dtPrev = 0.;
        
        std::ofstream data("detailed"+std::to_string(a)+"MSol.txt");
        
        data << "%adim" << TAB << "Myr" << TAB << "Myr" << TAB << "Rsol" << TAB << "adim" << TAB << "Msol" << TAB << "Msol" << TAB << "Msol" << TAB << "Msol" << TAB << "Msol" << TAB << "Msol" << TAB << "Msol" << TAB << "Msol" << TAB << "Msol" << TAB << "Msol" << TAB << "Msol" << TAB << "Msol" << TAB << "Rsol" << TAB << "Rsol" << TAB << "adim" << TAB << "adim" << TAB << "adim" << TAB << "adim" << TAB << "yr-1" << TAB << "yr-1" << TAB << "yr-1" << TAB << "yr-1" << TAB << "adim" << TAB << "adim" << TAB << "Myr" << TAB << "Myr" << TAB << "Lsol" << TAB << "Lsol" << TAB << "Tsol" << TAB << "Tsol" << TAB << "Msol*AU^2*yr-1" << TAB << "Msol*AU^2*yr-1" << TAB << "Myr" << TAB << "Myr" << TAB << "Myr" << TAB << "Myr" << TAB << "Myr" << TAB << "Myr" << TAB << "adim" << TAB << "adim" << TAB << "Myr-1" << TAB << "Myr-1" << TAB << "Msol" << TAB << "Msol" << TAB << "Msol" << TAB << "Msol" << TAB << "adim" << TAB << "Msol*AU^2*yr-1"<< TAB << "Msol*AU^2*yr-2" << TAB << "adim"<< TAB << "adim"<< TAB << "adim"<< TAB << "adim"<< TAB << "adim"<< TAB << "adim"<< TAB << "adim" << TAB << "adim" << TAB << "adim" << std::endl;
            data << "%SEED" << TAB << "dt" << TAB << "time" << TAB << "separation" << TAB << "eccentricity" << TAB << "mass1_0" << TAB << "mass2_0" << TAB << "mass1" << TAB << "mass2" << TAB << "massEnv1" << TAB << "massEnv2" << TAB << "massCore1" << TAB << "massCore2" << TAB << "massHeCore1" << TAB << "massHeCore2" << TAB << "massCOCore1" << TAB << "massCOCore2" << TAB << "Radius1" << TAB << "Radius2" << TAB << "RocheLobe1/a" << TAB << "RocheLobe2/a" << TAB << "Radius1/RL1" << TAB << "Radius2/RL2" << TAB << "omega1" << TAB << "omega2" << TAB << "omegaBreak1" << TAB << "omegaBreak2" << TAB << "S1type" << TAB << "S2type" << TAB << "Age1" << TAB << "Age2" << TAB << "Lum1" << TAB << "Lum2" << TAB << "Teff1" << TAB << "Teff2" << TAB << "AngMom1" << TAB << "AngMom2" << TAB << "tauDynamical1" << TAB << "tauDynamical2" << TAB << "tauThermal1" << TAB << "tauThermal2" << TAB << "tauNuclear1" << TAB << "tauNuclear2" << TAB << "ZThermal1" << TAB << "ZThermal2" << TAB << "ZNuclear1" << TAB << "ZNuclear2" << TAB << "dmWinds1" << TAB << "dmWinds2" << TAB << "dmMassTransfer1" << TAB << "dmMassTransfer2" << TAB << "MTtype" << TAB << "AngMomTotal"<< TAB << "EnergyTotal" << TAB << "Nanjing1" << TAB << "Nanjing2" << TAB << "Loveridge1" << TAB << "Loveridge2" << TAB << "Kruckow1Top" << TAB << "Kruckow2Top" << TAB << "Kruckow1Mid" << TAB << "Kruckow2Mid" << TAB << "Kruckow1Bot" << TAB << "Kruckow2Bot" << TAB << "Metallicity1" << TAB << "Metallicity2" << TAB << "massTransferTracker" << std::endl;
            data << testBinary.m_randomSeed << TAB << 0.0 << TAB << testBinary.m_time << TAB << testBinary.m_SemiMajorAxisPrime*AUToRsol << TAB << testBinary.m_EccentricityPrime << TAB << testBinary.star1.m_Mass0 << TAB << testBinary.star2.m_Mass0 << TAB << testBinary.star1.m_Mass << TAB << testBinary.star2.m_Mass << TAB << testBinary.star1.m_envMass << TAB << testBinary.star2.m_envMass << TAB << testBinary.star1.m_coreMass << TAB << testBinary.star2.m_coreMass << TAB << testBinary.star1.m_HeCoreMass << TAB << testBinary.star2.m_HeCoreMass << TAB << testBinary.star1.m_COCoreMass << TAB << testBinary.star2.m_COCoreMass << TAB << testBinary.star1.m_Radius << TAB << testBinary.star2.m_Radius << TAB << rocheLobeRadius(testBinary.star1.m_Mass, testBinary.star2.m_Mass) << TAB << rocheLobeRadius(testBinary.star2.m_Mass, testBinary.star1.m_Mass) << TAB << testBinary.rocheLobeTracker1 << TAB << testBinary.rocheLobeTracker2 << TAB << testBinary.star1.m_omega << TAB << testBinary.star2.m_omega << TAB << testBinary.star1.m_omegaBreak << TAB << testBinary.star2.m_omegaBreak << TAB << testBinary.star1.m_stellarType << TAB << testBinary.star2.m_stellarType << TAB << testBinary.star1.m_Age << TAB << testBinary.star2.m_Age << TAB << testBinary.star1.m_Luminosity << TAB << testBinary.star2.m_Luminosity << TAB << testBinary.star1.m_Temperature << TAB << testBinary.star2.m_Temperature << TAB << testBinary.star1.m_angularMomentum << TAB << testBinary.star2.m_angularMomentum << TAB << dynamicalTimescale(testBinary.star1.m_Mass, testBinary.star1.m_Radius) << TAB << dynamicalTimescale(testBinary.star2.m_Mass, testBinary.star2.m_Radius) << TAB << thermalTimescale(testBinary.star1.m_Mass, testBinary.star1.m_envMass, testBinary.star1.m_Radius, testBinary.star1.m_Luminosity, testBinary.star1.m_stellarType) << TAB << thermalTimescale(testBinary.star2.m_Mass, testBinary.star2.m_envMass, testBinary.star2.m_Radius, testBinary.star2.m_Luminosity, testBinary.star2.m_stellarType) << TAB << nuclearTimescale(testBinary.star1.m_Mass, testBinary.star1.m_Luminosity) << TAB << nuclearTimescale(testBinary.star2.m_Mass, testBinary.star2.m_Luminosity) << TAB << determineConvergedMassStepZetaThermal(testBinary.star1, *options) << TAB << determineConvergedMassStepZetaThermal(testBinary.star2, *options) << TAB << determineConvergedTimestepZetaNuclear(testBinary.star1, *options) << TAB << determineConvergedTimestepZetaNuclear(testBinary.star2, *options) << TAB << testBinary.star1.m_massLossDiff << TAB << testBinary.star2.m_massLossDiff << TAB << testBinary.star1.m_MassTransferDiff << TAB << testBinary.star2.m_MassTransferDiff << TAB << testBinary.m_massTransferType << TAB << testBinary.m_TotalAngularMomentumPrime << TAB << testBinary.m_TotalEnergyPrime << TAB << testBinary.star1.m_NanjingLambda << TAB << testBinary.star2.m_NanjingLambda << TAB << testBinary.star1.m_LoveridgeLambda << TAB << testBinary.star2.m_LoveridgeLambda << TAB << testBinary.star1.m_KruckowTopLambda << TAB << testBinary.star2.m_KruckowTopLambda << TAB << testBinary.star1.m_KruckowMidLambda << TAB << testBinary.star2.m_KruckowMidLambda << TAB << testBinary.star1.m_KruckowBotLambda << TAB << testBinary.star2.m_KruckowBotLambda << TAB << testBinary.star1.m_Metallicity << TAB << testBinary.star2.m_Metallicity << TAB << testBinary.m_massTransferTrackerHistory << std::endl;
        
        int j=0;
        
        do
        {   
            
            
            
            dt = std::min(testBinary.star1.makeTimestep(*options),testBinary.star2.makeTimestep(*options));
            dt=testBinary.chooseTimestep(*options,dt,dtPrev);
            
            testBinary.m_timePrev = testBinary.m_time;
            testBinary.m_time += dt;

            //evolve the system
            testBinary.star1.evolveOneTimestep(dt,false,*options);
            testBinary.star2.evolveOneTimestep(dt,false,*options);
            
            testBinary.evaluateBinary(*options, r, dt, dtPrev);
            
            dtPrev = dt;
            
         
            
            firstTimeFlag = false;
            
            data <<std::setprecision(10) << testBinary.m_randomSeed << TAB << dt << TAB << testBinary.m_time << TAB << testBinary.m_SemiMajorAxisPrime*AUToRsol << TAB << testBinary.m_EccentricityPrime << TAB << testBinary.star1.m_Mass0 << TAB << testBinary.star2.m_Mass0 << TAB << testBinary.star1.m_Mass << TAB << testBinary.star2.m_Mass << TAB << testBinary.star1.m_envMass << TAB << testBinary.star2.m_envMass << TAB << testBinary.star1.m_coreMass << TAB << testBinary.star2.m_coreMass << TAB << testBinary.star1.m_HeCoreMass << TAB << testBinary.star2.m_HeCoreMass << TAB << testBinary.star1.m_COCoreMass << TAB << testBinary.star2.m_COCoreMass << TAB << testBinary.star1.m_Radius << TAB << testBinary.star2.m_Radius << TAB << rocheLobeRadius(testBinary.star1.m_Mass, testBinary.star2.m_Mass) << TAB << rocheLobeRadius(testBinary.star2.m_Mass, testBinary.star1.m_Mass) << TAB << testBinary.rocheLobeTracker1 << TAB << testBinary.rocheLobeTracker2 << TAB << testBinary.star1.m_omega << TAB << testBinary.star2.m_omega << TAB << testBinary.star1.m_omegaBreak << TAB << testBinary.star2.m_omegaBreak << TAB << testBinary.star1.m_stellarType << TAB << testBinary.star2.m_stellarType << TAB << testBinary.star1.m_Age << TAB << testBinary.star2.m_Age << TAB << testBinary.star1.m_Luminosity << TAB << testBinary.star2.m_Luminosity << TAB << testBinary.star1.m_Temperature << TAB << testBinary.star2.m_Temperature << TAB << testBinary.star1.m_angularMomentum << TAB << testBinary.star2.m_angularMomentum << TAB << dynamicalTimescale(testBinary.star1.m_Mass, testBinary.star1.m_Radius) << TAB << dynamicalTimescale(testBinary.star2.m_Mass, testBinary.star2.m_Radius) << TAB << thermalTimescale(testBinary.star1.m_Mass, testBinary.star1.m_envMass, testBinary.star1.m_Radius, testBinary.star1.m_Luminosity, testBinary.star1.m_stellarType) << TAB << thermalTimescale(testBinary.star2.m_Mass, testBinary.star2.m_envMass, testBinary.star2.m_Radius, testBinary.star2.m_Luminosity, testBinary.star2.m_stellarType) << TAB << nuclearTimescale(testBinary.star1.m_Mass, testBinary.star1.m_Luminosity) << TAB << nuclearTimescale(testBinary.star2.m_Mass, testBinary.star2.m_Luminosity) << TAB << determineConvergedMassStepZetaThermal(testBinary.star1, *options) << TAB << determineConvergedMassStepZetaThermal(testBinary.star2, *options) << TAB << determineConvergedTimestepZetaNuclear(testBinary.star1, *options) << TAB << determineConvergedTimestepZetaNuclear(testBinary.star2, *options) << TAB << testBinary.star1.m_massLossDiff << TAB << testBinary.star2.m_massLossDiff << TAB << testBinary.star1.m_MassTransferDiff << TAB << testBinary.star2.m_MassTransferDiff << TAB << testBinary.m_massTransferType << TAB << testBinary.m_TotalAngularMomentumPrime << TAB << testBinary.m_TotalEnergyPrime << TAB << testBinary.star1.m_NanjingLambda << TAB << testBinary.star2.m_NanjingLambda << TAB << testBinary.star1.m_LoveridgeLambda << TAB << testBinary.star2.m_LoveridgeLambda << TAB << testBinary.star1.m_KruckowTopLambda << TAB << testBinary.star2.m_KruckowTopLambda << TAB << testBinary.star1.m_KruckowMidLambda << TAB << testBinary.star2.m_KruckowMidLambda << TAB << testBinary.star1.m_KruckowBotLambda << TAB << testBinary.star2.m_KruckowBotLambda << TAB << testBinary.star1.m_Metallicity << TAB << testBinary.star2.m_Metallicity << TAB << testBinary.m_massTransferTrackerHistory << std::endl;
            
        if (j++ > 50000)
        {
            std::cout<<"maximum number of timesteps"<<std::endl;
            break;
        }
        if (testBinary.star1.m_stellarType == MASSLESS_REMNANT or testBinary.star2.m_stellarType == MASSLESS_REMNANT)
            {
                break;
            }
        
        } while (testBinary.star1.m_stellarType < HELIUM_WHITE_DWARF or testBinary.star2.m_stellarType < HELIUM_WHITE_DWARF);
            
        std::cout<<testBinary.m_SemiMajorAxisInitial<<TAB<<testBinary.star1.m_stellarType<<TAB<<testBinary.star2.m_stellarType<<std::endl;
        
        data.close();
        
    }
     
}
 
BOOST_AUTO_TEST_SUITE_END()
