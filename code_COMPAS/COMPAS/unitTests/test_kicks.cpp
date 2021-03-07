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

 
BOOST_AUTO_TEST_SUITE(Kicks)

BOOST_AUTO_TEST_CASE(TestMaxwellianDistribution)
{
    //JWB - Draws kicks from our Maxwellian distribution and checks the RMS matches the chosen one
    //Draws 10000 samples per distribution for distributions with 10 randomly chosen sigmas
    
    gsl_rng* r = getGslRng();
    
    int nSigmas = 10;
    int nSamplesPerSigma = 10000;
    
    double sigmaMin = 0.;
    double sigmaMax = 1000.;

    double meanToleranceFactor = 3.0;
    double varianceToleranceFactor = 3.0;
    
    double pi = getPi();
    
    for (int i=0; i< nSigmas; i++)
    {
        //generate sigma randomly between [sigmaMin,sigmaMax]
        double sigma = (gsl_rng_uniform(r) * (sigmaMax-sigmaMin)) + sigmaMin;
        
        std::vector<double> drawnSamples;
        drawnSamples.resize(nSamplesPerSigma);
        
        //Draw from the distribution
        for (int j=0; j<nSamplesPerSigma; j++)
        {
            double sample =  kickVelocityDistributionMaxwell(r, sigma);
            drawnSamples[j] = sample;
        }
        
        
        //calculate the sample mean
        double mean = 0.;
        double rms = 0.;
        for (int j=0; j<nSamplesPerSigma; j++)
        {
            mean += drawnSamples[j];
            rms += drawnSamples[j]*drawnSamples[j];
        }
        
        mean /= nSamplesPerSigma;
        rms = std::sqrt(rms/nSamplesPerSigma);
        
        //analytic formula for mean
        double trueMean = 2.*sigma*std::sqrt(2./pi);
        
        std::cout<<"sigma\t"<<sigma<<"sample RMS\t"<<rms<<"\tsample mean\t"<<mean<<"\ttrue mean\t"<<trueMean<<std::endl;
        
        //calculate the variance
        double variance = 0;
        for (int j=0; j<nSamplesPerSigma; j++)
        {
            double deviant = (drawnSamples[j] - mean);
            variance += deviant*deviant;
        }
        variance /= nSamplesPerSigma;
        
        //analytic formula for the variance
        double trueVariance = sigma*sigma*((3. * pi) - 8.)/pi;
        
        std::cout<<"sigma\t"<<sigma<<"\tsample variance\t"<<variance<<"\ttrue variance\t"<<trueVariance<<std::endl;
        
        double meanTolerance = meanToleranceFactor * std::sqrt(trueVariance/nSamplesPerSigma);

        double varianceTolerance = varianceToleranceFactor * std::sqrt(trueVariance);

        BOOST_CHECK(std::abs(variance - trueVariance) < varianceTolerance);

        BOOST_CHECK(std::abs(mean - trueMean) < meanTolerance);
    }
        
}

//BOOST_AUTO_TEST_CASE(KickOpposingOrbitCollapsesBinary)
//{
//    //JWB - Tests that a binary which gets a kick which exactly cancels out the orbital motion of one of the
//    // stars subsequently collapses the binary. Repeats this 10 times
//    
//    programOptions options* = getDefaultCompasOptions();
//    
//    for (int i=0; i<10; ++i)
//    {
//    
//        Binary testBinary = generateRandomBinary();
//        
//        testBinary.star1.evolveOneTimestep(0.,false,*options);
//        testBinary.star2.evolveOneTimestep(0.,false,*options);
//        
//        Star starCopy;
//        
//        int whichStarExploded = 0;
//        
//        //evolve both stars until one of them goes supernova, keep a copy of the previous timestep so that we can correct the kick later
//        do {
//        
//            dt = std::min(testBinary.star1.makeTimestep(*options),testBinary.star1.makeTimestep(*options));
//            
//            starCopy = testBinary.star1;
//            
//            testBinary.star1.evolveOneTimeStep(dt,false,*options);
//            
//            if (testBinary.star1.m_stellarType == NEUTRON_STAR or testBinary.star1.m_stellarType == BLACK_HOLE)
//            {
//                whichStarExploded = 1;
//                break;
//            }
//            
//            starCopy = testBinary.star2;
//            
//            if (testBinary.star2.m_stellarType == NEUTRON_STAR or testBinary.star2.m_stellarType == BLACK_HOLE)
//            {
//                whichStarExploded = 2;
//                break;
//            }
//            
//        while (testBinary.star1.m_stellarType < HELIUM_WHITE_DWARF and testBinary.star2.m_stellarType < HELIUM_WHITE_DWARF);
//        
//        //if we've ended up with a binary white dwarf, skip it
//        if (testBinary.star1.m_stellarType != NEUTRON_STAR and testBinary.star2.m_stellarType != NEUTRON_STAR
//        and testBinary.star1.m_stellarType != BLACK_HOLE and testBinary.star2.m_stellarType != BLACK_HOLE)
//        {
//            --i;
//            continue;
//        }
//        
//        double kickVelocity;
//        double theta, phi;
//        
//        //get the kick from the exploded star
//        if (whichStarExploded == 1) kickVelocity = testBinary.star1.m_drawnKickVelocity;
//        else if (whichStarExploded == 2) kickVelocity = testBinary.star1.m_drawnKickVelocity;
//        else BOOST_FAIL ("made a compact object without a supernova");
//        
//        theta = testBinary.m_Theta;
//        phi = testBinary.m_Phi;
//        
//        
//    }
//}
 

 
BOOST_AUTO_TEST_SUITE_END()
