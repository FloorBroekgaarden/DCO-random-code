//
//  test_supernovae.cpp
//  
//  Created by Jim Barrett on 12/01/2017.
//
//    unit tests for the functionality of orbits and orbital calculations
//

#define BOOST_TEST_DYN_LINK

#include "../testingUtils.h"

#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(Orbits)
 
BOOST_AUTO_TEST_CASE(ZeroOrbitalVelocityLeadsToCollapse)
{
    //JWB - Tests that killing the orbital velocity of a binary collapses it on the next timestep
    //      evolves a Binary for a random number of timesteps before reducing its orbital velocity to zero
    //      repeats this for 20 binaries
    
    programOptions* options = getDefaultCompasOptions();

}

BOOST_AUTO_TEST_CASE(RemoveHalfMassBecomesUnbound)
{
    //JWB - Tests that removing half the mass of the system unbinds it. Takes 10 random systems, evolves them
    //      for a random number of timesteps and then drops half of the mass of the system (at a random ratio
    //      between the two stars)
    
//    programOptions* options = getDefaultCompasOptions();
//    
//    gsl_rng* r = getGslRng();
//    
//    int minimumTimesteps = 100;
//    int maximumTimesteps = 1000;
//    
//    //flag to skip the binary. Currently used if we form a massless remnant
//    bool skipBinary = false;
//    
//    for (int i=0; i<10; ++i)
//    {
//    
//        BinaryStar testBinary = generateRandomBinary();
//        
//        testBinary.star1.evolveOneTimestep(0.,false,*options);
//        testBinary.star2.evolveOneTimestep(0.,false,*options);
//        
//        double dt;
//        double dtPrev = 0.;
//        
//        int nTimesteps = (gsl_rng_uniform(r) * (maximumTimesteps - minimumTimesteps)) + minimumTimesteps;
//        
//        for (int j=0; j<nTimesteps; j++)
//        {
//            //calculate the relevant timescales
//            dt = std::min(testBinary.star1.makeTimestep(*options),testBinary.star2.makeTimestep(*options));

//            if (j!=0)
//            {
//                dt=testBinary.chooseTimestep(*options,dt,dtPrev);
//            }

//            //evolve the system
//            testBinary.star1.evolveOneTimestep(dt,false,*options);
//            testBinary.star2.evolveOneTimestep(dt,false,*options);
//            
//            testBinary.evaluateBinary(*options, r, dt, dtPrev);
//            
//            dtPrev = dt;
//            
//            //it doesn't make sense to unbound a massless remnant
//            if (testBinary.star1.m_stellarType == MASSLESS_REMNANT or testBinary.star2.m_stellarType == MASSLESS_REMNANT)
//            {
//                skipBinary = true;
//                break;
//            }
//            
//        }
//        
//        //if the binary is no use to us, then skip it
//        if (skipBinary)
//        {
//            skipBinary = false;
//            --i;
//            continue;
//        }

//        //drop half the mass of the system from the heavier star
//        double halfTotalMass = 0.5*(testBinary.star1.m_Mass + testBinary.star2.m_Mass);
//        
////        std::cout<<"first\t"<<testBinary.star1.m_Mass<<"\t"<<testBinary.star2.m_Mass<<"\t"<<halfTotalMass<<std::endl;
//        
//        if (testBinary.star1.m_Mass > testBinary.star2.m_Mass)
//        {
//            testBinary.star1.m_Mass -= halfTotalMass;
//        }
//        else
//        {
//            testBinary.star2.m_Mass -= halfTotalMass;
//        }
//        
//        testBinary.star1.evolveOneTimestep(0.,false,*options);
//        testBinary.star2.evolveOneTimestep(0.,false,*options);
//        
////        std::cout<<"second\t"<<testBinary.star1.m_Mass<<"\t"<<testBinary.star2.m_Mass<<"\t"<<halfTotalMass<<std::endl;
//        
//        //evolve it all one more timestep
//        dt = std::min(testBinary.star1.makeTimestep(*options),testBinary.star2.makeTimestep(*options));
//        dt = testBinary.chooseTimestep(*options,dt,dtPrev);
//        
////        std::cout<<"third\t"<<testBinary.star1.m_Mass<<"\t"<<testBinary.star2.m_Mass<<"\t"<<halfTotalMass<<std::endl;
////        std::cout<<"third (prevs)\t"<<testBinary.star1.m_MassPrev<<"\t"<<testBinary.star2.m_MassPrev<<"\t"<<halfTotalMass<<std::endl;
//        
//        testBinary.evaluateBinary(*options, r, dt, dtPrev);
//        
////        std::cout<<"fourth\t"<<testBinary.star1.m_Mass<<"\t"<<testBinary.star2.m_Mass<<"\t"<<halfTotalMass<<std::endl;
////        std::cout<<"fourth (prevs)\t"<<testBinary.star1.m_MassPrev<<"\t"<<testBinary.star2.m_MassPrev<<"\t"<<halfTotalMass<<std::endl;
//        
//        //check if it's bound
//        bool unboundSystem = -(testBinary.m_totalOrbitalEnergyPrime/testBinary.m_totalOrbitalEnergyPrev) > 0.0;
//        
//        double mu = (testBinary.star1.m_Mass * testBinary.star2.m_Mass)/(testBinary.star1.m_Mass + testBinary.star2.m_Mass);
//        double M = testBinary.star1.m_Mass + testBinary.star2.m_Mass;
//        
//        double E = testBinary.orbitalEnergy(mu,M,testBinary.m_SemiMajorAxis)
//        
//        BOOST_CHECK(E > 0.);
//              
//    }

}

BOOST_AUTO_TEST_SUITE_END()
