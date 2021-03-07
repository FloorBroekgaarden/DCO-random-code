//
//  test_supernovae.cpp
//  
//  Created by Jim Barrett on 12/01/2017.
//
//    unit tests for the functionality of supernovae
//

#define BOOST_TEST_DYN_LINK

#include "../testingUtils.h"


#include <iostream>
#include <boost/test/unit_test.hpp>
 
BOOST_AUTO_TEST_SUITE(Supernovae)
 
BOOST_AUTO_TEST_CASE(SupernovaReducesMass)
{
    //JWB - Tests that a star loses mass during supernova. Evolves 10 stars until they go supernova, keeping
    //      a copy 1 timestep behind, then compares the mass just before and after the supernova.
    
    programOptions* options = getDefaultCompasOptions();

    for (int i=0; i<10; ++i)
    {
          
        Star testStar = generateRandomStar();
        
        double prevMass=testStar.m_Mass;

        double dt = 0.;
        //evolve the star to supernova/white dwarf formation, keeping a copy of the mass from one step behind
        do 
        {
            prevMass=testStar.m_Mass;

            dt = testStar.makeTimestep(*options);
            
            testStar.evolveOneTimestep(dt, false, *options);
            
        } while (testStar.m_stellarType < HELIUM_WHITE_DWARF);

        // if no supernova then loop again
        if (testStar.m_stellarType != NEUTRON_STAR and testStar.m_stellarType != BLACK_HOLE)
        {
            i--;
            continue;
        }
        else
        {
            BOOST_CHECK(testStar.m_Mass < prevMass);
        }
    }
    
}
 
BOOST_AUTO_TEST_CASE(multiplication2)
{
    int f = 10;
    int a = 5;
 
    BOOST_CHECK(f*a == 50);
}
 
BOOST_AUTO_TEST_SUITE_END()
