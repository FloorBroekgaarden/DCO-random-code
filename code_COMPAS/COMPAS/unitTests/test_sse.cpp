//
//  test_sse.cpp
//  
//  Created by Jim Barrett on 23/01/2017.
//
//    unit tests for the functionality of single star evolution
//

#define BOOST_TEST_DYN_LINK

#include "../testingUtils.h"


#include <iostream>
#include <boost/test/unit_test.hpp>
 
BOOST_AUTO_TEST_SUITE(SSE)
 
BOOST_AUTO_TEST_CASE(NeutronStarsFormed)
{
    //JWB - Tests that stars with progenitor masses between 7 and 20 solar masses produce neutron stars at the end
    //      of their evolution. Evolves 20 stars with masses drawn randomly in this range
    
    double minProgMass = 7.;
    double maxProgMass = 20.;

    programOptions* options = getDefaultCompasOptions();

    gsl_rng* r = getGslRng();

    for (int i=0; i<20; ++i)
    {
        double mass = minProgMass + (maxProgMass-minProgMass)*gsl_rng_uniform(r); 
        Star testStar = generateStarWithSpecificMass(mass);

        double dt = 0.;
	testStar.evolveOneTimestep(dt,false,*options);

        //evolve the star to supernova/white dwarf formation
        do 
        {

            dt = testStar.makeTimestep(*options);
            
            testStar.evolveOneTimestep(dt, false, *options);
            
        } while (testStar.m_stellarType < HELIUM_WHITE_DWARF);
        
        
        BOOST_CHECK(testStar.m_stellarType == NEUTRON_STAR);
        
    }
    
}
 
BOOST_AUTO_TEST_SUITE_END()
