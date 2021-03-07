//
//  main.cpp
//  
//  Created by Jim Barrett on 12/01/2017.
//
//	Main file to ties together unit_tests
//

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN
#define BOOST_TEST_MODULE "COMPAS Tests"
#include <boost/test/included/unit_test.hpp>

#include "./test_supernovae.cpp"
#include "./test_kicks.cpp"
#include "./test_orbits.cpp"
#include "./test_sse.cpp"
#include "./test_commonEnvelopes.cpp"
