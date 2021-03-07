#include "AISclass.h"

// Default constructor
AISvariables::AISvariables(){

    DrawingFromAISdistributions = false;  // if true we sample from AIS distributions  
    RandomGaussianDraw = 99999;
    fractionSampled = 0;
    fexplAIS = 1;
    CounterDCOsAIS = 0; // CounterrDCOsAIS counts the number of DCOs of interest ("hits") in the exploratory phase AIS

}

