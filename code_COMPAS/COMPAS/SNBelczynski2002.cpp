//
//  finalBlackHoleMass.cpp


// Include header
#include "SNBelczynski2002.h"


void Belczynski2002Remnant(double *Mass, double *fallbackFraction, double COCoreMass){
    /*
     Formula used by Hurley code not given in Hurley et al 2000 but in Belczynski et al 2002
     
     Parameters
     -----------
     [Mass] :  Msol
     fallbackFraction : Fraction of mass falling back onto compact object [0,1]
     [COCoreMass]: ] Msol
     
     Sets the following pointers
     --------
     fallbackFraction : Fraction of mass falling back onto compact object [0,1]
     [Mass] : Remnant mass in Msol

     */
    
    // Calculate iron core mass and fallback
    double McFeNi = Belczynski2002IronCoreMass(COCoreMass);
    double fb     = Belczynski2002Fallback(*Mass, COCoreMass);
    
    *fallbackFraction = fb;
    *Mass = McFeNi + fb * (*Mass - McFeNi);
    
}







double Belczynski2002IronCoreMass(double COCoreMass){
    /*
     Calculate the iron-nickel core mass as a function of CO core mass
     
     Parameters
     -----------
     COCoreMass : double
     Carbon Oxygen core mass in Msol
     
     Returns
     --------
     McFeNi : double
     Iron-Nickel core mass in Msol

     */
    
    double McFeNi = 0.0;
    
    if(COCoreMass < 2.5){
        McFeNi = (0.161767 * COCoreMass) + 1.067055;
    }
    else{
        McFeNi = (0.314154 * COCoreMass) + 0.686088;
    }
    
    return McFeNi;
}

double Belczynski2002Fallback(double Mass, double COCoreMass){
    /*
     Calculate fallback using linear interpolation
     
     Parameters
     -----------
     Mass : double
     Mass in Msol
     COCoreMass : double
     CO core mass in Msol
     
     Returns
     --------
     fb : double
     Fallback fraction between 0 and 1

     */
    if(COCoreMass <= 5.0){
        return 0.0;
    }
    else if(COCoreMass > 5.0 and COCoreMass < 7.6){
        return (COCoreMass - 5.0)/2.6;
    }
    else{
        return 1.0;
    }
}










///////////////////////////////////////////////////////////////////////////////////////
//                          NOTES
///////////////////////////////////////////////////////////////////////////////////////



//double finalBlackHoleMass(double progenitorMass, programOptions options){
//
//    // This needs to be passed the options so that you can choose between using various prescriptions.
//
//    //if(options.)
//
//    // Should I assume that I have already tested for masses greater than the cutoff mass
//    // progenitorMass should be in solar masses
//    // Note also that progenitorMass is passed by reference
//    // Is using static variables here the correct thing to do? I only want to declare them once. I guess what I should do is have these be options/ in the constants file.
//
//    //static double minimumBlackHoleProgenitorMass    = 25.0;             // Minimum mass for progenitor to form black hole
//    static double initialMassAtSlopeChange          = 50.0;             // Initial mass at point where slope changes
//    static double finalMassAtSlopeChange            = 30.0;             // Final mass at point where slope changes
//    static double alpha                             = 3.9;              // Slope
//    static double constantOffset                    = 0.6;              // Constant offset after slope change
//    double finalMass                                = 0.0;                // Declare variable for final mass
//
//    if(progenitorMass < minimumBlackHoleProgenitorMass){
//        // I don't want it. This is really just like error checking since I should already have gotten rid of it.
//        finalMass = 0.0;
//    }
//    else if(progenitorMass >= minimumBlackHoleProgenitorMass and progenitorMass < initialMassAtSlopeChange){
//        // Power law
//        finalMass = finalMassAtSlopeChange * pow((progenitorMass/initialMassAtSlopeChange), alpha);
//    }
//    else if(progenitorMass >= initialMassAtSlopeChange){
//        // Constant offset
//        finalMass = constantOffset * progenitorMass;
//    }
//    else{
//        // Error checking
//        std::cerr << "ERROR: finalBlackHoleMass" << std::endl;
//    }
//
//    // Change the mass (from when this was a void function)
//    //progenitorMass = finalMass;
//
//    return finalMass;
//
//}
