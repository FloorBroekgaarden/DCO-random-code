//
//  magneticBraking.cpp


#include "magneticBraking.h"

// Derivation still in need of further testing

///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
const double gamma_MB = 5*pow(10.0,-29.0)*pow(100*Rsol,2.0)/(year*MyearToyear);  // [gamma_MB] = yr Rsol-2

// Functions used for derivation as in Repetto et. al 2014
double magneticBrakingSpinDownRepetto(double Radius, double Omega, double dt){
// Spin-down of the star due to MB, according to Repetto et al. 2013, eq. 1
    return -gamma_MB*pow(Radius,2.0)*pow(Omega,3.0)*dt;
}

double magneticBrakingAngularMomentumLossRepetto(double k, double Mass, double Radius, double Omega, double dt){
// Non-coupled method, according to Repetto et al. 2013, eq. 12
// Rate of angular momentum loss of the star in the stellar wind due to MB
    return k*Mass*pow(Radius,2.0)*magneticBrakingSpinDownRepetto(Radius, Omega, dt);
}

double magneticBrakingSemiMajorAxisRepetto(double k, double Mass, double Radius, double MassComp, double a, double dt){
// Non-coupled method, according to Repetto et al. 2013, eq. 13
// Semi-major axis decay rate of the binary due to MB of 1 star
    const double GMyear = G1*pow(10.0,12.0);
    double w = sqrt(GMyear*(Mass+MassComp)*pow(a,-3.0));
    std::cout <<  "KRep " << gamma_MB << std::endl;
    return -2*gamma_MB*k*pow(Radius,4.0)*(Mass+MassComp)*pow(MassComp,-1.0)*pow(w,2.0)*pow(a,-1.0)*dt;
    //    return  -2.0*gamma_MB*k*pow(Radius,4.0)*GMyear*pow(Mass+MassComp,2.0)*pow(MassComp,-1.0)*pow(a,-4.0)*dt;

}

///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////


// Function used for derivation as in Belczynski et ak. 2008
double magneticBrakingRappaport(double Mass, double Radius, double MassComp, double a, double dt){
    // Semi-major axis decay rate of the Binary as suggested by Rappaport et al. 1983
    const double KRappaport = 5.8*pow(10.0,-22.0);
    const double GMyear = G1*pow(10.0,12.0);
    double w = sqrt(GMyear*(Mass+MassComp)*pow(a,-3.0));
    //return  2.0*KRappaport*pow(Radius,2.0)*pow(Mass+MassComp,2.0)*GMyear*pow(MassComp,-1.0)*pow(a,-4.0)*dt;
//    return  -2.0*KRappaport*pow(Radius,2.0)*pow(Mass+MassComp,2.0)*pow(MassComp,-1.0)*pow(a,-1.0)*pow(w,2.0)*dt;
    std::cout << "KRap " << KRappaport << std::endl;
        return  -2.0*KRappaport*pow(Radius,4.0)*(Mass+MassComp)*pow(MassComp,-1.0)*pow(a,-1.0)*pow(w,2.0)*dt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////

double magneticBrakingIvanova(double Omega, double Mass, double Radius, double MassComp, double a, double dt){
    // Semi-major axis decay rate of the Binary as suggested by Ivanova et al. 1983
    const double GMyear = G1*pow(MyearToyear,2.0);
    const double OmegaX = 9.45*pow(10.0,8.0);   // [OmegaX] = Myr-1
    double KIvanova = 0.0;
    double w = sqrt(GMyear*(Mass+MassComp)*pow(a,-3.0));
    
    if(Omega <= OmegaX){
        KIvanova = 619.2*(pow(9.45*pow(10.0,7.0),-3.0));
//        return  -2.0*KIvanova*GMyear*pow(Radius,4.0)*pow(Mass+MassComp,2.0)*pow(Mass*MassComp,-1.0)*pow(a,-4.0)*dt;
        std::cout <<  "KIv " << KIvanova<< std::endl;
        return  -2.0*KIvanova*GMyear*pow(Radius,4.0)*(Mass+MassComp)*pow(Mass*MassComp,-1.0)*pow(a,-2.0)*pow(w,2.0)*dt;

    }
    else{
        KIvanova = 619.2*pow(10.0,1.7)*(pow(9.45*pow(10.0,7.0),-1.3));
        std::cout <<  "Kr " << KIvanova<< std::endl;
        return  -2.0*KIvanova*pow(Radius,4.0)*(Mass+MassComp)*pow(Mass*MassComp,-1.0)*pow(w,0.3)*dt;
    }
}