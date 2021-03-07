//
//  remnantStructure.cpp

#include "remnantStructure.h"

double SchwarzschildRadiusBlackHole(double Mass){
    /*
     Calculate the Schwarzschild radius of a black hole
     
     Given by Equation 94 in Hurley et al 2000
     
     Parameters
     ------------
     Mass : double
     Mass of BH in Msol
     
     Returns
     --------
     R_BH : double
     Radius of black hole in Rsol

     */
    // 2 G M / c^2 -> convert to nice units
    return 4.24E-6 * Mass;
}

double luminosityBlackHole(double Mass){
    /*
     Calculate the luminosity of a black hole - should be negligible - set small to avoid division by zero
     
     Given by Equation 95 in Hurley et al 2000
     
     Parameters
     -----------
     Mass : double
     Mass in Msol
     
     Returns
     --------
     L_BH : double
     Luminoisty in Lsol

     */
    return 1E-10;
}

double radiusNeutronStar(){
    /*
     Assign the radius of a NS as ~10km = 1.4E-5 Rsol
     
     Parameters
     -----------
     
     Returns
     --------
     R_NS : double
     Radius of a NS in Rsol

     */
    return (1.0 / 7.0) * 1E-4; // 10km in Rsol
}

double luminosityNeutronStar(double time, double Mass){
    /*
     Calculate the luminosity of a NS
     
     Given by Equation 93 in Hurley et al 2000
     
     Parameters
     -----------
     time : double
     Time since formation of the object in Myrs
     Mass : double
     Mass of the NS in Msol
     
     Returns
     --------
     L_NS : double
     Luminosity of a NS in Lsol

     */
    double top      = 0.02 * Mass;
    double bottom   = pow(std::max(time, 0.1), 2.0);
    return top / bottom;
}

double whiteDwarfA(int WD_TYPE){
    /*
     Calculate the effective baryon number for WD composition
     
     Given just after Equation 90 in Hurley et al 2000
     
     Parameters
     -----------
     WD_TYPE : int
     Which type (He/CO/ONe) of WD
     
     Returns
     --------
     A : double
     Effective baryon number

     */
    
    if(WD_TYPE == HELIUM_WHITE_DWARF){
        return 4.0;
    }
    else if(WD_TYPE == CARBON_OXYGEN_WHITE_DWARF){
        return 15.0;
    }
    else if(WD_TYPE == OXYGEN_NEON_WHITE_DWARF){
        return 17.0;
    }
    else{
        std::cout << "This isn't a white dwarf!" << std::endl;
        return 4.0;
    }
}

double luminosityWhiteDwarf(double time, double Mass, double Metallicity, int WD_TYPE){
    /*
     Calculate the luminosity of a white dwarf (WD) as it cools as a function of time
     
     Parameters
     -----------
     time : double
     Time since WD formation in MYRs
     Mass : double
     Mass in Msol
     WD_TYPE : int
     Type of WD (He/CO/ONe)
     
     Returns
     --------
     L_WD : double
     Luminosity of WD in Lsol

     */
    
    double top  = 635.0 * Mass * pow(Metallicity, 0.4);
    
    double A    = whiteDwarfA(WD_TYPE);
    
    double bottom = pow((A * (time + 0.1)), (1.4));
    
    return top / bottom;
}

double radiusWhiteDwarf(double Mass){
    /*
     Calculate the radius of a white dwarf
     
     Given by Equation 91 in Hurley et al 2000, taken from Tout et al 1997
     
     Parameters
     -----------
     Mass : double
     Mass in Msol
     
     Returns
     --------
     R_WD : double
     Radius of a WD in Rsol (since WD is ~ Earth sized, expect answer around 0.009)

     */
    
    double RNS = radiusNeutronStar(); // Taken to be 10km
    
    return std::max(RNS, 0.0115 * sqrt(pow((Mch/Mass), (2.0/3.0)) - pow((Mass/Mch), (2.0/3.0))));
}




double neutronStarRemnantMass(double coreMass){
    /*
     Calculate the mass of the neutron star remnant
     
     One of the functions that will be replaced by Fryer et al prescription?
     
     Given by Equation 92 in Hurley et al 2000
     
     Parameters
     -----------
     coreMass : double
     Core mass in Msol
     
     Returns
     --------
     MNS : double
     Neutron star mass in Msol

     */
    return 1.17 + 0.09 * coreMass;
}



double finalBlackHoleMass(double progenitorMass, programOptions options){
    /*
     Function which returns remnant mass based on the chosen prescription
     
     Parameters
     -----------
     options : programOptions
        Class containing all of the user specified options including the remnantMassPrescription
     
     Returns
     --------

     */
    
    if(boost::iequals(options.remnantMassPrescription, "postItNote")){
        std::cout << "Not yet implemented!" << std::endl;
    }
    else{
        std::cout << "Error, remnant mass prescription " << options.remnantMassPrescription << " not supported." << std::endl;
    }
    
    return 0;
}



double calculateTemperature(double luminosity, double radius){
    /*
     Calculate the effective temperature of the star in Tsol
     
     Parameters
     -----------
     [luminosity] =  Lsol
     [radius]     =  Rsol
     
     Returns
     --------
     [T] = Tsol

     */
    return pow(luminosity,(1.0/4.0))  *pow(radius,(-1.0/2.0));
}




double calculateTemperatureKelvin(double luminosity, double radius){
    /*
     Calculate the effective temperature of the star given the 
     Luminosity in Lsol and Radius in Rsol
     
     Parameters
     -----------
     [luminosity] =  Lsol
     [radius]     =  Rsol
     
     Returns
     --------
     [T] = K

     */
    return calculateTemperature(luminosity, radius)  *Tsol;
}



// For comparrison, include our post-it note prescription
double postItNotBlackHoleMass(double initialMass){
    /*
     Our post-it note prescription based on Woosley and Heger 2002
     
     Only applicable to zero metallicity stars with no mass loss.
     
     Parameters
     -----------
     initialMass : double
     Initial mass in Msol
     
     Returns
     --------
     finalMass : double
     Final mass in Msol
     
     */
    
    const double alpha = 3.9;       // Slope of initial/final mass relation
    //const double beta = 0.6;      // Constant fraction
    const double beta = 1.0;        // Constant fraction
    
    double finalMass = 0.0;
    
    if(initialMass <= 25.0){
        finalMass = 0.0;
    }
    else if(initialMass > 25.0 and initialMass <= 50.0){
        finalMass = 30.0 * pow((initialMass/50.0), (alpha));
    }
    else if(initialMass > 50.0){
        finalMass = beta * initialMass;
    }
    
    return finalMass;
    
}
