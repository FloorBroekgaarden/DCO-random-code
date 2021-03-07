//
//  fryerFallback.cpp

#include "SNFryer2012.h"

// Implementation
// Put everything together in the function you actually call!
void FryerRemnant(double *Mass, double COCoreMass, double *fallbackFraction, int engine, unsigned long randomSeed, const programOptions &options){
    /*
     Function to put together the various functions below based on users choice of SN engine.

     Parameters
     -----------
     POINTERS:
     Mass : Mass star in [Msol]  
     fallbackFraction : Fraction of mass [0,1] falling back onto proto compact object 

     VALUES
     COCoreMass : COCore mass at time of the SN in [Msol]
     engine : int Which supernova engine to use

     */
    bool   debugging = false;
    double Mproto    = 0.0;
    double fb        = 0.0;
    double Mfb       = 0.0;
    double Mrem_bary = 0.0;
    double Mrem_grav = 0.0;
    if(debugging){std::cout<<"Seed---------------"<<randomSeed<<"\n"
                             "In FryerRemnant() function"<<std::endl;}  

    if(engine == SN_DELAYED){
        if(debugging){std::cout<<"Engine = DELAYED"<<std::endl;}
        // Calculate Mproto
        Mproto = FryerCOCoreToProtoMassDelayed(COCoreMass);
        // Calculate fb, Mfb
        fb = FryerFallbackDelayed(*Mass, Mproto, COCoreMass);
        Mfb = FryerMassFallbackGeneral(*Mass, Mproto, fb);
        // Calculate Mremnant, baryonic, Mremnant, gravitational
        Mrem_bary = FryerBaryonicRemnantMassGeneral(Mproto, Mfb);
        Mrem_grav = FryerGravitationalRemnantMassGeneral(Mrem_bary, options);
    }


    else if (engine == SN_RAPID){
         if(debugging){std::cout<<"Engine = RAPID"<<std::endl;}
        // Calculate Mproto
        Mproto = FryerCOCoreToProtoMassRapid(COCoreMass);
        // Calculate fb, Mfb
        fb = FryerFallbackRapid(*Mass, Mproto, COCoreMass);
        Mfb = FryerMassFallbackGeneral(*Mass, Mproto, fb);
        // Calculate Mremnant, baryonic
        Mrem_bary = FryerBaryonicRemnantMassGeneral(Mproto, Mfb);
        Mrem_grav = FryerGravitationalRemnantMassGeneral(Mrem_bary, options);
    }


    else{
        std::cout << "Unrecognised Fryer supernova engine: choose RAPID or DELAYED (default?)" << std::endl;
    }


    if(debugging){
            std::cout << "MpreSN = " << *Mass << std::endl;
            std::cout << "COCore pre SN "<< COCoreMass <<std::endl;
            std::cout << "Mproto = " << Mproto << std::endl;
            std::cout << "Mfb = " << Mfb << std::endl;
            std::cout << "Mbary = " << Mrem_bary << std::endl;
            std::cout << "Mgrav = " << Mrem_grav << std::endl;
            std::cout << "exit - FryerRemnant()\n -----" << std::endl;
    }


    //Set Pointers
    *fallbackFraction = fb;
    *Mass  =  Mrem_grav;
}




double FryerSolveQuadratic(double a, double b, double c){
    /*
     Solve quadratic ax^2 + bx + c
     return either 0, 1, or 2 roots depending on discriminant

     */
    double discriminant = b*b - 4.0*a*c;
    double x = 0.0;
    
    if(discriminant > 0.0){
        //std::cout << "There are 2 real roots" << std::endl;
        double xplus  = (-b + sqrt(discriminant))/(2.0*a);
        double xminus = (-b - sqrt(discriminant))/(2.0*a);
        x = std::max(xplus, xminus);
    }
    else if(discriminant == 0){
        //std::cout << "There is 1 repeated root" << std::endl;
        x = -b/(2.0*a);
    }
    else{
        std::cerr << "Error. No real roots" << std::endl;
    }
    
    return x;
}

double FryerBaryonicRemnantMassGeneral(double Mproto, double Mfb){
    /*
     Equation 12 in Fryer et al 2012
     
     Parameters
     -----------
     Mproto : double
        Mass of proto compact object in Msol
     Mfb : double
        Mass falling back onto proto compact object Mfb = fb*(MpreSN - Mproto)
     
     Returns
     --------
     Mremnant, baryonic : double
        Baryonic mass of remnant

     */
    return Mproto + Mfb; //Mfb = fb*(MpreSN - Mproto);
}

double blackHoleFormationNeutrinoMassLoss(double Mbary, const programOptions &options){
    /*
     Parameters
     -----------
     Mbary : double
     Baryonic remnant mass in Msol
     
     Returns
     --------
     Mremnant, gravitational : double
     Gravitational mass of remnant in Msol
    */
    double Mgrav = 0.0;
    
    if(options.neutrinoMassLossAssumptionBH == NEUTRINO_MASS_LOSS_BH_FIXED_FRACTION){
        Mgrav = Mbary * (1.0 - options.neutrinoMassLossValueBH);
    }
    else if(options.neutrinoMassLossAssumptionBH == NEUTRINO_MASS_LOSS_BH_FIXED_MASS){
        Mgrav = Mbary - options.neutrinoMassLossValueBH;
    }
    else{
        std::cerr << "Error in blackHoleFormationNeutrinoMassLoss setting BH mass" << std::endl;
    }
    return Mgrav;
}

double FryerGravitationalRemnantMassGeneral(double Mbary, const programOptions &options){
    /*
     Equations 13 & 14 in Fryer et al 2012
     
     Parameters
     -----------
     Mbary : double
     Baryonic remnant mass in Msol
     
     Returns
     --------
     Mremnant, gravitational : double
     Gravitational mass of remnant in Msol

     */
    double Mbary_maximumNeutronStarMass = (0.075 * pow(options.maximumNeutronStarMass, 2)) + options.maximumNeutronStarMass;
    // this is the mass Mbary for which the gravitational (final) mass will be below the options.maximumNeutronStarMass limit given. 
    // see Eq. 13 in Fryer+2012
    if(Mbary < Mbary_maximumNeutronStarMass){
        // NS
        // Solve quadratic equation
        double tmp = FryerSolveQuadratic(0.075, 1.0, -Mbary);
        return tmp;
    }
    else if(Mbary >= options.maximumNeutronStarMass){
        // BH
        return blackHoleFormationNeutrinoMassLoss(Mbary, options);
    }
    else{
        std::cerr << "Error in FryerGravitationalRemnantMassGeneral" << std::endl;
        return 0.0;
    }
}

double FryerMassFallbackGeneral(double MpreSN, double Mproto, double fb){
    /*
     Equation 11 in Fryer et al 2012
     
     Parameters
     -----------
     MpreSN : double
     Pre supernova stellar mass in Msol
     COCoreMass : double
     Pre supernova Carbon Oxygen (CO) core mass in Msol
     fb : double
     Fallback [0,1]
     
     Returns
     --------
     Mfb : double
     Mass falling back onto proto object

     */
    return fb * (MpreSN - Mproto);
}


double FryerCOCoreToProtoMassRapid(double COCoreMass){
    /*
     Equation 15 in Fryer et al 2012, based on Woosley et al 2002
     
     Sets Mproto = 1.0 Msol regardless of progenitor

     Parameters
     ------------
     COCoreMass : double
     Carbon Oxygen (CO) core mass in Msol
     
     Returns
     --------
     Mproto : double
     Mass of Fe/Ni proto core in Msol

     */
    return 1.0;
}

double FryerFallbackRapid(double MpreSN, double Mproto, double COCoreMass){
    /*
     Calculate fallback using the rapid prescription
     
     Equations 15 & 16 from Fryer et al 2012
     
     Parameters
     -----------
     MpreSN : double
     Pre supernova stellar mass in Msol
     Mproto : double
     Fe/Ni proto compact object mass in Msol
     COCoreMass : double
     Pre supernova CO core mass in Msol
     
     Returns
     ---------
     fb : double
     Fraction of mass falling back onto proto object

     */
    double fb = 0.0;

    bool debugging = false;

    if(COCoreMass < 2.5){
        fb = 0.2/(MpreSN - Mproto);
    }
    else if(COCoreMass >=2.25 and COCoreMass < 6.0){
        fb = (0.286 * COCoreMass - 0.514)/(MpreSN - Mproto);
    }
    else if(COCoreMass >= 6.0 and COCoreMass < 7.0){
        fb = 1.0;
    }
    else if(COCoreMass >= 7.0 and COCoreMass < 11.0){
        double a1 = 0.25 - (1.275/(MpreSN - Mproto));
        double b1 = -11.0 * a1 + 1.0;
        fb = a1 * COCoreMass + b1;
    }
    else if(COCoreMass >= 11.0){
        fb = 1.0;
    }
    else{
        std::cerr << "Error in fallbackRapid" << std::endl;
        fb = 0.0;
    }

    if(fb > 1)
    {
        if(debugging){
            std::cout << "Fallback delayed gives fb > 1.0 so setting fb = 1" << std::endl;
        }
        fb = 1.0;
    }
    else if (fb < 0)
    {
        if(debugging){
            std::cout << "Fallback delayed gives fb < 0 so setting fb = 0" << std::endl;
        }
        fb = 0;
    }
  
  return fb;
}

double FryerCOCoreToProtoMassDelayed(double COCoreMass){
    /*
     Equation 18 in Fryer et al 2012
     
     Parameters
     -----------
     COCoreMass : double
     Mass of the CO core in Msol
     
     Returns
     --------
     Mproto : double
     Mass of the Fe/Ni proto core in Msol

     */
    
    if(COCoreMass < 3.5){
        return 1.2;
    }
    else if(COCoreMass >= 3.5 and COCoreMass < 6.0){
        return 1.3;
    }
    else if(COCoreMass >= 6.0 and COCoreMass < 11.0){
        return 1.4;
    }
    else if(COCoreMass >= 11.0){
        return 1.6;
    }
    else{
        std::cerr << "Error in COCoreToProtoMassDelayed" << std::endl;
        return 0;
    }
    
}

double FryerFallbackDelayed(double MpreSN, double Mproto, double COCoreMass){
    /*
     Calculate fallback using the delayed prescription
     
     Equation 19 of Fryer et al 2012
     
     Parameters
     -----------
     MpreSN : double
     Pre supernova stellar mass in Msol
     Mproto : double
     Fe/Ni proto compact object mass in Msol
     COCoreMass : double
     Pre supernova CO core mass in Msol
     
     Returns
     --------
     fb : double
     Fallback

     */
    double fb = 0;

    bool debugging = false;

    if(COCoreMass < 2.5){
        fb = 0.2/(MpreSN - Mproto);
    }
    else if(COCoreMass >= 2.5 and COCoreMass < 3.5){
        fb = (0.5 * COCoreMass - 1.05)/(MpreSN - Mproto);
    }
    else if(COCoreMass >= 3.5 and COCoreMass < 11.0){
        double a2 = 0.133 - (0.093)/(MpreSN - Mproto);
        double b2 = -11.0 * a2 + 1.0;
        fb = a2 * COCoreMass + b2;
    }
    else if(COCoreMass >= 11.0){
        fb = 1.0;
    }
    else{
        std::cerr << "Error in fallbackDelayed" << std::endl;
        fb = 0.0;
    }
  
    if(fb > 1)
    {
        if(debugging){
            std::cout << "Fallback delayed gives fb > 1.0 so setting fb = 1" << std::endl;
        }
        fb = 1.0;
    }
    else if (fb < 0)
    {
        if(debugging){
            std::cout << "Fallback delayed gives fb < 0 so setting fb = 0" << std::endl;
        }
        fb = 0;
    }
  return fb;

}




// old functions


// double fitToFigure2(double initialMass){
//     /*
//      Fit to figure 2 from Fryer and Kalogera 2001 for comparisson
     
//      Parameters
//      ------------
//      initialMass : double
//      Initial Mass in Msol
     
//      Returns
//      ---------
//      finalMass : double
//      Final Mass in Msol

//      */
//     const double y2 = 15.0;
//     const double y1 = 1.4;
//     const double x2 = 40.0;
//     const double x1 = 20.0;
//     const double gradient = (y2 - y1)/(x2 - x1);
//     const double intercept = y1 - (gradient * x1);
    
//     double finalMass = 0.0;
    
//     if(initialMass < 20.0){
//         finalMass = 1.4;
//     }
//     else if(initialMass >= 20.0 and initialMass <= 40.0){
//         finalMass = (gradient * initialMass) + intercept;
//     }
//     else if(initialMass > 40.0){
//         finalMass = initialMass;
//     }
//     else{
//         std::cerr << "Error in fitToFigure2" << std::endl;
//         finalMass = initialMass;
//     }
    
//     return finalMass;
    
// }

// double COCoreToProtoMassRapid(double COCoreMass){
//     /*
//      Equation 10 in Fryer et al 2012
     
//      Parameters
//      ------------
//      COCoreMass : double
//      Carbon Oxygen (CO) core mass in Msol
     
//      Returns
//      --------
//      Mproto : double
//      Mass of Fe/Ni proto core in Msol

//      */
//     if(COCoreMass < 4.82){
//         return 1.50;
//     }
//     else if(COCoreMass >= 4.82 and COCoreMass < 6.31){
//         return 2.11;
//     }
//     else if(COCoreMass >= 6.31 and COCoreMass < 6.75){
//         return 0.69 * COCoreMass - 2.26;
//     }
//     else if(COCoreMass >= 6.75){
//         return 0.37 * COCoreMass - 0.07;
//     }
//     else{
//         std::cerr << "Error in COCoreToProtoMass" << std::endl;
//         return 0;
//     }
    
// }

// double gravitationalRemnantMass(double MpreSN, double COCoreMass){
//     /*
//      Equation 13 in Fryer et al 2012
     
//      Parameters
//      ------------
//      MpreSN : double
//      Pre supernova stellar mass in Msol
//      COCoreMass : double
//      Pre supernova Carbon Oxygen (CO) core mass in Msol
     
//      Returns
//      --------
//      Mremnant, gravitational : double
//      Gravitational mass of remnant in Msol

//      */
    
//     double Mbary = baryonicRemnantMass(MpreSN, COCoreMass);
    
//     //std::cout << "MpreSN: " << MpreSN << std::endl;
//     //std::cout << "COCoreMass: " << COCoreMass << std::endl;
//     //std::cout << "Mbary: " << Mbary << std::endl;
    
//     //options.maximumNeutronStarMass

//     if(Mbary < MAXIMUM_NS_MASS){
//         // NS
//         // Solve quadratic equation
//         double tmp = solveQuadratic(0.075, 1.0, -Mbary);
//         return tmp;
        
//     }
//     else if(Mbary >= MAXIMUM_NS_MASS){
//         return 0.9 * Mbary;
//     }
//     else{
//         std::cerr << "Error in gravitationalRemnantMass" << std::endl;
//         return 0;
//     }
// }


// double lowerNSProgenitorMass(double Metallicity){
    
//      Equation 4 from Fryer et al 2012
     
//      Parameters
//      ------------
//      Metallicity : double
//      Metallicity (Z = 0.02 = Zsol)
     
//      Returns
//      --------
//      Lower limit on NS progenitor mass in Msol

     
//     if(Metallicity > 1E-3 * Zsol){
//         return 9.0 + 0.9 * log10(Metallicity / Zsol);
//     }
//     else{
//         return 6.3;
//     }
    
// }

// double intermediateStarRemnantMassesDelay(double MZAMS, double Metallicity){
    
//      Equation 5 from Fryer et al 2012
     
//      Parameters
//      ------------
//      MZAMS : double
//      Zero Age Main Sequence (ZAMS) mass of star (11.0 < Mstar < 30.0)
//      Metallicity : double
//      Metallicity Z (Z = 0.02 = Zsol)
     
//      Returns
//      --------
//      Mremnant, delay : double
//      Remnant mass in Msol

     
    
//     double Zmetal = Metallicity / Zsol;
    
//     return 1.1 + 0.2*exp((MZAMS-11.0)/4.0) - (2.0 + Zmetal) * exp(0.4*(MZAMS - 26.0));
    
// }

// double intermediateStarRemnantMassesRapid(double MZAMS, double Metallicity){
    
//      Equation 6 from Fryer et al 2012
     
//      Parameters
//      ------------
//      MZAMS : double
//      Zero Age Main Sequence (ZAMS) mass of star (11.0 < Mstar < 30.0)
//      Metallicity : double
//      Metallicity Z (Z = 0.02 = Zsol)
     
//      Returns
//      --------
//      Mremnant, rapid : double
//      Remnant mass in Msol

     
    
//     double Zmetal = Metallicity / Zsol;
    
//     if(MZAMS < 22.0){
//         return 1.1 + 0.2 * exp((MZAMS - 11.0)/7.5) + 10.0 * (1.0 + Zmetal) * exp(-(pow((MZAMS - 23.5),(2.0))/pow((1.0 + Zmetal), (2.0))));
//     }
//     else{
//         return intermediateStarRemnantMassesDelay(MZAMS, Metallicity) - 1.85 + 0.25 * Zmetal + 10.0 * (1.0 + Zmetal) * exp(-(pow((MZAMS - 23.5),(2.0))/pow((1.0 + Zmetal), (2.0))));
//     }
    
// }
// should put all together so that mass determines which procedure to use


// double massiveStarRemnantMassDelay(double MZAMS, double Metallicity){
    
//      Equation 7 from Fryer et al 2012
     
//      Parameters
//      -----------
//      MZAMS : double
//      Zero Age Main Sequence (ZAMS) mass of star in Msol (MZAMS > 30.0 Msol)
//      Metallicity : double
//      Metallicity Z (Z = 0.02 = Zsol)
     
//      Returns
//      --------
//      Mremnant, delay : double
//      Remnant mass in Msol

     
    
//     double Zmetal = Metallicity / Zsol;
    
//     double a = 33.5 + (4.75 + 1.25 * Zmetal)*(MZAMS - 34.0);
//     double b = MZAMS - sqrt(Zmetal)*(1.3 * MZAMS - 18.35);
    
//     return std::min(a, b);
// }

// double massiveStarRemnantMassRapid(double MZAMS, double Metallicity){
    
//      Equation 8 from Fryer et al 2012
     
//      Parameters
//      ------------
//      MZAMS : double
//      Zero Age Main Sequence (ZAMS) mass of star in Msol (MZAMS > 30.0 Msol)
//      Metallicity : double
//      Metallicity Z (Z = 0.02 = Zsol)
     
//      Returns
//      --------
//      Mremnant, rapid : double
//      Remnant mass in Msol

     
    
//     double Zmetal = Metallicity / Zsol;
    
//     return massiveStarRemnantMassDelay(MZAMS, Metallicity) - 1.85 + Zmetal*(75.0 - MZAMS)/(20.0);
// }

// double massiveStarRemnantMassWoosley(double MZAMS, double Metallicity){
    
//      Equation 9 from Fryer et al 2012
     
//      Is itself a fit to Woosley and Heger 2002, Heger et al 2003 for solar metallicities
     
//      For lower metallicities, use maximum of Equations 7 and 9
     
//      Parameters
//      -----------
//      MZAMS : double
//      Zero Age Main Sequence (ZAMS) mass of star in Msol (MZAMS > 30.0 Msol)
//      Metallicity : double
//      Metallicity Z (Z = 0.02 = Zsol)
     
//      Returns
//      --------
//      Mremnant : double
//      Remnant mass in Msol
     
     
    
//     if(MZAMS < 90.0){
//         return 1.8 + 0.04*(90.0 - MZAMS);
//     }
//     else{
//         return 1.8 + log10(MZAMS - 89.0);
//     }
    
// }

// double fryerRemnantMass2012(double MZAMS, double Metallicity, int engine){
//     /*
//      Return remnant mass for chosen supernova engine based on initial mass for a single star
     
//      Parameters
//      -----------
//      MZAMS : double
//      Zero Age Main Sequence (ZAMS) mass of star in Msol (MZAMS > 30.0 Msol)
//      Metallicity : double
//      Metallicity Z (Z = 0.02 = Zsol)
//      engine : int
//      Either 'delayed' or 'rapid' strings
     
//      Returns
//      ---------
//      Mremnant : double
//      Remnant mass in Msol

//      */
    
//     // TODO: Replace with string comparrison with boost iequals
    
//     //std::string RAPID = 'rapid';
//     //std::string DELAYED = 'delayed';
//     int DELAYED = 0;
//     int RAPID = 1;
    
//     if(engine == DELAYED){
//         if(MZAMS < 30.0){
//             return intermediateStarRemnantMassesDelay(MZAMS, Metallicity);
//         }
//         else if(MZAMS >= 30.0){
//             return std::max(massiveStarRemnantMassDelay(MZAMS, Metallicity), massiveStarRemnantMassWoosley(MZAMS, Metallicity));
//         }
//         else{
//             std::cerr << "Error in fryerRemnantMass2012" << std::endl;
//             return 0;
//         }
//     }
//     else if(engine == RAPID){
//         if(MZAMS < 30.0){
//             return intermediateStarRemnantMassesRapid(MZAMS, Metallicity);
//         }
//         else if(MZAMS >= 30.0){
//             return std::max(massiveStarRemnantMassRapid(MZAMS, Metallicity), massiveStarRemnantMassWoosley(MZAMS, Metallicity));
//         }
//         else{
//             std::cerr << "Error in fryerRemnantMass2012" << std::endl;
//             return 0;
//         }
//     }
//     else{
//         std::cerr << "Error in fryerRemnantMass2012. Undefined engine." << std::endl;
//         return 0;
//     }
// }

// // IMF shouldn't be necessary, but useful for testing
// double IMFSamples(double u, double mmin, double mmax){
    
//      Draw samples from an IMF between mmin, mmax
     
//      Parameters
//      -----------
//      u : double
//      Random number [0,1]
//      mmin : double
//      Minimum mass
//      mmax : double
//      Maximum mass
     
//      Returns
//      ---------
//      mass : double
//      Mass drawn from IMF

     
    
//     double n = -2.35;
    
//     return pow((u*(pow(mmax, n+1.0) - pow(mmin, n+1.0)) + pow(mmin, n+1.0)), (1.0/(n+1.0)));
// }

// double COCoreToProtoMass(double COCoreMass){
//     /*
//      Equation 10 in Fryer et al 2012
     
//      Parameters
//      ------------
//      COCoreMass : double
//      Carbon Oxygen (CO) core mass in Msol
     
//      Returns
//      --------
//      Mproto : double
//      Mass of Fe/Ni proto core in Msol

//      */
//     if(COCoreMass < 4.82){
//         return 1.50;
//     }
//     else if(COCoreMass >= 4.82 and COCoreMass < 6.31){
//         return 2.11;
//     }
//     else if(COCoreMass >= 6.31 and COCoreMass < 6.75){
//         return 0.69 * COCoreMass - 2.26;
//     }
//     else if(COCoreMass >= 6.75){
//         return 0.37 * COCoreMass - 0.07;
//     }
//     else{
//         std::cerr << "Error in COCoreToProtoMass" << std::endl;
//         return 0;
//     }
    
// }

// double fallbackFraction(double COCoreMass){
    
//      Equation 11 in Fryer et al 2012
     
//      Parameters
//      -----------
//      COCoreMass : double
//      Carbon Oxygen (CO) core mass in Msol
     
//      Returns
//      -------
//      fb : double
//      Fraction of mass falling back onto proto object

     
//     if(COCoreMass < 5.0){
//         return 0;
//     }
//     else if(COCoreMass >= 5.0 and COCoreMass < 7.6){
//         return (0.378 * COCoreMass) - 1.889;
//     }
//     else if(COCoreMass >= 7.6){
//         return 1.0;
//     }
//     else{
//         std::cerr << "Error in fallbackFraction" << std::endl;
//         return 0;
//     }
// }

// double massFallback(double MpreSN, double COCoreMass){
    
//      Equation 11 in Fryer et al 2012
     
//      Parameters
//      -----------
//      MpreSN : double
//      Pre supernova stellar mass in Msol
//      COCoreMass : double
//      Pre supernova Carbon Oxygen (CO) core mass in Msol
     
//      Returns
//      --------
//      Mfb : double
//      Mass falling back onto proto object

     
//     double fb = fallbackFraction(COCoreMass);
//     double Mproto = COCoreToProtoMass(COCoreMass);
//     return fb * (MpreSN - Mproto);
// }

// double baryonicRemnantMass(double MpreSN, double COCoreMass){
    
//      Equation 12 in Fryer et al 2012
     
//      Parameters
//      ------------
//      MpreSN : double
//      Pre supernova stellar mass in Msol
//      COCoreMass : double
//      Pre supernova Carbon Oxygen (CO) core mass in Msol
     
//      Returns
//      --------
//      Mremnant, baryonic : double
//      Baryonic mass of remnant in Msol

     
//     double fb = fallbackFraction(COCoreMass);
//     double Mproto = COCoreToProtoMass(COCoreMass);
//     return Mproto + fb*(MpreSN - Mproto);
// }
