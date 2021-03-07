//
//  massLoss.cpp
//
//  This file includes implementations of:
//  1) The original stellar wind prescription from Hurley et al 2000 (referred to as 'Hurley' winds, or 'old' winds in the LSC astro paper)
//  2) Wind prescriptions presented in Vink et al 2001 and Belczynski et al 2010 (referred to as 'Vink' winds, or 'new' winds in the LSC astro paper)
//

#include "massLoss.h"

// Mass loss
double massLossKudritzkiReimers(double Mass, double Luminosity, double Radius){
    /*
     Calculate mass loss rate on the GB and beyond based on a prescription taken from Kudritzki and Reimers 1978
     
     Given by Equation 106 in Hurley et al 2000
     
     Parameters
     -----------
     Mass : double
        Mass in Msol
     Luminosity : double
        Luminosity in Lsol
     Radius : double
        Radius in Rsol
     
     Returns
     --------
     Mdot_KR :
        Kudritzki and Reimers mass loss rate (in Msol yr^{-1})

     */
    return 4E-13 * (massLossEta * Luminosity * Radius / Mass); // Shouldn't be eta squared like in paper!
}

double massLossVassiliadisWood(double Mass, double Radius, double Luminosity){
    /*
     Calculate the mass loss on the AGB based on the Mira pulsation period (P0) from Vassiliadis and Wood 1993
     
     Given after Equation 106 in Hurley et al 2000
     
     Parameters
     -----------
     Mass : double
        Mass in Msol
     Radius : double
        Radius in Rsol
     Luminosity : double
        Luminosity in Lsol
     
     Returns
     --------
     Mdot_VW : double
        Mass loss rate on AGB in Msol per year

     */
    double logP0 = std::min(3.3, (-2.07 - 0.9*log10(Mass) + 1.94*log10(Radius)));
    double P0 = pow(10.0, (logP0));
    // In their fortran code, Hurley et al take P0 to be min(p0, 2000.0), implemented here as a minimum power
    double logMdot_VW = -11.4 + 0.0125 * (P0 - 100.0 * std::max((Mass - 2.5), 0.0));
    double Mdot_VW = pow(10.0, (logMdot_VW));
    // DEBUGGING
    //std::cout << "MdotVW, 1.36E-9*L = " << Mdot_VW << " " << 1.36E-9 * Luminosity << std::endl;
    return std::min(Mdot_VW, 1.36E-9*Luminosity);
}

double massLossNieuwenhuijzenDeJager(double Mass, double Radius, double Luminosity, double Metallicity){
    /*
     Calculate the mass loss rate for massive stars (L > 4000 L_sol) using the
     Nieuwenhuijzen & de Jager 1990 prescription, modified by a metallicity
     dependent factor (Kudritzki et al 1989).
     
     Given by Equation in Hurley et al 2000
     
     Parameters
     ----------
     Mass : double
        Mass in Msol
     Radius : double
        Radius in Rsol
     Luminosity : double
        Luminosity in Lsol
     Metallicity : double
        Metallicity Z (Z = 0.02 = Zsol)
     
     Returns
     --------
     Mdot_NdJ : double
        Nieuwenhuijzen & de Jager mass loss rate for massive stars (in Msol yr^-1)

     */
    double smoothTaper = std::min(1.0, (Luminosity-4000.0)/500.0); // Smooth taper between no mass loss and mass loss
    
    return sqrt((Metallicity / Zsol)) * smoothTaper * 9.6E-15 * pow(Radius, 0.81) * pow(Luminosity, 1.24) * pow(Mass, 0.16);
}

double massLossWolfRayetLike(double Luminosity, double Mu){
    /*
     Calculate the Wolf-Rayet like mass loss rate for small hydrogen-envelope mass (when mu < 1.0), taken from (Hamann, Koesterke & Wessolowski 1995, Hamann & Koesterke 1998).
     
     Given after Equation 106 in Hurley et al 2000
     
     Note that the reduction of this formula is imposed in order to match the observed number of black holes in binaries (Hurley et al 2000)
     
     Parameters
     ------------
     Luminosity : double
        Luminosity in Lsol
     Mu : double
        Small envelope parameter (see Equations 97,98)
     
     Returns
     ---------
     Mdot_WR : double
        Mass loss rate (in Msol yr^{-1})

     */
    // In the fortran code there is a parameter here hewind which by default is 1.0 - can be set to zero to disable this particular part of winds. We instead opt for all winds on or off.
    return 1E-13 * pow(Luminosity, 1.5) * (1.0 - Mu);
}

double massLossWolfRayet2(double Luminosity, double Metallicity, double Mu, double wolfRayetFactor){
    /*
     Calculate the Wolf-Rayet like mass loss rate for small hydrogen-envelope mass (when mu < 1.0), taken from (Hamann, Koesterke & Wessolowski 1995, Hamann & Koesterke 1998).
     
     Given by Equation 9 in Belczynski 2010
     
     Note that the reduction of this formula is imposed in order to match the observed number of black holes in binaries (Hurley et al 2000)
     
     Parameters
     ------------
     Luminosity : double
        Luminosity in Lsol
     Metallicity : double
        Metallicity Z (Z = 0.02 = Zsol)
     Mu : double
        Small envelope parameter (see Equations 97,98)
     wolfRayetFactor : double
        Multiplicative constant
     
     Returns
     ---------
     Mdot_WR : double
        Mass loss rate (in Msol yr^{-1})

     */
    // I think StarTrack may still do something different here, there are references to Hamann & Koesterke 1998 and Vink and de Koter 2005
    return wolfRayetFactor * 1E-13 * pow(Luminosity, 1.5) * pow(Metallicity/Zsol, 0.86) * (1.0 - Mu);
}

double massLossWolfRayet3(double Mass){
    /*
     Calculate the Wolf-Rayet like mass loss rate independent of WR star composition as given by Nugis & Lamers 2000
     
     Given by Equation 10 in Belczynski 2010. We do not use this equation by default.
     
     Parameters
     ------------
     Mass : double
        Mass in Msol
     
     Returns
     ---------
     Mdot_WR : double
        Mass loss rate (in Msol yr^{-1})

     */
    double log_Mdot_WR = -5.73 + 0.88*log(Mass);
    return exp(log_Mdot_WR);
}

double massLossLBV(double Radius, double Luminosity){
    /*
     Calculate LBV-like mass loss rate for stars beyond the Humphreys-Davidson limit (Humphreys & Davidson 1994)
     
     Parameters
     -----------
     Radius : double
        Radius in Rsol
     Luminosity : double
        Luminosity in Lsol
     
     Returns
     --------
     Mdot_LBV : double
        LBV-like mass loss rate (in Msol yr^{-1})
     */
    return 0.1 * pow(((1E-5 * Radius * sqrt(Luminosity)) - 1.0), (3.0)) * ((Luminosity / 6E5) - 1.0);
}

double massLossLBV2(double flbv){
    /*
     Calculate LBV-like mass loss rate for stars beyond the Humphreys-Davidson limit (Humphreys & Davidson 1994)
     
     Given by Equation 8 in Belczynski et al 2010
     
     Parameters
     -----------
     flbv : double
        Multiplicitive constant multiplying base rate of 1E-4 Msol yr^-1 for LBV mass loss (Belczynski et al 2010 sets this as 1.5, Mennekens & Vanbeveren 2014 set this as 10)
     
     Returns
     --------
     Mdot_LBV2 : double
        LBV-like mass loss rate (in Msol yr^{-1})

     */
    return flbv * 1E-4;
}

double logDensityAt50PercentVinf(double Z){
    /*
     Calculate the density at the point where v = 0.5 v_inf
     
     Given by Equation 14 in Vink et al 2001
     
     Parameters
     -----------
     Z : double
        Metallicity in Zsol
     
     Returns
     --------
     log_rho : double
        log of the density at the point where v = 0.5 v_inf

     */
    return -13.636 + 0.889 * log(Z);
}

double Tjump(double Z){
    /*
     Calculate the bi-stability jump temperature for a given metallicity
     
     Given by Equation 15 in Vink et al 2001
     
     Parameters
     -------------
     Z : double
        Metallicity in Zsol
     
     Returns
     ---------
     Tjump : double
        Jump temperature in K

     */
    return (61.2 + 2.59 * logDensityAt50PercentVinf(Z)) * 1E3;
}

double massLossOB(double Mass, double Radius, double Luminosity, double Teff, double Metallicity){
    /*
     Calculate mass loss rate for massive OB stars using the Vink et al 2001 prescription
     
     Equations 24 & 25 in Vink et al 2001, Equations 6 & 7 in Belczynski et al 2010
     
     Parameters
     -----------
     M : float
     Mass in Msol
     R : float
     Radius in Rsol
     L : float
     Luminosity in Lsol
     Teff : float
     Effective temperature in K
     Z : float
     Metallicity Z (Z = 0.02 = Zsol)
     
     Returns
     --------
     Mdot_OB : float
     Mass loss rate for hot OB stars in Msol yr^{-1}

     */
    // Declare variables
    double  V = 0;
    double  logMdotOB = 0;
    bool    debugging = false;
    
    if(debugging){
    std::cout << "Teff: " << Teff << std::endl;
    std::cout << "Luminosity: " << Luminosity << std::endl;
    std::cout << "Mass: " << Mass << std::endl;
    std::cout << "Z: " << Metallicity << std::endl;
    std::cout << "Test: " << 1E5*pow(10.0,-5.0) << std::endl;
    }
    
    if(Teff < 12500.0){

        if(debugging){std::cout << "Teff < 12500.0" << std::endl;}
        
        // Probably should just return 0?
        std::cerr << "Error : too cold - shouldn't get here" << std::endl;
        return 0;
    }
    else if(Teff >= 12500.0 and Teff <= 25000.0){
        if(debugging){std::cout << "Teff >= 12500.0 and Teff <= 25000.0" << std::endl;}
        V = 1.3; // v_inf/v_esc
        logMdotOB = -6.688 + 2.210*log10(Luminosity/1E5) - 1.339*log10(Mass/30.0) - 1.601*log10(V/2.0) + 0.85*log10(Metallicity/Zsol) + 1.07*log10(Teff/20000.0);
        return pow(10.0, logMdotOB);
    }
    else if(Teff >25000){
        if(debugging){std::cout << "Teff >25000" << std::endl;}
        if(Teff > 50000.0){
            if(debugging){
                std::cout << "Teff > 50000.0" << std::endl;
                std::cout << "Warning, winds (massLossOB) being used outside their comfort zone" << std::endl;
            }
        }
        V = 2.6; // v_inf/v_esc
        logMdotOB = -6.697 + 2.194*log10(Luminosity/1E5) - 1.313*log10(Mass/30.0) - 1.226*log10(V/2.0) + 0.85*log10(Metallicity/Zsol) + 0.933*log10(Teff/40000.0) - (10.92*log10(Teff/40000.0)*log10(Teff/40000.0));
        return pow(10.0, logMdotOB);
    }
    else{
        std::cerr << "Error: shouldn't be able to get here" << std::endl;
        return 0;
    }
}

double massLossRate(double Mass, double Luminosity, double Radius, double Mu, double Metallicity, int STELLAR_TYPE){
    /*
     Calculate the dominant mass loss mechanism and associated rate for the star at the current evolutionary phase
     
     Parameters
     -----------
     Mass : double          Mass in Msol
        Luminosity : double    Luminosity in Lsol
     Radius : double        Radius in Rsol
        Mu : double            Small envelope parameter
     Metallicity : double   Metallicity Z (Z = 0.02 = Zsol)
     
     STELLAR_TYPE : int     Phase of stellar evolution
     
     Returns
     --------
     Mdot : double          Mass loss rate in Msol per year

     */
    bool    debugging = false;
    double dms = 0.0;
    double dml = 0.0;
    double dmt = 0.0;
    
    if(debugging){std::cout << "Mass loss rate hurley " << std::endl;}
    
    if(Luminosity > 4E3 or STELLAR_TYPE == MS_MORE_THAN_07){
        dms = massLossNieuwenhuijzenDeJager(Mass, Radius, Luminosity, Metallicity);
    }
    else{
        dms = 0.0;
    }
    
    // If stellar type between 2 and 9
    if(STELLAR_TYPE >= HERTZSPRUNG_GAP and STELLAR_TYPE <= NAKED_HELIUM_STAR_GIANT_BRANCH){
        
        dml = massLossKudritzkiReimers(Mass, Luminosity, Radius);
        
        if(STELLAR_TYPE == EARLY_ASYMPTOTIC_GIANT_BRANCH or STELLAR_TYPE == THERMALLY_PULSING_ASYMPTOTIC_GIANT_BRANCH){
            
            dmt = massLossVassiliadisWood(Mass, Radius, Luminosity);
            
            double dms2 = std::max(dml, dmt);
            dms = std::max(dms, dms2);
        }
        else if(STELLAR_TYPE >= NAKED_HELIUM_STAR_MS){
            
            Mu = 0.0;
            dms = std::max(dml, massLossWolfRayetLike(Luminosity, Mu));
            
        }
        else{
            
            dms = std::max(dml, dms);
            
            if(Mu < 1.0){
                
                dml = massLossWolfRayetLike(Luminosity, Mu);
                
                dms = std::max(dml, dms);
                
            }
            
            double x = 1E-5 * Radius * sqrt(Luminosity);
            
            if(Luminosity > 6.0E5 and x > 1.0){
                
                dml = massLossLBV(Radius, Luminosity);
                
                dms = dms + dml;
                
            }
            
        }
        
    }
    
    // I think BSE/StarTrack have a multiplier factor that they use here.
    return dms;
}


// Based on implementation in StarTrack
double massLossRateVink(double Mass, double Luminosity, double Radius, double Mu, double Metallicity, double Teff, double flbv, double wolfRayetFactor, int STELLAR_TYPE, bool & LBVphaseFlag){
    /*
     Calculate the dominant mass loss mechanism and associated rate for the star at the current evolutionary phase
     
     Parameters
     -----------
     Mass : double
     Mass in Msol
     Luminosity : double
     Luminosity in Lsol
     Radius : double
     Radius in Rsol
     Mu : double
     Small envelope parameter
     Metallicity : double
     Metallicity Z (Z = 0.02 = Zsol)
     Teff : double
     Effective temperature in K
     flbv : double
     LBV mass loss multiplicative constant
     wolfRayetFactor : double
     WR mass loss multiplicative constant
     STELLAR_TYPE : int
     Phase of stellar evolution
     
     Returns
     --------
     Mdot : double
     Mass loss rate in Msol per year

     */
    
    double  dms = 0.0;
    double  dml = 0.0;
    bool    debugging = false;
    //debugging = true;
    //double dmt = 0.0;
    
    // What units is metallicity in? As I suspercted, is absolute
    //std::cout << "Metallicity = " << Metallicity << std::endl;
    
    double x = 1E-5 * Radius * sqrt(Luminosity);
    
    double  LBV_luminosity_limit = LBV_LUMINOSITY_LIMIT_STARTRACK;
//    double  LBV_luminosity_limit = LBV_LUMINOSITY_LIMIT_VANBEVEREN;
    
    // Investigating LBV mass loss
    if(debugging){
        std::cout << "Radius: " << Radius << std::endl;
        std::cout << "Luminosity: " << Luminosity << std::endl;
        std::cout << "x: " << x << std::endl;
    }
    
    Teff=Teff*Tsol; // Change Teff to Kelvin so it can be compared with values as stated in Vink prescription
    
    if(debugging){std::cout << "Mass loss rate vink " << std::endl;}
    
    if(STELLAR_TYPE >= MS_LESS_THAN_07 and STELLAR_TYPE < NAKED_HELIUM_STAR_MS and Teff < 12500.0){ // Cool stars, use Hurley et al 2000 winds

        double thisMassLossRate = massLossRate(Mass, Luminosity, Radius, Mu, Metallicity, STELLAR_TYPE);

        if(debugging){ 
            std::cout << "Hurley mass loss rate " << thisMassLossRate << std::endl;
        }

        dms = thisMassLossRate;
    }
    else if(STELLAR_TYPE >= MS_LESS_THAN_07 and STELLAR_TYPE < NAKED_HELIUM_STAR_MS and Teff >= 12500.0){ // Hot stars, use Vink et al 2001 winds (ignoring bistability jump)
        if(debugging){  std::cout << "massLossOB" << std::endl;}

        double thisMassLossRateOB = massLossOB(Mass, Radius, Luminosity, Teff, Metallicity);

        if(debugging){
            std::cout << "massLossOB " << thisMassLossRateOB << std::endl;
        }

        dms = thisMassLossRateOB;
    }
    else if(STELLAR_TYPE == NAKED_HELIUM_STAR_MS or STELLAR_TYPE == NAKED_HELIUM_STAR_HETZSPRUNG_GAP or STELLAR_TYPE == NAKED_HELIUM_STAR_GIANT_BRANCH){
        Mu = 0.0;

        double thisMassLossRateWR = massLossWolfRayet2(Luminosity, Metallicity, Mu, wolfRayetFactor);
        
        if(debugging){
            std::cout << "massLossWolfRayet2 " << thisMassLossRateWR << std::endl;
        }

        dms = thisMassLossRateWR;
    }
    else if(STELLAR_TYPE >= HELIUM_WHITE_DWARF){

        // No mass loss from remnants

        dms = 0.0;

    }
    else{

        std::cerr << "Error in massLossRateVink, shouldn't get here. No valid winds for stellar type:\t" << STELLAR_TYPE << std::endl;

    }
    
    if(STELLAR_TYPE >= MS_LESS_THAN_07 and STELLAR_TYPE < NAKED_HELIUM_STAR_MS and Luminosity > LBV_luminosity_limit  and x > 1.0){ // Luminous blue variable
		LBVphaseFlag = true;
		if(debugging){
			std::cout << "A star is going through the LBV phase. Mass, Stellartype: " << Mass << " " << STELLAR_TYPE << std::endl;
		}		
        double thisMassLossRateLBV = massLossLBV2(flbv);

        if(debugging){
            std::cout << "massLossLBV2: " << thisMassLossRateLBV << std::endl;
        }

        dms = thisMassLossRateLBV;

    }


    // BSE and StarTrack have some mulptilier they apply here
    return dms;

}

/////////////
//double massLossRateVink(double Mass, double Luminosity, double Radius, double Mu, double Metallicity, double Teff, int STELLAR_TYPE){
//    /*
//     Calculate the dominant mass loss mechanism and associated rate for the star at the current evolutionary phase
//     
//     Parameters
//     -----------
//     Mass : double
//     Mass in Msol
//     Luminosity : double
//     Luminosity in Lsol
//     Radius : double
//     Radius in Rsol
//     Mu : double
//     Small envelope parameter
//     Metallicity : double
//     Metallicity Z (Z = 0.02 = Zsol)
//     Teff : double
//     Effective temperature in K
//     STELLAR_TYPE : int
//     Phase of stellar evolution
//     
//     Returns
//     --------
//     Mdot : double
//     Mass loss rate in Msol per year

//     */
//    
//    double dms = 0.0;
//    double dml = 0.0;
//    double dmt = 0.0;
//    
//    if(STELLAR_TYPE == MS_LESS_THAN_07){
//        dml = massLossRate(Mass, Luminosity, Radius, Mu, Metallicity, STELLAR_TYPE);
//    }
//    else if(STELLAR_TYPE == MS_MORE_THAN_07){
//        dml = massLossOB(Mass, Radius, Luminosity, Teff, Metallicity);
//        
//        dms = dml;
//    }
//    else if(STELLAR_TYPE >= HERTZSPRUNG_GAP and STELLAR_TYPE <= NAKED_HELIUM_STAR_GIANT_BRANCH){
//        
//        dml = massLossKudritzkiReimers(Mass, Luminosity, Radius);
//        
//        if(STELLAR_TYPE == EARLY_ASYMPTOTIC_GIANT_BRANCH or STELLAR_TYPE == THERMALLY_PULSING_ASYMPTOTIC_GIANT_BRANCH){
//            
//            dmt = massLossVassiliadisWood(Mass, Radius, Luminosity);
//            
//            dms = std::max(dml, dmt);
//            
//        }
//        else if(STELLAR_TYPE >= NAKED_HELIUM_STAR_MS){
//            
//            Mu = 0.0;
//            dms = std::max(dml,  massLossWolfRayet2(Luminosity, Metallicity, Mu));//massLossWolfRayetLike(Luminosity, Mu));
//            
//        }
//        else{
//            
//            dms = std::max(dml, dms);
//            
//            if(Mu < 1.0){
//                
//                dml = massLossWolfRayet2(Luminosity, Metallicity, Mu);
//                
//                dms = std::max(dml, dms);
//                
//            }
//            
//            double x = 1E-5 * Radius * sqrt(Luminosity);
//            
//            if(Luminosity > 6.0E5 and x > 1.0){
//                
//                dml = massLossLBV2();
//                
//                dms = dml;
//                
//            }
//            
//        }
//        
//    }
//    
//    return dms;
//}

// probably using wrong log here
//double massLossOB(double Mass, double Radius, double Luminosity, double Teff, double Z){
//    /*
//     Calculate mass loss rate for massive OB stars using the Vink et al 2001 prescription
//     
//     Equations 24 & 25 in Vink et al 2001, Equations 6 & 7 in Belczynski et al 2010
//     
//     Parameters
//     -----------
//     M : float
//     Mass in Msol
//     R : float
//     Radius in Rsol
//     L : float
//     Luminosity in Lsol
//     Teff : float
//     Effective temperature in K
//     Z : float
//     Metallicity in Zsol
//     
//     Returns
//     --------
//     Mdot_OB : float
//     Mass loss rate for hot OB stars in Msol yr^{-1}
//     */
//    
//    // Declare variables
//    double V = 0;
//    double logMdotOB = 0;
//    
//    if(Teff < 12500.0){
//        // Probably should just return 0?
//        std::cout << "Error : too cold - shouldn't get here" << std::endl;
//        return 0;
//    }
//    else if(Teff >= 12500.0 and Teff <= 25000.0){
//        V = 1.3; // v_inf/v_esc
//        logMdotOB = -6.688 + 2.210*log(Luminosity/1E5) - 1.339*log(Mass/30.0) - 1.601*log(V/2.0) + 0.85*log(Z) + 1.07*log(Teff/20000.0);
//        return exp(logMdotOB);
//    }
//    else if(Teff >25000){
//        if(Teff > 50000.0){
//            std::cout << "Warning, winds (massLossOB) being used outside their comfort zone" << std::endl;
//        }
//        V = 2.6; // v_inf/v_esc
//        logMdotOB = -6.697 + 2.194*log(Luminosity/1E5) - 1.313*log(Mass/30.0) - 1.226*log(V/2.0) + 0.85*log(Z) + 0.933*log(Teff/40000.0) - 10.92*log(Teff/40000.0)*log(Teff/40000.0);
//        return exp(logMdotOB);
//    }
//    else{
//        std::cout << "Error: shouldn't be able to get here" << std::endl;
//    }
//    
//    return 0;
//}

