//
//  commonEnvelope.cpp
//

#include "commonEnvelope.hpp"

double  doubleCommonEnvelopePhase(double alphaCE, double lambda1, double lambda2, double a, double m1, double mfin1, double menv1, double rRLd1, double m2, double mfin2, double menv2, double rRLd2, int type1, int type2, bool & doubleCommonEnvelopeFlag){
    // Double common envelope phase prescription (Brown 1995) to calculate new semi-major axis due to the CEE, as described in Belczynsky et al. 2002, eq. (12)
    // alphaCE: CE efficiency parameter
    // lambda: parameter to measure "the central concentration of the donor"
    // a: semi-major axis before CEE, in AU units.
    // m1: mass of star 1 before CEE, in Msol units
    // mfin1: mass of star 1 after CEE, in Msol units
    // menv1: mass of star 1 envelope after CEE, in Msol units
    // rRLd1: Roche Lobe radius of star 1, in AU units
    // m2: mass of star 2 before CEE, in Msol units
    // mfin2: mass of star 2 after CEE, in Msol units
    // menv2: mass of star 2 envelope after CEE, in Msol units
    // rRLd2: Roche Lobe radius of star 2, in AU units
    bool    debugging = false;

    double  K1 = (2/(lambda1*alphaCE))*m1*menv1/rRLd1;
    double  K2 = (2/(lambda2*alphaCE))*m2*menv2/rRLd2;
    double  K3 = m1*m2/a;
    double  K4 = mfin1*mfin2;

	// Setting K_i to zero in case the star is a compact object, with envelope = 0.0
	if(type1 >= HELIUM_WHITE_DWARF){
		K1 = 0.0;
	}

	if(type2 >= HELIUM_WHITE_DWARF){
		K2 = 0.0;
	}

    if(debugging){
        std::cout << "lambda1 = " << lambda1 << std::endl;
        std::cout << "lambda2 = " << lambda2 << std::endl;
        std::cout << "alphaCE = " << alphaCE << std::endl;
        std::cout << "K1 = (2/(lambda*alphaCE))*m1*menv1/rRLd1 = " << K1 << std::endl;
        std::cout << "K2 = (2/(lambda*alphaCE))*m2*menv2/rRLd2 = " << K2 << std::endl;
        std::cout << "K3 = m1*m2/a = " << K3 << std::endl;
        std::cout << "mfin1 = " << mfin1 << std::endl;
        std::cout << "mfin2 = " << mfin2 << std::endl;
        std::cout << "K4 = mfin1*mfin2 = " << K4 << std::endl;
        std::cout << "afin = " << K4/(K1+K2+K3) << std::endl;
    }
    
    if((K1>0.0)and(K2>0.0)and(K3>0.0)and(K4>0.0)){
		if(debugging){std::cout<<"Double Common Envelope event."<<std::endl;}
        doubleCommonEnvelopeFlag = true;
    }

    return K4/(K1+K2+K3);
}

double calculateCommonEnvelopeLambdaParameter(const Star &star, const programOptions &options, unsigned long randomSeed){
    /*
    This function calculates the common envelope lambda parameter for a given star
    with a prescription chosen by the user 
    
    Parameters
    -----------
    star : Star
        The star to calculate lambda for
    options : programOptions    
        User specified program options
    randomSeed : unsigned long
        Random seed of the binary

    Returns
    --------
    lambda : double
        Common envelope lambda parameter
    */

    double menv = star.m_Mass - star.m_coreMass;    // Stars envelope mass, in Msol units
    double lambda = 0;

    bool debugging = true;

    if(options.commonEnvelopeLambdaPrescription == LAMBDA_FIXED){

        if(debugging){
            std::cout << "Using LAMBDA_FIXED" << std::endl;
        }

        lambda = options.commonEnvelopeLambda;
    }
    else if(options.commonEnvelopeLambdaPrescription == LAMBDA_LOVERIDGE){
        lambda = lambdaLoveridgeEnergyFormalism(star.m_Metallicity, star.m_MZAMS, star.m_Mass, menv, star.m_COCoreMass, star.m_Radius, false, randomSeed); // mass loss manually set to false here . Could use options.useMassLoss
    }
    else if(options.commonEnvelopeLambdaPrescription == LAMBDA_NANJING){

        if(debugging){
            std::cout << "Using LAMBDA_NANJING" << std::endl;
        }

        lambda = lambdaNanjing(star.m_MZAMS, star.m_Mass, star.m_coreMass, star.m_Radius, star.m_stellarType, star.m_Metallicity, options);
    }
    else if(options.commonEnvelopeLambdaPrescription == LAMBDA_KRUCKOW){
        lambda = lambdaKruckow(star.m_Radius, options.commonEnvelopeSlopeKruckow);
    }
    else if(options.commonEnvelopeLambdaPrescription == LAMBDA_DEWI){
        lambda = lambdaDewi(star.m_Mass, star.m_coreMass, star.m_Radius, star.m_RZAMS, star.m_Luminosity, star.m_stellarType, randomSeed);
    }
    else{
        std::cout << "Unrecognised common envelope lambda prescription. Shouldn't get here. Using lambda = 1.0" << std::endl;
        lambda = 1.0;
    }

    // Multiply calculated lambda by multiplier (default = 1) 
    lambda *= options.commonEnvelopeLambdaMultiplier;

    // If lambda is below some minimum value, set to zero (infinite binding energy) AVG.
    if(lambda < 0.00001){
        lambda = 0.0;
    }

    return lambda;
}



double massAccretionNeutronStar(const gsl_rng *r,  const programOptions &options, double M_comp, double R_comp ){
     /*
 *   documentation goes here
 *
 * */
     double Delta_M_NS = 0;

     if (options.commonEnvelopeMassAccretionPrescription == COMMON_ENVELOPE_ACCRETION_CONSTANT) {
          Delta_M_NS =  options.commonEnvelopeMassAccretionConstant ;  
     } 

     else if (options.commonEnvelopeMassAccretionPrescription == COMMON_ENVELOPE_ACCRETION_ZERO) {
          Delta_M_NS =  0.0 ; 
     }

     else if (options.commonEnvelopeMassAccretionPrescription == COMMON_ENVELOPE_ACCRETION_UNIFORM) {
         
         double rmin = options.commonEnvelopeMassAccretionMin;
         double rmax = options.commonEnvelopeMassAccretionMax; 
         Delta_M_NS = uniformDistribution( r , rmin, rmax) ;  //Oslowski+ (2011)
     }
 
     else if (options.commonEnvelopeMassAccretionPrescription == COMMON_ENVELOPE_ACCRETION_MACLEOD) {
        
         double mm =  -1.0714285714285712e-05   ;// gradient of the linear fit for gradient
         double cm =   0.00012057142857142856   ;// intercept of the linear fit for gradient
         double mc =   0.01588571428571428      ;// gradient of the linear fit for intercept 
         double cc =  -0.15462857142857137      ;// intercept of the linear fir for intercept
         double m = mm*M_comp + cm ; 
         double c = mc*M_comp + cc ;
         Delta_M_NS = m*R_comp  + c  ;                             //linear regression estimated from Macleod+ (2015)
         if (Delta_M_NS > options.commonEnvelopeMassAccretionMax){
             Delta_M_NS = options.commonEnvelopeMassAccretionMax ;
         } 
         else if (Delta_M_NS < options.commonEnvelopeMassAccretionMin){
              Delta_M_NS = options.commonEnvelopeMassAccretionMin ; 
         }
     }

     return Delta_M_NS;
}


