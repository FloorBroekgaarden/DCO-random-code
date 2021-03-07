#include "pulsarAccretionInducedEvolution.h"



//The function to calculate isolated spin down rate
double calculateSpinDownRate(double  omega, double momentOfInteria, double magneticField, double neutronStarRadius, double alpha){
   /*
    * Calculates the spin down rate for isolated neutron stars in SI
    *
    * See Equation X in arxiv.xxxx
    *
    * Parameters
    * --------------
    *
    *  omega : double
    *     angular velocity in SI
    *  momentOfInteria : double
    *     Moment of Interia of the neutron star calculated using equation of state  in kg m^2
    *  magneticField : double
    *     Magnetic field in Tesla
    *  neutronStarRadius : double
    *     Radius of the neutron star in metres
    *  alpha : double
    *     Angle between the spin and magnetic axes in radians
    *
    * Returns
    * ------------
    *  Spin down rate of an isolated neutron star in SI
    *
    */
   
   double omegaDotTop = 8.0 * pi * pow(omega, 3.0)* pow(neutronStarRadius, 6.0) * pow(magneticField, 2.0)* 1.0 ; //pow ( sin (alpha*0.0174533), 2.0);
   double omegaDotBottom = 3.0* mu_0 * pow(c, 3.0) * momentOfInteria;
   double omegaDot = - omegaDotTop / omegaDotBottom;
   return (omegaDot);
}


//The function to calculate isolated omega
double isolatedOmega(double  omega, double momentOfInteria, double initialMagneticField, double finalMagneticField, double neutronStarRadius, double alpha, double magneticFieldLowerLimit, double stepsize,double tau){
    /*
     * Calculates the spin down rate for isolated neutron stars in SI
     *
     * See Equation X in arxiv.xxxx
     *
     * Parameters
     * --------------
     *
     *  omega : double
     *     angular velocity in SI
     *  momentOfInteria : double
     *     Moment of Interia of the neutron star calculated using equation of state  in kg m^2
     *  magneticField : double
     *     Magnetic field in Tesla
     *  neutronStarRadius : double
     *     Radius of the neutron star in metres
     *  alpha : double
     *     Angle between the spin and magnetic axes in radians
     *
     * Returns
     * ------------
     *  Spin down rate of an isolated neutron star in SI
     *
     */
    
    double constantTop = 8.0 * pi * pow(neutronStarRadius, 6.0) * 1.0 ;      //pow ( sin(alpha*0.0174533), 2.0);
    double constantBottom = 3.0* mu_0 * pow(c, 3.0) * momentOfInteria;
    double constant = constantTop / constantBottom;
    double oneOverOmegaSquaredTermOne =  pow(magneticFieldLowerLimit, 2.0) * stepsize   ;
    double oneOverOmegaSquaredTermTwo =  tau * magneticFieldLowerLimit * (finalMagneticField - initialMagneticField) ;
    double oneOverOmegaSquaredTermThree = (tau/2.0) *  (pow(finalMagneticField, 2.0) - pow(initialMagneticField, 2.0)) ;
    double oneOverOmegaSquared = ((constant * 2.0) * ( oneOverOmegaSquaredTermOne - oneOverOmegaSquaredTermTwo - oneOverOmegaSquaredTermThree ) + 1.0/(pow(omega,2.0))) ;
    omega = pow(oneOverOmegaSquared, -0.5);
    return ( omega );
}


// The function to calculate isolated decay of magnnetic field of the neutron star
double isolatedMagneticFieldDecayFunction(double stepsize, double tau, double magneticFieldLowerLimit, double initialMagneticField){
    /*
     *  Calculates accretion induced magnetic field decay for an accreting neutron star
     *
     * See Equation X in arxiv.xxxx
     *
     * Parameters
     * -------------
     *
     *   tau : double
     *     exponenetial scaled down factor in seconds
     *   magneticFiledLowerLimit : double
     *     Lower cut-off limit of magnetic filed for a neutron star in Tesla
     *   initialMagneticField : double
     *     Initial magnetic field of the neutron star in Tesla
     *
     *  Returns
     * -----------
     * Magnetic filed of a neutron star due to accretion in Tesla
     */
    double isolatedMagneticFieldDecay = magneticFieldLowerLimit + exp(-stepsize/tau)*(initialMagneticField - magneticFieldLowerLimit);
    return (isolatedMagneticFieldDecay);
}



//The function to calculate Alfven Radius
double alfvenRadiusFunction(double mass, double mdot, double radius, double magneticField){
   /*
    *  Calculates the Alfven radius for an accreting neutron star
    *
    *  See Equation X in arxiv.xxxxxx
    *
    *  Parameters
    * -----------
    *  mass : double
    *      Mass of the neutron star in kg
    *  modt : double
    *      Mass transfer rate from the secondary in kg/s
    *  radius : double
    *      Radius of the neutron star in m
    *  magneticField : double
    *      Magnetic field of the neutron star in Tesla
    *
    *  Returns
    * ----------
    *  alfvenRadius : double
    *      Alfven radius in m
    */
   double p = ((pow(radius,6.0))/(pow(mass,0.5)*mdot));
   double q = pow(p, 2.0/7.0 );
   double alfvenRadiusWithoutMagneticField = constPartOfAlfvenRadius*q;
   double alfvenRadius = alfvenRadiusWithoutMagneticField * pow(magneticField,4.0/7.0);
   return (q);
}


// The function to calculate accretion induce decay of magnnetic field of the neutron star
double accretionInducedMagneticFieldDecayFunction(double massGainPerTimeStep, double kappa, double magneticFieldLowerLimit, double magneticField){
   /*
    *  Calculates accretion induced magnetic field decay for an accreting neutron star
    *
    * See Equation X in arxiv.xxxx
    *
    * Parameters
    * -------------
    *  massGainPerTimeStep : double
    *          Mass loss from the secondary for each iteration rate  in kg
    *  kappa : double
    *          Magnetic filed decay scale
    *  magneticFiledLowerLimit : double
    *          Lower cut-off limit of magnetic filed for a neutron star in Tesla
    *  magneticField : double
    *          Previuos magnetic field in Tesla
    *  Msol : constant
    *          Solar Mass in kg
    *  Returns
    * -----------
    * Magnetic filed of a neutron star due to accretion in Tesla
    */
   double accretionInducedMagneticField = magneticFieldLowerLimit + exp(-(1.0/kappa)*massGainPerTimeStep)*(magneticField-magneticFieldLowerLimit);
   //std::cout << "kappa = " << kappa << std::endl;
   // old exp(-(1.0/kappa)*massGainPerTimeStep/Msol)
  // std::cout << "argument of exponent = " << -(1.0/kappa)*massGainPerTimeStep << std::endl;
  // std::cout << "exponential term = " << exp(-(1.0/kappa)*massGainPerTimeStep) << std::endl;
  // std::cout << "part of acc ind mag fld   "  <<exp(-(1.0/kappa)*massGainPerTimeStep)*(magneticField-magneticFieldLowerLimit) << std::endl;
   return (accretionInducedMagneticField);
}


//The function to calculate the angular momentum of an accreting neutron star 
double angularMomentumDueToAccretionFunction (double epsilon, double alfvenRadius, double massGainPerTimeStep, double velocityDifference){
   /*
    *  Calculates the angular momentum for an accreting neutron star
    *  
    *  See Equation X in arxiv.xxxxxx
    *    
    *  Parameters
    *  -----------
    *   epsilon : double
    *         Uncertainty due to mass loss 
    *   alfvenRadius : double
    *         Alfven radius of the neutron star in m
    *   massGainPerTimeStep : double 
    *         Mass loss from the secondary for each iteration rate  in kg
    *   velocityDiffernce : double 
    *         Differnce in the keplerian angular velocity and surface angular velocity of the neutron star in m
    *  Returns
    *   --------
    *   angular momentum of the neutron star due to accretion in kg m^2 sec^-1
    *                      
    *
    */
      double angularMomentumDueToAccretion = 0.25*epsilon*pow(alfvenRadius,2.0)*massGainPerTimeStep*velocityDifference ;
      return (angularMomentumDueToAccretion);
}


// The function to calculate the angular velocity differnce 
double velocityDifferenceFunction (double surfaceAngularVelocity, double mass, double alfvenRadius){
   /*
    * Differnce in the keplerian angular velocity and surface angular velocity of the neutron star in m
    *
    * See Equation X in arxiv.xxxxxx
    *
    * Parameters
    * -------------
    *  surfaceAngularVelocity : double
    *          Surface angular velocity of the neutron star in SI (same as the "omega" we calculate)
    *  mass : double
    
    *  alfvenRadius : double
    *          Alfven radius of the neutron star in m
    *  Returns
    * -----------
    *   Differnce in the keplerian angular velocity and surface angular velocity of the neutron star in SI
    *
    */
  double keplarianVelocityAtAlfvenRadius = sqrt(G*mass)/sqrt(alfvenRadius/2.0);
  double keplarianAngularVelocityAtAlfvenRadius = 2.0* (keplarianVelocityAtAlfvenRadius/alfvenRadius);
  double velocityDifference = keplarianAngularVelocityAtAlfvenRadius - surfaceAngularVelocity ;
  return (velocityDifference);
}
//

//The function to update the magnetic field and spins of neutron stars
void updateMagneticFieldAndSpin(int booleanCheckForAccretion, bool commonEnvelopeFlag, double stepsize, double momentOfInertia, double mass,  double radius, double massGainPerTimeStep, double kappa, double epsilon, double tau, double magneticFieldLowerLimit, double &angularMomentum, double &omega, double &magneticField, double &omegaOverTime, double alpha){
   /*
    *  Description goes here
    *
    * Parameters
    * -----------
    *
    *  booleanCheckForAccretion : integer
    *     returns 1 when there is accretion, returns 0 when there is no accretion
    *  commonEnvelopeFlag : bool
    *     True if there is a common envelope
    *  stepsize : double
    *     timestep size for integartion in seconds
    *  momentOfInertia : double
    *     Moment of Inertia of the neutron star calculated using equation of state  in kg m^2
    *  mass : double
    *     Mass of the neutron star in kg
    *  modt : double
    *     Mass transfer rate from the secondary in kg/s
    *  radius : double
    *     Radius of the neutron star in m
    *  massGainPerTimeStep : double
    *     Mass loss from the secondary for each iteration rate  in kg
    *  kappa : double
    *     Magnetic filed decay scale
    *  epsilon : double
    *     Uncertainty due to mass loss
    *  tau : double
    *     exponenetial scaled down factor in seconds
    *  magneticFiledLowerLimit : double
    *     Lower cut-off limit of magnetic filed for a neutron star in Tesla
    *  angular momentum : double 
    *     Angular momentum of the neutron star in kg m^2 sec^-1
    *  magneticField : double
    *     Previuos magnetic field in Tesla
    *  omegaOverTime : double
    *     Rate of change of angular velocity in
    *  alpha : double
    *     Angle between the spin and magnetic axes in radians
    * Returns
    * ---------
    *  Angular Momentum of a neutron star
    */
    //stepsize = stepsize/5.0 ;
    double initialMagneticField = magneticField ;

    
    //Check point to decide if the neutron star is accreting
    if(booleanCheckForAccretion==false and commonEnvelopeFlag == false){
        
        
        magneticField = isolatedMagneticFieldDecayFunction(stepsize, tau, magneticFieldLowerLimit, initialMagneticField) ;
        omega = isolatedOmega(omega, momentOfInertia, initialMagneticField, magneticField, radius, alpha, magneticFieldLowerLimit, stepsize, tau);
        omegaOverTime = calculateSpinDownRate(omega, momentOfInertia, magneticField, radius, alpha) ; 
        angularMomentum = omega * momentOfInertia ;
        
   }
     else if((booleanCheckForAccretion==true) and (massGainPerTimeStep > 0.0 ))  {  

       // std:: cout << "booleanCheckForAccretion "  << booleanCheckForAccretion << std:: endl;
	//std::cout << "massGainPerTimestep " << massGainPerTimeStep << std::endl;        
        double mdot = massGainPerTimeStep/stepsize ;
        double alfvenRadius = alfvenRadiusFunction(mass, mdot, radius, magneticField);
        double velocityDifference = velocityDifferenceFunction(omega, mass, alfvenRadius);
        magneticField = accretionInducedMagneticFieldDecayFunction(massGainPerTimeStep, kappa, magneticFieldLowerLimit, magneticField) ;
        double deltaAngularMomentum = angularMomentumDueToAccretionFunction(epsilon, alfvenRadius, massGainPerTimeStep, velocityDifference);
        angularMomentum = angularMomentum + deltaAngularMomentum ;
        omega = angularMomentum/momentOfInertia;
        omegaOverTime = (deltaAngularMomentum/stepsize)/momentOfInertia ;
        
             
   }
 
     else if((commonEnvelopeFlag==true) and (massGainPerTimeStep > 0.0 ))  {

        //std:: cout << "commonEnvelopeFlag "  << commonEnvelopeFlag << std:: endl;
	//std::cout << "massGainPerTimestep " << massGainPerTimeStep << std::endl;
        double mdot = massGainPerTimeStep/stepsize ;
        double alfvenRadius = alfvenRadiusFunction(mass, mdot, radius, magneticField);
        double velocityDifference = velocityDifferenceFunction(omega, mass, alfvenRadius);
        //std::cout << "B before acc  " << magneticField << std::endl;
        magneticField = accretionInducedMagneticFieldDecayFunction(massGainPerTimeStep, kappa, magneticFieldLowerLimit, magneticField) ;
        //std::cout << "B after acc  " << magneticField << std::endl;
	//std::cout << "choose whether to allow spin up during CE. Currently allows by default" << std::endl;
        double deltaAngularMomentum = angularMomentumDueToAccretionFunction(epsilon, alfvenRadius, massGainPerTimeStep, velocityDifference);
        angularMomentum = angularMomentum + deltaAngularMomentum ;
        //std::cout << "omega before CE = " << omega << std::endl;
        omega = angularMomentum/momentOfInertia;
        omegaOverTime = (deltaAngularMomentum/stepsize)/momentOfInertia ;
	//std::cout << "omega after CE = " << omega << std::endl;


   }

    
    else {

        //std:: cout << "something is wrong "  << std:: endl;
        //std:: cout << "massGainPerTimeStep "  << massGainPerTimeStep << std:: endl;
        //std:: cout << "booleanCheckForAccretion "  << booleanCheckForAccretion << std:: endl;
   }
      


}






