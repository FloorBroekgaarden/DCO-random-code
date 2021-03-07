#include "KickAndOrbitalPlane.h"

#include <gsl/gsl_randist.h> 
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_erf.h>

double kickVelocityBrayEldridge(double Mejecta, double Mremnant, double alpha, double beta){
    /*
    Implements model for kick velocities from Bray & Eldridge 2016,2018:
    https://arxiv.org/abs/1605.09529
    https://arxiv.org/abs/1804.04414

    Parameters
    -----------
    Mejecta : double
        Mass of ejecta in solar masses
    Mremnant : double
        Mass of remnant in solar masses
    alpha : double
        Fitting coefficient
    beta : double
        Fitting coefficient
    
    Returns
    --------
    vkick : double
        Kick velocity in km s^-1
    */
    return alpha * (Mejecta / Mremnant) + beta;
}

double kickVelocityDistribution0(void){
    // Returns a kick velocity of 0 ms^-1
    return 0.0;
}

double maxwellCDF(double x, double sigma){
    /*
    We require inverse sampling from the Maxwell CDF for MCMC and importance sampling

    From https://en.wikipedia.org/wiki/Maxwell–Boltzmann_distribution , the Maxwell CDF is:

    erf(x/(sqrt(2.0)*sigma)) - (sqrt(2.0/pi) * x * exp(-(x*x)/(2.0*sigma*sigma))/sigma)

    where erf is the error function which has the form:

    erf(x) = (2/\sqrt(\pi)) \int_0^x dt \exp(-t^2).

    We use the gsl error function https://www.gnu.org/software/gsl/manual/html_node/Error-Function.html

    Parameters
    -----------
    x : double
    sigma : double

    Returns
    --------
    maxwellCDF : double
    */
    return gsl_sf_erf(x/(sqrt(2.0) * sigma)) - (sqrt(2.0/pi) * x * exp(-(x*x)/(2.0*sigma*sigma))/sigma);
}

struct my_f_params { 
    double y;       // Value of CDF, should be drawn as U(0,1)
    double sigma;   // sigma for kick distribution
};

double functionToRootFind(double x, void * p){
    /*
    Calculate the inverse of the Maxwell CDF
    
    Parameters
    -----------
    x : double
        Value of the kick vk which we want to find
    p : struct my_f_params 
        Structure containing y, the CDF draw U(0,1), and sigma which 
        sets the characteristic kick of the Maxwellian distribution

    Returns
    --------
    z : double
        Should be zero when x = vk, the value of the kick to draw
    */
    struct my_f_params * params = (struct my_f_params *)p;

    double y = (params->y);             // Value of CDF, should be drawn as U(0,1)
    double sigma = (params->sigma);     // sigma for kick distribution

    return maxwellCDF(x, sigma) - y;
}

double kickVelocityDistributionMaxwell(const gsl_rng *r, double sigma, double star_kick_magnitude_random_number){
    /*

    Draw a kick velocity in km s^-1 from a Maxwellian distribution of the form:

    sqrt(2/pi) * (x * x * exp(-(x * x)/(2.0 * a * a)) / (a*a*a))

    Roughly factor of 10 slower than old method, 10^-4 s vs 10^-5 s for old method.

    Implement option to only do it this way if using mcmc or importance sampling, otherwise use old method
    
    Inverse sampling done using root finding in GSL, adapted from the examples here:
    https://www.gnu.org/software/gsl/doc/html/roots.html

    Parameters
    -----------
    r: const gsl_rng
        Used for random number generator
    sigma : double
        Determines characteristic kick velocity in km s^-1 drawn from Maxwellian distribution
    star_kick_magnitude_random_number : double
        Random number between 0 and 1 used for drawing from the inverse CDF of the Maxwellian

    Returns
    --------
    vk : double
        Drawn kick velocity in km s^-1
    */

    bool old_method = false; // options.useMCMC

    if(old_method){

        // Old method
        double rms = 0.;
    
        //generate each component of a 3 dimensional velocity vector, sum the squares (because we need the magnitude)
        for (int i=0; i<3; i++){
            double vel = gsl_ran_gaussian(r,sigma);
            
            rms += vel*vel; 
        }
    
        //return the magnitude
        return std::sqrt(rms);

    }
    else{

        //double star_kick_magnitude_random_number = gsl_rng_uniform(r);

        bool debugging = false;

        if(debugging){
            std::cout << "Drawing kick" << std::endl;
        }

        int status;
        int iter = 0, max_iter = 100;
        const gsl_root_fsolver_type *T;
        gsl_root_fsolver *s;
        double result = 0;
        double x_lo = 0.0, x_hi = 5.0*sigma; // options.kickmax or whatever it is called.
        gsl_function F;

        double maximum_inverse = maxwellCDF(x_hi, sigma);

        if(debugging){
            std::cout << "maximum inverse is " << x_hi << " " << sigma << " " << maximum_inverse << std::endl;
            std::cout << "random number = " << star_kick_magnitude_random_number << std::endl;
        }

        while(star_kick_magnitude_random_number > maximum_inverse){

            x_hi *= 2.0;

            maximum_inverse = maxwellCDF(x_hi, sigma);

            if(debugging){

                std::cout << "random number past end point, increasing end point" << std::endl;
                std::cout << "new maximum inverse is " << maximum_inverse << std::endl;

            }

        }

        if(star_kick_magnitude_random_number > maximum_inverse){

            star_kick_magnitude_random_number = maximum_inverse;

            if(debugging){

                std::cout << "random number past end point, resetting to maximum" << std::endl;

            }

        }
        
        struct my_f_params params = {star_kick_magnitude_random_number, sigma}; // y, sigma

        F.function = &functionToRootFind;
        F.params = &params;

        // gsl_root_fsolver_brent
        // gsl_root_fsolver_bisection
        T = gsl_root_fsolver_brent;
        s = gsl_root_fsolver_alloc (T);

        gsl_root_fsolver_set (s, &F, x_lo, x_hi);

        //std::cout << "GSL_CONTINUE = " << GSL_CONTINUE << std::endl;
        //std::cout << "Status = " << status << std::endl;
        status = GSL_CONTINUE;

        while(status == GSL_CONTINUE and iter < max_iter){

            //std::cout << "in while loop" << std::endl;

            iter++;
            status = gsl_root_fsolver_iterate (s);
            result = gsl_root_fsolver_root (s);
            x_lo = gsl_root_fsolver_x_lower (s);
            x_hi = gsl_root_fsolver_x_upper (s);
            status = gsl_root_test_interval (x_lo, x_hi, 0, 0.001);

            //if (status == GSL_SUCCESS){
            //
            //    std::cout << "Converged:\n" << std::endl;
            //
            //    std::cout << iter << " " << x_lo << " " << x_hi << " " << r << " " << x_hi - x_lo << std::endl;
            //
            //}
        } 

        if(debugging){
            std::cout << "vk = " << result << std::endl;
        }
        
        // De-allocate memory for root solver
        gsl_root_fsolver_free (s);

        return result;

    }
    
}

double kickVelocityDistributionFlat(const gsl_rng *r, double kickVelocityDistributionMaximum, double star_kick_magnitude_random_number){
    /*
    Generates a kick velocity in km s^-1 from a uniform distribution between 0 and kickVelocityDistributionMaximum
    
    Parameters
    -----------
    r : const gsl_rng
        Used for gsl random number generation
    kickVelocityDistributionMaximum : double
        Maximum kick velocity in km s^-1 to draw
    star_kick_magnitude_random_number : double
        Random number uniform between 0 and 1
    
    Returns
    --------
    vk : double
        Drawn kick velocity in km s^-1
    
    */
    bool old_method = false; //true;

    if(old_method){
        double u = gsl_rng_uniform(r);
        return u * kickVelocityDistributionMaximum;
    }
    else{
        return star_kick_magnitude_random_number * kickVelocityDistributionMaximum;
    }
}

void blackHoleKicks(const programOptions &options, double &vk, double fallback, double blackHoleMass, int stellarType){
    /*
     Calculate the kick given to a black hole based on users chosen assumptions about black hole kicks, fallback and the magnitude of a kick drawn from a distribution
     
     Currently we can use:
     'full' : Black holes receive the same kicks as neutron stars
     'reduced' : Black holes receive the same momentum kick as a neutron star, but downweighted by the black hole mass
     'zero' : Black holes receive zero natal kick
     'fallback' : Black holes receive a kick downweighted by the amount of mass falling back onto them
     
     Parameters
     -----------
     options : programOptions
        Class contain user specified options relating to physical assumptions
     vk : double &
        Kick velocity that would otherwise be applied to a neutron star
     fallback : double
        Fraction of mass that falls back onto the proto-compact object
     blackHoleMass : double
        Mass of remnant (in Msol)
     stellarType : int
        Stellar type of remnant (to check if neutron star or black hole)
     
     Returns
     --------
     nothing, updates vk

     */
    
    bool    debugging = false;
    
    if(debugging){std::cout << "Using black hole kicks option " << options.blackHoleKicks << " " << options.blackHoleKicksOption << std::endl;}

    switch(options.blackHoleKicksOption){
            
        case FULL:
            vk = vk;                                        // BH recieves full kick - no adjustment necessary
            break;
        case REDUCED:
            vk *= neutronStarMass/blackHoleMass;            // Kick is reduced by the ratio of the black hole mass to neutron star mass i.e. v_bh = ns/bh  *v_ns
            break;
        case ZERO:
            vk = 0.0;                                       // BH Kicks are set to zero regardless of BH mass or kick velocity drawn.
            break;
        case FALLBACK:
            if(debugging){std::cout << "Reducing BH kicks by fallback = " << fallback << " from " << vk << std::endl;}
            vk *= (1.0 - fallback);                         // Using the so-called 'fallback' prescription for BH kicks
            if(debugging){std::cout << "Reducing BH kicks by fallback = " << fallback << " to " << vk << std::endl;}
            break;
    }
    
}

double radiusFunctionOfPsi(double a=1.0, double e=0.0, double psi=0.0){
    // Calculate orbital radius at some true anomaly psi
    // (Default is a circular orbit e=0, r = a for all psi)
    double r              = 0.0;
    double topFraction    = a*(1.0 - (e*e));
    double bottomFraction = 1.0 + e*cos(psi);
    r = topFraction/bottomFraction;
    return r;
}

double radiusFunctionOfEccentricAnomaly(double a, double e, double E){
    // Calculate radius using eccentric anomaly
    return a*(1.0 - e*cos(E));
}

double orbitalVelocity(double M=Msol, double r=AU, double a=AU){
    // Calculate orbital velocity at some true anomaly psi
    // Default is a circular orbit, V = sqrt(gm/a) = const.
    // Since this equation contains 'G', all other quantities must be in SI
    // to get answer in ms^-1
    double v        = 0.0;
    double vSquared = 0.0;
    double brackets = 0.0;
    brackets        = (2.0/r) - (1.0/a);
    vSquared        = G*M*brackets;
    v               = sqrt(vSquared);
    return v;
}


///// SOME EQUATIONS TO CALCULATE PRE-SUPERNOVA ORBITAL ENERGY ///////
double EPreSupernova(double mu, double M, double a){
    // Vis-Viva equation to calculate E
    return -(G*mu*M)/(2.0*a);
}

double EPreSupernovaLong(double m1, double m2, double r, double a){
    // Using the full expression summing kinetic energy of both stars and potential energy
    // Hopefully the most general expression, will apply in elliptical orbits
    double M = m1 + m2;
    double vRel = sqrt(G*M*((2.0/r)-(1.0/a)));
    double v1 = vRel*(m2/M);
    double v2 = vRel*(m1/M);
    return 0.5*m1*v1*v1 + 0.5*m2*v2*v2 -(G*m1*m2/r);
}

double EPreSupernovaAlt(double mu, double M, double vRel, double r){
    // Calculate the total energy prior to the supernova
    //double E = 0;
    //E = (0.5  *mu  *vRel  *vRel) - ((G  *M  *mu)/r);
    //return E;
    return (0.5  *mu  *vRel  *vRel) - ((G  *M  *mu)/r);
}

/////////////////////////////////////////////////////////////////////

// Calculate post-supernova orbital energy for an initially circular orbit
double EPostSupernovaInitiallyCircular(double a=1.0, double M1=0.5, double M2=0.5, double M1Prime=1.0, double M2Prime=1.0, double theta=0.0, double phi=0.0, double uk=0.0){
    // Calculate the post-supernova energy to see if the binary survives.
    // Note that this is in SI units (i.e. J) so need uK, a, M etc in SI
    // Will need to pass m1prime and m2prime as well as M/Mprime
    double M             = M1 + M2;
    double MPrime        = M1Prime + M2Prime;
    double EPrime        = 0.0;
    double MOverMPrime   = M/MPrime;
    double quadraticTerm = (1.0 + 2.0*uk*cos(theta)*cos(phi) + uk*uk);
    double Eterm         = (-G*M1Prime*M2Prime)/(2.0*a);
    double brackets      = (2.0 - MOverMPrime*quadraticTerm);
    //EPrime               = Eterm  *(2 - MOverMPrime*quadraticTerm);
    EPrime               = Eterm*brackets;
    return EPrime;
}

// Post supernova orbital eccentricity
double ePostSupernovaInitiallyCircular(double M, double MPrime, double theta, double phi, double uk){
    // Calculate ePrime for an initially circular orbit (e=0) using my equation/the one from Brandt & Podsiadlowski 1995 based on angular momentum arguments
    double MOverMPrime   = M/MPrime;
    double quadraticTerm = (1.0 + 2.0*uk*cos(theta)*cos(phi) + uk*uk);
    double aSquaredPlusBSquaredTerm = (1.0 + uk*cos(theta)*cos(phi))*(1.0 + uk*cos(theta)*cos(phi)) + (uk*sin(theta))*(uk*sin(theta));
    double brackets = (2.0 - MOverMPrime*quadraticTerm)*aSquaredPlusBSquaredTerm;
    double RHS = 1.0 - (MOverMPrime*brackets);
    double e = sqrt(RHS);
    return e;
}

double aPostSupernovaInitiallyCircular(double a, double M, double MPrime, double theta, double phi, double uk){
    // Calculate aPrime for an initially circular orbit (e=0)
    double aPrime         = 0.0;
    double bottomFraction = 1.0;                                        // On bottom so better initialise to 1. CPLB: Dividing by zero might be a useful failure mode?
    double MOverMPrime    = M/MPrime;
    double quadraticTerm  = 1.0 + 2.0*uk*cos(theta)*cos(phi) + uk*uk;
    bottomFraction        = 2.0 - MOverMPrime*quadraticTerm;
    aPrime                = a/bottomFraction;
    return aPrime;
}

// Initially circular versions
double iPostSupernovaInitiallyCircular(double uk, double theta, double phi){
    // Assumes that the spins of the stars are initially aligned with the orbital plane
    // Post-Supernova orbital inclination
    double cosi = (1.0 + uk*cos(theta)*cos(phi))/sqrt((uk*sin(theta))*(uk*sin(theta)) + (1.0 + uk*cos(theta)*cos(phi))*(1.0 + uk*cos(theta)*cos(phi)));
    return acos(cosi);
}

double cosiPostSupernovaInitiallyCircular(double uk, double theta, double phi){
    // Return cos(i) post supernova orbital inclination so you can choose to output cos(i) instead
    return (1.0 + uk*cos(theta)*cos(phi))/sqrt((uk*sin(theta))*(uk*sin(theta)) + (1.0 + uk*cos(theta)*cos(phi))*(1.0 + uk*cos(theta)*cos(phi)));
}

// Post supernova orbital parameters for general eccentric orbits
double ePostSupernova(double uk, double theta, double phi, double beta, double M, double MPrime, double r, double a,  unsigned long randomSeed){
    /*
     This is Equation (31) in my post-SN orbital characteristics 2 notes - simplifies to Equation (2.8) in Brandt & Podsiadlowski 1995 (http://arxiv.org/abs/astro-ph/9412023) when e = 0, beta = pi/2
     This is also given by Equation A.12 in Hurley et al 2002 (http://arxiv.org/pdf/astro-ph/0201220v1.pdf)
     
     Parameters
     -----------
     uk : double
     Dimensionless kick velocity vk/vrel
     theta : double
     Angle
     phi : double
     Angle
     beta : double
     Angle the preSN velocity vector makes to the radial position vector
     M : double
     Total system mass before the supernova
     Mprime : double
     Total system mass after the supernova
     r : double
     Current orbital separation in AU (equal to a for eccentricity of 0) at the time of the supernova (sampled randomly for an eccentric orbit)
     a : double
     Semi major axis in AU
     
     Returns
     --------
     ePrime : double
     Orbital eccentricity after a supernova

     */
    
    double MOverMPrime              = M/MPrime;
    double twoOverrMinusOneOvera    = 2.0/r - 1.0/a;
    double quadraticTerm            = 1.0 + 2.0*uk*cos(theta)*cos(phi) + uk*uk;
    double firstSquareBrackets      = uk*uk*sin(theta)*sin(theta) + (uk*cos(theta)*sin(phi)*cos(beta) - sin(beta)*(uk*cos(theta)*cos(phi) + 1))*(uk*cos(theta)*sin(phi)*cos(beta) - sin(beta)*(uk*cos(theta)*cos(phi) + 1.0));
    double secondSquareBrackets     = 2.0/r - MOverMPrime*twoOverrMinusOneOvera*quadraticTerm;
    double oneMinuseSquared         = r*r*MOverMPrime*twoOverrMinusOneOvera*firstSquareBrackets*secondSquareBrackets;
    double eSquared                 = 1.0 - oneMinuseSquared;
    
    // Deal with small number rounding problems
    if(eSquared < 1E-8){
        eSquared = 0.0;
    }
    
    double final_e                  = sqrt(eSquared);
    
//    std::cout << "MOverMPrime = " << MOverMPrime << std::endl;
//    std::cout << "twoOverrMinusOneOvera = " << twoOverrMinusOneOvera << std::endl;
//    std::cout << "quadraticTerm = " << quadraticTerm << std::endl;
//    std::cout << "firstSquareBrackets = " << firstSquareBrackets << std::endl;
//    std::cout << "secondSquareBrackets = " << secondSquareBrackets << std::endl;
//    std::cout << "oneMinuseSquared = " << oneMinuseSquared << std::endl;
//    std::cout << "eSquared = " << eSquared << std::endl;
//    std::cout << "final_e = " << final_e << std::endl;
    
    return final_e;
    
}

double aPostSupernova(double uk, double M, double MPrime, double r, double a, double theta, double phi,  unsigned long randomSeed){
    /*
     This function calculates the post-supernova semi major axis
     
     This is Equation (22) in the post-SN orbital characteristics 2 document
     
     Parameters
     ----------
     uk : double
     Dimensionless kick velocity in units of the preSN orbital velocity
     M : double
     Total system mass before the supernova
     MPrime : double
     Total system mass after the supernova
     r : double
     r is the instantaneous separation at the moment of the SN (accounts for eccentricity) in AU
     a : double
     Semi major axis of the orbit before the supernova
     theta : double
     theta is the kick direction angle out of the plane
     phi : double
     phi is the kick direction angle in the plane
     
     Returns
     -------
     aPrime : double
     Semi major axis of the orbit after the supernova

     */
    
    double MOverMPrime   = M/MPrime;
    double firstBrackets = (2.0/r - 1.0/a);
    double quadraticTerm = 1.0 + 2.0*uk*cos(theta)*cos(phi) + uk*uk;
    double bottom        = 2.0/r - MOverMPrime*firstBrackets*quadraticTerm;
    double final_a       = 1.0/bottom;
    
    return final_a;
    
}


// Versions for arbitrary initial eccentricity
double cosFinalPlaneTilt(double uk, double beta, double theta, double phi){
    // Calculate cos(i), the cos of the tilt between the pre and post SN orbital planes (as defined by the angular momentum)
    // uk is dimensionless kick velocity
    // beta is angle between velocity and radius vectors at moment of SN which encodes information about the location on the orbit (anomaly) and the eccentricity of the orbit
    // theta and phi are angles of kick velocity
    // Eq (40) in post-SN orbital characteristics 2
    
    // I think the top is okay
    double topFraction      = (sin(beta)  *(1.0 + uk  *cos(theta)  *cos(phi))) - (uk  *cos(theta)  *sin(phi)  *cos(beta));
    
    // Bottom may be wrong.
    //double bottomFraction   = sqrt( uk*uk*sin(beta)*sin(beta) + uk*uk*cos(theta)*cos(theta)*sin(phi)*sin(phi)*cos(beta)*cos(beta) + sin(beta)*sin(beta)*(1.0 + uk  *cos(theta)  *cos(phi)));
    
    double bottomFraction = sqrt(uk*uk*sin(theta)*sin(theta) + uk*uk*cos(theta)*cos(theta)*sin(phi)*sin(phi)*cos(beta)*cos(beta) + sin(beta)*sin(beta)*pow((uk*cos(theta)*cos(phi) + 1.0), 2.0) - 2.0*cos(theta)*sin(phi)*cos(beta)*sin(beta)*(uk*cos(theta)*cos(phi) + 1.0));
    
    // DEBUGGING
    //std::cout << "top " << topFraction << std::endl;
    //std::cout << "bottom " << bottomFraction << std::endl;
    
    return topFraction / bottomFraction;
}

double finalPlaneTilt(double uk, double beta, double theta, double phi){
    // arccos(above)
    double cosi = cosFinalPlaneTilt(uk, beta, theta, phi);
    return acos(cosi);
}

double finalPlaneTiltVector(double mu, double muPrime, double radius, double vrel, double vk, double theta, double phi, double beta){
    // Try using the linear algebra stuff to calculate the dot product for the tilt. Is it easier? Requires many more parameters (since need to define before and after)
    
    // Initialise vectors
    double LPrimex      = 0;
    double LPrimey      = 0;
    double LPrimez      = 0;
    double normLPrime   = 0;
    
    double Lx           = 0;
    double Ly           = 0;
    double Lz           = 0;
    double normL        = 0;

    // Calculate L (only has a z component)
    Lz = - mu  *radius  *vrel  *sin(beta);
    
    // Calculate L'
    LPrimex = muPrime  *radius  *vk  *sin(theta)  *sin(beta);
    LPrimey = - muPrime  *radius  *vk  *sin(theta)  *cos(beta);
    LPrimez = muPrime  *( radius  *vk  *cos(theta)  *sin(phi)  *cos(beta) - radius  *sin(beta)  *(vk  *cos(theta)  *cos(phi) + vrel));
    
    // Calculate L.L'/norm(L) norm(L')
    normL       = norm3(Lx, Ly, Lz);
    normLPrime  = norm3(LPrimex, LPrimey, LPrimez);
    
    double cosi = 0;
    
    cosi = dotProduct3(Lx, Ly, Lz, LPrimex, LPrimey, LPrimez) / (normL  *normLPrime);
    
    return acos(cosi);
}

// Some function here for fixing uK to a certain value -- where should it live?

void fixUK(const programOptions &options, double & uk){
    // determine if user wishes to set uK to a constant value, and change uk accordingly.
    if(options.useFixedUK){
        // Set uk to a predetermined value
        uk = options.fixedUK;
    }
    else{
        // Do nothing
        uk = uk;
    }
}


double postSNSystemicVelocity(double uk, double vorb, double Mtot, double MtotPrime, double M1Prime, double deltaM1, double M2, double theta, double phi, double beta){
    /*
    Calculate the systemic velocity (center-of-mass velocity) of the binary after the supernova

    See Equation 2.10 in Brandt & Podsiadlowski 1995 https://arxiv.org/pdf/astro-ph/9412023.pdf 
    or Equation A14 in Hurley et al 2002 https://arxiv.org/pdf/astro-ph/0201220.pdf

    Parameters
    -----------
    uk : double
        Dimensionless kick velocity vkick/vorb

    Returns
    --------
    vsys : double
        Post supernova systemic velocity
    
    */

    // Calculate some commonly used variables
    double uksquared = uk * uk;
    double M1Primesquared = M1Prime * M1Prime;
    
    double sintheta = sin(theta);
    double costheta = cos(theta);

    double cosphi = cos(phi);

    double sinbeta = sin(beta);
    double cosbeta = cos(beta);

    /////////////////////

    double prefactor = vorb / MtotPrime;

    double firstTerm = pow((deltaM1 * M2 / Mtot), 2);
    double secondTerm = M1Primesquared * uksquared;
    double thirdTerm = 2.0 * (M1Prime * deltaM1 * M2 / Mtot) * uk * sintheta * cosphi * cosbeta;
    double fourthTerm = 2.0 * (M1Prime * deltaM1 * M2 / Mtot) * uk * costheta * cosphi * sinbeta;

    return prefactor * sqrt(firstTerm + secondTerm + thirdTerm + fourthTerm);
}

