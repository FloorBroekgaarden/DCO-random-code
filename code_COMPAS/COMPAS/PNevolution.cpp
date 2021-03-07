//
//  PNevolution.cpp


#include "PNevolution.h"

// Implementation

double orbitalVelocity(const double &M, const double a){
    // This is Equation (38) in my precessionNotes document
    // Calculate orbital velocity for a circular orbit in G = c = 1 units
    // M is the total mass (in m)
    // a is the semi-major axis (in m)
    return sqrt(M / a); // CPLB: Average? SPS: updated above.
}

double orbitalFrequency(const double &M, const double a){
    // This is Equation (25) in my precessionNotes document
    // Calculate the orbital frequency Omega for a circular orbit in G = c = 1 units
    // M is the total mass (in m)
    // a is the semi-major axis (in m)
    double brackets = ( (M) / (a*a*a) ); // CPLB: Circular orbit? SPS: yup
    return sqrt(brackets);
}

double orbitalFrequencyFuncNu(const double &M, const double &nu){
    // This is Equation (39) in the precessionNotes document
    // Calculate the orbital frequency Omega for a circular orbit in G = c = 1 units
    // M is the total mass (in m)
    // nu is the orbital velocity (this is useful since this is the parameter we evolve in the differential equations)
    // CPLB: What are inputs? SPS: updated above.
    return (nu*nu*nu)/(M);
}

double gravitationalWaveFreqFuncNu(const double &M, const double &nu){
    // Calculate the gravitational wave frequency for a circular orbit as a function of orbital velocity nu
    // M is the total mass (in m)
    // Gravitational wave frequency = 2 * orbital frequency (dominant mode) CPLB: Approximately. Exact for circular.
    return 2.0 * orbitalFrequencyFuncNu(M, nu);
}

void dadt(const double &m1, const double &m2, const double &e, const double &a, double &dadt){
    // This is in Peters (1964 - http://journals.aps.org/pr/pdf/10.1103/PhysRev.136.B1224) Equation 5.6
    // Calculate the rate of change of the semi-major axis = dadt in G = c = 1 units
    // m1, m2 are the component masses (in m)
    // e is the eccentricity
    // a is the semi-major axis
    double M = m1 + m2;
    dadt = -(64.0/5.0)*(m1*m2*M/(a*a*a*pow((1.0 - e*e),3.5)))*(1.0 + (73.0/24.0)*e*e + (37.0/96.0)*e*e*e*e);
}

void dedt(const double &m1, const double &m2, const double &e, const double &a, double &dedt){
    // This is in Peters (1964 - http://journals.aps.org/pr/pdf/10.1103/PhysRev.136.B1224) Equation 5.7
    // Calculate the rate of change of eccentricity = dedt in G = c = 1 units
    // m1, m2 are the component masses (in m)
    // e is the eccentricity
    // a is the semi-major axis
    double M = m1 + m2;
    dedt = -(304.0/15.0)*e*(m1*m2*M/(a*a*a*a*pow((1.0 - e*e), 2.5)))*(1.0 + (121.0/304.0)*e*e);
}

void dLdt(const double &m1, const double &m2, const double &e, const double &a, double &dLdt){
    // This is in Peters (1964 - http://journals.aps.org/pr/pdf/10.1103/PhysRev.136.B1224) Equation 5.5
    // Peters equation for the evolution of the magnitude of L = dLdt due to radiation reaction
    // m1, m2 are the component masses (in m)
    // e is the eccentricity
    // a is the semi-major axis
    double M = m1 + m2;
    dLdt = -(32.0/5.0)*m1*m1*m2*m2*(sqrt(M)/(pow(a, 3.5) * (1.0 - e*e) * (1.0 - e*e)) ) * (1.0 + (7.0/8.0)*e*e);
}

double calculateDeltaPhi(const double &S1x, const double &S1y, const double &S1z, const double &S2x, const double &S2y, const double &S2z, const double &Lhatx, const double &Lhaty, const double &Lhatz){
    // This is Equation (7) in the precessionNotes
    // Calculate the angle between the projections of S1 and S2 on the orbital plane defined by L which is called deltaPhi
    // What if s1 // s2 // L i.e. aligned spins? Beware of divide by zero giving NAN. arbitrarily set to 0 to stop other functions breaking due to NAN.
    
    // Create unit spin vectors
    double uS1x = 0.0;
    double uS1y = 0.0;
    double uS1z = 0.0;
    
    double uS2x = 0.0;
    double uS2y = 0.0;
    double uS2z = 0.0;
    
    unit3(S1x, S1y, S1z, uS1x, uS1y, uS1z);
    unit3(S2x, S2y, S2z, uS2x, uS2y, uS2z);
    
    // Create temporary vectors
    
    double first_term_x         = 0.0;                                      // tmp for calculating dphi
    double first_term_y         = 0.0;                                      // tmp for calculating dphi
    double first_term_z         = 0.0;                                      // tmp for calculating dphi
    double first_term_denom     = 0.0;
    
    crossProduct3(uS1x, uS1y, uS1z, Lhatx, Lhaty, Lhatz, first_term_x, first_term_y, first_term_z);
    first_term_denom = norm3(first_term_x, first_term_y, first_term_z);
    
    double second_term_x        = 0.0;                                      // tmp for calculating dphi
    double second_term_y        = 0.0;                                      // tmp for calculating dphi
    double second_term_z        = 0.0;                                      // tmp for calculating dphi
    double second_term_denom    = 0.0;
    
    crossProduct3(uS2x, uS2y, uS2z, Lhatx, Lhaty, Lhatz, second_term_x, second_term_y, second_term_z);
    second_term_denom = norm3(second_term_x, second_term_y, second_term_z);
    
    double cosdphi = 0.0;
    
    if(first_term_denom != 0 and second_term_denom != 0){
        first_term_x /= first_term_denom;
        first_term_y /= first_term_denom;
        first_term_z /= first_term_denom;
    
        second_term_x /= second_term_denom;
        second_term_y /= second_term_denom;
        second_term_z /= second_term_denom;
    }
    
    cosdphi = dotProduct3(first_term_x, first_term_y, first_term_z, second_term_x, second_term_y, second_term_z);
    
    double dphi = 0.0;
    dphi = acos(cosdphi);
    
    return dphi;
}

double calculateTheta12(const double &S1x, const double &S1y, const double &S1z, const double &S2x, const double &S2y, const double &S2z){
    // This is Equation (8) in precessionNotes
    // Calculate the angle between the two spin vectors S1 and S2 called theta12 in the notes and in Gerosa et al
    
    // Calculate the unit vectors
    double uS1x = 0.0;
    double uS1y = 0.0;
    double uS1z = 0.0;
    
    double uS2x = 0.0;
    double uS2y = 0.0;
    double uS2z = 0.0;
    
    unit3(S1x, S1y, S1z, uS1x, uS1y, uS1z);
    unit3(S2x, S2y, S2z, uS2x, uS2y, uS2z);
    
    double costheta12   = 0.0;
    
    costheta12 = dotProduct3(uS1x, uS1y, uS1z, uS2x, uS2y, uS2z);
    
    return acos(costheta12);
}

double calculateTheta(const double &Sx, const double &Sy, const double &Sz, const double &Lhatx, const double &Lhaty, const double &Lhatz){
    // This is Equations (5) and (6) in the precession Notes
    // Calculate the angle between a spin vector S and the orbital angular momentum vector L
    
    // Calcluate the spin unit vector
    double uSx = 0.0;
    double uSy = 0.0;
    double uSz = 0.0;
    
    unit3(Sx, Sy, Sz, uSx, uSy, uSz);
    
    // Calcluate cos(theta) = uS dot L
    double costheta = dotProduct3(uSx, uSy, uSz, Lhatx, Lhaty, Lhatz);
    
    return acos(costheta);
}

double calculatePhi(const double &Sx, const double &Sy, const double &Sz, const double &Lhatx, const double &Lhaty, const double &Lhatz){
    // This is Equation (61) in the precession Notes
    // Calculate angle between projections of S and x on orbital plane
    
    // Create unit spin vector
    double uSx = 0.0;
    double uSy = 0.0;
    double uSz = 0.0;
    
    unit3(Sx, Sy, Sz, uSx, uSy, uSz);
    
    // Define x axis
    double uxx = 1.0;
    double uxy = 0.0;
    double uxz = 0.0;
    
    // Create temporary vectors
    
    double first_term_x         = 0.0;                                      // tmp for calculating phi
    double first_term_y         = 0.0;                                      // tmp for calculating phi
    double first_term_z         = 0.0;                                      // tmp for calculating phi
    double first_term_denom     = 0.0;
    
    crossProduct3(uSx, uSy, uSz, Lhatx, Lhaty, Lhatz, first_term_x, first_term_y, first_term_z);
    first_term_denom = norm3(first_term_x, first_term_y, first_term_z);
    
    double second_term_x        = 0.0;                                      // tmp for calculating phi
    double second_term_y        = 0.0;                                      // tmp for calculating phi
    double second_term_z        = 0.0;                                      // tmp for calculating phi
    double second_term_denom    = 0.0;
    
    crossProduct3(uxx, uxy, uxz, Lhatx, Lhaty, Lhatz, second_term_x, second_term_y, second_term_z);
    second_term_denom = norm3(second_term_x, second_term_y, second_term_z);
    
    double cosphi = 0.0;
    
    if(first_term_denom != 0 and second_term_denom != 0){
        first_term_x /= first_term_denom;
        first_term_y /= first_term_denom;
        first_term_z /= first_term_denom;
    
        second_term_x /= second_term_denom;
        second_term_y /= second_term_denom;
        second_term_z /= second_term_denom;
    }
    
    cosphi = dotProduct3(first_term_x, first_term_y, first_term_z, second_term_x, second_term_y, second_term_z);
    
    double phi = 0.0;
    phi = acos(cosphi);
    
    return phi;
    
}

void calculateS1(double &S1x, double &S1y, double &S1z, double chi1, double m1, double theta1, double deltaPhi){
    // Calculate S1
    // CPLB: chi1 * m1 * m1 * [output from calculate US1]
    
    S1x = chi1 * m1 * m1 * sin(theta1) * sin(deltaPhi);
    S1y = chi1 * m1 * m1 * sin(theta1) * cos(deltaPhi);
    S1z = chi1 * m1 * m1 * cos(theta1);
    
}

void calculateS2(double &S2x, double &S2y, double &S2z, double chi2, double m2, double theta2){
    // Calculate S2
    // CPLB: calculateS1(S2x, S2y, S2z, chi2, m2, theta2, 0.0)
    
    S2x = 0.0;
    S2y = chi2 * m2 * m2 * sin(theta2);
    S2z = chi2 * m2 * m2 * cos(theta2);
    
}

void calculateUS1(double &S1x, double &S1y, double &S1z, double theta1, double deltaPhi){
    // Calculate unit vector for S1
    
    S1x = sin(theta1) * sin(deltaPhi);
    S1y = sin(theta1) * cos(deltaPhi);
    S1z = cos(theta1);
    
}

void calculateUS2(double &S2x, double &S2y, double &S2z, double theta2){
    // Calculate unit vector for S2
    // CPLB: calculateUS1(S2x, S2y, S2z, theta2, 0.0)
    
    S2x = 0.0;
    S2y = sin(theta2);
    S2z = cos(theta2);
    
}

void Omega1(double &Omega1x, double &Omega1y, double &Omega1z, const double &m1, const double &m2, const double &uLx, const double &uLy, const double &uLz, const double &s1x, const double &s1y, const double &s1z, const double &s2x, const double &s2y, const double &s2z, const double &nu){
    
    // This is Equation (31) in precessionNotes and Equation (14) in Gerosa et al 2013 (http://arxiv.org/abs/1302.4442 )
    // Calculates the precession vector Omega1 for spin vector S1 in G = c = 1 units
    // m1, m2 are component masses (in m)
    // nu is the orbital velocity
    
    // Calculate some combinations of the components masses
    double q    = m2/m1;                    // Mass ratio
    double dm   = m1 - m2;                  // Difference in mass, delta m
    double Mtot = m1 + m2;                  // Total mass M
    double eta  = (m1*m2)/(Mtot*Mtot);      // Symmetric mass ratio eta
    
    // Precompute dot products
    double Ldots1 = dotProduct3(uLx, uLy, uLz, s1x, s1y, s1z);
    double Ldots2 = dotProduct3(uLx, uLy, uLz, s2x, s2y, s2z);
    
    // Calculate the components of Omega
    //CPLB: Gerose et al. only quote up to O(nu^6), where doees O(nu^7) term come from?
    Omega1x = eta*nu*nu*nu*nu*nu*(2.0 + (3.0*q)/(2.0))*uLx + ((nu*nu*nu*nu*nu*nu)/(2.0*Mtot*Mtot))*(s2x - 3.0*Ldots2*uLx - 3.0*q*Ldots1*uLx) + nu*nu*nu*nu*nu*nu*nu*( (9.0)/(16.0) + (5.0*eta)/(4.0) - (eta*eta)/(24.0) + ((dm)/(Mtot))*( (-9.0)/(16.0) + (5.0*eta)/(8.0) ) )*uLx;
    Omega1y = eta*nu*nu*nu*nu*nu*(2.0 + (3.0*q)/(2.0))*uLy + ((nu*nu*nu*nu*nu*nu)/(2.0*Mtot*Mtot))*(s2y - 3.0*Ldots2*uLy - 3.0*q*Ldots1*uLy) + nu*nu*nu*nu*nu*nu*nu*( (9.0)/(16.0) + (5.0*eta)/(4.0) - (eta*eta)/(24.0) + ((dm)/(Mtot))*( (-9.0)/(16.0) + (5.0*eta)/(8.0) ) )*uLy;
    Omega1z = eta*nu*nu*nu*nu*nu*(2.0 + (3.0*q)/(2.0))*uLz + ((nu*nu*nu*nu*nu*nu)/(2.0*Mtot*Mtot))*(s2z - 3.0*Ldots2*uLz - 3.0*q*Ldots1*uLz) + nu*nu*nu*nu*nu*nu*nu*( (9.0)/(16.0) + (5.0*eta)/(4.0) - (eta*eta)/(24.0) + ((dm)/(Mtot))*( (-9.0)/(16.0) + (5.0*eta)/(8.0) ) )*uLz;
    
    // Divide by Mtot
    Omega1x /= Mtot;
    Omega1y /= Mtot;
    Omega1z /= Mtot;
}

// Same (but different) for Omega2
void Omega2(double &Omega2x, double &Omega2y, double &Omega2z, const double &m1, const double &m2, const double &uLx, const double &uLy, const double &uLz, const double &s1x, const double &s1y, const double &s1z, const double &s2x, const double &s2y, const double &s2z, const double &nu){
    // This is Equation (32) in precessionNotes and Equation (15) in Gerosa et al 2013 (http://arxiv.org/abs/1302.4442 )
    // Calculates the precession vector Omega2 for spin vector S2 in G = c = 1 units
    // m1, m2 are component masses (in m)
    // nu is the orbital velocity
    
    // Calculate some combinations of the components masses
    double q    = m2/m1;                    // Mass ratio
    double dm   = m1 - m2;                  // Difference in mass, delta m
    double Mtot = m1 + m2;                  // Total mass M
    double eta  = (m1*m2)/(Mtot*Mtot);      // Symmetric mass ratio eta
    
    // Precompute dot products
    double Ldots1 = dotProduct3(uLx, uLy, uLz, s1x, s1y, s1z);
    double Ldots2 = dotProduct3(uLx, uLy, uLz, s2x, s2y, s2z);
    
    // Remember to switch all of the things around
    Omega2x = eta*nu*nu*nu*nu*nu*(2.0 + (3.0)/(2.0*q))*uLx + ((nu*nu*nu*nu*nu*nu)/(2.0*Mtot*Mtot))*(s1x - 3.0*Ldots1*uLx - ((3.0)/(q))*Ldots2*uLx) + nu*nu*nu*nu*nu*nu*nu*( (9.0)/(16.0) + (5.0*eta)/(4.0) - (eta*eta)/(24.0) - ((dm)/(Mtot))*( (-9.0)/(16.0) + (5.0*eta)/(8.0) ) )*uLx;
    Omega2y = eta*nu*nu*nu*nu*nu*(2.0 + (3.0)/(2.0*q))*uLy + ((nu*nu*nu*nu*nu*nu)/(2.0*Mtot*Mtot))*(s1y - 3.0*Ldots1*uLy - ((3.0)/(q))*Ldots2*uLy) + nu*nu*nu*nu*nu*nu*nu*( (9.0)/(16.0) + (5.0*eta)/(4.0) - (eta*eta)/(24.0) - ((dm)/(Mtot))*( (-9.0)/(16.0) + (5.0*eta)/(8.0) ) )*uLy;
    Omega2z = eta*nu*nu*nu*nu*nu*(2.0 + (3.0)/(2.0*q))*uLz + ((nu*nu*nu*nu*nu*nu)/(2.0*Mtot*Mtot))*(s1z - 3.0*Ldots1*uLz - ((3.0)/(q))*Ldots2*uLz) + nu*nu*nu*nu*nu*nu*nu*( (9.0)/(16.0) + (5.0*eta)/(4.0) - (eta*eta)/(24.0) - ((dm)/(Mtot))*( (-9.0)/(16.0) + (5.0*eta)/(8.0) ) )*uLz;
    
    // Divide by Mtot
    Omega2x /= Mtot;
    Omega2y /= Mtot;
    Omega2z /= Mtot;
}

void dS1dt(double &ds1xdt, double &ds1ydt, double &ds1zdt, const double &s1x, const double &s1y, const double &s1z, const double &s2x, const double &s2y, const double &s2z, const double &m1, const double &m2, const double &lhatx, const double&lhaty, const double &lhatz, const double &nu){
    // Calculate the time derivative of S1 - dS1/dt = Omega1 X S1
    // This is Equation (29) in precessionNotes and Equation (14) in Gerosa et al 2013 (http://arxiv.org/abs/1302.4442)
    // m1, m2 are component masses
    // nu is the orbital velocity
    
    // Declare precession vector Omega1
    double Omega1x = 0.0;
    double Omega1y = 0.0;
    double Omega1z = 0.0;
    
    // Calculate precession vector Omega1
    Omega1(Omega1x, Omega1y, Omega1z, m1, m2, lhatx, lhaty, lhatz, s1x, s1y, s1z, s2x, s2y, s2z, nu);
    
    // Calculate dS1/dt = Omega1 X S1
    crossProduct3(Omega1x, Omega1y, Omega1z, s1x, s1y, s1z, ds1xdt, ds1ydt, ds1zdt);
    
}

void dS2dt(double &ds2xdt, double &ds2ydt, double &ds2zdt, const double &s1x, const double &s1y, const double &s1z, const double &s2x, const double &s2y, const double &s2z, const double &m1, const double &m2, const double &lhatx, const double&lhaty, const double &lhatz, const double &nu){
    // Calculate the time derivative of S2 - dS2/dt = Omega2 X S2
    // This is Equation (30) in precessionNotes and Equation (15) in Gerosa et al 2013 (http://arxiv.org/abs/1302.4442)
    // m1, m2 are component masses
    // nu is the orbital velocity
    
    // Declare precession vector Omega2
    double Omega2x = 0.0;
    double Omega2y = 0.0;
    double Omega2z = 0.0;
    
    // Calculate precession vector Omega2
    Omega2(Omega2x, Omega2y, Omega2z, m1, m2, lhatx, lhaty, lhatz, s1x, s1y, s1z, s2x, s2y, s2z, nu);
    
    // Calculate dS2/dt = Omega2 X S2
    crossProduct3(Omega2x, Omega2y, Omega2z, s2x, s2y, s2z, ds2xdt, ds2ydt, ds2zdt);
    
}

void dnudt(const double &m1, const double &m2, double &dnudt, const double &S1x, const double &S1y, const double &S1z, const double &S2x, const double &S2y, const double &S2z, const double &uLx, const double &uLy, const double &uLz, const double &nu){
    // uS1x = x component of S1 unit vector
    // m1, m2 are component masses (in m)
    // nu = orbital velocity
    // S1, S2 = spin vectors
    // Everything is in G = c = 1 units
    // Note that L is already a unit vector here as indicated by the u
    // This is Equation (37) in precessionNotes and Equation (17) in Gerosa et al 2013 (http://arxiv.org/abs/1302.4442)
    
    // Calculate some combinations of the components masses
    double Mtot = m1 + m2;
    double eta  = (m1*m2)/(Mtot*Mtot);
    
    double prefactor = (32.0/5.0)*(eta/Mtot)*pow(nu,9.0);
    
    // Calculate spin unit vectors
    double uS1x = 0.0;
    double uS1y = 0.0;
    double uS1z = 0.0;
    double chi1 = 0.0;
    
    double uS2x = 0.0;
    double uS2y = 0.0;
    double uS2z = 0.0;
    double chi2 = 0.0;
    
    uS1x = S1x/norm3(S1x, S1y, S1z);
    uS1y = S1y/norm3(S1x, S1y, S1z);
    uS1z = S1z/norm3(S1x, S1y, S1z);
    chi1 = norm3(S1x, S1y, S1z)/(m1*m1);    // Dimensionless spin magnitude 1
    
    uS2x = S2x/norm3(S2x, S2y, S2z);
    uS2y = S2y/norm3(S2x, S2y, S2z);
    uS2z = S2z/norm3(S2x, S2y, S2z);
    chi2 = norm3(S2x, S2y, S2z)/(m2*m2);    // Dimensionless spin magnitude 2
    
    // Precompute vector calculations? -- only has dot products so basically nice
    double s1dotL   = dotProduct3(uS1x, uS1y, uS1z, uLx, uLy, uLz);
    double s2dotL   = dotProduct3(uS2x, uS2y, uS2z, uLx, uLy, uLz);
    double s1dots2  = dotProduct3(uS1x, uS1y, uS1z, uS2x, uS2y, uS2z);
    
    // Have each order of nu as a separate term. The term 'second' etc refers to power of nu, so second is term with nu^2
    double zeroth   = 1.0;

    double first    = 0.0;

    double second   = - nu*nu*((743.0 + 924.0*eta)/(336.0)); // CPLB: Corrected sign (double check!!) // SPS: Cheers!

    // This is the first awkward term.
    double third    = nu*nu*nu*( 4.0 * M_PI - (chi1 * s1dotL * ( (113.0*m1*m1)/(12.0*Mtot*Mtot) + (25.0*eta)/(4.0) ) ) - (chi2 * s2dotL * ( (113.0*m2*m2)/(12.0*Mtot*Mtot) + (25.0*eta)/(4.0))) );

    // This is the other awkward term. (of course awkward here means most important)
    double fourth   = nu*nu*nu*nu*( (34103.0)/(18144.0) + (13661.0*eta)/(2016.0) + (59.0*eta*eta)/(18.0) + ( (eta*chi1*chi2)/(48.0) )*(721.0*s1dotL*s2dotL - 247.0*s1dots2) + ((1.0)/(96.0))*( (pow((m1*chi1/Mtot),2.0)*(719.0*s1dotL*s1dotL - 233.0)) + (pow((m2*chi2/Mtot),2.0)*(719.0*s2dotL*s2dotL - 233.0)) ) );

    double fifth    = - nu*nu*nu*nu*nu*M_PI*( (4159.0 + 15876.0*eta)/(672.0) );

    double sixth    = nu*nu*nu*nu*nu*nu*( (16447322263.0)/(13970880.0) + (16.0*M_PI*M_PI)/(3.0) - (1712.0 * (gammaE + log(4.0*nu)))/(105.0) + ( (451.0*M_PI*M_PI)/(48.0) - (56198689.0)/(217728.0) )*eta + (541.0*eta*eta)/(896.0) - (5605.0*eta*eta*eta)/(2592.0) );

    double seventh  = nu*nu*nu*nu*nu*nu*nu*M_PI*( (-4415.0)/(4032.0) + (358675.0*eta)/(6048.0) + (91495.0*eta*eta)/(1512.0) );
    
    // Horrible eta equation goes here (have split up to emphasise PN expansion)
    dnudt = prefactor*(zeroth + first + second + third + fourth + fifth + sixth + seventh);
}

void dLhatdt(double &dLhatxdt, double &dLhatydt, double &dLhatzdt, const double &m1, const double &m2, const double &S1x, const double &S1y, const double &S1z, const double &S2x, const double &S2y, const double &S2z, const double &Lhatx, const double &Lhaty, const double &Lhatz, const double &nu, const double &ds1xdt, const double &ds1ydt, const double &ds1zdt, const double &ds2xdt, const double &ds2ydt, const double &ds2zdt){
    // Calculate the time derivative of the orbital angular momentum unit vector Lhat defined as dLhat_dt
    // This is Equation (23) in precessionNotes and Equation (16) in Gerosa et al 2013 (http://arxiv.org/abs/1302.4442)
    // m1, m2 are the component masses in G = c = 1 units (i.e. in m)
    // S1, S2 are the spin vectors
    // nu is the orbital velocity
    // Also depends on derivatives of individual spins ds1_dt and ds2_dt
    
    // Evolution of the direction of angular momentum vector L
    double Mtot = m1 + m2;
    double eta  = (m1*m2)/(Mtot*Mtot);
    
    // Precompute some cumbersome factors
    double prefactor            = (-nu)/(eta*Mtot*Mtot);
    double denominator          = 1.0 + nu*nu*((3.0/2.0) + (eta/6.0));
    double fractionPrefactor    = (2.0*eta*nu*nu*nu*nu*nu*nu*nu)/(Mtot);
    
    // Need to calculate a temporary vector quantity P1 = Lhat X S1
    double P1x = 0.0;
    double P1y = 0.0;
    double P1z = 0.0;
    
    crossProduct3(Lhatx, Lhaty, Lhatz, S1x, S1y, S1z, P1x, P1y, P1z);
    
    // Need to calculate a temporary vector quantity P2 = Lhat X S2
    double P2x = 0.0;
    double P2y = 0.0;
    double P2z = 0.0;
    
    crossProduct3(Lhatx, Lhaty, Lhatz, S2x, S2y, S2z, P2x, P2y, P2z);
    
    // Now can calculate dLhat/dt
    dLhatxdt = prefactor*( ds1xdt + ds2xdt - fractionPrefactor*(eta + (m2/Mtot)*(m2/Mtot)*(1.0 + (3.0*m2)/(16.0*m1)))*P1x - fractionPrefactor*(eta + (m1/Mtot)*(m1/Mtot)*(1.0 + (3.0*m1)/(16.0*m2)))*P2x)/denominator;
    dLhatydt = prefactor*( ds1ydt + ds2ydt - fractionPrefactor*(eta + (m2/Mtot)*(m2/Mtot)*(1.0 + (3.0*m2)/(16.0*m1)))*P1y - fractionPrefactor*(eta + (m1/Mtot)*(m1/Mtot)*(1.0 + (3.0*m1)/(16.0*m2)))*P2y)/denominator;
    dLhatzdt = prefactor*( ds1zdt + ds2zdt - fractionPrefactor*(eta + (m2/Mtot)*(m2/Mtot)*(1.0 + (3.0*m2)/(16.0*m1)))*P1z - fractionPrefactor*(eta + (m1/Mtot)*(m1/Mtot)*(1.0 + (3.0*m1)/(16.0*m2)))*P2z)/denominator;
}

double angMomentumNat2(double m1, double m2, double a, double e){
    // Takes total mass M = m1 + m2 and semi-major axis a and eccentricity e in natural units
    // Doesn't assume circular orbit
    // Returns the angular momentum in natural units (m^2)
    // This is Equation (24) in precessionNotes
    double M  = m1 + m2;
    double mu = m1*m2/(M);
    return mu * sqrt(M * a * (1.0 - e*e));
}

double angMomentumFuncNu(double m1, double m2, double nu){
    // Takes component masses m1 and m2 in natural units and calcluates the angular momentum
    // This is Equation (24) in precessionNotes
    return (m1 * m2) / (nu);
}

void calculateRotationMatrix(double iota, double (&rotationMatrix)[3][3]){
    
    // Calculate the rotation matrix for a rotation of iota counter clockwise around the y-axis
    // CPLB: assumes other elements zero? Might be worth explicitly ensuring this!
    
    rotationMatrix[0][0] = cos(iota);
    rotationMatrix[0][2] = sin(iota);
    rotationMatrix[1][1] = 1.0;
    rotationMatrix[2][0] = -sin(iota);
    rotationMatrix[2][2] = cos(iota);
    
}

void convertToLAL(double &S1x, double &S1y, double &S1z, double &S2x, double &S2y, double &S2z, double &Lhatx, double &Lhaty, double &Lhatz, double iota){
    // This is Equations (68) and (69) in precessionNotes
    // Note that S1, S2 aren't unit vectors
    
//    // DEBUG
//    cout << "S1 = " << endl;
//    cout << S1x << endl;
//    cout << S1y << endl;
//    cout << S1z << endl;
//    
//    cout << "S2 = " << endl;
//    cout << S2x << endl;
//    cout << S2y << endl;
//    cout << S2z << endl;
    
    // Calculate and store theta1, theta2, deltaPhi
    double theta1 = calculateTheta(S1x, S1y, S1z, Lhatx, Lhaty, Lhatz);
    double theta2 = calculateTheta(S2x, S2y, S2z, Lhatx, Lhaty, Lhatz);
    double deltaPhi = calculateDeltaPhi(S1x, S1y, S1z, S2x, S2y, S2z, Lhatx, Lhaty, Lhatz);
    
    double normS1 = norm3(S1x, S1y, S1z);
    double normS2 = norm3(S2x, S2y, S2z);
    
    // Firstly move to frame with Lhat = Zhat
    Lhatx = 0.0;
    Lhaty = 0.0;
    Lhatz = 1.0;
    
    // Then recalculate S1, S2 in this frame - norms should remain same so only need to worry about unit vectors (remember to multiply norm back in later)
    calculateUS1(S1x, S1y, S1z, theta1, deltaPhi);
    calculateUS2(S2x, S2y, S2z, theta2);
    
    // Calculate rotation matrix and apply it
    double rotationMatrix[3][3];
    initialiseArrayToZero(rotationMatrix);
    
    // Calculate the rotation matrix
    calculateRotationMatrix(iota, rotationMatrix);
    
    // Create vectors for S1, S2, L
    double S1vect[3][1];
    double S2vect[3][1];
    double Lvect[3][1];
    
    // Package S1, S2, L as arrays
    packageVectorAsArray(Lhatx, Lhaty, Lhatz, Lvect);
    packageVectorAsArray(S1x, S1y, S1z, S1vect);
    packageVectorAsArray(S2x, S2y, S2z, S2vect);
    
    // Will need output vectors for these
    double S1vectPrime[3][1];
    double S2vectPrime[3][1];
    double LvectPrime[3][1];
    
    // Apply rotation matrix to L, S1, S2 to rotate to new coordinate system
    matrixVectorMultiplication(rotationMatrix, S1vect, S1vectPrime);
    matrixVectorMultiplication(rotationMatrix, S2vect, S2vectPrime);
    matrixVectorMultiplication(rotationMatrix, Lvect, LvectPrime);
    
    // Unpack vectors so can use algebra on them (could make into a unpackVector function)
    S1x = S1vectPrime[0][0];
    S1y = S1vectPrime[1][0];
    S1z = S1vectPrime[2][0];
    
    S2x = S2vectPrime[0][0];
    S2y = S2vectPrime[1][0];
    S2z = S2vectPrime[2][0];
    
    Lhatx = LvectPrime[0][0];
    Lhaty = LvectPrime[1][0];
    Lhatz = LvectPrime[2][0];
    
    // Multiply the magnitude back in
    S1x *= normS1;
    S1y *= normS1;
    S1z *= normS1;
    
    S2x *= normS2;
    S2y *= normS2;
    S2z *= normS2;
    
    // DEBUG
    //    cout << "S1 = " << endl;
    //    cout << S1x << endl;
    //    cout << S1y << endl;
    //    cout << S1z << endl;
    //
    //    cout << "S2 = " << endl;
    //    cout << S2x << endl;
    //    cout << S2y << endl;
    //    cout << S2z << endl;
    
    
}





///////////////////////////////////////////////////////////////////
//                          NOTES
///////////////////////////////////////////////////////////////////
