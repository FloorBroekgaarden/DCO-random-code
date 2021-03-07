//
//  tides.cpp


#include "tides.h"

// Orbital angular momentum conservation solution
//void tidesConvservationAngularMomentum(double stuff, double aPrime){
//    // Calculate the final semi major axis a' given an initial system
//    //
//    
//    
//    //
//    // Solve an nth order polynomial of the form
//    //
//    // P(x) = a_0 + a_1 x + a_2 x^2 + ... + a_{n-1} x^{n-1}
//    
//    // Coefficients of the polynomial
//    //    // Coefficients of the cubic : Example for (x + 1)(2x + 1)(x - 7)
//    
//    // Number of coefficients (is order of polynomial + 1)
//    const size_t polynomialOrder = 4;
//    const size_t nCoefficients = polynomialOrder + 1;
//    const size_t nSolutions = 2 * polynomialOrder; // real part + imaginary part for each solution
//    
//    // Calculate the coefficients of the polynomial
//    
//    double m1 = 0.0;                                        // Initial primary mass in SI
//    double m2 = 0.0;                                        // Initial secondary mass in SI
//    double Mtot = m1 + m2;                                  // Initial total mass in SI
//    double mu = (m1 * m2)/(Mtot);                           // Initial reduced mass in SI
//    double a = 0.0;                                         // Initial semi-major axis in SI
//    double spin1 = 0.0;                                     // Initial spin angular momenta of object 1
//    double spin2 = 0.0;                                     // Initial spin angular momenta of object 2
//    
//    double m1prime = 0.0;                                   // Final primary mass in SI
//    double m2prime = 0.0;                                   // Final secondary mass in SI
//    double Mtotprime = m1prime + m2prime;                   // Final total mass in SI
//    double muprime = (m1prime * m2prime) / Mtotprime;       // Final reduced mass in SI
//    double aprime = 0.0;                                    // Final semi-major axis in SI
//    double I1prime = 0.0;                                   // Final moment of inertia of object 1
//    double I2prime = 0.0;                                   // Final moment of inertia of object 2
//    
//    // this is just a function of a, m1, m2
//    double omega_orb = 0.0;                                 // Initial orbital frequency
//    
//    double A = sqrt(G * Mtotprime) * muprime;
//    double B = (omega_orb * a*a * mu) + spin1 + spin2;
//    double C = (I1prime * sqrt(G * Mtotprime)) + (I2prime * sqrt(G * Mtotprime));
//    
//    double c_4 = A;
//    double c_3 = B;
//    double c_2 = 0.0;
//    double c_1 = 0.0;
//    double c_0 = C;
//    
//    // Package these coefficients into an array
//    double coefficients[nCoefficients] = {c_0, c_1, c_2, c_3, c_4};
//    double solutions[nSolutions];
//    
//    // Do GSL rootfinding magic
//    gsl_poly_complex_workspace *workspace = gsl_poly_complex_workspace_alloc (nCoefficients);
//    
//    gsl_poly_complex_solve(coefficients, nCoefficients, workspace, solutions);
//    
//    gsl_poly_complex_workspace_free(workspace);
//    
//    // DEBUG: Output the solutions
//    for (int i = 0; i < polynomialOrder; i++){
//        
//        std::cout << "z_" << i << " = " << solutions[2*i] << " + " << solutions[2*i+1] << "i" << std::endl;
//        
//    }
//    
//    // Pick the solution which is real and best conserves orbital energy
//    
//    double y = 0.0;
//    
//    // Return the new semi-major axis and synchronised rotation frequency
//    aPrime = sqrt(y);
//    
//}

//// NOTES


//void locked(){
//
//}

///
//// Prototype function
//void addFive(int &x);
//
//int main(int argc, const char * argv[])
//{
//    
//    int nValue = 5;
//    int *pnPtr = &nValue;
//    
//    cout << "The address of nValue is " << &nValue << endl; // print the address of variable nValue
//    cout << "The value of the pointer pnPtr is " << pnPtr << endl;
//    cout << "The actual value of the thing being pointed at by pnPtr is " << *pnPtr << endl;
//    
//    
//    int a = 0;
//    cout << "Before, a = " << a << endl;
//    addFive(a);
//    cout << "After, a = " << a << endl;}
//
////void f(int& answer);
//
//void addFive(int &x){
//    x += 5;
//}
