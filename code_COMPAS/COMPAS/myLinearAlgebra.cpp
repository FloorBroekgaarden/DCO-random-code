//
//  myLinearAlgebra.cpp

#include "myLinearAlgebra.h"


// Dot product
double dotProduct3(const double &ax, const double &ay, const double &az, const double &bx, const double &by, const double &bz){
    // Return the dot product a.b
    return ax*bx + ay*by + az*bz;
}

// Cross product
void crossProduct3(const double &ax, const double &ay, const double &az, const double &bx, const double &by, const double &bz, double &cx, double &cy, double &cz){
    // return c = a X b
    cx = ay*bz - az*by;
    cy = az*bx - ax*bz;
    cz = ax*by - ay*bx;
}

// Calculate norm of vector
double norm3(const double &ax, const double &ay, const double &az){
    // Return norm of vector |a| = (a.a)^(1/2)
    return sqrt(dotProduct3(ax, ay, az, ax, ay, az));
}

// Create a unit vector
void unit3(const double &ax, const double &ay, const double &az, double &uax, double &uay, double &uaz){
    
    double norma = norm3(ax, ay, az);
    
    uax = ax/norma;
    uay = ay/norma;
    uaz = az/norma;
    
}

void packageVectorAsArray(double &ax, double &ay, double &az, double (&array)[3][1]){
    
    // Take in the vector a in components and pacakge as an array
    // Can you return an array? Other wise also pass it the array by reference
    
    array[0][0] = ax;
    array[1][0] = ay;
    array[2][0] = az;
    
}


/////////////////////////////////////////////////////////////////////////////
//                                  NOTES
/////////////////////////////////////////////////////////////////////////////



