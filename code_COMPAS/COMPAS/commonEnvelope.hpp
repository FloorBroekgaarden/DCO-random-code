//
//  commonEnvelope.hpp
//  COMPAS
//
//  Created by Alejandro Vignal-Gomez on 17/02/2016.
//  Copyright © 2016 Simon Stevenson. All rights reserved.
//
//

#ifndef commonEnvelope_hpp
#define commonEnvelope_hpp

// Standard includes
#include <stdio.h>
#include <cmath>

// COMPAS includes
#include "loveridgeCoefficients.h"
#include "programOptions.h"
#include "star.h"
#include "lambdaNanjing.h"
#include "lambdaDewi.h"
#include "lambdaKruckow.h"
#include "lambdaLoveridge.h"
#include "uniformDistribution.h"

// Prototypes
double  doubleCommonEnvelopePhase(double alphaCE, double lambda1, double lambda2, double a, double m1, double mfin1, double menv1, double rRLd1, double m2, double mfin2, double menv2, double rRLd2, int type1, int type2, bool & doubleCommonEnvelopeFlag);
double calculateCommonEnvelopeLambdaParameter(const Star &star, const programOptions &options, unsigned long randomSeed);
double massAccretionNeutronStar(const gsl_rng *r, const programOptions &options, double M_comp, double R_comp) ;

#endif /* commonEnvelope_hpp */
