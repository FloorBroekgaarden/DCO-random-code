//
//  printingHistory.hpp
//  binaries
//
//  Created by Alejandro Vignal-Gomez on 11/07/2016.
//  Copyright © 2016 Simon Stevenson. All rights reserved.
//

#ifndef printingHistory_hpp
#define printingHistory_hpp

#include <stdio.h>

// My includes
#include "constants.h"                  // A header file containing mathematical and physical constants
#include "programOptions.h"
#include "star.h"

// Prototypes
//std::string chooseHistorySNe(bool flagSNe1, bool flagSNe2, bool flagECSN1, bool flagECSN2);
std::string chooseHistorySNe(int stellarType1, int stellarTypePrev1, int stellarType2, int stellarTypePrev2);
std::string chooseHistoryMT(bool flagRLOF1, bool flagRLOF2);
std::string chooseHistoryCEE(bool commonEnvelopeFlag, bool wetMergerFlag, bool failedCommonEnvelope);
std::string chooseHistoryMasslessRemnant(int stellarType1, int stellarType2);
std::string chooseHistoryCoalescence(double a, double R1, double R2, bool failedCommonEnvelope);
std::string chooseHistorySNeUnbound(double a, double R1, double R2);
std::string chooseHistoryGravitationallyUnbound(double E);

#endif /* printingHistory_hpp */
