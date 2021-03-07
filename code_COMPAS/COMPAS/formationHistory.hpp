//
//  Created by Coen Neijssel on 14/11/2016.
//  Copyright © 2016 Simon Stevenson. All rights reserved.
//

#ifndef formationHistory_hpp
#define formationHistory_hpp

#include <fstream>
#include <stdio.h>

#include "BinaryStar.h"


void updateFormationHistory(bool flagRLOFprimary, bool flagRLOFsecondary, bool commonEnvelopeFlag, bool same_RLOF_loop, bool flagExperiencedCCSNPrimary, bool flagExperiencedCCSNSecondary, bool flagExperiencedECSNPrimary, bool flagExperiencedECSNSecondary,int stellar_type_primary, int stellar_type_secondary, int &eventCounter,  int &mt_primary_counter, int &mt_secondary_counter, int &mt_primary_ep1, int &mt_primary_ep1_K1,  int &mt_primary_ep1_K2, int &mt_primary_ep2, int &mt_primary_ep2_K1,  int &mt_primary_ep2_K2, int &mt_primary_ep3, int &mt_primary_ep3_K1, int &mt_primary_ep3_K2, int &mt_secondary_ep1, int &mt_secondary_ep1_K1, int &mt_secondary_ep1_K2, int &mt_secondary_ep2, int &mt_secondary_ep2_K1, int &mt_secondary_ep2_K2, int &mt_secondary_ep3, int &mt_secondary_ep3_K1, int &mt_secondary_ep3_K2,  int &SN_primary_type_1, int &SN_primary_type_2, int &SN_primary_type_3, int &SN_secondary_type_1, int &SN_secondary_type_2, int &SN_secondary_type_3, int &CEE, int &CEE_instigator, int &CEE_failed, int &CEE_failed_instigator, int &CEE_wet, int &CEE_wet_instigator, int &stellar_type_K1,  int &stellar_type_K2, bool m_stellarMerger, int primaryTypePrev, int secondaryTypePrev);

#endif 
