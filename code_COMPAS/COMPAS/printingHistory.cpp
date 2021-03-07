//
//  printingHistory.cpp

#include "printingHistory.hpp"

/*
std::string chooseHistorySNe(bool flagSNe1, bool flagSNe2, bool flagECSN1, bool flagECSN2){
    std::string temp;
    bool    debugging = false;
    
    std::cout << "went in choose history SNe" << std::endl;
    
    // If Sne, record event
    if(flagSNe1 == true or flagSNe2 == true){
        if(flagSNe1 == true and flagSNe2 == true){
            temp = "->SNe:1&2";
        }
        else if(flagSNe1 == true){
            temp = "->SNe:1";
        }
        else if(flagSNe2 == true){
            temp = "->SNe:2";
        }
        else{
            temp = "";
            std::cerr << "Error in printing supernova history." << std::endl;
        }
    }
    
    if(flagECSN1 == true or flagECSN2 == true){
        
        std::cout << "Went in ecs bit" << std::endl;
        
        if(flagECSN1 == true and flagECSN2 == true){
            temp = "->ECSNe:1&2";
        }
        else if(flagECSN1 == true){
            temp = "->ECSNe:1";
            std::cout << "got here1" << std::endl;
        }
        else if(flagECSN2 == true){
            temp = "->ECSNe:2";
            std::cout << "got here2" << std::endl;
        }
        else{
            temp = "";
            std::cerr << "Error in printing electron capture supernova history." << std::endl;
        }
        
    }
    
    return temp;
}
 */

std::string chooseHistorySNe(int stellarType1, int stellarTypePrev1, int stellarType2, int stellarTypePrev2){
    std::string temp;
    bool    debugging = false;
    bool    flagSNe1 = false;
    bool    flagSNe2 = false;
    
//    debugging = true;
    
    if(debugging){
        std::cout << "went in choose history SNe" << std::endl;
        std::cout << "stellarTypePrev1: " << stellarTypePrev1 << std::endl;
        std::cout << "stellarType1: " << stellarType1 << std::endl;
        std::cout << "stellarTypePrev2: " << stellarTypePrev2 << std::endl;
        std::cout << "stellarType2: " << stellarType2 << std::endl;

    }
    
    if(((stellarType1 == NEUTRON_STAR)and(stellarTypePrev1 <= NAKED_HELIUM_STAR_GIANT_BRANCH))or((stellarType1 == BLACK_HOLE)and(stellarTypePrev1 <= NAKED_HELIUM_STAR_GIANT_BRANCH))){
        flagSNe1 = true;
    }
    
    if(((stellarType2 == NEUTRON_STAR)and(stellarTypePrev2 <= NAKED_HELIUM_STAR_GIANT_BRANCH))or((stellarType2 == BLACK_HOLE)and(stellarTypePrev2 <= NAKED_HELIUM_STAR_GIANT_BRANCH))){
        flagSNe2 = true;
    }
    
    // If Sne, record event
    if(flagSNe1 == true or flagSNe2 == true){
        if(flagSNe1 == true and flagSNe2 == true){
            temp = "->SNe:1&2";
        }
        else if(flagSNe1 == true){
            temp = "->SN:1";
        }
        else if(flagSNe2 == true){
            temp = "->SN:2";
        }
        else{
            temp = "";
            std::cerr << "Error in printing supernova history." << std::endl;
        }
    }
    /*
    if(flagECSN1 == true or flagECSN2 == true){
        
        std::cout << "Went in ecs bit" << std::endl;
        
        if(flagECSN1 == true and flagECSN2 == true){
            temp = "->ECSNe:1&2";
        }
        else if(flagECSN1 == true){
            temp = "->ECSNe:1";
            std::cout << "got here1" << std::endl;
        }
        else if(flagECSN2 == true){
            temp = "->ECSNe:2";
            std::cout << "got here2" << std::endl;
        }
        else{
            temp = "";
            std::cerr << "Error in printing electron capture supernova history." << std::endl;
        }
        
    }
     */
    
    return temp;
}

std::string chooseHistoryMT(bool flagRLOF1, bool flagRLOF2){
    std::string temp;
    bool    debugging = false;

    
    if(flagRLOF1 == true or flagRLOF2 == true){
        if(flagRLOF1 == true and flagRLOF2 == true){
            if(debugging){
                std::cout << "Contact system has been formed -- currently ends evolution" << std::endl;
            }
            temp = "->RLOF:1&2";
        }
        else if(flagRLOF1 == true){
            temp = "->MT:1";
        }
        else if(flagRLOF2 == true){
            temp = "->MT:2";
        }
        else{
            temp = "";
        }
    }
    return temp;
}


std::string chooseHistoryCEE(bool commonEnvelopeFlag, bool wetMergerFlag, bool failedCommonEnvelope){
    std::string temp;
    bool    debugging = false;
    
    // Check if there was a common envelope
    if(commonEnvelopeFlag){
        if(debugging){
            std::cout << "CE flag";
        }
        temp = "->CEE";
    }
    
    // Check if there was a wet merger
    if(wetMergerFlag){
        temp += "->WetMerger";
    }
    
    if(failedCommonEnvelope){
        temp += "->failedCEE";
    }
    
    return temp;
}


std::string chooseHistoryMasslessRemnant(int stellarType1, int stellarType2){
    std::string temp;
    bool    debugging = false;
    
    if(debugging){
        std::cout << "Massless remnant. Stop evolution." << std::endl;
    }
    
    if(stellarType1 == MASSLESS_REMNANT and stellarType2 == MASSLESS_REMNANT){ //Does this even makes sense to ask?
        temp = "->MasslessRemnant:1&2";
    }
    else if(stellarType1 == MASSLESS_REMNANT){
        temp = "->MasslessRemnant:1";
    }
    else if(stellarType2 == MASSLESS_REMNANT){
        temp = "->MasslessRemnant:2";
    }
    else{
        temp = "";
    }
    
    
    return temp;
}


std::string chooseHistoryCoalescence(double a, double R1, double R2, bool failedCommonEnvelope){
// Check if binary components are touching -- print if this happens (should usually be avoided as MT or CE should happen prior to this)
    std::string temp;
    bool    debugging = false;
    
    if(!failedCommonEnvelope){
        if(a <= RsolToAU*(R1+R2) and a > 0.0){
        
            if(debugging){
                std::cout << "a: " << a*AUToRsol << std::endl;
                std::cout << "R1: " << R1 << std::endl;
                std::cout << "R2: " << R2 << std::endl;
                std::cout << "R1+R2: " << R1+R2 << std::endl;
                std::cout << "Too close: coalescence" << std::endl;
            }
            
            temp = "->Coalescence";
        }
    }
    
    return temp;
}


std::string chooseHistorySNeUnbound(double a, double R1, double R2){
    // Check if binary semi major axis is still positive (negative means unbound, presumably by SN)
    std::string temp;
    bool    debugging = false;
    
    if(a <= RsolToAU*(R1+R2) and a > 0.0){
        
        if(debugging){
            std::cout << "a: " << a*AUToRsol << std::endl;
            std::cout << "R1: " << R1 << std::endl;
            std::cout << "R2: " << R2 << std::endl;
            std::cout << "R1+R2: " << R1+R2 << std::endl;
            std::cout << "Binary unbound by supernova" << std::endl;
        }
        
        temp = "->SNe:unbound";
    }
    
    return temp;
}


std::string chooseHistoryGravitationallyUnbound(double E){
    // Check if system remains gravitationally bound
    std::string temp;
    bool    debugging = false;
    
    if(E > 0.0){
        if(debugging){
            std::cout << "System unbound" << std::endl;
            std::cout << E << std::endl;
        }
        
        temp = "->Grav:unbound";
    }
    
    return temp;
}
