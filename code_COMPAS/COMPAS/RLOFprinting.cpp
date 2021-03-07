
#include "RLOFprinting.h"

// Default constructor
RLOFprinting::RLOFprinting(){

    // Initialise variables
    cur_props.m_ID =0;
    cur_props.m_randomSeed=0;
    cur_props.m_mass1=0.0;
    cur_props.m_mass2=0.0;
    cur_props.m_radius1=0.0;
    cur_props.m_radius2=0.0;
    cur_props.m_type1=0;
    cur_props.m_type2=0;
    cur_props.m_separation=0.0;
    cur_props.m_eventCounter=0;
    cur_props.m_time=0.0;
    cur_props.m_flagRLOF1=false;
    cur_props.m_flagRLOF2=false;
    cur_props.m_printMS1=false;
    cur_props.m_printMS2=false;
    cur_props.m_printHeMS1=false;
    cur_props.m_printHeMS2=false;
    cur_props.m_CEEflag=false;

}

void RLOFprinting::setProps(const BinaryStar &currentBinary, bool Current){
    /*
    Set the properties for the current binary
    */
    Properties *props;
    bool setParams= true;
    bool dummy1= false;
    bool dummy2= false;
    if(Current){
        props = &cur_props;
        dummy1 = monitorStartCaseA(currentBinary);
    }
    else{
        props = &prev_props;
        if((cur_props.m_printMS1)or(cur_props.m_printMS2)or
           (cur_props.m_printHeMS1)or(cur_props.m_printHeMS2)){
            setParams=false;
        }//dont set if we are monitoring caseA
    }
    if(setParams){
        props->m_ID           =currentBinary.m_ID;
        props->m_randomSeed   =currentBinary.m_randomSeed;
        props->m_mass1        =currentBinary.star1.m_Mass;
        props->m_mass2        =currentBinary.star2.m_Mass;
        props->m_radius1      =currentBinary.star1.m_Radius;
        props->m_radius2      =currentBinary.star2.m_Radius;
        props->m_type1        =currentBinary.star1.m_stellarType;
        props->m_type2        =currentBinary.star2.m_stellarType;
        props->m_separation   =currentBinary.m_SemiMajorAxisPrime*AUToRsol;
        props->m_time         =currentBinary.m_time;
        props->m_flagRLOF1    =currentBinary.star1.flagRLOF;
        props->m_flagRLOF2    =currentBinary.star2.flagRLOF;
        props->m_CEEflag      =currentBinary.m_commonEnvelopeFlag;
    }
}


bool RLOFprinting::monitorStartCaseA(const BinaryStar &currentBinary){
    
    bool setPrev = false;

    bool isMS = false;

    if(currentBinary.star1.m_stellarType == MS_LESS_THAN_07 or currentBinary.star1.m_stellarType==MS_MORE_THAN_07){
        isMS = true;
    }

    //TODO reduce copy code by pointer star1 vs star2 similar to setProps function
    if(currentBinary.star1.flagRLOF){
        if (isMS and cur_props.m_printMS1 == false){
            cur_props.m_printMS1=true;
            setPrev = true;
        }
        else if(currentBinary.star1.m_stellarType == NAKED_HELIUM_STAR_MS and cur_props.m_printHeMS1 == false){
            cur_props.m_printHeMS1=true;
            setPrev = true;
        }
    }
    if(currentBinary.star2.flagRLOF){
        if (isMS and cur_props.m_printMS2 == false){
            cur_props.m_printMS2=true;
            setPrev = true;
        }
        else if(currentBinary.star2.m_stellarType == NAKED_HELIUM_STAR_MS and cur_props.m_printHeMS2 == false){
            cur_props.m_printHeMS2=true;
            setPrev = true;
        }
    }

    return setPrev;
}


bool RLOFprinting::monitorEndCaseA(const BinaryStar &currentBinary){

    bool setEnd = false;
    
    if(((currentBinary.star1.m_stellarType != MS_LESS_THAN_07) and
        (currentBinary.star1.m_stellarType != MS_MORE_THAN_07)) and
        (cur_props.m_printMS1==true)){
        cur_props.m_printMS1=false;
        setEnd = true;
    }
    if(((currentBinary.star2.m_stellarType != MS_LESS_THAN_07) and
        (currentBinary.star2.m_stellarType != MS_MORE_THAN_07)) and
        (cur_props.m_printMS2==true)){
        cur_props.m_printMS2=false;
        setEnd = true;
    }   
    if((currentBinary.star1.m_stellarType != NAKED_HELIUM_STAR_MS) and
        (cur_props.m_printHeMS1==true)){
        cur_props.m_printHeMS1=false;
        setEnd = true;
    }
    if((currentBinary.star2.m_stellarType != NAKED_HELIUM_STAR_MS) and
        (cur_props.m_printHeMS2==true)){
        cur_props.m_printHeMS2=false;
        setEnd = true;
    }
    return setEnd;

}



bool RLOFprinting::printParameters(const BinaryStar &currentBinary){

    bool print = false;
    
    //If not MS mass transfer just print
    //TODO currently anti bug is also if stellartype < WD
    if(currentBinary.star1.flagRLOF){
        if ((prev_props.m_type1!=MS_LESS_THAN_07) and
            (prev_props.m_type1!=MS_MORE_THAN_07) and
            (prev_props.m_type1!=NAKED_HELIUM_STAR_MS) and
            (prev_props.m_type1<HELIUM_WHITE_DWARF)){
            print = true;
            cur_props.m_printMS1=false;
            cur_props.m_printHeMS1=false;
        }
    }
    if(currentBinary.star2.flagRLOF){
        if ((prev_props.m_type2!=MS_LESS_THAN_07) and
            (prev_props.m_type2!=MS_MORE_THAN_07) and
            (prev_props.m_type2!=NAKED_HELIUM_STAR_MS) and
            (prev_props.m_type2<HELIUM_WHITE_DWARF)){
            print = true;
            cur_props.m_printMS2=false;
            cur_props.m_printHeMS2=false;
        }
    }

    bool printEndCaseA = false;
    printEndCaseA = monitorEndCaseA(currentBinary);
    if(printEndCaseA){
        print = true;}
    
    return print;

}

void RLOFprinting::printRLOFParameters(const boost::filesystem::path &outputPath, std::string filename, const BinaryStar &currentBinary){
    
    // std::ios_base::app appends to the end of the file.
    std::ofstream RLOF((outputPath/filename).string(), std::ios_base::app);
    
    RLOF << cur_props.m_ID          <<TAB<<   cur_props.m_randomSeed    <<TAB<<
            cur_props.m_mass1       <<TAB<<   cur_props.m_mass2         <<TAB<<
            cur_props.m_radius1     <<TAB<<   cur_props.m_radius2       <<TAB<<
            cur_props.m_type1       <<TAB<<   cur_props.m_type2         <<TAB<<
            cur_props.m_separation  <<TAB<<   cur_props.m_eventCounter  <<TAB<<
            cur_props.m_time        <<TAB<<   cur_props.m_flagRLOF1     <<TAB<<
            cur_props.m_flagRLOF2   <<TAB<<   cur_props.m_CEEflag       <<TAB<<
            prev_props.m_mass1      <<TAB<<   prev_props.m_mass2        <<TAB<<
            prev_props.m_radius1    <<TAB<<   prev_props.m_radius2      <<TAB<<
            prev_props.m_type1      <<TAB<<   prev_props.m_type2        <<TAB<<
            prev_props.m_separation <<TAB<<   prev_props.m_eventCounter <<TAB<<
            prev_props.m_time       <<TAB<<   prev_props.m_flagRLOF1    <<TAB<<
            prev_props.m_flagRLOF2  <<TAB<<   

            currentBinary.star1.m_zetaThermal       <<TAB<<
            currentBinary.star1.m_zetaNuclear       <<TAB<<
            currentBinary.star1.m_zetaSoberman      <<TAB<<
            currentBinary.star1.m_zetaSobermanHelium<<TAB<<
            currentBinary.star1.m_zetaHurley        <<TAB<<
            currentBinary.star1.m_zetaHurleyHelium  <<TAB<<
            currentBinary.star1.m_zetaSimple        <<TAB<<

            currentBinary.star2.m_zetaThermal       <<TAB<<
            currentBinary.star2.m_zetaNuclear       <<TAB<<
            currentBinary.star2.m_zetaSoberman      <<TAB<<
            currentBinary.star2.m_zetaSobermanHelium<<TAB<<
            currentBinary.star2.m_zetaHurley        <<TAB<<
            currentBinary.star2.m_zetaHurleyHelium  <<TAB<<
            currentBinary.star2.m_zetaSimple        <<TAB<<

            currentBinary.m_zetaRLOFanalytic        <<TAB<<
            currentBinary.m_zetaRLOFnumerical       << std::endl;

    //every time we print a MT event happend
    cur_props.m_eventCounter +=1;

    RLOF.close();  

}
