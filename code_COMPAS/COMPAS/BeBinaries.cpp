
#include "BeBinaries.h"

// Default constructor
BeBinaries::BeBinaries(){

    // Initialise variables
    cur_props.ID                = -1;
    cur_props.randomSeed        = -1;

    cur_props.dt                = 0;
    cur_props.total_time        = 0;

    cur_props.NS_mass           = NEVER_SET;
    cur_props.NS_spin           = NEVER_SET;

    cur_props.companion_mass    = NEVER_SET;
    cur_props.companion_Lum     = NEVER_SET;
    cur_props.companion_Teff    = NEVER_SET;
    cur_props.companion_R       = NEVER_SET;

    cur_props.companion_accreted= false;

    cur_props.separation        = NEVER_SET;
    cur_props.eccentricity      = NEVER_SET;
}

void BeBinaries::setProps(const BinaryStar &currentBinary){
    /*
    Set the properties for the current binary
    */
    Properties *props;
    prev_props = cur_props;
    props = &cur_props;


    props->ID             = currentBinary.m_ID;
    props->randomSeed     = currentBinary.m_randomSeed;
    props->separation     = currentBinary.m_SemiMajorAxisPrime*AUToRsol;
    props->dt             = currentBinary.m_dt;
    props->total_time    += currentBinary.m_dt;
    props->eccentricity   = currentBinary.m_Eccentricity;


    //It ain't pretty this duplicate code, planned to assign with pointer but didnt
    //work

    //if(currentBinary.star1.m_stellarType == NEUTRON_STAR &
    //  currentBinary.star2.m_stellarType == MAIN_SEQUENCE){
    //   NS        = *currentBinary.star1
    //   companion = *currentBinary.star2
    //}

    if(currentBinary.star1.m_stellarType == NEUTRON_STAR &
       (currentBinary.star2.m_stellarType == MS_MORE_THAN_07 or
        currentBinary.star2.m_stellarType == MS_LESS_THAN_07)){
        props->NS_mass            = currentBinary.star1.m_Mass;
        //props->NS_spin            = currentBinary.star1.m_Mass;
        props->companion_mass     = currentBinary.star2.m_Mass;
        props->companion_Lum      = currentBinary.star2.m_Luminosity;
        props->companion_Teff     = currentBinary.star2.m_Temperature;
        props->companion_R        = currentBinary.star2.m_Radius;
        //props->companion_accreted = false;
    }


    else if( currentBinary.star2.m_stellarType == NEUTRON_STAR &
             (currentBinary.star1.m_stellarType == MS_MORE_THAN_07 or
              currentBinary.star1.m_stellarType == MS_LESS_THAN_07)){
        props->NS_mass            = currentBinary.star2.m_Mass;
        //props->NS_spin            = currentBinary.star1.m_Mass;
        props->companion_mass     = currentBinary.star1.m_Mass;
        props->companion_Lum      = currentBinary.star1.m_Luminosity;
        props->companion_Teff     = currentBinary.star1.m_Temperature;
        props->companion_R        = currentBinary.star1.m_Radius;
        //props->companion_accreted = false;
    }

    else{
        std::cerr << currentBinary.m_randomSeed << 
        "\tError in Be Binaries printing stellar types should not end here (No NS+MS) " << endl;
    }

}




void BeBinaries::averageParameters(){


}



void BeBinaries::printBeBinariesParameters(const boost::filesystem::path &outputPath, std::string filename){
    
    // std::ios_base::app appends to the end of the file.
    std::ofstream BeB((outputPath/filename).string(), std::ios_base::app);

    BeB  <<  cur_props.ID                  <<TAB<<     cur_props.randomSeed     <<TAB<<
             cur_props.dt                  <<TAB<<     cur_props.total_time     <<TAB<<
             cur_props.NS_mass             <<TAB<<     cur_props.NS_spin        <<TAB<<
             cur_props.companion_mass      <<TAB<<     cur_props.companion_Lum  <<TAB<<
             cur_props.companion_Teff      <<TAB<<     cur_props.companion_R    <<TAB<<
             cur_props.companion_accreted  <<TAB<<     cur_props.separation     <<TAB<<
             cur_props.eccentricity        <<  std::endl;
    

    BeB.close();  

}
