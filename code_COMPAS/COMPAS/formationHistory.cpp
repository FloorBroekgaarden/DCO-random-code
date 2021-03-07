#include "formationHistory.hpp"
#include <fstream>

void updateFormationHistory(bool flagRLOFprimary, bool flagRLOFsecondary,bool same_RLOF_loop, bool commonEnvelopeFlag, bool flagExperiencedCCSNPrimary, bool flagExperiencedCCSNSecondary, bool flagExperiencedECSNPrimary, bool flagExperiencedECSNSecondary,int stellar_type_primary, int stellar_type_secondary,  int &eventCounter,    int &mt_primary_counter, int &mt_secondary_counter, int &mt_primary_ep1, int &mt_primary_ep1_K1,  int &mt_primary_ep1_K2, int &mt_primary_ep2, int &mt_primary_ep2_K1,  int &mt_primary_ep2_K2, int &mt_primary_ep3, int &mt_primary_ep3_K1, int &mt_primary_ep3_K2, int &mt_secondary_ep1, int &mt_secondary_ep1_K1, int &mt_secondary_ep1_K2, int &mt_secondary_ep2, int &mt_secondary_ep2_K1, int &mt_secondary_ep2_K2, int &mt_secondary_ep3, int &mt_secondary_ep3_K1, int &mt_secondary_ep3_K2,  int &SN_primary_type_1, int &SN_primary_type_2, int &SN_primary_type_3, int &SN_secondary_type_1, int &SN_secondary_type_2, int &SN_secondary_type_3, int &CEE, int &CEE_instigator, int &CEE_failed, int &CEE_failed_instigator, int &CEE_wet, int &CEE_wet_instigator, int &stellar_type_K1,  int &stellar_type_K2, bool m_stellarMerger, int primaryTypePrev, int secondaryTypePrev){


if((flagRLOFprimary==true)||(flagRLOFsecondary==true)){

    if (same_RLOF_loop ==true){
    }
    else if ((stellar_type_primary !=stellar_type_K1)||(stellar_type_secondary!=stellar_type_K2)){
        if(flagRLOFprimary==true){
            stellar_type_K1 = stellar_type_primary ; 
            stellar_type_K2 = stellar_type_secondary; 
            eventCounter +=1;
            mt_primary_counter +=1;
            if (mt_primary_counter ==1){
               mt_primary_ep1 = eventCounter;
               mt_primary_ep1_K1=  stellar_type_K1;
               mt_primary_ep1_K2=  stellar_type_K2;
            }
            if (mt_primary_counter ==2){
               mt_primary_ep2 = eventCounter;
               mt_primary_ep2_K1=  stellar_type_K1;
               mt_primary_ep2_K2=  stellar_type_K2;
            }
            if (mt_primary_counter ==3){
               mt_primary_ep3 = eventCounter;
               mt_primary_ep3_K1=  stellar_type_K1;
               mt_primary_ep3_K2=  stellar_type_K2;
           }
       } //close check RLOF star1
        
    else if(flagRLOFsecondary==true){
            stellar_type_K1 = stellar_type_primary ; 
            stellar_type_K2 = stellar_type_secondary; 
            eventCounter +=1;
            mt_secondary_counter +=1;
            if (mt_secondary_counter ==1){
               mt_secondary_ep1 = eventCounter;
               mt_secondary_ep1_K1 = stellar_type_K1;
               mt_secondary_ep1_K2 = stellar_type_K2;
            }
            if (mt_secondary_counter ==2){
               mt_secondary_ep2 = eventCounter;
               mt_secondary_ep2_K1 = stellar_type_K1;
               mt_secondary_ep2_K2 = stellar_type_K2;
            }
            if (mt_secondary_counter ==3){
               mt_secondary_ep3 = eventCounter;
               mt_secondary_ep3_K1 = stellar_type_K1;
               mt_secondary_ep3_K2 = stellar_type_K2;
            }
        } //close check RLOF star1
    }//close check stellar type change
}//close RLOF check


if(commonEnvelopeFlag == true or flagExperiencedCCSNPrimary == true or flagExperiencedCCSNSecondary or
	flagExperiencedECSNPrimary == true or flagExperiencedECSNSecondary){
//TODO COEN why not break up CEE and SN in two if loops also in BinaryStar.cpp?
    if (m_stellarMerger){
        eventCounter +=1;
        CEE_failed = eventCounter;
        if (flagRLOFprimary && flagRLOFsecondary){
            CEE_failed_instigator =3;}
        else if (flagRLOFprimary){
            CEE_failed_instigator =1;}
        else if (flagRLOFsecondary){
            CEE_failed_instigator =2;}
    }

    else if(commonEnvelopeFlag){
        eventCounter +=1;
        CEE =eventCounter;
        if (flagRLOFprimary && flagRLOFsecondary){
            CEE_instigator =3;}
        else if (flagRLOFprimary){
            CEE_instigator =1;}
        else if (flagRLOFsecondary){
            CEE_instigator =2;}
    }
    else if((flagExperiencedCCSNPrimary)&&(SN_primary_type_1==false)&&(SN_primary_type_2==false)&&(SN_primary_type_3==false)){
        eventCounter +=1;
        SN_primary_type_1 =eventCounter;
    }
    else if((flagExperiencedCCSNSecondary)&&(SN_secondary_type_1==false)&&(SN_secondary_type_2==false)&&(SN_secondary_type_3==false)){
        eventCounter +=1;
        SN_secondary_type_1 =eventCounter;
    }
    else if((flagExperiencedECSNPrimary)&&(SN_primary_type_1==false)&&(SN_primary_type_2==false)&&(SN_primary_type_3==false)){
        eventCounter +=1;
        SN_primary_type_2 =eventCounter;
    }
    else if((flagExperiencedECSNSecondary)&&(SN_secondary_type_1==false)&&(SN_secondary_type_2==false)&&(SN_secondary_type_3==false)){
        eventCounter +=1;
        SN_secondary_type_2 =eventCounter;
    }
}




}




