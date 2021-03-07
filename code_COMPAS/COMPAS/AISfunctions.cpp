//
//  AISfunctions.cpp
//  COMPAS
//
//  Created by Floor Broekgaarden on 24/04/2018
//  
//

#include "AISfunctions.h"


void printAISexploratorySettings(const programOptions &options){
	// function that prints the selections of the 'hits' of the exploratory phase of the Adaptive Importance sampling
	std::cout << std::endl;
	std::cout << " ------------------------------------------------------------- " << std::endl;
    std::cout << " Running the Adaptive Importance Sampling exploratory phase " << std::endl; 
    
    // SIMON: I replaced options.AISDCOtypeString == "ALL" with boost::iequals(options.AISDCOtypeString, "ALL") which should be safer/less error prone
    if(boost::iequals(options.AISDCOtypeString, "ALL")){
    	std::cout << " - Selecting all DCOs    " << std::endl;
    }
    else if(boost::iequals(options.AISDCOtypeString, "BBH")){
    	std::cout << " - Selecting all BBH mergers    " << std::endl;
    }
    else if(boost::iequals(options.AISDCOtypeString, "BNS")){
    	std::cout << " - Selecting all BNS mergers    " << std::endl;
    }
    else if(boost::iequals(options.AISDCOtypeString, "BHNS")){
    	std::cout << " - Selecting all BHNS mergers    " << std::endl;
    }

    if(options.AISHubble == 1){
        std::cout << " - Selecting only binaries that merge in Hubble time    " << std::endl;
    }	    
    else{
    	std::cout << " - Selecting binaries that both do and do not merge in Hubble time   " << std::endl;
    }

    if(options.AISRLOF == 1){
        std::cout << " - Excluding binaries that have RLOFSecondaryAfterCEE == 1    " << std::endl;
    }
    else{
    	std::cout << " - Not excluding binaries that have RLOFSecondaryAfterCEE == 1   " << std::endl;
    }

    if(options.AISPessimistic){
        std::cout << " - Selecting only Pessimistic CE binaries    " << std::endl;
    }	    
    else{
    	std::cout << " - Selecting Optimistic CE binaries   " << std::endl;
    }

    std::cout << " ------------------------------------------------------------- " << std::endl;
    std::cout << std::endl;	
}


int calculateDCOhit(const programOptions &options, const BinaryStar &currentBinary){
	// Function that checks if we have found a hit, it returns DCOhit which equals 0 if there is no hit and 1 if there is a hit
	// it checks it with the corresponding settings of the eploratory phase in Adaptive Importance Sampling that are defined in pythonSubmit

	int DCOhit = 0;

// Check if we have found a DCO of the DCO type of interest and test for other filters

    if(boost::iequals(options.AISDCOtypeString, "ALL")){
        // we have a possible hit if binary is DCO 

        DCOhit = 1;
        // Only if we found a DCO of interest, check the other flags
        // If we want to exclude DCOs that do not merge in Hubble, set DCOhit = 0 for mergers not merging in Hubble
        if(options.AISHubble){
            if (currentBinary.m_mergesInHubbleTimeFlag ==0){
            	DCOhit = 0;
            } 
        }
        // If we want to exclude DCOs that have RLOFafterSecondaryAfterCEE, set DCOhit = 0 for mergers RLOFafterSecondaryAfterCEE =1
        if(options.AISRLOF){
            if (currentBinary.m_RLOFSecondaryAfterCEE ==1){
                DCOhit = 0; 
            }
        }
        // If we want to exclude DCOs that have optimisticCEFlag, set DCOhit = 0 for mergers optimisticCEFlag =1
        if(options.AISPessimistic){
            if (currentBinary.m_optimisticCommonEnvelopeFlag ==1){
                DCOhit = 0; 
            }
        }     
    }

    else if(boost::iequals(options.AISDCOtypeString,"BBH")){
        // we have a possible hit if binary is BBH 
        if (currentBinary.star1.m_stellarType  == BLACK_HOLE and currentBinary.star2.m_stellarType == BLACK_HOLE){
            DCOhit = 1;
            // Only if we found a DCO of interest, check the other flags
            // If we want to exclude DCOs that do not merge in Hubble, set DCOhit = 0 for mergers not merging in Hubble
            if(options.AISHubble){
                if (currentBinary.m_mergesInHubbleTimeFlag ==0){
                    DCOhit = 0;
                } 
            }
            // If we want to exclude DCOs that have RLOFafterSecondaryAfterCEE, set DCOhit = 0 for mergers RLOFafterSecondaryAfterCEE =1
            if(options.AISRLOF){
                if (currentBinary.m_RLOFSecondaryAfterCEE ==1){
                    DCOhit = 0; 
                }
            }
            // If we want to exclude DCOs that have optimisticCEFlag, set DCOhit = 0 for mergers optimisticCEFlag =1
            if(options.AISPessimistic){
                if (currentBinary.m_optimisticCommonEnvelopeFlag ==1){
                    DCOhit = 0; 
                }
            }     
        }
    }
    else if(boost::iequals(options.AISDCOtypeString,"BNS")){
        // we have a possible hit if binary is BNS 
        if (currentBinary.star1.m_stellarType == NEUTRON_STAR  and currentBinary.star2.m_stellarType == NEUTRON_STAR){
            DCOhit = 1;
            // Only if we found a DCO of interest, check the other flags
            // If we want to exclude DCOs that do not merge in Hubble, set DCOhit = 0 for mergers not merging in Hubble
            if(options.AISHubble){
                if (currentBinary.m_mergesInHubbleTimeFlag ==0){
                    DCOhit = 0;
                } 
            }
            // If we want to exclude DCOs that have RLOFafterSecondaryAfterCEE, set DCOhit = 0 for mergers RLOFafterSecondaryAfterCEE =1
            if(options.AISRLOF){
                if (currentBinary.m_RLOFSecondaryAfterCEE ==1){
                    DCOhit = 0; 
                }
            }
            // If we want to exclude DCOs that have optimisticCEFlag, set DCOhit = 0 for mergers optimisticCEFlag =1
            if(options.AISPessimistic){
                if (currentBinary.m_optimisticCommonEnvelopeFlag ==1){
                    DCOhit = 0; 
                }
            }     
        }
    }
    else if(boost::iequals(options.AISDCOtypeString,"BHNS")){
        // we have a possible hit if binary is BHNS 
        if ((currentBinary.star1.m_stellarType == NEUTRON_STAR and currentBinary.star2.m_stellarType == BLACK_HOLE) or (currentBinary.star1.m_stellarType == BLACK_HOLE and currentBinary.star2.m_stellarType == NEUTRON_STAR)){
            DCOhit = 1;
                        
            // Only if we found a DCO of interest, check the other flags
            // If we want to exclude DCOs that do not merge in Hubble, set DCOhit = 0 for mergers not merging in Hubble
            if(options.AISHubble){
                if (currentBinary.m_mergesInHubbleTimeFlag ==0){
                    DCOhit = 0;
                } 
            }
            // If we want to exclude DCOs that have RLOFafterSecondaryAfterCEE, set DCOhit = 0 for mergers RLOFafterSecondaryAfterCEE =1
            if(options.AISRLOF){
                if (currentBinary.m_RLOFSecondaryAfterCEE ==1){
                    DCOhit = 0; 
                }
            }
            // If we want to exclude DCOs that have optimisticCEFlag, set DCOhit = 0 for mergers optimisticCEFlag =1
            if(options.AISPessimistic){
                if (currentBinary.m_optimisticCommonEnvelopeFlag ==1){
                    DCOhit = 0; 
                }
            }     
        }
    }
    return DCOhit;
}


void defineAISgaussians(AISvariables & aisvariables){ // const AISvariables
    // Function that reads in the gaussian means and covariances  from STROOPWAFELgaussians.txt and returns vectors mu_x and cov_x 
    // that will be given to the distribution functions in BinaryStar such that one can draw random samples for Adaptive Importance Sampling (step 2)

    double  mu1, mu2, mu3, cov1, cov2, cov3; // temp storage of data
    double n_Gaussians; // number of Gaussians in STROOPWAFELgaussians.txt

    std::ifstream infile;
    std::string line;
    int num = 0;
    
	if ( !boost::filesystem::exists( "STROOPWAFELgaussians.txt" ) )
	{
	  std::cout << "Can't find the STROOPWAFELgaussians.txt file!" << std::endl;
	}
	else{
	    infile.open("STROOPWAFELgaussians.txt");// file containing means and covariances of Gaussians in 6 columns
	    if(infile.fail()) // checks to see if file opended
	    { 
	        std::cout << "Error opening STROOPWAFELgaussians.txt" << std::endl;
	    }

	    while(getline(infile, line))  // reads file to end of file
	    { 
	        if(num==0){
	            // read in the first column which contains the nr of Gaussians (i.e. hits) used
	            ++num;
	            infile >>  n_Gaussians ; //
                std::cout << n_Gaussians << "  nGaussians   " << std::endl;
	            // set size of vectors to the found length

                aisvariables.mu_M1.resize(n_Gaussians);
                aisvariables.mu_loga.resize(n_Gaussians); 
                aisvariables.mu_q.resize(n_Gaussians); 
                aisvariables.cov_M1.resize(n_Gaussians);
                aisvariables.cov_loga.resize(n_Gaussians);
                aisvariables.cov_q.resize(n_Gaussians);
	        }   
	        
            if(num==n_Gaussians+1){
                // this is a hack since while loop will continue +1 beyond nr of gaussians
                // this makes sure we dont read in last line twice (or get segmentation fault)
            }
  
	        else{
	            // read colums of STROOPWAFELgaussians.txt into vectors

	            infile >> mu1; 
	            infile >> mu2; 
	            infile >> mu3; 
	            infile >> cov1;
	            infile >> cov2;
	            infile >> cov3;

	            aisvariables.mu_M1[num-1] = mu1;   
	            aisvariables.mu_loga[num-1] = mu2;
	            aisvariables.mu_q[num-1] = mu3;
	            aisvariables.cov_M1[num-1] = cov1;
	            aisvariables.cov_loga[num-1] = cov2;
	            aisvariables.cov_q[num-1] = cov3; 

	            ++num;
	        }
	    }

	    std::cout << " ------------------------------------------------------------- " << std::endl;
	    std::cout << " Started Adaptive Importance Sampling (step 2)  "                 << std::endl;
	    std::cout << " COMPAS will use an instrumental distribution based on " << aisvariables.cov_q.size() << " Gaussian distributions, i.e. 'hits'.  " <<std::endl;
	    std::cout << " ------------------------------------------------------------- " << std::endl;
	}
}


void updateFexpl(const programOptions &options, AISvariables & aisvariables, int i){
    // Function that updates the fraction that should be spend on the exploratory phase fexpl 
    // this is Equation 16 in Broekgaarden+18  

    // estimated weight of the unidentified target population region 
    // uses previously found fexplAIS to estimate z2

    if(options.nBatchesUsed > 0){
        double z2 = 1.0 / (aisvariables.fexplAIS * options.nBatchesUsed * float(options.nBinaries)); // estimated weight of unidentified region
        double z1 = float(aisvariables.CounterDCOsAIS) / float(i); // the estimated rate of the target population region

        double numerator = z1* (sqrt(1. - z1) - sqrt(z2));  // numerator of eq 16
        double denominator = sqrt(1. - z1) * ( sqrt(z2 * (1. - z1)) + z1); // dominator of Eq 16
        
        aisvariables.fexplAIS = 1. - (float(numerator) / denominator); 
    }
    else{
        std::cerr << "\tError in nBatches. Can't be negative.\n";
    }
}

