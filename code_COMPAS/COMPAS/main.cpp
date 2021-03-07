//
//  main.cpp
//

// Standard Includes
#include <iostream>
#include <fstream>
#include <ctime>
#include <string>
#include <stdlib.h>
#include <vector>

// COMPAS Includes
#include "BinaryStar.h"                 // A header file containing the class for the binary
#include "constants.h"                  // A header file containing mathematical and physical constants
#include "eIntegral.h"                  // A header file containing the time to coalescence eccentricity integral
#include "timeToCoalescence.h"          // A header file containing the time to coalescence calculation using interpolated values
#include "programOptions.h"             // A header file containing the class to store program options
#include "programOptionsSorter.h"       // A header file containing a function that gets options from commmand line and assigns them to programOptions
#include "rocheLobeRadius.hpp"          // A header file containing a function that computes roche lobe radius
#include "printingHistory.hpp"
#include "printing.h"                   //A header file containing all different printing functions.
#include "RLOFprinting.h"
#include "BeBinaries.h"
#include "AISfunctions.h"                // file containing functions for Adaptive Importance Sampling  // Floor 
#include "AISclass.h"

// GSL Includes
#include <gsl/gsl_rng.h>                // GSL random number generator
#include <gsl/gsl_spline.h>             // GSL splines for interpolation
#include <gsl/gsl_errno.h>              // GSL error handling
#include <gsl/gsl_interp.h>             // GSL interpolation

// Boost Includes
#include <boost/program_options.hpp>    // Boost command line options tools
#include <boost/filesystem.hpp>         // Boost filesystem tools for handling paths etc.

using namespace std;

// Should put these in a header file somewhere
void COMPASSingle(const programOptions &options, const gsl_rng *r);
void COMPASBinary(programOptions &options, gsl_rng *r);

int main(int argc, char * argv[]){
    // Create instance of class to hold program options
    programOptions options = programOptions();
    // Integer to hold status of program after reading command line - decides whether to run main program
    int programStatus = 0;              // Whether or not to continue running program after sorting command line options

    // Get command line options
    commandLineSorter(argc, argv, programStatus, options);
    if(programStatus == CONTINUE){

        ////////////////////////////////////////////////////////////////////
        //              SETUP GSL RNG
        ////////////////////////////////////////////////////////////////////

        // Set the generator and seed using environment variables (use GSL_RNG_SEED to seed)
        gsl_rng_env_setup();

        // Seed the random number generator
        if(!getenv("GSL_RNG_SEED")){
            gsl_rng_default_seed = time(0);

        }

//        // Alternatively could call gsl_rng_set(r, time(0)) after gsl_rng_alloc if you don't want to allow the user to override the seed with their environment variable.
//        if(!options.fixedRandomSeed){
//            cout << "default seed: " << gsl_rng_default_seed << endl;
//        }

        // Declare Variables
        static gsl_rng *r=NULL;

        // Allocate memory to RNG
        if(!r){
            if (!options.quiet) cout << "Allocating rng." << endl;
            r = gsl_rng_alloc (gsl_rng_default);
        }

        // if single do single else do binary
        if(options.singleStar){
            COMPASSingle(options, r);
        }
        else{
            COMPASBinary(options, r);
        }


        // De-allocate memory to rng 'r'
        gsl_rng_free (r);

        // Make me happy -- Always do, baby ;)
        if (!options.quiet) cout << endl << "Success!" << endl;

        programStatus = SUCCESS;

    }

    return programStatus;

}

void COMPASSingle(const programOptions &options, const gsl_rng *r){
    /*
    A function to evolve ONLY single stars
    */
    //cout << "Does nothing at the moment" << endl;
    
    ofstream outs;
    
    // Declare variables -- Loop over masses to evolve
    int nStepsMass = 100;                                       // Make into an option
    double MassMin = 5.0;                                       // Make into an option
    double MassMax = 100.0;                                     // Make into an option
    double stepMass = (MassMax - MassMin)/nStepsMass;  
    double Mass_to_evolve = 0.0;

    // Which metallicity to use
    double Metallicity_to_evolve = options.metallicity; // 0.001 // 0.02
    
    bool useFakeTimestep = false;
    
    // double dt = 0.005;
    double dt = 0.0;
    
    // double radius_proxy_mass_transfer = 200.0;
    double timestep_reduction_factor=1.0; // Should put this as an option in COMPAS
    
    for(int i =0; i<=nStepsMass; i++){



        Mass_to_evolve = MassMin + (i*stepMass);

        cout << "this mass = " << Mass_to_evolve << endl;
        
        string base_filename = options.outputPathString + "/COMPAS_vink_winds_five_percent_solar_mass_";

        string extension = ".txt";
        
        string Number;                          // string which will contain the result
        ostringstream convert;                  // stream used for the conversion
        convert << i;                           // insert the textual representation of the i-th binary
        Number = convert.str();
        
        string this_filename = base_filename + Number + extension;
        
        cout << "this filename = " << this_filename << endl;
        
        // Open this output
        outs.open(this_filename);
        
        // Header for the output file
        outs << "#Age\tdt\tTime\tType\tMass\tMass0\tRadius\tRZAMS\tLuminosity\tTemperature\tCoreMass\tHeCoreMass\tCOCoreMass\tMdot\ttMS" << endl;
        
        // Create a star
        Star star_to_evolve = Star(Mass_to_evolve, Metallicity_to_evolve, options, r);
        
        // Initialise the star
        star_to_evolve.evolveOneTimestep(dt, useFakeTimestep, options, r); // fixes bug with fact that timescales aren't calculated when star is initialised (should fix)
        
        // Evolve star until end of its life
        while(star_to_evolve.m_stellarType < HELIUM_WHITE_DWARF){
            
            dt = star_to_evolve.makeTimestep(options, r);
            
            if(star_to_evolve.m_stellarType > MS_MORE_THAN_07){
                dt = dt/timestep_reduction_factor;
            }
            
            star_to_evolve.evolveOneTimestep(dt, useFakeTimestep, options, r); //dt
            
            bool debugging = false;
            
            // Apply wind mass loss
            if(options.useMassLoss){

                star_to_evolve.m_Mass = star_to_evolve.massLossCaller(star_to_evolve.m_Mass, star_to_evolve.m_Mass0, star_to_evolve.m_coreMass, star_to_evolve.m_Luminosity, star_to_evolve.m_Radius, star_to_evolve.m_Mu, star_to_evolve.m_Metallicity, star_to_evolve.m_Temperature, star_to_evolve.m_logMetallicityXi, star_to_evolve.m_Mdot, star_to_evolve.m_dt, star_to_evolve.m_Age, options.luminousBlueVariableFactor, options.wolfRayetFactor, star_to_evolve.m_stellarType, options.useMassLoss, false, options.massLossPrescription, debugging, star_to_evolve.m_an_coefficients, star_to_evolve.m_timescales);

                double agePrime1 = star_to_evolve.m_Age;
            
                double f_rej_1 = massTransferRejuvenationFactor(star_to_evolve, options);
            
                star_to_evolve.m_Age = agePrime1 * f_rej_1;
                
            }
            
            // Write output
            cout << star_to_evolve.m_Age << "\t" << dt << "\t" << star_to_evolve.m_time << "\t" << star_to_evolve.m_stellarType << "\t" << star_to_evolve.m_Mass << "\t" << star_to_evolve.m_Mass0 << "\t" << star_to_evolve.m_Radius << "\t" << star_to_evolve.m_RZAMS << "\t" << star_to_evolve.m_Luminosity << "\t" << star_to_evolve.m_Temperature << "\t" << star_to_evolve.m_coreMass << "\t" << star_to_evolve.m_HeCoreMass << "\t" << star_to_evolve.m_COCoreMass << "\t" <<star_to_evolve.m_Mdot << "\t" << star_to_evolve.m_timescales[0] << endl;
            
            outs << star_to_evolve.m_Age << "\t" << dt << "\t" << star_to_evolve.m_time << "\t" << star_to_evolve.m_stellarType << "\t" << star_to_evolve.m_Mass << "\t" << star_to_evolve.m_Mass0 << "\t" << star_to_evolve.m_Radius << "\t" << star_to_evolve.m_RZAMS << "\t" << star_to_evolve.m_Luminosity << "\t" << star_to_evolve.m_Temperature << "\t" << star_to_evolve.m_coreMass << "\t" << star_to_evolve.m_HeCoreMass << "\t" << star_to_evolve.m_COCoreMass << "\t" <<star_to_evolve.m_Mdot << "\t" << star_to_evolve.m_timescales[0] << endl;
            
        }
    
        // Close the output file
        outs.flush();
        outs.close();

        
    }

}

    
void COMPASBinary(programOptions &options, gsl_rng *r){
    /*
    A function to evolve ONLY binaries
    */

    // Declare variables
    int binaryStatus = 0;               // Whether or not to evolve binary

    // SETUP ASCII OUTPUT FILE STREAMS
    // Create string containing filename
    string systemParametersFile         = "systemParameters.txt";
    string doubleCompactObjectsFile     = "doubleCompactObjects.txt";
    string formationHistoryFile         = "formationHistory.txt";
    string supernovaFile                = "supernovae.txt";
    string commonEnvelopesFile          = "commonEnvelopes.txt";
    string RLOFFile                     = "RLOF.txt";
    string BeBinariesFile               = "BeBinaries.txt";
    string pulsarEvolutionFile          = "pulsarEvolution.txt";
    string rejectedSamplesFile          = "rejectedSamples.txt";  // rejectedSampling


    // Print headers for output files -- Update these when you add new parameters to the output! 
	//(we should also make this a function since used in multiple place)
	printSystemHeader(options.outputPath, systemParametersFile);
	printFormationHeader(options.outputPath, formationHistoryFile);
    printSupernovaeHeader(options.outputPath, supernovaFile);
    printCommonEnvelopeHeader(options.outputPath, commonEnvelopesFile);
	printDoubleCompactObjectsHeader(options.outputPath, doubleCompactObjectsFile);
       
    
    // File streams
    // std::ios_base::app appends to the end of the file after the headers.
    //ofstream systems((options.outputPath/systemParametersFile).string(), std::ios_base::app); 
    ofstream doubleCompactObjects((options.outputPath/doubleCompactObjectsFile).string(), std::ios_base::app);
    //ofstream data;

    if(options.BeBinaries){
		printBeBinariesHeader(options.outputPath, "BeBinaries.txt");
	}
    if(options.CHEvolution){
        printCHEvolutionHeader(options.outputPath, "CHEvolutionParameters.txt");
    }
    if(options.RLOFPrinting){
        printRLOFHeader(options.outputPath, RLOFFile);
    }
    if(options.evolvePulsars){
        printpulsarEvolutionHeader(options.outputPath, pulsarEvolutionFile);
    }
    if(options.RejectedSamplesPrinting){
        printRejectedSamplesHeader(options.outputPath, rejectedSamplesFile); 
    }

    // GENERATE AND EVOLVE BINARIES
    if (!options.quiet) cout << "Now generating binaries." << endl;



    // Floor 24/04/2018 - print the selected options for Adaptive Importance Sampling Exploratory phase in the beginning of the run
    if(options.AISexploratoryphase){
        printAISexploratorySettings(options);   
    }

    //Floor 24/05/2018 - Adaptive Importance Sampling step 2: class AISvariables contains parameters for sampling from Gaussians. 
    // In particular it contains vectors that define the gaussians for AIS step 2 aisvariables 
    AISvariables aisvariables = AISvariables();
    if(options.AISrefinementPhase){
        //  If we are sampling using Adaptive Importance Sampling (step 2):read in gaussians 
        defineAISgaussians(aisvariables); 
    }
    
    // Time the program
    clock_t clockStart;
    clockStart  = clock();

    // Generate and evolve binaries
    for(int i = 0; i < options.nBinaries; i++){
        // By default we will evolve this binary
        binaryStatus = CONTINUE;

        // Generate a random binary according to the user options
        // BinaryStar currentBinary = BinaryStar();
        // Floor : added means and covariances as variables in case using Adaptive Importance Sampling
        BinaryStar currentBinary = BinaryStar(r, options, aisvariables); // Initialise current binary using the constructor 
        

        RLOFprinting RLOFp = RLOFprinting();
        BeBinaries BeB     = BeBinaries();
        // ALEJANDRO - initialized each star random seed to the binary random seed value, to keep track of errors in populations.
        currentBinary.star1.m_randomSeed    = currentBinary.m_randomSeed;
        currentBinary.star2.m_randomSeed    = currentBinary.m_randomSeed;

        if(options.populationDataPrinting){
            cout << "\nGenerating a new binary - " << i << endl;
            cout << "Binary has masses " << currentBinary.star1.m_Mass << " & " << currentBinary.star2.m_Mass << endl;
            cout << "Binary has initial separation " << currentBinary.m_SemiMajorAxisPrime << endl;
            cout << "randomSeed " << currentBinary.m_randomSeed << endl;
        }

        // Initialise time elapsed and timestep to 0
        double dt = 0.0;

        //---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
        // For now, just evolve the stars independently.
        //Evolve each star 1 and 2 timestep and compute using the mimimum time
        //bool supernova_star1=false;
        //bool supernova_star2=false;
        int j;
        double  t_MS1 = currentBinary.star1.gettMS();
        double  t_MS2 = currentBinary.star2.gettMS();
        double  dtPrev = 0;
        bool    debugging = false;
        //////////////////////////////////////////

        if(debugging){
            cout << "a (AU)" << TAB << TAB << currentBinary.m_SemiMajorAxisPrime << endl;
            cout << "e" << TAB << TAB << currentBinary.getEccentricity() << endl;
            cout << "w (yr-1)" << TAB << currentBinary.m_orbitalVelocityPrime << endl;
            cout << "w1 (yr-1)" << TAB << currentBinary.star1.m_omega << endl;
            cout << "w2 (yr-1)" << TAB << currentBinary.star2.m_omega << endl;
            cout << "wb1 (yr-1)" << TAB << currentBinary.star1.m_omegaBreak << endl;
            cout << "wb2 (yr-1)" << TAB << currentBinary.star2.m_omegaBreak << endl;
            cout << "wZ1 (yr-1)" << TAB << currentBinary.star1.m_omegaZAMS << endl;
            cout << "wZ2 (yr-1)" << TAB << currentBinary.star2.m_omegaZAMS << endl;
            cout << "m1 (Msol)" << TAB << currentBinary.star1.m_Mass << endl;
            cout << "m2 (Msol)" << TAB << currentBinary.star2.m_Mass << endl;
            cout << "R1 (Rsol)" << TAB << currentBinary.star1.m_Radius << endl;
            cout << "R2 (Rsol)" << TAB << currentBinary.star2.m_Radius << endl;
            cout << "Z1 " << TAB << currentBinary.star1.m_Metallicity << endl;
            cout << "Z2" << TAB << currentBinary.star2.m_Metallicity << endl;
            cout << "L" << TAB << currentBinary.m_TotalAngularMomentumPrime << endl;
            cout << "t_MS1 (Myr)" << TAB << t_MS1 << endl;
            cout << "t_MS2 (Myr)" << TAB << t_MS2 << endl<<endl;
        }

        if(options.detailedOutput){

            printDetailedOutputHeader(options.outputPath, i);
            //data.open((options.outputPath/dataOutputFile).string(),std::ios_base::app);

            printDetailedOutput(options.outputPath, i, currentBinary);
		}
        // Check if the secondary is massive enough to bother evolving (i.e. the system has some possibility of forming a double compact object)
        if(options.onlyDoubleCompactObjects){

            if(currentBinary.star2.m_Mass < MINIMUM_MASS_SECONDARY){
                binaryStatus = BINARY_NOT_OF_INTEREST;

                //specifically print formation channel here too because we want all channels
                //also the ones that stop so we can see why
                //printFormationChannel(options.outputPath, "formationHistory.txt", currentBinary);
            }
            else{
                binaryStatus = CONTINUE;
            }

        }

        // If the current binary is one of interest, evolve it
        if(binaryStatus == CONTINUE){

            // Set up strings to keep track of evolutionary history of a binary
            string  history = "";
            string  tempHistory = "";
            string  pastHistoryMT = "";
            string  tempHistoryMT = "";
            const   string  voidHistoryMT = "";


            // Evolve the current binary up to the maximum evolution time
            for( j = 0 ; currentBinary.m_time <= options.maxEvolutionTime ; j++){

                if(j>=options.maxNumberOfTimestepIterations){
                    std::cerr << currentBinary.m_randomSeed << "\tMax number timesteps reached: " << j << endl;
                    break;
                }
                //-------------------------------------------------------------------------------------------------------------------------------------------------//
                // Begins timestep calculation
                dt=min(currentBinary.star1.makeTimestep(options, r), currentBinary.star2.makeTimestep(options, r));

                if(debugging){
                    std::cout << "dt1, dt2, dt: " << currentBinary.star1.makeTimestep(options, r) << " " << currentBinary.star2.makeTimestep(options, r) << " " << dt << std::endl;
                }
// 				currentBinary.printingBinaryVariables();

                if(j==0){
                    dt=dt/1000;
                    dtPrev=dt;
                    currentBinary.m_timePrev = currentBinary.m_time;
                    currentBinary.m_time += dt; // maybe put this inside chooseTimestep function?
                }
                else{
//                        if(currentBinary.m_massTransferZeroFlag==false){
                    dt=currentBinary.chooseTimestep(options,dt,dtPrev);
                    currentBinary.m_timePrev = currentBinary.m_time;
                    currentBinary.m_time += dt; // maybe put this inside chooseTimestep function?
					currentBinary.m_dt = dt;
                }

                if(debugging){
                    std::cout << "After .chooseTimestep, dt: " << dt << std::endl;
                }
                dtPrev=dt;
                currentBinary.m_TotalAngularMomentumPrev=currentBinary.m_TotalAngularMomentumPrime; // Is this line ok here?
                // Ends timestep calculation
                //-------------------------------------------------------------------------------------------------------------------------//
                //-------------------------------------------------------------------------------------------------------------------------//
                // Begins SSE
                

                if(debugging){cout << "Before SSE"<<endl;}
                currentBinary.star1.evolveOneTimestep(dt, false, options, r);
                currentBinary.star2.evolveOneTimestep(dt, false, options, r);
                if(debugging){cout << "After SSE"<<endl;}
                // currentBinary.comparingLambdas(options);

                // history += chooseHistorySNe(currentBinary.star1.flagSN, currentBinary.star2.flagSN);    // In some cases, double print Sne (here and after binary)
                if(options.RLOFPrinting){
                    //change paramters if RLOF
                    RLOFp.setProps(currentBinary, false);
                }

                // Check if either star has experienced an error, and end the simulation if they have.
                if(currentBinary.star1.m_error){
                    if(debugging){std::cerr << currentBinary.m_randomSeed << "\tError flag in star1" << endl;}
                    break;
                }
                if(currentBinary.star2.m_error){
                    if(debugging){std::cerr << currentBinary.m_randomSeed << "\tError flag in star2" << endl;}
                    break;
                }

                // Check if a massless remnant has been formed
                if(currentBinary.star1.m_stellarType == MASSLESS_REMNANT or currentBinary.star2.m_stellarType == MASSLESS_REMNANT){break;}

                // Check if binary components are touching -- end simulation if this happens (should usually be avoided as MT or CE should happen prior to this)
                if(currentBinary.m_SemiMajorAxisPrime <= RsolToAU*(currentBinary.star1.m_Radius+currentBinary.star2.m_Radius) and currentBinary.m_SemiMajorAxisPrime > 0.0){
					if(debugging){
						std::cout << "Stars are touching" << std::endl;
					}
					break;
				}
                // Check if binary semi major axis is still positive (negative means unbound, presumably by SN)
				if(options.evolveUnboundSystems == false)
					if(currentBinary.m_SemiMajorAxisPrime <= 0.0 or currentBinary.m_Eccentricity > 1.0){break;}

                // Check if system remains gravitationally bound
                //  epsilon = -EPrime / E;
				if(options.evolveUnboundSystems == false){
					if(-(currentBinary.m_totalOrbitalEnergyPrime/currentBinary.m_totalOrbitalEnergyPrev) > 0.0){
						if(debugging){
							std::cout << "System unbound, DeltaE_{orb} = " << -(currentBinary.m_totalOrbitalEnergyPrime/currentBinary.m_totalOrbitalEnergyPrev) << std::endl;
						}
						break;
					}
				}
                //-------------------------------------------------------------------------------------------------------------------------//
                // Beging evaluate binary
                if(debugging){cout << "Begin evaluate binary" << endl;}
                currentBinary.evaluateBinary(options, r, dt, dtPrev);
                if(debugging){cout << "End evaluate binary" << endl << endl;}

                if(options.RLOFPrinting){
                    //change paramters if RLOF
                    RLOFp.setProps(currentBinary, true);
                    bool printRL =  RLOFp.printParameters( currentBinary);
                    if(printRL){
                        RLOFp.printRLOFParameters(options.outputPath, "RLOF.txt", currentBinary);
                    }
                }
                if(options.BeBinaries){
                    if((currentBinary.star1.m_stellarType == NEUTRON_STAR &
                       (currentBinary.star2.m_stellarType == MS_MORE_THAN_07 or
                        currentBinary.star2.m_stellarType == MS_LESS_THAN_07))
                       or
                        (currentBinary.star1.m_stellarType == NEUTRON_STAR &
                        (currentBinary.star2.m_stellarType == MS_MORE_THAN_07 or
                         currentBinary.star2.m_stellarType == MS_LESS_THAN_07))){
                        BeB.setProps(currentBinary);
                        BeB.printBeBinariesParameters(options.outputPath,BeBinariesFile);
                    }
                }


                // Check if binary evolution caused the stars in the binary to touch.
                if(currentBinary.m_SemiMajorAxisPrime <= RsolToAU*(currentBinary.star1.m_Radius+currentBinary.star2.m_Radius) and currentBinary.m_SemiMajorAxisPrime > 0.0){
                    if(debugging){
                        std::cout << "Stars are touching" << std::endl;
                    }
                    break;
                }

                // Check if binary semi major axis is still positive (negative means unbound, presumably by SN)
				if(options.evolveUnboundSystems == false)
					if(currentBinary.m_SemiMajorAxisPrime <= 0 or currentBinary.m_Eccentricity > 1.0){break;}

                // End evaluate binary
                //-------------------------------------------------------------------------------------------------------------------------//
                if(options.detailedOutput){
                    printDetailedOutput(options.outputPath, i, currentBinary);
                }

                //-------------------------------------------------------------------------------------------------------------------------//<< TAB << currentBinary.m_massTransferTrackerHistory

                if(currentBinary.m_stellarMerger == true){
                    if(debugging){cout << "Stellar merger. New semi-major axis:" << currentBinary.m_SemiMajorAxisPrime << endl;}
                    break;
                }

                if(currentBinary.m_errorFlag){
                    std::cerr << currentBinary.m_randomSeed << "\tError in Binary evolution." << endl;
                    break;
                }

                if(options.CHEvolution){
                    if(currentBinary.m_flagFirstTimeStep){
                        printCHEvolutionParameters(options.outputPath, "CHEvolutionParameters.txt", currentBinary);
                    }
                }
                if(options.evolvePulsars) {
                    if ((currentBinary.star1.m_stellarType == NEUTRON_STAR) or (currentBinary.star2.m_stellarType == NEUTRON_STAR)){
                    	printpulsarEvolutionParameters(options.outputPath, pulsarEvolutionFile, currentBinary);
                    }
                }
                // ALEJANDRO - 16/03/2017 - Check for double WD systems so we stop evolving them as they take time to process and we are currently not interested in dealing with them.
				if	((currentBinary.star1.m_stellarType >= HELIUM_WHITE_DWARF) and (currentBinary.star1.m_stellarType <= OXYGEN_NEON_WHITE_DWARF) and 
					(currentBinary.star2.m_stellarType >= HELIUM_WHITE_DWARF) and (currentBinary.star2.m_stellarType <= OXYGEN_NEON_WHITE_DWARF)){
					if(debugging){
						std::cout << currentBinary.m_randomSeed << "\tDouble White Dwarf system formed. Type1, Type2:\t" << currentBinary.star1.m_stellarType << ", " << currentBinary.star2.m_stellarType << std::endl;
					}
					break;
				}
				
				
			    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                // Check merging
                if((currentBinary.star1.m_stellarType == NEUTRON_STAR) and (currentBinary.star2.m_stellarType == NEUTRON_STAR)){
                    if(debugging){cout << "NSNS" << endl;}
                }

                if((currentBinary.star1.m_stellarType == NEUTRON_STAR or currentBinary.star1.m_stellarType == BLACK_HOLE) and (currentBinary.star2.m_stellarType == NEUTRON_STAR or currentBinary.star2.m_stellarType == BLACK_HOLE)){
                    bool    useWinds = false;
                    if(debugging){
                        cout << "S1, S2: " << currentBinary.star1.m_stellarType << "\t" << currentBinary.star2.m_stellarType << endl;
                    }
					
					// Check CBC
                    currentBinary.coalesce(options, doubleCompactObjects);

                    // Floor 10/11/2018 -- check if we have found a hit in exploratory phase of Adaptive Importance Sampling
                    // based on settings . And add  1 to the counter of DCOhits if we found a hit
                    if(options.AISexploratoryphase){ 
                        // track if we have a hit
                        int DCOhit = calculateDCOhit(options, currentBinary);
                        // DCOhit = 0 if no hit, DCOhit = 1 if we found a hit
                        aisvariables.CounterDCOsAIS = aisvariables.CounterDCOsAIS + DCOhit;
                    }
					break;
                }

                ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                if(currentBinary.m_errorFlag == true){break;}

                //After everything is evaluated in the first timestep set the first time step flag to false
                if(currentBinary.m_flagFirstTimeStep){currentBinary.m_flagFirstTimeStep = false;}

            }

			//This location is for printing parameters after evolving the entire binary 
            printFormationChannel(options.outputPath, formationHistoryFile, currentBinary);
			printSystemParameters(options.outputPath, systemParametersFile, currentBinary);

            
            // Printing for debbuging
            if(debugging){
                cout << endl << endl;
                cout << "a (AU)" << TAB << TAB << currentBinary.m_SemiMajorAxisPrime << endl;
                cout << "e" << TAB << TAB << currentBinary.getEccentricity() << endl;
                cout << "w (yr-1)" << TAB << currentBinary.m_orbitalVelocityPrime << endl;
                cout << "w1 (yr-1)" << TAB << currentBinary.star1.m_omega << endl;
                cout << "w2 (yr-1)" << TAB << currentBinary.star2.m_omega << endl;
                cout << "wZ1 (yr-1)" << TAB << currentBinary.star1.m_omegaZAMS << endl;
                cout << "wZ2 (yr-1)" << TAB << currentBinary.star2.m_omegaZAMS << endl;
                cout << "m1 (Msol)" << TAB << currentBinary.star1.m_Mass << endl;
                cout << "m2 (Msol)" << TAB << currentBinary.star2.m_Mass << endl;
                cout << "R1 (Rsol)" << TAB << currentBinary.star1.m_Radius << endl;
                cout << "R2 (Rsol)" << TAB << currentBinary.star2.m_Radius << endl;
                cout << "L1 ZAMS (Lsol)" << TAB << currentBinary.star1.m_LZAMS << endl;
                cout << "L2 ZAMS (Lsol)" << TAB << currentBinary.star1.m_LZAMS << endl;
                cout << "T1 ZAMS (Tsol)" << TAB << currentBinary.star1.m_TZAMS << endl;
                cout << "T2 ZAMS (Tsol)" << TAB << currentBinary.star1.m_TZAMS << endl;
                cout << "S1" << TAB << currentBinary.star1.m_stellarType << endl;
                cout << "S2" << TAB << currentBinary.star2.m_stellarType << endl;
                cout << "Met1 " << TAB << currentBinary.star1.m_Metallicity << endl;
                cout << "Met2" << TAB << currentBinary.star2.m_Metallicity << endl;
                cout << "t_MS1 (Myr)" << TAB << t_MS1 << endl;
                cout << "t_MS2 (Myr)" << TAB << t_MS2 << endl<<endl;
                cout << "time (Myr)" << TAB << currentBinary.star1.m_time << endl<<endl;

                cout << "Ei: " << TAB << currentBinary.getTotalEnergy() << TAB << "Li: " << TAB << currentBinary.getTotalAngularMomentum() << endl;
                cout << "Ef: " << TAB << currentBinary.m_TotalEnergyPrime << TAB << "Lf: " << TAB << currentBinary.m_TotalAngularMomentumPrime << endl;

                cout << "CE flag: " << TAB << currentBinary.m_commonEnvelopeFlag << endl;
                cout << endl << "number of timesteps: " << j << endl;
            }


        if(options.populationDataPrinting)
            cout << history << endl;
        }


        if(options.AISexploratoryphase){
            //  Floor 14/11/2018 - update stopping criteria expl phase. 
            if(aisvariables.CounterDCOsAIS >= 2){  // if we have found at least 2 hits in total (to overcome possibility of 100% of samples are hits)    

                //  update fexplAIS to estimate how long we should be spending on sampling from exploratory phase. 
                updateFexpl(options, aisvariables, i); // floor: we could update this only every other 10 runs in the future..
            } 

            aisvariables.fractionSampled = float(i) / options.nBinaries; // calculate fraction so far spend on expl phase
            if (aisvariables.fexplAIS!=1){  
                // if the fraction of total samples that we spend is larger than fexpl we should stop and switch to refinement phase
                if(aisvariables.fractionSampled >= aisvariables.fexplAIS){ 
                    cout << " We are stopping this simulation since we have spend a fraction fexpl on exploratory phase  \n " ;   // floor 
                    break;         
                }                  
            }       
        }
    }


    double duration;
    duration = ( clock() - clockStart ) / (double) CLOCKS_PER_SEC;
    if (!options.quiet) cout << "Time taken was " << duration << "s." << endl;


    // Close the outputs

    doubleCompactObjects.flush();
    doubleCompactObjects.close();
}
