#pathToData = './COMPASOutput.h5'
from __future__ import print_function
from __future__ import division

import sys

import h5py  as h5   #for handling data format
import numpy as np  #for array handling
import os           #For checking existence data
import WriteH5File

import h5py as h5
import sys

sys.path.append('/Users/floorbroekgaarden/Projects/BHNS_project/Scripts')
import ClassCosmicIntegrator  as CI #Given settings and redshifts returns rates (2D arrays) Loads the data
import coencodeVarious        as CV
from PostProcessingScripts import * 

def maskTargetDCOsSTROOPWAFEL(DCOtype, boolDCOmask, f): #, otherSelection, otherparam):
    """returns mask of DCOs of interest
    fxparam  is hdf5 keyname of file where variable for which you want to mask DCOs is in 
    DCOtype = 'BBH' / 'ALL' / 'BHNS' or 'BNS' 
    boolDCOmask = [Hubble, RLOF, Optimistic] # boolean values whether to mask mergers in a HUbble time, 
    binaries that have RLOFSecondaryAfterCEE = True, and Optimistic binaries (i.e. optimisticCEFlag == 0)
    pathToDirectory is pathname to Directory where _oratory & _sampling directories are
    """
    
    Hubble, RLOF, Optimistic = boolDCOmask
    

 
    
    fDCO = f['doubleCompactObjects']
    
        
    
    # mask binaries of given DCO type
    if (DCOtype == 'BNS') | (DCOtype=='NSNS'):
        mask0 = ((fDCO['stellarType1'][...] == 13) & (fDCO['stellarType2'][...] == 13))
    elif (DCOtype == 'BHNS') | (DCOtype == 'NSBH'):
        mask0 = ((fDCO['stellarType1'][...] == 13) & (fDCO['stellarType2'][...] == 14)) | \
            ((fDCO['stellarType1'][...] == 14) & (fDCO['stellarType2'][...] == 13) )          
    elif (DCOtype == 'BBH') | (DCOtype=='BHBH'):
        mask0 = ((fDCO['stellarType1'][...] == 14) & (fDCO['stellarType2'][...] == 14))
    elif (DCOtype == 'all') | (DCOtype == 'ALL') :
        mask0 = ((fDCO['stellarType1'][...] == 14) | (fDCO['stellarType1'][...] == 13))
    else:
        print('error: DCO type not known')
        
    # Hubble mask
    if Hubble:
        mask1 = (fDCO['mergesInHubbleTimeFlag'][...]==True) 
    elif not Hubble:
        mask1 = (fDCO['mergesInHubbleTimeFlag'][...]==True) |  (fDCO['mergesInHubbleTimeFlag'][...]==False) 

    # RLOF mask
    if RLOF:
        mask2 = np.ones_like(mask1)
    elif not RLOF:
        mask2 = np.ones_like(mask1)
    # Optimistic mask :  if True mask systems that have optimistic CE flag ==1
    if Optimistic:
        mask3 = np.logical_not(fDCO["optimisticCEFlag"][...] == 1) + \
        np.logical_not(fDCO["optimisticCEFlag"][...] == 0) 

    elif not Optimistic:
        mask3 = np.logical_not(fDCO["optimisticCEFlag"][...] == 1)  
    print('333')
    # combine the different masks and the oratory and refinement masks
    # tempcombinedmask = mask0*mask1*mask3
    combinedmask = mask0 * mask1 * mask2 * mask3
    combinedmask = combinedmask.squeeze()
    # if otherSelection =='UFD':
    #     KpcToKM = 3.086 * 10**(16) # kpc to km  
    #     MyrToYr = 1E6 # yrs
    #     YrToSec = 3.154 *1E7 #sec        
    #     UFD_epsilon = otherparam[0]
    #     UFD_Rvir = otherparam[1]
    #     Xbh1 = otherparam[2]
    #     Rns = otherparam[3]

    #     fSN = f['supernovae']
    #     seedsOfIntererst = fDCO['seed'][...].squeeze()
    #     seedsSN = fSN['randomSeed'][...].squeeze()
    #     bools = np.in1d(seedsSN, seedsOfIntererst)        
        
    #     tc  = fDCO['tc'][...].squeeze()
    #     vsys = fSN['systemicVelocity'][...].squeeze()[bools]
    #     vsysSN2 = vsys[1:][::2]
    #     traveldistance = tc * vsysSN2 *  MyrToYr * YrToSec
    #     radiusUFDgalaxy = UFD_epsilon * UFD_Rvir * KpcToKM
    #     maskCandidatesUFD = (traveldistance <= radiusUFDgalaxy) | ((vsysSN2 <= 44) & (tc * MyrToYr *YrToSec<= radiusUFDgalaxy)) 
        
    #     combinedmask = maskCandidatesUFD*combinedmask
    

    return combinedmask


def reduceH5file(pathToData, pathToDataWithoutCOMPASname, DCOtype, optimistic, BPSmodelName):
    """

    DCOtype in 'BBH', 'BNS', 'BHNS'
     """

    ff  = h5.File(pathToData) # + 'COMPASOutput.h5')
    # print(pathToData + 'COMPASOutput.h5')
    print("The main files I have at my disposal are:\n",list(ff.keys()))


    # print("The main files I have at my disposal are:\n",list(Data['formationChannels'].keys()))

	# Which Files do I want?
	# options: ['RLOF', 'XRayBinaries', 'commonEnvelopes', 'cppSource', 'doubleCompactObjects', 'formationChannels', 'pulsarEvolution', 'runtimes', 'supernovae', 'systems'])
    filesOfInterest   = {1:'doubleCompactObjects',2:'systems', 3:'supernovae', 4:'weights_detected',  5:'weights_intrinsic'}

	# #Give a list of columns you want, if you want all, say ['All']
	# columnsOfInterest = {1:['All'],\
	#                      2:['All'],\
	#                      3:['SEED', 'MZAMS_1', 'MZAMS_2']}
    columnsOfInterest =	   {1:['All'],\
							2:['All'],\
							3:['All'],\
                            4:['All'],\
                            5:['All']\
                                  }

	# #example of the seeds dictionary the actual one will be defined later
	# seedsOfInterest   = {1:None,\
	#                      2:None,\
	#                      3:None}
    
    # seedsDCO = ff['doubleCompactObjects']['seed'][()]
    seedsSystems = ff['systems']['SEED'][...].squeeze()
    seedsSN = ff['supernovae']['randomSeed'][...].squeeze()


    boolDCOmask = []
    boolDCOmask.append(1)
    boolDCOmask.append(1)
    boolDCOmask.append(optimistic)
    print('boolDCOmask', boolDCOmask)

    maskDCO = maskTargetDCOsSTROOPWAFEL(DCOtype=DCOtype, boolDCOmask=boolDCOmask, f=ff) #,\
                                            # otherSelection=None, otherparam=None)


    #what are the seeds of the systems that form BBHs
    seedsDCO = ff['doubleCompactObjects']['seed'][...].squeeze()[maskDCO]
    seedsSystemsMask = np.in1d(seedsSystems,seedsDCO)
    seedsSNMask = np.in1d(seedsSN, seedsDCO)
    seedsSystems = seedsSystems[seedsSystemsMask]
    seedsSN = seedsSN[seedsSNMask]
    seedsWeights = ff['weights_detected']['SEED'][...].squeeze()
    seedsWeightsMask = np.in1d(seedsWeights, seedsDCO)
    seedsWeights = seedsWeights[seedsWeightsMask]
    print('This should be true:')
    print(len(seedsWeights)==len(seedsDCO))
    print(len(seedsWeights),len(seedsDCO))

    # print((seedsDCO))
    # print((seedsSystems))
    # print((seedsSN))   
    # print(len(seedsDCO))
    # print(len(seedsSystems))
    # print(len(seedsSN))

    seedsOfInterest   = {1:seedsDCO,\
                          2:seedsSystems,\
                          3:seedsSN,\
                          4:seedsWeights,\
                          5:seedsWeights,\
                          }



    pathToNewData = pathToDataWithoutCOMPASname + 'COMPASCompactOutput_' + DCOtype + '_' + BPSmodelName +'temp.h5'


    WriteH5File.reduceH5(pathToOld = pathToData, pathToNew = pathToNewData,\
                     dictFiles=filesOfInterest, dictColumns=columnsOfInterest, dictSeeds=seedsOfInterest)






    ff.close() 




    print('-----------------------------')

    print()
    print('-----------------------------------------------')
    print('completed')
    print('SUCCESS!!!! YES!')




# for DCOtype in ['BHNS']:

#     BPSmodelName = 'B'
#     pathToData = '/Volumes/Andromeda/DATA/AllDCO/fiducial/' + 'COMPASCompactOutput_' + DCOtype + '_' + BPSmodelName +'.h5'
#     pathToDataWithoutCOMPASname = '/Volumes/Andromeda/DATA/AllDCO/fiducial/'
#     optimistic=int(1)
    
#     print()
#     print('-------------------------------------------------------')
#     print(' now at ', DCOtype, 'running for modelname', BPSmodelName)
#     print('which is in this directory', pathToDataWithoutCOMPASname)
#     print('optimistic = ', optimistic)

#     reduceH5file(pathToData=pathToData, pathToDataWithoutCOMPASname=pathToDataWithoutCOMPASname, DCOtype=DCOtype, optimistic=optimistic, BPSmodelName=BPSmodelName)



# for DCOtype in ['BHNS']:
#     BPSmodelName = 'A'
#     pathToData = '/Volumes/Andromeda/DATA/AllDCO/fiducial/' + 'COMPASCompactOutput_' + DCOtype + '_' + BPSmodelName +'.h5'
#     pathToDataWithoutCOMPASname = '/Volumes/Andromeda/DATA/AllDCO/fiducial/'
#     optimistic=int(0)
#     # BPSmodelName = 'A'
#     print()
#     print('-------------------------------------------------------')
#     print(' now at ', DCOtype, 'running for modelname', BPSmodelName)
#     print('which is in this directory', pathToDataWithoutCOMPASname)
#     print('optimistic = ', optimistic)

#     reduceH5file(pathToData=pathToData, pathToDataWithoutCOMPASname=pathToDataWithoutCOMPASname, DCOtype=DCOtype, optimistic=optimistic, BPSmodelName=BPSmodelName)



for DCOtype in ['BHNS']:
    BPSmodelName = 'G'
    pathToData = '/Volumes/Andromeda/DATA/AllDCO/zeroBHkick/' + 'COMPASCompactOutput_' + DCOtype + '_' + BPSmodelName +'.h5'
    pathToDataWithoutCOMPASname = '/Volumes/Andromeda/DATA/AllDCO/zeroBHkick/'
    optimistic=int(0)
    # BPSmodelName = 'G'
    print()
    print('-------------------------------------------------------')
    print(' now at ', DCOtype, 'running for modelname', BPSmodelName)
    print('which is in this directory', pathToDataWithoutCOMPASname)
    print('optimistic = ', optimistic)

    reduceH5file(pathToData=pathToData, pathToDataWithoutCOMPASname=pathToDataWithoutCOMPASname, DCOtype=DCOtype, optimistic=optimistic, BPSmodelName=BPSmodelName)


for DCOtype in ['BHNS']:
    BPSmodelName = 'L'
    pathToData = '/Volumes/Andromeda/DATA/AllDCO/ccSNkick_30km_s/' + 'COMPASCompactOutput_' + DCOtype + '_' + BPSmodelName +'.h5'
    pathToDataWithoutCOMPASname = '/Volumes/Andromeda/DATA/AllDCO/ccSNkick_30km_s/'
    optimistic=int(0)
    # BPSmodelName = 'L'
    print()
    print('-------------------------------------------------------')
    print(' now at ', DCOtype, 'running for modelname', BPSmodelName)
    print('which is in this directory', pathToDataWithoutCOMPASname)
    print('optimistic = ', optimistic)

    reduceH5file(pathToData=pathToData, pathToDataWithoutCOMPASname=pathToDataWithoutCOMPASname, DCOtype=DCOtype, optimistic=optimistic, BPSmodelName=BPSmodelName)


for DCOtype in ['BHNS']:
    BPSmodelName = 'K'
    pathToData = '/Volumes/Andromeda/DATA/AllDCO/ccSNkick_100km_s/' + 'COMPASCompactOutput_' + DCOtype + '_' + BPSmodelName +'.h5'
    # pathToData = '/Volumes/Andromeda/DATA/AllDCO/ccSNkick_100km_s/COMPASOutput.h5'
    pathToDataWithoutCOMPASname = '/Volumes/Andromeda/DATA/AllDCO/ccSNkick_100km_s/'
    optimistic=int(0)
    # BPSmodelName = 'K'
    print()
    print('-------------------------------------------------------')
    print(' now at ', DCOtype, 'running for modelname', BPSmodelName)
    print('which is in this directory', pathToDataWithoutCOMPASname)
    print('optimistic = ', optimistic)

    reduceH5file(pathToData=pathToData, pathToDataWithoutCOMPASname=pathToDataWithoutCOMPASname, DCOtype=DCOtype, optimistic=optimistic, BPSmodelName=BPSmodelName)



for DCOtype in ['BHNS']:
    BPSmodelName = 'D'
    pathToData = '/Volumes/Andromeda/DATA/AllDCO/alpha2_0/' + 'COMPASCompactOutput_' + DCOtype + '_' + BPSmodelName +'.h5'
    pathToDataWithoutCOMPASname = '/Volumes/Andromeda/DATA/AllDCO/alpha2_0/'
    optimistic=int(0)
    # BPSmodelName = 'D'
    print()
    print('-------------------------------------------------------')
    print(' now at ', DCOtype, 'running for modelname', BPSmodelName)
    print('which is in this directory', pathToDataWithoutCOMPASname)
    print('optimistic = ', optimistic)

    reduceH5file(pathToData=pathToData, pathToDataWithoutCOMPASname=pathToDataWithoutCOMPASname, DCOtype=DCOtype, optimistic=optimistic, BPSmodelName=BPSmodelName)


for DCOtype in ['BHNS']:
    BPSmodelName = 'C'
    pathToData = '/Volumes/Andromeda/DATA/AllDCO/alpha0_5/' + 'COMPASCompactOutput_' + DCOtype + '_' + BPSmodelName +'.h5'
    # pathToData = '/Volumes/Andromeda/DATA/AllDCO/alpha0_5/COMPASOutput.h5'
    pathToDataWithoutCOMPASname = '/Volumes/Andromeda/DATA/AllDCO/alpha0_5/'
    optimistic=int(0)
    # BPSmodelName = 'C'
    print()
    print('-------------------------------------------------------')
    print(' now at ', DCOtype, 'running for modelname', BPSmodelName)
    print('which is in this directory', pathToDataWithoutCOMPASname)
    print('optimistic = ', optimistic)

    reduceH5file(pathToData=pathToData, pathToDataWithoutCOMPASname=pathToDataWithoutCOMPASname, DCOtype=DCOtype, optimistic=optimistic, BPSmodelName=BPSmodelName)


# for DCOtype in ['BHNS']:
#     pathToData = '/Volumes/Virgo/DATA/AllDCO/rapid/COMPASOutput.h5'
#     pathToDataWithoutCOMPASname = '/Volumes/Virgo/DATA/AllDCO/rapid/'
#     optimistic=int(0)
#     BPSmodelName = 'F'
#     print()
#     print('-------------------------------------------------------')
#     print(' now at ', DCOtype, 'running for modelname', BPSmodelName)
#     print('which is in this directory', pathToDataWithoutCOMPASname)
#     print('optimistic = ', optimistic)

#     reduceH5file(pathToData=pathToData, pathToDataWithoutCOMPASname=pathToDataWithoutCOMPASname, DCOtype=DCOtype, optimistic=optimistic, BPSmodelName=BPSmodelName)



###########




# for DCOtype in ['BBH']:

#     pathToData = '/Volumes/Andromeda/DATA/AllDCO/fiducial/COMPASOutput.h5'
#     pathToDataWithoutCOMPASname = '/Volumes/Andromeda/DATA/AllDCO/fiducial/'
#     optimistic=int(1)
#     BPSmodelName = 'B'
#     print()
#     print('-------------------------------------------------------')
#     print(' now at ', DCOtype, 'running for modelname', BPSmodelName)
#     print('which is in this directory', pathToDataWithoutCOMPASname)
#     print('optimistic = ', optimistic)

#     reduceH5file(pathToData=pathToData, pathToDataWithoutCOMPASname=pathToDataWithoutCOMPASname, DCOtype=DCOtype, optimistic=optimistic, BPSmodelName=BPSmodelName)



# for DCOtype in ['BBH', 'BNS']:
#     pathToData = '/Volumes/Andromeda/DATA/AllDCO/fiducial/COMPASOutput.h5'
#     pathToDataWithoutCOMPASname = '/Volumes/Andromeda/DATA/AllDCO/fiducial/'
#     optimistic=int(0)
#     BPSmodelName = 'A'
#     print()
#     print('-------------------------------------------------------')
#     print(' now at ', DCOtype, 'running for modelname', BPSmodelName)
#     print('which is in this directory', pathToDataWithoutCOMPASname)
#     print('optimistic = ', optimistic)

#     reduceH5file(pathToData=pathToData, pathToDataWithoutCOMPASname=pathToDataWithoutCOMPASname, DCOtype=DCOtype, optimistic=optimistic, BPSmodelName=BPSmodelName)



# for DCOtype in ['BBH', 'BNS', 'BHNS']:
#     pathToData = '/Volumes/Andromeda/DATA/AllDCO/zeroBHkick/COMPASOutput.h5'
#     pathToDataWithoutCOMPASname = '/Volumes/Andromeda/DATA/AllDCO/zeroBHkick/'
#     optimistic=int(0)
#     BPSmodelName = 'G'
#     print()
#     print('-------------------------------------------------------')
#     print(' now at ', DCOtype, 'running for modelname', BPSmodelName)
#     print('which is in this directory', pathToDataWithoutCOMPASname)
#     print('optimistic = ', optimistic)

#     reduceH5file(pathToData=pathToData, pathToDataWithoutCOMPASname=pathToDataWithoutCOMPASname, DCOtype=DCOtype, optimistic=optimistic, BPSmodelName=BPSmodelName)


# for DCOtype in ['BBH', 'BNS', 'BHNS']:
#     pathToData = '/Volumes/Andromeda/DATA/AllDCO/ccSNkick_30km_s/COMPASOutput.h5'
#     pathToDataWithoutCOMPASname = '/Volumes/Andromeda/DATA/AllDCO/ccSNkick_30km_s/'
#     optimistic=int(0)
#     BPSmodelName = 'L'
#     print()
#     print('-------------------------------------------------------')
#     print(' now at ', DCOtype, 'running for modelname', BPSmodelName)
#     print('which is in this directory', pathToDataWithoutCOMPASname)
#     print('optimistic = ', optimistic)

#     reduceH5file(pathToData=pathToData, pathToDataWithoutCOMPASname=pathToDataWithoutCOMPASname, DCOtype=DCOtype, optimistic=optimistic, BPSmodelName=BPSmodelName)


# for DCOtype in ['BBH', 'BNS', 'BHNS']:
#     pathToData = '/Volumes/Andromeda/DATA/AllDCO/ccSNkick_100km_s/COMPASOutput.h5'
#     pathToDataWithoutCOMPASname = '/Volumes/Andromeda/DATA/AllDCO/ccSNkick_100km_s/'
#     optimistic=int(0)
#     BPSmodelName = 'K'
#     print()
#     print('-------------------------------------------------------')
#     print(' now at ', DCOtype, 'running for modelname', BPSmodelName)
#     print('which is in this directory', pathToDataWithoutCOMPASname)
#     print('optimistic = ', optimistic)

#     reduceH5file(pathToData=pathToData, pathToDataWithoutCOMPASname=pathToDataWithoutCOMPASname, DCOtype=DCOtype, optimistic=optimistic, BPSmodelName=BPSmodelName)



# for DCOtype in ['BBH', 'BNS', 'BHNS']:
#     pathToData = '/Volumes/Andromeda/DATA/AllDCO/alpha2_0/COMPASOutput.h5'
#     pathToDataWithoutCOMPASname = '/Volumes/Andromeda/DATA/AllDCO/alpha2_0/'
#     optimistic=int(0)
#     BPSmodelName = 'D'
#     print()
#     print('-------------------------------------------------------')
#     print(' now at ', DCOtype, 'running for modelname', BPSmodelName)
#     print('which is in this directory', pathToDataWithoutCOMPASname)
#     print('optimistic = ', optimistic)

#     reduceH5file(pathToData=pathToData, pathToDataWithoutCOMPASname=pathToDataWithoutCOMPASname, DCOtype=DCOtype, optimistic=optimistic, BPSmodelName=BPSmodelName)


# for DCOtype in ['BBH', 'BNS', 'BHNS']:
#     pathToData = '/Volumes/Andromeda/DATA/AllDCO/alpha0_5/COMPASOutput.h5'
#     pathToDataWithoutCOMPASname = '/Volumes/Andromeda/DATA/AllDCO/alpha0_5/'
#     optimistic=int(0)
#     BPSmodelName = 'C'
#     print()
#     print('-------------------------------------------------------')
#     print(' now at ', DCOtype, 'running for modelname', BPSmodelName)
#     print('which is in this directory', pathToDataWithoutCOMPASname)
#     print('optimistic = ', optimistic)

#     reduceH5file(pathToData=pathToData, pathToDataWithoutCOMPASname=pathToDataWithoutCOMPASname, DCOtype=DCOtype, optimistic=optimistic, BPSmodelName=BPSmodelName)














# if __name__ == "__main__":
#     pathToData = (sys.argv[1])
#     pathToDataWithoutCOMPASname  = (sys.argv[2])
#     DCOtype = (sys.argv[3])
#     optimistic = int(sys.argv[4])
#     BPSmodelName = (sys.argv[5])

#     print(optimistic)
    
# #    print('test')
#     reduceH5file(pathToData,  pathToDataWithoutCOMPASname,  DCOtype, optimistic, BPSmodelName)


### EXAMPLE:
#python reduceHdf5DataMinimal.py '/Volumes/Virgo/DATA/BHNS/alpha_10/COMPASOutput.h5' '/Volumes/Virgo/DATA/BHNS/alpha_10/'  'BNS' 0 'A'
# python rewriteToCompactFileWithCIweights.py '/Volumes/Andromeda/DATA/AllDCO/fiducial/COMPASOutput.h5' '/Volumes/Andromeda/DATA/AllDCO/fiducial/'  'BNS' 1 'B'
