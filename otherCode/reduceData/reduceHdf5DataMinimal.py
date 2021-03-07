#pathToData = './COMPASOutput.h5'


import sys

import h5py  as h5   #for handling data format
import numpy as np  #for array handling
import os           #For checking existence data
import WriteH5File



def reduceH5file(pathToData, pathToNewData):

    # read data
    Data  = h5.File(pathToData)
    print("The main files I have at my disposal are:\n",list(Data.keys()))


    # print("The main files I have at my disposal are:\n",list(Data['formationChannels'].keys()))

	# Which Files do I want?
	# options: ['RLOF', 'XRayBinaries', 'commonEnvelopes', 'cppSource', 'doubleCompactObjects', 'formationChannels', 'pulsarEvolution', 'runtimes', 'supernovae', 'systems'])
    filesOfInterest   = {1:'doubleCompactObjects',2:'systems'}

	# #Give a list of columns you want, if you want all, say ['All']
	# columnsOfInterest = {1:['All'],\
	#                      2:['All'],\
	#                      3:['SEED', 'MZAMS_1', 'MZAMS_2']}
    columnsOfInterest =	   {1:[ 'ID', 'M1', 'M1ZAMS', 'M2', 'M2ZAMS', 'Metallicity1', 'Metallicity2',  \
							 	 'mergesInHubbleTimeFlag', 'optimisticCEFlag',   'seed', \
							    'separationDCOFormation', 'separationInitial',  'stellarType1', 'stellarType2', 'tc', 'tform',  'weight'],\
							2:['ID', 'Metallicity1', 'Metallicity2', 'SEED',  'weight'],\
							}

	# #example of the seeds dictionary the actual one will be defined later
	# seedsOfInterest   = {1:None,\
	#                      2:None,\
	#                      3:None}
    
    seedsDCO = Data['doubleCompactObjects']['seed'][()]
    seedsSystems = Data['systems']['SEED'][()]


    seedsOfInterest   = {1:seedsDCO,\
                          2:seedsSystems\
                          }

    WriteH5File.reduceH5(pathToOld = pathToData, pathToNew = pathToNewData,\
                     dictFiles=filesOfInterest, dictColumns=columnsOfInterest, dictSeeds=seedsOfInterest)



if __name__ == "__main__":
    pathToData = (sys.argv[1])
    pathToNewData = (sys.argv[2])
    
#    print('test')
    reduceH5file(pathToData, pathToNewData)




