import numpy as np
import os
import gc
import time
import subprocess as sp
import h5py as h5

###Options####
collectAndCombineSystems = True
collectAndCombineDoubleCompactObjects = True
collectAndCombineErrors = True
collectAndCombineFormation = True
collectAndCombineXRayBinaries = True
collectAndCombineCommonEnvelopes = True
collectAndCombineSupernovae = True
collectAndCombineRuntimes = True
collectAndCombineRLOF = True
collectAndCombinepulsarEvolution = True
runPlottingRoutines = False #See plottingRoutines.py for options
cleanUp = False

#-- Base directory where results are
baseDirectory = os.getcwd()

#-- Where to output post-processed results
rootOutputDir = os.getcwd()

# If weights files exist, add them to the HDF5 file
if os.path.isfile(os.path.join(baseDirectory, "allDoubleCompactObjectsWeights.txt")):
    addHdfGroupDCOWeights = True
else: 
    addHdfGroupDCOWeights = False
if os.path.isfile(os.path.join(baseDirectory, "allSystemsWeights.txt")):
    addHdfGroupAllWeights = True
else: 
    addHdfGroupAllWeights = False

# only if rejected samples file exist in output0 (means in all outputs):
if os.path.isfile(os.path.join(baseDirectory+'/output0', "rejectedSamples.txt")):
    collectAndCombineRejectedSamples = True
else: 
    collectAndCombineRejectedSamples = False

##############

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--masterFolderDir", help="Directory to masterFolder")
args = parser.parse_args()

if(runPlottingRoutines):
    from plottingRoutines import plottingRoutines

compasRootDir = os.environ.get('COMPAS_ROOT_DIR')

sourceDirectory = os.path.join(compasRootDir, 'COMPAS')

#make the HDF5 file
h5File = h5.File(os.path.join(rootOutputDir, 'COMPASOutput.h5'),'w')





def copyPythonSubmit(baseDirectory="./"):
    """
    """
    #put the contents of the python submit in there
    ps = open(os.path.join(baseDirectory, "pythonSubmit.py")) #os.path.join(baseDirectory, 'masterFolder/pythonSubmit.py'))
    #makes one long string and stores it, with \n linebreaks included
    pythonSubmit = ''
    for l in ps:
        pythonSubmit += l
    ps.close()

    #psDset = h5File.create_dataset('pythonSubmit', dtype="S" + str(len(pythonSubmit)), shape=(1,))
    #psDset[0] = pythonSubmit 
    h5File.attrs.create('pythonSubmit', pythonSubmit, dtype=h5.special_dtype(vlen=str))

    return

def collectSource():
    """
    """
    #put the c++ source in there
    sourceGroup = h5File.create_group('cppSource')

    for root,dirs,files in os.walk(sourceDirectory):
        for f in files:
            if f[-2:] == '.h' or f[-4:] == '.cpp' or f[-4:] == '.hpp':

                #makes one long string and stores it, with \n linebreaks included
                src = ''
                path = os.path.join(root,f)
                fil = open(path)
                for l in fil:
                    src += l
                fil.close()
                #srcDset = sourceGroup.create_dataset(f,dtype="S"+str(len(src)),shape=(1,))
                #srcDset[0] = src
                sourceGroup.attrs.create(f, src, dtype=h5.special_dtype(vlen=str))
    return

def collectAndCombine(collected_file, individual_output_file, nHeaders=3, baseDirectory='.'):
    """
    There was lots of duplicated code in this python function for collecting and combining different text files. 
    This function aims to generalise that code (originally written by Jim) and make it reusable.

    Parameters
    -----------
    collected_file : string
        Filename of the collected output (e.g. allXRayBinaries.txt etc)
    individual_output_file : string
        Filename of the individual output (e.g. XRayBinaries.txt)
    nHeaders : int
        Number of lines of headers in individual output. These will be copied to the collected output only once.
    baseDirectory : string
        The base directory which contains all of the files to be collected

    Returns
    --------
    """
    allOutputsFile = open(os.path.join(rootOutputDir, collected_file), 'w')
    headersWritten = False

    if(nHeaders < 0):
        raise ValueError("nHeaders must be >= 0, you gave " + str(nHeaders))

    if(np.logical_not(type(nHeaders) == int)):
        raise TypeError("nHeaders should be an int. Try int(nHeaders)")

    for root,dirs,files in os.walk(baseDirectory):

        for f in files:
            
            if f == individual_output_file:

                currentLine = 0

                path = os.path.join(root,f)
                
                outputFile = open(path)

                if not headersWritten:
                    for i in range(nHeaders):
                        line = outputFile.readline()
                        nCols = len(line.split('\t'))
                        allOutputsFile.write(line)
                        headersWritten = True
                else:
                    [outputFile.readline() for i in range(nHeaders)]

                for line in outputFile:

                    thisNCols = len(line.split('\t'))
                    
                    if(nHeaders > 0):
                        if(thisNCols == nCols):
                            allOutputsFile.write(line)
                            #print("Line OK!")
                        else:
                            print("Line incorrect length in ", path)
                            print("length ", thisNCols)
                    else:
                        allOutputsFile.write(line)
                                   
                outputFile.close()

    allOutputsFile.close()

def addHdf5Group(groupName, asciiFileName):

    #too slow to go line by line, so load in a modest (in term sof memory) amount at a time
    chunkSize = 500000

    group = h5File.create_group(groupName)

    f = open(os.path.join(rootOutputDir, asciiFileName))

    units = f.readline()[:-1].split()
    dtypeStrs = f.readline()[:-1].split()
    headers = f.readline()[:-1].split()

    nColumns = len(headers)

    dtypes = []
    for dts in dtypeStrs:
        if dts == 'int':
            dtypes.append(np.int64)
        elif dts == 'float':
            dtypes.append(np.float64)
        elif dts == 'bool':
            dtypes.append(bool)
        else:
            raise ValueError("Unrecognised datatype : " + dts)


    #get the length of the file (minus headers)
    fileLength = int(sp.check_output('wc -l ' + os.path.join(rootOutputDir, asciiFileName), shell=True).split()[0]) - 3

    #create the dataset for each column (and add the units as an attribute)
    for header,dtype,unit in zip(headers,dtypes,units):

        dset = group.create_dataset(header,dtype=dtype,shape=(fileLength,1))
        dset.attrs['units'] = unit


    chunkBegin = 0
    chunkEnd = 0
    while chunkEnd < fileLength:

        data = []

        chunkEnd = chunkBegin + chunkSize

        #dont try to load in more data than you've got
        if chunkEnd > fileLength:
            chunkEnd = fileLength

        #read in a modest number of lines
        for i in range(chunkEnd-chunkBegin):

            data.append(f.readline()[:-1].split())

        data = np.array(data)

        #cast the columns to the correct type en masse
        for en, d in enumerate(data.T):

            (group[headers[en]])[chunkBegin:chunkEnd] = np.array(d,dtype=dtypes[en]).reshape((-1,1))

        chunkBegin = chunkEnd

def cleanUpOutputFolders():
    """
    """
    os.chdir(rootOutputDir)
    
    tarCommand = "tar czf tarredFolders.tar.gz output*"
    
    print(tarCommand)
    
    os.system(tarCommand)

    rmCommand = "rm -r ./output*"
    
    print(rmCommand)

    os.system(rmCommand)

    return

copyPythonSubmit(baseDirectory=args.masterFolderDir)
collectSource()



if collectAndCombineSystems:
    collectAndCombine('allSystems.dat','systemParameters.txt', baseDirectory=baseDirectory)
    addHdf5Group('systems', 'allSystems.dat')

if collectAndCombineDoubleCompactObjects:
    collectAndCombine('allDoubleCompactObjects.dat', 'doubleCompactObjects.txt', baseDirectory=baseDirectory)
    addHdf5Group('doubleCompactObjects', 'allDoubleCompactObjects.dat')

if(collectAndCombineSupernovae):
    collectAndCombine('allSupernovae.txt', 'supernovae.txt', baseDirectory=baseDirectory)
    addHdf5Group('supernovae', 'allSupernovae.txt')

if(collectAndCombineRuntimes):
    collectAndCombine('allRuntimes.txt', 'runtimes.txt', baseDirectory=baseDirectory)
    addHdf5Group('runtimes', 'allRuntimes.txt')

if(collectAndCombineXRayBinaries):
    collectAndCombine('allXRayBinaries.txt', 'XRayBinaries.txt', baseDirectory=baseDirectory)
    addHdf5Group('XRayBinaries', 'allXRayBinaries.txt')

if(collectAndCombineCommonEnvelopes):
    collectAndCombine('allCommonEnvelopes.txt', 'commonEnvelopes.txt', baseDirectory=baseDirectory)
    addHdf5Group('commonEnvelopes', 'allCommonEnvelopes.txt')

if(collectAndCombinepulsarEvolution):
    collectAndCombine('allpulsarEvolution.txt', 'pulsarEvolution.txt', baseDirectory=baseDirectory)
    addHdf5Group('pulsarEvolution', 'allpulsarEvolution.txt')

if(collectAndCombineFormation):
    collectAndCombine('allFormation.txt', 'formationHistory.txt', baseDirectory=baseDirectory)
    addHdf5Group('formationChannels', 'allFormation.txt')


if(addHdfGroupDCOWeights):
    addHdf5Group('doubleCompactObjectsWeights', 'allDoubleCompactObjectsWeights.txt')
if(addHdfGroupAllWeights):
    addHdf5Group('systemsWeights', 'allSystemsWeights.txt')

if(collectAndCombineErrors):
    collectAndCombine('allErrors.txt', 'error.txt', nHeaders=0, baseDirectory=baseDirectory)

if(collectAndCombineRLOF):
    collectAndCombine('allRLOF.txt', 'RLOF.txt', baseDirectory=baseDirectory)
    addHdf5Group('RLOF', 'allRLOF.txt')

if(collectAndCombineRejectedSamples):
    collectAndCombine('allRejectedSamples.txt', 'rejectedSamples.txt', baseDirectory=baseDirectory)
    addHdf5Group('rejectedSamples', 'allRejectedSamples.txt')   

if(runPlottingRoutines):
    plottingRoutines()

if(cleanUp):
    cleanUpOutputFolders()

print("post-processing completed successfully")
