import numpy as np
import os
import time
from subprocess import Popen, PIPE
import sys
import pickle
import compas_hpc_functions as chpc





# def jobFileName(pythonNameJob, dirname):
# 	if pythonNameJob=='calculateRate':

# 		JobPythonName = 'test.py'
# 		pythonNameJob = JobPythonName #os.path.join(dirname, JobPythonName)

# 	elif pythonNameJob=='B':
# 		pythonNameJob =='B'


# 	return pythonNameJob



def runAndWrite(pythonNameJob, pathToData,  dirname,  DCOtype, optimistic, BPSmodelName, rootdir):
	"""
	pythonNameJob = pythonScript name should end with .py 
	dirname, 
	BPSmodelName, 
	rootdir

	 """
	pathToDataWithoutCOMPASname = dirname +'/'
	print('running', pythonNameJob)
	print('for dirname', pathToDataWithoutCOMPASname)
	print('model name = ', BPSmodelName)
	
	

	number_of_nodes		='1'
	number_of_cores 	= '1'	
	JobWalltime			= "100:00:00"
	JobMemory			= "50000"
	send_email = True # Whether you want to recieve an email when your jobs are done
	user_email = 'floor.broekgaarden@cfa.harvard.edu'




	# copy python script that I want to run
	bashCommand = 'cp ' + os.path.join(rootdir, pythonNameJob) + " " + os.path.join(dirname, '.')
	chpc.runBashCommand(bashCommand, verbose=True)


	# what will be the command to run 
	JobCommand  = 'python ' + pythonNameJob   + " " + pathToData + " " +  pathToDataWithoutCOMPASname + " " +   DCOtype + " " +  'Optimistic' + str(optimistic) + " " +  BPSmodelName #+   " " + input1 + " " + input2


	os.chdir(dirname)
	JobOutputFile = os.path.join(dirname, pythonNameJob +'_' +DCOtype + '_' + BPSmodelName  +  '.out')
	JobErrorFile  = os.path.join(dirname, pythonNameJob +'_' +DCOtype + '_' + BPSmodelName  + '.err')
	JobLogFile    = os.path.join(dirname, pythonNameJob +'_' +DCOtype + '_' + BPSmodelName  + '.log')

	JobCommand += " > " + JobLogFile


	JobString= chpc.GenerateSlurmOdysseyJobString(pythonNameJob, number_of_nodes, number_of_cores, JobOutputFile, JobErrorFile, JobWalltime, JobMemory, send_email, user_email=user_email, run_directory=dirname, command=JobCommand)

	print(JobString)

	# sbatchAISCommand = 'sbatch --dependency=afterok:' + str(AISUpdateWeightsJobID)

	# print(sbatchAISCommand)

	# Added to print bash script
	# Save AIS - Combine job to a file
	sbatchJobFile = open(dirname+'/' + pythonNameJob +'.bash' ,'w')
	sbatchJobFile.write(JobString)
	sbatchJobFile.close()

	bashCommand = 'sbatch ' + os.path.join(dirname, pythonNameJob+'.bash') 
	chpc.runBashCommand(bashCommand, verbose=True)

if __name__ == "__main__":
    # z = float((sys.argv[1]))
    # pythonNameJob = (sys.argv[1])
    # dirname = (sys.argv[2])
    # modelname = (sys.argv[3])

    rootdir = '/n/holystore01/LABS/berger_lab/Lab/fbroekgaarden/DATA/all_dco_legacy/PostProcessing'


    pathToData = (sys.argv[1])
    dirname  = (sys.argv[2])
    DCOtype = (sys.argv[3])
    optimistic = int(sys.argv[4])
    BPSmodelName = (sys.argv[5])

    pythonNameJob = 'rewriteToCompactFileWithCIweights.py'

#    print('test')
    runAndWrite(pythonNameJob, pathToData,  dirname,  DCOtype, optimistic, BPSmodelName, rootdir)


# example:
#runAndWritePP.py   'test.py' '/n/holystore01/LABS/berger_lab/Lab/fbroekgaarden/DATA/all_dco_legacy/unstableCaseBB' 'E' 



# pythonNameJob, dirname, modelname, rootdir
#python runAndWritePP.py '/Volumes/Virgo/DATA/BHNS/alpha_10/COMPASOutput.h5' '/Volumes/Virgo/DATA/BHNS/alpha_10'  'BNS' 0 'A'
# python runAndWritePP.py '/n/holystore01/LABS/berger_lab/Lab/fbroekgaarden/DATA/all_dco_legacy/fiducial/COMPASOutput.h5' '/n/holystore01/LABS/berger_lab/Lab/fbroekgaarden/DATA/all_dco_legacy/fiducial/'  'BNS' 0 'A'


