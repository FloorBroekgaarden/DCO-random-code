# from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import h5py as h5
import time
import sys
import copy
#Quick fudge to make import from ../Scripts work
import sys
sys.path.append('/Users/floorbroekgaarden/Projects/BHNS_project/Scripts')
import string

import ClassCosmicIntegrator  as CI #Given settings and redshifts returns rates (2D arrays) Loads the data
import coencodeVarious        as CV
from PostProcessingScripts import * 



from IPython.core.display import display, HTML
display(HTML("<style>.container { width:100% !important; }</style>"))


import pandas as pd




####################################
#path to the data


def writeToRatesFile(modelname='Z', pathCOMPASOutput='home', DCOtype='BHNS', Optmistic=False ):


	print('writing and calculating rate for model ', modelname)
	#Will only look at BBHs so might as well set everything
	minz = 0.
	maxz = 2.
	resz = 100
	Data = CI.CosmicIntegrator(COMPASpath = pathCOMPASOutput, DCOtypes=DCOtype,\
	       minRedshift=minz,   maxRedshift=maxz, nrRedshiftBins=resz, optimistic=Optmistic, Cosmology='WMAP')

	#I use the custom cosmology because this was the flatlambda prescription used before WMAP Stevenson et al 2019
	#Doesnt matter to much (between WMAP and 
	#this it is 22, and 22.7 per year) but to prevent redoing all the numbers in the tex for referee

	print(Data.COMPAS.mass1)
	print(len(Data.COMPAS.mass1))


	#######




	rates    = []
	totals   = []
	ratesIntrinsic_z0 = []
	labelslist = []

	GSMFs = ['Panter et al. (2004) Single', 'Furlong et al. (2015) Single', 'Furlong et al. (2015) Double']
	MZs   = [ 'Langer et al. (2006)'      , 'Langer et al. +offset (2006)', 'Ma et al. (2015)']
	SFRs  = ['Madau et al. (2014)'         ,'Strolger et al. (2004)',     'Madau et al. (2017)']


	# Neijssel:

	Data.MSSFR.Zprescription         = 'logNormal' 
	Data.MSSFR.SFRprescription       = 'Neijssel et al. (2019)'
	Data.MSSFR.logNormalPrescription = 'Neijssel Phenomenological'
	Data.MSSFR.GSMFprescription      = None
	Data.MSSFR.ZMprescription        = None
	Data.cosmologicalIntegration()
	weightSTROOPWAFEL = Data.COMPAS.weight # //floor weight
	Row        =  np.sum(Data.PerSystemPerRedshift_ratesObserved*weightSTROOPWAFEL, axis=0) # //floor weight
	#     Row        =np.sum(Data.PerSystemPerRedshift_ratesObserved, axis=0)
	ratesPerSystem_z0 = Data.PerSystemPerRedshift_ratesIntrinsic[0,:] * Data.COMPAS.weight
	ratesIntrinsic_z0.append(np.sum(ratesPerSystem_z0))
	print(ratesIntrinsic_z0)
	rates.append(Row)
	totals.extend([np.sum(Row)])
	labelslist.append('Neijssel 2019')
	iii=1
	for ind_GSMF, GSMF in enumerate(GSMFs):
	         for ind_MZ, MZ in enumerate(MZs):
	             for ind_SFR, SFR in enumerate(SFRs):
	                iii+=1
	                print('now at ', iii)
	                Data.MSSFR.Zprescription         = 'MZ_GSMF'
	                Data.MSSFR.SFRprescription       = SFR
	#                 Data.MSSFR.logNormalPrescription = logNormal[nrL]
	                Data.MSSFR.GSMFprescription      = GSMF
	                Data.MSSFR.ZMprescription        = MZ
	                Data.cosmologicalIntegration()
	                weightSTROOPWAFEL = Data.COMPAS.weight # //floor weight
	                Row        =  np.sum(Data.PerSystemPerRedshift_ratesObserved*weightSTROOPWAFEL, axis=0) # //floor weight
	                ratesPerSystem_z0 = Data.PerSystemPerRedshift_ratesIntrinsic[0,:] * Data.COMPAS.weight

	                ratesIntrinsic_z0.append(np.sum(ratesPerSystem_z0))
	                rates.append(Row)
	                totals.extend([np.sum(Row)])

	                labelslist.append(SFR + r'$\_+\_' + GSMF + r'$\_+\_' + MZ )





	if DCOtype=='BBH':
		DCOname = 'BHBH'
	elif DCOtype=='BNS':
		DCOname = 'NSNS'
	elif DCOtype=='BHNS':
		DCOname = 'BHNS'

       

	df = pd.read_csv('rates_MSSFR_Models_'+DCOname+'.csv', index_col=0)
	namez0 = modelname + ' intrinsic (z=0) [Gpc^{-3} yr^{-1}]'
	nameObs = modelname + ' observed (design LVK) [yr^{-1}]'
	# rates0 =  df[name0]

	df[namez0] = ratesIntrinsic_z0
	df[nameObs] = totals 

	df.to_csv('rates_MSSFR_Models_'+DCOtype+'.csv')
	return






# Models to RUN 

# May 20: I am updating my data with the AllDCO focused runs :-) 

# this is an overwrite with better data (old ones are in BHNS copy)


# for DCOtype in ['BHNS', 'BBH', 'BNS']:
# 	print('at DCOtype =', DCOtype)
# 	pathCOMPASOutput = '/Volumes/Andromeda/DATA/AllDCO/fiducial/'
# 	modelname = 'A'
# 	writeToRatesFile(modelname=modelname, pathCOMPASOutput=pathCOMPASOutput, DCOtype=DCOtype, Optmistic=False)


# 	pathCOMPASOutput = '/Volumes/Andromeda/DATA/AllDCO/fiducial/'
# 	modelname = 'B'
# 	writeToRatesFile(modelname=modelname, pathCOMPASOutput=pathCOMPASOutput, DCOtype=DCOtype, Optmistic=True)


# 	pathCOMPASOutput = '/Volumes/Andromeda/DATA/AllDCO/zeroBHkick/'
# 	modelname = 'G'
# 	writeToRatesFile(modelname=modelname, pathCOMPASOutput=pathCOMPASOutput, DCOtype=DCOtype, Optmistic=True)


for DCOtype in ['BBH', 'BNS']:

	pathCOMPASOutput = '/Volumes/Andromeda/DATA/AllDCO/zeroBHkick/'
	modelname = 'G'
	print('modelname')
	writeToRatesFile(modelname=modelname, pathCOMPASOutput=pathCOMPASOutput, DCOtype=DCOtype, Optmistic=True)

	print('at DCOtype =', DCOtype)
	pathCOMPASOutput = '/Volumes/Andromeda/DATA/AllDCO/fiducial/'
	modelname = 'A'
	print('modelname')
	writeToRatesFile(modelname=modelname, pathCOMPASOutput=pathCOMPASOutput, DCOtype=DCOtype, Optmistic=False)


	pathCOMPASOutput = '/Volumes/Andromeda/DATA/AllDCO/fiducial/'
	modelname = 'B'
	print('modelname')
	writeToRatesFile(modelname=modelname, pathCOMPASOutput=pathCOMPASOutput, DCOtype=DCOtype, Optmistic=True)





# pathCOMPASOutput = '/Volumes/Andromeda/DATA/AllDCO/alpha0_5/'
# modelname, DCOtype = 'M', 'BNS'
# writeToRatesFile(modelname=modelname, pathCOMPASOutput=pathCOMPASOutput, DCOtype=DCOtype, Optmistic=False)

# DCOtype='BHNS'
# writeToRatesFile(modelname=modelname, pathCOMPASOutput=pathCOMPASOutput, DCOtype=DCOtype, Optmistic=False)

# DCOtype='BBH'
# writeToRatesFile(modelname=modelname, pathCOMPASOutput=pathCOMPASOutput, DCOtype=DCOtype, Optmistic=False)



# pathCOMPASOutput = '/Volumes/Andromeda/DATA/AllDCO/alpha2_0/'
# modelname, DCOtype = 'N', 'BNS'
# writeToRatesFile(modelname=modelname, pathCOMPASOutput=pathCOMPASOutput, DCOtype=DCOtype, Optmistic=False)

# DCOtype='BHNS'
# writeToRatesFile(modelname=modelname, pathCOMPASOutput=pathCOMPASOutput, DCOtype=DCOtype, Optmistic=False)

# DCOtype='BBH'
# writeToRatesFile(modelname=modelname, pathCOMPASOutput=pathCOMPASOutput, DCOtype=DCOtype, Optmistic=False)



