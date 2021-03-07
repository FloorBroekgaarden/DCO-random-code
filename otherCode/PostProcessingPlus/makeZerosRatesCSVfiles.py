# CREATE PANDAS FILE 
 
import h5py as h5
import numpy as np
import string
import pandas as pd


nModels=26
BPSnameslist = list(string.ascii_uppercase)[0:nModels]


GSMFs = ['Panter et al. (2004) Single', 'Furlong et al. (2015) Single', 'Furlong et al. (2015) Double']
MZs   = [ 'Langer et al. (2006)'      , 'Langer et al. +offset (2006)', 'Ma et al. (2015)']
SFRs  = ['Madau et al. (2014)'         ,'Strolger et al. (2004)',     'Madau et al. (2017)']




MSSFRnameslist = []
MSSFRnameslist.append('000') # add phenomenological 

for ind_GSMF, GSMF in enumerate(GSMFs):
    ind_y = ind_GSMF + 1
    for ind_MZ, MZ in enumerate(MZs):
        ind_z = ind_MZ +1
        for ind_SFR, SFR in enumerate(SFRs):
            ind_x = ind_SFR+1


NAMES = []

for ind_l, L in enumerate(BPSnameslist):
    str_z0 = str(L + ' intrinsic (z=0) [Gpc^{-3} yr^{-1}]')
    str_obs = str(L + ' observed (design LVK) [yr^{-1}]')
    NAMES.append(str_z0)
    NAMES.append(str_obs)
    
    
print(NAMES)

datas=[]

for i in range(len(BPSnameslist)):
    datas.append(np.zeros_like(MSSFRnameslist))
    datas.append(np.zeros_like(MSSFRnameslist))


for ind_dco, DCOtype in enumerate(['BHBH', 'BHNS', 'NSNS']):    
    
	df = pd.DataFrame(data=datas, index=NAMES).T
	df.index =   MSSFRnameslist
	df.index.names = ['xyz']
	df.columns.names = ['mu']



	df.to_csv('RatesMSSFRandBPS_' + DCOtype + '.csv')