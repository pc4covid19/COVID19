#!/usr/bin/env python
# coding: utf-8

# In[1]:

import scipy.io as sio
import numpy as np
import os
import sys
os.chdir('../../')
sys.path.append('.')
from pyMCDS import pyMCDS


# In[ ]:

data = []
    
temp1_DM = []  # 3
temp1_TC = []  # 3
temp1_TH1 = []  # 3
temp1_TH2 = []  # 3
temp1_BC = []  # 3
temp1_PS = []  # 3
temp1_DL = []  # 3

for j in range(12): 

    file_name = 'dm_tc.dat'
    path = 'output_R'+str("%02d"%j)
    
    os.chdir(path)
    d = np.loadtxt(file_name)
    os.chdir('../')

    DM = d[:,0]
    TC = d[:,1]
    TH1 = d[:,2]
    TH2 = d[:,3]
    BC = d[:,6]
    PS = d[:,7]
    DL = d[:,12]
        
    temp1_DM.append( DM )
    temp1_TC.append( TC )
    temp1_TH1.append( TH1 )
    temp1_TH2.append( TH2 )
    temp1_BC.append( BC )
    temp1_PS.append( PS )
    temp1_DL.append( DL )
 
aDM = np.asarray(temp1_DM)
aTC = np.asarray(temp1_TC)
aTH1 = np.asarray(temp1_TH1)
aTH2 = np.asarray(temp1_TH2)
aBC = np.asarray(temp1_BC)
aPS = np.asarray(temp1_PS)
aDL = np.asarray(temp1_DL)

meanDM = np.mean(aDM, axis=0)
meanTC = np.mean(aTC, axis=0)
meanTH1 = np.mean(aTH1, axis=0)
meanTH2 = np.mean(aTH2, axis=0)
meanBC = np.mean(aBC, axis=0)
meanPS = np.mean(aPS, axis=0)
meanDL = np.mean(aDL, axis=0)

stdDM = np.std(aDM, axis=0)
stdTC = np.std(aTC, axis=0)
stdTH1 = np.std(aTH1, axis=0)
stdTH2 = np.std(aTH2, axis=0)
stdBC = np.std(aBC, axis=0)
stdPS = np.std(aPS, axis=0)
stdDL = np.std(aDL, axis=0)

data.append( np.vstack((meanDM, meanTC, meanTH1, meanTH2, meanBC, meanPS, meanDL, stdDM, stdTC, stdTH1, stdTH2, stdBC, stdPS, stdDL)) )

timedata = np.asarray(data)

sio.savemat('timeReplicateDat.mat', {'timedata':timedata})

