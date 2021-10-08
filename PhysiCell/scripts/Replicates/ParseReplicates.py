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
    
temp2_CD8 = []  # 3
temp2_mac = []  # 4
temp2_neut= []  # 5
temp2_DC = []   # 6
temp2_CD4 = []   # 7
temp2_fib = []   # 8
temp2_vir = []   # env
temp2_IFN = []   # env
temp2_Ig = []   # env
temp2_pI = []   # env
temp2_aI = []   # env
temp2_col = []   # env

for j in range(12): 

    temp1_CD8 = []  # 3
    temp1_mac = []  # 4
    temp1_neut= []  # 5
    temp1_DC = []   # 6
    temp1_CD4 = []   # 7
    temp1_fib = []   # 8
    temp1_vir = []   # env
    temp1_IFN = []   # env
    temp1_Ig = []   # env
    temp1_pI = []   # env
    temp1_aI = []   # env
    temp1_col = []   # env

    for k in range(361):
        str_name = 'output{:08d}.xml'.format(k)
        path = 'output_R'+str("%02d"%j)
        mcds = pyMCDS(str_name, path)  # /case1/run3/output

        cycle = mcds.data['discrete_cells']['cycle_model']
        cycle = cycle.astype(int)
        cell_type = mcds.data['discrete_cells']['cell_type']
        cell_type = cell_type.astype(int)

        CD8 = np.where((cell_type ==3) & (cycle < 100))
        mac = np.where((cell_type ==4) & (cycle < 100))
        neut = np.where((cell_type ==5) & (cycle < 100))
        DC = np.where((cell_type ==6) & (cycle < 100))
        CD4 = np.where((cell_type ==7) & (cycle < 100))
        fib = np.where((cell_type ==8) & (cycle < 100))
        virI = np.sum(mcds.data['discrete_cells']['assembled_virion'])
        virm = mcds.data['continuum_variables']['virion']
        virs = np.sum(virm['data'])*8000
        vir = np.add(virI,virs)
        IFNm = mcds.data['continuum_variables']['interferon 1']
        IFN = np.sum(IFNm['data'])*8000
        Igm = mcds.data['continuum_variables']['Ig']
        Ig = np.sum(Igm['data'])*8000
        pIm = mcds.data['continuum_variables']['pro-inflammatory cytokine']
        pI = np.sum(pIm['data'])*8000
        aIm = mcds.data['continuum_variables']['anti-inflammatory cytokine']
        aI = np.sum(aIm['data'])*8000
        colm = mcds.data['continuum_variables']['collagen']
        col = np.sum(colm['data'])*8000

        temp1_CD8.append( len(CD8[0]) )
        temp1_mac.append( len(mac[0]) ) 
        temp1_neut.append( len(neut[0]) ) 
        temp1_DC.append( len(DC[0]) ) 
        temp1_CD4.append( len(CD4[0]) ) 
        temp1_fib.append( len(fib[0]) ) 
        temp1_vir.append( vir ) 
        temp1_IFN.append( IFN )
        temp1_Ig.append( Ig ) 
        temp1_pI.append( pI ) 
        temp1_aI.append( aI ) 
        temp1_col.append( col ) 
        
    temp2_CD8.append(temp1_CD8)
    temp2_mac.append(temp1_mac)
    temp2_neut.append(temp1_neut)
    temp2_DC.append(temp1_DC)
    temp2_CD4.append(temp1_CD4)
    temp2_fib.append(temp1_fib)
    temp2_vir.append( temp1_vir )
    temp2_IFN.append( temp1_IFN )
    temp2_Ig.append( temp1_Ig )
    temp2_pI.append(temp1_pI)
    temp2_aI.append( temp1_aI ) 
    temp2_col.append( temp1_col )
 
aCD8 = np.asarray(temp2_CD8)
amac = np.asarray(temp2_mac)
aneut = np.asarray(temp2_neut)
aDC = np.asarray(temp2_DC)
aCD4 = np.asarray(temp2_CD4)
afib = np.asarray(temp2_fib)
avir = np.asarray(temp2_vir)
aIFN = np.asarray(temp2_IFN)
aIg = np.asarray(temp2_Ig)
apI = np.asarray(temp2_pI)
aaI = np.asarray(temp2_aI)
acol = np.asarray(temp2_col)

meanCD8 = np.mean(aCD8, axis=0)
meanmac = np.mean(amac, axis=0)
meanneut = np.mean(aneut, axis=0)
meanDC = np.mean(aDC, axis=0)
meanCD4 = np.mean(aCD4, axis=0)
meanfib = np.mean(afib, axis=0)
meanvir = np.mean(avir, axis=0)
meanIFN = np.mean(aIFN, axis=0)
meanIg = np.mean(aIg, axis=0)
meanpI = np.mean(apI, axis=0)
meanaI = np.mean(aaI, axis=0)
meancol = np.mean(acol, axis=0)

stdCD8 = np.std(aCD8, axis=0)
stdmac = np.std(amac, axis=0)
stdneut = np.std(aneut, axis=0)
stdDC = np.std(aDC, axis=0)
stdCD4 = np.std(aCD4, axis=0)
stdfib = np.std(afib, axis=0)
stdvir = np.std(avir, axis=0)
stdIFN = np.std(aIFN, axis=0)
stdIg = np.std(aIg, axis=0)
stdpI = np.std(apI, axis=0)
stdaI = np.std(aaI, axis=0)
stdcol = np.std(acol, axis=0)

data.append( np.vstack((meanCD8, meanmac, meanneut, meanDC, meanCD4, meanfib, meanvir, meanIFN, meanIg, meanpI, meanaI, meancol, stdCD8, stdmac, stdneut, stdDC, stdCD4, stdfib, stdvir, stdIFN, stdIg,stdpI,stdaI,stdcol)) )


# In[10]:

timedata = np.asarray(data)

sio.savemat('timeReplicate.mat', {'timedata':timedata})

