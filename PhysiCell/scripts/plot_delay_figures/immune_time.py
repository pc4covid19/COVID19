#!/usr/bin/env python
# coding: utf-8

# In[1]:


import scipy.io as sio
import numpy as np
from pyMCDS import pyMCDS


# In[ ]:


data = []
for i in range(1, 5):
    
    temp2_CD8 = []  # 3
    temp2_mac = []  # 4
    temp2_neut= []  # 5
    temp2_DC = []   # 6
    temp2_CD4 = []   # 7
    temp2_fib = []   # 8

    for j in range(1, 11): 

        temp1_CD8 = []  # 3
        temp1_mac = []  # 4
        temp1_neut= []  # 5
        temp1_DC = []   # 6
        temp1_CD4 = []   # 7
        temp1_fib = []   # 8

        for k in range(181):
            str_name = 'output{:08d}.xml'.format(k)
            path = '../run{:d}/run{:d}_{:d}/output'.format(i, i, j)
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

            temp1_CD8.append( len(CD8[0]) )
            temp1_mac.append( len(mac[0]) ) 
            temp1_neut.append( len(neut[0]) ) 
            temp1_DC.append( len(DC[0]) ) 
            temp1_CD4.append( len(CD4[0]) ) 
            temp1_fib.append( len(fib[0]) ) 
            
        temp2_CD8.append(temp1_CD8)
        temp2_mac.append(temp1_mac)
        temp2_neut.append(temp1_neut)
        temp2_DC.append(temp1_DC)
        temp2_CD4.append(temp1_CD4)
        temp2_fib.append(temp1_fib)
     
    aCD8 = np.asarray(temp2_CD8)
    amac = np.asarray(temp2_mac)
    aneut = np.asarray(temp2_neut)
    aDC = np.asarray(temp2_DC)
    aCD4 = np.asarray(temp2_CD4)
    afib = np.asarray(temp2_fib)

    meanCD8 = np.mean(aCD8, axis=0)
    meanmac = np.mean(amac, axis=0)
    meanneut = np.mean(aneut, axis=0)
    meanDC = np.mean(aDC, axis=0)
    meanCD4 = np.mean(aCD4, axis=0)
    meanfib = np.mean(afib, axis=0)

    stdCD8 = np.std(aCD8, axis=0)
    stdmac = np.std(amac, axis=0)
    stdneut = np.std(aneut, axis=0)
    stdDC = np.std(aDC, axis=0)
    stdCD4 = np.std(aCD4, axis=0)
    stdfib = np.std(afib, axis=0)

    data.append( np.vstack((meanCD8, meanmac, meanneut, meanDC, meanCD4, meanfib, stdCD8, stdmac, stdneut, stdDC, stdCD4, stdfib)) )


# In[10]:

timedata = np.asarray(data)

sio.savemat('time_1_4.mat', {'timedata':timedata})

