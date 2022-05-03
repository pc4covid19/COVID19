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


for q in range(10):
    data = []
        
    temp2_CD8 = []  # 3
    temp2_mac = []  # 4
    temp2_macM2 = []  # 4
    temp2_maci = []  # 4
    temp2_mach = []  # 4
    temp2_macex = []  # 4
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
    temp2_epi = []  # 1

    for j in range(12): 

        temp1_CD8 = []  # 3
        temp1_mac = []  # 4
        temp1_macM2 = []  # 4
        temp1_maci = []  # 4
        temp1_mach = []  # 4
        temp1_macex = []  # 4
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
        temp1_epi = []  # 1

        for k in range(13):
            str_name = 'output{:08d}.xml'.format(k)
            path = 'output_S'+str("%06d"%q)+'_R'+str("%02d"%j)
            mcds = pyMCDS(str_name, path)  # /case1/run3/output

            cycle = mcds.data['discrete_cells']['cycle_model']
            cycle = cycle.astype(int)
            phase = mcds.data['discrete_cells']['ability_to_phagocytose_infected_cell']
            phase = phase.astype(int)
            active = mcds.data['discrete_cells']['activated_immune_cell']
            active = active.astype(int)
            ex = mcds.data['discrete_cells']['M2_phase']
            ex = ex.astype(int)
            ex2 = mcds.data['discrete_cells']['total_volume']
            ex2 = ex.astype(int)
            cell_type = mcds.data['discrete_cells']['cell_type']
            cell_type = cell_type.astype(int)

            CD8 = np.where((cell_type ==3) & (cycle < 100))
            mac1 = np.where((cell_type ==4) & (cycle < 100) & (active == 1) & (ex == 0))
            mac2 = np.where((cell_type ==4) & (cycle < 100) & (ex == 1))
            mac3 = np.where((cell_type ==4) & (cycle < 100) & (active == 0))
            mac4 = np.where((cell_type ==4) & (cycle < 100) & (phase == 1))
            mac5 = np.where((cell_type ==4) & (cycle < 100) & (ex2 > 6500))
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
            epi = np.where((cell_type ==1) & (cycle < 100))

            temp1_CD8.append( len(CD8[0]) )
            temp1_mac.append( len(mac1[0]) ) 
            temp1_macM2.append( len(mac2[0]) ) 
            temp1_maci.append( len(mac3[0]) ) 
            temp1_mach.append( len(mac4[0]) ) 
            temp1_macex.append( len(mac5[0]) ) 
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
            temp1_epi.append( len(epi[0]) ) 
            
        temp2_CD8.append(temp1_CD8)
        temp2_mac.append(temp1_mac)
        temp2_macM2.append(temp1_macM2)
        temp2_maci.append(temp1_maci)
        temp2_mach.append(temp1_mach)
        temp2_macex.append(temp1_macex)
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
        temp2_epi.append( temp1_epi )
     
    aCD8 = np.asarray(temp2_CD8)
    amac = np.asarray(temp2_mac)
    amacM2 = np.asarray(temp2_macM2)
    amaci = np.asarray(temp2_maci)
    amach = np.asarray(temp2_mach)
    amacex = np.asarray(temp2_macex)
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
    aepi = np.asarray(temp2_epi)

    meanCD8 = np.mean(aCD8, axis=0)
    meanmac = np.mean(amac, axis=0)
    meanmacM2 = np.mean(amacM2, axis=0)
    meanmaci = np.mean(amaci, axis=0)
    meanmach = np.mean(amach, axis=0)
    meanmacex = np.mean(amacex, axis=0)
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
    meanepi = np.mean(aepi, axis=0)

    stdCD8 = np.std(aCD8, axis=0)
    stdmac = np.std(amac, axis=0)
    stdmacM2 = np.std(amacM2, axis=0)
    stdmaci = np.std(amaci, axis=0)
    stdmach = np.std(amach, axis=0)
    stdmacex = np.std(amacex, axis=0)
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
    stdepi = np.std(aepi, axis=0)

    data.append( np.vstack((meanCD8, meanmac, meanmacM2, meanmaci, meanmach, meanmacex, meanneut, meanDC, meanCD4, meanfib, meanvir, meanIFN, meanIg, meanpI, meanaI, meancol,meanepi, stdCD8, stdmac, stdmacM2, stdmaci,stdmach,stdmacex, stdneut, stdDC, stdCD4, stdfib, stdvir, stdIFN, stdIg,stdpI,stdaI,stdcol,stdepi)) )


    # In[10]:

    timedata = np.asarray(data)
    sio.savemat('timeReplicate{:06d}.mat'.format(q), {'timedata':timedata})