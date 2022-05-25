#!/usr/bin/env python
# coding: utf-8

import scipy.io as sio
import numpy as np
import matplotlib.pyplot as plt
import os
os.chdir('../../')

for q in range(10):
    a = sio.loadmat('timeReplicate{:06d}.mat'.format(q))
    cells = a['timedata']
    cells.shape

    mac = np.squeeze(np.sum([cells[:,1:5,:]], axis=2))
    mac_v = np.squeeze(np.sum([cells[:,1+17:5+17,:]], axis=2))

    t = np.linspace(0, 12, 13)


    immune_cells = ['CD8 T', 'Mac', 'M2Mac', 'Maci', 'Mach', 'Macexh', 'Neut', 'DC', 'CD4 T', 'Fib', 'virion', 'IFN', 'Ig', 'pro-I', 'anti-I', 'collagen', 'epi']
    innate = ['totalMac', 'Mac', 'M2Mac', 'Maci', 'Mach', 'Neut']
    set1 = ['CD8 T', 'DC', 'CD4 T','Ig']
    counter1 = ['CD8 T', 'Mac', 'M2Mac', 'Maci', 'Mach', 'Macexh', 'Neut', 'DC', 'CD4 T', 'Fib']
    set2 = ['vir', 'IFN']
    set3 = ['Ig']
    set4 = ['pI', 'aI']
    repair = ['col', 'Fib']
    set6 = ['epi']

    innate_r = np.vstack((mac, np.squeeze(cells[:,1,:]), np.squeeze(cells[:,2,:]),np.squeeze(cells[:,3,:]),np.squeeze(cells[:,4,:]), np.squeeze(cells[:,6,:])))
    innate_v = np.vstack((mac_v, np.squeeze(cells[:,1+17,:]),np.squeeze(cells[:,2+17,:]),np.squeeze(cells[:,3+17,:]),np.squeeze(cells[:,4+17,:]), np.squeeze(cells[:,6+17,:])))

    immun_r = np.vstack((np.squeeze(cells[:,0,:]), np.squeeze(cells[:,7,:]),np.squeeze(cells[:,8,:]),np.squeeze(cells[:,12,:])))
    immun_v = np.vstack((np.squeeze(cells[:,0+17,:]), np.squeeze(cells[:,7+17,:]),np.squeeze(cells[:,8+17,:]),np.squeeze(cells[:,12+17,:])))

    repair_r = np.vstack((np.squeeze(cells[:,15,:]), np.squeeze(cells[:,9,:])))
    repair_v = np.vstack((np.squeeze(cells[:,15+17,:]), np.squeeze(cells[:,9+17,:])))

    Dat = ['DM', 'TC', 'TH1', 'TH2', 'BC', 'PS', 'DL']

    colorv = ['black' , 'lime', 'darkgreen', '#238b45', '#a8ddb5', 'cyan'] #change style on avg so better show variation darken dull color
    colorv2 = ['red' , 'deeppink', 'orange', 'blueviolet']


    fig, ax = plt.subplots(figsize=(6, 4))

    # plot physicell files

    for i in range(1):
        
        for j in range( len(innate_r) ):
            ax.plot(t, innate_r[j,:], color= colorv[j], linewidth=2)
            ax.fill_between(t, innate_r[j,:]-innate_v[j,:], innate_r[j,:]+innate_v[j,:], 
                               label=innate[j], color= colorv[j], alpha=0.35)   
        
        ax.legend(fontsize=12, loc=2, ncol=2)
              
        ax.set_xlabel('Time (days)', fontsize=16)

        ax.set_ylabel('# of immune cells', fontsize=16)
            
        ax.set_ylim([-50, 400])
        ax.tick_params(axis='both', which='major', labelsize=16)

        # Customize the major grid
        ax.grid(which='major', linestyle='solid', linewidth='2', color='w')
        ax.set_facecolor("#EEEEEE")  # #E6E6E6, #D3D3D3

    plt.tight_layout()
    fig.savefig('innateimmune{:06d}.png'.format(q), dpi=600, pad_inches=0.1, bbox_inches='tight')  # dpi=600,
    plt.cla()  
    plt.clf()

    fig, ax = plt.subplots(figsize=(6, 4))

    for i in range(1):
        
        for j in range( len(immun_r) ):
            ax.plot(t, immun_r[j,:], color= colorv2[j], linewidth=2)
            ax.fill_between(t, immun_r[j,:]-immun_v[j,:], immun_r[j,:]+immun_v[j,:], 
                               label=set1[j], color= colorv2[j], alpha=0.35)   
        
        ax.legend(fontsize=12, loc=2, ncol=2)
              
        ax.set_xlabel('Time (days)', fontsize=16)

        ax.set_ylabel('# of immune cells', fontsize=16)
            
        ax.set_ylim([-50, 400])
        ax.tick_params(axis='both', which='major', labelsize=16)

        # Customize the major grid
        ax.grid(which='major', linestyle='solid', linewidth='2', color='w')
        ax.set_facecolor("#EEEEEE")  # #E6E6E6, #D3D3D3

    plt.tight_layout()
    fig.savefig('immune{:06d}.png'.format(q), dpi=600, pad_inches=0.1, bbox_inches='tight')  # dpi=600,
    plt.cla()  
    plt.clf()

    for i in range(1):
        
        fig, ax1 = plt.subplots(figsize=(6, 4)) 
      
        color = 'C0'
        ax1.set_xlabel('Time (days)', fontsize=16) 
        ax1.set_ylabel('virion', color = color, fontsize=16) 
        for j in range( len(counter1), len(counter1)+1, 1):
            ax1.plot(t, cells[i,j,:], color = color, linewidth=2)
            ax1.fill_between(t, cells[i,j,:]-cells[i,j+17,:], cells[i,j,:]+cells[i,j+17,:], 
                               label=immune_cells[j], color='C0', alpha=0.35)  
        ax1.tick_params(axis ='y', which='major', labelsize=16, labelcolor = color) 
        ax1.tick_params(axis ='x', which='major', labelsize=16, labelcolor = color) 
      
        # Adding Twin Axes to plot using dataset_2
        ax2 = ax1.twinx() 
      
        color = 'C1'
        ax2.set_ylabel('IFN', color = color, fontsize=16) 
        for j in range( len(counter1)+1, len(counter1)+len(set2), 1):
            ax2.plot(t, cells[i,j,:], color = color, linewidth=2)
            ax2.fill_between(t, cells[i,j,:]-cells[i,j+17,:], cells[i,j,:]+cells[i,j+17,:], 
                               label=immune_cells[j], color='C1', alpha=0.35)  
        ax2.tick_params(axis ='y', which='major', labelsize=16, labelcolor = color) 
      
        # Adding title
        #plt.title('Use different y-axes on the left and right of a Matplotlib plot', fontweight ="bold") 
        #ax.grid(which='major', linestyle='solid', linewidth='2', color='w')
        #ax.set_facecolor("#EEEEEE")  # #E6E6E6, #D3D3D3

    plt.tight_layout()
    fig.savefig('populationvir{:06d}.png'.format(q), dpi=600, pad_inches=0.1, bbox_inches='tight')  # dpi=600, 
    plt.clf()

    fig, ax = plt.subplots(figsize=(6, 4))
    for i in range(1):
        # Creating plot with dataset_1
        fig, ax1 = plt.subplots(figsize=(6, 4)) 
      
        color = 'C0'
        ax1.set_xlabel('Time (days)', fontsize=16) 
        ax1.set_ylabel('pro-inflammatory cytokine', color = color, fontsize=16) 
        for j in range( len(counter1)+len(set2)+len(set3), len(counter1)+len(set2)+len(set3)+1 , 1):
            ax1.plot(t, cells[i,j,:], color = color, linewidth=2)
            ax1.fill_between(t, cells[i,j,:]-cells[i,j+17,:], cells[i,j,:]+cells[i,j+17,:], 
                               label=immune_cells[j], color='C0', alpha=0.35)  
        ax1.tick_params(axis ='y', which='major', labelsize=16, labelcolor = color) 
        ax1.tick_params(axis ='x', which='major', labelsize=16, labelcolor = color) 
      
        # Adding Twin Axes to plot using dataset_2
        ax2 = ax1.twinx() 
      
        color = 'C1'
        ax2.set_ylabel('anti-inflammatory cytokine', color = color, fontsize=16) 
        for j in range( len(counter1)+len(set2)+len(set3)+1, len(counter1)+len(set2)+len(set3)+2 , 1):
            ax2.plot(t, cells[i,j,:], color = color, linewidth=2)
            ax2.fill_between(t, cells[i,j,:]-cells[i,j+17,:], cells[i,j,:]+cells[i,j+17,:], 
                               label=immune_cells[j], color='C1', alpha=0.35)  
        ax2.tick_params(axis ='y', which='major', labelsize=16, labelcolor = color) 
      
        # Adding title
        #plt.title('Use different y-axes on the left and right of a Matplotlib plot', fontweight ="bold") 
        #ax.grid(which='major', linestyle='solid', linewidth='2', color='w')
        #ax.set_facecolor("#EEEEEE")  # #E6E6E6, #D3D3D3

    plt.tight_layout()
    fig.savefig('populationaIpI{:06d}.png'.format(q), dpi=600, pad_inches=0.1, bbox_inches='tight')  # dpi=600, 
    plt.clf()

    fig, ax = plt.subplots(figsize=(6, 4))
    for i in range(1):
        
        for j in range( len(counter1)+len(set2)+len(set4)+len(repair), len(counter1)+len(set2)+len(set4)+len(repair)+len(set6) , 1):
            ax.plot(t, cells[i,j,:], linewidth=2)
            ax.fill_between(t, cells[i,j,:]-cells[i,j+17,:], cells[i,j,:]+cells[i,j+17,:], 
                               label=immune_cells[j], alpha=0.35)   
        
        ax.legend(fontsize=12, ncol=2)
              
        ax.set_xlabel('Time (days)', fontsize=16)

        ax.set_ylabel('epi [Count]', fontsize=16)
            
        ax.set_ylim([0, 3000])
        ax.tick_params(axis='both', which='major', labelsize=16)

        # Customize the major grid
        ax.grid(which='major', linestyle='solid', linewidth='2', color='w')
        ax.set_facecolor("#EEEEEE")  # #E6E6E6, #D3D3D3

    plt.tight_layout()
    fig.savefig('populationepi{:06d}.png'.format(q), dpi=600, pad_inches=0.1, bbox_inches='tight')  # dpi=600, 
