#!/usr/bin/env python
# coding: utf-8

import scipy.io as sio
import numpy as np
import matplotlib.pyplot as plt
import os
os.chdir('../../')

a = sio.loadmat('timeReplicate.mat')
cells = a['timedata']

cells.shape

t = np.linspace(0, 12, 145)


immune_cells = ['CD8 T', 'Mac', 'M2Mac', 'Maci', 'Mach', 'Macexh', 'Neut', 'DC', 'CD4 T', 'Fib', 'virion', 'IFN', 'Ig', 'pro-I', 'anti-I', 'collagen', 'epi']
set1 = ['CD8 T', 'Mac', 'M2Mac', 'Maci', 'Mach', 'Macexh', 'Neut', 'DC', 'CD4 T', 'Fib']
set2 = ['vir', 'IFN']
set3 = ['Ig']
set4 = ['pI', 'aI']
set5 = ['col']
set6 = ['epi']

colorv = ['red' , 'blue', 'darkgreen', 'orange', 'yellow', 'cyan', 'white', 'white', 'white', 'white']
# [ 1 2 3 4 5 6 7 8 9 10 11 ]

fig, ax = plt.subplots(figsize=(6, 4))

for i in range(1):
    
    for j in range( len(set1) ):
        ax.plot(t, cells[i,j,:], color= colorv[j], linewidth=2)
        ax.fill_between(t, cells[i,j,:]-cells[i,j+17,:], cells[i,j,:]+cells[i,j+17,:], 
                           label=immune_cells[j], color= colorv[j], alpha=0.35)   
    
    ax.legend(fontsize=12, loc=2, ncol=2)
          
    ax.set_xlabel('Time (days)', fontsize=16)

    ax.set_ylabel('# of immune cells', fontsize=16)
        
    ax.set_ylim([-50, 400])
    ax.tick_params(axis='both', which='major', labelsize=16)

    # Customize the major grid
    ax.grid(which='major', linestyle='solid', linewidth='2', color='w')
    ax.set_facecolor("#EEEEEE")  # #E6E6E6, #D3D3D3

plt.tight_layout()
fig.savefig('populationimmune.png', dpi=600, pad_inches=0.1, bbox_inches='tight')  # dpi=600,
plt.cla()  
plt.clf()

for i in range(1):
    
    fig, ax1 = plt.subplots(figsize=(6, 4)) 
  
    color = 'C0'
    ax1.set_xlabel('Time (days)', fontsize=16) 
    ax1.set_ylabel('virion', color = color, fontsize=16) 
    for j in range( len(set1), len(set1)+1, 1):
        ax1.plot(t, cells[i,j,:], color = color, linewidth=2)
        ax1.fill_between(t, cells[i,j,:]-cells[i,j+17,:], cells[i,j,:]+cells[i,j+17,:], 
                           label=immune_cells[j], color='C0', alpha=0.35)  
    ax1.tick_params(axis ='y', which='major', labelsize=16, labelcolor = color) 
    ax1.tick_params(axis ='x', which='major', labelsize=16, labelcolor = color) 
  
    # Adding Twin Axes to plot using dataset_2
    ax2 = ax1.twinx() 
  
    color = 'C1'
    ax2.set_ylabel('IFN', color = color, fontsize=16) 
    for j in range( len(set1)+1, len(set1)+len(set2), 1):
        ax2.plot(t, cells[i,j,:], color = color, linewidth=2)
        ax2.fill_between(t, cells[i,j,:]-cells[i,j+17,:], cells[i,j,:]+cells[i,j+17,:], 
                           label=immune_cells[j], color='C1', alpha=0.35)  
    ax2.tick_params(axis ='y', which='major', labelsize=16, labelcolor = color) 
  
    # Adding title
    #plt.title('Use different y-axes on the left and right of a Matplotlib plot', fontweight ="bold") 
    #ax.grid(which='major', linestyle='solid', linewidth='2', color='w')
    #ax.set_facecolor("#EEEEEE")  # #E6E6E6, #D3D3D3

plt.tight_layout()
fig.savefig('populationvir.png', dpi=600, pad_inches=0.1, bbox_inches='tight')  # dpi=600, 
plt.clf()

fig, ax = plt.subplots(figsize=(6, 4))
for i in range(1):
    
    for j in range( len(set1)+len(set2), len(set1)+len(set2)+len(set3) , 1):
        ax.plot(t, cells[i,j,:], linewidth=2)
        ax.fill_between(t, cells[i,j,:]-cells[i,j+17,:], cells[i,j,:]+cells[i,j+17,:], 
                           label=immune_cells[j], alpha=0.35)   
    
    ax.legend(fontsize=12, ncol=2)
          
    ax.set_xlabel('Time (days)', fontsize=16)

    ax.set_ylabel('Antibody [Count]', fontsize=16)
        
    #ax.set_ylim([-50, 400])
    ax.tick_params(axis='both', which='major', labelsize=16)

    # Customize the major grid
    ax.grid(which='major', linestyle='solid', linewidth='2', color='w')
    ax.set_facecolor("#EEEEEE")  # #E6E6E6, #D3D3D3

plt.tight_layout()
fig.savefig('populationIg.png', dpi=600, pad_inches=0.1, bbox_inches='tight')  # dpi=600, 
plt.clf()

for i in range(1):
    # Creating plot with dataset_1
    fig, ax1 = plt.subplots(figsize=(6, 4)) 
  
    color = 'C0'
    ax1.set_xlabel('Time (days)', fontsize=16) 
    ax1.set_ylabel('pro-inflammatory cytokine', color = color, fontsize=16) 
    for j in range( len(set1)+len(set2)+len(set3), len(set1)+len(set2)+len(set3)+1 , 1):
        ax1.plot(t, cells[i,j,:], color = color, linewidth=2)
        ax1.fill_between(t, cells[i,j,:]-cells[i,j+17,:], cells[i,j,:]+cells[i,j+17,:], 
                           label=immune_cells[j], color='C0', alpha=0.35)  
    ax1.tick_params(axis ='y', which='major', labelsize=16, labelcolor = color) 
    ax1.tick_params(axis ='x', which='major', labelsize=16, labelcolor = color) 
  
    # Adding Twin Axes to plot using dataset_2
    ax2 = ax1.twinx() 
  
    color = 'C1'
    ax2.set_ylabel('anti-inflammatory cytokine', color = color, fontsize=16) 
    for j in range( len(set1)+len(set2)+len(set3)+1, len(set1)+len(set2)+len(set3)+2 , 1):
        ax2.plot(t, cells[i,j,:], color = color, linewidth=2)
        ax2.fill_between(t, cells[i,j,:]-cells[i,j+17,:], cells[i,j,:]+cells[i,j+17,:], 
                           label=immune_cells[j], color='C1', alpha=0.35)  
    ax2.tick_params(axis ='y', which='major', labelsize=16, labelcolor = color) 
  
    # Adding title
    #plt.title('Use different y-axes on the left and right of a Matplotlib plot', fontweight ="bold") 
    #ax.grid(which='major', linestyle='solid', linewidth='2', color='w')
    #ax.set_facecolor("#EEEEEE")  # #E6E6E6, #D3D3D3

plt.tight_layout()
fig.savefig('populationaIpI.png', dpi=600, pad_inches=0.1, bbox_inches='tight')  # dpi=600, 
plt.clf()


for i in range(1):
    
    for j in range( len(set1)+len(set2)+len(set3)+len(set4), len(set1)+len(set2)+len(set3)+len(set4)+len(set5) , 1):
        ax.plot(t, cells[i,j,:], linewidth=2)
        ax.fill_between(t, cells[i,j,:]-cells[i,j+17,:], cells[i,j,:]+cells[i,j+17,:], 
                           label=immune_cells[j], alpha=0.35)   
    
    ax.legend(fontsize=12, ncol=2)
          
    ax.set_xlabel('Time (days)', fontsize=16)

    ax.set_ylabel('Collagen [Count]', fontsize=16)
        
    #ax.set_ylim([-50, 400])
    ax.tick_params(axis='both', which='major', labelsize=16)

    # Customize the major grid
    ax.grid(which='major', linestyle='solid', linewidth='2', color='w')
    ax.set_facecolor("#EEEEEE")  # #E6E6E6, #D3D3D3

plt.tight_layout()
fig.savefig('populationcol.png', dpi=600, pad_inches=0.1, bbox_inches='tight')  # dpi=600, 
plt.clf()

fig, ax = plt.subplots(figsize=(6, 4))
for i in range(1):
    
    for j in range( len(set1)+len(set2)+len(set3)+len(set4)+len(set5), len(set1)+len(set2)+len(set3)+len(set4)+len(set5)+len(set6) , 1):
        ax.plot(t, cells[i,j,:], linewidth=2)
        ax.fill_between(t, cells[i,j,:]-cells[i,j+17,:], cells[i,j,:]+cells[i,j+17,:], 
                           label=immune_cells[j], alpha=0.35)   
    
    ax.legend(fontsize=12, ncol=2)
          
    ax.set_xlabel('Time (days)', fontsize=16)

    ax.set_ylabel('epi [Count]', fontsize=16)
        
    #ax.set_ylim([-50, 400])
    ax.tick_params(axis='both', which='major', labelsize=16)

    # Customize the major grid
    ax.grid(which='major', linestyle='solid', linewidth='2', color='w')
    ax.set_facecolor("#EEEEEE")  # #E6E6E6, #D3D3D3

plt.tight_layout()
fig.savefig('populationepi.png', dpi=600, pad_inches=0.1, bbox_inches='tight')  # dpi=600, 

