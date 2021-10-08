#!/usr/bin/env python
# coding: utf-8

# In[1]:


import scipy.io as sio
import numpy as np
import matplotlib.pyplot as plt


# In[2]:


a = sio.loadmat('time_1_4.mat')
cells = a['timedata']


# In[3]:


cells.shape


# In[4]:


t = np.linspace(0, 180/12, 181)
title_list = ['No delay', 'Half day delay (ddT=0)', 'Half day delay (ddT=0.25 d)', 'One day delay']


# In[5]:


# temp2_CD8 = []  # 3
# temp2_mac = []  # 4
# temp2_neut= []  # 5
# temp2_DC = []   # 6
# temp2_CD4 = []   # 7
# temp2_fib = []   # 8


# In[6]:


immune_cells = ['CD8 T', 'Mac', 'Neut', 'DC', 'CD4 T', 'Fib']

colorv = ['red' , 'seagreen', 'cyan', '#810f7c', 'orange', 'blueviolet']
# In[15]:


fig, ax = plt.subplots(nrows=1, ncols=4, figsize=(22, 4))

for i in range(4):
    
    for j in range( len(immune_cells) ):
        ax[i].plot(t, cells[i,j,:], color= colorv[j], linewidth=2)
        ax[i].fill_between(t, cells[i,j,:]-cells[i,j+6,:], cells[i,j,:]+cells[i,j+6,:], 
                           label=immune_cells[j], color= colorv[j], alpha=0.35)    
    # ax[i].plot(t, meanR, 'green')
    # ax[i].fill_between(t, meanR-stdR, meanR+stdR, alpha=0.2, label='resistant', facecolor='green')
    
    if i==3:  
        ax[i].legend(fontsize=16, loc=9, bbox_to_anchor=(1.2, 1.))
        
    title_name = title_list[i]    
    ax[i].set_title(title_name, fontsize=18)  
    ax[i].set_xlabel('Time (days)', fontsize=16)
    if i==0: 
        ax[i].set_ylabel('# of immune cells', fontsize=16)
        
    ax[i].set_ylim([-50, 400])
    ax[i].tick_params(axis='both', which='major', labelsize=16)

    # Customize the major grid
    ax[i].grid(which='major', linestyle='solid', linewidth='2', color='w')
    ax[i].set_facecolor("#EEEEEE")  # #E6E6E6, #D3D3D3

# fig.suptitle("Comparison of different time delay", fontsize=22)
# plt.subplots_adjust(left=0.1, top=0.8, wspace = 0.2, hspace = 0.4)
plt.tight_layout()
fig.savefig('population-dynamics.png', dpi=600, pad_inches=0.1, bbox_inches='tight')  # dpi=600, 


# In[ ]:





# ## plt.subplots practice

# In[76]:


fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(10, 5))
rows, cols = axes.shape

for i in range(rows):
    for j in range(cols):
        name = str(i) + ' ' + 'and' + ' ' + str(j) 
        axes[i,j].plot(range(0, 10))
        axes[i,j].set_title(name)

plt.subplots_adjust(left=0.1, top=0.9, wspace = 0.2, hspace = 0.4)
# plt.tight_layout()


# In[69]:


axes.shape

