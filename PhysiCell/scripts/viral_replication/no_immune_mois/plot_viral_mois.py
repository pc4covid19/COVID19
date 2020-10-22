import sys
sys.path.append('../../../')
from pyMCDS import pyMCDS
import pickle, os
import numpy as np
import matplotlib.pyplot as plt
import seaborn

class DataLoader:
    def __init__(self, path):
        self.path = path
        self.cur_path = os.getcwd()
        os.chdir(self.path)
        self.data = {}
        self.load_data(path)
        os.chdir(self.cur_path)

    def load_data(self, path):
        for moi in [10, 1, 0.1, 0.01]:
            with open("viral_test_{}/data_arr_{}.pickle".format(moi,moi), "rb") as f:
                dat = pickle.load(f)
            virion_dat, disc_vir, cell_dat, immune_dat, cont_immune, frac_inf, assembled_vir, inf_cells = dat
            self.data[moi] = {}
            self.data[moi]['env_virion'] = np.array(virion_dat)
            self.data[moi]['discrete_virion'] = np.array(disc_vir)
            self.data[moi]['cells'] = np.array(cell_dat)
            self.data[moi]['inf_cells'] = np.array(inf_cells)
            self.data[moi]['assembled_vir'] = np.array(assembled_vir)
            self.data[moi]['uninf_cells'] = self.data[moi]['cells'] - self.data[moi]['inf_cells']
            init_cell_cnt = cell_dat[0]
            self.data[moi]['dead_cells'] = init_cell_cnt - self.data[moi]['cells']

    def get_data(self, moi=10, key="env_virion"):
        return self.data[moi][key]


keys = ['cells', 'env_virion', 'discrete_virion', 'assembled_vir'] 
orig_data = DataLoader('.')

fig, axs = plt.subplots(nrows=4, ncols=1, sharex=True, sharey=True, figsize=(10,10))

for ir, row_ax in enumerate(axs):
    key = keys[ir]
    orig_ax = row_ax
    orig_ax.set_ylim((1e-2,1e8))
    orig_ax.set(yscale="log")
    orig_ax.set_ylabel("{} count".format(key))

    for moi in [10, 1, 0.1, 0.01]:
        orig_ax.plot(orig_data.get_data(moi=moi, key=key), label=moi)
    if ir == 0:
        orig_ax.legend(frameon=False)
    if ir == 3:
        orig_ax.set_xlabel("time (hours)")

plt.savefig("moi_plots.png")
plt.close()

keys = ['uninf_cells', 'inf_cells', 'dead_cells']
fig, axs = plt.subplots(nrows=3, ncols=1, sharex=True, sharey=True, figsize=(10,10))

for ir, row_ax in enumerate(axs):
    key = keys[ir]
    orig_ax = row_ax
    orig_ax.set_ylim((0,3e3))
    orig_ax.set_ylabel("{} count".format(key))
    for moi in [10, 1, 0.1, 0.01]:
        orig_ax.plot(orig_data.get_data(moi=moi, key=key), label=moi)
    if ir == 0:
        orig_ax.legend(frameon=False)
    if ir == 2:
        orig_ax.set_xlabel("time (hours)")

plt.savefig("cell_plots.png")
plt.close()


# keys = ['uninf_cells', 'inf_cells', 'dead_cells']
# colors = ["red", "green", "blue"]
# for moi in [10,1,0.01, 0.1]:
#     fig, axs = plt.subplots(figsize=(10,10))
#     for ik,key in enumerate(keys):
#         axs.set_ylim((0,3e3))
#         axs.set_ylabel("count")
#         axs.set_xlabel("time (hours)")
#         axs.set_title("MOI = {}".format(moi))
#         axs.plot(orig_data.get_data(moi=moi, key=key), label="or_{}".format(key), color=colors[ik])
#         axs.plot(custom_data.get_data(moi=moi, key=key), label="ch_{}".format(key), color=colors[ik], linestyle="--")
#         axs.legend(frameon=False)
#     
#     plt.savefig("cell_per_moi_{}.png".format(moi))
#     plt.close()
# 
# keys = ['cells', 'env_virion', 'discrete_virion', 'assembled_vir'] 
# colors = ["red", "green", "blue", "black"]
# for moi in [10,1,0.01, 0.1]:
#     fig, axs = plt.subplots(figsize=(10,10))
#     for ik,key in enumerate(keys):
#         axs.set_ylim((1e-2,1e8))
#         axs.set(yscale="log")
#         axs.set_ylabel("count")
#         axs.set_xlabel("time (hours)")
#         axs.set_title("MOI = {}".format(moi))
#         od = orig_data.get_data(moi=moi, key=key)
#         if len(np.where(od==0)[0]) > 1:
#             if np.where(od==0)[0][0] > 0:
#                 od = od[:np.where(od==0)[0][0]]
#             else:
#                 od = od[:np.where(od==0)[0][1]]
#         cd = custom_data.get_data(moi=moi, key=key)
#         if len(np.where(cd==0)[0]) > 1:
#             if np.where(cd==0)[0][0] > 0:
#                 cd = cd[:np.where(cd==0)[0][0]]
#             else:
#                 cd = cd[:np.where(cd==0)[0][1]]
#         axs.plot(od, label="or_{}".format(key), color=colors[ik])
#         axs.plot(cd, label="ch_{}".format(key), color=colors[ik], linestyle="--")
#         axs.legend(frameon=False)
#     
#     plt.savefig("vir_per_moi_{}.png".format(moi))
#     plt.close()
