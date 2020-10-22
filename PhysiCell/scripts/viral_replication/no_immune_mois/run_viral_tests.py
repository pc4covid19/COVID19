import xmltodict, shutil, sys, os, subprocess, pickle
sys.path.append('../../../')
from pyMCDS import pyMCDS
import numpy as np

conf_file = "Viral_replication_no_immune.xml"

class Sim:
    def __init__(self, name, config):
        self.name = name
        self.conf_file = config
        with open(self.conf_file, "r") as f:
            lines = f.read()
            self.xmldict = xmltodict.parse(lines)

    @property
    def moi(self):
        return self._moi

    @moi.setter
    def moi(self, moi):
        # set xml option for moi
        self.xmldict['PhysiCell_settings']['user_parameters']['multiplicity_of_infection']['#text'] = moi
        self._moi = moi

    @property
    def output(self):
        return self._output

    @output.setter
    def output(self, output):
        # folder setup
        self.xmldict['PhysiCell_settings']['save']['folder'] = output
        self._output = output

    def run(self):
        # make folder
        fold_name = "viral_test_{}".format(self.name)
        os.mkdir(fold_name)
        # copy executable
        shutil.copy("../../../COVID19", fold_name)
        # get in there and prep the rest
        os.chdir(fold_name)
        os.mkdir(self.output)
        # write current xml 
        self.mod_conf = "{}_conf.xml".format(self.name)
        with open(self.mod_conf, "w") as f:
            f.write(xmltodict.unparse(self.xmldict))
        # run simulator
        rc = subprocess.run([os.path.abspath("COVID19"), self.mod_conf])
        assert rc.returncode == 0, "Simulator failed to run {}".format(os.getcwd())
        # analyze
        self.analyze()
        # go back to our main folder
        os.chdir(parent_dir)

    def analyze(self):
        # process output
        self.virion_dat = []
        self.cell_dat = []
        self.disc_vir = []
        self.immune_dat = []
        self.cont_immune = []
        self.frac_inf = []
        self.inf_cells = []
        self.assembled_vir = []
        for i in range(0,144):
            mcds = pyMCDS("output{:08d}.xml".format(i), "output/")
            total_volume = 8000 
            virion_dens = mcds.data['continuum_variables']['virion']['data']
            virions = (virion_dens*total_volume).sum()
            self.virion_dat.append(virions)

            vir_arr = mcds.data['discrete_cells']['virion']

            infected = vir_arr >= 1
            frac_inf_tp = len(np.where(infected)[0])/len(vir_arr) if len(vir_arr) != 0 else 0
            self.frac_inf.append(frac_inf_tp)
            self.inf_cells.append(len(np.where(infected)[0]))

            disc_vir_cnt = vir_arr.sum()
            self.disc_vir.append(disc_vir_cnt)

            assembl_vir_cnt = mcds.data['discrete_cells']['assembled_virion'].sum()
            self.assembled_vir.append(assembl_vir_cnt)

            live_cell_cnt = (mcds.data['discrete_cells']['cell_type']==1).sum()
            self.cell_dat.append(live_cell_cnt)

            immune_cnt = (mcds.data['discrete_cells']['cell_type']>1).sum()
            self.immune_dat.append(immune_cnt)

            cont_immune_cnt = mcds.data['continuum_variables']['interferon 1']['data'].sum()
            cont_immune_cnt += mcds.data['continuum_variables']['pro-inflammatory cytokine']['data'].sum()
            cont_immune_cnt += mcds.data['continuum_variables']['chemokine']['data'].sum()
            cont_immune_cnt += mcds.data['continuum_variables']['debris']['data'].sum()
            cont_immune_cnt = cont_immune_cnt * total_volume
            self.cont_immune.append(cont_immune_cnt)

        self.save_dat()

    def save_dat(self):
        dat = (self.virion_dat, self.disc_vir, self.cell_dat, \
               self.immune_dat, self.cont_immune, self.frac_inf,\
               self.assembled_vir, self.inf_cells)
        with open("data_arr_{}.pickle".format(self.name), "wb") as f:
            pickle.dump(dat, f)


parent_dir = os.getcwd()
sims = {}
for moi in ['10', '1', '0.1', '0.01']:
    sim_obj = Sim(moi, conf_file)
    sim_obj.moi = moi
    sim_obj.output = "output"
    sim_obj.run()
    sims[moi] = sim_obj
