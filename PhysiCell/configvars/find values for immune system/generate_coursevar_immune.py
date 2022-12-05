# This script provides parameter samples from latin hypercube. The script creates
# a new folder (subdirectory) for each set of parameters, makes changes to a default
# configuration (.xml) file using specified parameter values (in an accompanying .txt file),
# copies the new config file into the new folder

import xml.etree.ElementTree as ET
from shutil import copyfile
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import random
import time

os.chdir('../../')

def generate_parSamples(parameters, default_value, variation, Samples_number,Replicas_number, fileOut):

    file = open(fileOut, "w")
    samples = np.logspace(-4, 0, num=Samples_number)
    #define max values
    macH=84    
    IgH=0.19
    CD8H=250
    CD4H=9
    DCH=48
    DMH=0.1
    TCH=1.1E5
    TH1H=13
    TH2H=3.5E3
    BcH=5.2e03
    PsH=2.5e04
    DLH=0.084
    #define base values
    macL=50  
    IgL=0
    CD8L=0
    CD4L=0
    DCL=28
    DML=0
    TCL=0
    TH1L=0
    TH2L=0
    BcL=0
    PsL=0
    DLL=0
    #set the samples to supply
    mac_s=macH*samples+(1-samples)*macL
    Ig_s=IgH*samples+(1-samples)*IgL
    CD8_s=CD8H*samples+(1-samples)*CD8L
    CD4_s=CD4H*samples+(1-samples)*CD4L
    DC_s=DCH*samples+(1-samples)*DCL
    DM_s=DMH*samples+(1-samples)*DML
    TC_s=TCH*samples+(1-samples)*TCL
    TH1_s=TH1H*samples+(1-samples)*TH1L
    TH2_s=TH2H*samples+(1-samples)*TH2L
    Bc_s=BcH*samples+(1-samples)*BcL
    Ps_s=PsH*samples+(1-samples)*PsL
    DL_s=DLH*samples+(1-samples)*DLL
    #Write file with samples
    for sample_id in range(Samples_number):
        for replica_id in range(Replicas_number):
            folder = 'output_S'+str("%06d"%sample_id)+'_R'+str("%02d"%replica_id)
            file.write("folder"+" "+folder+"\n")
            # set system time as seed
            # create a seed
            seed_value = round(random.randrange(round(sys.maxsize/100000000000000)))
            omp_num_threads = 8
            mac_si=round(mac_s[sample_id])
            CD8_si=round(CD8_s[sample_id])
            CD4_si=round(CD4_s[sample_id])
            DC_si=round(DC_s[sample_id])
            # set reduced size for data and output
            # save this seed somewhere. So if you like the result you can use this seed to reproduce it
            # Set of parameters
            for id_par in range(0, len(parameters)):
                file.write('number_of_macrophages'+" "+str(mac_si)+"\n")                
                file.write('microenvironment_setup/variable[@name="Ig"]/initial_condition'+" "+str(Ig_s[sample_id])+"\n")
                file.write('number_of_CD8_Tcells'+" "+str(CD8_si)+"\n")
                file.write('number_of_CD4_Tcells'+" "+str(CD4_si)+"\n")
                file.write('number_of_DCs'+" "+str(DC_si)+"\n")
                file.write('DM_init'+" "+str(DM_s[sample_id])+"\n")
                file.write('TC_init'+" "+str(TC_s[sample_id])+"\n")
                file.write('TH1_init'+" "+str(TH1_s[sample_id])+"\n")
                file.write('TH2_init'+" "+str(TH2_s[sample_id])+"\n")
                file.write('Bc_init'+" "+str(Bc_s[sample_id])+"\n")
                file.write('Ps_init'+" "+str(Ps_s[sample_id])+"\n")
                file.write('DL_init'+" "+str(DL_s[sample_id])+"\n")
                file.write("random_seed"+" "+str(seed_value)+"\n")
                file.write("omp_num_threads"+" "+str(omp_num_threads)+"\n")
            file.write("#"+"\n")
    file.close()

def generate_configXML(params_file):
    xml_file_in = 'config/PhysiCell_settings.xml'
    xml_file_out = 'config/tmp.xml'
    copyfile(xml_file_in,xml_file_out)
    tree = ET.parse(xml_file_out)
    xml_root = tree.getroot()
    first_time = True
    output_dirs = []
    #reduce save intervals
    SVG_interval = "2880"
    full_data_interval = "720"
    
    with open(params_file) as f:
        for line in f:
            print(len(line),line)
            if (line[0] == '#'):
                continue
            (key, val) = line.split()
            print(key,val)
            print('test')
            if (key == 'folder'):
                if first_time:  # we've read the 1st 'folder'
                    first_time = False
                else:  # we've read  additional 'folder's
                    # write the config file to the previous folder (output) dir and start a simulation
                    print('---write (previous) config file and start its sim')
                    tree.write(xml_file_out)
                xml_file_out = val + '/config.xml'  # copy config file into the output dir
                output_dirs.append(val)
            if ('.' in key):
                k = key.split('.')
                uep = xml_root
                for idx in range(len(k)):
                    uep = uep.find('.//' + k[idx])  # unique entry point (uep) into xml
                uep.text = val
            else:
                if (key == 'folder' and not os.path.exists(val)):
                    print('creating ' + val)
                    os.makedirs(val)
                xml_root.find('.//' + key).text = val
            xml_root.find('save/full_data/interval').text = full_data_interval
            xml_root.find('save/SVG/interval').text = SVG_interval
    tree.write(xml_file_out)
    print(output_dirs)

if __name__ == '__main__':
    parameters = np.array(['total_variation']) #not needed
    default_value = np.array([1])
    
    variation = 1 #percentage variation wanted
    file = "ParameterSamples.txt"
    Samples_number = 6
    Replicas_number = 10
    # Generate samples from Latin Hypercube
    generate_parSamples(parameters, default_value, variation, Samples_number,Replicas_number, file)
    # Create .xml and folder to each simulation
    generate_configXML(file)