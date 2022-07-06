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
    min = default_value - variation*default_value
    max = default_value + variation*default_value
    
    #generate variation from 0 to 1 for num samples
    samples = np.linspace(0, 1, num=Samples_number)
    samples = min + samples*(max-min)

    #Write file with samples
    for sample_id in range(Samples_number):
        for replica_id in range(Replicas_number):
            folder = 'output_S'+str("%06d"%sample_id)+'_R'+str("%02d"%replica_id)
            file.write("folder"+" "+folder+"\n")
            # set system time as seed
            # create a seed
            seed_value = round(random.randrange(sys.maxsize)/100000000)
            omp_num_threads = 1
            # set reduced size for data and output
            # save this seed somewhere. So if you like the result you can use this seed to reproduce it
            # Set of parameters
            for id_par in range(0, len(parameters)):
                file.write(parameters[id_par]+" "+str(samples[sample_id])+"\n")
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
    parameters = np.array(["macrophage_max_recruitment_rate"])
    default_value = np.array([2e-8])
    
    variation = 1
    file = "ParameterSamples.txt"
    Samples_number = 10
    Replicas_number = 6
    # Generate samples from Latin Hypercube
    generate_parSamples(parameters, default_value, variation, Samples_number,Replicas_number, file)
    # Create .xml and folder to each simulation
    generate_configXML(file)