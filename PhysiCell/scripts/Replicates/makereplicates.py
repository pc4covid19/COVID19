# This script provides parameter samples from latin hypercube. The script creates
# a new folder (subdirectory) for each set of parameters, makes changes to a default
# configuration (.xml) file using specified parameter values (in an accompanying .txt file),
# copies the new config file into the new folder

import xml.etree.ElementTree as ET
from shutil import copyfile
import numpy as np
import os
import sys

Replicas_number= 12

os.chdir('../../')
fileOut='test.txt'
file = open(fileOut, "w")
for replica_id in range(Replicas_number):
    folder = 'output_R'+str("%02d"%replica_id)
    file.write("folder"+" "+folder+"\n")
    file.write("#"+"\n")
    #os.mkdir(folder)
file.close()

xml_file_in = 'config/PhysiCell_settings.xml'
xml_file_out = 'config/tmp.xml'
copyfile(xml_file_in,xml_file_out)
tree = ET.parse(xml_file_out)
xml_root = tree.getroot()
first_time = True
output_dirs = []
with open(fileOut) as f:
    for line in f:
        print(len(line),line)
        if (line[0] == '#'):
            continue
        (key, val) = line.split()
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

tree.write(xml_file_out)
print(output_dirs)