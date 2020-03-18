import argparse
import os
import cooler
import pandas as pd
import numpy as np

# Parse file name
parser = argparse.ArgumentParser()
parser.add_argument('path', type = str, help = 'path to files')
parser.add_argument('cool_file', type = str, help = 'cooler to convert')
parser.add_argument('genome', type = str, help = 'genome')
args = parser.parse_args()

# Read arguments
f = args.cool_file
path = args.path
genome = args.genome

# Create folder for maps in different formats
if not os.path.exists(path + 'formatted_maps/'):
    os.mkdir(path + 'formatted_maps/')

# GZIPPED AND TAB SEPARATED
# Read matrix and metadata
replicate = f.split('/')[-1].replace('.cool', '')
c = cooler.Cooler(f)
m = c.matrix(balance = True)[:]
m = pd.DataFrame(m).interpolate(method = 'linear', axis = 0).as_matrix()
m = pd.DataFrame(m).interpolate(method = 'linear', axis = 1).as_matrix()
resolution = c.info['bin-size']

# Create tab separated map
if not os.path.exists(path + 'formatted_maps/tab/'):
    os.mkdir(path + 'formatted_maps/tab/')
pd.DataFrame(m).to_csv(path + 'formatted_maps/tab/' + replicate + '.tab', sep = '\t', index = 0, header = 0)
os.system('gzip ' + path + 'formatted_maps/tab/' + replicate + '.tab')
pd.DataFrame(m).to_csv(path + 'formatted_maps/tab/' + replicate + '.tab', sep = '\t', index = 0, header = 0)

# RAW TAB SEPARATED
# Read matrix and metadata
replicate = f.split('/')[-1].replace('.cool', '')
c = cooler.Cooler(f)
m = c.matrix(balance = False)[:]
m = pd.DataFrame(m).interpolate(method = 'linear', axis = 0).as_matrix()
m = pd.DataFrame(m).interpolate(method = 'linear', axis = 1).as_matrix()
resolution = c.info['bin-size']

# Create tab separated map
if not os.path.exists(path + 'formatted_maps/tab/'):
    os.mkdir(path + 'formatted_maps/tab/')
pd.DataFrame(m).to_csv(path + 'formatted_maps/tab/' + replicate + '_raw.tab', sep = '\t', index = 0, header = 0)

# HIC
# Create .hic map
if not os.path.exists(path + 'formatted_maps/hic/'):
    os.mkdir(path + 'formatted_maps/hic/')
with open(path + 'formatted_maps/hic/' + replicate + '.pre', 'w') as f:
    start1 = 0
    for i in range(0, len(m)):
        start2 = 0
        for j in range(0, len(m[i])):
            f.write(' '.join(['0', '1', str(start1), '0', '0', '1', str(start2), '1', str(m[i][j])]) + '\n')
            start2+=resolution
        start1+=resolution
os.system('java -jar /home/magnitov/Software/juicer_tools_1.11.09_jcuda.0.8.jar pre -r ' + str(resolution) + ' -n ' + path + 'formatted_maps/hic/' + replicate + '.pre ' + path + 'formatted_maps/hic/' + replicate + '.hic ' + path + genome + '_chromsizes.txt')