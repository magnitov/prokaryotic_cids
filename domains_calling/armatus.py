import argparse
import subprocess
import os
import pandas as pd
import numpy as np
import cooler

# Function that returns command to run armatus
def run_armatus(replicate, out_dir, gamma, resolution):
    armatus_path = '/home/magnitov/Software/armatus-2.2/bin/armatus'
    input_matrix = '/home/magnitov/analysis/MICROBES/' + out_dir + '/formatted_maps/tab/' + replicate + '.tab.gz'
    out_pref = replicate + '.TADs.' + str(gamma)
    return('{} --gammaMax {} -j -r {} -c 1 --input {} --output {}'.format(armatus_path, gamma, resolution, input_matrix, out_pref))

# --------------------------- MAIN ----------------------------------- #
# Parsing input arguments
parser = argparse.ArgumentParser()
parser.add_argument('tab_path', type = str, help = 'File name')
parser.add_argument('resolution', type = int, help = 'Resolution')
args = parser.parse_args()

# Path to the input files
tab_path = args.tab_path
resolution = args.resolution

# Create folders to store output files
out_dir = tab_path.split('/')[1]
replicate = tab_path.split('/')[-1].replace('.tab.gz', '')
if not os.path.exists('./' + out_dir + '/TADs/armatus/'):
    os.mkdir('./' + out_dir + '/TADs/armatus/')

# Run TAD caller
os.chdir('./' + out_dir + '/TADs/armatus/')
for gamma in np.arange(0.01, 5.01, 0.01):
    gamma = round(gamma, 2)
    # Call TADs
    print('Calling TADs at gamma=' + str(gamma))
    process = subprocess.Popen(run_armatus(replicate, out_dir, gamma, resolution), shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = process.communicate()
    # Filter small TADs
    if os.stat(replicate + '.TADs.' + str(gamma) + '.consensus.txt').st_size != 0:
        tads = pd.read_csv(replicate + '.TADs.' + str(gamma) + '.consensus.txt', sep = '\t', header = None)
        tads[2] = tads[2].values+1
        index_keep = []
        for i in range(0, len(tads)):
            if tads[2].values[i]-tads[1].values[i] > 3*resolution:
                index_keep.append(i)
        tads = tads.iloc[index_keep]
        tads.to_csv(replicate + '.TADs.' + str(gamma) + '.consensus.txt', sep = '\t', header = 0, index = 0)
