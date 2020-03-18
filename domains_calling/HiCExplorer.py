import argparse
import subprocess
import os
import pandas as pd
import numpy as np
import cooler

# Function that returns command to run armatus
def run_hicexplorer(replicate, out_dir, ws, delta):
    input_matrix = '/home/magnitov/analysis/MICROBES/' + out_dir + '/formatted_maps/h5/' + replicate + '.h5'
    out_pref = replicate + '.ws_' + str(ws) + '.delta_' + str(delta)
    return('hicFindTADs --matrix {} --outPrefix {} --minDepth {} --delta {} --correctForMultipleTesting fdr'.format(input_matrix, out_pref, ws, delta))

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
replicate = tab_path.split('/')[-1].replace('.h5', '')
if not os.path.exists('./' + out_dir + '/TADs/HiCExplorer/'):
    os.mkdir('./' + out_dir + '/TADs/HiCExplorer/')

# Run TAD caller
os.chdir('./' + out_dir + '/TADs/HiCExplorer/')

for ws in np.arange(3*resolution, 21*resolution, resolution):
    for delta in np.arange(0.005, 0.105, 0.005):
        delta = np.round(delta, 3)
        # Call TADs
        print('Calling TADs at ws=' + str(ws) + ', delta=' + str(delta))
        process = subprocess.Popen(run_hicexplorer(replicate, out_dir, ws, delta), shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = process.communicate()
        # Filter small TADs
        if os.stat(replicate + '.ws_' + str(ws) + '.delta_' + str(delta) + '_domains.bed').st_size != 0:
            tads = pd.read_csv(replicate + '.ws_' + str(ws) + '.delta_' + str(delta) + '_domains.bed', sep = '\t', header = None)
            index_keep = []
            for i in range(0, len(tads)):
                if tads[2].values[i]-tads[1].values[i] > 3*resolution:
                    index_keep.append(i)
            tads = tads.iloc[index_keep][[0, 1, 2]]
            tads.to_csv(replicate + '.ws_' + str(ws) + '.delta_' + str(delta) + '_domains.bed', sep = '\t', header = 0, index = 0)

        os.system('rm *.bm')

os.system('rm *.gff *.h5 *.bedgraph *boundaries.bed')
