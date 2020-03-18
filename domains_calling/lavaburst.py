import lavaburst
import os
import numpy as np
import pandas as pd
import argparse

def call_tads(matrix, gamma, method, resolution):
    # Create TADs caller and do segmentation
    if method == 'armatus':
        caller = lavaburst.scoring.armatus_score(matrix, gamma = gamma)
    elif method == 'corner':
        caller = lavaburst.scoring.corner_score(matrix, gamma = gamma)
    elif method == 'modularity':
        caller = lavaburst.scoring.modularity_score(matrix, gamma = gamma)
    elif method == 'variance':
        caller = lavaburst.scoring.variance_score(matrix, gamma = gamma)
    else:
        print('Incorrect method!')
    model = lavaburst.SegModel(caller)
    segments = model.optimal_segmentation()
    segments = [seg for seg in segments if (seg[1]-seg[0])*resolution > 3*resolution]
    return(segments)

# Parsing input arguments
parser = argparse.ArgumentParser()
parser.add_argument('tab_path', type = str, help = 'File name')
parser.add_argument('resolution', type = int, help = 'Resolution')
args = parser.parse_args()

# Path to the input files
tab_path = args.tab_path
resolution = args.resolution
matrix = pd.read_table(tab_path, sep = '\t', header = None).as_matrix()

# Create folders to store output files
out_dir = tab_path.split('/')[1]
replicate = tab_path.split('/')[-1].replace('.tab', '')
if not os.path.exists('./' + out_dir + '/TADs/lavaburst/'):
    os.mkdir('./' + out_dir + '/TADs/lavaburst/')

# Run TAD caller
for gamma in np.arange(0.01, 5.01, 0.01):
    gamma = round(gamma, 2)
    print('Calling TADs at gamma=' + str(gamma))
    for scoring_function in ['armatus', 'corner', 'modularity', 'variance']:
        print('\t', scoring_function)
        segments = call_tads(matrix, gamma, scoring_function, resolution)
        tads = pd.DataFrame({'Chromosome' : ['1']*len(segments), 
                             'Start' : [seg[0]*resolution for seg in segments], 
                             'End' : [seg[1]*resolution for seg in segments]})
        tads.to_csv('./' + out_dir + '/TADs/lavaburst/' + replicate + '.' + scoring_function + '.' + str(gamma) + '', 
                    sep = '\t', header = 0, index = 0)

