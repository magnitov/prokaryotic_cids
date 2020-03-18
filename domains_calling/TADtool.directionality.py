import numpy as np
import pandas as pd
import os
import argparse

# Parsing input arguments
parser = argparse.ArgumentParser()
parser.add_argument('path_input', type = str, help = 'Input maps', default = '../RESOLUTION/Ecoli_20M_10kb_rep2.tab')
parser.add_argument('path_input_bed', type = str, help = 'Input bed', default = '../RESOLUTION/ec1_10kb_windows.bed')
parser.add_argument('path_output', type = str, help = 'Output folder for bed files', default = '../RESOLUTION/TADtool/test_direct_Ecoli_20M_10kb_rep1')
parser.add_argument('resolution', type = int, help = 'Resolution of map')
args = parser.parse_args()

path_input = args.path_input
path_input_bed = args.path_input_bed
path_output = args.path_output
resolution = args.resolution
b = []
for i in np.arange(0,0.0201,0.0001):
    b.append(round(i,4))
for ws in np.arange(2, 21)*resolution:
    for cutoff in b:
        os.system(' '.join([
            'tadtool tads',
            path_input,
            path_input_bed,
            '-a directionality',
            str(ws),
            str(cutoff),
            path_output + '_ws_' + str(ws) + '_co_' + str(cutoff)
        ]))


