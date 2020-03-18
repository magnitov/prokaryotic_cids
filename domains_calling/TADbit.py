import pandas as pd
import os
import argparse
from pytadbit import Chromosome

# Parsing input arguments
parser = argparse.ArgumentParser()
parser.add_argument('path_input_raw', type = str, help = 'Input maps raw', default = 'RESOLUTION/Ecoli_20M_3kb_rep2_raw.tab')
parser.add_argument('path_input_norm', type = str, help = 'Input maps normalized', default = 'RESOLUTION/Ecoli_20M_3kb_rep2.tab')
parser.add_argument('path_output', type = str, help = 'Output folder for bed files', default = 'RESOLUTION/Tadbit/Ecoli_20M_3kb_rep1_final')
parser.add_argument('restictase', type = str, help = 'Restrictase', default = 'HpaII')
parser.add_argument('resolution', type = int, help = 'Resolution of map')
args = parser.parse_args()

path_input_raw = args.path_input_raw
path_input_norm = args.path_input_norm
path_output = args.path_output
restrictase = args.restictase
resolution = args.resolution

my_chrom = Chromosome(name = '1',centromere_search = False)
my_chrom.add_experiment(restrictase + '1_stat', resolution = resolution, hic_data = path_input_raw,\
                        norm_data = path_input_norm, enzyme = restrictase)
exp = my_chrom.experiments[restrictase + '1_stat']
my_chrom.find_tad([restrictase + '1_stat'], verbose = True, batch_mode = False)
my_chrom.experiments[restrictase + '1_stat']
exp.write_tad_borders(density = True, savedata = 'tmp.txt', normalized = False)
data = pd.read_csv('tmp.txt', sep = '\t')
data['start'] -= 1
data *= resolution
data = data[ data['end'] - data['start'] >= resolution * 4 ]
data['ix'] = 'chr1'
data = data[['ix', 'start', 'end']]
data.to_csv(path_output, sep = '\t', header = False, index = False)
