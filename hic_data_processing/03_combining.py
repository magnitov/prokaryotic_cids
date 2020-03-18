import argparse
import os
import logging
import pandas as pd
from hiclib import mapping
from mirnylib import h5dict, genome
logging.basicConfig(level = logging.DEBUG)


def parse_bams(chromosome_names, cell_line, path, genome_version, enzyme):

    if not os.path.exists(path + 'maps/' + cell_line):
        os.mkdir(path + 'maps/' + cell_line)

    for chrm_list in chromosome_names:

        if len(chrm_list) > 1:
            mapped_reads = h5dict.h5dict(path + 'maps/' + cell_line +  '/mapped_reads_full.hdf5')
        else:
            mapped_reads = h5dict.h5dict(path + 'maps/' + cell_line +  '/mapped_reads_' + chrm_list[0] + '.hdf5')
        
        genome_db = genome.Genome('/home/magnitov/data/genomes/' + genome_version, gapFile = 'gap.txt' , readChrms = chrm_list, forceOrder = True)

        mapping.parse_sam(
            sam_basename1 = path + 'bam/' + cell_line + '/' + cell_line + '_R1.bam',
            sam_basename2 = path + 'bam/' + cell_line + '/' + cell_line + '_R2.bam',
            out_dict = mapped_reads,
            genome_db = genome_db,
            enzyme_name = enzyme)


# ----------- MAIN ------------------------------------- #

# Parsing input arguments
parser = argparse.ArgumentParser()
parser.add_argument('path', type = str, help = 'Path to the files')
parser.add_argument('filename', type = str, help = 'File name')
parser.add_argument('genome', type = str, help = 'Genome version')
parser.add_argument('enzyme', type = str, help = 'Restriction enzyme')
args = parser.parse_args()

# Path to the input files, genome version and restriction enzyme
path = args.path
cell_line = args.filename
genome_version = args.genome
enzyme = args.enzyme

# Create folders to store output files
if not os.path.exists(path + 'maps/'):
    os.mkdir(path + 'maps/')

# Get chromosome names
chrom_info = pd.read_table('/home/magnitov/data/genomes/' + genome_version + '_chromsizes.txt', header = None, sep = ',')
chromosome_names = [str(name).replace('chr', '') for name in chrom_info[0].values]
if len(chromosome_names) > 1:
    chromosome_names = chromosome_names + [chromosome_names]
# Parsing
parse_bams(chromosome_names, cell_line, path, genome_version, enzyme)
