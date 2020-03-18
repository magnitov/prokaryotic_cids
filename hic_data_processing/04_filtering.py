import argparse
import os
import pandas as pd
import logging
logging.basicConfig(level = logging.DEBUG)
from hiclib import fragmentHiC, binnedData
from mirnylib import h5dict, genome
import cooler

def filtration(chromosome_names, cell_line, path, genome_version, enzyme, resolution_list):
    for chrm_list in chromosome_names:
        genome_db = genome.Genome('/home/magnitov/data/genomes/' + genome_version, gapFile = 'gap.txt', readChrms = chrm_list, forceOrder = True)
        # Read mapped reads        
        if len(chrm_list) > 1:
            fragments = fragmentHiC.HiCdataset(
                filename = path + 'filtered_maps/' + cell_line + '/fragment_dataset_full.hdf5',
                genome = genome_db,
                enzymeName = enzyme,
                mode = 'w')
            fragments.parseInputData(dictLike = path + 'maps/' + cell_line + '/mapped_reads_full.hdf5')
        else:
            fragments = fragmentHiC.HiCdataset(
                filename = path + 'filtered_maps/' + cell_line + '/fragment_dataset_' + chrm_list[0] + '.hdf5',
                genome = genome_db,
                enzymeName = enzyme,
                mode = 'w')
            fragments.parseInputData(dictLike = path + 'maps/' + cell_line + '/mapped_reads_' + chrm_list[0] + '.hdf5')
        # Apply filters
        fragments.filterDuplicates()
        # Save statistics
        if len(chrm_list) > 1 or len(chromosome_names) == 1:
            fragments.writeFilteringStats()
            fragments.printMetadata(saveTo = path + 'processing_stats/' + cell_line + '/processing_stats_' + cell_line + '.txt')
        if len(chrm_list) == 1 and len(chromosome_names) > 1:
            fragments.writeFilteringStats()
            fragments.printMetadata(saveTo = path + 'processing_stats/' + cell_line + '/processing_stats_' + cell_line + '_' + chrm_list[0] + '.txt')

        # Sort reads and calculate contact probability (both normalized and not)
        fragments._sortData()
        if len(chrm_list) > 1:
            contact_probs = fragments.plotScaling(normalize = True, plot = False)
            pd.DataFrame({'Distance' : contact_probs[0], 'Probability' : contact_probs[1]}).to_csv(path + 'contact_probs/' + cell_line + '/contact_probs_' + cell_line + '_full_norm.txt', header = 1, index = 0, sep = '\t')
            contact_probs = fragments.plotScaling(normalize = False, plot = False)
            pd.DataFrame({'Distance' : contact_probs[0], 'Probability' : contact_probs[1]}).to_csv(path + 'contact_probs/' + cell_line + '/contact_probs_' + cell_line + '_full.txt', header = 1, index = 0, sep = '\t')

        if len(chrm_list) == 1:
            contact_probs = fragments.plotScaling(normalize = True, plot = False)
            pd.DataFrame({'Distance' : contact_probs[0], 'Probability' : contact_probs[1]}).to_csv(path + 'contact_probs/' + cell_line + '/contact_probs_' + cell_line + '_' + chrm_list[0] + '_norm.txt', header = 1, index = 0, sep = '\t')
            contact_probs = fragments.plotScaling(normalize = False, plot = False)
            pd.DataFrame({'Distance' : contact_probs[0], 'Probability' : contact_probs[1]}).to_csv(path + 'contact_probs/' + cell_line + '/contact_probs_' + cell_line + '_' + chrm_list[0] + '.txt', header = 1, index = 0, sep = '\t')

        # Save into .cool and .hdf5 files
        for resolution in resolution_list:
            fragments.saveCooler(filename = path + 'filtered_maps/' + cell_line + '/heatmap-' + chrm_list[0] + '-' + str(resolution/1000) + 'K.cool', resolution = resolution)

# --------------------------- MAIN ----------------------------------- #

# Parsing input arguments
parser = argparse.ArgumentParser()
parser.add_argument('path', type = str, help = 'Path to the files')
parser.add_argument('filename', type = str, help = 'File name')
parser.add_argument('genome', type = str, help = 'Genome version')
parser.add_argument('enzyme', type = str, help = 'Restriction enzyme')
parser.add_argument('resolution_list', type = str, help = 'List of map resolutions')
args = parser.parse_args()

# Path to the input files, genome_version and restriction enzyme
path = args.path
cell_line = args.filename
genome_version = args.genome
enzyme = args.enzyme
resolution_list = [int(x) for x in args.resolution_list.split(',')]

# Create folders to store output files
if not os.path.exists(path + 'processing_stats/'):
    os.mkdir(path + 'processing_stats/')
if not os.path.exists(path + 'processing_stats/' + cell_line + '/'):
        os.mkdir(path + 'processing_stats/' + cell_line + '/')
if not os.path.exists(path + 'contact_probs/'):
    os.mkdir(path + 'contact_probs/')
if not os.path.exists(path + 'contact_probs/' + cell_line + '/'):
        os.mkdir(path + 'contact_probs/' + cell_line + '/')
if not os.path.exists(path + 'filtered_maps/'):
    os.mkdir(path + 'filtered_maps/')
if not os.path.exists(path + 'filtered_maps/' + cell_line + '/'):
        os.mkdir(path + 'filtered_maps/' + cell_line + '/')

# Get chromosome names
chrom_info = pd.read_table('/home/magnitov/data/genomes/' + genome_version + '_chromsizes.txt', header = None, sep = ',')
chromosome_names = [str(name).replace('chr', '') for name in chrom_info[0].values]
if len(chromosome_names) > 1:
    chromosome_names = chromosome_names + [chromosome_names]

# Filter fragments
filtration(chromosome_names, cell_line, path, genome_version, enzyme, resolution_list)
