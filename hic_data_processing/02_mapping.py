import argparse
import os
import logging
from hiclib import mapping
from mirnylib import h5dict, genome
logging.basicConfig(level = logging.DEBUG)

# Function to perform mapping of Hi-C reads to reference genome
def map_reads(cell_line, path_input, path_output, genome_version):
    if not os.path.exists(path_output + 'bam/' + cell_line):
        os.mkdir(path_output + 'bam/' + cell_line)

    mapping.iterative_mapping(
        bowtie_path = '/usr/bin/bowtie2',
        bowtie_index_path = path_output + 'index_' + cell_line + '/' + genome_version,
        fastq_path = path_input + cell_line + '_R1.fastq.gz',
        out_sam_path = path_output + 'bam/' + cell_line + '/' + cell_line + '_R1.bam',
        min_seq_len = 25,
        seq_start = 4,
        len_step = 3,
        nthreads = 8,
        temp_dir = path_output + 'tmp_' + cell_line,
        bowtie_flags = '--very-sensitive')

    mapping.iterative_mapping(
        bowtie_path = '/usr/bin/bowtie2',
        bowtie_index_path = path_output + 'index_' + cell_line + '/' + genome_version,
        fastq_path = path_input + cell_line + '_R2.fastq.gz',
        out_sam_path = path_output + 'bam/' + cell_line + '/' + cell_line + '_R2.bam',
        min_seq_len = 25,
        seq_start = 4,
        len_step = 3,
        nthreads = 8,
        temp_dir = path_output + 'tmp_' + cell_line,
        bowtie_flags = '--very-sensitive')

# ----------- MAIN ------------------------------------- #

# Parsing input arguments
parser = argparse.ArgumentParser()
parser.add_argument('path_input', type = str, help = 'Path to the input files')
parser.add_argument('path_output', type = str, help = 'Path to the output files')
parser.add_argument('filename', type = str, help = 'File name')
parser.add_argument('genome', type = str, help = 'Genome version')
args = parser.parse_args()

# Input and output folders, file and genome version
path_input = args.path_input
path_output = args.path_output
cell_line = args.filename
genome_version = args.genome

# Create folders to store temporary and output files
if not os.path.exists(path_output):
    os.makedirs(path_output)
if not os.path.exists(path_output + 'tmp_' + cell_line + '/'):
    os.mkdir(path_output + 'tmp_' + cell_line + '/')
if not os.path.exists(path_output + 'bam/'):
    os.mkdir(path_output + 'bam/')

# Index genome with Bowtie2
if not os.path.exists(path_output + 'index_' + cell_line + '/'):
    os.mkdir(path_output + 'index_' + cell_line + '/')
os.system('bowtie2-build /home/magnitov/data/genomes/' + genome_version + '.fa ' + path_output + 'index_' + cell_line + '/' + genome_version)

# Mapping
map_reads(cell_line, path_input, path_output, genome_version)

# Remove temporary folder we don't need any more
os.system('rm -r ' + path_output + 'tmp_' + cell_line + '/')
os.system('rm -r ' + path_output + 'index_' + cell_line + '/')
