import os
import argparse

# Parsing input arguments
parser = argparse.ArgumentParser()
parser.add_argument('path_input', type = str, help = 'Path to the input files')
parser.add_argument('path_output', type = str, help = 'Path to the output files')
args = parser.parse_args()

# Input and output folders
path_input = args.path_input
path_output = args.path_output
if not os.path.exists(path_output):
    os.mkdir(path_output)
if not os.path.exists(path_output + 'quality/'):
    os.mkdir(path_output + 'quality/')

# Get names of all fastq files from the input folder
files = [f for f in os.listdir(path_input) if '.fastq.gz' in f]

# Run FastQC
os.chdir(path_input)
os.system('/home/magnitov/Software/FastQC/./fastqc --extract -o ' + path_output + 'quality/ ' + ' '.join(files))
os.chdir(path_output + 'quality/')
os.system('rm *.zip')
os.system('multiqc --filename ' + path_input.split('/')[-2] + '_QC .')
